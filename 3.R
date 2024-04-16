library("infercnv")
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

PCSC <- readRDS("./2.1.harmony_res2/PCSC_annoted.Rds")

# doublet ------------------------------------------------------------

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

PCSC_split <- SplitObject(PCSC, split.by = "orig.ident")

for (i in 1:length(PCSC_split)) {
  # pK Identification
  sweep.res.list <- paramSweep_v3(PCSC_split[[i]], PCs = 1:10, sct = FALSE)
  
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- PCSC_split[[i]]@meta.data$celtype
  homotypic.prop <- modelHomotypic(annotations)  
  DoubletRate = ncol(PCSC_split[[i]])*8*1e-6     #按每增加1000个细胞，双细胞比率增加千分之8来计算
  
  # 估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。
  nExp_poi <- round(DoubletRate*nrow(PCSC_split[[i]]@meta.data))  
  
  # 计算双细胞比例
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # doubletFinder
  PCSC_split[[i]] <- doubletFinder_v3(PCSC_split[[i]], PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  PCSC_split[[i]] <- doubletFinder_v3(PCSC_split[[i]], PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  # judge
  PCSC_split[[i]]@meta.data[,"Doublet"] <- PCSC_split[[i]]@meta.data[,length(colnames(PCSC_split[[i]]@meta.data))]

}

#merge seurat
PCSC$doublet <- "Singlet"
for (i in 1:length(PCSC_split)) {
  PCSC@meta.data$doublet[rownames(PCSC@meta.data) %in% 
                           rownames(PCSC_split[[i]]@meta.data)[PCSC_split[[i]]@meta.data$Doublet == "Doublet"]] <- "Doublet"
}

PCSC <- subset(PCSC, doublet != "Doublet")

saveRDS(PCSC, "./4.3.Doublet/PCSC_doublet_Removed.Rds")

pdf("./4.3.Doublet/DimPlot_doublet.pdf", width = 8.5, height = 7)
DimPlot(PCSC, group.by="doublet") + scale_color_npg()
dev.off()

# inferCNV ------------------------------------------

folder <- "4.1.inferCNV"
if(!dir.exists(folder)){
  dir.create(folder)
}

# Immune as reference ------------------------------------------------

PCSC <- readRDS("./2.1.harmony/PCSC_annoted.Rds")
PCSC@meta.data$cellGroup <- "Immune_cells"
PCSC@meta.data$cellGroup[ PCSC@meta.data$celltype %in% c("Endothelial", "Mesenchymal")] <- "Stromal_cells"
PCSC@meta.data$cellGroup[ PCSC@meta.data$celltype == "Epithelial" ] <- "Epithelial_cells"
Idents(PCSC) <- "cellGroup"
PCSC_inf <- subset(PCSC, cellGroup %in% c("Immune_cells", "Epithelial_cells"))

folder <- "4.1.inferCNV/Imm_Epi"
if(!dir.exists(folder)){
  dir.create(folder)
}

##### Create Infercnv Object #####

load("./geneordering.RData")
list <- rownames(table(PCSC@meta.data$orig.ident))

folder <- "4.1.inferCNV/Imm_Epi/objective"
if(!dir.exists(folder)){
  dir.create(folder)
}

for (i in list){
  sam.name <- paste0("infer_",i)
  rtmp <- subset(PCSC_inf, subset = orig.ident == i)
  matrixtmp <- as.matrix(rtmp@assays$RNA@counts)
  annotmp <- as.matrix(rtmp@active.ident)
  inferobject <- CreateInfercnvObject(raw_counts_matrix = matrixtmp,
                                      annotations_file = annotmp,
                                      gene_order_file = geneOrderingfile,
                                      ref_group_names = c("Immune_cells"))
  save(inferobject,file=paste0(folder, "/",sam.name,".RData"))
  rm(rtmp,matrixtmp,inferobject)
}

##### Infer CNV #####

folderresult <- "4.1.inferCNV/Imm_Epi/inferResult"
if(!dir.exists(folderresult)){
  dir.create(folderresult)
}

for (i in list){
  sam.name <- paste0("infer_",i)
  load(file=paste0(folder, "/",sam.name,".RData"))
  foldertmp <- paste0(folderresult, "/", "infer_", i)
  if(!dir.exists(foldertmp)){
    dir.create(foldertmp)
  }
  inferResult = infercnv::run(inferobject,
                              cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir=foldertmp,  # 输出文件夹
                              cluster_by_groups=T,   # 聚类
                              denoise=T,    #去噪
                              HMM=T,
                              num_threads = 20)
  save(inferResult,file=paste0(foldertmp, "/", sam.name, "_result.RData"))
  rm(inferobject,inferResult)
}

##### Infer cancer cells #####

for (i in list) {
  
  infercnv_obj = readRDS(paste0("/data_analysis_2/csl/cyf/SC/4.1.inferCNV/Imm_Epi/inferResult/infer_", i, "/run.final.infercnv_obj"))
  expr <- infercnv_obj@expr.data
  
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc <- normal_loc$Immune_cells
  
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  test_loc <- test_loc$Epithelial_cells
  
  anno.df=data.frame(
    CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
    class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
  )
  head(anno.df)
  
  gn <- rownames(expr)
  geneFile <- geneOrderingfile
  sub_geneFile <- geneFile[intersect(gn,rownames(geneFile)),]
  expr=expr[intersect(gn,rownames(geneFile)),]
  head(sub_geneFile,4)
  expr[1:4,1:4]
  
  #聚类，8类，提取结果
  set.seed(20220624)
  kmeans.result <- kmeans(t(expr), 8)
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  kmeans_df$CB <- rownames(kmeans_df)
  kmeans_df <- kmeans_df %>% inner_join(anno.df, by="CB")   #合并
  kmeans_df_s <- arrange(kmeans_df, kmeans_class)   #排序
  rownames(kmeans_df_s) <- kmeans_df_s$CB
  kmeans_df_s$CB <- NULL
  kmeans_df_s$kmeans_class <- as.factor(kmeans_df_s$kmeans_class)   #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:8排序
  head(kmeans_df_s)
  
  #定义热图的注释，及配色
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col = "NA"), labels = 1:22, labels_gp = gpar(cex = 1.5)))
  color_v <- RColorBrewer::brewer.pal(8, "Dark2")[1:8] #类别数
  names(color_v) <- as.character(1:8)
  left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="#A52A2A","normal" = "#053E7A"),kmeans_class=color_v))
  
  #####  1. heatmap  #####
  #绘图
  pdf(paste0("/data_analysis_2/csl/cyf/SC/4.1.inferCNV/Imm_Epi/inferResult/infer_", i, "/inferCancer_heatmap_N.pdf"),width = 7.5,height = 5)
  ht <-  Heatmap(t(expr)[rownames(kmeans_df_s),],   #绘图数据的CB顺序和注释CB顺序保持一致
                 col = colorRamp2(c(0.8,1,1.2), c("#053E7A","#F0F0F0","#A52A2A")),  #这里的刻度会有所变化
                 cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = F,
                 column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")),  #控制染色体顺序，即使你的基因排序文件顺序是错的
                 column_gap = unit(2, "mm"),
                 heatmap_legend_param = list(title = "Modified expression", direction = "vertical", 
                                             title_position = "leftcenter-rot", at=c(0.8,1,1.2), 
                                             legend_height = unit(3, "cm")),
                 top_annotation = top_anno, left_annotation = left_anno, #添加注释
                 row_title = NULL,column_title = NULL)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  #每一类对应的CB保存在kmeans_df_s数据框中
  write.table(kmeans_df_s, file = paste0("/data_analysis_2/csl/cyf/SC/4.1.inferCNV/Imm_Epi/inferResult/infer_", i, "/kmeans_df_s.txt"), quote = FALSE, sep = '\t', row.names = T, col.names = T)
  
  #####  1. CNV score  #####
  expr2 <- expr-1
  expr2 <- expr2 ^ 2
  CNV_score <- as.data.frame(colMeans(expr2))
  colnames(CNV_score) <- "CNV_score"
  CNV_score$CB <- rownames(CNV_score)
  kmeans_df_s$CB <- rownames(kmeans_df_s)
  CNV_score <- CNV_score %>% inner_join(kmeans_df_s,by="CB")
  
  CNV_score %>% ggplot(aes(kmeans_class,CNV_score)) + geom_violin(aes(fill=kmeans_class),color="NA") +
    scale_fill_manual(values = color) +
    theme_bw()
  ggsave(paste0("/data_analysis_2/csl/cyf/SC/4.1.inferCNV/Imm_Epi/inferResult/infer_", i, "/vlnPlot_CNV_level_N.pdf"),width = 20,height = 12,units = "cm")
}

#LOPC6

PCSC_inf <- subset(PCSC, cellGroup %in% c("Immune_cells"))
LOPC1 <- subset(PCSC_inf, subset = orig.ident == "LOPC1")
LOPC2 <- subset(PCSC_inf, subset = orig.ident == "LOPC2")
LOPC3 <- subset(PCSC_inf, subset = orig.ident == "LOPC3")
LOPC4 <- subset(PCSC_inf, subset = orig.ident == "LOPC4")
LOPC5 <- subset(PCSC_inf, subset = orig.ident == "LOPC5")

matrix1 <- as.matrix(LOPC1@assays$RNA@counts)
matrix2 <- as.matrix(LOPC1@assays$RNA@counts)
matrix3 <- as.matrix(LOPC1@assays$RNA@counts)
matrix4 <- as.matrix(LOPC1@assays$RNA@counts)
matrix5 <- as.matrix(LOPC1@assays$RNA@counts)

matrixtmp <- cbind(matrix1, matrix2, matrix3, matrix4, matrix5)
annotmp <- matrix(nrow = ncol(matrixtmp), ncol = 1)
rownames(annotmp) <- colnames(matrixtmp)
annotmp[,1] <- "Immune_cells_ref"

LOPC6 <- subset(PCSC, subset = orig.ident == "LOPC6")
matrix6 <- as.matrix(LOPC6@assays$RNA@counts)
anno <- as.matrix(LOPC6@active.ident)
matrixtmp <- cbind(matrixtmp, matrix6)
annotmp <- rbind(annotmp, anno)

i="LOPC6"
sam.name <- paste0("infer_",i)
inferobject <- CreateInfercnvObject(raw_counts_matrix = matrixtmp,
                                    annotations_file = annotmp,
                                    gene_order_file = geneOrderingfile,
                                    ref_group_names = c("Immune_cells_ref"))
save(inferobject,file=paste0(folder, "/",sam.name,".RData"))

load(file=paste0(folder, "/",sam.name,".RData"))
foldertmp <- paste0(folderresult, "/", "infer_", i)
inferResult = infercnv::run(inferobject,
                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                            out_dir=foldertmp,  # 输出文件夹
                            cluster_by_groups=T,   # 聚类
                            denoise=T,    #去噪
                            HMM=T,
                            num_threads = 20)
save(inferResult, file=paste0(foldertmp, "/", sam.name, "_result.RData"))
rm(inferobject,inferResult)

#infer cancer cell
infercnv_obj = readRDS(paste0("/data_analysis_2/csl/cyf/SC/4.1.inferCNV/Imm_Epi/inferResult/infer_", i, "/run.final.infercnv_obj"))
expr <- infercnv_obj@expr.data

loc <- infercnv_obj@observation_grouped_cell_indices
normal_loc <- c(loc$Immune_cells, loc$Stromal_cells)
test_loc <- loc$Epithelial_cells

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)

gn <- rownames(expr)
geneFile <- geneOrderingfile
sub_geneFile <- geneFile[intersect(gn,rownames(geneFile)),]
expr=expr[intersect(gn,rownames(geneFile)),]
head(sub_geneFile,4)
expr[1:4,1:4]

#聚类，8类，提取结果
set.seed(20220624)
kmeans.result <- kmeans(t(expr), 8)
kmeans_df <- data.frame(as.matrix(kmeans.result$cluster))
colnames(kmeans_df) <- "kmeans_class"
kmeans_df$CB <- rownames(kmeans_df)
kmeans_df <- kmeans_df %>% inner_join(anno.df, by="CB")   #合并
kmeans_df_s <- arrange(kmeans_df, kmeans_class)   #排序
rownames(kmeans_df_s) <- kmeans_df_s$CB
kmeans_df_s$CB <- NULL
kmeans_df_s$kmeans_class <- as.factor(kmeans_df_s$kmeans_class)   #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:8排序
head(kmeans_df_s)

#定义热图的注释，及配色
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col = "NA"), labels = 1:22, labels_gp = gpar(cex = 1.5)))
color_v <- RColorBrewer::brewer.pal(8, "Dark2")[1:8] #类别数
names(color_v) <- as.character(1:8)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="#A52A2A","normal" = "#053E7A"),kmeans_class=color_v))

#####  1. heatmap  #####
#绘图
pdf(paste0("/data_analysis_2/csl/cyf/SC/4.1.inferCNV/Imm_Epi/inferResult/infer_", i, "/inferCancer_heatmap_N.pdf"),width = 7.5,height = 5)
ht <-  Heatmap(t(expr)[rownames(kmeans_df_s),],   #绘图数据的CB顺序和注释CB顺序保持一致
               col = colorRamp2(c(0.8,1,1.2), c("#053E7A","#F0F0F0","#A52A2A")),  #这里的刻度会有所变化
               cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = F,
               column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")),  #控制染色体顺序，即使你的基因排序文件顺序是错的
               column_gap = unit(2, "mm"),
               heatmap_legend_param = list(title = "Modified expression", direction = "vertical", 
                                           title_position = "leftcenter-rot", at=c(0.8,1,1.2), 
                                           legend_height = unit(3, "cm")),
               top_annotation = top_anno, left_annotation = left_anno, #添加注释
               row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

#每一类对应的CB保存在kmeans_df_s数据框中
write.table(kmeans_df_s, file = paste0("/data_analysis_2/csl/cyf/SC/4.1.inferCNV/Imm_Epi/inferResult/infer_", i, "/kmeans_df_s.txt"), quote = FALSE, sep = '\t', row.names = T, col.names = T)

#####  1. CNV score  #####
expr2 <- expr-1
expr2 <- expr2 ^ 2
CNV_score <- as.data.frame(as.matrix(colMeans(expr2)))
colnames(CNV_score) <- "CNV_score"
CNV_score$CB <- rownames(CNV_score)
kmeans_df_s$CB <- rownames(kmeans_df_s)
CNV_score <- CNV_score %>% inner_join(kmeans_df_s,by="CB")

CNV_score %>% ggplot(aes(kmeans_class,CNV_score)) + geom_violin(aes(fill=kmeans_class),color="NA") +
  scale_fill_manual(values = color) +
  theme_bw()
ggsave(paste0("/data_analysis_2/csl/cyf/SC/4.1.inferCNV/Imm_Epi/inferResult/infer_", i, "/vlnPlot_CNV_level_N.pdf"),width = 20,height = 12,units = "cm")

##### integrate cancer group #####

cluster <- c()

for (i in list) {
  clustertmp <- read.table(paste0("/data_analysis_2/csl/cyf/SC/4.1.inferCNV/Imm_Epi/inferResult/infer_", i, "/kmeans_df_s.txt"), sep = '\t', row.names = 1, header = T)
  clustertmp$kmeans_class <- paste0(i, "_", clustertmp$kmeans_class)
  cluster <- rbind(cluster, clustertmp)
}

Epi <- subset(PCSC, celltype %in% c("Epithelial"))
cluster <- cluster[rownames(cluster) %in% colnames(Epi),]

cluster$CancerType <- "Normal"
cluster$CancerType[cluster$kmeans_class %in% c(
  "EOPC1_1", "EOPC1_2", "EOPC1_5",
  "EOPC2_1", "EOPC2_2", "EOPC2_3", "EOPC2_4", "EOPC2_8", 
  "EOPC3_1", "EOPC3_4", "EOPC3_7",
  "EOPC4_2",
  "LOPC1_1", "LOPC1_2", "LOPC1_4", "LOPC1_5", "LOPC1_7", "LOPC1_3",
  "LOPC2_5",
  "LOPC3_1", "LOPC3_4", "LOPC3_8", "LOPC3_6", 
  "LOPC4_2", "LOPC4_4", "LOPC4_5", "LOPC4_6", 
  "LOPC5_2", "LOPC5_4", "LOPC5_5", "LOPC5_6", "LOPC5_8",
  "LOPC6_1", "LOPC6_2", "LOPC6_4", "LOPC6_6", "LOPC6_8" 
)] <- "Cancer"
cluster$CancerType[cluster$kmeans_class %in% c(
  "EOPC1_8", 
  "EOPC2_6",
  "EOPC3_8", 
  "EOPC4_7",
  "LOPC1_8",
  "LOPC2_3",
  "LOPC4_1"
)] <- "Other_Epi"

PCSC@meta.data$inferCNVCellStat <- PCSC@meta.data$celltype
PCSC@meta.data$inferCNVCellStat <- as.character(PCSC@meta.data$inferCNVCellStat)
PCSC@meta.data$inferCNVCellStat[which(rownames(PCSC@meta.data) %in% rownames(cluster[cluster$CancerType == "Cancer",]))] <- "Cancer"
PCSC@meta.data$inferCNVCellStat[which(rownames(PCSC@meta.data) %in% rownames(cluster[cluster$CancerType == "Other_Epi",]))] <- "Other_Epi"
PCSC@meta.data$inferCNVCellStat[which(rownames(PCSC@meta.data) %in% rownames(cluster[cluster$CancerType == "Normal",]))] <- "Normal"

Idents(PCSC) <- "inferCNVCellStat"

pal_for_celltype <- color[c(1:13)]
names(pal_for_celltype) <- c("Epithelial","Normal","Other_Epi","T-cell","B-cell",
                             "Myeloid","Mast","Endothelial","Fibroblast", "Non1", "Non2", "SMC")

pdf("./4.1.inferCNV/InferCancer/UMAP_cellType_inferCancer_res0.2.pdf")
DimPlot(PCSC, label = T, pt.size=0.5, shuffle = T, cols = pal_for_celltype, label.size = 5) + 
  NoLegend()
dev.off()

saveRDS(PCSC, "./4.1.inferCNV/InferCancer/PCSC_inferred.Rds")
PCSC <- readRDS("./4.1.inferCNV/InferCancer/PCSC_inferred.Rds")

# All marker bubble ------------------------------------------

Allmarker <- c("EPCAM","CDH1", 
               "KLK3", "AR", "AMACR",
               "TP63", "KRT5", "FGFR2",
               "PTPRC", 
               "CD3E", "CD3D", "CD2", 
               "MS4A1", "CD79A", 
               "CD14", "FCGR3A", "LYZ", 
               "KIT", "MS4A2", "TPSAB1", 
               "VIM", 
               "PECAM1", "ENG", "VWF", 
               "FAP", "PDGFRA", "THY1", 
               "COL3A1", "ACTA2", "MYL9")

PCSC$inferCNVCellStat <- factor(PCSC$inferCNVCellStat, 
                                levels = c("Cancer", "Other_Epi", "Normal", "T-cell", "B-cell",
                                           "Myeloid", "Mast", "Endothelial", "Fibroblast", "SMC"))

pdf("./4.1.inferCNV/InferCancer/AllmarkerBubble.pdf",width = 15, height = 4.6)
jjDotPlot(PCSC,
          gene = Allmarker,
          id = "inferCNVCellStat",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0,0,0,0),
          dot.min = 0,
          dot.max = 6,
          x.text.vjust = 0.5)
dev.off()

# Epi clustering ------------------------------------------

Epi <- subset(PCSC, cellGroup %in% c("Epithelial_cells"))
Idents(Epi) <- "inferCNVCellStat"

Epi <- NormalizeData(Epi)

Epi <- FindVariableFeatures(Epi,
                            selection.method = "vst",
                            nfeatures = 2000,
                            verbose = FALSE) %>%
  ScaleData(vars.to.regress = c("percent.mt"))
Epi <- RunPCA(Epi)

# harmony 

Epi <- RunHarmony(Epi, "orig.ident", plot_convergence = T)

Epi <- RunUMAP(Epi, reduction = "harmony", dims = 1:40)

pal_for_celltype <- color[c(1:9)]
names(pal_for_celltype) <- c("Cancer","Other_Epi","Normal","T-cell","B-cell",
                             "Myeloid","Mast","Endothelial","Mesenchymal")

pdf("./4.1.inferCNV/InferCancer/DimPlot_Epi_celltype_umap.pdf", width = 5, height = 5)
DimPlot(Epi, group.by = "inferCNVCellStat" ,cols = pal_for_celltype, pt.size=0.5, shuffle = T) +
  ggtitle("")
dev.off()

pdf("./4.1.inferCNV/InferCancer/DimPlot_Epi_celltype_umap_NoLegend.pdf", width = 5, height = 5)
DimPlot(Epi, group.by = "inferCNVCellStat" ,cols = pal_for_celltype, pt.size=0.5, shuffle = T) +
  ggtitle("") + NoLegend()
dev.off()

# Feature VlnPlot --------------------------------------------

pdf("./4.1.inferCNV/InferCancer/inferCancer_feature_VlnPlot.pdf", height = 6, width = 3)
VlnPlot(object = Epi, features = c("KLK3", "AR", 
                                   "KRT5", "TP63", 
                                   "PTPRC", "CD3D"), 
        pt.size=0, cols = c("#A52A2ACC", "#053E7ACC", "#2F4F4FCC"), ncol = 2)
dev.off()

# Epi score (ssGSEA) ------------------------------------

CR_Epi <- readRDS("./4.1.inferCNV/InferCancer/CR_Epi.Rds")
expr <- as.matrix(Epi@assays$RNA@data)

Score <- gsva(expr, CR_Epi, verbose =TRUE, parallel.sz=10, method="ssgsea")

saveRDS(Score, "./4.1.inferCNV/InferCancer/ssGSEA_CRset_matrix.Rds")

#integrate

Epi_CR <- CreateSeuratObject(Score,
                             project = "Epi_CR", 
                             names.field = 1,
                             names.delim = "_")
Epi_CR@meta.data <- Epi@meta.data
Epi_CR@reductions <- Epi@reductions
Idents(Epi_CR) <- "inferCNVCellStat"

pdf("./4.1.inferCNV/InferCancer/inferCancer_ssGSEA_CRsets_VlnPlot.pdf", height = 6, width = 3)
VlnPlot(object = Epi_CR, features = c("Basal", "Luminal", "Neuroendocrine", "Club", "Hillock"), 
        pt.size=0, ncol = 2, cols = c("#A52A2ACC", "#053E7ACC", "#2F4F4FCC"))
dev.off()

saveRDS(Epi_hallmark,"./4.1.inferCNV/InferCancer/Epi_CRset_seuratObj.Rds")

# Epi score (AUCell) ------------------------------------

CR_Epi <- readRDS("./4.1.inferCNV/InferCancer/CR_Epi.Rds")
expr <- as.matrix(Epi@assays$RNA@data)

cells_rankings <- AUCell_buildRankings(expr, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(CR_Epi, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./4.1.inferCNV/InferCancer/Epi_CRset_AUC.Rds")

Score <- cells_AUC@assays@data$AUC %>% as.data.frame()

#integrate

Epi_CR <- CreateSeuratObject(Score,
                             project = "Epi_CR", 
                             names.field = 1,
                             names.delim = "_")
Epi_CR@meta.data <- Epi@meta.data
Epi_CR@reductions <- Epi@reductions
Idents(Epi_CR) <- "inferCNVCellStat"

pdf("./4.1.inferCNV/InferCancer/inferCancer_AUCell_CRsets_VlnPlot.pdf", height = 6, width = 3)
VlnPlot(object = Epi_CR, features = c("Basal", "Luminal", "Neuroendocrine", "Club", "Hillock"), 
        pt.size=0, ncol = 2, cols = c("#A52A2ACC", "#053E7ACC", "#E46B4FCC"))
dev.off()
