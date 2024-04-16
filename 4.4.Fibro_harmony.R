library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library("SingleR")
library(limma)
library(monocle)
library("clusterProfiler")
library(org.Hs.eg.db)
library("enrichplot")
library(harmony)
library(irGSEA)
library(reshape2)
library(ggpubr)

# subset analysis ---------------------------------------------------------------------

PCSC <- readRDS("./4.3.Doublet/PCSC_inferred_doublet.Rds")

# Mes --------------------------------------------------------------------

color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")
                    

folder <- "5.5.Fibr_harmony"
if(!dir.exists(folder)){
  dir.create(folder)
}

# 1. subset, harmony and clustering ----------------------------------------------------------------

Fibr <- subset(PCSC, celltype %in% c("Smooth_muscle", "Fibroblast"))

Fibr <- NormalizeData(Fibr)
Fibr <- FindVariableFeatures(Fibr,
                             selection.method = "vst",
                             nfeatures = 2000,
                             verbose = FALSE) %>%
  ScaleData(vars.to.regress = c("percent.mt"))
Fibr <- RunPCA(Fibr)

pdf("./5.5.Fibr_harmony/harmony_convergence.pdf")
Fibr <- RunHarmony(Fibr, "orig.ident", plot_convergence = T)
dev.off()

Fibr <- RunUMAP(Fibr, reduction = "harmony", dims = 1:20)

Fibr <- FindNeighbors(Fibr, reduction = "harmony", 
                      dims = 1:20)
for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1)) {
  Fibr <- FindClusters(Fibr, resolution = res)}

apply(Fibr@meta.data[, grep("RNA_snn_res", colnames(Fibr@meta.data))],2,table)

library(clustree)

p2_tree=clustree(Fibr@meta.data, prefix = "RNA_snn_res.") +
  scale_color_manual(values = color)
ggsave(plot=p2_tree, filename="./5.5.Fibr_harmony/Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(Fibr) <- "RNA_snn_res.0.8"

# 2.Cell cycle scoring ---------------------------------------------------------

folder <- "./5.5.Fibr_harmony/cellcycle"
if(!dir.exists(folder)){
  dir.create(folder)
}

str(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Fibr <- CellCycleScoring(Fibr,
                         g2m.features = g2m.genes,
                         s.features = s.genes)

pdf("./5.5.Fibr_harmony/cellcycle/cellcycle.pdf",width = 12, height = 4)
DimPlot(Fibr,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase") +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/cellcycle/cellcycle_SampleSplit.pdf",width = 12, height = 4)
DimPlot(Fibr, 
        reduction = "pca", 
        group.by = "orig.ident",
        split.by = "Phase",
        shuffle = T) +
  scale_color_manual(values = color) 
dev.off()

# 3.plotting --------------------------------------------------------

pdf("./5.5.Fibr_harmony/CellCluster-UMAPPlot_res0.8.pdf",width = 4,height = 4)
DimPlot(Fibr, reduction = "umap", pt.size = 0.3, label = T, shuffle = T) + NoLegend() +
  scale_color_manual(values = color)
dev.off()

# 4.AllMarkers and remove other celltypes -----------------------------------------

Features <- c("EPCAM","CDH1", "PTPRC", 
              "CD3E", "CD3G", "CD3D", "CD2", "MS4A1", "CD79A", "CD14", 
              "FCGR3A", "CD68", "CD163", "LYZ", 
              "KIT", "MS4A2", "TPSAB1", "TPSB2", "VIM", "PECAM1", "ENG", "VWF", 
              "FAP", "THY1", "COL3A1", "ACTA2")

library(scRNAtoolVis)

jjDotPlot(Fibr,
          gene = Features,
          id = "RNA_snn_res.0.8",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0.5,
          dot.max = 5)

## plotting -------------------------------------------------

pdf("./5.5.Fibr_harmony/CellCluster-UMAPPlot_res0.8.pdf",width = 4,height = 4)
DimPlot(Fibr, reduction = "umap", pt.size = 0.3, label = T, shuffle = T) + NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/CellCluster-UMAPPlot_SampleType.pdf",width = 4,height = 4)
DimPlot(Fibr, reduction = "umap", group.by = "SampleType", pt.size = 0.3, shuffle = T) +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/CellCluster-UMAPPlot_SampleType_noLegend.pdf",width = 4,height = 4)
DimPlot(Fibr, reduction = "umap", group.by = "SampleType", pt.size = 0.3, shuffle = T) +
  NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/CellCluster-UMAPPlot_Sample.pdf",width = 4,height = 4)
DimPlot(Fibr, reduction = "umap", group.by = "orig.ident", pt.size = 0.3, label = F, shuffle = T) +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/CellCluster-UMAPPlot_Sample_noLegend.pdf",width = 4,height = 4)
DimPlot(Fibr, reduction = "umap", group.by = "orig.ident", pt.size = 0.3, label = F, shuffle = T) + 
  NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/CellCluster-UMAPPlot_SampleTypeSplit.pdf",width = 8,height = 4)
DimPlot(Fibr, reduction = "umap", split.by = "SampleType", pt.size = 0.3, label = T) + 
  NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/cellcycle/CellCluster-UMAPPlot_cellcycle_res0.8.pdf",width = 4,height = 4)
DimPlot(Fibr, reduction = "umap", group.by = "Phase", pt.size = 0.3, label = T) +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/cellcycle/CellCluster-UMAPPlot_SampleTypeSplit_cellcycle_res0.8.pdf",width = 8,height = 4)
DimPlot(Fibr, reduction = "umap", group.by = "Phase", split.by = "SampleType", pt.size = 0.3, label = T) +
  scale_color_manual(values = color)
dev.off()

saveRDS(Fibr, "./5.5.Fibr_harmony/Fibr_harmonied.Rds")

# 5.FindMarkers ------------------------------------

Fibr_nomix <- subset(Fibr, RNA_snn_res.0.8 %in% c(0,1,2,3,4,6,7,10,12))
Fibr_nomix$RNA_snn_res.0.8 <- factor(Fibr_nomix$RNA_snn_res.0.8, 
                                     levels = c(0,1,2,3,4,6,7,10,12))
Idents(Fibr_nomix) <- "RNA_snn_res.0.8"

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)

Allmarkers <- FindAllMarkers(Fibr_nomix,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = T)
saveRDS(Allmarkers, "./5.5.Fibr_harmony/Markers/cluster_markers.Rds")

topmarkers <- Allmarkers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

#Heatmap by cluster

Cluster <- data.frame(ID = colnames(Fibr_nomix), Cluster = Idents(Fibr_nomix))
Cluster$Cluster <- factor(Cluster$Cluster, levels = c("12","4","7","3","2","1","0","10","6"))
Cluster <- Cluster[order(Cluster$Cluster),]
annotation <- data.frame(Cluster = paste0("c", Cluster$Cluster))
rownames(annotation) <- Cluster$ID

top <- data.frame(gene = topmarkers$gene, Cluster = topmarkers$cluster)
top$Cluster <- factor(top$Cluster, levels = c("12","4","7","3","2","1","0","10","6"))
top <- top[order(top$Cluster),]

data <- as.matrix(Fibr_nomix@assays$RNA@data)
genes <- top$gene
data <- data[rownames(data) %in% genes,]
data <- data[match(genes, rownames(data)),]
data <- data[, Cluster$ID]

Cluster_colors <- color[c(1:9)]
names(Cluster_colors) <- paste0("c", c(0,1,2,3,4,6,7,10,12))
ann_colors <- list(Cluster = Cluster_colors)

pdf("./5.5.Fibr_harmony/Markers/cluster_heatmap.pdf", width = 5, height = 8)
bk <- c(seq(-3,0,by=0.01),seq(0.01,3,by=0.01))
pheatmap(data,
         cluster_cols = F,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = F,
         show_rownames = T,
         scale="row",  #矫正
         breaks=bk,
         legend_breaks=seq(-3,3,1),
         fontsize = 8,
         fontsize_row = 5,
         fontsize_col = 12,
         main = "",
         angle_col = "90",
         #border_color = "grey",
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cutree_rows = 5)
dev.off()

# 6.CytoTrace ------------------------------------

exp <- as.matrix(Fibr@assays$RNA@counts)
celltype <- Fibr@meta.data$RNA_snn_res.0.8
names(celltype) <- rownames(Fibr@meta.data)

Sample <- Fibr@meta.data$orig.ident
names(Sample) <- rownames(Fibr@meta.data)
exp <- exp[-grep("^RP[SL]", rownames(exp)),]
Fibr_CytoTrace_result <- CytoTRACE(exp, enableFast = TRUE, batch = Sample, ncores = 20, subsamplesize = 1000)
saveRDS(Fibr_CytoTrace_result, "./5.5.Fibr_harmony/CytoTRACE/Fibro_CytoTrace_result.Rds")

emb <- Fibr@reductions$umap@cell.embeddings
plotCytoTRACE(Fibr_CytoTrace_result, emb = emb, outputDir = "./5.5.Fibr_harmony/CytoTRACE/")
plotCytoGenes(Fibr_CytoTrace_result, numOfGenes = 20, outputDir = "./5.5.Fibr_harmony/CytoTRACE/")

# 7.Monocle ------------------------------------

#构造表达及注释数据
exp.matrix <- as(as.matrix(Fibr@assays$RNA@data), "sparseMatrix")
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(exp.matrix)
exp_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- Fibr@meta.data
rownames(sample_ann) <- colnames(exp.matrix)
exp_pd <- new("AnnotatedDataFrame", data =sample_ann)

#生成monocle对象
exp.monocle<-newCellDataSet(exp.matrix, phenoData = exp_pd,featureData = exp_fd, expressionFamily = negbinomial.size())
head(pData(exp.monocle))
head(fData(exp.monocle))

#计算sizefactor
exp.monocle <- estimateSizeFactors(exp.monocle)
exp.monocle <- estimateDispersions(exp.monocle)

#根据seurat cluster计算差异表达基因并挑选用于构建拟时序轨迹的基因
diff_test_res<-differentialGeneTest(exp.monocle,fullModelFormulaStr = "~RNA_snn_res.0.8") 
ordering_genes<-row.names (subset(diff_test_res, qval < 0.01))
saveRDS(ordering_genes, "./5.5.Fibr_harmony/Monocle/ordering_genes.Rds")
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
#plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

#修改monocle对象中的列名示例
saveRDS(exp.monocle, "./5.5.Fibr_harmony/Monocle/exp.monocle.Rds")

#将不同分组情况的拟时序轨迹图画到一起

plot1<-plot_cell_trajectory(exp.monocle, color_by = "RNA_snn_res.0.8",cell_size = 0.2, shuffle = T)
plot2<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size = 0.2, shuffle = T)
plot3<-plot_cell_trajectory(exp.monocle, color_by = "SampleType",cell_size = 0.2, shuffle = T)
pdf("./5.5.Fibr_harmony/Monocle/trajectory_plot.pdf",width = 3,height = 9)
CombinePlots(plots = list(plot1,plot2,plot3),legend = NULL,ncol=1)
dev.off()
rm(plot1,plot2,plot3)

#split by celltype

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

celltype_name <- as.character(levels(exp.monocle$RNA_snn_res.0.8))
Plist <- list()

#cluster
for (i in 1:length(celltype_name)) {
  a <- rownames(Fibr@meta.data)[ Fibr@meta.data$RNA_snn_res.0.8 == celltype_name[i] ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$RNA_snn_res.0.8 == celltype_name[i] ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(celltype_name[i])
  Plist[[i]] <- plot
  names(Plist) <- celltype_name[i]
}

Plist[[length(Plist)+1]] <- plot2
pdf("./5.5.Fibr_harmony/Monocle/trajectory_plot_celltypeSplit.pdf",width = 9,height = 6)
CombinePlots(plots = Plist,legend = NULL,ncol=4)
dev.off()

# 8.Annotation ------------------------------------

Features <- c("VIM", "FGF1", "FAP", "S100A4", "TNC", "POSTN", "DES", "PDGFRA", "PDGFRB", 
              "THY1", "PDPN", "ITGB1", "CAV1", "ACTA2", "SCRG1", "CNN1", "EPCAM", "CD3E",
              "LYZ", "PECAM1")

pdf("./5.5.Fibr_harmony/Scatter/Marker_scatter.pdf", width = 13, height = 9)
FeaturePlot(Fibr, features = Features, cols = c("grey", "#A52A2A"), ncol = 5)
dev.off()


## annotation --------------------------------------

celltype=data.frame(ClusterID=0:12,
                    celltype='SMC-1')   
celltype[celltype$ClusterID %in% c(3),2]='SMC-2'  
celltype[celltype$ClusterID %in% c(4),2]='myCAF-1' 
celltype[celltype$ClusterID %in% c(7),2]='myCAF-2' 
celltype[celltype$ClusterID %in% c(0,10),2]='iCAF-1' 
celltype[celltype$ClusterID %in% c(6),2]='iCAF-2'
celltype[celltype$ClusterID %in% c(12),2]='Myoblast'
celltype[celltype$ClusterID %in% c(5),2]='Epi_mix'
celltype[celltype$ClusterID %in% c(9),2]='Tcell_mix'
celltype[celltype$ClusterID %in% c(11),2]='Mye_mix'
celltype[celltype$ClusterID %in% c(8),2]='Endo_mix'
head(celltype)

Fibr@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  Fibr@meta.data[which(Fibr@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(Fibr@meta.data$celltype)

Fibr$celltype <- factor(Fibr$celltype, levels = c("SMC-1", "SMC-2", "myCAF-1", "myCAF-2", "iCAF-1", "iCAF-2", 
                                                  "Myoblast", "Epi_mix", "Endo_mix", "Tcell_mix", "Mye_mix"))
Idents(Fibr) <- "celltype"

## Findmarkers ----------------------------------------------

Fibr_nomix <- subset(Fibr, RNA_snn_res.0.8 %in% c(0,1,2,3,4,6,7,10,12))
Fibr_nomix$celltype <- factor(Fibr_nomix$celltype, 
                              levels = c("SMC-1", "SMC-2", "myCAF-1", "myCAF-2", 
                                         "iCAF-1", "iCAF-2", "Myoblast"))
Idents(Fibr_nomix) <- "celltype"

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)

Allmarkers <- FindAllMarkers(Fibr_nomix,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = T)
saveRDS(Allmarkers, "./5.5.Fibr_harmony/Markers/celltype_markers.Rds")

topmarkers <- Allmarkers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

#Heatmap by cluster

Cluster <- data.frame(ID = colnames(Fibr_nomix), Cluster = Idents(Fibr_nomix))
Cluster <- Cluster[order(Cluster$Cluster),]
annotation <- data.frame(Cluster = Cluster$Cluster)
rownames(annotation) <- Cluster$ID

top <- data.frame(gene = topmarkers$gene, Cluster = topmarkers$cluster)
top <- top[order(top$Cluster),]

data <- as.matrix(Fibr_nomix@assays$RNA@data)
genes <- top$gene
data <- data[rownames(data) %in% genes,]
data <- data[match(genes, rownames(data)),]
data <- data[, Cluster$ID]

Cluster_colors <- color[c(1:7)]
names(Cluster_colors) <- as.character(levels(Idents(Fibr_nomix)))
ann_colors <- list(Cluster = Cluster_colors)

pdf("./5.5.Fibr_harmony/Markers/celltype_heatmap.pdf", width = 5, height = 8)
bk <- c(seq(-3,0,by=0.01),seq(0.01,3,by=0.01))
pheatmap(data,
         cluster_cols = F,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = F,
         show_rownames = T,
         scale="row",  #矫正
         breaks=bk,
         legend_breaks=seq(-3,3,1),
         fontsize = 8,
         fontsize_row = 8,
         fontsize_col = 12,
         main = "",
         angle_col = "90",
         #border_color = "grey",
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cutree_rows = 5)
dev.off()

# plotting

pdf("./5.5.Fibr_harmony/CellType-UMAPPlot_res0.8.pdf", width = 5, height = 5)
DimPlot(Fibr, label = TRUE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color) + 
  NoLegend() + 
  theme(plot.title = element_blank())
dev.off()

Features <- c("ACTA2", "CNN1", "TNC", "DES", "ATF3",
              "COL1A1", "COL1A2", "THY1", "RGS5", "CCL21",
              "FAP", "PDGFRA", "PDPN", "IGF1", "CCL2",
              "MYF5", "EPCAM", "PECAM1", "CD3E", "CD68")

pdf("./5.5.Fibr_harmony/AllmarkerBubble.pdf",width = 5, height = 4.8)
jjDotPlot(Fibr,
          gene = Features,
          id = "celltype",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 5,
          x.text.vjust = 0.5)
dev.off()

# vioheatmap
Cluster <- data.frame(Cluster = Idents(Fibr))
data <- as.matrix(Fibr@assays$RNA@data)
data <- data[rownames(data) %in% Features,]
data <- as.data.frame(t(data))
data <- merge(Cluster, data, by='row.names')
data <- column_to_rownames(data, "Row.names")
library(reshape2)
data <- melt(data)
colnames(data) <- c("Cluster", "gene", "Exp")

# add median expression to group per gene
map_df(unique(data$Cluster),function(x){
  tmp <- data %>% filter(Cluster == x)
  map_df(unique(tmp$gene),function(j){
    tmp1 <- tmp %>% filter(gene == j)
    # calculate median expressions
    tmp1$Avg_Exp <- median(tmp1$Exp)
    return(tmp1)
  }) -> res
  return(res)
}) -> mean

mean$gene <- factor(mean$gene, levels = Features)

library(ggthemes)
pdf("./5.5.Fibr_harmony/VlnHeatmapPlot_markers.pdf")
ggplot(mean, aes(x = gene, y = Cluster)) +
  geom_jjviomap(aes(val = Exp, fill = Cluster), width = 1) +
  coord_fixed() +
  theme_base() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color) +
  NoLegend()
dev.off()

saveRDS(Fibr, "./5.5.Fibr_harmony/Fibr_annoted.Rds")

## Roe -------------------------------------------

tab <- table(Fibr$celltype, Fibr$SampleType)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.5.Fibr_harmony/celltype_Roe_SampleType.pdf", width = 2, height = 3)
pheatmap(Roe,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(colors = c("white","#A52A2A"))(100),
         show_colnames = T,
         show_rownames = T,
         scale="none",  #矫正
         #breaks=bk,
         #legend_breaks=seq(-0.5,0.5,0.5),
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Ro/e",
         angle_col = "0",
         display_numbers = T,
         legend = F)
dev.off()

tab <- table(Fibr_nomix$celltype, Fibr_nomix$SampleType)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.5.Fibr_harmony/celltype_Roe_SampleType_nomix.pdf", width = 2, height = 3)
pheatmap(Roe,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(colors = c("white","#A52A2A"))(100),
         show_colnames = T,
         show_rownames = T,
         scale="none",  #矫正
         #breaks=bk,
         #legend_breaks=seq(-0.5,0.5,0.5),
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Ro/e",
         angle_col = "0",
         display_numbers = T,
         legend = F)
dev.off()

pdf("./5.5.Fibr_harmony/celltype_Roe_SampleType_nomix_invert.pdf", width = 5, height = 1.5)
pheatmap(Roe,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(colors = c("white","#A52A2A"))(100),
         show_colnames = T,
         show_rownames = T,
         scale="none",  #矫正
         #breaks=bk,
         #legend_breaks=seq(-0.5,0.5,0.5),
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Ro/e",
         angle_col = "45",
         display_numbers = T,
         legend = F)
dev.off()

# 9.nomix -------------------------------------------

# 9.1.clustering ---------------------------------------

Fibr <- readRDS("./5.5.Fibr_harmony/Fibr_annoted.Rds")
Fibr_nomix <- subset(Fibr, celltype != "Epi_mix" & celltype != "Endo_mix" & celltype != "Tcell_mix" & celltype != "Mye_mix")

Fibr_nomix <- NormalizeData(Fibr_nomix)
Fibr_nomix <- FindVariableFeatures(Fibr_nomix,
                             selection.method = "vst",
                             nfeatures = 2000,
                             verbose = FALSE) %>%
  ScaleData(vars.to.regress = c("percent.mt"))
Fibr_nomix <- RunPCA(Fibr_nomix)

pdf("./5.5.Fibr_harmony/Nomix_harmony_convergence.pdf")
Fibr_nomix <- RunHarmony(Fibr_nomix, "orig.ident", plot_convergence = T)
dev.off()

Fibr_nomix <- RunUMAP(Fibr_nomix, reduction = "harmony", dims = 1:20)

Fibr_nomix <- FindNeighbors(Fibr_nomix, reduction = "harmony", 
                      dims = 1:20)
for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1)) {
  Fibr_nomix <- FindClusters(Fibr_nomix, resolution = res)}

apply(Fibr_nomix@meta.data[, grep("RNA_snn_res", colnames(Fibr_nomix@meta.data))],2,table)

library(clustree)

p2_tree=clustree(Fibr_nomix@meta.data, prefix = "RNA_snn_res.") +
  scale_color_manual(values = color)
ggsave(plot=p2_tree, filename="./5.5.Fibr_harmony/Nomix_Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(Fibr_nomix) <- "RNA_snn_res.0.8"

pdf("./5.5.Fibr_harmony/Nomix_CellCluster-UMAPPlot_res0.8.pdf",width = 4,height = 4)
DimPlot(Fibr_nomix, label = TRUE, pt.size=0.3, shuffle = T, label.box = T,
        label.size = 4, label.color = "black", cols = color, repel = FALSE) + NoLegend() +
  theme(plot.title = element_blank())
dev.off()

# 9.2.Other Markers -----------------------------------------

Features <- c("EPCAM","CDH1", "PTPRC", 
              "CD3E", "CD3G", "CD3D", "CD2", "MS4A1", "CD79A", "CD14", 
              "FCGR3A", "CD68", "CD163", "LYZ", 
              "KIT", "MS4A2", "TPSAB1", "TPSB2", "VIM", "PECAM1", "ENG", "VWF", 
              "FAP", "THY1", "COL3A1", "ACTA2")

library(scRNAtoolVis)

jjDotPlot(Fibr_nomix,
          gene = Features,
          id = "RNA_snn_res.0.8",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0.5,
          dot.max = 5)

## plotting 

pdf("./5.5.Fibr_harmony/Nomix_CellCluster-UMAPPlot_res0.8.pdf",width = 4,height = 4)
DimPlot(Fibr_nomix, reduction = "umap", pt.size = 0.3, label = T, shuffle = T) + NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/Nomix_CellCluster-UMAPPlot_SampleType.pdf",width = 4,height = 4)
DimPlot(Fibr_nomix, reduction = "umap", group.by = "SampleType", pt.size = 0.3, shuffle = T) +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/Nomix_CellCluster-UMAPPlot_SampleType_noLegend.pdf",width = 4,height = 4)
DimPlot(Fibr_nomix, reduction = "umap", group.by = "SampleType", pt.size = 0.3, shuffle = T) +
  NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/Nomix_CellCluster-UMAPPlot_Sample.pdf",width = 4,height = 4)
DimPlot(Fibr_nomix, reduction = "umap", group.by = "orig.ident", pt.size = 0.3, label = F, shuffle = T) +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/Nomix_CellCluster-UMAPPlot_Sample_noLegend.pdf",width = 4,height = 4)
DimPlot(Fibr_nomix, reduction = "umap", group.by = "orig.ident", pt.size = 0.3, label = F, shuffle = T) + 
  NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/Nomix_CellCluster-UMAPPlot_SampleTypeSplit.pdf",width = 8,height = 4)
DimPlot(Fibr_nomix, reduction = "umap", split.by = "SampleType", pt.size = 0.3, label = T) + 
  NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/cellcycle/Nomix_CellCluster-UMAPPlot_cellcycle_res0.8.pdf",width = 4,height = 4)
DimPlot(Fibr_nomix, reduction = "umap", group.by = "Phase", pt.size = 0.3, label = T) +
  scale_color_manual(values = color)
dev.off()

pdf("./5.5.Fibr_harmony/cellcycle/Nomix_CellCluster-UMAPPlot_SampleTypeSplit_cellcycle_res0.8.pdf",width = 8,height = 4)
DimPlot(Fibr_nomix, reduction = "umap", group.by = "Phase", split.by = "SampleType", pt.size = 0.3, label = T) +
  scale_color_manual(values = color)
dev.off()

saveRDS(Fibr_nomix, "./5.5.Fibr_harmony/Nomix_Fibr_harmonied.Rds")

# 9.3.FindMarkers ------------------------------------

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)

Allmarkers <- FindAllMarkers(Fibr_nomix,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = T)
saveRDS(Allmarkers, "./5.5.Fibr_harmony/Markers/Nomix_cluster_markers.Rds")

topmarkers <- Allmarkers %>% group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

#Heatmap by cluster

Cluster <- data.frame(ID = colnames(Fibr_nomix), Cluster = Idents(Fibr_nomix))
Cluster <- Cluster[order(Cluster$Cluster),]
annotation <- data.frame(Cluster = paste0("c", Cluster$Cluster))
rownames(annotation) <- Cluster$ID

top <- data.frame(gene = topmarkers$gene, Cluster = topmarkers$cluster)
top <- top[order(top$Cluster),]

data <- as.matrix(Fibr_nomix@assays$RNA@data)
genes <- top$gene
data <- data[rownames(data) %in% genes,]
data <- data[match(genes, rownames(data)),]
data <- data[, Cluster$ID]

Cluster_colors <- color[c(1:11)]
names(Cluster_colors) <- paste0("c", c(0:10))
ann_colors <- list(Cluster = Cluster_colors)

pdf("./5.5.Fibr_harmony/Markers/Nomix_cluster_heatmap.pdf", width = 5, height = 8)
bk <- c(seq(-2,0,by=0.01),seq(0.01,2,by=0.01))
pheatmap(data,
         cluster_cols = F,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = F,
         show_rownames = T,
         scale="row",  #矫正
         breaks=bk,
         legend_breaks=seq(-2,2,1),
         fontsize = 8,
         fontsize_row = 5,
         fontsize_col = 12,
         main = "",
         angle_col = "90",
         #border_color = "grey",
         annotation_col = annotation,
         annotation_colors = ann_colors)
dev.off()

pdf("./5.5.Fibr_harmony/Markers/Nomix_cluster_heatmap_top5.pdf", width = 5, height = 6)
bk <- c(seq(-2,0,by=0.01),seq(0.01,2,by=0.01))
pheatmap(data,
         cluster_cols = F,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = F,
         show_rownames = T,
         scale="row",  #矫正
         breaks=bk,
         legend_breaks=seq(-2,2,1),
         fontsize = 8,
         fontsize_row = 8,
         fontsize_col = 12,
         main = "",
         angle_col = "90",
         #border_color = "grey",
         annotation_col = annotation,
         annotation_colors = ann_colors)
dev.off()

# 9.4.Monocle ------------------------------------

#构造表达及注释数据
exp.matrix <- as(as.matrix(Fibr@assays$RNA@data), "sparseMatrix")
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(exp.matrix)
exp_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- Fibr@meta.data
rownames(sample_ann) <- colnames(exp.matrix)
exp_pd <- new("AnnotatedDataFrame", data =sample_ann)

#生成monocle对象
exp.monocle<-newCellDataSet(exp.matrix, phenoData = exp_pd,featureData = exp_fd, expressionFamily = negbinomial.size())
head(pData(exp.monocle))
head(fData(exp.monocle))

#计算sizefactor
exp.monocle <- estimateSizeFactors(exp.monocle)
exp.monocle <- estimateDispersions(exp.monocle)

#根据seurat cluster计算差异表达基因并挑选用于构建拟时序轨迹的基因
diff_test_res<-differentialGeneTest(exp.monocle,fullModelFormulaStr = "~RNA_snn_res.0.8") 
ordering_genes<-row.names (subset(diff_test_res, qval < 0.01))
saveRDS(ordering_genes, "./5.5.Fibr_harmony/Monocle/ordering_genes.Rds")
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
#plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

#修改monocle对象中的列名示例
saveRDS(exp.monocle, "./5.5.Fibr_harmony/Monocle/exp.monocle.Rds")

#将不同分组情况的拟时序轨迹图画到一起

plot1<-plot_cell_trajectory(exp.monocle, color_by = "RNA_snn_res.0.8",cell_size = 0.2, shuffle = T)
plot2<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size = 0.2, shuffle = T)
plot3<-plot_cell_trajectory(exp.monocle, color_by = "SampleType",cell_size = 0.2, shuffle = T)
pdf("./5.5.Fibr_harmony/Monocle/trajectory_plot.pdf",width = 3,height = 9)
CombinePlots(plots = list(plot1,plot2,plot3),legend = NULL,ncol=1)
dev.off()
rm(plot1,plot2,plot3)

#split by celltype

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

celltype_name <- as.character(levels(exp.monocle$RNA_snn_res.0.8))
Plist <- list()

#cluster
for (i in 1:length(celltype_name)) {
  a <- rownames(Fibr@meta.data)[ Fibr@meta.data$RNA_snn_res.0.8 == celltype_name[i] ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$RNA_snn_res.0.8 == celltype_name[i] ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(celltype_name[i])
  Plist[[i]] <- plot
  names(Plist) <- celltype_name[i]
}

Plist[[length(Plist)+1]] <- plot2
pdf("./5.5.Fibr_harmony/Monocle/trajectory_plot_celltypeSplit.pdf",width = 9,height = 6)
CombinePlots(plots = Plist,legend = NULL,ncol=4)
dev.off()

# 9.5.Annotation -----------------------------------------

Features <- c("VIM", "FGF1", "FAP", "S100A4", "TNC", "POSTN", "DES", "PDGFRA", "PDGFRB", 
              "THY1", "PDPN", "ITGB1", "CAV1", "ACTA2", "CNN1")

pdf("./5.5.Fibr_harmony/Scatter/Nomix_Marker_scatter.pdf", width = 13, height = 9)
FeaturePlot(Fibr_nomix, features = Features, cols = c("grey", "#A52A2A"), ncol = 5)
dev.off()

Plist <- list()
for (i in 1:length(Features)) {
  P <- FeaturePlot(Fibr_nomix, features = Features[i], cols = c("grey", "#A52A2A")) + NoLegend()
  Plist[[i]] <- P
}

png("./5.5.Fibr_harmony/Scatter/Nomix_Marker_scatter.png", res = 300, units = 'in', width = 6, height = 9)
grid.arrange(grobs = Plist, ncol = 3)
dev.off()

## annotation --------------------------------------

celltype=data.frame(ClusterID=0:12,
                    celltype='SMC_KLF2')   
celltype[celltype$ClusterID %in% c(2),2]='SMC_HOPX'  
celltype[celltype$ClusterID %in% c(3,8),2]='myCAF_RGS5' 
celltype[celltype$ClusterID %in% c(6),2]='myCAF_CCL21' 
celltype[celltype$ClusterID %in% c(0,9),2]='iCAF_APOD' 
celltype[celltype$ClusterID %in% c(5),2]='iCAF_CCL2'
celltype[celltype$ClusterID %in% c(10),2]='Myoblast_MYF5'
head(celltype)

Fibr_nomix@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  Fibr_nomix@meta.data[which(Fibr_nomix@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(Fibr_nomix@meta.data$celltype)

Fibr_nomix$celltype <- factor(Fibr_nomix$celltype, levels = c("SMC_KLF2", "SMC_HOPX", "myCAF_RGS5", "myCAF_CCL21", 
                                                              "iCAF_APOD", "iCAF_CCL2", "Myoblast_MYF5"))
Idents(Fibr_nomix) <- "celltype"

## Findmarkers ----------------------------------------------

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)

Allmarkers <- FindAllMarkers(Fibr_nomix,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = T)
saveRDS(Allmarkers, "./5.5.Fibr_harmony/Markers/Nomix_celltype_markers.Rds")

sigmarkers <- Allmarkers[which(Allmarkers$avg_log2FC > 0.25 & 
                                 Allmarkers$p_val_adj < 0.01 &
                                 Allmarkers$pct.1 > 0.3),]

sigmarkers_order <- c()
for (i in as.character(unique(Allmarkers$cluster))) {
  
  tmp <- sigmarkers[which(sigmarkers$cluster == i),]
  tmp <- tmp[order(tmp$avg_log2FC, decreasing = T),]
  sigmarkers_order <- rbind(sigmarkers_order, tmp)
  
}

write.table(sigmarkers_order, "./5.5.Fibr_harmony/Markers/Nomix_celltype_sigmarkers.txt",
            quote = F, sep = "\t", row.names = F)

topmarkers <- Allmarkers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

#Heatmap by cluster

Cluster <- data.frame(ID = colnames(Fibr_nomix), Cluster = Idents(Fibr_nomix))
Cluster <- Cluster[order(Cluster$Cluster),]
annotation <- data.frame(Cluster = Cluster$Cluster)
rownames(annotation) <- Cluster$ID

top <- data.frame(gene = topmarkers$gene, Cluster = topmarkers$cluster)
top <- top[order(top$Cluster),]

data <- as.matrix(Fibr_nomix@assays$RNA@data)
genes <- top$gene
data <- data[rownames(data) %in% genes,]
data <- data[match(genes, rownames(data)),]
data <- data[, Cluster$ID]

Cluster_colors <- color[c(1:7)]
names(Cluster_colors) <- as.character(levels(Idents(Fibr_nomix)))
ann_colors <- list(Cluster = Cluster_colors)

pdf("./5.5.Fibr_harmony/Markers/Nomix_celltype_heatmap.pdf", width = 5, height = 8)
bk <- c(seq(-2,0,by=0.01),seq(0.01,2,by=0.01))
pheatmap(data,
         cluster_cols = F,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = F,
         show_rownames = T,
         scale="row",  #矫正
         breaks=bk,
         legend_breaks=seq(-2,2,1),
         fontsize = 8,
         fontsize_row = 8,
         fontsize_col = 12,
         main = "",
         angle_col = "90",
         #border_color = "grey",
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cutree_rows = 5)
dev.off()

# plotting

pdf("./5.5.Fibr_harmony/Nomix_CellType-UMAPPlot_res0.8.pdf", width = 5, height = 5)
DimPlot(Fibr_nomix, label = TRUE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color, repel = T) + 
  NoLegend() + 
  theme(plot.title = element_blank())
dev.off()

pdf("./5.5.Fibr_harmony/Nomix_CellType-UMAPPlot_res0.8_noLebel.pdf", width = 3, height = 3)
DimPlot(Fibr_nomix, label = FALSE, pt.size=0.3, group.by = 'celltype', 
        shuffle = T , cols = color) + 
  NoLegend() + 
  theme(plot.title = element_blank())
dev.off()

png("./5.5.Fibr_harmony/Nomix_CellType-UMAPPlot_res0.8.png", width = 360, height = 360)
DimPlot(Fibr_nomix, label = TRUE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color, repel = T) + 
  NoLegend() + 
  theme(plot.title = element_blank())
dev.off()

Features <- c("CNN1", "MYL9", "KLF2", "HOPX",
              "COL1A1", "COL1A2", "DCN",
              "THY1", "ACTA2", "TAGLN", "RGS5", "CCL21",
              "FAP", "PDGFRA", "PDPN", "APOD", "CCL2",
              "MYF5")

pdf("./5.5.Fibr_harmony/Nomix_Rotated_AllmarkerBubble.pdf",width = 2, height = 4)
jjDotPlot(Fibr_nomix,
          gene = Features,
          id = "celltype",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 3
          ) +
  coord_flip() + theme(legend.text = element_text(size = 8),
                       legend.title = element_text(size = 10),
                       axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 9),
                       axis.title.x = element_blank(),
                       axis.text.y = element_text(size = 8),
                       axis.title.y = element_blank())
dev.off()

# vioheatmap
Cluster <- data.frame(Cluster = Idents(Fibr_nomix))
data <- as.matrix(Fibr_nomix@assays$RNA@data)
data <- data[rownames(data) %in% Features,]
data <- as.data.frame(t(data))
data <- merge(Cluster, data, by='row.names')
data <- column_to_rownames(data, "Row.names")
library(reshape2)
data <- melt(data)
colnames(data) <- c("Cluster", "gene", "Exp")

# add median expression to group per gene
map_df(unique(data$Cluster),function(x){
  tmp <- data %>% filter(Cluster == x)
  map_df(unique(tmp$gene),function(j){
    tmp1 <- tmp %>% filter(gene == j)
    # calculate median expressions
    tmp1$Avg_Exp <- median(tmp1$Exp)
    return(tmp1)
  }) -> res
  return(res)
}) -> mean

mean$gene <- factor(mean$gene, levels = Features)

library(ggthemes)
pdf("./5.5.Fibr_harmony/VlnHeatmapPlot_markers.pdf")
ggplot(mean, aes(x = gene, y = Cluster)) +
  geom_jjviomap(aes(val = Exp, fill = Cluster), width = 1) +
  coord_fixed() +
  theme_base() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color) +
  NoLegend()
dev.off()

saveRDS(Fibr_nomix, "./5.5.Fibr_harmony/Nomix_Fibr_annoted.Rds")

## Roe -------------------------------------------

tab <- table(Fibr_nomix$celltype, Fibr_nomix$SampleType)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.5.Fibr_harmony/Nomix_celltype_Roe_SampleType.pdf", width = 5, height = 2)
pheatmap(t(Roe),
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(colors = c("white","#A52A2A"))(100),
         show_colnames = T,
         show_rownames = T,
         scale="none",  #矫正
         #breaks=bk,
         #legend_breaks=seq(-0.5,0.5,0.5),
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Ro/e",
         angle_col = "45",
         display_numbers = T,
         legend = F)
dev.off()

pdf("./5.5.Fibr_harmony/celltype_Roe_SampleType_nomix_invert.pdf", width = 2, height = 3)
pheatmap(Roe,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(colors = c("white","#A52A2A"))(100),
         show_colnames = T,
         show_rownames = T,
         scale="none",  #矫正
         #breaks=bk,
         #legend_breaks=seq(-0.5,0.5,0.5),
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Ro/e",
         angle_col = "45",
         display_numbers = T,
         legend = F)
dev.off()

# 10.celltype ecrichment -------------------------------------

Allmarkers <- readRDS("./5.5.Fibr_harmony/Markers/Nomix_celltype_markers.Rds")

# 10.1.enrichment KEGG (GSEA) -----------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
#remotes::install_github("YuLab-SMU/createKEGGdb")
#createKEGGdb::create_kegg_db('hsa')
#install.packages("KEGG.db_1.0.tar.gz",repos=NULL,type="source")
library(KEGG.db)

Allmarkers$cluster <- as.character(Allmarkers$cluster)

KEGG_gse_all <- list()
for (i in 1:length(unique(Allmarkers$cluster))) {
  gene <- Allmarkers[Allmarkers$cluster == unique(Allmarkers$cluster)[i],]
  gene <- data.frame(Gene = gene$gene, logFC = gene$avg_log2FC)
  genename <- as.character(gene[,1])
  gene_map <- AnnotationDbi::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
  colnames(gene_map)[1]<-"Gene"
  
  gene <- inner_join(gene_map, gene, by = "Gene")
  gene <- gene[,-1]
  gene <- na.omit(gene)
  gene <- gene[order(gene$logFC, decreasing = T),]
  
  geneList = gene[,2]
  names(geneList) = as.character(gene[,1])
  
  KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1, eps = 1e-100)
  saveRDS(KEGG_gseresult, paste0("5.5.Fibr_harmony/Markers/Nomix_", unique(Allmarkers$cluster)[i], "_KEGG_gseresult.Rds"))
  
  KEGG_gse_all[[i]] <- KEGG_gseresult
  names(KEGG_gse_all)[i] <- unique(Allmarkers$cluster)[i]
}

saveRDS(KEGG_gse_all, "5.5.Fibr_harmony/Markers/Nomix_all_KEGG_gseresult.Rds")

library(enrichplot)
pdf("5.5.Fibr_harmony/Markers/Nomix_KEGG_gseplot_SMC1.pdf",width = 12)
gseaplot2(KEGG_gse_all[[1]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=,title="SMC-1")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_KEGG_gseplot_SMC2.pdf",width = 12)
gseaplot2(KEGG_gse_all[[2]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="SMC-2")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_KEGG_gseplot_MyoFib1.pdf",width = 12)
gseaplot2(KEGG_gse_all[[1]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="SMC-1")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_KEGG_gseplot_MyoFib2.pdf",width = 12)
gseaplot2(KEGG_gse_all[[2]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="SMC-2")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_KEGG_gseplot_CAF1.pdf",width = 12)
gseaplot2(KEGG_gse_all[[5]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="CAF-1")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_KEGG_gseplot_CAF2.pdf",width = 12)
gseaplot2(KEGG_gse_all[[6]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="CAF2")
dev.off()

# 10.1.2.enrichment GOBP (GSEA) -----------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)

Allmarkers$cluster <- as.character(Allmarkers$cluster)

GOBP_gse_all <- list()
for (i in 1:length(unique(Allmarkers$cluster))) {
  gene <- Allmarkers[Allmarkers$cluster == unique(Allmarkers$cluster)[i],]
  gene <- data.frame(Gene = gene$gene, logFC = gene$avg_log2FC)
  genename <- as.character(gene[,1])
  gene_map <- AnnotationDbi::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
  colnames(gene_map)[1]<-"Gene"
  
  gene <- inner_join(gene_map, gene, by = "Gene")
  gene <- gene[,-1]
  gene <- na.omit(gene)
  gene <- gene[order(gene$logFC, decreasing = T),]
  
  geneList = gene[,2]
  names(geneList) = as.character(gene[,1])
  
  GOBP_gseresult <- gseGO(geneList, OrgDb = "org.Hs.eg.db", nPerm = 1000, minGSSize = 10, 
                          maxGSSize = 1000, pvalueCutoff=1, eps = 1e-100, ont = "BP")
  saveRDS(GOBP_gseresult, paste0("5.5.Fibr_harmony/Markers/Nomix_", unique(Allmarkers$cluster)[i], "_GOBP_gseresult.Rds"))
  
  GOBP_gse_all[[i]] <- GOBP_gseresult
  names(GOBP_gse_all)[i] <- unique(Allmarkers$cluster)[i]
}

saveRDS(GOBP_gse_all, "5.5.Fibr_harmony/Markers/Nomix_all_GOBP_gseresult.Rds")

library(enrichplot)
pdf("5.5.Fibr_harmony/Markers/Nomix_GOBP_gseplot_SMC1.pdf",width = 5,height = 3)
gseaplot2(GOBP_gse_all[[1]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=,title="SMC-1")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_GOBP_gseplot_SMC2.pdf",width = 5,height = 3)
gseaplot2(GOBP_gse_all[[2]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="SMC-2")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_GOBP_gseplot_MyoFib1.pdf",width = 5,height = 3)
gseaplot2(GOBP_gse_all[[3]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="MyoFibro-1")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_GOBP_gseplot_MyoFib2.pdf",width = 5,height = 3)
gseaplot2(GOBP_gse_all[[4]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="MyoFibro-2")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_GOBP_gseplot_CAF1.pdf",width = 5,height = 3)
gseaplot2(GOBP_gse_all[[5]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="CAF-1")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_GOBP_gseplot_CAF2.pdf",width = 5,height = 3)
gseaplot2(GOBP_gse_all[[6]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="CAF2")
dev.off()
pdf("5.5.Fibr_harmony/Markers/Nomix_GOBP_gseplot_Myoblast.pdf",width = 5,height = 3)
gseaplot2(GOBP_gse_all[[7]],1:5,pvalue_table = T,rel_heights=c(1.5,0,0),color=color,title="Myoblast")
dev.off()


# 11.2.enrichment KEGG (clusterProfiler) -----------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)

Allmarkers$cluster <- as.character(Allmarkers$cluster)
Upmarkers <- Allmarkers[which(Allmarkers$avg_log2FC > 0.25 & Allmarkers$p_val_adj < 0.05),]
topmarkers <- Upmarkers %>% group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

KEGG_enrich_all <- list()
for (i in 1:length(unique(topmarkers$cluster))) {
  gene <- topmarkers[topmarkers$cluster == unique(topmarkers$cluster)[i],]
  genename <- as.character(gene$gene)
  gene_map <- AnnotationDbi::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
  gene_map <- gene_map[is.na(gene_map[,"ENTREZID"])==F,]
  gene <- gene_map$ENTREZID
  
  KEGG_enrichresult <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
  saveRDS(KEGG_enrichresult, paste0("5.5.Fibr_harmony/Markers/Nomix_", unique(topmarkers$cluster)[i], "_KEGG_enrichresult.Rds"))
  
  KEGG_enrich_all[[i]] <- KEGG_enrichresult
  names(KEGG_enrich_all)[i] <- unique(topmarkers$cluster)[i]
}

saveRDS(KEGG_enrich_all, "5.5.Fibr_harmony/Markers/Nomix_all_KEGG_enrichresult.Rds")

# 12.GSVA -----------------------------------------
library(msigdbr)
KEGG <- msigdbr(species = "Homo sapiens",
                category = "C2",
                subcategory = "CP:KEGG") %>%
  select(gs_name, gene_symbol) %>%
  as.data.frame()
KEGG_list <- split(KEGG$gene_symbol, KEGG$gs_name)
names(KEGG_list) <- gsub("KEGG_", "", names(KEGG_list))
names(KEGG_list) <- gsub("_", " ", names(KEGG_list))
names(KEGG_list) <- str_to_title(names(KEGG_list))

Fibr_nomix <- readRDS("./5.5.Fibr_harmony/Nomix_Fibr_annoted.Rds")

##12.1.SMC 1 vs 2 -----------------------------------

SMC <- subset(Fibr_nomix, celltype %in% c("SMC_KLF2", "SMC_HOPX"))
mtx <- as.matrix(SMC@assays$RNA@data)

GSVA_result <- gsva(expr=mtx, 
                    gset.idx.list=KEGG_list, 
                    kcdf="Gaussian",
                    verbose=T, 
                    parallel.sz = 20)
saveRDS(GSVA_result, "5.5.Fibr_harmony/celltype_GSVA/Nomix_SMC_GSVA_result.Rds")


type <- SMC$celltype
type <- sort(type)

GSVA_SMC <- GSVA_result[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("SMC_KLF2", table(type)[1]), rep("SMC_HOPX", table(type)[2])), levels = c('SMC_KLF2', 'SMC_HOPX'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_SMC)
head(design)

# Tunor VS Normal
compare <- makeContrasts(SMC_HOPX - SMC_KLF2, levels=design)
fit <- lmFit(GSVA_SMC, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 10, wt = t), top_n(dat_plot, -10, wt = t))
dat_plot$threshold = factor(ifelse(dat_plot$t  > 0, 'Up', 'Down'), levels=c('Up','Down','NoSignifi'))

# 排序
dat_plot <- dat_plot %>% arrange(t)

# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

# 绘制
library(ggplot2)
library(ggprism)

p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#A52A2ACC', 'Down'='#053E7ACC')) +
  xlab('') + 
  ylab('t value of GSVA score, SMC_HOPX versus SMC_KLF2') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_SMC <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black')
ggsave("5.5.Fibr_harmony/celltype_GSVA/Nomix_gsva_bar_SMC_top10.pdf", p_SMC, width = 8, height  = 8)

##12.2.myoFibr 1 vs 2 -----------------------------------

MyoFib <- subset(Fibr_nomix, celltype %in% c("MyoFib_RGS5", "MyoFib_CCL21"))
mtx <- as.matrix(MyoFib@assays$RNA@data)

GSVA_result <- gsva(expr=mtx, 
                    gset.idx.list=KEGG_list, 
                    kcdf="Gaussian",
                    verbose=T, 
                    parallel.sz = 20)
saveRDS(GSVA_result, "5.5.Fibr_harmony/celltype_GSVA/Nomix_MyoFib_GSVA_result.Rds")


type <- MyoFib$celltype
type <- sort(type)

GSVA_MyoFib <- GSVA_result[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("MyoFib_RGS5", table(type)[3]), rep("MyoFib_CCL21", table(type)[4])), levels = c('MyoFib_RGS5', 'MyoFib_CCL21'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_MyoFib)
head(design)

# Tunor VS Normal
compare <- makeContrasts(MyoFib_CCL21 - MyoFib_RGS5, levels=design)
fit <- lmFit(GSVA_MyoFib, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 10, wt = t), top_n(dat_plot, -10, wt = t))
dat_plot$threshold = factor(ifelse(dat_plot$t  > 0, 'Up', 'Down'), levels=c('Up','Down','NoSignifi'))

# 排序
dat_plot <- dat_plot %>% arrange(t)

# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

# 绘制
library(ggplot2)
library(ggprism)

p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#A52A2ACC', 'Down'='#053E7ACC')) +
  xlab('') + 
  ylab('t value of GSVA score, MyoFib_CCL21 versus MyoFib_RGS5') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_MyoFib <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black')
ggsave("5.5.Fibr_harmony/celltype_GSVA/Nomix_gsva_bar_MyoFib_top10.pdf", p_MyoFib, width = 8, height  = 8)

##12.3.CAF 1 vs 2 -----------------------------------
"SMC_KLF2", "SMC_HOPX", "MyoFib_RGS5", "MyoFib_CCL21", 
"CAF_APOD", "CAF_CCL2", "Myoblast_MYF5"
CAF <- subset(Fibr_nomix, celltype %in% c("CAF_APOD", "CAF_CCL2"))
mtx <- as.matrix(CAF@assays$RNA@data)

GSVA_result <- gsva(expr=mtx, 
                    gset.idx.list=KEGG_list, 
                    kcdf="Gaussian",
                    verbose=T, 
                    parallel.sz = 20)
saveRDS(GSVA_result, "5.5.Fibr_harmony/celltype_GSVA/Nomix_CAF_GSVA_result.Rds")


type <- CAF$celltype
type <- sort(type)

GSVA_CAF <- GSVA_result[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("CAF_APOD", table(type)[5]), rep("CAF_CCL2", table(type)[6])), levels = c('CAF_APOD', 'CAF_CCL2'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_CAF)
head(design)

# Tunor VS Normal
compare <- makeContrasts(CAF_CCL2 - CAF_APOD, levels=design)
fit <- lmFit(GSVA_CAF, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 10, wt = t), top_n(dat_plot, -10, wt = t))
dat_plot$threshold = factor(ifelse(dat_plot$t  > 0, 'Up', 'Down'), levels=c('Up','Down','NoSignifi'))

# 排序
dat_plot <- dat_plot %>% arrange(t)

# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

# 绘制
library(ggplot2)
library(ggprism)

p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#A52A2ACC', 'Down'='#053E7ACC')) +
  xlab('') + 
  ylab('t value of GSVA score, CAF-2 versus CAF-1') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_y_continuous(limits = c (-15, 25))
p

# 添加标签

p_CAF <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black')
ggsave("5.5.Fibr_harmony/celltype_GSVA/Nomix_gsva_bar_CAF_top10.pdf", p_CAF, width = 8, height  = 8)

# 12.4. all plot ---------------------------------------

pdf("5.5.Fibr_harmony/celltype_GSVA/Nomix_gsva_bar_all_top10.pdf", width = 24, height = 8)
gridExtra::grid.arrange(grobs = list(p_SMC,p_MyoFib,p_CAF), ncol = 3)
dev.off()


# 13. CAF_APOD score in TCGA --------------------------------------

#### read gene matrix and MP score --------------------
expMtx <- read.table("./0.TCGA/mRNA_Matrix_TPM.txt", row.names = 1, header = T, sep = "\t")
expMtx <- as.matrix(expMtx)

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_TCGA.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)

clinical <- read.table("./0.TCGA/Clinical.txt", row.names = 1, header = T, sep = "\t")
age <- as.data.frame(clinical[,1])
colnames(age) <- "age"
rownames(age) <- gsub("_", ".", rownames(clinical))
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"

survival <- read.table("./0.TCGA/TCGA-PRAD.survival.tsv", row.names = 1, header = T, sep = "\t")
rownames(survival) <- gsub("-", ".", rownames(survival))
survival <- survival[,c(1,3)]

#### calculate CAF_APOD score ------------------------------

Allmarkers <- readRDS("./5.5.Fibr_harmony/Markers/Nomix_celltype_markers.Rds")

APOD <- Allmarkers[Allmarkers$cluster == "CAF_APOD",]
APOD_markers <- APOD[which(APOD$p_val_adj < 0.01 & APOD$pct.1 > 0.3 & APOD$avg_log2FC > 1),]

CCL2 <- Allmarkers[Allmarkers$cluster == "CAF_CCL2",]
CCL2_markers <- CCL2[which(CCL2$p_val_adj < 0.01 & CCL2$pct.1 > 0.3 & CCL2$avg_log2FC > 1),]

library(AUCell)
gs <- list(CAF_APOD = APOD_markers$gene,
           CAF_CCL2 = CCL2_markers$gene)

cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.5.Fibr_harmony/Bulk_CAF_score/CAF_gs_AUC.Rds")

CAF <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

##### MP1 and APOD --------------------------------------------------------

MP1_CAF <- merge(t(AUCell_mtx), CAF, by = 0) %>% column_to_rownames("Row.names")
MP1_CAF <- MP1_CAF[-grep(".11$|.06$", rownames(MP1_CAF)),]

for (i in (7 : (ncol(MP1_CAF)))) {
  
  P <- ggplot(MP1_CAF, aes(x=MP1_CAF[,1], y=MP1_CAF[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = MP1_CAF[,1], y = MP1_CAF[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(MP1_CAF)[1]," AUCell")) + ylab(paste0(colnames(MP1_CAF)[i]," AUCell"))
  
  ggsave(paste0("./5.5.Fibr_harmony/Bulk_CAF_score/CorPlot_MP1_",colnames(MP1_CAF)[i],".pdf"), P, width = 2, height = 2)
}

##### EO vs LO and CAF --------------------------------------------------------

Age_CAF <- merge(age, CAF, by = 0) %>% column_to_rownames("Row.names")

# correlation

for (i in (3 : (ncol(Age_CAF)))) {
  
  P <- ggplot(Age_CAF, aes(x=Age_CAF[,1], y=Age_CAF[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = Age_CAF[,1], y = Age_CAF[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_CAF)[i]," AUCell"))
  
  ggsave(paste0("./5.5.Fibr_harmony/Bulk_CAF_score/CorPlot_Age_",colnames(Age_CAF)[i],".pdf"), P, width = 2, height = 2)
}

for (i in (3 : (ncol(Age_CAF)))) {
  
  P <- ggplot(Age_CAF, aes(x=Age_CAF[,1], y=Age_CAF[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'pearson', aes(x = Age_CAF[,1], y = Age_CAF[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_CAF)[i]," AUCell"))
  
  ggsave(paste0("./5.5.Fibr_harmony/Bulk_CAF_score/CorPlot_pearson_Age_",colnames(Age_CAF)[i],".pdf"), P, width = 2, height = 2)
}

# boxplot

data <- melt(Age_CAF[,c(2:ncol(Age_CAF))])
P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("AUCell score") +
  xlab("") + 
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("./5.5.Fibr_harmony/Bulk_CAF_score/BoxPlot_CAF_EOLO.pdf", width = 2, height = 2)
P + NoLegend()
dev.off()

#### CAF and survival --------------------------------------

survival_CAF <- merge(CAF, survival, by = 0) %>% column_to_rownames("Row.names")
survival_CAF <- survival_CAF[-grep(".11$|.06$", rownames(survival_CAF)),]
survival_CAF$OS.time <-survival_CAF$OS.time / 365

coxph(Surv(OS.time, OS) ~ CAF_CAF, survival_CAF)

library(survival)
library("survminer")

# APOD

quar <- quantile(survival_CAF$CAF_APOD, probs = seq(0, 1, 0.25))

survival_CAF$APOD_HL <- 0
survival_CAF$APOD_HL[which(survival_CAF$CAF_APOD > quar[2])] <- 1

diff <- survdiff(Surv(OS.time, OS) ~ APOD_HL, data = survival_CAF)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
P1 <- pValue
fit <- survfit(Surv(OS.time, OS) ~ APOD_HL, data = survival_CAF)

ggsurvplot(fit, 
           data=survival_CAF,
           conf.int=F,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=F,
           size = 1,
           ylim=c(0,1.00),
           legend.title= "",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("#000080", "#A00000"),
           risk.table.height=.25)

# CCL2

quar <- quantile(survival_CAF$CAF_CCL2, probs = seq(0, 1, 0.25))

survival_CAF$CCL2_HL <- 0
survival_CAF$CCL2_HL[which(survival_CAF$CAF_CCL2 > quar[3])] <- 1
diff <- survdiff(Surv(OS.time, OS) ~ CCL2_HL, data = survival_CAF)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
P1 <- pValue
fit <- survfit(Surv(OS.time, OS) ~ CCL2_HL, data = survival_CAF)

ggsurvplot(fit, 
           data=survival_CAF,
           conf.int=F,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=F,
           size = 1,
           ylim=c(0,1.00),
           legend.title= "",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("#000080", "#A00000"),
           risk.table.height=.25)

#### CAF and survival in LOPC --------------------------------------

survival_Age <- merge(survival, age, by = 0) %>% column_to_rownames("Row.names")
survival_EOPC <- survival_Age[survival_Age$Type == "LOPC",]
survival_CAF <- merge(CAF, survival_EOPC, by = 0) %>% column_to_rownames("Row.names")
survival_CAF$OS.time <-survival_CAF$OS.time / 365


library(survival)
library("survminer")

# APOD

quar <- quantile(survival_CAF$CAF_APOD, probs = seq(0, 1, 0.25))

survival_CAF$APOD_HL <- 0
survival_CAF$APOD_HL[which(survival_CAF$CAF_APOD > quar[2])] <- 1

diff <- survdiff(Surv(OS.time, OS) ~ APOD_HL, data = survival_CAF)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
P1 <- pValue
fit <- survfit(Surv(OS.time, OS) ~ APOD_HL, data = survival_CAF)

ggsurvplot(fit, 
           data=survival_CAF,
           conf.int=F,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=F,
           size = 1,
           ylim=c(0,1.00),
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("#053E7A", "#A52A2A"),
           legend = "none",
           font.x = c(16, "plain", "black"),
           font.y = c(16, "plain", "black"))

# CCL2

quar <- quantile(survival_CAF$CAF_CCL2, probs = seq(0, 1, 0.25))

survival_CAF$CCL2_HL <- 0
survival_CAF$CCL2_HL[which(survival_CAF$CAF_CCL2 > quar[2])] <- 1
diff <- survdiff(Surv(OS.time, OS) ~ CCL2_HL, data = survival_CAF)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
P1 <- pValue
fit <- survfit(Surv(OS.time, OS) ~ CCL2_HL, data = survival_CAF)

ggsurvplot(fit, 
           data=survival_CAF,
           conf.int=F,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=F,
           size = 1,
           ylim=c(0,1.00),
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("#053E7A", "#A52A2A"),
           legend = "none",
           font.x = c(16, "plain", "black"),
           font.y = c(16, "plain", "black"))

#### CAF and TNM --------------------------------------------

#M
clinical$clinical_M[grep("M1", clinical$clinical_M)] <- "M1"

#GS
clinical$gleason_score <- factor(clinical$gleason_score, 
                                 levels = c("6","7","8","9","10"))

rownames(clinical) <- gsub("_", ".", rownames(clinical))

mtx <- merge(CAF, clinical, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,c(1,2,6,8,9,10)]

pairwise.wilcox.test(mtx$CAF_APOD, mtx$clinical_M, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_APOD, mtx$pathologic_T, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_APOD, mtx$pathologic_N, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_CCL2, mtx$clinical_M, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_CCL2, mtx$pathologic_T, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_CCL2, mtx$pathologic_N, p.adjust.method = "none")

##### plot ----------------------------------------

for (i in 3:ncol(mtx)) {
  data <- mtx[, c(1,i)]
  NAs <- which(is.na(data[,2]))
  if(length(NAs) > 0) {data <- data[-NAs,]} else {data <- data}
  levels <- length(table(data[,2]))
  col <- paste0(color,"CC")[1:levels]
  
  P = ggplot(data, aes(x = data[,2], y = data[,1], fill = data[,2])) +
    scale_fill_manual(values = col) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(colnames(mtx)[i]) +
    ylab(paste0(colnames(data)[1], " AUCell")) +
    xlab("") +
    NoLegend()
  
  ggsave(filename = paste0("5.5.Fibr_harmony/Bulk_CAF_score/CAF_APOD_CliBox_",  colnames(mtx)[i], ".pdf"),
         P,
         height = 3,
         width = 1 + 0.3*1*length(table(mtx[,i])))
  
}

for (i in 3:ncol(mtx)) {
  data <- mtx[, c(2,i)]
  NAs <- which(is.na(data[,2]))
  if(length(NAs) > 0) {data <- data[-NAs,]} else {data <- data}
  levels <- length(table(data[,2]))
  col <- paste0(color,"CC")[1:levels]
  
  P = ggplot(data, aes(x = data[,2], y = data[,1], fill = data[,2])) +
    scale_fill_manual(values = col) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(colnames(mtx)[i]) +
    ylab(paste0(colnames(data)[1], " AUCell")) +
    xlab("") +
    NoLegend()
  
  ggsave(filename = paste0("5.5.Fibr_harmony/Bulk_CAF_score/CAF_CCL2_CliBox_",  colnames(mtx)[i], ".pdf"),
         P,
         height = 3,
         width = 1 + 0.3*1*length(table(mtx[,i])))
  
}

# Ptrend

GS <- mtx[,c(1,4)]
GS$gleason_score <- as.numeric(GS$gleason_score)
cor.test(GS$CAF_APOD, GS$gleason_score, method = "pearson")

pT <- mtx[,c(1,6)]
pT$pathologic_T <- factor(pT$pathologic_T, levels = names(table(pT$pathologic_T)))
pT$pathologic_T <- as.numeric(pT$pathologic_T)
cor.test(pT$CAF_APOD, pT$pathologic_T, method = "pearson")

GS <- mtx[,c(2,4)]
GS$gleason_score <- as.numeric(GS$gleason_score)
cor.test(GS$CAF_CCL2, GS$gleason_score, method = "pearson")

pT <- mtx[,c(2,6)]
pT$pathologic_T <- factor(pT$pathologic_T, levels = names(table(pT$pathologic_T)))
pT$pathologic_T <- as.numeric(pT$pathologic_T)
cor.test(pT$CAF_CCL2, pT$pathologic_T, method = "pearson")

##### CAF avg plot -------------------------------

mtx$CAF <- (mtx$CAF_APOD + mtx$CAF_CCL2) / 2

for (i in 3:(ncol(mtx)-1)) {
  data <- mtx[, c(7,i)]
  NAs <- which(is.na(data[,2]))
  if(length(NAs) > 0) {data <- data[-NAs,]} else {data <- data}
  levels <- length(table(data[,2]))
  col <- paste0(color,"CC")[1:levels]
  
  P = ggplot(data, aes(x = data[,2], y = data[,1], fill = data[,2])) +
    scale_fill_manual(values = col) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(colnames(mtx)[i]) +
    ylab("Avg. CAF AUCell") +
    xlab("") +
    NoLegend()
  
  ggsave(filename = paste0("5.5.Fibr_harmony/Bulk_CAF_score/CAFs_avg_CliBox_",  colnames(mtx)[i], ".pdf"),
         P,
         height = 3,
         width = 1 + 0.3*1*length(table(mtx[,i])))
  
}

# Ptrend

GS <- mtx[,c(7,4)]
GS$gleason_score <- as.numeric(GS$gleason_score)
cor.test(GS$CAF, GS$gleason_score, method = "pearson")

pT <- mtx[,c(7,6)]
pT$pathologic_T <- factor(pT$pathologic_T, levels = names(table(pT$pathologic_T)))
pT$pathologic_T <- as.numeric(pT$pathologic_T)
cor.test(pT$CAF, pT$pathologic_T, method = "pearson")

#### BMPs and clinical --------------------------

BMPs <- expMtx[c("BMP4","BMP5","BMP7"),]
BMPs <- log2(BMPs + 1)
mtx <- merge(t(BMPs), clinical, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,c(1,2,3,7,9,10,11)]

for (i in 4:ncol(mtx)) {
  data <- mtx[, c(1,i)]
  NAs <- which(is.na(data[,2]))
  if(length(NAs) > 0) {data <- data[-NAs,]} else {data <- data}
  levels <- length(table(data[,2]))
  col <- paste0(color,"CC")[1:levels]
  
  P = ggplot(data, aes(x = data[,2], y = data[,1], fill = data[,2])) +
    scale_fill_manual(values = col) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(colnames(mtx)[i]) +
    ylab(paste0(colnames(data)[1], " exp.")) +
    xlab("") +
    NoLegend()
  
  ggsave(filename = paste0("5.5.Fibr_harmony/Bulk_CAF_score/BMP4_",  colnames(mtx)[i], ".pdf"),
         P,
         height = 3,
         width = 1 + 0.3*1*length(table(mtx[,i])))
  
}

mtx$BMP <- (mtx$BMP4 + mtx$BMP5 + mtx$BMP7) / 3

for (i in 4:(ncol(mtx)-1)) {
  data <- mtx[, c(8,i)]
  NAs <- which(is.na(data[,2]))
  if(length(NAs) > 0) {data <- data[-NAs,]} else {data <- data}
  levels <- length(table(data[,2]))
  col <- paste0(color,"CC")[1:levels]
  
  P = ggplot(data, aes(x = data[,2], y = data[,1], fill = data[,2])) +
    scale_fill_manual(values = col) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(colnames(mtx)[i]) +
    ylab("BMPs Avg. Exp.") +
    xlab("") +
    NoLegend()
  
  ggsave(filename = paste0("5.5.Fibr_harmony/Bulk_CAF_score/BMPs_",  colnames(mtx)[i], ".pdf"),
         P,
         height = 3,
         width = 1 + 0.3*1*length(table(mtx[,i])))
  
}

pT <- mtx[,c(1,7)]
pT$pathologic_T <- factor(pT$pathologic_T, levels = names(table(pT$pathologic_T)))
pT$pathologic_T <- as.numeric(pT$pathologic_T)
cor.test(pT$BMP4, pT$pathologic_T, method = "pearson")

GS <- mtx[,c(1,5)]
GS$gleason_score <- as.numeric(GS$gleason_score)
cor.test(GS$BMP4, GS$gleason_score, method = "pearson")

#### CAF and EMT score --------------------------------------------

EMT <- readRDS("5.1.1.epithelial_NMF/process_k3/Bulk validation/Hypo_EMT_GSVA.Rds")
CAF_EMT <- merge(CAF, t(EMT), by=0) %>% column_to_rownames("Row.names")
CAF_EMT <- CAF_EMT[-grep(".11$|.06$", rownames(CAF_EMT)),]

for (i in c(1 : 2)) {
  
  P <- ggplot(CAF_EMT, aes(x=CAF_EMT[,3], y=CAF_EMT[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = CAF_EMT[,3], y = CAF_EMT[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("EMT GSVA") + ylab(paste0(colnames(CAF_EMT)[i]," AUCell"))
  
  ggsave(paste0("./5.5.Fibr_harmony/Bulk_CAF_score/CorPlot_EMT_",colnames(CAF_EMT)[i],".pdf"), P, width = 2, height = 2)
}

APOD_EMT <- CAF_EMT[,c(1,3)]
APOD_EMT$Type <- "CAF_APOD"
colnames(APOD_EMT)[1] <- "AUCell"
CCL2_EMT <- CAF_EMT[,c(2,3)]
CCL2_EMT$Type <- "CAF_CCL2"
colnames(CCL2_EMT)[1] <- "AUCell"

data <- rbind(APOD_EMT,CCL2_EMT)

pdf("./5.5.Fibr_harmony/Bulk_CAF_score/CorPlot_EMT_CAF.pdf", width = 2, height = 2)
pal <- c(CAF_APOD = "#CD7054CC", CAF_CCL2 = "#4A7088CC")
ggplot(data = data, aes(x=AUCell, y=EMT, color=Type)) +
  geom_point(size = 0.5) +
  theme_bw() + 
  scale_color_manual(values = pal) +
  geom_smooth(method = lm, formula = y ~ x) +
  theme(axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        axis.title.x=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
  NoLegend()
dev.off()

### BMPs and EMT score ------------------------------------

BMPs <- expMtx[c("BMP4","BMP5","BMP7"),]
BMPs <- log2(BMPs + 1)
BMPs_EMT <- merge(t(BMPs), t(EMT), by=0) %>% column_to_rownames("Row.names")
BMPs_EMT <- BMPs_EMT[-grep(".11$|.06$", rownames(BMPs_EMT)),]

for (i in c(1 : 3)) {
  
  P <- ggplot(BMPs_EMT, aes(x=BMPs_EMT[,4], y=BMPs_EMT[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = BMPs_EMT[,4], y = BMPs_EMT[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("EMT GSVA") + ylab(paste0(colnames(BMPs_EMT)[i]," exp."))
  
  ggsave(paste0("./5.5.Fibr_harmony/Bulk_CAF_score/CorPlot_EMT_",colnames(BMPs_EMT)[i],".pdf"), P, width = 2, height = 2)
}

BMP4_EMT <- BMPs_EMT[,c(1,4)]
BMP4_EMT$Type <- "BMP4"
colnames(BMP4_EMT)[1] <- "Norm.Exp."
BMP5_EMT <- BMPs_EMT[,c(2,4)]
BMP5_EMT$Type <- "BMP5"
colnames(BMP5_EMT)[1] <- "Norm.Exp."
BMP7_EMT <- BMPs_EMT[,c(3,4)]
BMP7_EMT$Type <- "BMP7"
colnames(BMP7_EMT)[1] <- "Norm.Exp."

data <- rbind(BMP4_EMT,BMP5_EMT,BMP7_EMT)

pdf("./5.5.Fibr_harmony/Bulk_CAF_score/CorPlot_EMT_BMPs.pdf", width = 2, height = 2)
pal <- c(BMP4 = "#CD7054CC", BMP5 = "#CDAD00CC", BMP7 = "#4A7088CC")
ggplot(data = data, aes(x=Norm.Exp., y=EMT, color=Type)) +
  geom_point(size = 0.5) +
  theme_bw() + 
  scale_color_manual(values = pal) +
  geom_smooth(method = lm, formula = y ~ x) +
  theme(axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        axis.title.x=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
  NoLegend()
dev.off()

### CAF EMT Age and clinical -------------------------------------

#T
clinical$pathologic_T[grep("T2", clinical$pathologic_T)] <- "T2"
clinical$pathologic_T[grep("T3", clinical$pathologic_T)] <- "T3"

#M
clinical$clinical_M[grep("M1", clinical$clinical_M)] <- "M1"

#GS
clinical$gleason_score <- factor(clinical$gleason_score, 
                                 levels = c("6","7","8","9","10"))

rownames(clinical) <- gsub("_", ".", rownames(clinical))

mtx <- merge(CAF_EMT, clinical, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,c(1,2,3,8,10,11,12)]

mtx <- merge(mtx, age, by = 0) %>% column_to_rownames("Row.names")
mtx_LO <- mtx[mtx$Type == "LOPC",]

pairwise.wilcox.test(mtx$EMT, mtx$clinical_M, p.adjust.method = "none")
pairwise.wilcox.test(mtx$EMT, mtx$pathologic_T, p.adjust.method = "none")
pairwise.wilcox.test(mtx$EMT, mtx$pathologic_N, p.adjust.method = "none")

pairwise.wilcox.test(mtx_LO$EMT, mtx_LO$clinical_M, p.adjust.method = "none")
pairwise.wilcox.test(mtx_LO$EMT, mtx_LO$pathologic_T, p.adjust.method = "none")
pairwise.wilcox.test(mtx_LO$EMT, mtx_LO$pathologic_N, p.adjust.method = "none")

pairwise.wilcox.test(mtx$CAF_APOD, mtx$clinical_M, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_APOD, mtx$pathologic_T, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_APOD, mtx$pathologic_N, p.adjust.method = "none")

pairwise.wilcox.test(mtx_LO$CAF_APOD, mtx_LO$clinical_M, p.adjust.method = "none")
pairwise.wilcox.test(mtx_LO$CAF_APOD, mtx_LO$pathologic_T, p.adjust.method = "none")
pairwise.wilcox.test(mtx_LO$CAF_APOD, mtx_LO$pathologic_N, p.adjust.method = "none")

pairwise.wilcox.test(mtx$CAF_CCL2, mtx$clinical_M, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_CCL2, mtx$pathologic_T, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_CCL2, mtx$pathologic_N, p.adjust.method = "none")

pairwise.wilcox.test(mtx_LO$CAF_CCL2, mtx_LO$clinical_M, p.adjust.method = "none")
pairwise.wilcox.test(mtx_LO$CAF_CCL2, mtx_LO$pathologic_T, p.adjust.method = "none")
pairwise.wilcox.test(mtx_LO$CAF_CCL2, mtx_LO$pathologic_N, p.adjust.method = "none")

for (i in 4:7) {
  data <- mtx[, c(3,i)]
  NAs <- which(is.na(data[,2]))
  if(length(NAs) > 0) {data <- data[-NAs,]} else {data <- data}
  levels <- length(table(data[,2]))
  col <- paste0(color,"CC")[1:levels]
  
  P = ggplot(data, aes(x = data[,2], y = data[,1], fill = data[,2])) +
    scale_fill_manual(values = col) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(colnames(mtx)[i]) +
    ylab("EMT GSVA") +
    xlab("") +
    NoLegend()
  
  ggsave(filename = paste0("5.5.Fibr_harmony/Bulk_CAF_score/EMT_",  colnames(mtx)[i], ".pdf"),
         P,
         height = 3,
         width = 1 + 0.3*1*length(table(mtx[,i])))
  
}


# Ptrend

GS <- mtx[,c(3,5)]
GS$gleason_score <- as.numeric(GS$gleason_score)
cor.test(GS$EMT, GS$gleason_score, method = "pearson")

pT <- mtx[,c(3,7)]
pT$pathologic_T <- factor(pT$pathologic_T, levels = names(table(pT$pathologic_T)))
pT$pathologic_T <- as.numeric(pT$pathologic_T)
cor.test(pT$EMT, pT$pathologic_T, method = "pearson")

# 13.2. GSE21034 CAF EMT Age and clinical -------------------------------------

expMtx <- read.table("./0.TCGA/GSE21034/MergeExpro_contrib1-GPL5188.txt", header = T, sep = "\t", check.names = F)
expMtx <- as.matrix(expMtx)
rownames(expMtx) <- expMtx[,1]
expMtx <- expMtx[,2:ncol(expMtx)]
dimnames <- list(rownames(expMtx),colnames(expMtx))
expMtx <- matrix(as.numeric(as.matrix(expMtx)),nrow=nrow(expMtx),dimnames=dimnames)
expMtx <- avereps(expMtx)
expMtx <- log2(expMtx + 1)

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE21034.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)

clinical <- read.table("./0.TCGA/GSE21034/GSE21034_clinical_simplified_untreated.txt", row.names = 1, header = T, sep = "\t")
age <- as.data.frame(clinical[,2])
colnames(age) <- "age"
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"

### calculate CAF_APOD score ------------------------------

Allmarkers <- readRDS("./5.5.Fibr_harmony/Markers/Nomix_celltype_markers.Rds")

APOD <- Allmarkers[Allmarkers$cluster == "CAF_APOD",]
APOD_markers <- APOD[which(APOD$p_val_adj < 0.01 & APOD$pct.1 > 0.3 & APOD$avg_log2FC > 1),]

CCL2 <- Allmarkers[Allmarkers$cluster == "CAF_CCL2",]
CCL2_markers <- CCL2[which(CCL2$p_val_adj < 0.01 & CCL2$pct.1 > 0.3 & CCL2$avg_log2FC > 1),]

library(AUCell)
gs <- list(CAF_APOD = APOD_markers$gene,
           CAF_CCL2 = CCL2_markers$gene)

cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.5.Fibr_harmony/Bulk_CAF_score/GSE21034/CAF_gs_AUC.Rds")

CAF <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

### calculate EMT score ------------------------------

cancersea <- read.gmt("hallmark_cancersea.gmt")

select_gs <- list(EPITHELIAL_MESENCHYMAL_TRANSITION = cancersea$gene[cancersea$term == "EPITHELIAL_MESENCHYMAL_TRANSITION"])

library(GSVA)

GSVA <- gsva(expMtx, select_gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "./5.5.Fibr_harmony/Bulk_CAF_score/GSE21034/EMT_GSVA.Rds")

### CAF and clinical --------------------------------------------

mtx <- merge(CAF, clinical, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,c(1,2,6,9,10)]

pairwise.wilcox.test(mtx$CAF_APOD, mtx$SMS, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_APOD, mtx$ECE, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_APOD, mtx$SVI, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_APOD, mtx$LNI, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_CCL2, mtx$SMS, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_CCL2, mtx$ECE, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_CCL2, mtx$SVI, p.adjust.method = "none")
pairwise.wilcox.test(mtx$CAF_CCL2, mtx$LNI, p.adjust.method = "none")

### EMT and clinical --------------------------------------------

EMT <- t(GSVA) %>% as.data.frame()
colnames(EMT) <- "EMT"
mtx <- merge(EMT, clinical, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,c(1,9,10,11,12)]

pairwise.wilcox.test(mtx$EMT, mtx$SMS, p.adjust.method = "none")
pairwise.wilcox.test(mtx$EMT, mtx$ECE, p.adjust.method = "none")
pairwise.wilcox.test(mtx$EMT, mtx$SVI, p.adjust.method = "none")
pairwise.wilcox.test(mtx$EMT, mtx$LNI, p.adjust.method = "none")

for (i in 2:5) {
  data <- mtx[, c(1,i)]
  NAs <- which(is.na(data[,2]))
  if(length(NAs) > 0) {data <- data[-NAs,]} else {data <- data}
  levels <- length(table(data[,2]))
  col <- paste0(color,"CC")[1:levels]
  
  P = ggplot(data, aes(x = data[,2], y = data[,1], fill = data[,2])) +
    scale_fill_manual(values = col) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(colnames(mtx)[i]) +
    ylab("EMT GSVA") +
    xlab("") +
    NoLegend()
  
  ggsave(filename = paste0("./5.5.Fibr_harmony/Bulk_CAF_score/GSE21034/EMT_",  colnames(mtx)[i], ".pdf"),
         P,
         height = 3,
         width = 1 + 0.3*1*length(table(mtx[,i])))
  
}


# Ptrend

GS <- mtx[,c(3,5)]
GS$gleason_score <- as.numeric(GS$gleason_score)
cor.test(GS$EMT, GS$gleason_score, method = "pearson")

pT <- mtx[,c(3,7)]
pT$pathologic_T <- factor(pT$pathologic_T, levels = names(table(pT$pathologic_T)))
pT$pathologic_T <- as.numeric(pT$pathologic_T)
cor.test(pT$EMT, pT$pathologic_T, method = "pearson")

# 14. Bulk validation: CRPC vs CSPC ----------------------------------------

Allmarkers <- readRDS("./5.5.Fibr_harmony/Markers/Nomix_celltype_markers.Rds")

APOD <- Allmarkers[Allmarkers$cluster == "CAF_APOD",]
APOD_markers <- APOD[which(APOD$p_val_adj < 0.01 & APOD$pct.1 > 0.3 & APOD$avg_log2FC > 1),]

CCL2 <- Allmarkers[Allmarkers$cluster == "CAF_CCL2",]
CCL2_markers <- CCL2[which(CCL2$p_val_adj < 0.01 & CCL2$pct.1 > 0.3 & CCL2$avg_log2FC > 1),]

Markers <- readRDS("./5.1.1.epithelial_NMF/process_k3/MP_markers.Rds")
ARMP <- Markers[[1]]

gs <- list(ARMP = ARMP,
           iCAF_APOD = APOD_markers$gene,
           iCAF_CCL2 = CCL2_markers$gene)

## GSE189343 (为转移组织，只能看ARMP) -----------------------------------------------

exp <- read.table("0.TCGA/GSE189343/GSE189343_mRNA.txt", header = T, row.names = 1, sep = "\t")
exp <- as.matrix(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)

saveRDS(cells_AUC, "5.5.Fibr_harmony/Bulk_CAF_score/GSE189343/ARMP_iCAF_gs_AUC.Rds")

score <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

clinical <- read.table("0.TCGA/GSE189343/GSE189343_SampleInfo.txt", header = T, row.names = 1, sep = "\t")

mtx <- merge(clinical, score, by=0) %>% column_to_rownames("Row.names")
mtx <- mtx[,c("X1.TREATMENT.GROUP","ARMP")]
data <- melt(mtx)
colnames(data)[c(1,2)] <- c("Group", "ARMP")
data$Group[data$Group == "Castration resistant"] <- "CR"
data$Group[data$Group == "Short-term castrated"] <- "Castrated"
data$Group[data$Group == "Hormone-naïve"] <- "HN"
data$Group <- factor(data$Group, levels = c("CR","Castrated","HN"))

#plot

P = ggplot(data, aes(x = Group, y = value, fill = Group)) +
  scale_fill_manual(values = c("#A52A2ACC", "#E46B4FCC", "#053E7ACC")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  ylab("AR-MP AUCell score") +
  xlab("") + 
  ggtitle("GSE189343 \n (metastasis tissues)") +
  stat_compare_means(aes(group=Group),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif") +
  NoLegend()

pdf("5.5.Fibr_harmony/Bulk_CAF_score/GSE189343/boxplot_ARMP_AUC_diff.pdf", width=2, height=3)
print(P)
dev.off()

## GSE32269(CRPC为转移，只看ARMP) --------------------------------------------

exp <- read.table("0.TCGA/GSE32269/GSE32269_exp.txt", header = T, row.names = 1, sep = "\t")
exp <- as.matrix(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)

saveRDS(cells_AUC, "5.5.Fibr_harmony/Bulk_CAF_score/GSE32269/ARMP_iCAF_gs_AUC.Rds")

score <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

clinical <- read.table("0.TCGA/GSE32269/GSE32269_SampleInfo.txt", header = T, row.names = 1, sep = "\t")

mtx <- merge(clinical, score, by=0) %>% column_to_rownames("Row.names")
mtx <- mtx[,c("X1.METASTASIS","ARMP")]
data <- melt(mtx)
colnames(data)[c(1,2)] <- c("Group", "ARMP")
data <- data[which(data$Group != "normal"),]
data$Group[data$Group == "bone metastases (CRPC)"] <- "CRPC Bone Met"
data$Group[data$Group == "localized prostate cancer"] <- "HNPC Localized"
data$Group <- factor(data$Group, levels = c("CRPC Bone Met","HNPC Localized"))

#plot

P = ggplot(data, aes(x = Group, y = value, fill = Group)) +
  scale_fill_manual(values = c("#A52A2ACC", "#053E7ACC")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  ylab("AR-MP AUCell score") +
  xlab("") + 
  ggtitle("GSE32269") +
  stat_compare_means(aes(group=Group),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif") +
  NoLegend()

pdf("5.5.Fibr_harmony/Bulk_CAF_score/GSE32269/boxplot_ARMP_AUC_diff.pdf", width=1.5, height=3)
print(P)
dev.off()

## GSE70768 ------------------------------------------------

library(limma)
expMtx <- read.table("./0.TCGA/GSE70768/matrix_GSE70768.txt", header = T, sep = "\t", check.names = F)
expMtx <- as.matrix(expMtx)
rownames(expMtx) <- expMtx[,1]
expMtx <- expMtx[,2:ncol(expMtx)]
dimnames <- list(rownames(expMtx),colnames(expMtx))
expMtx <- matrix(as.numeric(as.matrix(expMtx)),nrow=nrow(expMtx),dimnames=dimnames)
expMtx <- avereps(expMtx)

clinical <- read.table("./0.TCGA/GSE70768/clinical.txt", row.names = 1, header = T, sep = "\t")
clinical$Type[clinical$Age <= 55] <- "EOPC"
clinical$Type[clinical$Age > 55] <- "LOPC"

clinical <- clinical[clinical$Sample != "Benign",]
expMtx <- expMtx[,which(colnames(expMtx) %in% rownames(clinical))]

### calculate CAF_APOD score ------------------------------

cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.5.Fibr_harmony/Bulk_CAF_score/GSE70768/ARMP_iCAF_gs_AUC.Rds")

score <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()
mtx <- merge(score, clinical,by=0) %>% column_to_rownames("Row.names")

### ARMP and CRPC ------------------------------------

data <- mtx[,c(1,16)]
data <- melt(data)

P = ggplot(data, aes(x = Sample, y = value, fill = Sample)) +
  scale_fill_manual(values = c("#A52A2ACC", "#053E7ACC")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  ylab("AR-MP AUCell score") +
  xlab("") + 
  ggtitle("GSE70768") +
  stat_compare_means(aes(group=Sample),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif") +
  NoLegend()

pdf("5.5.Fibr_harmony/Bulk_CAF_score/GSE70768/boxplot_ARMP_AUC_diff.pdf", width=1.5, height=3)
print(P)
dev.off()

### CAF and CRPC ------------------------------------

data <- mtx[,c(2,3,16)]
data <- melt(data)

P = ggplot(data, aes(x = variable, y = value, fill = Sample)) +
  scale_fill_manual(values = c("#A52A2ACC", "#053E7ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_blank()) +
  ylab("AUCell") +
  stat_compare_means(aes(group=Sample),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

pdf("5.5.Fibr_harmony/Bulk_CAF_score/GSE70768/boxplot_iCAF_AUC_diff.pdf", width=3, height=3)
print(P)
dev.off()

### CAF and EO/LO ------------------------------------

data <- mtx[,c(2,3,17)]
data <- melt(data)

P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#A52A2ACC", "#053E7ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_blank()) +
  ylab("AUCell") +
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

pdf("5.5.Fibr_harmony/Bulk_CAF_score/GSE70768/boxplot_iCAF_AUC_diff.pdf", width=3, height=3)
print(P)
dev.off()

### calculate EMT score ------------------------------

cancersea <- read.gmt("hallmark_cancersea.gmt")

select_gs <- list(EPITHELIAL_MESENCHYMAL_TRANSITION = cancersea$gene[cancersea$term == "EPITHELIAL_MESENCHYMAL_TRANSITION"])

library(GSVA)

GSVA <- gsva(expMtx, select_gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "./5.5.Fibr_harmony/Bulk_CAF_score/GSE70768/EMT_GSVA.Rds")

### EMT and ECE ---------------------------------------

mtx <- merge(t(GSVA), clinical,by=0) %>% column_to_rownames("Row.names")

data <- mtx[, c(1,12)]
data <- data[-which(data$ECE == "unknown"),]
levels <- length(table(data[,2]))
col <- paste0(color,"CC")[1:levels]

P = ggplot(data, aes(x = data[,2], y = data[,1], fill = data[,2])) +
  scale_fill_manual(values = col) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("ECE") +
  ylab("EMT GSVA") +
  xlab("") +
  NoLegend()

## Vo et al. (CRPC为转移，只看ARMP) --------------------------------------------

exp <- read.table("0.TCGA/Cell/ExpMtx.txt", header = T, row.names = 1, sep = "\t")
exp <- as.matrix(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)

saveRDS(cells_AUC, "5.5.Fibr_harmony/Bulk_CAF_score/Cell/ARMP_iCAF_gs_AUC.Rds")

score <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

clinical <- read.table("0.TCGA/Cell/sample.txt", header = T, row.names = 1, sep = "\t")

mtx <- merge(clinical, score, by=0) %>% column_to_rownames("Row.names")
mtx <- mtx[which(mtx$Type %in% c("CRPC","CSPC")),]

### ARMP and CRPC ------------------------------------

data <- mtx[,c(1,2)]
data <- melt(data)

P = ggplot(data, aes(x = Type, y = value, fill = Type)) +
  scale_fill_manual(values = c("#A52A2ACC", "#053E7ACC")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  ylab("AR-MP AUCell score") +
  xlab("") + 
  ggtitle("Vo et al.") +
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif") +
  NoLegend()

pdf("5.5.Fibr_harmony/Bulk_CAF_score/Cell/boxplot_ARMP_AUC_diff.pdf", width=1.5, height=3)
print(P)
dev.off()

### CAF and CRPC ------------------------------------

data <- mtx[,c(1,3,4)]
data <- melt(data)

P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#A52A2ACC", "#053E7ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_blank()) +
  ylab("AR-MP AUCell score") + 
  ggtitle("Vo et al.") +
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

pdf("5.5.Fibr_harmony/Bulk_CAF_score/Cell/boxplot_iCAF_AUC_diff.pdf", width=3, height=3)
print(P)
dev.off()

# 15. meta-analysis EO vs LO iCAF ---------------------------------------

Allmarkers <- readRDS("./5.5.Fibr_harmony/Markers/Nomix_celltype_markers.Rds")

APOD <- Allmarkers[Allmarkers$cluster == "CAF_APOD",]
APOD_markers <- APOD[which(APOD$p_val_adj < 0.01 & APOD$pct.1 > 0.3 & APOD$avg_log2FC > 1),]

CCL2 <- Allmarkers[Allmarkers$cluster == "CAF_CCL2",]
CCL2_markers <- CCL2[which(CCL2$p_val_adj < 0.01 & CCL2$pct.1 > 0.3 & CCL2$avg_log2FC > 1),]

gs <- list(iCAF_APOD = APOD_markers$gene,
           iCAF_CCL2 = CCL2_markers$gene)

## TCGA ----------------------------

expMtx <- read.table("./0.TCGA/mRNA_Matrix_TPM.txt", row.names = 1, header = T, sep = "\t")
expMtx <- as.matrix(expMtx)
cells_rankings <- AUCell_buildRankings(expMtx, nCores=10, plotStats=TRUE)
cells_AUC_TCGA <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC_TCGA, "./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_TCGA.Rds")

cells_AUC_TCGA <- readRDS("./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_TCGA.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/Clinical.txt", row.names = 1, header = T, sep = "\t")
age <- as.data.frame(clinical[,1])
colnames(age) <- "age"
rownames(age) <- gsub("_", ".", rownames(clinical))
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]

mtx_TCGA_logit <- mtx
mtx_TCGA_logit$Type <- ifelse(mtx_TCGA_logit$Type == "EOPC", 1, 0)

## GSE141551 ---------------------------------

exp <- read.table("./0.TCGA/GSE141551/MergeExpro_contrib1-GPL14951.txt", row.names = 1, header = T, sep = "\t")
exp <- as.matrix(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=20, plotStats=TRUE)
cells_AUC_GSE141551 <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE141551, "./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE141551.Rds")

cells_AUC_GSE141551 <- readRDS("./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE141551.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE141551@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE141551/SampleInfo_contrib1-GPL14951.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(agegroup = clinical[,"X1.AGE_GROUP"])
rownames(age) <- gsub("_", ".", rownames(clinical))
age$ageup <- as.numeric(substr(age$agegroup,4,5))
age$agelow <- as.numeric(substr(age$agegroup,1,2))
age$avg <- apply(age[,c(2,3)], 1, median)
age$Type[age$ageup < 55] <- "EOPC"
age$Type[age$agelow >= 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-c(1:4)]

mtx_GSE141551_logit <- mtx
mtx_GSE141551_logit$Type <- ifelse(mtx_GSE141551_logit$Type == "EOPC", 1, 0)

## GSE157547 ----------------------------------------

exp <- read.table("./0.TCGA/GSE157547/MergeExpro_contrib1-GPL5175.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
exp <- 2^exp
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE157547 <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE157547, "./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE157547.Rds")

cells_AUC_GSE157547 <- readRDS("./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE157547.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE157547@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE157547/SampleInfo_contrib1-GPL5175.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]

mtx_GSE157547_logit <- mtx
mtx_GSE157547_logit$Type <- ifelse(mtx_GSE157547_logit$Type == "EOPC", 1, 0)

## GSE153352 --------------------------------------------------

exp <- read.table("./0.TCGA/GSE153352/MergeExpro_contrib1-GPL5175.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
exp <- 2^exp
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE153352 <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE153352, "./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE153352.Rds")

cells_AUC_GSE153352 <- readRDS("./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE153352.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE153352@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE153352/SampleInfo_contrib1-GPL5175.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]

mtx_GSE153352_logit <- mtx
mtx_GSE153352_logit$Type <- ifelse(mtx_GSE153352_logit$Type == "EOPC", 1, 0)

## GSE88808 -----------------------------------------

exp <- read.table("./0.TCGA/GSE88808/MergeExpro_contrib1-GPL22571.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE88808 <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE88808, "./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE88808.Rds")

cells_AUC_GSE88808 <- readRDS("./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE88808.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE88808@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE88808/SampleInfo_contrib1-GPL22571.txt", row.names = 1, header = T, sep = "\t")
clinical <- clinical[clinical$"X1.TISSUE" == "tumor",]
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]

mtx_GSE88808_logit <- mtx
mtx_GSE88808_logit$Type <- ifelse(mtx_GSE88808_logit$Type == "EOPC", 1, 0)

## GSE62116 ----------------------------------------

exp <- read.table("./0.TCGA/GSE62116/MergeExpro_contrib1-GPL5188.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
exp <- 2^exp
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE62116 <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE62116, "./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE62116.Rds")

cells_AUC_GSE62116 <- readRDS("./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE62116.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE62116@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE62116/SampleInfo_contrib1-GPL5188.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]

mtx_GSE62116_logit <- mtx
mtx_GSE62116_logit$Type <- ifelse(mtx_GSE62116_logit$Type == "EOPC", 1, 0)

## GSE21034 ----------------------------------------

exp <- read.table("./0.TCGA/GSE21034/MergeExpro_contrib1-GPL5188.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE21034 <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE21034, "./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE21034.Rds")

cells_AUC_GSE21034 <- readRDS("./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE21034.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE21034@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE21034/GSE21034_clinical.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"DxAge"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]

mtx_GSE21034_logit <- mtx
mtx_GSE21034_logit$Type <- ifelse(mtx_GSE21034_logit$Type == "EOPC", 1, 0)

## GSE183019 ----------------------------------------

exp <- read.table("./0.TCGA/GSE183019/GSE183019_mRNA_processed_TPM.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE183019 <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE183019, "./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE183019.Rds")

cells_AUC_GSE183019 <- readRDS("./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE183019.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE183019@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE183019/clinical.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"Age"])
rownames(age) <- gsub("-", ".", clinical$Sample)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]

mtx_GSE183019_logit <- mtx
mtx_GSE183019_logit$Type <- ifelse(mtx_GSE183019_logit$Type == "EOPC", 1, 0)

## GSE201284 ----------------------------------------

exp <- read.table("./0.TCGA/GSE201284/GSE201284_mRNA_processed_TPM.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE201284 <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE201284, "./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE201284.Rds")

cells_AUC_GSE201284 <- readRDS("./5.5.Fibr_harmony/Bulk_CAF_score/iCAF_AUC_GSE201284.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE201284@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE201284/clinical.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"Age"])
rownames(age) <- gsub("-", ".", clinical$Sample)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]

mtx_GSE201284_logit <- mtx
mtx_GSE201284_logit$Type <- ifelse(mtx_GSE201284_logit$Type == "EOPC", 1, 0)

## logistic regression and meta-analysis ---------------------------------------

mtx_GSE141551_logit_1 <- rbind(mtx_TCGA_logit[1:300,], mtx_GSE183019_logit[1:203,])

MP1_logit_list <- list(TCGA = mtx_TCGA_logit, GSE183019 = mtx_GSE183019_logit, GSE201284 = mtx_GSE201284_logit, 
                       GSE141551 = mtx_GSE141551_logit_1, GSE153352=mtx_GSE153352_logit, GSE157547=mtx_GSE157547_logit, 
                       GSE88808 = mtx_GSE88808_logit, GSE62116 = mtx_GSE62116_logit, GSE21034 = mtx_GSE21034_logit)

result <- c()
for (i in 1:length(MP1_logit_list)) {
  data <- MP1_logit_list[[i]]
  data$iCAF_APOD <- scale(data$iCAF_APOD)
  logit <- glm(Type ~ iCAF_APOD, data = data, family = binomial(link = "logit"))
  OR <- exp(summary(logit)$coefficients)[2,1]
  LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2]
  UCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2]
  P <- summary(logit)$coefficients[2,4]
  TE <- summary(logit)$coefficients[2,1]
  seTE <- summary(logit)$coefficients[2,2]
  
  res <- data.frame(Dataset = names(MP1_logit_list)[i], OR = OR, LCI = LCI, 
                    UCI = UCI, P = P, TE = TE, seTE = seTE)
  
  result <- rbind(result, res)
}

library(meta)

str(result)

meta <- metagen(TE = TE, seTE = seTE, studlab = Dataset, data = result, sm = "OR", fixed = F)

pdf("./5.5.Fibr_harmony/Bulk_CAF_score/bulk_APOD_meta.pdf", width = 20, height = 7)
forest(meta, text.random = "", text.w.random = "", print.subgroup.name = F, lwd = 3, hetlab = "",
       test.subgroup = F, col.study = "#556B2F", leftcols = "studlab", leftlabs = "Dataset",
       print.tau2 = F)
dev.off()


mtx_GSE153352_logit_1 <- mtx_GSE153352_logit
mtx_GSE153352_logit_1$Type[mtx_GSE153352_logit_1$Type == 1] <- 2
mtx_GSE153352_logit_1$Type[mtx_GSE153352_logit_1$Type == 0] <- 1
mtx_GSE153352_logit_1$Type[mtx_GSE153352_logit_1$Type == 2] <- 0
mtx_GSE157547_logit_1 <- mtx_GSE157547_logit
mtx_GSE157547_logit_1$Type[mtx_GSE157547_logit_1$Type == 1] <- 2
mtx_GSE157547_logit_1$Type[mtx_GSE157547_logit_1$Type == 0] <- 1
mtx_GSE157547_logit_1$Type[mtx_GSE157547_logit_1$Type == 2] <- 0

MP1_logit_list <- list(TCGA = mtx_TCGA_logit, GSE183019 = mtx_GSE183019_logit, GSE201284 = mtx_GSE201284_logit, 
                       GSE141551 = mtx_GSE141551_logit, GSE153352=mtx_GSE153352_logit_1, GSE157547=mtx_GSE157547_logit_1, 
                       GSE88808 = mtx_GSE88808_logit, GSE62116 = mtx_GSE62116_logit, GSE21034 = mtx_GSE21034_logit)

result <- c()
for (i in 1:length(MP1_logit_list)) {
  data <- MP1_logit_list[[i]]
  data$iCAF_CCL2 <- scale(data$iCAF_CCL2)
  logit <- glm(Type ~ iCAF_CCL2, data = data, family = binomial(link = "logit"))
  OR <- exp(summary(logit)$coefficients)[2,1]
  LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2]
  UCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2]
  P <- summary(logit)$coefficients[2,4]
  TE <- summary(logit)$coefficients[2,1]
  seTE <- summary(logit)$coefficients[2,2]
  
  res <- data.frame(Dataset = names(MP1_logit_list)[i], OR = OR, LCI = LCI, 
                    UCI = UCI, P = P, TE = TE, seTE = seTE)
  
  result <- rbind(result, res)
}

library(meta)

str(result)

meta <- metagen(TE = TE, seTE = seTE, studlab = Dataset, data = result, sm = "OR", fixed = F)

pdf("./5.5.Fibr_harmony/Bulk_CAF_score/bulk_CCL2_meta.pdf", width = 20, height = 7)
forest(meta, text.random = "", text.w.random = "", print.subgroup.name = F, lwd = 3, hetlab = "",
       test.subgroup = F, col.study = "#556B2F", leftcols = "studlab", leftlabs = "Dataset",
       print.tau2 = F)
dev.off()
