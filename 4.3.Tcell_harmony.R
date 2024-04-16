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

# subset analysis ---------------------------------------------------------------------

PCSC <- readRDS("./4.3.Doublet/PCSC_inferred_doublet.Rds")

# Tcell --------------------------------------------------------------------

folder <- "5.2.Tcell_harmony"
if(!dir.exists(folder)){
  dir.create(folder)
}

# 1. subset ----------------------------------------------------------------

Tcell <- subset(PCSC, celltype == "T-cell")

# 2.harmony ---------------------------------------------------------

Tcell <- NormalizeData(Tcell)

Tcell <- FindVariableFeatures(Tcell,
                              selection.method = "vst",
                              nfeatures = 2000,
                              verbose = FALSE) %>%
  ScaleData(vars.to.regress = c("percent.mt"))
Tcell <- RunPCA(Tcell)

pdf("./5.2.Tcell_harmony/harmony_convergence.pdf")
Tcell <- RunHarmony(Tcell, "orig.ident", plot_convergence = T)
dev.off()

Tcell <- RunUMAP(Tcell, reduction = "harmony", dims = 1:20)

# 3.clustering -------------------------------------------------------

Tcell <- FindNeighbors(Tcell, reduction = "harmony", 
                       dims = 1:20)
for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1)) {
  Tcell <- FindClusters(Tcell,  
                        resolution = res)}

apply(Tcell@meta.data[, grep("RNA_snn_res", colnames(Tcell@meta.data))],2,table)

library(clustree)

p2_tree=clustree(Tcell@meta.data, prefix = "RNA_snn_res.") +
  scale_color_npg()
ggsave(plot=p2_tree, filename="./5.2.Tcell_harmony/Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(Tcell) <- "RNA_snn_res.0.8"

pdf("./5.2.Tcell_harmony/CellCluster-UMAPPlot_res0.8.pdf",width = 7,height = 7)
DimPlot(Tcell, label.size = 5, reduction = "umap", pt.size = 0.1, label = T, cols = color) + NoLegend()
dev.off()
png("./5.2.Tcell_harmony/CellCluster-UMAPPlot_res0.8.png", width = 360, height = 360)
DimPlot(Tcell, label.size = 5, reduction = "umap", pt.size = 0.1, label = T, cols = color) + NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/CellCluster-UMAPPlot_SampleType_res0.8.pdf",width = 7,height = 7)
DimPlot(Tcell, reduction = "umap", group.by = "SampleType", pt.size = 0.1, label = T)
dev.off()

pdf("./5.2.Tcell_harmony/CellCluster-UMAPPlot_Sample_res0.8.pdf",width = 7,height = 7)
DimPlot(Tcell, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F, shuffle = T)
dev.off()

pdf("./5.2.Tcell_harmony/CellCluster-UMAPPlot_SampleTypeSplit_res0.8.pdf",width = 14,height = 7)
DimPlot(Tcell, reduction = "umap", split.by = "SampleType", pt.size = 0.1, label = T) + NoLegend()
dev.off()

#2.Cell cycle scoring

folder <- "./5.2.Tcell_harmony/cellcycle"
if(!dir.exists(folder)){
  dir.create(folder)
}

str(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Tcell <- CellCycleScoring(Tcell,
                          g2m.features = g2m.genes,
                          s.features = s.genes)
View(Tcell@meta.data)

pdf("./5.2.Tcell_harmony/cellcycle/cellcycle.pdf",width = 12, height = 4)
DimPlot(Tcell,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()

pdf("./5.2.Tcell_harmony/cellcycle/cellcycle_SampleSplit.pdf",width = 12, height = 4)
DimPlot(Tcell, 
        reduction = "pca", 
        group.by = "orig.ident",
        split.by = "Phase",
        shuffle = T)
dev.off()

pdf("./5.2.Tcell_harmony/cellcycle/CellCluster-UMAPPlot_cellcycle_res0.8.pdf",width = 7,height = 7)
DimPlot(Tcell, reduction = "umap", group.by = "Phase", pt.size = 0.1, label = T)
dev.off()

pdf("./5.2.Tcell_harmony/cellcycle/CellCluster-UMAPPlot_SampleTypeSplit_cellcycle_res0.8.pdf",width = 13,height = 7)
DimPlot(Tcell, reduction = "umap", group.by = "Phase", split.by = "SampleType", pt.size = 0.1, label = T)
dev.off()

table(Tcell@meta.data$orig.ident,Tcell@meta.data$RNA_snn_res.0.8)
write.table(table(Tcell@meta.data$SampleType,Tcell@meta.data$RNA_snn_res.0.8),
            "./5.2.Tcell_harmony/Sampletype/TypeSplit_cellcluster_proportion.txt",sep="\t",quote = F)

saveRDS(Tcell, "./5.2.Tcell_harmony/Tcell_harmonied.Rds")

# 5.find markers cluster------------------------------------

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 1000 * 1024^2)

Idents(Tcell) <- "RNA_snn_res.0.8"
Allmarkers_cl <- FindAllMarkers(Tcell,
                                min.pct = 0.1, 
                                logfc.threshold = 0,
                                verbose = FALSE)
saveRDS(Allmarkers_cl, "./5.2.Tcell_harmony/Markers/cluster_markers.Rds")

topmarkers <- Allmarkers_cl %>% group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)


Tcell_scale <- Tcell
Tcell_scale <- ScaleData(Tcell_scale, features = rownames(Tcell_scale))

#Heatmap
Tcell_reorder <- Tcell_scale
Tcell_reorder$RNA_snn_res.0.8 <- factor(Tcell_reorder$RNA_snn_res.0.8,
                                           levels = c("10","5","9","6","12","3","4","0","2","1","7","13","8","11"))

pdf("./5.2.Tcell_harmony/Markers/Heatmap_clusters.pdf", width = 10, height = 10)
DoHeatmap(Tcell_reorder, features = topmarkers$gene, group.by = "RNA_snn_res.0.8", group.colors=color)
dev.off()

#Heatmap by cluster

Cluster <- data.frame(ID = colnames(Tcell_reorder), Cluster = Tcell_reorder$RNA_snn_res.0.8)
Cluster <- Cluster[order(Cluster$Cluster),]
annotation <- data.frame(Cluster = Cluster$Cluster)
rownames(annotation) <- Cluster$ID

top <- data.frame(gene = topmarkers$gene, Cluster = topmarkers$cluster)
top$Cluster <- factor(top$Cluster, 
                      levels = c("10","5","9","6","12","3","4","0","2","1","7","13","8","11"))
top <- top[order(top$Cluster),]

data <- as.matrix(Tcell_reorder@assays$RNA@data)
genes <- top$gene
data <- data[rownames(data) %in% genes,]
data <- data[match(genes, rownames(data)),]
data <- data[, Cluster$ID]

Cluster_colors <- color[c(1:14)]
names(Cluster_colors) <- as.character(levels(Tcell_reorder$RNA_snn_res.0.8))
ann_colors <- list(Cluster = Cluster_colors)

pdf("./5.2.Tcell_harmony/Markers/celltype_heatmap.pdf", width = 5, height = 8)
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

# 6. Annotation ---------------------------------------------------------------

Tcell <- readRDS("./5.2.Tcell_harmony/Tcell_harmonied.Rds")
DefaultAssay(Tcell) <- "RNA"

# 6.1.公认marker注释 ---------------------------------------------

pdf("./5.2.Tcell_harmony/markerPlot/CD4.pdf", width = 4, height = 4)
FeaturePlot(object = Tcell, features = c("CD4"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2)
dev.off()
png("./5.2.Tcell_harmony/markerPlot/CD4.png", width = 240, height = 240)
FeaturePlot(object = Tcell, features = c("CD4"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2)
dev.off()

pdf("./5.2.Tcell_harmony/markerPlot/CD8.pdf", width = 4, height = 4)
FeaturePlot(object = Tcell, features = c("CD8A"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2)
dev.off()
png("./5.2.Tcell_harmony/markerPlot/CD8.png", width = 240, height = 240)
FeaturePlot(object = Tcell, features = c("CD8A"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2)
dev.off()

pdf("./5.2.Tcell_harmony/markerPlot/NKT.pdf", width = 8.5, height = 8.5)
FeaturePlot(object = Tcell, features = c("KLRB1", "FCGR3A", "KLRK1", "CD160"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 2)
dev.off()

pdf("./5.2.Tcell_harmony/markerPlot/Naive.pdf", width = 13, height = 8.5)
FeaturePlot(object = Tcell, features = c("CCR7", "SELL", "TCF7", "LEF1",
                                         "IL7R"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 3)
dev.off()

pdf("./5.2.Tcell_harmony/markerPlot/Effector.pdf", width = 13, height = 8.5)
FeaturePlot(object = Tcell, features = c("GZMA", "LAG3", "FASLG", "CD44",
                                         "FAS"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 3)
dev.off()

pdf("./5.2.Tcell_harmony/markerPlot/Memory.pdf", width = 8.5, height = 8.5)
FeaturePlot(object = Tcell, features = c("CD44", "IFNG", "S100A4", "GPR183"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 2)
dev.off()

pdf("./5.2.Tcell_harmony/markerPlot/Treg.pdf", width = 8.5, height = 8.5)
FeaturePlot(object = Tcell, features = c("FOXP3", "CTLA4", "IL2RA", "ENTPD1"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 2)
dev.off()

pdf("./5.2.Tcell_harmony/markerPlot/cytotoxic.pdf", width = 13, height = 13)
FeaturePlot(object = Tcell, features = c("RPF1", "GZMB", "NKG7", "CCL4",
                                         "CST7", "GZMA", "IFNG", "CCL3",
                                         "GZMK"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 3)
dev.off()

pdf("./5.2.Tcell_harmony/markerPlot/dysfunction.pdf", width = 17.5, height = 8.5)
FeaturePlot(object = Tcell, features = c("TIGIT", "PDCD1", "CXCL13", "LAG3",
                                         "HAVCR2", "CD38", "ENTPD1", "CTLA4"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 4)
dev.off()

pdf("./5.2.Tcell_harmony/markerPlot/Naive_Effector_or_memory.pdf", width = 17.5, height = 8.5)
FeaturePlot(object = Tcell, features = c("TIGIT", "PDCD1", "CXCL13", "LAG3",
                                         "HAVCR2", "CD38", "ENTPD1", "CTLA4"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 4)
dev.off()

## VlnHeatPlot ------------------------------------------------------

Features <- c("CD3D", "CD3E", "CD4", "CD8A", "CCR7", "TCF7", "CD27", "IL7R", "GZMA", "GZMB", 
              "CCL5", "FCGR3A", "FOXP3", "CTLA4", "TRAV1-2", "TRDC", "TRDV2", "CLIC3", "NKG7", "NCAM1",
              "KLK3")

Cluster <- data.frame(Cluster = Idents(Tcell))
data <- as.matrix(Tcell@assays$RNA@data)
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
pdf("./5.2.Tcell_harmony/markerPlot/VlnHeatmapPlot_markers.pdf")
ggplot(mean, aes(x = gene, y = Cluster)) +
  geom_jjviomap(aes(val = Exp, fill = Cluster), width = 1) +
  coord_fixed() +
  theme_base() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

# 6.2.Monocel -- CD4 cells ------------------------------------------

folder <- "5.2.Tcell_harmony/Monocle/CD4_Tcell"
if(!dir.exists(folder)){
  dir.create(folder)
}

CD4_Tcell <- subset(Tcell, subset = RNA_snn_res.0.8 %in% c(3,5,6,9,12))

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/CD4_T_UMAP.pdf",width = 4.5,height = 4)
DimPlot(CD4_Tcell, reduction = "umap", group.by = "RNA_snn_res.0.8", pt.size = 0.1, label = T) +
  scale_color_npg()
dev.off()

#构造表达及注释数据
exp.matrix <- as(as.matrix(CD4_Tcell@assays$RNA@data), "sparseMatrix")
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(exp.matrix)
exp_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- CD4_Tcell@meta.data
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
saveRDS(ordering_genes, "./5.2.Tcell_harmony/Monocle/CD4_Tcell/ordering_genes.Rds")
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
#plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

#修改monocle对象中的列名示例
names(pData(exp.monocle))[names(pData(exp.monocle))=="RNA_snn_res.0.8"]="Cluster"
saveRDS(exp.monocle, "./5.2.Tcell_harmony/Monocle/CD4_Tcell/exp.monocle.Rds")

#将不同分组情况的拟时序轨迹图画到一起

plot1<-plot_cell_trajectory(exp.monocle, color_by = "Cluster",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color[c(4,6,7,10,13)])
plot2<-plot_cell_trajectory(exp.monocle, color_by = "State",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
plot3<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size = 0.2, shuffle = T) + 
  scale_color_gradient(low = "#101D44", high = "#9FB2ED")
plot4<-plot_cell_trajectory(exp.monocle, color_by = "SampleType",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color[c(2,1)])
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/trajectory_plot.pdf",width = 12,height = 10)
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = NULL,ncol=2)
dev.off()

#split by clusters

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

#cluster
for (i in c(3,5,6,9,12)) {
  a <- rownames(CD4_Tcell@meta.data)[ CD4_Tcell@meta.data$RNA_snn_res.0.8 == i ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$Cluster == i ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(paste0("Cluster ", i))
  assign(paste0("P", i), plot)
}

Plist <- list(P3,P5,P6,P9,P12)
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/trajectory_plot_clusterSplit.pdf",width = 9,height = 5)
CombinePlots(plots = Plist,legend = NULL,ncol=3)
dev.off()

#split by sampletype

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

#cluster
for (i in c("EOPC", "LOPC")) {
  a <- rownames(CD4_Tcell@meta.data)[ CD4_Tcell@meta.data$SampleType == i ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, a] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(i)
  assign(paste0("P_", i), plot)
}

Plist <- list(P_EOPC,P_LOPC)
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/trajectory_plot_SampleTypeSplit.pdf",width = 6, height = 2.5)
CombinePlots(plots = Plist,legend = NULL,ncol=2)
dev.off()

#Staes_SampleType_Roe--------------------------------------

tab <- table(exp.monocle$SampleType, exp.monocle$State)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/State_Roe_SampleType.pdf", width = 4, height = 1)
pheatmap(Roe,
         cluster_cols = T,
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

#features plots

gene <- c("TIGIT", "FOXP3", "CTLA4", "LAG3", "HAVCR2",
          "CCR7", "SELL", "TCF7", "IL7R", "LEF1")

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/feature_naive_dysfunction_plot.pdf", width = 5, height = 5)
plot_genes_in_pseudotime(exp.monocle[gene,], color_by = "Pseudotime", ncol =2)
dev.off()

#Stimulatory and Inhibitory phenotype score (irGSEA) ------------------------------

library("irGSEA")

Stim_Inhi <- list()
Stim_Inhi$Stimulatory <- c("BTN3A1", "BTN3A2", "CCL5", "CD27", "CD28", "CD40", 
                           "CD40LG", "CD70", "CD80", "CX3CL1", "CXCL10", "CXCL9", 
                           "ENTPD1", "HMGB1", "ICAM1", "ICOS", "ICOSLG", "IFNA1", 
                           "IFNA2", "IFNG", "IL1A", "IL1B", "IL2", "IL2RA", "ITGB2", 
                           "PRF1", "SELP", "TLR4", "TNF", "TNFRSF14", "TNFRSF18", 
                           "TNFRSF4", "TNFRSF9", "TNFSF4", "TNFSF9", "C10orf54", 
                           "CD276", "CD48", "CD86", "HHLA2", "IL6", "IL6R", "KLRC1", 
                           "KLRK1", "LTA", "MICB", "NT5E", "PVR", "RAET1E", "TMEM173", 
                           "TMIGD2", "TNFRSF13B", "TNFRSF13C", "TNFRSF17", "TNFRSF25", 
                           "TNFRSF8", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", 
                           "TNFSF18", "ULBP1")
Stim_Inhi$Inhibitory <- c("ADORA2A", "ARG1", "BTLA", "CD274", "CD276", "CTLA4", 
                          "EDNRB", "HAVCR2", "IDO1", "IL10", "IL13", "IL4", "KIR2DL1", 
                          "KIR2DL2", "KIR2DL3", "LAG3", "PDCD1", "SLAMF7", "TGFB1", 
                          "TIGIT", "VEGFA", "VEGFB", "C10orf54", "VTCN1", "CD160", 
                          "CD244", "CD96", "CSF1R", "IL10RB", "KDR", "LGALS9", 
                          "PDCD1LG2", "PVRL2", "TGFBR1", "CD24", "CD47", "NKG2D", 
                          "B7H3", "TIM3", "GALECTIN9", "OX40", "OX40L", "GITR", 
                          "GITRL", "B7H4")


# calculate score
CD4_Tcell_irGSEA <- irGSEA.score(object = CD4_Tcell, assay = "RNA", slot = "data", 
                                 seeds = 123, ncores = 20, min.cells = 3, 
                                 min.feature = 0, custom = T, geneset = Stim_Inhi, 
                                 msigdb = F, method = c("AUCell", "UCell", "singscore", "ssgsea"),
                                 aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                 kcdf = 'Gaussian')

AUCell <- as.data.frame(t(as.data.frame(CD4_Tcell_irGSEA@assays$AUCell@data)))
AUCell$SampleType <- CD4_Tcell_irGSEA$SampleType

UCell <- as.data.frame(t(as.data.frame(CD4_Tcell_irGSEA@assays$UCell@data)))
UCell$SampleType <- CD4_Tcell_irGSEA$SampleType

ssgsea <- as.data.frame(t(as.data.frame(CD4_Tcell_irGSEA@assays$ssgsea@data)))
ssgsea$SampleType <- CD4_Tcell_irGSEA$SampleType

singscore <- as.data.frame(t(as.data.frame(CD4_Tcell_irGSEA@assays$singscore@data)))
singscore$SampleType <- CD4_Tcell_irGSEA$SampleType

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Scatter_Stim_Inhi_AUCell.pdf")
pal_for_sampletype <- c(EOPC = "#E64B35", LOPC = "#3C5488")
names(pal_for_sampletype) <- c("EOPC", "LOPC")
ggplot(data = AUCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("AUCell")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Scatter_Stim_Inhi_UCell_nolegend.pdf")
pal_for_sampletype <- c(EOPC = "#E64B35", LOPC = "#3C5488")
names(pal_for_sampletype) <- c("EOPC", "LOPC")
ggplot(data = UCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("UCell") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Scatter_Stim_Inhi_ssgsea_nolegend.pdf")
pal_for_sampletype <- c(EOPC = "#E64B35", LOPC = "#3C5488")
names(pal_for_sampletype) <- c("EOPC", "LOPC")
ggplot(data = ssgsea, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("ssGSEA") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Scatter_Stim_Inhi_singscore_nolegend.pdf")
pal_for_sampletype <- c(EOPC = "#E64B35", LOPC = "#3C5488")
names(pal_for_sampletype) <- c("EOPC", "LOPC")
ggplot(data = singscore, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("singscore") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Scatter_Stim_Inhi_AUCell_nolegend.pdf")
pal_for_sampletype <- c(EOPC = "#E64B35", LOPC = "#3C5488")
names(pal_for_sampletype) <- c("EOPC", "LOPC")
ggplot(data = AUCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("AUCell") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

#Violin plot
library(reshape2)
AUCell_plot = melt(AUCell)
colnames(AUCell_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Violin_Stim_Inhi_AUCell.pdf")
ggviolin(AUCell_plot, x="Status", y="value", fill = "SampleType", 
         ylab="AUCell Score",
         xlab="",
         palette = c("#E64B35", "#3C5488"),
         width=0.6, add = c("boxplot"), 
         error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()

UCell_plot = melt(UCell)
colnames(UCell_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Violin_Stim_Inhi_UCell.pdf")
ggviolin(UCell_plot, x="Status", y="value", fill = "SampleType", 
         ylab="UCell Score",
         xlab="",
         palette = c("#E64B35", "#3C5488"),
         width=0.6, add = c("boxplot"), 
         error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()

ssgsea_plot = melt(ssgsea)
colnames(ssgsea_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Violin_Stim_Inhi_ssgsea.pdf")
ggviolin(AUCell_plot, x="Status", y="value", fill = "SampleType", 
         ylab="ssGSEA Score",
         xlab="",
         palette = c("#E64B35", "#3C5488"),
         width=0.6, add = c("boxplot"), 
         error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()

singscore_plot = melt(singscore)
colnames(singscore_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Violin_Stim_Inhi_singscore.pdf")
ggviolin(AUCell_plot, x="Status", y="value", fill = "SampleType", 
         ylab="singscore",
         xlab="",
         palette = c("#E64B35", "#3C5488"),
         width=0.6, add = c("boxplot"), 
         error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()

## Differences between sampletype ----------------------------------------

fre <- table(CD4_Tcell$RNA_snn_res.0.8, CD4_Tcell$orig.ident)
sum <- table(CD4_Tcell$orig.ident)

pct <- c()
for (i in 1:nrow(fre)) {
  pct_tmp <- fre[i,]/sum
  pct <- rbind(pct, pct_tmp)
}
rownames(pct) <- paste0("c",c(0:13))

pct <- pct[c("c3", "c5", "c6", "c9", "c12"),]

pct <- t(pct)
pct <- as.data.frame(pct)

pct$type <- substr(rownames(pct), 1, 4)

library(ggpubr)

#把数据转换成gglpot2输入文件
data <- melt(pct, id.vars = c("type"))
colnames(data) <- c("Type","Group","Pct")

#绘制boxplot
p=ggviolin(data, x="Group", y="Pct", fill = "Type", 
           ylab="Pct. of each cluster in CD4+ T-cells",
           xlab="",
           palette = c("#3C5488","#E64B35"),
           width=0.6, add = c("dotplot"), 
           error.plot = "crossbar", add.params = list(size=0.5)) 
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.05, 1), symbols = c("P < 0.05", " ")),
                        label = "p.signif")

#输出图片
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/CD4_difference.pdf", width=5, height=4)
print(p1)
dev.off()

#FOXP3 & CTLA4
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/DotPlot_FOXP3_CTLA4.pdf", width=5, height=4)
DotPlot(CD4_Tcell, features=c("FOXP3", "CTLA4"), group.by="SampleType")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/VlnPlot_FOXP3_CTLA4.pdf", width=3, height=4)
VlnPlot(CD4_Tcell, features=c("FOXP3", "CTLA4"), group.by="SampleType")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/FOXP3.pdf", width=3, height=4)
VlnPlot(CD4_Tcell, features=c("FOXP3", "CTLA4"))
dev.off()

#FOXP3_proportion

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/FOXP3_positive.pdf")
CD4_Tcell$FOXP3_positive <- ifelse(CD4_Tcell@assays$RNA@data["FOXP3",] > 0,"Positive","Negative")
DimPlot(CD4_Tcell, split.by="SampleType", group.by="FOXP3_positive", pt.size=0.5, cols=c("grey", "#E64B35"))
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/FOXP3_positive_bar.pdf", width=3, height=4)
FOXP3_pos <- table(CD4_Tcell$SampleType, CD4_Tcell$FOXP3_positive)
sum <- apply(FOXP3_pos, 1, sum)
FOXP3_pos_pct <- data.frame()
for (n in 1:2) {
  FOXP3_pos_pct = rbind( FOXP3_pos_pct, FOXP3_pos[n,] / sum[n] )  
} 
colnames(FOXP3_pos_pct) <- colnames(FOXP3_pos)
FOXP3_pos_pct$group <- rownames(FOXP3_pos)

library(reshape2)
FOXP3_pos_plot = melt(FOXP3_pos_pct)
colnames(FOXP3_pos_plot) = c('group','positive','percent')

ggplot( FOXP3_pos_plot, aes( x = group, weight = percent, fill = positive))+
  geom_bar( position = "stack") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
dev.off()

#CTLA4
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/CTLA4_positive.pdf")
CD4_Tcell$CTLA4_positive <- ifelse(CD4_Tcell@assays$RNA@data["CTLA4",] > 0,"Positive","Negative")
DimPlot(CD4_Tcell, split.by="SampleType", group.by="CTLA4_positive", pt.size=0.5, cols=c("grey", "#E64B35"))
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/CTLA4_positive_bar.pdf", width=3, height=4)
CTLA4_pos <- table(CD4_Tcell$SampleType, CD4_Tcell$CTLA4_positive)
sum <- apply(CTLA4_pos, 1, sum)
CTLA4_pos_pct <- data.frame()
for (n in 1:2) {
  CTLA4_pos_pct = rbind( CTLA4_pos_pct, CTLA4_pos[n,] / sum[n] )  
} 
colnames(CTLA4_pos_pct) <- colnames(CTLA4_pos)
CTLA4_pos_pct$group <- rownames(CTLA4_pos)

library(reshape2)
CTLA4_pos_plot = melt(CTLA4_pos_pct)
colnames(CTLA4_pos_plot) = c('group','positive','percent')

ggplot( CTLA4_pos_plot, aes( x = group, weight = percent, fill = positive))+
  geom_bar( position = "stack") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
dev.off()

# 6.3.Monocel -- CD8 cells ------------------------------------------

folder <- "5.2.Tcell_harmony/Monocle/CD8_Tcell"
if(!dir.exists(folder)){
  dir.create(folder)
}

CD8_Tcell <- subset(Tcell, subset = RNA_snn_res.0.8 %in% c(0,1,2,4,7,13))

#构造表达及注释数据
exp.matrix <- as(as.matrix(CD8_Tcell@assays$RNA@data), "sparseMatrix")
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(exp.matrix)
exp_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- CD8_Tcell@meta.data
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
saveRDS(ordering_genes, "./5.2.Tcell_harmony/Monocle/CD8_Tcell/ordering_genes.Rds")
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
#plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

#修改monocle对象中的列名示例
names(pData(exp.monocle))[names(pData(exp.monocle))=="RNA_snn_res.0.8"]="Cluster"
saveRDS(exp.monocle, "./5.2.Tcell_harmony/Monocle/CD8_Tcell/exp.monocle.Rds")

#将不同分组情况的拟时序轨迹图画到一起

plot1<-plot_cell_trajectory(exp.monocle, color_by = "Cluster",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color[c(1,2,3,5,8,14)])
plot2<-plot_cell_trajectory(exp.monocle, color_by = "State",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
plot3<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size = 0.2, shuffle = T) + 
  scale_color_gradient(low = "#101D44", high = "#9FB2ED")
plot4<-plot_cell_trajectory(exp.monocle, color_by = "SampleType",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color[c(2,1)])
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/trajectory_plot.pdf",width = 12,height = 10)
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = NULL,ncol=2)
dev.off()

#split by clusters

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

#cluster
for (i in c(0,1,2,4,7,13)) {
  a <- rownames(CD8_Tcell@meta.data)[ CD8_Tcell@meta.data$RNA_snn_res.0.8 == i ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$Cluster == i ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(paste0("Cluster ", i))
  assign(paste0("P", i), plot)
}

Plist <- list(P3,P5,P6,P9,P12)
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/trajectory_plot_clusterSplit.pdf",width = 9,height = 5)
CombinePlots(plots = Plist,legend = NULL,ncol=3)
dev.off()

#split by sampletype

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

#cluster
for (i in c("EOPC", "LOPC")) {
  a <- rownames(CD8_Tcell@meta.data)[ CD8_Tcell@meta.data$SampleType == i ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$SampleType == i ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(i)
  assign(paste0("P_", i), plot)
}

Plist <- list(P_EOPC,P_LOPC)
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/trajectory_plot_SampleTypeSplit.pdf",width = 6, height = 2.5)
CombinePlots(plots = Plist,legend = NULL,ncol=2)
dev.off()

#Staes_SampleType_Roe--------------------------------------

tab <- table(exp.monocle$SampleType, exp.monocle$State)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/State_Roe_SampleType.pdf", width = 4, height = 2)
bk <- c(seq(0.5,1.5,by=0.01))
pheatmap(Roe,
         cluster_cols = T,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = T,
         show_rownames = T,
         scale="none",  #矫正
         breaks=bk,
         legend_breaks=seq(-0.5,0.5,0.5),
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Ro/e",
         angle_col = "0",
         display_numbers = T,
         legend = F)
dev.off()


#Stimulatory and Inhibitory phenotype score (irGSEA) ------------------------------

library("AUCell")

Stim_Inhi <- list()
Stim_Inhi$Stimulatory <- c("BTN3A1", "BTN3A2", "CCL5", "CD27", "CD28", "CD40", 
                           "CD40LG", "CD70", "CD80", "CX3CL1", "CXCL10", "CXCL9", 
                           "ENTPD1", "HMGB1", "ICAM1", "ICOS", "ICOSLG", "IFNA1", 
                           "IFNA2", "IFNG", "IL1A", "IL1B", "IL2", "IL2RA", "ITGB2", 
                           "PRF1", "SELP", "TLR4", "TNF", "TNFRSF14", "TNFRSF18", 
                           "TNFRSF4", "TNFRSF9", "TNFSF4", "TNFSF9", "C10orf54", 
                           "CD276", "CD48", "CD86", "HHLA2", "IL6", "IL6R", "KLRC1", 
                           "KLRK1", "LTA", "MICB", "NT5E", "PVR", "RAET1E", "TMEM173", 
                           "TMIGD2", "TNFRSF13B", "TNFRSF13C", "TNFRSF17", "TNFRSF25", 
                           "TNFRSF8", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", 
                           "TNFSF18", "ULBP1")
Stim_Inhi$Inhibitory <- c("ADORA2A", "ARG1", "BTLA", "CD274", "CD276", "CTLA4", 
                          "EDNRB", "HAVCR2", "IDO1", "IL10", "IL13", "IL4", "KIR2DL1", 
                          "KIR2DL2", "KIR2DL3", "LAG3", "PDCD1", "SLAMF7", "TGFB1", 
                          "TIGIT", "VEGFA", "VEGFB", "C10orf54", "VTCN1", "CD160", 
                          "CD244", "CD96", "CSF1R", "IL10RB", "KDR", "LGALS9", 
                          "PDCD1LG2", "PVRL2", "TGFBR1", "CD24", "CD47", "NKG2D", 
                          "B7H3", "TIM3", "GALECTIN9", "OX40", "OX40L", "GITR", 
                          "GITRL", "B7H4")


# calculate score
CD8_Tcell_irGSEA <- irGSEA.score(object = CD8_Tcell, assay = "RNA", slot = "data", 
                                 seeds = 123, ncores = 20, min.cells = 3, 
                                 min.feature = 0, custom = T, geneset = Stim_Inhi, 
                                 msigdb = F, method = c("AUCell", "UCell", "singscore", "ssgsea"),
                                 aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                 kcdf = 'Gaussian')

saveRDS(CD8_Tcell_irGSEA, "5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Stimu_Inhib_irGSEA.Rds")

AUCell <- as.data.frame(t(as.data.frame(CD8_Tcell_irGSEA@assays$AUCell@data)))
AUCell$SampleType <- CD8_Tcell_irGSEA$SampleType

UCell <- as.data.frame(t(as.data.frame(CD8_Tcell_irGSEA@assays$UCell@data)))
UCell$SampleType <- CD8_Tcell_irGSEA$SampleType

ssgsea <- as.data.frame(t(as.data.frame(CD8_Tcell_irGSEA@assays$ssgsea@data)))
ssgsea$SampleType <- CD8_Tcell_irGSEA$SampleType

singscore <- as.data.frame(t(as.data.frame(CD8_Tcell_irGSEA@assays$singscore@data)))
singscore$SampleType <- CD8_Tcell_irGSEA$SampleType

#Scatter Plot 

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Scatter_Stim_Inhi_AUCell.pdf")
pal_for_sampletype <- c(EOPC = "#053E7A", LOPC = "#A52A2A")
names(pal_for_sampletype) <- c("EOPC", "LOPC")
ggplot(data = AUCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("AUCell") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Scatter_Stim_Inhi_UCell_nolegend.pdf")
ggplot(data = UCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("UCell") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Scatter_Stim_Inhi_ssgsea_nolegend.pdf")
ggplot(data = ssgsea, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("ssGSEA") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Scatter_Stim_Inhi_singscore_nolegend.pdf")
ggplot(data = singscore, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("singscore") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Scatter_Stim_Inhi_AUCell_nolegend.pdf")
ggplot(data = AUCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("AUCell") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

#Violin Plot 

library(reshape2)
AUCell_plot = melt(AUCell)
colnames(AUCell_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Violin_Stim_Inhi_AUCell.pdf")
ggviolin(AUCell_plot, x="Status", y="value", fill = "SampleType", 
           ylab="AUCell Score",
           xlab="",
           palette = c("#053E7ACC", "#A52A2ACC"),
           width=0.6, add = c("boxplot"), 
           error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()

UCell_plot = melt(UCell)
colnames(UCell_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Violin_Stim_Inhi_UCell.pdf")
ggviolin(UCell_plot, x="Status", y="value", fill = "SampleType", 
         ylab="UCell Score",
         xlab="",
         palette = c("#053E7ACC", "#A52A2ACC"),
         width=0.6, add = c("boxplot"), 
         error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()

ssgsea_plot = melt(ssgsea)
colnames(ssgsea_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Violin_Stim_Inhi_ssgsea.pdf")
ggviolin(AUCell_plot, x="Status", y="value", fill = "SampleType", 
         ylab="ssGSEA Score",
         xlab="",
         palette = c("#053E7ACC", "#A52A2ACC"),
         width=0.6, add = c("boxplot"), 
         error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()

singscore_plot = melt(singscore)
colnames(singscore_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Violin_Stim_Inhi_singscore.pdf")
ggviolin(AUCell_plot, x="Status", y="value", fill = "SampleType", 
         ylab="singscore",
         xlab="",
         palette = c("#053E7ACC", "#A52A2ACC"),
         width=0.6, add = c("boxplot"), 
         error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()

## Differences between sampletype ----------------------------------------

fre <- table(CD8_Tcell$RNA_snn_res.0.8, CD8_Tcell$orig.ident)
sum <- table(CD8_Tcell$orig.ident)

pct <- c()
for (i in 1:nrow(fre)) {
  pct_tmp <- fre[i,]/sum
  pct <- rbind(pct, pct_tmp)
}
rownames(pct) <- paste0("c",c(0:13))

pct <- pct[c("c0", "c1", "c2", "c4", "c7", "c13"),]

pct <- t(pct)
pct <- as.data.frame(pct)

pct$type <- substr(rownames(pct), 1, 4)

library(ggpubr)

#把数据转换成gglpot2输入文件
data <- melt(pct, id.vars = c("type"))
colnames(data) <- c("Type","Group","Pct")

#绘制boxplot
p=ggviolin(data, x="Group", y="Pct", fill = "Type", 
           ylab="Pct. of each cluster in CD8+ T-cells",
           xlab="",
           palette = c("#E64B35", "#3C5488"),
           width=0.6, add = c("dotplot"), 
           error.plot = "crossbar", add.params = list(size=0.5)) 
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.05, 1), symbols = c("P < 0.05", " ")),
                        label = "p.signif")

#输出图片
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/CD8_difference.pdf", width=5, height=4)
print(p1)
dev.off()

# select markers TIGIT LAG3 (dysfunction) IL7R and CCR7 (naive) GZMA IFNG (cytotoxic) KLRG1 B3GAT1 (senescence)---------------------------------

# exp

exp <- CD8_Tcell@assays$RNA@data[c("TIGIT","LAG3", "IL7R", "CCR7", "IFNG", "GZMA", "KLRG1", "B3GAT1"),] %>% as.data.frame() %>% t()
SampleType <- as.data.frame(CD8_Tcell@meta.data$SampleType)
colnames(SampleType) <- "SampleType"
exp <- cbind(exp, SampleType)
exp <- melt(exp)

P = ggplot(exp, aes(x = variable, y = value, fill = SampleType)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  ylab("Norm. Exp.") +
  xlab("") + 
  ggtitle("") +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

pdf("5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Boxplot_selected_genes.pdf", width = 6, height = 3)
print(P)
dev.off()

# proportion

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/LAG3_positive.pdf")
Plist <- list()
genes <- c("TIGIT","LAG3", "IL7R", "CCR7", "IFNG", "GZMA", "KLRG1", "B3GAT1")
for (i in genes) {
  CD8_Tcell$positive <- ifelse(CD8_Tcell@assays$RNA@data[i,] > 0, "Positive", "Negative")
  
  pos <- table(CD8_Tcell$SampleType, CD8_Tcell$positive)
  chi <- chisq.test(pos)
  p <- signif(chi$p.value,3)
  pos <- melt(pos)
  colnames(pos) <- c("SampleType", "Status", "nCells")
  
  barplot <- ggplot(pos, aes(x=SampleType, y=nCells,fill=Status)) +
    geom_bar(stat="identity", position = 'fill') +
    theme_bw() +
    scale_fill_manual(values = c("#BBBBBB", "#556B2FCC")) +
    annotate(geom="text", x=1.5, y=1.05, label=paste0("Chisq P-value = ", p), size = 2) +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ylab(paste0(i, " Pos. Pct.")) +
    xlab("") + 
    ggtitle("") +
    NoLegend()
  
  Plist[[i]] <- barplot
}

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Bar_selected_genes_positive.pdf", width=5, height=12)
grid.arrange(grobs = Plist, ncol = 2)
dev.off()

#dysfunction markers -------------------------------------

exp <- CD8_Tcell@assays$RNA@data[c("TIGIT","LAG3"),] %>% as.data.frame() %>% t()
SampleType <- as.data.frame(CD8_Tcell@meta.data$SampleType)
colnames(SampleType) <- "SampleType"
exp <- cbind(exp, SampleType)
exp <- melt(exp)

P = ggplot(exp, aes(x = variable, y = value, fill = SampleType)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  ylab("AUCell score") +
  xlab("") + 
  ggtitle("") +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

pdf("5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Boxplot_TIGIT_LAG3.pdf", width = 4, height = 3)
print(P)
dev.off()

#LAG3
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/LAG3_positive.pdf")
CD8_Tcell$LAG3_positive <- ifelse(CD8_Tcell@assays$RNA@data["LAG3",] > 0,"Positive","Negative")
DimPlot(CD8_Tcell, split.by="SampleType", group.by="LAG3_positive", pt.size=0.5, cols=c("grey", "#E46B4F"))
dev.off()

LAG3_pos <- table(CD8_Tcell$SampleType, CD8_Tcell$LAG3_positive)
chi <- chisq.test(LAG3_pos)
p <- signif(chi$p.value,3)
pos <- melt(LAG3_pos)
colnames(pos) <- c("SampleType", "Status", "nCells")

barplot <- ggplot(pos, aes(x=SampleType, y=nCells,fill=Status)) +
  geom_bar(stat="identity", position = 'fill') +
  theme_bw() +
  scale_fill_manual(values = c("#BBBBBB", "#556B2FCC")) +
  annotate(geom="text", x=1.5, y=1.05, label=paste0("Chisq P-value = ", p), size = 2) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("LAG3 Pos. Pct.") +
  xlab("") + 
  ggtitle("") +
  NoLegend()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/LAG3_positive_bar.pdf", width=2, height=3)
barplot
dev.off()

#TIGIT

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/TIGIT_positive.pdf")
CD8_Tcell$TIGIT_positive <- ifelse(CD8_Tcell@assays$RNA@data["TIGIT",] > 0,"Positive","Negative")
DimPlot(CD8_Tcell, split.by="SampleType", group.by="TIGIT_positive", pt.size=0.5, cols=c("grey", "#E46B4F"))
dev.off()

TIGIT_pos <- table(CD8_Tcell$SampleType, CD8_Tcell$TIGIT_positive)
chi <- chisq.test(TIGIT_pos)
p <- signif(chi$p.value,3)
pos <- melt(TIGIT_pos)
colnames(pos) <- c("SampleType", "Status", "nCells")

barplot <- ggplot(pos, aes(x=SampleType, y=nCells,fill=Status)) +
  geom_bar(stat="identity", position = 'fill') +
  theme_bw() +
  scale_fill_manual(values = c("#BBBBBB", "#556B2FCC")) +
  annotate(geom="text", x=1.5, y=1.05, label=paste0("Chisq P-value = ", p), size = 2) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("TIGIT Pos. Pct.") +
  xlab("") + 
  ggtitle("") +
  NoLegend()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/TIGIT_positive_bar.pdf", width=2, height=3)
barplot
dev.off()

#Naive markers IL7R and SELL ----------------------------------------

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/DotPlot_IL7R_SELL.pdf")
DotPlot(CD8_Tcell, features=c("IL7R", "SELL"), group.by="SampleType")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/VlnPlot_IL7R_SELL.pdf")
VlnPlot(CD8_Tcell, features=c("IL7R", "SELL"), group.by="SampleType")
dev.off()

#senescence KLRG1 and B3GAT1  -------------------------------------

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/DotPlot_KLRG1_B3GAT1.pdf")
DotPlot(CD8_Tcell, features=c("KLRG1", "B3GAT1"), group.by="SampleType")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/VlnPlot_KLRG1_B3GAT1.pdf")
VlnPlot(CD8_Tcell, features=c("KLRG1", "B3GAT1"), group.by="SampleType")
dev.off()

# 6.4. SingleR注释 ----------------------------------------

library(SingleR)
#celldex::DatabaseImmuneCellExpressionData()
#An external refdata is also acceptable 
load("./DatabaseImmuneCellExpressionData.Rdata")
testdata <- GetAssayData(Tcell, slot="data")
#testdata <- PCSC@assays$SCT@data
clusters <- Tcell@meta.data$RNA_snn_res.0.8

#按cluster进行注释
cellpred <- SingleR(test = testdata, ref = refdata, 
                    labels = refdata$label.main, 
                    clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype <- data.frame(ClusterID=rownames(cellpred), 
                       celltype=cellpred$labels, stringsAsFactors = F)

Tcell@meta.data$singleR_celltype = "Unknown"
for(i in 1:nrow(celltype)){
  Tcell@meta.data[which(Tcell@meta.data$seurat_clusters == celltype$ClusterID[i]),'singleR_celltype'] <- celltype$celltype[i]
}

plot_anno <- DimPlot(Tcell, group.by="singleR_celltype", pt.size=0.1, label=TRUE, label.size=5, reduction='umap', shuffle = T)
ggsave(filename = "./5.2.Tcell_harmony/CellCluster-UMAPPlot_res0.8_SingleR.pdf",
       plot = plot_anno)

# 6.5. scType注释 ---------------------------------------------

library("HGNChelper")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# __secondary----------------------------------------------------------
db_ = "./4.annotation/ScTypeDB_short.xlsx"
tissue = "Tcell" 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix

es.max = sctype_score(scRNAseqData = Tcell[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(Tcell@meta.data$RNA_snn_res.0.8), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(Tcell@meta.data[Tcell@meta.data$RNA_snn_res.0.8==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Tcell@meta.data$RNA_snn_res.0.8==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "Mixed"
print(sctype_scores[,1:3])

Tcell@meta.data$scType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  Tcell@meta.data$scType[Tcell@meta.data$RNA_snn_res.0.8 == j] = as.character(cl_type$type[1])
}

pdf("./5.2.Tcell_harmony/CellCluster-UMAPPlot_res0.8_ScType.pdf", )
DimPlot(Tcell, reduction = "umap", label = TRUE, repel = TRUE, 
        pt.size=0.1, group.by = 'scType', shuffle = T) + NoLegend()
dev.off()

# 6.6. rename -----------------------------------------

celltype=data.frame(ClusterID=0:13,
                    celltype='CD8_Tem')   
celltype[celltype$ClusterID %in% c(0),2]='CD8_Tn'  
celltype[celltype$ClusterID %in% c(3),2]='Tcm'
celltype[celltype$ClusterID %in% c(4),2]='Tcm'
celltype[celltype$ClusterID %in% c(5),2]='Treg'
celltype[celltype$ClusterID %in% c(6),2]='CD4_Tn'
celltype[celltype$ClusterID %in% c(8),2]='NKT'
celltype[celltype$ClusterID %in% c(9),2]='CD4_Tcon'
celltype[celltype$ClusterID %in% c(10),2]='Epi_mix'
celltype[celltype$ClusterID %in% c(11),2]='NK'
celltype[celltype$ClusterID %in% c(12),2]='CD4_Tn'
head(celltype)

Tcell@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  Tcell@meta.data[which(Tcell@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(Tcell@meta.data$celltype)

Idents(Tcell) <- "celltype"
Tcell <- RenameIdents(Tcell,
                      "Tcm" = "Tcm",
                      "CD4_Tn" = "CD4_Tn",
                      "CD4_Tcon" = "CD4_Tcon",
                      "Treg" = "Treg",
                      "CD8_Tn" = "CD8_Tn",
                      "CD8_Tem" = "CD8_Tem",
                      "NK" = "NK",
                      "NKT" = "NKT",
                      "Epi_mix" = "Epi_mix")
Tcell$celltype <- Idents(Tcell)

pdf("./5.2.Tcell_harmony/Tcell_celltype_umap.pdf", width = 5, height = 5)
DimPlot(Tcell, label = TRUE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color) + NoLegend() + 
  ggtitle("")
dev.off()
png("./5.2.Tcell_harmony/Tcell_celltype_umap.png", width = 360, height = 360)
DimPlot(Tcell, label = TRUE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color) + NoLegend() + 
  ggtitle("")
dev.off()

saveRDS(Tcell, "./5.2.Tcell_harmony/Tcell_annotated.Rds")

# Roe

tab <- table(Tcell@meta.data$celltype, Tcell@meta.data$SampleType)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.2.Tcell_harmony/Roe_celltype_SampleType.pdf", width = 2.5, height = 4)
bk <- c(seq(-0.5,0.5,by=0.01))
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

# Allmarker bubble

Features <- c("CD3E", "CD4", "CD8A", "CCR7", "TCF7", "CD27", "IL7R", "GZMA", "GZMB", 
              "CCL5", "FCGR3A", "FOXP3", "CTLA4", "TRAV1-2", "TRDC", "TRDV2", "CLIC3", "NKG7", "NCAM1",
              "KLK3")

pdf("./5.2.Tcell_harmony/AllmarkerBubble.pdf",width = 6, height = 4.8)
jjDotPlot(Tcell,
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

# 7. exclude Epi-mix -----------------------------------------------

Tcell <- readRDS("./5.2.Tcell_harmony/Tcell_annotated.Rds")
Tcell_nomix <- subset(Tcell, celltype != "Epi_mix")

# 7.1. clustering ----------------------------------------------

Tcell_nomix <- NormalizeData(Tcell_nomix)

Tcell_nomix <- FindVariableFeatures(Tcell_nomix,
                              selection.method = "vst",
                              nfeatures = 2000,
                              verbose = FALSE) %>%
  ScaleData(vars.to.regress = c("percent.mt"))
Tcell_nomix <- RunPCA(Tcell_nomix)

pdf("./5.2.Tcell_harmony/Nomix_harmony_convergence.pdf")
Tcell_nomix <- RunHarmony(Tcell_nomix, "orig.ident", plot_convergence = T)
dev.off()

Tcell_nomix <- RunUMAP(Tcell_nomix, reduction = "harmony", dims = 1:20)

# clustering

Tcell_nomix <- FindNeighbors(Tcell_nomix, reduction = "harmony", 
                       dims = 1:20)
for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1)) {
  Tcell_nomix <- FindClusters(Tcell_nomix,  
                        resolution = res)}

apply(Tcell_nomix@meta.data[, grep("RNA_snn_res", colnames(Tcell_nomix@meta.data))],2,table)

library(clustree)

p2_tree=clustree(Tcell_nomix@meta.data, prefix = "RNA_snn_res.") +
  scale_color_manual(values = color)
ggsave(plot=p2_tree, filename="./5.2.Tcell_harmony/Nomix_Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(Tcell_nomix) <- "RNA_snn_res.0.8"

Tcell_nomix <- subset(Tcell_nomix, RNA_snn_res.0.8 != 13)
Tcell_nomix$RNA_snn_res.0.8 <- factor(Tcell_nomix$RNA_snn_res.0.8,
                                      levels = c(0:12))

pdf("./5.2.Tcell_harmony/Nomix_CellCluster-UMAPPlot_res0.8.pdf",width = 4,height = 4)
DimPlot(Tcell_nomix, label.size = 5, reduction = "umap", pt.size = 0.1, label = T, cols = color) + NoLegend()
dev.off()
png("./5.2.Tcell_harmony/Nomix_CellCluster-UMAPPlot_res0.8.png", width = 360, height = 360)
DimPlot(Tcell_nomix, label.size = 5, reduction = "umap", pt.size = 0.1, label = T, cols = color) + NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_CellCluster-UMAPPlot_SampleType_res0.8.pdf",width = 7,height = 7)
DimPlot(Tcell_nomix, reduction = "umap", group.by = "SampleType", pt.size = 0.1, label = T)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_CellCluster-UMAPPlot_Sample_res0.8.pdf",width = 7,height = 7)
DimPlot(Tcell_nomix, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F, shuffle = T)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_CellCluster-UMAPPlot_SampleTypeSplit_res0.8.pdf",width = 14,height = 7)
DimPlot(Tcell_nomix, reduction = "umap", split.by = "SampleType", pt.size = 0.1, label = T) + NoLegend()
dev.off()

saveRDS(Tcell_nomix, "./5.2.Tcell_harmony/Nomix_Tcell_harmonied.Rds")

# 7.2.find markers cluster------------------------------------

Tcell_nomix <- readRDS("./5.2.Tcell_harmony/Nomix_Tcell_annotated.Rds")

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 1000 * 1024^2)

Idents(Tcell_nomix) <- "RNA_snn_res.0.8"
Allmarkers_cl <- FindAllMarkers(Tcell_nomix,
                                min.pct = 0.1, 
                                logfc.threshold = 0,
                                verbose = FALSE)
saveRDS(Allmarkers_cl, "./5.2.Tcell_harmony/Markers/Nomix_cluster_markers.Rds")

topmarkers <- Allmarkers_cl %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

#Heatmap by cluster

Cluster <- data.frame(ID = colnames(Tcell_nomix), Cluster = Idents(Tcell_nomix))
#cluster$Cluster <- factor(cluster$Cluster, levels = c(1,2,..))
Cluster <- Cluster[order(Cluster$Cluster),]
annotation <- data.frame(Cluster = Cluster$Cluster)
rownames(annotation) <- Cluster$ID

top <- data.frame(gene = topmarkers$gene, Cluster = topmarkers$cluster)
top <- top[order(top$Cluster),]

data <- as.matrix(Tcell_nomix@assays$RNA@data)
genes <- top$gene
data <- data[rownames(data) %in% genes,]
data <- data[match(genes, rownames(data)),]
data <- data[, Cluster$ID]

pdf("./5.2.Tcell_harmony/Markers/Nomix_cluster_heatmap.pdf", width = 5, height = 15)
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
         fontsize_row = 7,
         fontsize_col = 12,
         main = "",
         angle_col = "90",
         #border_color = "grey",
         annotation_col = annotation,
         #annotation_colors = c(color, "#776644"),
         cutree_rows = 5)
dev.off()

# 7.3. Annotation ---------------------------------------------------------------

DefaultAssay(Tcell_nomix) <- "RNA"

# 7.3.1.公认marker注释 ---------------------------------------------

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/CD4.pdf", width = 4, height = 4)
FeaturePlot(object = Tcell_nomix, features = c("CD4"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2)
dev.off()
png("./5.2.Tcell_harmony/Nomix_markerPlot/CD4.png", width = 240, height = 240)
FeaturePlot(object = Tcell_nomix, features = c("CD4"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/CD8.pdf", width = 4, height = 4)
FeaturePlot(object = Tcell_nomix, features = c("CD8A"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2)
dev.off()
png("./5.2.Tcell_harmony/Nomix_markerPlot/CD8.png", width = 240, height = 240)
FeaturePlot(object = Tcell_nomix, features = c("CD8A"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/NKT.pdf", width = 8.5, height = 8.5)
FeaturePlot(object = Tcell_nomix, features = c("KLRB1", "FCGR3A", "KLRK1", "CD160"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 2)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/NK.pdf", width = 4, height = 4)
FeaturePlot(object = Tcell_nomix, features = c("NCAM1"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/Naive.pdf", width = 13, height = 8.5)
FeaturePlot(object = Tcell_nomix, features = c("CCR7", "SELL", "TCF7", "LEF1",
                                         "IL7R"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 3)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/Effector.pdf", width = 13, height = 8.5)
FeaturePlot(object = Tcell_nomix, features = c("GZMA", "LAG3", "FASLG", "CD44",
                                         "FAS"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 3)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/Memory.pdf", width = 8.5, height = 8.5)
FeaturePlot(object = Tcell_nomix, features = c("CD44", "IFNG", "S100A4", "GPR183"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 2)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/Treg.pdf", width = 8.5, height = 8.5)
FeaturePlot(object = Tcell_nomix, features = c("FOXP3", "CTLA4", "IL2RA", "ENTPD1"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 2)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/cytotoxic.pdf", width = 13, height = 13)
FeaturePlot(object = Tcell_nomix, features = c("RPF1", "GZMB", "NKG7", "CCL4",
                                         "CST7", "GZMA", "IFNG", "CCL3",
                                         "GZMK"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 3)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/dysfunction.pdf", width = 17.5, height = 8.5)
FeaturePlot(object = Tcell_nomix, features = c("TIGIT", "PDCD1", "CXCL13", "LAG3",
                                         "HAVCR2", "CD38", "ENTPD1", "CTLA4"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 4)
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_markerPlot/Naive_Effector_or_memory.pdf", width = 17.5, height = 8.5)
FeaturePlot(object = Tcell_nomix, features = c("TIGIT", "PDCD1", "CXCL13", "LAG3",
                                         "HAVCR2", "CD38", "ENTPD1", "CTLA4"),cols = c("grey", "#A52A2A"), 
            reduction = "umap", pt.size = 0.2, ncol = 4)
dev.off()

# 7.3.2.Monocel -- CD4 cells ------------------------------------------

folder <- "5.2.Tcell_harmony/Monocle/CD4_Tcell"
if(!dir.exists(folder)){
  dir.create(folder)
}

CD4_Tcell <- subset(Tcell_nomix, subset = RNA_snn_res.0.8 %in% c(2,4,6,7,10))

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_CD4_T_UMAP.pdf",width = 4.5,height = 4)
DimPlot(CD4_Tcell, reduction = "umap", group.by = "RNA_snn_res.0.8", 
        pt.size = 0.1, label = T, cols = color) 
dev.off()

#构造表达及注释数据
exp.matrix <- as(as.matrix(CD4_Tcell@assays$RNA@data), "sparseMatrix")
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(exp.matrix)
exp_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- CD4_Tcell@meta.data
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
saveRDS(ordering_genes, "./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_ordering_genes.Rds")
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
#plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

#修改monocle对象中的列名示例
saveRDS(exp.monocle, "./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_exp.monocle.Rds")

#将不同分组情况的拟时序轨迹图画到一起

plot1<-plot_cell_trajectory(exp.monocle, color_by = "RNA_snn_res.0.8",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color[c(3,5,7,8,11)])
plot2<-plot_cell_trajectory(exp.monocle, color_by = "State",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
plot3<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size = 0.2, shuffle = T) + 
  scale_color_gradient(low = "#101D44", high = "#9FB2ED")
plot4<-plot_cell_trajectory(exp.monocle, color_by = "SampleType",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color[c(2,1)])
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_trajectory_plot.pdf",width = 12,height = 10)
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = NULL,ncol=2)
dev.off()

#split by clusters

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

#cluster
for (i in c(2,4,6,7,10)) {
  a <- rownames(CD4_Tcell@meta.data)[ CD4_Tcell@meta.data$RNA_snn_res.0.8 == i ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$RNA_snn_res.0.8 == i ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(paste0("Cluster ", i)) + 
    scale_color_gradient(low = "#101D44", high = "#9FB2ED")
  assign(paste0("P", i), plot)
}

Plist <- list(P2,P4,P6,P7,P10)
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_trajectory_plot_clusterSplit.pdf",width = 9,height = 5)
CombinePlots(plots = Plist,legend = NULL,ncol=3)
dev.off()

#split by sampletype

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

#cluster
for (i in c("EOPC", "LOPC")) {
  a <- rownames(CD4_Tcell@meta.data)[ CD4_Tcell@meta.data$SampleType == i ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, a] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(i)
  assign(paste0("P_", i), plot)
}

Plist <- list(P_EOPC,P_LOPC)
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_trajectory_plot_SampleTypeSplit.pdf",width = 6, height = 2.5)
CombinePlots(plots = Plist,legend = NULL,ncol=2)
dev.off()

#Staes_SampleType_Roe--------------------------------------

tab <- table(exp.monocle$SampleType, exp.monocle$State)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_State_Roe_SampleType.pdf", width = 4, height = 2)
pheatmap(Roe,
         cluster_cols = T,
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

#features plots

gene <- c("TIGIT", "FOXP3", "CTLA4", "LAG3", "HAVCR2",
          "CCR7", "SELL", "TCF7", "IL7R", "LEF1")

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_feature_naive_dysfunction_plot.pdf", width = 5, height = 5)
plot_genes_in_pseudotime(exp.monocle[gene,], color_by = "Pseudotime", ncol =2)
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_feature_naive_dysfunction_plot_ClusterGroup.pdf", width = 5, height = 5)
plot_genes_in_pseudotime(exp.monocle[gene,], color_by = "RNA_snn_res.0.8", ncol =2)
dev.off()

#Stimulatory and Inhibitory phenotype score (irGSEA) ------------------------------

library("irGSEA")

Stim_Inhi <- list()
Stim_Inhi$Stimulatory <- c("BTN3A1", "BTN3A2", "CCL5", "CD27", "CD28", "CD40", 
                           "CD40LG", "CD70", "CD80", "CX3CL1", "CXCL10", "CXCL9", 
                           "ENTPD1", "HMGB1", "ICAM1", "ICOS", "ICOSLG", "IFNA1", 
                           "IFNA2", "IFNG", "IL1A", "IL1B", "IL2", "IL2RA", "ITGB2", 
                           "PRF1", "SELP", "TLR4", "TNF", "TNFRSF14", "TNFRSF18", 
                           "TNFRSF4", "TNFRSF9", "TNFSF4", "TNFSF9", "C10orf54", 
                           "CD276", "CD48", "CD86", "HHLA2", "IL6", "IL6R", "KLRC1", 
                           "KLRK1", "LTA", "MICB", "NT5E", "PVR", "RAET1E", "TMEM173", 
                           "TMIGD2", "TNFRSF13B", "TNFRSF13C", "TNFRSF17", "TNFRSF25", 
                           "TNFRSF8", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", 
                           "TNFSF18", "ULBP1")
Stim_Inhi$Inhibitory <- c("ADORA2A", "ARG1", "BTLA", "CD274", "CD276", "CTLA4", 
                          "EDNRB", "HAVCR2", "IDO1", "IL10", "IL13", "IL4", "KIR2DL1", 
                          "KIR2DL2", "KIR2DL3", "LAG3", "PDCD1", "SLAMF7", "TGFB1", 
                          "TIGIT", "VEGFA", "VEGFB", "C10orf54", "VTCN1", "CD160", 
                          "CD244", "CD96", "CSF1R", "IL10RB", "KDR", "LGALS9", 
                          "PDCD1LG2", "PVRL2", "TGFBR1", "CD24", "CD47", "NKG2D", 
                          "B7H3", "TIM3", "GALECTIN9", "OX40", "OX40L", "GITR", 
                          "GITRL", "B7H4")


# calculate score
CD4_Tcell_irGSEA <- irGSEA.score(object = CD4_Tcell, assay = "RNA", slot = "data", 
                                 seeds = 123, ncores = 20, min.cells = 3, 
                                 min.feature = 0, custom = T, geneset = Stim_Inhi, 
                                 msigdb = F, method = c("AUCell", "UCell", "singscore", "ssgsea"),
                                 aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                 kcdf = 'Gaussian')

AUCell <- as.data.frame(t(as.data.frame(CD4_Tcell_irGSEA@assays$AUCell@data)))
AUCell$SampleType <- CD4_Tcell_irGSEA$SampleType

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_Scatter_Stim_Inhi_AUCell.pdf")
pal_for_sampletype <- c(EOPC = "#053E7ACC", LOPC = "#A52A2ACC")
names(pal_for_sampletype) <- c("EOPC", "LOPC")
ggplot(data = AUCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("AUCell")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_Scatter_Stim_Inhi_AUCell_nolegend.pdf")
pal_for_sampletype <- c(EOPC = "#053E7ACC", LOPC = "#A52A2ACC")
names(pal_for_sampletype) <- c("EOPC", "LOPC")
ggplot(data = AUCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("AUCell") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

#Violin plot
library(reshape2)
AUCell_plot = melt(AUCell)
colnames(AUCell_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/Nomix_Violin_Stim_Inhi_AUCell.pdf")
ggviolin(AUCell_plot, x="Status", y="value", fill = "SampleType", 
         ylab="AUCell Score",
         xlab="",
         palette = c("#053E7ACC", "#A52A2ACC"),
         width=0.6, add = c("boxplot"), 
         error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()


## Differences between sampletype (没跑) ----------------------------------------

fre <- table(CD4_Tcell$RNA_snn_res.0.8, CD4_Tcell$orig.ident)
sum <- table(CD4_Tcell$orig.ident)

pct <- c()
for (i in 1:nrow(fre)) {
  pct_tmp <- fre[i,]/sum
  pct <- rbind(pct, pct_tmp)
}
rownames(pct) <- paste0("c",c(0:13))

pct <- pct[c("c3", "c5", "c6", "c9", "c12"),]

pct <- t(pct)
pct <- as.data.frame(pct)

pct$type <- substr(rownames(pct), 1, 4)

library(ggpubr)

#把数据转换成gglpot2输入文件
data <- melt(pct, id.vars = c("type"))
colnames(data) <- c("Type","Group","Pct")

#绘制boxplot
p=ggviolin(data, x="Group", y="Pct", fill = "Type", 
           ylab="Pct. of each cluster in CD4+ T-cells",
           xlab="",
           palette = c("#3C5488","#E64B35"),
           width=0.6, add = c("dotplot"), 
           error.plot = "crossbar", add.params = list(size=0.5)) 
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.05, 1), symbols = c("P < 0.05", " ")),
                        label = "p.signif")

#输出图片
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/CD4_difference.pdf", width=5, height=4)
print(p1)
dev.off()

#FOXP3 & CTLA4
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/DotPlot_FOXP3_CTLA4.pdf", width=5, height=4)
DotPlot(CD4_Tcell, features=c("FOXP3", "CTLA4"), group.by="SampleType")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/VlnPlot_FOXP3_CTLA4.pdf", width=3, height=4)
VlnPlot(CD4_Tcell, features=c("FOXP3", "CTLA4"), group.by="SampleType")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/FOXP3.pdf", width=3, height=4)
VlnPlot(CD4_Tcell, features=c("FOXP3", "CTLA4"))
dev.off()

#FOXP3_proportion

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/FOXP3_positive.pdf")
CD4_Tcell$FOXP3_positive <- ifelse(CD4_Tcell@assays$RNA@data["FOXP3",] > 0,"Positive","Negative")
DimPlot(CD4_Tcell, split.by="SampleType", group.by="FOXP3_positive", pt.size=0.5, cols=c("grey", "#E64B35"))
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/FOXP3_positive_bar.pdf", width=3, height=4)
FOXP3_pos <- table(CD4_Tcell$SampleType, CD4_Tcell$FOXP3_positive)
sum <- apply(FOXP3_pos, 1, sum)
FOXP3_pos_pct <- data.frame()
for (n in 1:2) {
  FOXP3_pos_pct = rbind( FOXP3_pos_pct, FOXP3_pos[n,] / sum[n] )  
} 
colnames(FOXP3_pos_pct) <- colnames(FOXP3_pos)
FOXP3_pos_pct$group <- rownames(FOXP3_pos)

library(reshape2)
FOXP3_pos_plot = melt(FOXP3_pos_pct)
colnames(FOXP3_pos_plot) = c('group','positive','percent')

ggplot( FOXP3_pos_plot, aes( x = group, weight = percent, fill = positive))+
  geom_bar( position = "stack") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
dev.off()

#CTLA4
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/CTLA4_positive.pdf")
CD4_Tcell$CTLA4_positive <- ifelse(CD4_Tcell@assays$RNA@data["CTLA4",] > 0,"Positive","Negative")
DimPlot(CD4_Tcell, split.by="SampleType", group.by="CTLA4_positive", pt.size=0.5, cols=c("grey", "#E64B35"))
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/CTLA4_positive_bar.pdf", width=3, height=4)
CTLA4_pos <- table(CD4_Tcell$SampleType, CD4_Tcell$CTLA4_positive)
sum <- apply(CTLA4_pos, 1, sum)
CTLA4_pos_pct <- data.frame()
for (n in 1:2) {
  CTLA4_pos_pct = rbind( CTLA4_pos_pct, CTLA4_pos[n,] / sum[n] )  
} 
colnames(CTLA4_pos_pct) <- colnames(CTLA4_pos)
CTLA4_pos_pct$group <- rownames(CTLA4_pos)

library(reshape2)
CTLA4_pos_plot = melt(CTLA4_pos_pct)
colnames(CTLA4_pos_plot) = c('group','positive','percent')

ggplot( CTLA4_pos_plot, aes( x = group, weight = percent, fill = positive))+
  geom_bar( position = "stack") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
dev.off()

# 7.3.3.Monocel -- CD8 cells ------------------------------------------

folder <- "5.2.Tcell_harmony/Monocle/CD8_Tcell"
if(!dir.exists(folder)){
  dir.create(folder)
}

CD8_Tcell <- subset(Tcell_nomix, subset = RNA_snn_res.0.8 %in% c(0,1,3,5,9,12))

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Nomix_CD8_T_UMAP.pdf",width = 4.5,height = 4)
DimPlot(CD8_Tcell, reduction = "umap", group.by = "RNA_snn_res.0.8", 
        pt.size = 0.1, label = T, cols = color) 
dev.off()

#构造表达及注释数据
exp.matrix <- as(as.matrix(CD8_Tcell@assays$RNA@data), "sparseMatrix")
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(exp.matrix)
exp_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- CD8_Tcell@meta.data
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
saveRDS(ordering_genes, "./5.2.Tcell_harmony/Monocle/CD8_Tcell/Nomix_ordering_genes.Rds")
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
#plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

#修改monocle对象中的列名示例
saveRDS(exp.monocle, "./5.2.Tcell_harmony/Monocle/CD8_Tcell/Nomix_exp.monocle.Rds")

#将不同分组情况的拟时序轨迹图画到一起

plot1<-plot_cell_trajectory(exp.monocle, color_by = "RNA_snn_res.0.8",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color[c(1,2,4,6,10,13)])
plot2<-plot_cell_trajectory(exp.monocle, color_by = "State",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
plot3<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size = 0.2, shuffle = T) + 
  scale_color_gradient(low = "#101D44", high = "#9FB2ED")
plot4<-plot_cell_trajectory(exp.monocle, color_by = "SampleType",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color[c(2,1)])
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Nomix_trajectory_plot.pdf",width = 12,height = 10)
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = NULL,ncol=2)
dev.off()

#split by clusters

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

#cluster
for (i in c(0,1,3,5,9,12)) {
  a <- rownames(CD8_Tcell@meta.data)[ CD8_Tcell@meta.data$RNA_snn_res.0.8 == i ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$RNA_snn_res.0.8 == i ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(paste0("Cluster ", i)) + 
    scale_color_gradient(low = "#101D44", high = "#9FB2ED")
  assign(paste0("P", i), plot)
}

Plist <- list(P0,P1,P3,P5,P9,P12)
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Nomix_trajectory_plot_clusterSplit.pdf",width = 9,height = 5)
CombinePlots(plots = Plist,legend = NULL,ncol=3)
dev.off()

#split by sampletype

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

#cluster
for (i in c("EOPC", "LOPC")) {
  a <- rownames(CD8_Tcell@meta.data)[ CD8_Tcell@meta.data$SampleType == i ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$SampleType == i ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(i)
  assign(paste0("P_", i), plot)
}

Plist <- list(P_EOPC,P_LOPC)
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Nomix_trajectory_plot_SampleTypeSplit.pdf",width = 6, height = 2.5)
CombinePlots(plots = Plist,legend = NULL,ncol=2)
dev.off()

#Staes_SampleType_Roe--------------------------------------

tab <- table(exp.monocle$SampleType, exp.monocle$State)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Nomix_State_Roe_SampleType.pdf", width = 4, height = 2)
bk <- c(seq(0.5,1.5,by=0.01))
pheatmap(Roe,
         cluster_cols = T,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = T,
         show_rownames = T,
         scale="none",  #矫正
         breaks=bk,
         legend_breaks=seq(-0.5,0.5,0.5),
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Ro/e",
         angle_col = "0",
         display_numbers = T,
         legend = F)
dev.off()

gene <- c("TIGIT", "PDCD1", "LAG3", 
          "GZMA", "GZMK", "IFNG",
          "SELL", "TCF7", "IL7R")

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Nomix_feature_naive_dysfunction_plot.pdf", width = 5, height = 5)
plot_genes_in_pseudotime(exp.monocle[gene,], color_by = "RNA_snn_res.0.8", ncol =2)
dev.off()

#Stimulatory and Inhibitory phenotype score (irGSEA) ------------------------------

library("AUCell")

Stim_Inhi <- list()
Stim_Inhi$Stimulatory <- c("BTN3A1", "BTN3A2", "CCL5", "CD27", "CD28", "CD40", 
                           "CD40LG", "CD70", "CD80", "CX3CL1", "CXCL10", "CXCL9", 
                           "ENTPD1", "HMGB1", "ICAM1", "ICOS", "ICOSLG", "IFNA1", 
                           "IFNA2", "IFNG", "IL1A", "IL1B", "IL2", "IL2RA", "ITGB2", 
                           "PRF1", "SELP", "TLR4", "TNF", "TNFRSF14", "TNFRSF18", 
                           "TNFRSF4", "TNFRSF9", "TNFSF4", "TNFSF9", "C10orf54", 
                           "CD276", "CD48", "CD86", "HHLA2", "IL6", "IL6R", "KLRC1", 
                           "KLRK1", "LTA", "MICB", "NT5E", "PVR", "RAET1E", "TMEM173", 
                           "TMIGD2", "TNFRSF13B", "TNFRSF13C", "TNFRSF17", "TNFRSF25", 
                           "TNFRSF8", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", 
                           "TNFSF18", "ULBP1")
Stim_Inhi$Inhibitory <- c("ADORA2A", "ARG1", "BTLA", "CD274", "CD276", "CTLA4", 
                          "EDNRB", "HAVCR2", "IDO1", "IL10", "IL13", "IL4", "KIR2DL1", 
                          "KIR2DL2", "KIR2DL3", "LAG3", "PDCD1", "SLAMF7", "TGFB1", 
                          "TIGIT", "VEGFA", "VEGFB", "C10orf54", "VTCN1", "CD160", 
                          "CD244", "CD96", "CSF1R", "IL10RB", "KDR", "LGALS9", 
                          "PDCD1LG2", "PVRL2", "TGFBR1", "CD24", "CD47", "NKG2D", 
                          "B7H3", "TIM3", "GALECTIN9", "OX40", "OX40L", "GITR", 
                          "GITRL", "B7H4")


# calculate score
CD8_Tcell_irGSEA <- irGSEA.score(object = CD8_Tcell, assay = "RNA", slot = "data", 
                                 seeds = 123, ncores = 20, min.cells = 3, 
                                 min.feature = 0, custom = T, geneset = Stim_Inhi, 
                                 msigdb = F, method = c("AUCell", "UCell", "singscore", "ssgsea"),
                                 aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                 kcdf = 'Gaussian')

saveRDS(CD8_Tcell_irGSEA, "5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Nomix_Stimu_Inhib_irGSEA.Rds")

AUCell <- as.data.frame(t(as.data.frame(CD8_Tcell_irGSEA@assays$AUCell@data)))
AUCell$SampleType <- CD8_Tcell_irGSEA$SampleType

#Scatter Plot 

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Nomix_Scatter_Stim_Inhi_AUCell.pdf")
pal_for_sampletype <- c(EOPC = "#053E7ACC", LOPC = "#A52A2ACC")
names(pal_for_sampletype) <- c("EOPC", "LOPC")
ggplot(data = AUCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("AUCell") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Nomix_Scatter_Stim_Inhi_AUCell_nolegend.pdf")
ggplot(data = AUCell, aes(x=Stimulatory, y=Inhibitory, color=SampleType)) +
  geom_point() +
  scale_color_manual(values = pal_for_sampletype) +
  geom_smooth(method = lm) +
  theme_classic() +
  ggtitle("AUCell") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
dev.off()

#Violin Plot 

library(reshape2)
AUCell_plot = melt(AUCell)
colnames(AUCell_plot) = c('SampleType','Status','value')
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Nomix_Violin_Stim_Inhi_AUCell.pdf")
ggviolin(AUCell_plot, x="Status", y="value", fill = "SampleType", 
         ylab="AUCell Score",
         xlab="",
         palette = c("#053E7ACC", "#A52A2ACC"),
         width=0.6, add = c("boxplot"), 
         error.plot = "crossbar") + rotate_x_text(60) +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")
dev.off()


## Differences between sampletype （没跑） ----------------------------------------

fre <- table(CD8_Tcell$RNA_snn_res.0.8, CD8_Tcell$orig.ident)
sum <- table(CD8_Tcell$orig.ident)

pct <- c()
for (i in 1:nrow(fre)) {
  pct_tmp <- fre[i,]/sum
  pct <- rbind(pct, pct_tmp)
}
rownames(pct) <- paste0("c",c(0:13))

pct <- pct[c("c0", "c1", "c2", "c4", "c7", "c13"),]

pct <- t(pct)
pct <- as.data.frame(pct)

pct$type <- substr(rownames(pct), 1, 4)

library(ggpubr)

#把数据转换成gglpot2输入文件
data <- melt(pct, id.vars = c("type"))
colnames(data) <- c("Type","Group","Pct")

#绘制boxplot
p=ggviolin(data, x="Group", y="Pct", fill = "Type", 
           ylab="Pct. of each cluster in CD8+ T-cells",
           xlab="",
           palette = c("#E64B35", "#3C5488"),
           width=0.6, add = c("dotplot"), 
           error.plot = "crossbar", add.params = list(size=0.5)) 
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.05, 1), symbols = c("P < 0.05", " ")),
                        label = "p.signif")

#输出图片
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/CD8_difference.pdf", width=5, height=4)
print(p1)
dev.off()

# select markers TIGIT LAG3 (dysfunction) IL7R and CCR7 (naive) GZMA IFNG (cytotoxic) KLRG1 B3GAT1 (senescence)---------------------------------

# exp

exp <- CD8_Tcell@assays$RNA@data[c("TIGIT","LAG3", "IL7R", "CCR7", "IFNG", "GZMA", "KLRG1", "B3GAT1"),] %>% as.data.frame() %>% t()
SampleType <- as.data.frame(CD8_Tcell@meta.data$SampleType)
colnames(SampleType) <- "SampleType"
exp <- cbind(exp, SampleType)
exp <- melt(exp)

P = ggplot(exp, aes(x = variable, y = value, fill = SampleType)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  ylab("Norm. Exp.") +
  xlab("") + 
  ggtitle("") +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

pdf("5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Boxplot_selected_genes.pdf", width = 6, height = 3)
print(P)
dev.off()

# proportion

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/LAG3_positive.pdf")
Plist <- list()
genes <- c("TIGIT","LAG3", "IL7R", "CCR7", "IFNG", "GZMA", "KLRG1", "B3GAT1")
for (i in genes) {
  CD8_Tcell$positive <- ifelse(CD8_Tcell@assays$RNA@data[i,] > 0, "Positive", "Negative")
  
  pos <- table(CD8_Tcell$SampleType, CD8_Tcell$positive)
  chi <- chisq.test(pos)
  p <- signif(chi$p.value,3)
  pos <- melt(pos)
  colnames(pos) <- c("SampleType", "Status", "nCells")
  
  barplot <- ggplot(pos, aes(x=SampleType, y=nCells,fill=Status)) +
    geom_bar(stat="identity", position = 'fill') +
    theme_bw() +
    scale_fill_manual(values = c("#BBBBBB", "#556B2FCC")) +
    annotate(geom="text", x=1.5, y=1.05, label=paste0("Chisq P-value = ", p), size = 2) +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ylab(paste0(i, " Pos. Pct.")) +
    xlab("") + 
    ggtitle("") +
    NoLegend()
  
  Plist[[i]] <- barplot
}

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Bar_selected_genes_positive.pdf", width=5, height=12)
grid.arrange(grobs = Plist, ncol = 2)
dev.off()

#dysfunction markers -------------------------------------

exp <- CD8_Tcell@assays$RNA@data[c("TIGIT","LAG3"),] %>% as.data.frame() %>% t()
SampleType <- as.data.frame(CD8_Tcell@meta.data$SampleType)
colnames(SampleType) <- "SampleType"
exp <- cbind(exp, SampleType)
exp <- melt(exp)

P = ggplot(exp, aes(x = variable, y = value, fill = SampleType)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  ylab("AUCell score") +
  xlab("") + 
  ggtitle("") +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

pdf("5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/Boxplot_TIGIT_LAG3.pdf", width = 4, height = 3)
print(P)
dev.off()

#LAG3
pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/LAG3_positive.pdf")
CD8_Tcell$LAG3_positive <- ifelse(CD8_Tcell@assays$RNA@data["LAG3",] > 0,"Positive","Negative")
DimPlot(CD8_Tcell, split.by="SampleType", group.by="LAG3_positive", pt.size=0.5, cols=c("grey", "#E46B4F"))
dev.off()

LAG3_pos <- table(CD8_Tcell$SampleType, CD8_Tcell$LAG3_positive)
chi <- chisq.test(LAG3_pos)
p <- signif(chi$p.value,3)
pos <- melt(LAG3_pos)
colnames(pos) <- c("SampleType", "Status", "nCells")

barplot <- ggplot(pos, aes(x=SampleType, y=nCells,fill=Status)) +
  geom_bar(stat="identity", position = 'fill') +
  theme_bw() +
  scale_fill_manual(values = c("#BBBBBB", "#556B2FCC")) +
  annotate(geom="text", x=1.5, y=1.05, label=paste0("Chisq P-value = ", p), size = 2) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("LAG3 Pos. Pct.") +
  xlab("") + 
  ggtitle("") +
  NoLegend()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/LAG3_positive_bar.pdf", width=2, height=3)
barplot
dev.off()

#TIGIT

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/TIGIT_positive.pdf")
CD8_Tcell$TIGIT_positive <- ifelse(CD8_Tcell@assays$RNA@data["TIGIT",] > 0,"Positive","Negative")
DimPlot(CD8_Tcell, split.by="SampleType", group.by="TIGIT_positive", pt.size=0.5, cols=c("grey", "#E46B4F"))
dev.off()

TIGIT_pos <- table(CD8_Tcell$SampleType, CD8_Tcell$TIGIT_positive)
chi <- chisq.test(TIGIT_pos)
p <- signif(chi$p.value,3)
pos <- melt(TIGIT_pos)
colnames(pos) <- c("SampleType", "Status", "nCells")

barplot <- ggplot(pos, aes(x=SampleType, y=nCells,fill=Status)) +
  geom_bar(stat="identity", position = 'fill') +
  theme_bw() +
  scale_fill_manual(values = c("#BBBBBB", "#556B2FCC")) +
  annotate(geom="text", x=1.5, y=1.05, label=paste0("Chisq P-value = ", p), size = 2) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("TIGIT Pos. Pct.") +
  xlab("") + 
  ggtitle("") +
  NoLegend()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/Stimu_Inhib/TIGIT_positive_bar.pdf", width=2, height=3)
barplot
dev.off()

#Naive markers IL7R and SELL ----------------------------------------

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/DotPlot_IL7R_SELL.pdf")
DotPlot(CD8_Tcell, features=c("IL7R", "SELL"), group.by="SampleType")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/VlnPlot_IL7R_SELL.pdf")
VlnPlot(CD8_Tcell, features=c("IL7R", "SELL"), group.by="SampleType")
dev.off()

#senescence KLRG1 and B3GAT1  -------------------------------------

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/DotPlot_KLRG1_B3GAT1.pdf")
DotPlot(CD8_Tcell, features=c("KLRG1", "B3GAT1"), group.by="SampleType")
dev.off()

pdf("./5.2.Tcell_harmony/Monocle/CD8_Tcell/VlnPlot_KLRG1_B3GAT1.pdf")
VlnPlot(CD8_Tcell, features=c("KLRG1", "B3GAT1"), group.by="SampleType")
dev.off()

# 7.3.4. SingleR注释 ----------------------------------------

library(SingleR)
#celldex::DatabaseImmuneCellExpressionData()
#An external refdata is also acceptable 
load("./DatabaseImmuneCellExpressionData.Rdata")
testdata <- GetAssayData(Tcell, slot="data")
#testdata <- PCSC@assays$SCT@data
clusters <- Tcell@meta.data$RNA_snn_res.0.8

#按cluster进行注释
cellpred <- SingleR(test = testdata, ref = refdata, 
                    labels = refdata$label.main, 
                    clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype <- data.frame(ClusterID=rownames(cellpred), 
                       celltype=cellpred$labels, stringsAsFactors = F)

Tcell@meta.data$singleR_celltype = "Unknown"
for(i in 1:nrow(celltype)){
  Tcell@meta.data[which(Tcell@meta.data$seurat_clusters == celltype$ClusterID[i]),'singleR_celltype'] <- celltype$celltype[i]
}

plot_anno <- DimPlot(Tcell, group.by="singleR_celltype", pt.size=0.1, label=TRUE, label.size=5, reduction='umap', shuffle = T)
ggsave(filename = "./5.2.Tcell_harmony/CellCluster-UMAPPlot_res0.8_SingleR.pdf",
       plot = plot_anno)

# 7.3.5. rename -----------------------------------------

celltype=data.frame(ClusterID=0:13,
                    celltype='CD8_Tem')   
celltype[celltype$ClusterID %in% c(0),2]='CD8_Tn'  
celltype[celltype$ClusterID %in% c(2),2]='CD4_Tn'
celltype[celltype$ClusterID %in% c(4),2]='Tcm'
celltype[celltype$ClusterID %in% c(6),2]='Treg'
celltype[celltype$ClusterID %in% c(7),2]='Tcm'
celltype[celltype$ClusterID %in% c(8),2]='NKT'
celltype[celltype$ClusterID %in% c(10),2]='CD4_Tn'
celltype[celltype$ClusterID %in% c(11),2]='NK'
head(celltype)

Tcell_nomix@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  Tcell_nomix@meta.data[which(Tcell_nomix@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(Tcell_nomix@meta.data$celltype)

Idents(Tcell_nomix) <- "celltype"
Tcell_nomix <- RenameIdents(Tcell_nomix,
                      "Tcm" = "Tcm",
                      "CD4_Tn" = "CD4_Tn",
                      "Treg" = "Treg",
                      "CD8_Tn" = "CD8_Tn",
                      "CD8_Tem" = "CD8_Tem",
                      "NKT" = "NKT",
                      "NK" = "NK")
Tcell_nomix$celltype <- Idents(Tcell_nomix)

pdf("./5.2.Tcell_harmony/Nomix_Tcell_celltype_umap.pdf", width = 5, height = 5)
DimPlot(Tcell_nomix, label = TRUE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color, repel = T) + 
  NoLegend() + 
  theme(plot.title = element_blank())
dev.off()

pdf("./5.2.Tcell_harmony/Nomix_Tcell_celltype_umap_noLabel.pdf", width = 3, height = 3)
DimPlot(Tcell_nomix, label = FALSE, pt.size=0.3, group.by = 'celltype', 
        shuffle = T , cols = color) + 
  NoLegend() + 
  theme(plot.title = element_blank())
dev.off()

png("./5.2.Tcell_harmony/Nomix_Tcell_celltype_umap.png", width = 360, height = 360)
DimPlot(Tcell_nomix, label = TRUE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color, repel = T) + 
  NoLegend() + 
  theme(plot.title = element_blank())
dev.off()

saveRDS(Tcell_nomix, "./5.2.Tcell_harmony/Nomix_Tcell_annotated.Rds")

# Roe

tab <- table(Tcell_nomix@meta.data$celltype, Tcell_nomix@meta.data$SampleType)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.2.Tcell_harmony/Nomix_Roe_celltype_SampleType.pdf", width = 2.5, height = 4)
bk <- c(seq(-0.5,0.5,by=0.01))
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

# Allmarker bubble

Features <- c("CD3E", "CD4", "CD8A", 
              "IFNG", "GPR183",
              "TCF7", "IL7R", 
              "CTLA4", "FOXP3",
              "LAG3", "GZMB", 
              "FCGR3A", "KLRK1",
              "NCAM1","KLRB1")

pdf("./5.2.Tcell_harmony/Nomix_AllmarkerBubble.pdf",width = 6, height = 4.8)
jjDotPlot(Tcell_nomix,
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

pdf("./5.2.Tcell_harmony/Nomix_Rotated_AllmarkerBubble.pdf",width = 2.3, height = 3.8)
jjDotPlot(Tcell_nomix,
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


# 8.findmarkers celltype ---------
plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 1000 * 1024^2)

Idents(Tcell_nomix) <- "celltype"
Allmarkers <- FindAllMarkers(Tcell_nomix,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = FALSE)
saveRDS(Allmarkers, "./5.2.Tcell_harmony/Markers/Nomix_celltype_markers.Rds")

sigmarkers <- Allmarkers[which(Allmarkers$avg_log2FC > 0.25 & 
                                 Allmarkers$p_val_adj < 0.01 &
                                 Allmarkers$pct.1 > 0.3),]

sigmarkers_order <- c()
for (i in as.character(unique(Allmarkers$cluster))) {
  
  tmp <- sigmarkers[which(sigmarkers$cluster == i),]
  tmp <- tmp[order(tmp$avg_log2FC, decreasing = T),]
  sigmarkers_order <- rbind(sigmarkers_order, tmp)
  
}

write.table(sigmarkers_order, "./5.2.Tcell_harmony/Markers/Nomix_celltype_sigmarkers.txt",
            quote = F, sep = "\t", row.names = F)

topmarkers <- Allmarkers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

#Heatmap
DoHeatmap(Tcell_nomix, features = topmarkers$gene, group.by = "celltype")

#Heatmap by cluster

Cluster <- data.frame(Cluster = Idents(Endo))
data <- as.matrix(Endo@assays$RNA@scale.data)
data <- data[rownames(data) %in% topmarkers$gene,]
data <- as.data.frame(t(data))
data <- merge(Cluster, data, by='row.names')

data_list <- split(data, data$Cluster)
mean <- c()
for (i in 1:length(data_list)) {
  mean_tmp <- apply(data_list[[i]][,3:ncol(data_list[[i]])], 2, mean)
  mean <- cbind(mean, mean_tmp)
}
colnames(mean) <- paste0("c", c(0:8))

library(RColorBrewer)
pdf("./5.3.Endo_harmony/Markers/Nomix_heatmap_top10Marker.pdf")
heatmap(mean, Colv = NA, cexCol = 1, col = colorRampPalette(brewer.pal(9, "OrRd"))(50))
dev.off()

# 7. functional annotation -------------------------------------------

# 7.1.irGSEA Hallmark ---------------------------------------------

# calculate score
Idents(Tcell) <- "SampleType"
Tcell <- irGSEA.score(object = Tcell, assay = "RNA", slot = "data", 
                      seeds = 123, ncores = 20, min.cells = 3, 
                      min.feature = 0, custom = F, geneset = NULL, 
                      msigdb = T, species = "Homo sapiens", category = "H",  
                      subcategory = NULL, geneid = "symbol",
                      method = c("AUCell", "UCell", "singscore", "ssgsea"),
                      aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                      kcdf = 'Gaussian')

saveRDS(Tcell, "./5.2.Tcell_harmony/irGSEA/HALLMARK/Tcell_irGSEA_sampletyep.Rds")

#VlnPlot

library(gridExtra)

#GLYCOLYSIS
plot1 <- halfvlnplot.AUCell <- irGSEA.halfvlnplot(object = Tcell,
                                                 method = "AUCell",
                                                 color.cluster = c("#E64B35", "#3C5488"),
                                                 show.geneset = "HALLMARK-GLYCOLYSIS") + NoLegend() + labs(title="")
plot2 <- halfvlnplot.AUCell <- irGSEA.halfvlnplot(object = Tcell,
                                                 method = "UCell",
                                                 color.cluster = c("#E64B35", "#3C5488"),
                                                 show.geneset = "HALLMARK-GLYCOLYSIS") + NoLegend() + labs(title="")
plot3 <- halfvlnplot.AUCell <- irGSEA.halfvlnplot(object = Tcell,
                                                 method = "ssgsea",
                                                 color.cluster = c("#E64B35", "#3C5488"),
                                                 show.geneset = "HALLMARK-GLYCOLYSIS") + NoLegend() + labs(title="")
plot4 <- halfvlnplot.AUCell <- irGSEA.halfvlnplot(object = Tcell,
                                                 method = "singscore",
                                                 color.cluster = c("#E64B35", "#3C5488"),
                                                 show.geneset = "HALLMARK-GLYCOLYSIS") + NoLegend() + labs(title="")
plot <- grid.arrange(grobs = list(plot1,plot2,plot3,plot4), ncol=2, top = "HALLMARK-GLYCOLYSIS")
ggsave("./5.2.Tcell_harmony/irGSEA/HALLMARK/VlnPlot_Glycolysis.pdf", plot)

#hypoxia
plot1 <- halfvlnplot.AUCell <- irGSEA.halfvlnplot(object = Tcell,
                                                  method = "AUCell",
                                                  color.cluster = c("#E64B35", "#3C5488"),
                                                  show.geneset = "HALLMARK-HYPOXIA") + NoLegend() + labs(title="")
plot2 <- halfvlnplot.AUCell <- irGSEA.halfvlnplot(object = Tcell,
                                                  method = "UCell",
                                                  color.cluster = c("#E64B35", "#3C5488"),
                                                  show.geneset = "HALLMARK-HYPOXIA") + NoLegend() + labs(title="")
plot3 <- halfvlnplot.AUCell <- irGSEA.halfvlnplot(object = Tcell,
                                                  method = "ssgsea",
                                                  color.cluster = c("#E64B35", "#3C5488"),
                                                  show.geneset = "HALLMARK-HYPOXIA") + NoLegend() + labs(title="")
plot4 <- halfvlnplot.AUCell <- irGSEA.halfvlnplot(object = Tcell,
                                                  method = "singscore",
                                                  color.cluster = c("#E64B35", "#3C5488"),
                                                  show.geneset = "HALLMARK-HYPOXIA") + NoLegend() + labs(title="")
plot <- grid.arrange(grobs = list(plot1,plot2,plot3,plot4), ncol=2, top = "HALLMARK-HYPOXIA")
ggsave("./5.2.Tcell_harmony/irGSEA/HALLMARK/VlnPlot_HYPOXIA.pdf", plot)

pdf("./5.2.Tcell_harmony/irGSEA/HALLMARK/VlnPlot_ssgsea_Glycolysis_celltypeSplit.pdf")
DefaultAssay(Tcell) <- "ssgsea"
VlnPlot(Tcell, "HALLMARK-GLYCOLYSIS",group.by="celltype",split.by="SampleType",
        cols = c("#E64B35", "#3C5488"), pt.size=0)
dev.off()

pdf("./5.2.Tcell_harmony/irGSEA/HALLMARK/BoxPlot_ssgsea_Hypoxia_celltypeSplit.pdf")
exp <- data.frame(SampleType=Tcell$SampleType, celltype=Tcell$celltype, HYPOXIA=Tcell@assays$ssgsea@data["HALLMARK-HYPOXIA",])
rownames(exp) <- colnames(Tcell)
p=ggboxplot(exp, x="celltype", y="HYPOXIA", fill = "SampleType", 
            ylab="ssGSEA HYPOXIA score",
            xlab="",
            palette = c("#E64B35", "#3C5488"))
p=p + rotate_x_text(45) + theme(axis.title.y = element_text(vjust=2, size=18))
p1=p+stat_compare_means(aes(group=SampleType),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","")),
                        label = "p.signif")
p1
dev.off()

# 7.2.scFEA analyses -----------------------------------------------

scFEA <- read.csv("scFEA/output/balance_20221010-150304.csv", row.names =1)
scFEA <- t(scFEA)
Tcell[["FEA"]] <- CreateAssayObject(counts = scFEA)

DefaultAssay(Tcell) <- 'FEA'

Idents(Tcell) <- "SampleType"
FEA <- scale(scFEA)
cells <- data.frame(group = as.character(Idents(Tcell)))
rownames(cells) <- colnames(Tcell)
FEA <- merge(cells, FEA, by = 0)
FEA <- column_to_rownames(FEA, "Row.names")
FEA <- FEA[,which(colnames(FEA) %in% c("Glucose", "Citrate", "Pyruvate", "Acetyl.CoA", 
                                "Malate", "Lactate", "Succinate", "Oxaloacetate",
                                "Fumarate", "Succinyl.CoA", "group"))]

FEA_list <- split(FEA, FEA$group)

mean <- c()
for (i in 1:length(FEA_list)) {
  mean_tmp <- apply(FEA_list[[i]][,2:ncol(FEA_list[[i]])], 2, mean)
  mean <- cbind(mean, mean_tmp)
}
colnames(mean) <- names(FEA_list)

pdf("./5.2.Tcell_harmony/scFEA/scFEA_heatmap.pdf", height=10, width = 6)
pheatmap(mean, show_colnames =T, show_rownames = T, cluster_cols = FALSE)
dev.off()

saveRDS(Tcell, "./5.2.Tcell_harmony/Tcell_annoted_FEA.Rds")

#Lactate

Tcell <- ScaleData(Tcell)

exp <- data.frame(SampleType=Tcell$SampleType, celltype=Tcell$celltype, Lactate=Tcell@assays$FEA@scale.data["Lactate",])
rownames(exp) <- colnames(Tcell)
exp <- exp[which(exp$celltype != "Epi_contaminate"),]

library(ggpubr)

p=ggboxplot(exp, x="celltype", y="Lactate", fill = "SampleType", 
           ylab="Scaled Lactate Flux",
           xlab="",
           palette = c("#E64B35", "#3C5488"))  + scale_y_continuous(limits = c(-2.5,1.5))
p=p + rotate_x_text(45) + theme(axis.title.y = element_text(vjust=2, size=18))
p1=p+stat_compare_means(aes(group=SampleType),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","")),
                        label = "p.signif")

#输出图片
pdf("./5.2.Tcell_harmony/scFEA/Lactate.pdf", width=6, height=5)
print(p1)
dev.off()

#CTLA4 boxplot

DefaultAssay(Tcell) <- 'RNA'

exp <- data.frame(SampleType=Tcell$SampleType, celltype=Tcell$celltype, CTLA4=Tcell@assays$RNA@scale.data["CTLA4",])
rownames(exp) <- colnames(Tcell)
exp <- exp[which(exp$celltype != "Epi_contaminate"),]

p=ggboxplot(exp, x="celltype", y="CTLA4", fill = "SampleType", 
            ylab="CTLA4 Expression",
            xlab="",
            palette = c("#E64B35", "#3C5488"))
p=p + rotate_x_text(45) + theme(axis.title.y = element_text(vjust=2, size=18))
p1=p+stat_compare_means(aes(group=SampleType),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","")),
                        label = "p.signif")
pdf("./5.2.Tcell_harmony/scFEA/CTLA4.pdf", width=6, height=5)
print(p1)
dev.off()

exp <- data.frame(SampleType=Tcell$SampleType, celltype=Tcell_scale$celltype, CD28=Tcell_scale@assays$RNA@scale.data["CD28",])
rownames(exp) <- colnames(Tcell)
exp <- exp[which(exp$celltype != "Epi_contaminate"),]

p=ggboxplot(exp, x="celltype", y="CD28", fill = "SampleType", 
            ylab="CD28 Expression",
            xlab="",
            palette = c("#E64B35", "#3C5488"))
p=p + rotate_x_text(45) + theme(axis.title.y = element_text(vjust=2, size=18))
p1=p+stat_compare_means(aes(group=SampleType),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","")),
                        label = "p.signif")
pdf("./5.2.Tcell_harmony/scFEA/CD28.pdf", width=6, height=5)
print(p1)
dev.off()


# 7.3.TGFB ----------------------------------------------
exp <- data.frame(SampleType=Tcell$SampleType, celltype=Tcell$celltype, TGFB1=Tcell@assays$RNA@data["TGFB1",])
rownames(exp) <- colnames(Tcell)
exp <- exp[which(exp$celltype != "Epi_contaminate"),]

p=ggboxplot(exp, x="celltype", y="TGFB1", fill = "SampleType", 
            ylab="TGFB1 Expression",
            xlab="",
            palette = c("#E64B35", "#3C5488"))
p=p + rotate_x_text(45) + theme(axis.title.y = element_text(vjust=2, size=18))
p1=p+stat_compare_means(aes(group=SampleType),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","")),
                        label = "p.signif")




#glycolysis genes

saveRDS(Tcell_scale, "./5.2.Tcell_harmony/Tcell_scaleall.Rds")
genes <- c("CTLA4", "SLC16A3", "PFKP", "PFKL", "PGAM2", "ENO3", "SLC2A1", "HK2",
           "ENO2", "PKM", "BSG", "PKLR", "PGAM1", "ALDOA", "PGK1", "SLC16A1", "LDHA",
           "PFKM", "TPI1", "LDHB", "ENO1", "GAPDH", "HK1", "GPI")
mtx <- as.data.frame(t(as.matrix(Tcell@assays$RNA@data[rownames(Tcell@assays$RNA@data) %in% genes,])))
mtx$celltype <- Tcell$celltype

library(corrplot)
Treg_mtx <- mtx[mtx$celltype == "Treg",][, -ncol(mtx)]
Treg_cor <- cor.mtest(Treg_mtx, conf.level = 0.95)
pdf("correlation.pdf",height=16,width=16)
corrplot(corr=cor(Treg_mtx),
         method = "circle",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         p.mat = Treg_cor$p,
         sig.level = 0.001,
         insig = "pch",
         number.cex = 1,
         type = "upper",
         col=colorRampPalette(c("blue", "white", "red"))(50),
)
dev.off()
# 6.5.MEBOCOST --------------------------------------


# 6.