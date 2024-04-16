library(Matrix)
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

# Myeloid --------------------------------------------------------------------

color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")
                    

folder <- "5.4.Mono_harmony"
if(!dir.exists(folder)){
  dir.create(folder)
}

# 1. subset, harmony and clustering ----------------------------------------------------------------

Mono <- NormalizeData(Mono)
Mono <- FindVariableFeatures(Mono,
                                   selection.method = "vst",
                                   nfeatures = 2000,
                                   verbose = FALSE) %>%
  ScaleData(vars.to.regress = c("percent.mt"))
Mono <- RunPCA(Mono)

pdf("./5.4.Mono_harmony/Nomix_harmony_convergence.pdf")
Mono <- RunHarmony(Mono, "orig.ident", plot_convergence = T)
dev.off()

Mono <- RunUMAP(Mono, reduction = "harmony", dims = 1:20)

Mono <- FindNeighbors(Mono, reduction = "harmony", 
                            dims = 1:20)
for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1)) {
  Mono <- FindClusters(Mono, resolution = res)}

apply(Mono@meta.data[, grep("RNA_snn_res", colnames(Mono@meta.data))],2,table)

library(clustree)

p2_tree=clustree(Mono@meta.data, prefix = "RNA_snn_res.") +
  scale_color_manual(values = color)
ggsave(plot=p2_tree, filename="./5.4.Mono_harmony/Nomix_Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(Mono) <- "RNA_snn_res.0.8"

## plotting ------------------------------------

pdf("./5.4.Mono_harmony/Nomix_CellCluster-UMAPPlot_res0.8.pdf",width = 4,height = 4)
DimPlot(Mono, reduction = "umap", pt.size = 0.3, label = T, shuffle = T) + NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.4.Mono_harmony/Nomix_CellCluster-UMAPPlot_SampleType.pdf",width = 4,height = 4)
DimPlot(Mono, reduction = "umap", group.by = "SampleType", pt.size = 0.3, shuffle = T) +
  scale_color_manual(values = color)
dev.off()

pdf("./5.4.Mono_harmony/Nomix_CellCluster-UMAPPlot_SampleType_noLegend.pdf",width = 4,height = 4)
DimPlot(Mono, reduction = "umap", group.by = "SampleType", pt.size = 0.3, shuffle = T) +
  NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.4.Mono_harmony/Nomix_CellCluster-UMAPPlot_Sample.pdf",width = 4,height = 4)
DimPlot(Mono, reduction = "umap", group.by = "orig.ident", pt.size = 0.3, label = F, shuffle = T) +
  scale_color_manual(values = color)
dev.off()

pdf("./5.4.Mono_harmony/Nomix_CellCluster-UMAPPlot_Sample_noLegend.pdf",width = 4,height = 4)
DimPlot(Mono, reduction = "umap", group.by = "orig.ident", pt.size = 0.3, label = F, shuffle = T) + 
  NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.4.Mono_harmony/Nomix_CellCluster-UMAPPlot_SampleTypeSplit.pdf",width = 8,height = 4)
DimPlot(Mono, reduction = "umap", split.by = "SampleType", pt.size = 0.3, label = T) + 
  NoLegend() +
  scale_color_manual(values = color)
dev.off()

pdf("./5.4.Mono_harmony/cellcycle/Nomix_CellCluster-UMAPPlot_cellcycle_res0.3.pdf",width = 4,height = 4)
DimPlot(Mono, reduction = "umap", group.by = "Phase", pt.size = 0.3, label = T) +
  scale_color_manual(values = color)
dev.off()

pdf("./5.4.Mono_harmony/cellcycle/Nomix_CellCluster-UMAPPlot_SampleTypeSplit_cellcycle_res0.3.pdf",width = 8,height = 4)
DimPlot(Mono, reduction = "umap", group.by = "Phase", split.by = "SampleType", pt.size = 0.3, label = T) +
  scale_color_manual(values = color)
dev.off()

saveRDS(Mono, "./5.4.Mono_harmony/Nomix_Myeloid_harmonied.Rds")

## FindMarkers ------------------------------------

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)

Allmarkers <- FindAllMarkers(Mono,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = FALSE)
saveRDS(Allmarkers, "./5.4.Mono_harmony/Markers/Nomix_cluster_markers_nomix.Rds")

topmarkers <- Allmarkers %>% group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

#Heatmap by cluster

Cluster <- data.frame(ID = colnames(Mono), Cluster = Mono$RNA_snn_res.0.8)
Cluster <- Cluster[order(Cluster$Cluster),]
annotation <- data.frame(Cluster = Cluster$Cluster)
rownames(annotation) <- Cluster$ID

top <- data.frame(gene = topmarkers$gene, Cluster = topmarkers$cluster)
top <- top[order(top$Cluster),]

data <- as.matrix(Mono@assays$RNA@data)
genes <- top$gene
data <- data[rownames(data) %in% genes,]
data <- data[match(genes, rownames(data)),]
data <- data[, Cluster$ID]

Cluster_colors <- color
names(Cluster_colors) <- as.character(levels(Mono$RNA_snn_res.0.8))
ann_colors <- list(Cluster = Cluster_colors)

pdf("./5.4.Mono_harmony/Markers/Nomix_cluster_heatmap_top10_nomix.pdf", width = 5, height = 15)
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

pdf("./5.4.Mono_harmony/Markers/Nomix_cluster_heatmap_top8_nomix.pdf", width = 5, height = 10)
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

png("./5.4.Mono_harmony/Markers/Nomix_cluster_heatmap_top8_nomix.png", units = "in", width = 5, height = 11.4, res = 288)
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
         fontsize = 7.6,
         fontsize_row = 7.6,
         fontsize_col = 12,
         main = "",
         angle_col = "90",
         #border_color = "grey",
         annotation_col = annotation,
         annotation_colors = ann_colors,
         legend = F)
dev.off()

png("./5.4.Mono_harmony/Markers/Nomix_cluster_heatmap_top5_nomix.png", units = "in", width = 5, height = 6, res = 288)
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
         fontsize = 7.6,
         fontsize_row = 7.6,
         fontsize_col = 12,
         main = "",
         angle_col = "90",
         #border_color = "grey",
         annotation_col = annotation,
         annotation_colors = ann_colors,
         legend = F)
dev.off()

## Pan markers -----------------------------------

Marker <- c("JCHAIN", "CD1C",  ##DC
            "VCAN", "S100A8",  #Mono
            "CD86", "CD163"  #Macro
)

png("./5.4.Mono_harmony/Markers/Nomix_Pan_Markers.png", res = 300, units = 'in', width = 7, height = 9)
FeaturePlot(Mono, features = Marker, ncol = 2, cols = c("grey88", "#A52A2A"))
dev.off()

## Markers -----------------------------------------

Features <- c("CD68", "LYZ", 
              "CD86", "CD163", "C1QC",  #Macrophage
              "CD14", "S100A8", "FCGR3A",  #Monocyte
              "ITGAX", "CD83",  #Dendritic
              "JCHAIN", "LILRA4",  #pDC
              "CLEC9A", "IDO1",  #cDC1
              "CD1C", "FCER1A",  #cDC2
              "ACPP", "KLK3", "PECAM1", "LTB", "CD3E")

library(scRNAtoolVis)

pdf("./5.4.Mono_harmony/Nomix_AllmarkerBubble_cluster.pdf",width = 10, height = 6)
jjDotPlot(Mono,
          gene = Features,
          id = "RNA_snn_res.0.8",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0.5,
          dot.max = 5,
          x.text.vjust = 0.5)
dev.off()

## annotation --------------------------------------

celltype=data.frame(ClusterID=0:14,
                    celltype='Macro_APOE')   
celltype[celltype$ClusterID %in% c(0),2]='Macro_FTL'
celltype[celltype$ClusterID %in% c(1,7),2]='cDC2_CD1C'
celltype[celltype$ClusterID %in% c(4),2]='Macro_FOLR2'
celltype[celltype$ClusterID %in% c(5),2]='Mono_FCN1'
celltype[celltype$ClusterID %in% c(6),2]='Mono_IL1B'
celltype[celltype$ClusterID %in% c(8),2]='Macro_C3'
celltype[celltype$ClusterID %in% c(9),2]='Macro_CXCL11'
celltype[celltype$ClusterID %in% c(10),2]='Macro_RGS1'
celltype[celltype$ClusterID %in% c(11),2]='Macro_MT1A'
celltype[celltype$ClusterID %in% c(12),2]='cDC1_CLEC9A'
celltype[celltype$ClusterID %in% c(13),2]='pDC_JCHAIN'
celltype[celltype$ClusterID %in% c(14),2]='cDC3_LAMP3'
head(celltype)

Mono@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  Mono@meta.data[which(Mono@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(Mono@meta.data$celltype)

Mono$celltype <- factor(Mono$celltype, levels = 
                                c('pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
                                  'Mono_IL1B', 'Mono_FCN1',
                                  'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                  'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL'))
Idents(Mono) <- "celltype"

pdf("./5.4.Mono_harmony/Nomix_CellType-UMAPPlot_nomix_res0.8.pdf", width = 5, height = 5)
DimPlot(Mono, label = TRUE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color, repel = T) + 
  NoLegend() + 
  ggtitle("")
dev.off()

pdf("./5.4.Mono_harmony/Nomix_CellType-UMAPPlot_nomix_NoLebel_res0.8.pdf", width = 3.5, height = 3.5)
DimPlot(Mono, label = FALSE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color, repel = T) + 
  NoLegend() + 
  ggtitle("")
dev.off()

pdf("./5.4.Mono_harmony/Nomix_CellType-UMAPPlot_nomix_Legend_res0.8.pdf", width = 5, height = 5)
DimPlot(Mono, label = FALSE, pt.size=0.3, group.by = 'celltype', shuffle = T, 
        label.size = 3, label.color = "black", label.box = T, cols = color, repel = T) +
  ggtitle("")
dev.off()

## FindAllMarkers ------------------------------------

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)

Allmarkers <- FindAllMarkers(Mono,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = FALSE)
saveRDS(Allmarkers, "./5.4.Mono_harmony/Markers/cluster_markers_nomix.Rds")

## DotPlot ---------------------------------------

# pDC_JCHAIN  cDC1_CLEC9A    cDC2_CD1C   cDC3_LAMP3    Mono_IL1B    Mono_FCN1
# 30           36          814           25          355          362
# Macro_CXCL11     Macro_C3   Macro_RGS1  Macro_FOLR2   Macro_APOE   Macro_MT1A
# 228          303          190          396          839          115
# Macro_FTL
# 530



Features <- c("JCHAIN", "GZMB", "IL3RA",
              "CLEC9A","CLNK","IDO1",
              "CD1C", "FCER1A","CD1E",
              "LAMP3", "CCR7", "FSCN1",
              "IL1B", "EREG", "VCAN",
              "DNAJB1", "HSPH1", "BAG3",
              "CXCL11", "CXCL9", "CXCL10",
              "C3", "CSF1R", "CX3CR1",
              "FTL", "YBX1", "PRDX1",
              "FOS", "KLF2", "EGR1",
              "FOLR2", "MRC1", "MAF",
              "APOE", "C1QC", "C1QA",
              "MT1A", "MT1F", "SPP1",
              "CD3D", "CD2", "TRAC")

pdf("./5.4.Mono_harmony/Nomix_AllmarkerBubble.pdf",width = 10, height = 7)
jjDotPlot(Mono,
          gene = Features,
          id = "celltype",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 4,
          x.text.vjust = 0.5,
          col.min = -2,
          col.max = 2)
dev.off()

saveRDS(Mono, "./5.4.Mono_harmony/Nomix_Myeloid_annoted.Rds")

## Zhang Myeloid score ----------------------------------

Zhang_gs <- read.table("Zhang_Myeloid.txt", header = T, sep = "\t")

#cDC1-BATF3     cDC2-CD1C    Macro-IL1B   Macro-NLRP3    Macro-PLTP
#124           110           127            96           175
#Mono-CD14 Mono-CD14CD16     Mono-CD16 Monolike-FCN1    pDC-LILRA4
#49            14           120            37           153
#TAM-C1QC      TAM-SPP1
#160           101

gs <- list()
for (i in c("pDC-LILRA4","cDC1-BATF3","cDC2-CD1C","Mono-CD14","Mono-CD16",
            "Monolike-FCN1","Macro-IL1B","Macro-NLRP3","Macro-PLTP","TAM-C1QC","TAM-SPP1")) {
  gs[[i]] <- Zhang_gs[Zhang_gs$term == i,]$gene
}

library(AUCell)

expMtx <- GetAssayData(Mono, assay = "RNA", slot = "data")
cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.4.Mono_harmony/Zhang_gs/Zhang_gs_AUC.Rds")

AUCell <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

celltype <- Mono$celltype %>% as.character() %>% as.data.frame()
colnames(celltype) <- "celltype"
rownames(celltype) <- Cells(Mono)

mtx <- cbind(AUCell, celltype)

mtx <- split(mtx, mtx$celltype)
median <- c()
for (i in 1:length(mtx)) {
  median_tmp <- mtx[[i]][,-grep("celltype",colnames(mtx[[i]]))] %>% apply(., 2, median) 
  median <- rbind(median, median_tmp)
  rownames(median)[i] <- names(mtx)[i]
}

## pheatmap

pdf("./5.4.Mono_harmony/Zhang_gs/Zhang_gs_heatmap.pdf", width = 4.8, height = 4)
bk <- c(seq(-2,0,by=0.01),seq(0.01,2,by=0.01))
pheatmap(median,
         cluster_cols = F,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = T,
         show_rownames = T,
         scale="column",  #矫正
         breaks=bk,
         legend_breaks=seq(-2,2,1),
         fontsize = 8,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "",
         angle_col = "45",
         border_color = "grey",
         display_numbers = T,
         number_color = "black",
         legend = F
)
dev.off()

### SPP1 score in Macro ---------------------------------

SPP1 <- as.data.frame(AUCell$"TAM-SPP1")
colnames(SPP1) <- "SPP1_Score"
rownames(SPP1) <- rownames(AUCell)

celltype <- Mono$celltype %>% as.character() %>% as.data.frame()
colnames(celltype) <- "celltype"
rownames(celltype) <- Cells(Mono)

SPP1 <- as.data.frame(cbind(SPP1, celltype))

SPP1 <- SPP1[which(SPP1$celltype %in% c('Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                        'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL')),]
SPP1$celltype <- factor(SPP1$celltype, levels = c('Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                                  'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL'))

library(reshape2)
SPP1 <- melt(SPP1)

P1 = ggplot(SPP1, aes(x = celltype, y = value, fill = celltype)) +
  scale_fill_manual(values = color[c(7:13)]) +
  geom_boxplot(width = 0.6, outlier.colour = "grey40", position=position_dodge(0.8), outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_blank()) +
  ylab("TAM_SPP1 AUCell") +
  NoLegend()

Macro <- subset(Mono, celltype %in% c('Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                            'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL'))
Macro$celltype <- factor(Macro$celltype, levels = c('Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                                    'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL'))

P2 <- VlnPlot(Macro, "SPP1", cols = color[c(7:13)]) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_blank()) +
  ylab("SPP1 exp") +
  NoLegend()

P3 <- VlnPlot(Macro, "FN1", cols = color[c(7:13)]) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_blank()) +
  ylab("FN1 exp") +
NoLegend()

library(gridExtra)
pdf("./5.4.Mono_harmony/Zhang_gs/SPP1_anno.pdf", 
    width = 3, height = 7)
grid.arrange(grobs = list(P1, P2, P3), ncol = 1)
dev.off()


## Macrophage score ----------------------------------

Macrp_gs <- read.table("Macrophage_gs.txt", header = T, sep = "\t")

gs <- list()
for (i in c("M1","M2","MDSC","Angiogenesis","Phagocytosis")) {
  gs[[i]] <- Macrp_gs[Macrp_gs$term == i,]$gene
}

library(AUCell)

expMtx <- GetAssayData(Mono, assay = "RNA", slot = "data")
cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.4.Mono_harmony/Macro_gs/Macro_gs_AUC.Rds")

AUCell <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

celltype <- Mono$celltype %>% as.character() %>% as.data.frame()
colnames(celltype) <- "celltype"
rownames(celltype) <- Cells(Mono)

mtx <- as.data.frame(cbind(AUCell, celltype))

mtx <- mtx[which(mtx$celltype %in% c('Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                     'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL')),]

mtx <- split(mtx, mtx$celltype)
mean <- c()
for (i in 1:length(mtx)) {
  mean_tmp <- mtx[[i]][,-grep("celltype",colnames(mtx[[i]]))] %>% colMeans()
  mean <- rbind(mean, mean_tmp)
  rownames(mean)[i] <- names(mtx)[i]
}

median <- c()
for (i in 1:length(mtx)) {
  median_tmp <- mtx[[i]][,-grep("celltype",colnames(mtx[[i]]))] %>% apply(., 2, median) 
  median <- rbind(median, median_tmp)
  rownames(median)[i] <- names(mtx)[i]
}

mean <- mean[c('Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
               'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL'),]

median <- median[c('Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                   'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL'),]

## pheatmap

pdf("./5.4.Mono_harmony/Macro_gs/Macro_gs_heatmap.pdf", width = 2.5, height = 3)
bk <- c(seq(-2,0,by=0.01),seq(0.01,2,by=0.01))
pheatmap(mean,
         cluster_cols = F,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = T,
         show_rownames = T,
         scale="column",  #矫正
         breaks=bk,
         legend_breaks=seq(-2,2,1),
         fontsize = 8,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "",
         angle_col = "45",
         border_color = "grey",
         display_numbers = T,
         number_color = "black",
         legend = F
)
dev.off()

### Macro_APOE M2 and MDSC score (EO vs LO) -------------------------------------

Mono <- readRDS("5.4.Mono_harmony/Nomix_Myeloid_annoted.Rds")
cells_AUC <- readRDS( "./5.4.Mono_harmony/Macro_gs/Macro_gs_AUC.Rds")
AUCell <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

Cell_group <- data.frame(Cell = Cells(Mono), CellType = Mono$celltype, SampleType = Mono$SampleType)
rownames(Cell_group) <- NULL
Cell_group <- column_to_rownames(Cell_group, "Cell")
APOE_group <- Cell_group[Cell_group$CellType == "Macro_APOE",]

APOE_mtx <- merge(APOE_group, AUCell, by = 0) %>% column_to_rownames("Row.names")
APOE_M2_MDSC <- APOE_mtx[,c(2,4,5)]
data <- melt(APOE_M2_MDSC)

P = ggplot(data, aes(x = variable, y = value, fill = SampleType)) +
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
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("./5.4.Mono_harmony/Macro_gs/BoxPlot_APOE_M2_MDSC_EOLO.pdf", width = 3, height = 3)
P + NoLegend()
dev.off()

APOE_M2 <- APOE_M2_MDSC[,c(1,2)]
P = ggplot(APOE_M2, aes(x = SampleType, y = M2, fill = SampleType)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        axis.title.x=element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_blank()) + 
  ylab("M2 AUCell") +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("./5.4.Mono_harmony/Macro_gs/BoxPlot_APOE_M2_EOLO.pdf", width = 1, height = 1.5)
P + NoLegend()
dev.off()

APOE_MDSC <- APOE_M2_MDSC[,c(1,3)]
P = ggplot(APOE_MDSC, aes(x = SampleType, y = MDSC, fill = SampleType)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        axis.title.x=element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_blank()) + 
  ylab("MDSC AUCell") +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("./5.4.Mono_harmony/Macro_gs/BoxPlot_APOE_MDSC_EOLO.pdf", width = 1, height = 1.5)
P + NoLegend()
dev.off()

## enrichment analysis --------------------------------

library(SCP)

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)

Allmarkers <- FindAllMarkers(Mono,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = T)
saveRDS(Allmarkers, "./5.4.Mono_harmony/Markers/Nomix_celltype_markers_nomix.Rds")

Allmarkers_sig <- Allmarkers[which(Allmarkers$p_val_adj < 0.01 & 
                                     Allmarkers$pct.1 > 0.3 & 
                                     Allmarkers$avg_log2FC > 0.25), ]

Allmarkers_order <- c()
for (i in unique(Allmarkers_sig$cluster)) {
  
  tmp <- Allmarkers_sig[Allmarkers_sig$cluster == i,]
  tmp <- tmp[order(tmp$avg_log2FC, decreasing = T),]
  Allmarkers_order <- rbind(Allmarkers_order, tmp)
  
}

topmarkers <- Allmarkers %>% group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

ht <- FeatureHeatmap(
  srt = Mono, 
  group.by = "celltype", 
  features = topmarkers$gene, 
  feature_split = topmarkers$cluster,
  species = "Homo_sapiens", 
  db = c("GO_BP"), 
  anno_terms = TRUE,
  nlabel = 0,
  height = 9, 
  width = 8,
  group_palcolor = color,
  feature_split_palcolor = color,
)

saveRDS(ht, "5.4.Mono_harmony/enrichment/Heatmap_celltype_top50gene_topGOBP.Rds")

pdf("5.4.Mono_harmony/enrichment/Heatmap_celltype_top50gene_topGOBP.pdf", width = 20, height = 20)
ht$plot
dev.off()

# plot genes

topgenes <- Allmarkers %>% group_by(cluster) %>%
  top_n(n = 9, wt = avg_log2FC)
plot_gene <- data.frame(celltype = topgenes$cluster, gene = topgenes$gene)
pn <- length(unique(plot_gene$celltype))
dim1 <- 3*pn
dim2 <- 3
gene.text=t(matrix(plot_gene$gene,nrow = dim2,ncol = dim1))
rownames(gene.text)= seq(-0.5,-(dim1-0.5),-1)
colnames(gene.text)=seq(0.5,(dim2-0.5),1)
gene.text=as.data.frame(gene.text)
gene.text$dim1=rownames(gene.text)
gene.text=reshape2::melt(gene.text,id="dim1")
colnames(gene.text)[2:3]=c("dim2","gene")
gene.text$dim1=as.numeric(as.character(gene.text$dim1))
gene.text$dim2=as.numeric(as.character(gene.text$dim2))

plot <- gene.text %>% 
  ggplot(aes(x=dim2,y=dim1)) + geom_text(aes(label=gene)) +
  #geom_hline(yintercept = seq(-dim1,0,5)[-c(1,length(seq(-dim1,0,5)))],color="black",linetype=5)+
  scale_x_continuous("",expand = c(0,0),limits = c(0,10))+
  scale_y_continuous("",expand = c(0,0),limits = c(-dim1,0))+
  labs(title = paste(unique(plot_gene$program),collapse = "; "))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 20)
  )

pdf("5.4.Mono_harmony/enrichment/Gene_box.pdf", width = 7, height = 10)
plot
dev.off()

## Roe -------------------------------------------

tab <- table(Mono$celltype, Mono$SampleType)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.4.Mono_harmony/Nomix_celltype_Roe_SampleType.pdf", width = 3, height = 5)
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

pdf("./5.4.Mono_harmony/Nomix_celltype_Roe_SampleType_invert.pdf", width = 6, height = 2)
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

## Proportion -------------------------------------------

tab <- table(Mono$celltype, Mono$SampleType)
Group_type <- melt(tab)
colnames(Group_type) <- c("CellType", "Group", "nCells")

pdf("./5.4.Mono_harmony/Nomix_celltype_proportion_SampleType.pdf", height = 2, width = 14)
Group_type %>%
  ggplot(aes(x=Group, y=nCells,fill=CellType)) +
  geom_bar(stat="identity", position = 'fill') +
  theme_bw() +
  scale_fill_manual(values = color) +
  ylab("Pct. Cells") +
  xlab("") +
  NoLegend() +
  rotate()
dev.off()

# Monocle2 --------------------------------------------

#构造表达及注释数据
exp.matrix <- as(as.matrix(Mono@assays$RNA@data), "sparseMatrix")
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(exp.matrix)
exp_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- Mono@meta.data
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
diff_test_res<-differentialGeneTest(exp.monocle,fullModelFormulaStr = "~celltype") 
ordering_genes<-row.names (subset(diff_test_res, qval < 0.01))
saveRDS(ordering_genes, "./5.4.Mono_harmony/Monocle/Nomix_ordering_genes.Rds")
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
#plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

#修改monocle对象中的列名示例
saveRDS(exp.monocle, "./5.4.Mono_harmony/Monocle/Nomix_exp.monocle.Rds")

#将不同分组情况的拟时序轨迹图画到一起

plot1<-plot_cell_trajectory(exp.monocle, color_by = "celltype",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
plot2<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size = 0.2, shuffle = T) + 
  scale_color_gradient(low = "#101D44", high = "#9FB2ED")
plot3<-plot_cell_trajectory(exp.monocle, color_by = "SampleType",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
pdf("./5.4.Mono_harmony/Monocle/Nomix_trajectory_plot.pdf",width = 3,height = 9)
CombinePlots(plots = list(plot1,plot2,plot3),legend = NULL,ncol=1)
dev.off()

#split by celltype

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

celltype_name <- as.character(unique(exp.monocle$celltype))
Plist <- list()

#cluster
for (i in 1:length(celltype_name)) {
  a <- rownames(Mono@meta.data)[ Mono@meta.data$celltype == celltype_name[i] ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$celltype == celltype_name[i] ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend() + ggtitle(celltype_name[i])
  Plist[[i]] <- plot
  names(Plist)[i] <- celltype_name[i]
}

Plist[[length(Plist)+1]] <- plot2+NoLegend()+ggtitle("Pseudotime")
pdf("./5.4.Mono_harmony/Monocle/Nomix_trajectory_plot_celltypeSplit.pdf",width = 12,height = 12)
CombinePlots(plots = Plist,legend = NULL,ncol=3)
dev.off()

#split by sampletype

monocle <- exp.monocle
dimS <- exp.monocle@reducedDimS

#cluster
for (i in c("active_Smoker", "non_Smoker")) {
  a <- rownames(CD4_Tcell@meta.data)[ CD4_Tcell@meta.data$SampleType == i ]
  S <- dimS[,colnames(dimS) %in% a]
  mon <- monocle[, monocle$SampleType == i ] 
  mon@reducedDimS <- S
  plot <- plot_cell_trajectory(mon, color_by = "Pseudotime",cell_size = 0.05, shuffle = T)  + 
    NoLegend()
  assign(paste0("P_", i), plot)
}

Plist <- list(P_active_Smoker,P_non_Smoker)
pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/trajectory_plot_SampleTypeSplit.pdf",width = 6, height = 10)
CombinePlots(plots = Plist,legend = NULL,ncol=1)
dev.off()

#features plots

gene <- c("TIGIT", "FOXP3", "CTLA4", "LAG3", "HAVCR2",
          "CCR7", "SELL", "TCF7", "IL7R", "LEF1")

pdf("./5.2.Tcell_harmony/Monocle/CD4_Tcell/feature_naive_dysfunction_plot.pdf", width = 6, height = 5)
plot_genes_in_pseudotime(exp.monocle[gene,], color_by = "Pseudotime", ncol =2)
dev.off()

# slingshot ---------------------------------------

library(SCP)

Mono <- RunSlingshot(
  srt = Mono, 
  group.by = "celltype", 
  reduction = "UMAP")
dev.off()

pdf("./5.4.Mono_harmony/Slingshot/Trajectory.pdf", width = 10, height = 5)
CellDimPlot(Mono, 
            group.by = "celltype", 
            reduction = "UMAP",
            dims = c(1, 2), 
            lineages = paste0("Lineage", 1:3),
            palcolor = color)
dev.off()

pdf("./5.4.Mono_harmony/Slingshot/Lineage_seperate.pdf", width = 7, height = 2)
FeatureDimPlot(
  Mono, 
  features = paste0("Lineage", 1:3), 
  reduction = "UMAP", 
  theme_use = "theme_blank")
dev.off()

# scVelo -----------------------------------------

library(SCP)

## cell ID
Cells <- Cells(Mono)
## change ID
Cells_orig <- substr(Cells, nchar(Cells)-2, nchar(Cells))
Cells <- substr(Cells, 7, nchar(Cells)-4)
Cells <- paste0("PCSC", Cells_orig, "_RNA_10X:", Cells)

# MP + SampleType & write h5ad
Mono <- RenameCells(Mono, new.names=Cells) 
library(SCP)
adata <- srt_to_adata(Mono)
adata$write_h5ad("./5.4.Mono_harmony/scVelo/Myeloid.h5ad")

#python

c
import scvelo as scv
import pandas as pd
import numpy as np
#scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

# 读取文件

ldata = scv.read("/data_njmu/Lab_Seq/PCa_singlecell/data_tmp/loom/Loom_integration.h5ad", cache=False)
## 把AnnData里barcode中的'v3:|x'去掉
obs_name_rep = ldata.obs_names.str.replace('x','')
ldata.obs_names = obs_name_rep
obs_name_rep = ldata.obs_names.str.replace('PCSC004','PCSC005')
ldata.obs_names = obs_name_rep
ldata.var_names_make_unique()

adata = scv.read("./5.4.Mono_harmony/scVelo/Myeloid.h5ad")
adata.var_names_make_unique()

adata_m = scv.utils.merge(adata, ldata)

# MP1_OnOff剪切和未剪切的比例
scv.pl.proportions(adata_m, 'celltype')

# 筛选基因 及 在PCA空间中最近的邻居之间计算的一阶矩和二阶矩
scv.pp.filter_and_normalize(adata_m, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_m, n_pcs=30, n_neighbors=30)

# 估计RNA速率
scv.tl.velocity(adata_m)

# 细胞-细胞间的转换概率
scv.tl.velocity_graph(adata_m)
scv.pl.velocity_embedding_stream(adata_m, basis='umap', color="celltype", palette = ["#A52A2A","#053E7A","#E46B4F","#8B8B00","#B8860B","#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B","#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054"])
scv.pl.velocity_embedding(adata_m, arrow_length=8, arrow_size=3, dpi=120, color="celltype", palette = ["#A52A2A","#053E7A","#E46B4F","#8B8B00","#B8860B","#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B","#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054"])
scv.pl.velocity_embedding_grid(adata_m, arrow_length=3, arrow_size=2, dpi=120, color="celltype", palette = ["#A52A2A","#053E7A","#E46B4F","#8B8B00","#B8860B","#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B","#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054"])

scv.tl.score_genes_cell_cycle(adata_m)
scv.pl.scatter(adata_m, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])

# cellphoneDB MP1_On all -------------------------------------------------------------

## prepare data ------------------------------

Mono <- readRDS("./5.4.Mono_harmony/Nomix_Myeloid_annoted.Rds")
Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
Epi_MP1On <- subset(Epi, MP1_OnOff == "On")
Mye_Epi <- merge(Epi_MP1On, Mono)
Mye_Epi@meta.data$celltype[Mye_Epi@meta.data$celltype == "Epithelial"] <- "Epi_MP1_On"
Mye_Epi$Sample_Celltype <- paste(Mye_Epi$orig.ident, Mye_Epi$celltype, sep = "_")
Mye_Epi$Group_Celltype <- paste(Mye_Epi$SampleType, Mye_Epi$celltype, sep = "_")

expr <- as.matrix(Mye_Epi@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("./5.4.Mono_harmony/CellPhoneDB/input_Sample/expr.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
write.table(expr, paste0("./5.4.Mono_harmony/CellPhoneDB/input_Group/expr.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(Mye_Epi@meta.data), celltype = Mye_Epi@meta.data$Sample_Celltype)  
write.table(celltype, paste0("./5.4.Mono_harmony/CellPhoneDB/input_Sample/celltype.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype
celltype <- data.frame(Cell = rownames(Mye_Epi@meta.data), celltype = Mye_Epi@meta.data$Group_Celltype)  
write.table(celltype, paste0("./5.4.Mono_harmony/CellPhoneDB/input_Group/celltype.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

## cellphoneDB analysis ----------------------------------------------

cellphonedb method statistical_analysis --counts-data gene_name --output-path ./5.4.Mono_harmony/CellPhoneDB/output_Sample --threads 20 ./5.4.Mono_harmony/CellPhoneDB/input_Sample/celltype.txt ./5.4.Mono_harmony/CellPhoneDB/input_Sample/expr.txt;
cellphonedb plot heatmap_plot --pvalues-path ./5.4.Mono_harmony/CellPhoneDB/output_Sample/pvalues.txt --output-path ./5.4.Mono_harmony/CellPhoneDB/output_Sample ./5.4.Mono_harmony/CellPhoneDB/input_Sample/celltype.txt;

cellphonedb method statistical_analysis --counts-data gene_name --output-path ./5.4.Mono_harmony/CellPhoneDB/output_Group --threads 20 ./5.4.Mono_harmony/CellPhoneDB/input_Group/celltype.txt ./5.4.Mono_harmony/CellPhoneDB/input_Group/expr.txt;
cellphonedb plot heatmap_plot --pvalues-path ./5.4.Mono_harmony/CellPhoneDB/output_Group/pvalues.txt --output-path ./5.4.Mono_harmony/CellPhoneDB/output_Group ./5.4.Mono_harmony/CellPhoneDB/input_Group/celltype.txt;

## interpret - APOE MP1 EO vs LO sig counts ---------------------------------------

library(ktplots)

pvals_EO <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Group/pvalues.txt", check.names = FALSE)
means_EO <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Group/means.txt", check.names = FALSE)
decon_EO <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Group/deconvoluted.txt", check.names = FALSE)
sigmeans_EO <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Group/significant_means.txt", check.names = FALSE)

count_net_EO <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Group/count_network.txt", check.names = FALSE)
inter_net_EO <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Group/interaction_count.txt", check.names = FALSE)

library(dplyr)
library(CellChat)
library(tidyr)
count_inter <- count_net_EO
#count_inter$count <- count_inter$count/50
count_inter<-spread(count_inter, TARGET, count)
rownames(count_inter) <- count_inter$SOURCE
EOPC <- count_inter[grep("EOPC",rownames(count_inter)),grep("EOPC", colnames(count_inter))]
LOPC <- count_inter[grep("LOPC",rownames(count_inter)),grep("LOPC", colnames(count_inter))]

rownames(EOPC) <- substr(rownames(EOPC), 6, nchar(rownames(EOPC)))
colnames(EOPC) <- substr(colnames(EOPC), 6, nchar(colnames(EOPC)))
EOPC <- as.matrix(EOPC)
EOPC <- EOPC[c('pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
               'Mono_IL1B', 'Mono_FCN1',
               'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
               'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL',
               'Epi_MP1_On'),]
rownames(LOPC) <- substr(rownames(LOPC), 6, nchar(rownames(LOPC)))
colnames(LOPC) <- substr(colnames(LOPC), 6, nchar(colnames(LOPC)))
LOPC <- as.matrix(LOPC)
LOPC <- LOPC[c('pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
               'Mono_IL1B', 'Mono_FCN1',
               'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
               'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL',
               'Epi_MP1_On'),]

color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")

library(gridExtra)
pdf("./5.4.Mono_harmony/CellPhoneDB/output_Group/chord_Mye_MP1.pdf", width = 8, height = 4)
par(mfrow = c(1,2))
netVisual_circle(EOPC, 
                 targets.use = "Epi_MP1_On",
                 weight.scale = T, 
                 color.use = color,
                 edge.weight.max = 200,
                 arrow.size=0,
                 label.edge = T,
                 vertex.label.cex = 0.8,
                 edge.label.cex = 0.6,
                 title.name = "EOPC")
netVisual_circle(LOPC, 
                 targets.use = "Epi_MP1_On",
                 weight.scale = T, 
                 color.use = color,
                 edge.weight.max = 200,
                 arrow.size=0,
                 label.edge = T,
                 vertex.label.cex = 0.8,
                 edge.label.cex = 0.6,
                 title.name = "EOPC")
dev.off()

pdf("./5.4.Mono_harmony/CellPhoneDB/output_Group/chord_APOE_MP1.pdf", width = 8, height = 4)
par(mfrow = c(1,2))
netVisual_circle(EOPC, 
                 sources.use = "Macro_APOE",
                 targets.use = "Epi_MP1_On",
                 weight.scale = T, 
                 color.use = color,
                 edge.weight.max = 200,
                 arrow.size=0,
                 label.edge = T,
                 vertex.label.cex = 0.8,
                 edge.label.cex = 0.6,
                 title.name = "EOPC")
netVisual_circle(LOPC, 
                 sources.use = "Macro_APOE",
                 targets.use = "Epi_MP1_On",
                 weight.scale = T, 
                 color.use = color,
                 edge.weight.max = 200,
                 arrow.size=0,
                 label.edge = T,
                 vertex.label.cex = 0.8,
                 edge.label.cex = 0.6,
                 title.name = "EOPC")
dev.off()

## interpret - APOE MP1 interactions ---------------------------------------

library(ktplots)

pvals <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Sample/pvalues.txt", check.names = FALSE)
means <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Sample/means.txt", check.names = FALSE)
decon <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Sample/deconvoluted.txt", check.names = FALSE)
sigmeans <- read.delim("./5.4.Mono_harmony/CellPhoneDB/output_Sample/significant_means.txt", check.names = FALSE)

#APOE to other - EO sig

pvals_APOE_to_other <- pvals[, grep(".*Macro_APOE\\|.*Epi_MP1_On$", colnames(pvals))]
pvals_APOE_to_other <- pvals_APOE_to_other[, which(substr(colnames(pvals_APOE_to_other),1,5) == substr(colnames(pvals_APOE_to_other),18,22))]
pvals_APOE_to_other <- cbind(pvals[,1:11], pvals_APOE_to_other)
pvals_APOE_to_other <- pvals_APOE_to_other[which(pvals_APOE_to_other$'EOPC1_Macro_APOE|EOPC1_Epi_MP1_On' < 0.05 &
                                                   pvals_APOE_to_other$'EOPC2_Macro_APOE|EOPC2_Epi_MP1_On' < 0.05 &
                                                   pvals_APOE_to_other$'EOPC3_Macro_APOE|EOPC3_Epi_MP1_On' < 0.05 &
                                                   pvals_APOE_to_other$'EOPC4_Macro_APOE|EOPC4_Epi_MP1_On' < 0.05),]

means_APOE_to_other <- means[, grep(".*Macro_APOE\\|.*Epi_MP1_On$", colnames(means))]
means_APOE_to_other <- cbind(pvals[,1:11], means_APOE_to_other)

pdf("./5.4.Mono_harmony/CellPhoneDB/output_Sample/APOE_MP1_interaction.pdf", width = 7, height = 12)
plot_cpdb(cell_type1 = 'Macro_APOE', cell_type2 = "Epi_MP1_On", scdata = Mye_Epi,
          idents = 'celltype', split.by = "orig.ident", means = means_APOE_to_other, 
          pvals = pvals_APOE_to_other, highlight = "#A52A2A", keep_significant_only = T,
          highlight_size = NULL) +
  theme(axis.text  = element_text(size = 10, color = 'black'))
dev.off()

# seperate analysis ---------------------------------------

for (i in c("EOPC1", "EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  pvals <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/pvalues.txt"), check.names = FALSE)
  assign(paste0("pvals_", i), pvals)
  means <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/means.txt"), check.names = FALSE)
  assign(paste0("means_", i), means)
  decon <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/deconvoluted.txt"), check.names = FALSE)
  assign(paste0("decon_", i), decon)
  mynet <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/count_network.txt"), check.names = FALSE)
  assign(paste0("mynet_", i), mynet)
  sig <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/significant_means.txt"), check.names = FALSE)
  assign(paste0("sig_", i), sig)
}

### MP1_On and Macro_APOE ---------------------------------------------------------------

MtoC_pval <- data.frame(pair = pvals_EOPC1$interacting_pair, EOPC1 = pvals_EOPC1$'Macro_APOE|Epi_MP1_On')
MtoC_mean <- data.frame(pair = means_EOPC1$interacting_pair, EOPC1 = means_EOPC1$'Macro_APOE|Epi_MP1_On')
MtoC_mean <- MtoC_mean[match(MtoC_pval$pair,MtoC_mean$pair),]

for (i in c("EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  
  pvals <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/pvalues.txt"), check.names = FALSE)
  means <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/means.txt"), check.names = FALSE)
  
  MtoC_pval_tmp <- data.frame(pair = pvals$interacting_pair, pvals = pvals$'Macro_APOE|Epi_MP1_On')
  MtoC_means_tmp <- data.frame(pair = pvals$interacting_pair, means = means$'Macro_APOE|Epi_MP1_On')
  MtoC_means_tmp <- MtoC_means_tmp[match(MtoC_pval_tmp$pair,MtoC_means_tmp$pair),]
  
  colnames(MtoC_pval_tmp)[2] <- i
  colnames(MtoC_means_tmp)[2] <- i
  
  MtoC_pval <- merge(MtoC_pval, MtoC_pval_tmp, by = "pair")
  MtoC_mean <- merge(MtoC_mean, MtoC_means_tmp, by = "pair")
}

for (i in c("LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  
  pvals <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/pvalues.txt"), check.names = FALSE)
  means <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/means.txt"), check.names = FALSE)
  
  MtoC_pval_tmp <- data.frame(pair = pvals$interacting_pair, pvals = pvals$'Macro_APOE|Epi_MP1_On')
  MtoC_means_tmp <- data.frame(pair = pvals$interacting_pair, means = means$'Macro_APOE|Epi_MP1_On')
  MtoC_means_tmp <- MtoC_means_tmp[match(MtoC_pval_tmp$pair,MtoC_means_tmp$pair),]
  
  colnames(MtoC_pval_tmp)[2] <- i
  colnames(MtoC_means_tmp)[2] <- i
  
  MtoC_pval <- merge(MtoC_pval, MtoC_pval_tmp, by = "pair", all = T)
  MtoC_mean <- merge(MtoC_mean, MtoC_means_tmp, by = "pair", all = T)
  
}

#### shared in EOPC ----------------------------------

# APOE_to_MP1

EOPC_Share_pval <- MtoC_pval[which(MtoC_pval$EOPC1 < 0.05 & MtoC_pval$EOPC2 < 0.05 & MtoC_pval$EOPC3 < 0.05 & MtoC_pval$EOPC4 < 0.05),]
EOPC_Share_mean <- MtoC_mean[match(EOPC_Share_pval$pair,MtoC_mean$pair),]

data_EOPC_Share_pval <- melt(EOPC_Share_pval)
colnames(data_EOPC_Share_pval)[3] <- "pvals"
data_EOPC_Share_mean <- melt(EOPC_Share_mean)
colnames(data_EOPC_Share_mean)[3] <- "means"

data <- merge(data_EOPC_Share_pval, data_EOPC_Share_mean, by = c("pair", "variable"))
data$pvals[data$pvals == 0] <- 0.001
data$pvals <- -log10(data$pvals)

P <- ggplot(data, aes(x=variable, y=pair, size=means)) + 
  geom_point(aes(colour = pvals)) +
  theme_bw()+theme(axis.text.x = element_text(angle=90,size=10)) +
  scale_x_discrete(position = "bottom") +
  scale_colour_distiller(type = "seq", palette = "PuBu", direction = 1) +
  theme(panel.grid.major = element_blank(),
        legend.key.size = unit(7,'pt'),
        legend.key.height = unit(7,'pt'),
        legend.key.width = unit(7,'pt'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y = element_text(size=7, colour="black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_blank())
ggsave('./5.4.Mono_harmony/CellPhoneDB/output_Sample/APOE_MP1_interaction.pdf', P, width = 5, height = 5)

data_select <- data[data$pair %in% c("FN1_integrin_aVb5_complex", "FN1_integrin_aVb1_complex",
                                     "FN1_integrin_a3b1_complex", "FN1_integrin_a2b1_complex",
                                     "EREG_EGFR"),]
P <- ggplot(data_select, aes(x=variable, y=pair, size=means)) + 
  geom_point(aes(colour = pvals)) +
  scale_size_continuous(range=c(1,4)) +
  theme_bw() +
  scale_x_discrete(position = "bottom") +
  scale_colour_distiller(type = "seq", palette = "PuBu", direction = 1) +
  theme(panel.grid.major = element_blank(),
        legend.key.size = unit(7,'pt'),
        legend.key.height = unit(7,'pt'),
        legend.key.width = unit(7,'pt'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y = element_text(size=7, colour="black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_blank())
ggsave('./5.4.Mono_harmony/CellPhoneDB/output_Sample/APOE_MP1_interaction_selected.pdf', P, width = 4, height = 1.2)

# cellphoneDB MP1_On-------------------------------------------------------------

## prepare data ------------------------------------------------

Mono <- readRDS("./5.4.Mono_harmony/Nomix_Myeloid_annoted.Rds")
Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
Epi_MP1On <- subset(Epi, MP1_OnOff == "On")
Mye_Epi <- merge(Epi_MP1On, Mono)
Mye_Epi@meta.data$celltype[Mye_Epi@meta.data$celltype == "Epithelial"] <- "Epi_MP1_On"

for (i in unique(Mye_Epi$SampleType)) {
  Obj <- subset(Mye_Epi, SampleType == i)
  expr <- as.matrix(Obj@assays$RNA@counts)
  expr <- data.frame(Gene = rownames(expr), expr)
  write.table(expr, paste0("./5.4.Mono_harmony/CellPhoneDB/input/expr_", i, ".txt"),
              sep='\t', quote=F, row.names = F)  #表达谱
  celltype <- data.frame(Cell = rownames(Obj@meta.data), celltype = Obj@meta.data$celltype)  
  write.table(celltype, paste0("./5.4.Mono_harmony/CellPhoneDB/input/celltype_", i, ".txt"), 
              sep='\t', quote=F, row.names = F)  #celltype
}

for (i in unique(Mye_Epi$orig.ident)) {
  Obj <- subset(Mye_Epi, orig.ident == i)
  expr <- as.matrix(Obj@assays$RNA@counts)
  expr <- data.frame(Gene = rownames(expr), expr)
  write.table(expr, paste0("./5.4.Mono_harmony/CellPhoneDB/input/expr_", i, ".txt"),
              sep='\t', quote=F, row.names = F)  #表达谱
  celltype <- data.frame(Cell = rownames(Obj@meta.data), celltype = Obj@meta.data$celltype)  
  write.table(celltype, paste0("./5.4.Mono_harmony/CellPhoneDB/input/celltype_", i, ".txt"), 
              sep='\t', quote=F, row.names = F)  #celltype
}

## cellphoneDB analysis ----------------------------------------------

for i in EOPC LOPC;do
mkdir ./5.4.Mono_harmony/CellPhoneDB/output_$i;
done

for i in EOPC LOPC;do
cellphonedb method statistical_analysis --counts-data gene_name --output-path ./5.4.Mono_harmony/CellPhoneDB//output_$i --threads 20 ./5.4.Mono_harmony/CellPhoneDB//input/celltype_$i.txt ./5.4.Mono_harmony/CellPhoneDB//input/expr_$i.txt;
done

for i in EOPC LOPC;do
cellphonedb plot dot_plot --means-path ./5.4.Mono_harmony/CellPhoneDB//output_$i/means.txt --pvalues-path ./5.4.Mono_harmony/CellPhoneDB//output_$i/pvalues.txt --output-path ./5.4.Mono_harmony/CellPhoneDB//output_$i;
cellphonedb plot heatmap_plot --pvalues-path ./5.4.Mono_harmony/CellPhoneDB//output_$i/pvalues.txt --output-path ./5.4.Mono_harmony/CellPhoneDB//output_$i ./5.4.Mono_harmony/CellPhoneDB//input/celltype_$i.txt;
done

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
mkdir ./5.4.Mono_harmony/CellPhoneDB/output_$i;
done

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
cellphonedb method statistical_analysis --counts-data gene_name --output-path ./5.4.Mono_harmony/CellPhoneDB//output_$i --threads 20 ./5.4.Mono_harmony/CellPhoneDB//input/celltype_$i.txt ./5.4.Mono_harmony/CellPhoneDB//input/expr_$i.txt;
done

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
cellphonedb plot dot_plot --means-path ./5.4.Mono_harmony/CellPhoneDB//output_$i/means.txt --pvalues-path ./5.4.Mono_harmony/CellPhoneDB//output_$i/pvalues.txt --output-path ./5.4.Mono_harmony/CellPhoneDB//output_$i;
cellphonedb plot heatmap_plot --pvalues-path ./5.4.Mono_harmony/CellPhoneDB//output_$i/pvalues.txt --output-path ./5.4.Mono_harmony/CellPhoneDB//output_$i ./5.4.Mono_harmony/CellPhoneDB//input/celltype_$i.txt;
done

## plotting -----------------------------------------------------

library(ktplots)

for (i in c("EOPC", "LOPC")) {
  pvals <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/pvalues.txt"), check.names = FALSE)
  assign(paste0("pvals_", i), pvals)
  means <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/means.txt"), check.names = FALSE)
  assign(paste0("means_", i), means)
  decon <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/deconvoluted.txt"), check.names = FALSE)
  assign(paste0("decon_", i), decon)
  mynet <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/count_network.txt"), check.names = FALSE)
  assign(paste0("mynet_", i), mynet)
}

for (i in c("EOPC1", "EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  pvals <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/pvalues.txt"), check.names = FALSE)
  assign(paste0("pvals_", i), pvals)
  means <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/means.txt"), check.names = FALSE)
  assign(paste0("means_", i), means)
  decon <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/deconvoluted.txt"), check.names = FALSE)
  assign(paste0("decon_", i), decon)
  mynet <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/count_network.txt"), check.names = FALSE)
  assign(paste0("mynet_", i), mynet)
}

### overview - EO vs LO counts -----------------------------------------------------

counts_all <- data.frame(counts = mynet_EOPC1[c(2:14),3])
rownames(counts_all) <- mynet_EOPC1[c(2:14),2]
for (i in list(mynet_EOPC2, mynet_EOPC3, mynet_EOPC4, mynet_LOPC1, mynet_LOPC2, mynet_LOPC3, mynet_LOPC4, mynet_LOPC5, mynet_LOPC6)) {
  n = grep("Epi_MP1_On", i[,1])[length(grep("Epi_MP1_On", i[,1]))]
  counts <- data.frame(counts = i[c(2:n),3])
  rownames(counts) <- i[c(2:n),2]
  counts_all <- merge(counts_all, counts, by=0, all = T) %>% column_to_rownames("Row.names")
}
colnames(counts_all) <- c("EOPC1", "EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")
counts_all[is.na(counts_all)] <- 0
counts_all <- t(counts_all) %>% as.data.frame()
counts_all$Group <- substr(rownames(counts_all), 1, 4)

data_plot <- melt(counts_all)
data_plot$variable <- factor(data_plot$variable, levels = c('pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
                                                            'Mono_IL1B', 'Mono_FCN1',
                                                            'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                                            'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL'))

P = ggplot(data_plot, aes(x = variable, y = value, fill = Group)) +
  scale_fill_manual(values = c("#053E7ACC", "#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = "grey40", position=position_dodge(0.8), outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("MP1_On Cancer & Myeloid Cell interactions") +
  ylab("Interaction Pair Counts") +
  xlab("") + 
  stat_compare_means(method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif",
                     label.x.npc = "center")

pdf("./5.4.Mono_harmony/CellPhoneDB/Epi_Mye_interaction_pairs.pdf", 
    width = 5, height = 4)
P
dev.off()

### MP1_On and Macro_APOE ---------------------------------------------------------------

result <- plot_cpdb("Epi_MP1_On", "Macro_APOE", 
                    subset(Mye_Epi, orig.ident == "EOPC1"), 
                    "celltype", 
                    means_EOPC1, 
                    pvals_EOPC1, 
                    cluster_rows = FALSE, scale = F,
                    keep_significant_only = T, return_table = T)

CtoM_pval <- data.frame(pair = result[result$Var2 == "Epi_MP1_On-Macro_APOE",][,1], EOPC1 = result[result$Var2 == "Epi_MP1_On-Macro_APOE",][,4])
CtoM_pval <- CtoM_pval[-which(is.na(CtoM_pval[,2])),]
rownames(CtoM_pval) <- NULL
CtoM_pval <- column_to_rownames(CtoM_pval, "pair")

MtoC_pval <- data.frame(pair = result[result$Var2 == "Macro_APOE-Epi_MP1_On",][,1], EOPC1 = result[result$Var2 == "Macro_APOE-Epi_MP1_On",][,4])
MtoC_pval <- MtoC_pval[-which(is.na(MtoC_pval[,2])),]
rownames(MtoC_pval) <- NULL
MtoC_pval <- column_to_rownames(MtoC_pval, "pair")

for (i in c("EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  
  result <- plot_cpdb("Epi_MP1_On", "Macro_APOE", 
                      subset(Mye_Epi, orig.ident == i), 
                      "celltype", 
                      read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/means.txt"), check.names = FALSE),
                      read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_", i, "/pvalues.txt"), check.names = FALSE),
                      cluster_rows = FALSE, scale = F,
                      keep_significant_only = T, return_table = T)
  
  CtoM_pval_tmp <- data.frame(pair = result[result$Var2 == "Epi_MP1_On-Macro_APOE",][,1], pvals = result[result$Var2 == "Epi_MP1_On-Macro_APOE",][,4])
  colnames(CtoM_pval_tmp)[2] <- i
  CtoM_pval_tmp <- CtoM_pval_tmp[-which(is.na(CtoM_pval_tmp[,2])),]
  rownames(CtoM_pval_tmp) <- NULL
  CtoM_pval_tmp <- column_to_rownames(CtoM_pval_tmp, "pair")
  CtoM_pval <- merge(CtoM_pval, CtoM_pval_tmp, by = 0, all = T) %>% column_to_rownames("Row.names")
  
  MtoC_pval_tmp <- data.frame(pair = result[result$Var2 == "Macro_APOE-Epi_MP1_On",][,1], pvals = result[result$Var2 == "Macro_APOE-Epi_MP1_On",][,4])
  colnames(MtoC_pval_tmp)[2] <- i
  MtoC_pval_tmp <- MtoC_pval_tmp[-which(is.na(MtoC_pval_tmp[,2])),]
  rownames(MtoC_pval_tmp) <- NULL
  MtoC_pval_tmp <- column_to_rownames(MtoC_pval_tmp, "pair")
  MtoC_pval <- merge(MtoC_pval, MtoC_pval_tmp, by = 0, all = T) %>% column_to_rownames("Row.names")
  
}

# shared in EOPC

# MP1_to_Mye

EOPC_Share <- CtoM_pval[which(CtoM_pval$EOPC1 != "NA" & CtoM_pval$EOPC2 != "NA" & CtoM_pval$EOPC3 != "NA" & CtoM_pval$EOPC4 != "NA"),]
EOPC_Share_sig <- EOPC_Share

ss <- EOPC_Share_sig < 0.01
EOPC_Share_sig[ss] <- '**'
s <- EOPC_Share_sig >= 0.01 & EOPC_Share_sig < 0.05
EOPC_Share_sig[s] <- '*'
ns <- is.na(EOPC_Share_sig)
EOPC_Share_sig[ns] <- 'n.s'

ss <- EOPC_Share < 0.01
EOPC_Share[ss] <- 3
s <- EOPC_Share >= 0.01 & EOPC_Share < 0.05
EOPC_Share[s] <- 1.5
ns <- is.na(EOPC_Share)
EOPC_Share[ns] <- 0

EOPC_Share <- as.matrix(EOPC_Share)
EOPC_Share_sig <- as.matrix(EOPC_Share_sig)

pdf("./5.4.Mono_harmony/CellPhoneDB/Heatmap_MP1_to_Mye_sig.pdf", width = 6, height = 6)
pheatmap(EOPC_Share,
         scale = "none",
         cluster_row = F, 
         cluster_col = F, 
         border=NA,
         display_numbers = EOPC_Share_sig,
         fontsize_number = 6, 
         number_color = "black",
         cellwidth = 10, 
         cellheight = 8,
         color=c(colorRampPalette(colors = c("white","#E46B4F"))(50)),
         border_color="grey30",
         legend = F,
         fontsize = 8
)
dev.off()



# MP1_to_Mye

MtoC_pval[which(MtoC_pval$EOPC1 != "NA" & MtoC_pval$EOPC2 != "NA" & MtoC_pval$EOPC3 != "NA" & MtoC_pval$EOPC4 != "NA"),]
EOPC_Share <- MtoC_pval[which(MtoC_pval$EOPC1 != "NA" & MtoC_pval$EOPC2 != "NA" & MtoC_pval$EOPC3 != "NA" & MtoC_pval$EOPC4 != "NA"),]
EOPC_Share_sig <- EOPC_Share

ss <- EOPC_Share_sig < 0.01
EOPC_Share_sig[ss] <- '**'
s <- EOPC_Share_sig >= 0.01 & EOPC_Share_sig < 0.05
EOPC_Share_sig[s] <- '*'
ns <- is.na(EOPC_Share_sig)
EOPC_Share_sig[ns] <- 'n.s'

ss <- EOPC_Share < 0.01
EOPC_Share[ss] <- 3
s <- EOPC_Share >= 0.01 & EOPC_Share < 0.05
EOPC_Share[s] <- 1.5
ns <- is.na(EOPC_Share)
EOPC_Share[ns] <- 0

EOPC_Share <- as.matrix(EOPC_Share)
EOPC_Share_sig <- as.matrix(EOPC_Share_sig)

pdf("./5.4.Mono_harmony/CellPhoneDB/Heatmap_Mye_to_MP1_sig.pdf", width = 6, height = 6)
pheatmap(EOPC_Share,
         scale = "none",
         cluster_row = F, 
         cluster_col = F, 
         border=NA,
         display_numbers = EOPC_Share_sig,
         fontsize_number = 6, 
         number_color = "black",
         cellwidth = 10, 
         cellheight = 8,
         color=c(colorRampPalette(colors = c("white","#E46B4F"))(50)),
         border_color="grey30",
         legend = F,
         fontsize = 8
)
dev.off()


# cellphoneDB MP1_High-------------------------------------------------------------

## prepare data ------------------------------------------------

Mono <- readRDS("./5.4.Mono_harmony/Nomix_Myeloid_annoted.Rds")
Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
Epi_MP1H <- subset(Epi, MP1_HL == "High")
Mye_Epi <- merge(Epi_MP1H, Mono)
Mye_Epi@meta.data$celltype[Mye_Epi@meta.data$celltype == "Epithelial"] <- "Epi_MP1_High"

for (i in unique(Mye_Epi$orig.ident)) {
  Obj <- subset(Mye_Epi, orig.ident == i)
  expr <- as.matrix(Obj@assays$RNA@counts)
  expr <- data.frame(Gene = rownames(expr), expr)
  write.table(expr, paste0("./5.4.Mono_harmony/CellPhoneDB/input_HL/expr_", i, ".txt"),
              sep='\t', quote=F, row.names = F)  #表达谱
  celltype <- data.frame(Cell = rownames(Obj@meta.data), celltype = Obj@meta.data$celltype)  
  write.table(celltype, paste0("./5.4.Mono_harmony/CellPhoneDB/input_HL/celltype_", i, ".txt"), 
              sep='\t', quote=F, row.names = F)  #celltype
}

## cellphoneDB analysis ----------------------------------------------

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
mkdir ./5.4.Mono_harmony/CellPhoneDB/output_HL_$i;
done

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
cellphonedb method statistical_analysis --counts-data gene_name --output-path ./5.4.Mono_harmony/CellPhoneDB//output_HL_$i --threads 20 ./5.4.Mono_harmony/CellPhoneDB//input_HL/celltype_$i.txt ./5.4.Mono_harmony/CellPhoneDB//input_HL/expr_$i.txt;
done

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
cellphonedb plot dot_plot --means-path ./5.4.Mono_harmony/CellPhoneDB//output_HL_$i/means.txt --pvalues-path ./5.4.Mono_harmony/CellPhoneDB//output_HL_$i/pvalues.txt --output-path ./5.4.Mono_harmony/CellPhoneDB//output_HL_$i;
cellphonedb plot heatmap_plot --pvalues-path ./5.4.Mono_harmony/CellPhoneDB//output_HL_$i/pvalues.txt --output-path ./5.4.Mono_harmony/CellPhoneDB//output_HL_$i ./5.4.Mono_harmony/CellPhoneDB//input_HL/celltype_$i.txt;
done

## plotting -----------------------------------------------------

library(ktplots)

for (i in c("EOPC1", "EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  pvals <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_HL_", i, "/pvalues.txt"), check.names = FALSE)
  assign(paste0("pvals_", i), pvals)
  means <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_HL_", i, "/means.txt"), check.names = FALSE)
  assign(paste0("means_", i), means)
  decon <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_HL_", i, "/deconvoluted.txt"), check.names = FALSE)
  assign(paste0("decon_", i), decon)
  mynet <- read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_HL_", i, "/count_network.txt"), check.names = FALSE)
  assign(paste0("mynet_", i), mynet)
}

### overview - EO vs LO counts -----------------------------------------------------

counts_all <- data.frame(counts = mynet_EOPC1[c(2:14),3])
rownames(counts_all) <- mynet_EOPC1[c(2:14),2]
for (i in list(mynet_EOPC2, mynet_EOPC3, mynet_EOPC4, mynet_LOPC1, mynet_LOPC2, mynet_LOPC3, mynet_LOPC4, mynet_LOPC5, mynet_LOPC6)) {
  n = grep("Epi_MP1_High", i[,1])[length(grep("Epi_MP1_High", i[,1]))]
  counts <- data.frame(counts = i[c(2:n),3])
  rownames(counts) <- i[c(2:n),2]
  counts_all <- merge(counts_all, counts, by=0, all = T) %>% column_to_rownames("Row.names")
}
colnames(counts_all) <- c("EOPC1", "EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")
counts_all[is.na(counts_all)] <- 0
counts_all <- t(counts_all) %>% as.data.frame()
counts_all$Group <- substr(rownames(counts_all), 1, 4)

data_plot <- melt(counts_all)
data_plot$variable <- factor(data_plot$variable, levels = c('pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
                                                            'Mono_IL1B', 'Mono_FCN1',
                                                            'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                                            'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL'))

P = ggplot(data_plot, aes(x = variable, y = value, fill = Group)) +
  scale_fill_manual(values = c("#053E7ACC", "#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = "grey40", position=position_dodge(0.8), outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("MP1_On Cancer & Myeloid Cell interactions") +
  ylab("Interaction Pair Counts") +
  xlab("") + 
  stat_compare_means(method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif",
                     label.x.npc = "center")

pdf("./5.4.Mono_harmony/CellPhoneDB/Epi_HL_Mye_interaction_pairs.pdf", 
    width = 5, height = 4)
P
dev.off()

### MP1_On and Macro_APOE ---------------------------------------------------------------

result <- plot_cpdb("Epi_MP1_High", "Macro_APOE", 
                    subset(Mye_Epi, orig.ident == "EOPC1"), 
                    "celltype", 
                    means_EOPC1, 
                    pvals_EOPC1, 
                    cluster_rows = FALSE, scale = F,
                    keep_significant_only = T, return_table = T)

CtoM_pval <- data.frame(pair = result[result$Var2 == "Epi_MP1_High-Macro_APOE",][,1], EOPC1 = result[result$Var2 == "Epi_MP1_High-Macro_APOE",][,4])
CtoM_pval <- CtoM_pval[-which(is.na(CtoM_pval[,2])),]
rownames(CtoM_pval) <- NULL
CtoM_pval <- column_to_rownames(CtoM_pval, "pair")

MtoC_pval <- data.frame(pair = result[result$Var2 == "Macro_APOE-Epi_MP1_High",][,1], EOPC1 = result[result$Var2 == "Macro_APOE-Epi_MP1_High",][,4])
MtoC_pval <- MtoC_pval[-which(is.na(MtoC_pval[,2])),]
rownames(MtoC_pval) <- NULL
MtoC_pval <- column_to_rownames(MtoC_pval, "pair")

for (i in c("EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  
  result <- plot_cpdb("Epi_MP1_High", "Macro_APOE", 
                      subset(Mye_Epi, orig.ident == i), 
                      "celltype", 
                      read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_HL_", i, "/means.txt"), check.names = FALSE),
                      read.delim(paste0("5.4.Mono_harmony/CellPhoneDB/output_HL_", i, "/pvalues.txt"), check.names = FALSE),
                      cluster_rows = FALSE, scale = F,
                      keep_significant_only = T, return_table = T)
  
  CtoM_pval_tmp <- data.frame(pair = result[result$Var2 == "Epi_MP1_High-Macro_APOE",][,1], pvals = result[result$Var2 == "Epi_MP1_High-Macro_APOE",][,4])
  colnames(CtoM_pval_tmp)[2] <- i
  CtoM_pval_tmp <- CtoM_pval_tmp[-which(is.na(CtoM_pval_tmp[,2])),]
  rownames(CtoM_pval_tmp) <- NULL
  CtoM_pval_tmp <- column_to_rownames(CtoM_pval_tmp, "pair")
  CtoM_pval <- merge(CtoM_pval, CtoM_pval_tmp, by = 0, all = T) %>% column_to_rownames("Row.names")
  
  MtoC_pval_tmp <- data.frame(pair = result[result$Var2 == "Macro_APOE-Epi_MP1_High",][,1], pvals = result[result$Var2 == "Macro_APOE-Epi_MP1_High",][,4])
  colnames(MtoC_pval_tmp)[2] <- i
  MtoC_pval_tmp <- MtoC_pval_tmp[-which(is.na(MtoC_pval_tmp[,2])),]
  rownames(MtoC_pval_tmp) <- NULL
  MtoC_pval_tmp <- column_to_rownames(MtoC_pval_tmp, "pair")
  MtoC_pval <- merge(MtoC_pval, MtoC_pval_tmp, by = 0, all = T) %>% column_to_rownames("Row.names")
  
}

# shared in EOPC

# MP1_to_Mye

EOPC_Share <- CtoM_pval[which(CtoM_pval$EOPC1 != "NA" & CtoM_pval$EOPC2 != "NA" & CtoM_pval$EOPC3 != "NA" & CtoM_pval$EOPC4 != "NA"),]
EOPC_Share_sig <- EOPC_Share

ss <- EOPC_Share_sig < 0.01
EOPC_Share_sig[ss] <- '**'
s <- EOPC_Share_sig >= 0.01 & EOPC_Share_sig < 0.05
EOPC_Share_sig[s] <- '*'
ns <- is.na(EOPC_Share_sig)
EOPC_Share_sig[ns] <- 'n.s'

ss <- EOPC_Share < 0.01
EOPC_Share[ss] <- 3
s <- EOPC_Share >= 0.01 & EOPC_Share < 0.05
EOPC_Share[s] <- 1.5
ns <- is.na(EOPC_Share)
EOPC_Share[ns] <- 0

EOPC_Share <- as.matrix(EOPC_Share)
EOPC_Share_sig <- as.matrix(EOPC_Share_sig)

pdf("./5.4.Mono_harmony/CellPhoneDB/Heatmap_MP1_High_to_Mye_sig.pdf", width = 6, height = 6)
pheatmap(EOPC_Share,
         scale = "none",
         cluster_row = F, 
         cluster_col = F, 
         border=NA,
         display_numbers = EOPC_Share_sig,
         fontsize_number = 6, 
         number_color = "black",
         cellwidth = 10, 
         cellheight = 8,
         color=c(colorRampPalette(colors = c("white","#E46B4F"))(50)),
         border_color="grey30",
         legend = F,
         fontsize = 8
)
dev.off()



# MP1_to_Mye

EOPC_Share <- MtoC_pval[which(MtoC_pval$EOPC1 != "NA" & MtoC_pval$EOPC2 != "NA" & MtoC_pval$EOPC3 != "NA" & MtoC_pval$EOPC4 != "NA"),]
EOPC_Share_sig <- EOPC_Share

ss <- EOPC_Share_sig < 0.01
EOPC_Share_sig[ss] <- '**'
s <- EOPC_Share_sig >= 0.01 & EOPC_Share_sig < 0.05
EOPC_Share_sig[s] <- '*'
ns <- is.na(EOPC_Share_sig)
EOPC_Share_sig[ns] <- 'n.s'

ss <- EOPC_Share < 0.01
EOPC_Share[ss] <- 3
s <- EOPC_Share >= 0.01 & EOPC_Share < 0.05
EOPC_Share[s] <- 1.5
ns <- is.na(EOPC_Share)
EOPC_Share[ns] <- 0

EOPC_Share <- as.matrix(EOPC_Share)
EOPC_Share_sig <- as.matrix(EOPC_Share_sig)

pdf("./5.4.Mono_harmony/CellPhoneDB/Heatmap_Mye_High_to_MP1_sig.pdf", width = 6, height = 6)
pheatmap(EOPC_Share,
         scale = "none",
         cluster_row = F, 
         cluster_col = F, 
         border=NA,
         display_numbers = EOPC_Share_sig,
         fontsize_number = 6, 
         number_color = "black",
         cellwidth = 10, 
         cellheight = 8,
         color=c(colorRampPalette(colors = c("white","#E46B4F"))(50)),
         border_color="grey30",
         legend = F,
         fontsize = 8
)
dev.off()

# FA and Cholesterol metablism ----------------------------------

KEGG <- readRDS("kegg_hsa_gmt.Rds")

gs_selected <- KEGG[KEGG$term %in% c("Fatty acid biosynthesis",
                                     "Fatty acid degradation",
                                     "Biosynthesis of unsaturated fatty acids",
                                     "Fatty acid metabolism",
                                     "Cholesterol metabolism"),]
gs <- list()
for (i in c("Fatty acid biosynthesis",
            "Fatty acid degradation",
            "Biosynthesis of unsaturated fatty acids",
            "Fatty acid metabolism",
            "Cholesterol metabolism")) {
  gs[[i]] <- gs_selected[which(gs_selected$term == i),]$gene
}

library(GSVA)

expMtx <- GetAssayData(Mono, assay = "RNA", slot = "data")
GSVA <- gsva(expMtx, gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "./5.4.Mono_harmony/Lipid/5gs_GSVA.Rds")

#MP1 OnOff
celltype <- Mono$celltype %>% as.data.frame()
colnames(celltype) <- "celltype"

GSVA <- GSVA %>% as.data.frame() %>% t()
mtx <- cbind(GSVA, celltype)
mtx$celltype <- factor(mtx$celltype, levels = c('pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
                                                'Mono_IL1B', 'Mono_FCN1',
                                                'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                                'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL'))
mtx <- mtx[order(mtx$celltype),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = celltype)) +
  scale_fill_manual(values = paste0(color, "CC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 10),
        axis.text.y=element_text(size=10, colour="black"), 
        axis.title.y=element_text(size = 12),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("GSVA score") +
  xlab("")
pdf("./5.4.Mono_harmony/Lipid/BoxPlot_5gs_GSVA.pdf", width = 8, height = 4)
P
dev.off()

# CellphoneDB -----------------------------------------------

## Macro_APOE and Tcells (EOPC samples) ---------------------------------------

PCSC <- readRDS("PCSC_nomix_new_SCT.RDS")
APOE_T <- subset(PCSC, celltype %in% c("Macro_APOE", "Tcm", "CD4_Tn", "Treg",
                                       "CD8_Tn", "CD8_Tem", "NK", "NKT"))
APOE_T_EO <- subset(APOE_T, SampleType == "EOPC")

expr <- as.matrix(APOE_T_EO@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("5.4.Mono_harmony/cpdb_EO/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(APOE_T_EO@meta.data), celltype = APOE_T_EO@meta.data$celltype)  
write.table(celltype, paste0("5.4.Mono_harmony/cpdb_EO/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

cellphonedb method statistical_analysis --counts-data gene_name --output-path 5.4.Mono_harmony/cpdb_EO/output_all --threads 20 5.4.Mono_harmony/cpdb_EO/input/celltype_all.txt 5.4.Mono_harmony/cpdb_EO/input/expr_all.txt
cellphonedb plot dot_plot --means-path 5.4.Mono_harmony/cpdb_EO/output_all/means.txt --pvalues-path 5.4.Mono_harmony/cpdb_EO/output_all/pvalues.txt --output-path 5.4.Mono_harmony/cpdb_EO/output_all
cellphonedb plot heatmap_plot --pvalues-path 5.4.Mono_harmony/cpdb_EO/output_all/pvalues.txt --output-path 5.4.Mono_harmony/cpdb_EO/output_all 5.4.Mono_harmony/cpdb_EO/input/celltype_all.txt

## Macro_APOE and Tcells (LOPC samples) ---------------------------------------

PCSC <- readRDS("PCSC_nomix_new_SCT.RDS")
APOE_T <- subset(PCSC, celltype %in% c("Macro_APOE", "Tcm", "CD4_Tn", "Treg",
                                       "CD8_Tn", "CD8_Tem", "NK", "NKT"))
APOE_T_LO <- subset(APOE_T, SampleType == "LOPC")

expr <- as.matrix(APOE_T_LO@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("5.4.Mono_harmony/cpdb_LO/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(APOE_T_LO@meta.data), celltype = APOE_T_LO@meta.data$celltype)  
write.table(celltype, paste0("5.4.Mono_harmony/cpdb_LO/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

cellphonedb method statistical_analysis --counts-data gene_name --output-path 5.4.Mono_harmony/cpdb_LO/output_all --threads 20 5.4.Mono_harmony/cpdb_LO/input/celltype_all.txt 5.4.Mono_harmony/cpdb_LO/input/expr_all.txt
cellphonedb plot dot_plot --means-path 5.4.Mono_harmony/cpdb_LO/output_all/means.txt --pvalues-path 5.4.Mono_harmony/cpdb_LO/output_all/pvalues.txt --output-path 5.4.Mono_harmony/cpdb_LO/output_all
cellphonedb plot heatmap_plot --pvalues-path 5.4.Mono_harmony/cpdb_LO/output_all/pvalues.txt --output-path 5.4.Mono_harmony/cpdb_LO/output_all 5.4.Mono_harmony/cpdb_LO/input/celltype_all.txt

## Macro_APOE and Tcells (sample split) ---------------------------------------

PCSC <- readRDS("PCSC_nomix_new_SCT.RDS")
APOE_T <- subset(PCSC, celltype %in% c("Macro_APOE", "Tcm", "CD4_Tn", "Treg",
                                       "CD8_Tn", "CD8_Tem", "NK", "NKT"))
APOE_T$SampleCelltype <- paste0(APOE_T$orig.ident, "_", APOE_T$celltype)

expr <- as.matrix(APOE_T@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("5.4.Mono_harmony/cpdb/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
write.table(expr, paste0("5.4.Mono_harmony/cpdb/input_nosplit/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(APOE_T@meta.data), celltype = APOE_T@meta.data$SampleCelltype)  
write.table(celltype, paste0("5.4.Mono_harmony/cpdb/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype
celltype <- data.frame(Cell = rownames(APOE_T@meta.data), celltype = APOE_T@meta.data$celltype)  
write.table(celltype, paste0("5.4.Mono_harmony/cpdb/input_nosplit/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

cellphonedb method statistical_analysis --counts-data gene_name --output-path 5.4.Mono_harmony/cpdb/output_all --threads 20 5.4.Mono_harmony/cpdb/input/celltype_all.txt 5.4.Mono_harmony/cpdb/input/expr_all.txt
cellphonedb plot dot_plot --means-path 5.4.Mono_harmony/cpdb/output_all/means.txt --pvalues-path 5.4.Mono_harmony/cpdb/output_all/pvalues.txt --output-path 5.4.Mono_harmony/cpdb/output_all
cellphonedb plot heatmap_plot --pvalues-path 5.4.Mono_harmony/cpdb/output_all/pvalues.txt --output-path 5.4.Mono_harmony/cpdb/output_all 5.4.Mono_harmony/cpdb/input/celltype_all.txt

cellphonedb method statistical_analysis --counts-data gene_name --output-path 5.4.Mono_harmony/cpdb/output_nosplit --threads 20 5.4.Mono_harmony/cpdb/input_nosplit/celltype_all.txt 5.4.Mono_harmony/cpdb/input_nosplit/expr_all.txt
cellphonedb plot dot_plot --means-path 5.4.Mono_harmony/cpdb/output_nosplit/means.txt --pvalues-path 5.4.Mono_harmony/cpdb/output_nosplit/pvalues.txt --output-path 5.4.Mono_harmony/cpdb/output_nosplit
cellphonedb plot heatmap_plot --pvalues-path 5.4.Mono_harmony/cpdb/output_nosplit/pvalues.txt --output-path 5.4.Mono_harmony/cpdb/output_nosplit 5.4.Mono_harmony/cpdb/input_nosplit/celltype_all.txt

## Macro_APOE and Tcells (sampletype split) ---------------------------------------

PCSC <- readRDS("PCSC_nomix_new_SCT.RDS")
APOE_T <- subset(PCSC, celltype %in% c("Macro_APOE", "Tcm", "CD4_Tn", "Treg",
                                       "CD8_Tn", "CD8_Tem", "NK", "NKT"))
APOE_T$TypeCelltype <- paste0(APOE_T$SampleType, "_", APOE_T$celltype)

expr <- as.matrix(APOE_T@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("5.4.Mono_harmony/cpdb_SampleType/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(APOE_T@meta.data), celltype = APOE_T@meta.data$TypeCelltype)  
write.table(celltype, paste0("5.4.Mono_harmony/cpdb_SampleType/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

cellphonedb method statistical_analysis --counts-data gene_name --output-path 5.4.Mono_harmony/cpdb_SampleType/output_all --threads 20 5.4.Mono_harmony/cpdb_SampleType/input/celltype_all.txt 5.4.Mono_harmony/cpdb_SampleType/input/expr_all.txt
cellphonedb plot dot_plot --means-path 5.4.Mono_harmony/cpdb_SampleType/output_all/means.txt --pvalues-path 5.4.Mono_harmony/cpdb_SampleType/output_all/pvalues.txt --output-path 5.4.Mono_harmony/cpdb_SampleType/output_all
cellphonedb plot heatmap_plot --pvalues-path 5.4.Mono_harmony/cpdb_SampleType/output_all/pvalues.txt --output-path 5.4.Mono_harmony/cpdb_SampleType/output_all 5.4.Mono_harmony/cpdb_SampleType/input/celltype_all.txt

## intepret ----------------------------------------

library(ktplots)
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
library(CellChat)
library(tidyr)

color <- c("#A52A2A", "#CD7054", "#B8860B", "#CDAD00", 
           "#556B2F", "#008B8B", "#4A7088", "#7A378B")

### chord plot EOPC ----------------------------

pvals_EO <- read.delim("5.4.Mono_harmony/cpdb_EO/output_all/pvalues.txt", check.names = FALSE)
means_EO <- read.delim("5.4.Mono_harmony/cpdb_EO/output_all/means.txt", check.names = FALSE)
decon_EO <- read.delim("5.4.Mono_harmony/cpdb_EO/output_all/deconvoluted.txt", check.names = FALSE)
sigmeans_EO <- read.delim("5.4.Mono_harmony/cpdb_EO/output_all/significant_means.txt", check.names = FALSE)

count_net_EO <- read.delim("5.4.Mono_harmony/cpdb_EO/output_all/count_network.txt", check.names = FALSE)
inter_net_EO <- read.delim("5.4.Mono_harmony/cpdb_EO/output_all/interaction_count.txt", check.names = FALSE)

count_inter <- count_net_EO
#count_inter$count <- count_inter$count/50
count_inter<-spread(count_inter, TARGET, count)
rownames(count_inter) <- count_inter$SOURCE
count_inter <- count_inter[, -1]
count_inter <- as.matrix(count_inter)

color <- c("#A52A2A", "#2F4F4F", "#B8860B", "#4A7088", 
           "#556B2F", "#008B8B", "#CDAD00", "#7A378B")
pdf("5.4.Mono_harmony/cpdb_EO/chord_APOE_T.pdf", width = 4, height = 4)
netVisual_circle(count_inter, 
                 sources.use = "Macro_APOE",
                 weight.scale = T, 
                 color.use = color,
                 edge.weight.max = 200,
                 arrow.size=0,
                 label.edge = T,
                 vertex.label.cex = 0.8,
                 edge.label.cex = 0.6,
                 title.name = "EOPC")
dev.off()

### DotPlot Macro_APOE to T -------------------------------------------

MtoC_pval <- data.frame(pair = pvals_EOPC1$interacting_pair, EOPC1 = pvals_EOPC1$'Macro_APOE|Epi_MP1_On')
MtoC_mean <- data.frame(pair = means_EOPC1$interacting_pair, EOPC1 = means_EOPC1$'Macro_APOE|Epi_MP1_On')
MtoC_mean <- MtoC_mean[match(MtoC_pval$pair,MtoC_mean$pair),]

APOEtoT_pval <- pvals_EO[,which(substr(colnames(pvals_EO),1,10) == "Macro_APOE" & substr(colnames(pvals_EO),nchar(colnames(pvals_EO))-3,nchar(colnames(pvals_EO))) != "APOE")]
APOEtoT_pval <- data.frame(pair = pvals_EO$interacting_pair, APOEtoT_pval)
APOEtoT_mean <- means_EO[,which(substr(colnames(pvals_EO),1,10) == "Macro_APOE" & substr(colnames(pvals_EO),nchar(colnames(pvals_EO))-3,nchar(colnames(pvals_EO))) != "APOE")]
APOEtoT_mean <- data.frame(pair = means_EO$interacting_pair, APOEtoT_mean)
APOEtoT_mean <- APOEtoT_mean[match(APOEtoT_pval$pair,APOEtoT_mean$pair),]

colnames(APOEtoT_pval)[2:ncol(APOEtoT_pval)] <- substr(colnames(APOEtoT_pval)[2:ncol(APOEtoT_pval)],12,nchar(colnames(APOEtoT_pval)[2:ncol(APOEtoT_pval)]))
colnames(APOEtoT_mean)[2:ncol(APOEtoT_mean)] <- substr(colnames(APOEtoT_mean)[2:ncol(APOEtoT_mean)],12,nchar(colnames(APOEtoT_mean)[2:ncol(APOEtoT_mean)]))

sig_pvals <- APOEtoT_pval[-which(APOEtoT_pval$CD4_Tn > 0.05 &
                                   APOEtoT_pval$CD8_Tem > 0.05 &
                                   APOEtoT_pval$CD8_Tn > 0.05 &
                                   APOEtoT_pval$NK > 0.05 &
                                   APOEtoT_pval$NKT > 0.05 &
                                   APOEtoT_pval$Tcm > 0.05 &
                                   APOEtoT_pval$Treg > 0.05),]
sig_means <- APOEtoT_mean[match(sig_pvals$pair,APOEtoT_mean$pair),]

data_sig_pvals <- melt(sig_pvals)
colnames(data_sig_pvals)[3] <- "pvals"
data_sig_means <- melt(sig_means)
colnames(data_sig_means)[3] <- "means"

data <- merge(data_sig_pvals, data_sig_means, by = c("pair", "variable"))
data$pvals[data$pvals == 0] <- 0.001
data$pvals <- -log10(data$pvals)

P <- ggplot(data, aes(x=variable, y=pair, size=means)) + 
  geom_point(aes(colour = pvals)) +
  scale_size_continuous(range=c(1,4)) +
  theme_bw()+theme(axis.text.x = element_text(angle=90,size=10)) +
  scale_x_discrete(position = "bottom") +
  scale_colour_distiller(type = "seq", palette = "PuBu", direction = 1) +
  theme(panel.grid.major = element_blank(),
        legend.key.size = unit(7,'pt'),
        legend.key.height = unit(7,'pt'),
        legend.key.width = unit(7,'pt'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y = element_text(size=7, colour="black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_blank())
ggsave('./5.4.Mono_harmony/cpdb_EO/APOE_T_interaction.pdf', P, width = 4, height = 12)

### chord plot LOPC ----------------------------

pvals_LO <- read.delim("5.4.Mono_harmony/cpdb_LO/output_all/pvalues.txt", check.names = FALSE)
means_LO <- read.delim("5.4.Mono_harmony/cpdb_LO/output_all/means.txt", check.names = FALSE)
decon_LO <- read.delim("5.4.Mono_harmony/cpdb_LO/output_all/deconvoluted.txt", check.names = FALSE)
sigmeans_LO <- read.delim("5.4.Mono_harmony/cpdb_LO/output_all/significant_means.txt", check.names = FALSE)

count_net_LO <- read.delim("5.4.Mono_harmony/cpdb_LO/output_all/count_network.txt", check.names = FALSE)
inter_net_LO <- read.delim("5.4.Mono_harmony/cpdb_LO/output_all/interaction_count.txt", check.names = FALSE)

count_inter <- count_net_LO
#count_inter$count <- count_inter$count/50
count_inter<-spread(count_inter, TARGET, count)
rownames(count_inter) <- count_inter$SOURCE
count_inter <- count_inter[, -1]
count_inter <- as.matrix(count_inter)

color <- c("#A52A2A", "#2F4F4F", "#B8860B", "#CD7054", 
           "#556B2F", "#008B8B", "#CDAD00", "#7A378B")
pdf("5.4.Mono_harmony/cpdb_LO/chord_APOE_T.pdf", width = 4, height = 4)
netVisual_circle(count_inter, 
                 sources.use = "Macro_APOE",
                 weight.scale = T, 
                 color.use = color,
                 edge.weight.max = 200,
                 arrow.size=0,
                 label.edge = T,
                 vertex.label.cex = 0.8,
                 edge.label.cex = 0.6,
                 title.name = "EOPC")
dev.off()

### interaction comparison ------------------------------

count_net <- read.delim("/data_analysis_2/csl/cyf/SC/5.4.Mono_harmony/cpdb/output_all/count_network.txt")
count_net <- count_net[which(substr(count_net$SOURCE,1,5) == substr(count_net$TARGET,1,5)),]

EO <- count_net[which(substr(count_net$SOURCE,1,4) == "EOPC"),]
EO <- EO[which(substr(EO$SOURCE,7,nchar(EO$SOURCE)) == "Macro_APOE" & substr(EO$SOURCE,7,nchar(EO$TARGET)) != "Macro_APOE"),]
EO$Type <- "EOPC"
EO <- EO[,c(2,3,4)]
LO <- count_net[which(substr(count_net$SOURCE,1,4) == "LOPC"),]
LO <- LO[which(substr(LO$SOURCE,7,nchar(LO$SOURCE)) == "Macro_APOE" & substr(LO$SOURCE,7,nchar(LO$TARGET)) != "Macro_APOE"),]
LO$Type <- "LOPC"
LO <- LO[,c(2,3,4)]

com <- rbind(EO,LO)
com$cell <- substr(com$TARGET, 7, nchar(com$TARGET))
rownames(com) <- NULL
com <- column_to_rownames(com, "TARGET")

P = ggplot(com, aes(x = cell, y = count, fill = Type)) +
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
  ylab("Interaction Pairs") +
  xlab("") + 
  ggtitle("Interaction between Macro_APOE and T-cells") +
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

### dotplot EO vs LO --------------------------------------

pvals <- read.delim("5.4.Mono_harmony/cpdb_SampleType/output_all/pvalues.txt", check.names = FALSE)
means <- read.delim("5.4.Mono_harmony/cpdb_SampleType/output_all/means.txt", check.names = FALSE)
decon <- read.delim("5.4.Mono_harmony/cpdb_SampleType/output_all/deconvoluted.txt", check.names = FALSE)
sigmeans <- read.delim("5.4.Mono_harmony/cpdb_SampleType/output_all/significant_means.txt", check.names = FALSE)

#APOE to other
pvals_APOE_to_other <- pvals[, grep(".*Macro_APOE\\|.*", colnames(pvals))]
pvals_APOE_to_other <- pvals_APOE_to_other[, -grep(".*Macro_APOE$", colnames(pvals_APOE_to_other))]
pvals_APOE_to_other <- pvals_APOE_to_other[,grep("EOPC.*EOPC.*|LOPC.*LOPC.*", colnames(pvals_APOE_to_other))]
pvals_APOE_to_other <- cbind(pvals[,1:11], pvals_APOE_to_other)

means_APOE_to_other <- means[, grep(".*Macro_APOE\\|.*", colnames(means))]
means_APOE_to_other <- means_APOE_to_other[, -grep(".*Macro_APOE$", colnames(means_APOE_to_other))]
means_APOE_to_other <- means_APOE_to_other[,grep("EOPC.*EOPC.*|LOPC.*LOPC.*", colnames(means_APOE_to_other))]
means_APOE_to_other <- cbind(pvals[,1:11], means_APOE_to_other)

plot_cpdb(cell_type1 = 'Macro_APOE', cell_type2 = "", scdata = APOE_T,
          idents = 'celltype', split.by = "SampleType", means = means, 
          pvals = pvals, highlight = "#A52A2A", keep_significant_only = T,
          gene.family = c('costimulatory','coinhibitory'),
          highlight_size = NULL) +
  theme(axis.text  = element_text(size = 10, color = 'black'))


# CellChat "Secreted Signaling" and "Cell-Cell Contact"----------------------------------------------------------------------

library(CellChat)
library(patchwork)

PCSC <- readRDS("PCSC_nomix_new_SCT.RDS")
APOE_T <- subset(PCSC, celltype %in% c("Macro_APOE", "Tcm", "CD4_Tn", "Treg",
                                       "CD8_Tn", "CD8_Tem", "NK", "NKT"))
APOE_T$celltype <- as.character(APOE_T$celltype)

## EOPC ------------------------------------------------

### prepare ---------------------------------------

EOPC <- subset(APOE_T, SampleType == "EOPC")
data.input <- GetAssayData(EOPC, slot = "data")
meta = EOPC@meta.data
unique(meta$celltype)

## create CellChat

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

## database
CellChatDB <- CellChatDB.human
#showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","Cell-Cell Contact")) # use Secreted Signaling
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)   # This step is necessary even if using the whole database
future::plan("multicore", workers = 10)     # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

#### infer cell communications ------------------------------------------

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "5.4.Mono_harmony/CellChat/celltype_compare/net_lr_EOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "5.4.Mono_harmony/CellChat/celltype_compare/net_pathway_EOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
EOPC <- cellchat
saveRDS(EOPC, "5.4.Mono_harmony/CellChat/celltype_compare/cellchat_obj_EOPC.Rds")

### LOPC ------------------------------------------------

#### prepare ---------------------------------------

LOPC <- subset(APOE_T, SampleType == "LOPC")
data.input <- GetAssayData(LOPC, slot = "data")
meta = LOPC@meta.data
unique(meta$celltype)

## create CellChat

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

## database
CellChatDB <- CellChatDB.human
#showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","Cell-Cell Contact")) # use Secreted Signaling
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)   # This step is necessary even if using the whole database
future::plan("multicore", workers = 10)     # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

#### infer cell communications ------------------------------------------

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, "5.4.Mono_harmony/CellChat/celltype_compare/net_lr_LOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "5.4.Mono_harmony/CellChat/celltype_compare/net_pathway_LOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
LOPC <- cellchat
saveRDS(LOPC, "5.4.Mono_harmony/CellChat/celltype_compare/cellchat_obj_LOPC.Rds")

### integrating objects -------------
PCSC.list <- list(EOPC=EOPC, LOPC=LOPC)
cellchat <- mergeCellChat(PCSC.list, add.names = names(PCSC.list), cell.prefix = TRUE)
saveRDS(cellchat, "5.4.Mono_harmony/CellChat/celltype_compare/cellchat_obj_merge.Rds")

#### compare two groups ----------------------

#compare weight in particular cell groups

s.cell <- c("Macro_APOE")
r.cell <- c("CD4_Tn","CD8_Tem","CD8_Tn","NK","NKT","Tcm","Treg")
weight1 <- PCSC.list[[1]]@net$weight[s.cell,r.cell]
weight1 <- data.frame(Macro_APOE = weight1) %>% t
weight2 <- PCSC.list[[2]]@net$weight[s.cell,r.cell]
weight2 <- data.frame(Macro_APOE = weight2) %>% t
weight.max <- max(max(weight1), max(weight2))
pdf("5.4.Mono_harmony/CellChat/celltype_compare/EO_vs_LO_APOE_to_T.pdf")
for (i in 1:length(PCSC.list)) {
  par(mfrow = c(1,2))
  netVisual_circle(PCSC.list[[i]]@net$weight, 
                   weight.scale = T,
                   edge.weight.max = weight.max, 
                   edge.width.max = 12, 
                   sources.use = "Macro_APOE",
                   label.edge = F,
                   vertex.label.cex = 0.8,
                   title.name = paste0("Weight of interactions - ", names(PCSC.list)[i]))
}
dev.off()

#### 配体-受体对比分析 -----------------------------------
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(1,2,3,5,6,7,8),  comparison = c(1, 2), angle.x = 45)
ggsave("5.4.Mono_harmony/CellChat/celltype_compare/Compare_APOE_to_T_bubble.pdf", p, width = 5, height = 7)

# CellChat all interactions  ----------------------------------------------------------------------

library(CellChat)
library(patchwork)

PCSC <- readRDS("PCSC_nomix_new_SCT.RDS")
APOE_T <- subset(PCSC, celltype %in% c("Macro_APOE", "Tcm", "CD4_Tn", "Treg",
                                       "CD8_Tn", "CD8_Tem", "NK", "NKT"))
APOE_T$celltype <- as.character(APOE_T$celltype)

## EOPC ------------------------------------------------

### prepare ---------------------------------------

EOPC <- subset(APOE_T, SampleType == "EOPC")
data.input <- GetAssayData(EOPC, slot = "data")
meta = EOPC@meta.data
unique(meta$celltype)

## create CellChat

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

## database
CellChatDB <- CellChatDB.human
#showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
#CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","Cell-Cell Contact")) # use Secreted Signaling
cellchat@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)   # This step is necessary even if using the whole database
future::plan("multicore", workers = 10)     # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

#### infer cell communications ------------------------------------------

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "5.4.Mono_harmony/CellChat/all_interactions/net_lr_EOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "5.4.Mono_harmony/CellChat/all_interactions/net_pathway_EOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
EOPC <- cellchat
saveRDS(EOPC, "5.4.Mono_harmony/CellChat/all_interactions/cellchat_obj_EOPC.Rds")

### LOPC ------------------------------------------------

#### prepare ---------------------------------------

LOPC <- subset(APOE_T, SampleType == "LOPC")
data.input <- GetAssayData(LOPC, slot = "data")
meta = LOPC@meta.data
unique(meta$celltype)

## create CellChat

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

## database
CellChatDB <- CellChatDB.human
#showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
#CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","Cell-Cell Contact")) # use Secreted Signaling
cellchat@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)   # This step is necessary even if using the whole database
future::plan("multicore", workers = 10)     # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

#### infer cell communications ------------------------------------------

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, "5.4.Mono_harmony/CellChat/all_interactions/net_lr_LOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "5.4.Mono_harmony/CellChat/all_interactions/net_pathway_LOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
LOPC <- cellchat
saveRDS(LOPC, "5.4.Mono_harmony/CellChat/all_interactions/cellchat_obj_LOPC.Rds")

### integrating objects -------------
PCSC.list <- list(EOPC=EOPC, LOPC=LOPC)
cellchat <- mergeCellChat(PCSC.list, add.names = names(PCSC.list), cell.prefix = TRUE)
saveRDS(cellchat, "5.4.Mono_harmony/CellChat/all_interactions/cellchat_obj_merge.Rds")

#### EOPC, Mye to T, visualization ----------------------------------------------------------------------

EOPC <- readRDS("5.4.Mono_harmony/CellChat/all_interactions/cellchat_obj_EOPC.Rds")
levels(EOPC@idents)
p <- netVisual_bubble(EOPC, sources.use = c(4), targets.use = c(1,2,3,5,6,7,8), angle.x = 45)

#### compare two groups ----------------------

#compare weight in particular cell groups

s.cell <- c("Macro_APOE")
r.cell <- c("CD4_Tn","CD8_Tem","CD8_Tn","NK","NKT","Tcm","Treg")
weight1 <- PCSC.list[[1]]@net$weight[s.cell,r.cell]
weight1 <- data.frame(Macro_APOE = weight1) %>% t
weight2 <- PCSC.list[[2]]@net$weight[s.cell,r.cell]
weight2 <- data.frame(Macro_APOE = weight2) %>% t
weight.max <- max(max(weight1), max(weight2))
pdf("5.4.Mono_harmony/CellChat/all_interactions/EO_vs_LO_APOE_to_T.pdf")
for (i in 1:length(PCSC.list)) {
  par(mfrow = c(1,2))
  netVisual_circle(PCSC.list[[i]]@net$weight, 
                   weight.scale = T,
                   edge.weight.max = weight.max, 
                   edge.width.max = 12, 
                   sources.use = "Macro_APOE",
                   label.edge = F,
                   vertex.label.cex = 0.8,
                   title.name = paste0("Weight of interactions - ", names(PCSC.list)[i]))
}
dev.off()

#### 配体-受体对比分析 -----------------------------------
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(1,2,3,5,6,7,8),  comparison = c(1, 2), angle.x = 45)
ggsave("5.4.Mono_harmony/CellChat/all_interactions/Compare_APOE_to_T_bubble.pdf", p, width = 7, height = 7)

# Macro_APOE in TCGA EO vs LO / MP1 / survival ---------------------------------

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

#### calculate Macro_APOE score ------------------------------

Allmarkers <- readRDS("./5.4.Mono_harmony/Markers/Nomix_celltype_markers_nomix.Rds")
APOE <- Allmarkers[Allmarkers$cluster == "Macro_APOE",]
APOE_markers <- APOE[which(APOE$p_val_adj < 0.01 & APOE$pct.1 > 0.3 & APOE$avg_log2FC > 0.25),]

library(AUCell)
gs <- list(Macro_APOE = APOE_markers$gene)

cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.4.Mono_harmony/Bulk_APOE_score/APOE_gs_AUC.Rds")

APOE <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

#### calculate Macro_gs score ------------------------------

Macrp_gs <- read.table("Macrophage_gs.txt", header = T, sep = "\t")

gs <- list()
for (i in c("M1","M2","MDSC","Angiogenesis","Phagocytosis")) {
  gs[[i]] <- Macrp_gs[Macrp_gs$term == i,]$gene
}

library(AUCell)

cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.4.Mono_harmony/Bulk_APOE_score/Macro_gs_AUC.Rds")

Macro <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

##### MP1 and APOE --------------------------------------------------------

MP1_APOE <- merge(t(AUCell_mtx), APOE, by = 0) %>% column_to_rownames("Row.names")
MP1_APOE <- MP1_APOE[-grep(".11$|.06$", rownames(MP1_APOE)),]

for (i in (7 : (ncol(MP1_APOE)))) {
  
  P <- ggplot(MP1_APOE, aes(x=MP1_APOE[,1], y=MP1_APOE[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = MP1_APOE[,1], y = MP1_APOE[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(MP1_APOE)[1]," AUCell")) + ylab(paste0(colnames(MP1_APOE)[i]," AUCell"))
  
  ggsave(paste0("./5.4.Mono_harmony/Bulk_APOE_score/CorPlot_MP1_",colnames(MP1_APOE)[i],".pdf"), P, width = 2, height = 2)
}

##### M2/MDSC and APOE --------------------------------------------------------

APOE_M2_MDSC <- merge(Macro, APOE, by = 0) %>% column_to_rownames("Row.names")
APOE_M2_MDSC <- APOE_M2_MDSC[-grep(".11$|.06$", rownames(APOE_M2_MDSC)),]
APOE_M2_MDSC <- APOE_M2_MDSC[,c(2,3,6)]

cor.test(APOE_M2_MDSC$M2, APOE_M2_MDSC$Macro_APOE, method = "spearman")
cor.test(APOE_M2_MDSC$MDSC, APOE_M2_MDSC$Macro_APOE, method = "spearman")

for (i in c(1,2)) {
  
  P <- ggplot(APOE_M2_MDSC, aes(x=APOE_M2_MDSC[,i], y=APOE_M2_MDSC[,3])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = APOE_M2_MDSC[,i], y = APOE_M2_MDSC[,3]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(APOE_M2_MDSC)[i]," AUCell")) + ylab(paste0(colnames(APOE_M2_MDSC)[3]," AUCell"))
  
  ggsave(paste0("./5.4.Mono_harmony/Bulk_APOE_score/CorPlot_APOE_",colnames(APOE_M2_MDSC)[i],".pdf"), P, width = 2, height = 2)
}

##### EO vs LO and APOE --------------------------------------------------------

Age_APOE <- merge(age, APOE, by = 0) %>% column_to_rownames("Row.names")

# correlation

for (i in (3 : (ncol(Age_APOE)))) {
  
  P <- ggplot(Age_APOE, aes(x=Age_APOE[,1], y=Age_APOE[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = Age_APOE[,1], y = Age_APOE[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_APOE)[i]," AUCell"))
  
  ggsave(paste0("./5.4.Mono_harmony/Bulk_APOE_score/CorPlot_Age_",colnames(Age_APOE)[i],".pdf"), P, width = 2, height = 2)
}

for (i in (3 : (ncol(Age_APOE)))) {
  
  P <- ggplot(Age_APOE, aes(x=Age_APOE[,1], y=Age_APOE[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'pearson', aes(x = Age_APOE[,1], y = Age_APOE[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_APOE)[i]," AUCell"))
  
  ggsave(paste0("./5.4.Mono_harmony/Bulk_APOE_score/CorPlot_pearson_Age_",colnames(Age_APOE)[i],".pdf"), P, width = 2, height = 2)
}

# boxplot

data <- melt(Age_APOE[,c(2:ncol(Age_APOE))])
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
pdf("./5.4.Mono_harmony/Bulk_APOE_score/BoxPlot_APOE_EOLO.pdf", width = 2, height = 2.5)
P + NoLegend()
dev.off()

##### EO vs LO and M2/MDSC --------------------------------------------------------

Age_Macro <- merge(age, Macro, by = 0) %>% column_to_rownames("Row.names")
Age_Macro <- Age_Macro[, c(1,2,4,5)]

# correlation

for (i in (3 : (ncol(Age_Macro)))) {
  
  P <- ggplot(Age_Macro, aes(x=Age_Macro[,1], y=Age_Macro[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = Age_Macro[,1], y = Age_Macro[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_Macro)[i]," AUCell"))
  
  ggsave(paste0("./5.4.Mono_harmony/Bulk_APOE_score/CorPlot_Age_",colnames(Age_Macro)[i],".pdf"), P, width = 2, height = 2)
}

for (i in (3 : (ncol(Age_Macro)))) {
  
  P <- ggplot(Age_Macro, aes(x=Age_Macro[,1], y=Age_Macro[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'pearson', aes(x = Age_Macro[,1], y = Age_Macro[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_Macro)[i]," AUCell"))
  
  ggsave(paste0("./5.4.Mono_harmony/Bulk_APOE_score/CorPlot_pearson_Age_",colnames(Age_Macro)[i],".pdf"), P, width = 2, height = 2)
}

# boxplot

data <- melt(Age_Macro[,c(2:ncol(Age_Macro))])
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
pdf("./5.4.Mono_harmony/Bulk_APOE_score/BoxPlot_APOE_EOLO.pdf", width = 2, height = 2.5)
P + NoLegend()
dev.off()

#### APOE and survival --------------------------------------

survival_APOE <- merge(APOE, survival, by = 0) %>% column_to_rownames("Row.names")
survival_APOE <- survival_APOE[-grep(".11$|.06$", rownames(survival_APOE)),]
survival_APOE$OS.time <-survival_APOE$OS.time / 365

coxph(Surv(OS.time, OS) ~ Macro_APOE, survival_APOE)

library(survival)
library("survminer")

quar <- quantile(survival_APOE$Macro_APOE, probs = seq(0, 1, 0.25))

survival_APOE$APOE_HL <- 0
survival_APOE$APOE_HL[which(survival_APOE$Macro_APOE > quar[2])] <- 1

diff <- survdiff(Surv(OS.time, OS) ~ APOE_HL, data = survival_APOE)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
P1 <- pValue
fit <- survfit(Surv(OS.time, OS) ~ APOE_HL, data = survival_APOE)

ggsurvplot(fit, 
           data=survival_APOE,
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
dev.off()

#### APOE and survival in EOPC --------------------------------------

survival_Age <- merge(survival, age, by = 0) %>% column_to_rownames("Row.names")
survival_EOPC <- survival_Age[survival_Age$Type == "EOPC",]
survival_APOE <- merge(APOE, survival_EOPC, by = 0) %>% column_to_rownames("Row.names")
survival_APOE$OS.time <-survival_APOE$OS.time / 365

coxph(Surv(OS.time, OS) ~ Macro_APOE, survival_APOE)

library(survival)
library("survminer")

quar <- quantile(survival_APOE$Macro_APOE, probs = seq(0, 1, 0.25))

survival_APOE$APOE_HL <- 0
survival_APOE$APOE_HL[which(survival_APOE$Macro_APOE > quar[3])] <- 1

diff <- survdiff(Surv(OS.time, OS) ~ APOE_HL, data = survival_APOE)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
P1 <- pValue
fit <- survfit(Surv(OS.time, OS) ~ APOE_HL, data = survival_APOE)

pdf("./5.4.Mono_harmony/Bulk_APOE_score/APOE_survival_in_EOPC.pdf", width = 4, height = 4)
ggsurvplot(fit, 
           data=survival_APOE,
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
dev.off()
