library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)

color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")

# SC data -------------------------------------------

PCSC <- readRDS("./6.0.cellchat_Epi_all/Epi_all_annotated_new.Rds")
PCSC <- SCTransform(PCSC, assay = "RNA", vars.to.regress = c("percent.mt"))
saveRDS(PCSC, "PCSC_nomix_new_SCT.RDS")
Idents(PCSC) <- "inferCNVCellStat"
PCSC <- readRDS("PCSC_nomix_SCT.RDS")

##variable genes

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)
cluster_markers_all <- FindAllMarkers(object = PCSC, 
                                      assay = "SCT",
                                      slot = "data",
                                      verbose = TRUE, 
                                      only.pos = TRUE)
saveRDS(cluster_markers_all, "./12.ST_LOPC/3.anchoring/markers_sc_celltype.RDS")
saveRDS(cluster_markers_all, "./12.ST_LOPC/3.anchoring/markers_sc_celltype.RDS")

saveRDS(cluster_markers_all, "./12.ST_LOPC/3.anchoring/markers_sc_nomix_celltype.RDS")
saveRDS(cluster_markers_all, "./12.ST_LOPC/3.anchoring/markers_sc_nomix_celltype.RDS")

#random 1/10 cells

cells <- PCSC$inferCNVCellStat %>% as.data.frame()
colnames(cells) <- "celltype"

celltype <- c("Cancer","T-cell","B-cell", "Myeloid",
              "Mast","Endothelial","Fibroblast", "SMC")

set.seed(20230327)
RandomCell <- c()
for (i in celltype) {
  cell_tmp <- rownames(cells)[cells$celltype == i]
  n = length(cell_tmp) / 10
  if (n > 1000) {n = 1000} else
    if (n < 200) {n = 200} else
    {n = n}
  cell_tmp <- sample(cell_tmp, n)
  RandomCell <- c(RandomCell, cell_tmp)
}

PCSC_RandomCell <- PCSC[, RandomCell]
Idents(PCSC_RandomCell) <- "inferCNVCellStat"
saveRDS(PCSC_RandomCell, "PCSC_RandomCell_new.Rds")

# 1.1.Data processing - EOPC4-K -------------------------------------------------------

# load data from "filtered_feature_bc_matrix.h5" & "spatial"

EOPC4_ST <- Load10X_Spatial(data.dir = "DATA/PCSC008-K/", 
                            assay = "Spatial", 
                            slice = "slice", 
                            filter.matrix = TRUE)

EOPC4_ST@meta.data$orig.ident <- "EOPC"

# 1.2.ST ---------------------------------------------------------------------------------

## Data visualization ----------------------------

plot1 <- VlnPlot(EOPC4_ST, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(EOPC4_ST, features = "nCount_Spatial") + theme(legend.position = "right")

pdf("./11.ST_EOPC/1.QC/nCounts.pdf", height = 5, width = 12)
wrap_plots(plot1, plot2)
dev.off()

plot1 <- VlnPlot(EOPC4_ST, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(EOPC4_ST, features = "nFeature_Spatial") + theme(legend.position = "right")

pdf("./11.ST_EOPC/1.QC/nFeature.pdf", height = 5, width = 12)
wrap_plots(plot1, plot2)
dev.off()

## Normalization -----------------------------

EOPC4_ST[["percent.mt"]] <- PercentageFeatureSet(EOPC4_ST, pattern = "^MT-")
EOPC4_ST <- NormalizeData(EOPC4_ST, verbose = TRUE)
str(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
EOPC4_ST <- CellCycleScoring(EOPC4_ST, g2m.features=g2m.genes, s.features=s.genes)
EOPC4_ST <- SCTransform(EOPC4_ST, assay = "Spatial", vars.to.regress = c("percent.mt"))

## PCA ---------------------------------------

EOPC4_ST <- RunPCA(EOPC4_ST, verbose = FALSE)   %>%
      FindNeighbors(reduction = "pca", dims = 1:30) 

for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 0.9,1)) {
  EOPC4_ST <- FindClusters(EOPC4_ST, resolution = res)}

apply(EOPC4_ST@meta.data[, grep("SCT_snn_res", colnames(EOPC4_ST@meta.data))],2,table)

library(clustree)

p2_tree=clustree(EOPC4_ST@meta.data, prefix = "SCT_snn_res.") 
ggsave(plot=p2_tree, filename="./11.ST_EOPC/2.clustering/Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(EOPC4_ST) <- "SCT_snn_res.0.9"

#tsne & UMAP

EOPC4_ST <- RunUMAP(EOPC4_ST, reduction = "pca", dims = 1:30)
saveRDS(EOPC4_ST, "11.ST_EOPC/EOPC4_ST_clustered.Rds")

#plot

p1 <- DimPlot(EOPC4_ST, label = TRUE, shuffle = T, pt.size = 1, label.size = 4, label.box = T) + 
  scale_color_manual(values = color) + 
  scale_fill_manual(values = color) + NoLegend()
p2 <- SpatialDimPlot(EOPC4_ST, label = TRUE, label.size = 3, pt.size.factor = 1.2) + 
  scale_fill_manual(values = color) + NoLegend()

pdf("./11.ST_EOPC/2.clustering/umap_plot.pdf", height = 5, width = 10)
wrap_plots(p1, p2)
dev.off()

pdf("./11.ST_EOPC/2.clustering/spatial_plot_split_by_clusters.pdf", width = 14, height = 14)
p4 <- SpatialDimPlot(EOPC4_ST, 
                     cells.highlight = CellsByIdentities(object = EOPC4_ST),
                     facet.highlight = TRUE, ncol = 4) + NoLegend()
p4
dev.off()

#interactive plot

#SpatialDimPlot(EOPC4_ST, interactive = TRUE)
#SpatialFeaturePlot(EOPC4_ST, features = c("EPCAM", "MLH1", "MSH2", 
#                                    "MSH6", "PMS2"))
#LinkedDimPlot(EOPC4_ST)

#pdf("./2.clustering/cluster0.pdf", height = 5, width = 10)
#LinkedDimPlot(EOPC4_ST)
#dev.off()

## Find spatially variable genes ---------------------------------------

EOPC4_ST <- FindSpatiallyVariableFeatures(EOPC4_ST, assay = "SCT",
                                          features = VariableFeatures(EOPC4_ST)[1:1000], 
                                          selection.method = "markvariogram")
saveRDS(EOPC4_ST, "./11.ST_EOPC/2.clustering/EOPC4_clustered.Rds")

top.features <- head(SpatiallyVariableFeatures(EOPC4_ST, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(EOPC4_ST, features = top.features, ncol = 3, alpha = c(0.1, 1))

SpatialFeaturePlot(EOPC4_ST, features = c("LYZ"), interactive = FALSE)

# 1.3.1.anchored (Seurat) ---------------------------------

anchors <- FindTransferAnchors(reference = PCSC, query = EOPC4_ST, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = PCSC$inferCNVCellStat, prediction.assay = TRUE, 
                                  weight.reduction = EOPC4_ST[["pca"]], dims = 1:30)
EOPC4_ST[["predictions"]] <- predictions.assay
saveRDS(EOPC4_ST, "./11.ST_EOPC/3.anchoring/EOPC4_ST_celltype_SeuratAnchored.Rds")

DefaultAssay(EOPC4_ST) <- "predictions"
pdf("./11.ST_EOPC/3.anchoring/anchor_ST_celltype_Seurat.pdf")
SpatialFeaturePlot(EOPC4_ST, rownames(EOPC4_ST), pt.size.factor = 1.6, ncol = 4, crop = TRUE)
dev.off()

# 1.3.2.anchored (spacexr) ---------------------------------

library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)

#prepare SC data
PCSC_RandomCell <- readRDS("PCSC_RandomCell_new.Rds")
counts <- as.matrix(PCSC_RandomCell@assays$RNA@counts)
meta_data <- PCSC_RandomCell@meta.data
cell_types <- meta_data$inferCNVCellStat
names(cell_types) <- rownames(meta_data)
cell_types <- as.factor(cell_types)
nUMI <- meta_data$nCount_RNA
names(nUMI) <- rownames(meta_data)
### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)

#prepare ST data
coords <- GetTissueCoordinates(EOPC4_ST, scale = NULL) 
spatialcounts <- as.matrix(EOPC4_ST@assays$Spatial@counts)
nUMI <- EOPC4_ST@meta.data$nCount_Spatial
names(nUMI) <- rownames(EOPC4_ST@meta.data)
### Create SpatialRNA object
puck <- SpatialRNA(coords, spatialcounts, nUMI)

#RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') #doublet_mode = 'doublet'

saveRDS(myRCTD, "./11.ST_EOPC/3.anchoring/spacexr_anchored.Rds")

#1.3.3. annotation ---------------------------------------------

Features=c("EPCAM","CDH1", 
           "AR", "AMACR", 
           "TP63", "KRT5",
           "PTPRC", 
           "CD3E", "CD3G", "CD2", 
           "MS4A1", "CD79A", 
           "CD14", "FCGR3A", "LYZ", 
           "KIT", "MS4A2", "TPSAB1", 
           "VIM", 
           "PECAM1", "ENG", "VWF", 
           "FAP", "PDGFRA", "THY1", 
           "COL3A1", "ACTA2", "MYL9")

pdf("./11.ST_EOPC/2.clustering/AllmarkerBubble.pdf",width = 10, height = 6)
jjDotPlot(EOPC4_ST,
          gene = Features,
          id = "SCT_snn_res.0.9",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 6,
          x.text.vjust = 0.5)
dev.off()

for (i in Features) {
  Scatter <- FeaturePlot(EOPC4_ST, features = i, cols = c("grey", "#A52A2A"),
                         pt.size = 0.3, ncol = ) + NoLegend()
  ggsave(plot = Scatter, device = "tiff", filename = paste0("./11.ST_EOPC/2.clustering/Scatter/", i, ".tiff"), width = 2, height = 2)
}


celltype <- data.frame(cluster = levels(EOPC4_ST@meta.data$SCT_snn_res.0.9), celltype = "Cancer-1")
celltype$celltype[celltype$cluster %in% c(5,6)] <- "Cancer-2"
celltype$celltype[celltype$cluster %in% c(8)] <- "Cancer-3"
celltype$celltype[celltype$cluster %in% c(1)] <- "Myeloid"
celltype$celltype[celltype$cluster %in% c(0)] <- "SMC-1"
celltype$celltype[celltype$cluster %in% c(2,9,10,13,14)] <- "SMC-2"
celltype$celltype[celltype$cluster %in% c(12)] <- "Fibroblast"
celltype$celltype[celltype$cluster %in% c(7)] <- "Endothelial"

EOPC4_ST@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  EOPC4_ST@meta.data[which(EOPC4_ST@meta.data$SCT_snn_res.0.9 == celltype$cluster[i]),'celltype'] <- celltype$celltype[i]}
table(EOPC4_ST@meta.data$celltype)

Idents(EOPC4_ST) <- "celltype"
Idents(EOPC4_ST) <- factor(Idents(EOPC4_ST), levels = c("Cancer-1", "Cancer-2", "Cancer-3", "Fibroblast", 
                                                        "Myeloid", "SMC-1", "SMC-2", "Endothelial"))

saveRDS(EOPC4_ST, "./11.ST_EOPC/EOPC4_ST_annotated.Rds")

p1 <- DimPlot(EOPC4_ST, label = TRUE, shuffle = T, pt.size = 1, label.size = 4, label.box = T) + 
  scale_color_manual(values = color[c(1,5,3,7,2,4,6,14)]) + 
  scale_fill_manual(values = color[c(1,5,3,7,2,4,6,14)]) + NoLegend()
p2 <- SpatialDimPlot(EOPC4_ST, label = FALSE, pt.size.factor = 1.5) + 
  scale_fill_manual(values = color[c(1,5,3,7,2,4,6,14)]) + NoLegend()

pdf("./11.ST_EOPC/2.clustering/celltype_umap_plot.pdf", height = 5, width = 10)
wrap_plots(p1, p2)
dev.off()

pdf("./11.ST_EOPC/2.clustering/celltype_Spatial_plot.pdf", height = 4, width = 4)
SpatialDimPlot(EOPC4_ST, label = FALSE, pt.size.factor = 1.5) + 
  scale_fill_manual(values = color[c(1,5,3,7,2,4,6,14)])
dev.off()

## DotPlot

EOPC4_ST$celltype <- factor(EOPC4_ST$celltype, 
                            levels = c("Cancer-1","Cancer-2","Cancer-3","Myeloid",
                                       "Endothelial","Fibroblast","SMC-1","SMC-2"))

pdf("./11.ST_EOPC/2.clustering/AllmarkerBubble.pdf",width = 10, height = 6)
jjDotPlot(EOPC4_ST,
          gene = Features,
          id = "celltype",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 6,
          x.text.vjust = 0.5)
dev.off()

# 1.3.4.Hypoxia score ----------------------------------------

hsets <- read.gmt("./5.1.1.epithelial_NMF/hallmark_cancersea.gmt")
hsets$term <- gsub("_", " ", hsets$term)
hsets$term <- str_to_title(hsets$term)

gs_list <- list(Hypoxia = hsets$gene[hsets$term == "Hypoxia"])

expMtx <- GetAssayData(EOPC4_ST, assay = "Spatial", slot = "data")
GSVA <- gsva(expMtx, gs_list, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "11.ST_EOPC/Hypoxia_GSVA.Rds")

expMtx <- as.matrix(EOPC4_ST@assays$Spatial@data)
cells_rankings <- AUCell_buildRankings(expMtx, nCores=10, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs_list, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./11.ST_EOPC/HypoxiaScore_AUC.Rds")

EOPC4_ST[["Hypoxia"]] <- CreateAssayObject(counts = cells_AUC@assays@data$AUC)
DefaultAssay(EOPC4_ST) <- "Hypoxia"

pdf("./11.ST_EOPC/SpatialPlot_Hypoxia_AUCell.pdf")
SpatialFeaturePlot(EOPC4_ST, features = "Hypoxia")
dev.off()

# 1.3.5.MP1 score -----------------------------------------

Epi <- subset(EOPC4_ST, celltype %in% c("Cancer-1","Cancer-2","Cancer-3"))

library(AUCell)
Markers <- readRDS("./5.1.1.epithelial_NMF/process_k3/MP_markers.Rds")
expMtx <- as.matrix(Epi@assays$Spatial@data)
cells_rankings <- AUCell_buildRankings(expMtx, nCores=10, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./11.ST_EOPC/4.Epi/MPScore_AUC.Rds")

Epi[["MP_AUCell"]] <- CreateAssayObject(counts = cells_AUC@assays@data$AUC)
DefaultAssay(Epi) <- "MP_AUCell"

pdf("./11.ST_EOPC/4.Epi/SpatialPlot_MP1_score.pdf", width = 5, height = 5)
SpatialFeaturePlot(Epi, features = "MP1")
dev.off()

#boxplot

exp <- data.frame(celltype=Epi$celltype, MP1=Epi@assays$MP_AUCell@data["MP1",])
rownames(exp) <- colnames(Epi)
exp$celltype <- factor(exp$celltype, levels = c("Cancer-1", "Cancer-2", "Cancer-3"))
p=ggboxplot(exp, x="celltype", y="MP1",fill= "celltype",
            ylab="MP1 AUCell score",
            xlab="",
            palette = color[c(1,5,3)])
p=p + rotate_x_text(45) + theme(axis.title.y = element_text(vjust=2, size=18))
p1=p+stat_compare_means(aes(group=celltype),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","")),
                        label = "p.signif")
p1

pdf("./11.ST_EOPC/4.Epi/BoxPlot_MP1Score_CancerGroup.pdf", width = 3, height = 5)
p1 + NoLegend()
dev.off()

wilcox.test(exp$MP1[exp$celltype == "Cancer-1"], exp$MP1[exp$celltype == "Cancer-2"])
wilcox.test(exp$MP1[exp$celltype == "Cancer-1"], exp$MP1[exp$celltype == "Cancer-3"])
wilcox.test(exp$MP1[exp$celltype == "Cancer-2"], exp$MP1[exp$celltype == "Cancer-3"])


### fatty acid metabolism ---------------------------------------

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

expMtx <- GetAssayData(Epi, assay = "Spatial", slot = "data")
GSVA <- gsva(expMtx, gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "11.ST_EOPC/4.Epi/Lipid_5gs_GSVA.Rds")

Epi[["lipid"]] <- CreateAssayObject(counts = GSVA)
DefaultAssay(Epi) <- "lipid"

pdf("11.ST_EOPC/4.Epi/lipid/Vlnplot_5gs.pdf")
VlnPlot(Epi, features = rownames(Epi), cols = color[c(1,5,3)], pt.size = 0)
dev.off()

pdf("11.ST_EOPC/4.Epi/lipid/Boxplot_5gs.pdf")

celltype <- Epi$celltype %>% as.data.frame()
colnames(celltype) <- "celltype"

mtx <- cbind(t(as.data.frame(Epi@assays$lipid@data)), celltype)
mtx <- melt(mtx)

ggplot(mtx, aes(x = variable, y = value, fill = celltype)) +
  scale_fill_manual(values = paste0(color, "CC")) +
  geom_boxplot(width = 0.8, outlier.colour = NA, position=position_dodge(0.9), notch = T) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("GSVA Score") +
  xlab("") +
  scale_fill_manual(values = color[c(1,5,3)]) +
  rotate_x_text(angle = 45) +
  stat_compare_means(aes(group=celltype),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","")),
                     label = "p.signif")

dev.off()

# 1.3.6.Myeloid cells -----------------------------------------

Myeloid <- subset(EOPC4_ST, celltype == "Myeloid")

Myeloid <- SCTransform(Myeloid, assay = "Spatial", vars.to.regress = c("percent.mt"))

## PCA ---------------------------------------

Myeloid <- RunPCA(Myeloid, verbose = FALSE)   %>%
  FindNeighbors(reduction = "pca", dims = 1:30) 

for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 0.9,1)) {
  Myeloid <- FindClusters(Myeloid, resolution = res)}

apply(Myeloid@meta.data[, grep("SCT_snn_res", colnames(Myeloid@meta.data))],2,table)

library(clustree)

p2_tree=clustree(Myeloid@meta.data, prefix = "SCT_snn_res.") 
ggsave(plot=p2_tree, filename="./11.ST_EOPC/5.Myeloid/Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(Myeloid) <- "SCT_snn_res.0.8"

#tsne & UMAP

Myeloid <- RunUMAP(Myeloid, reduction = "pca", dims = 1:30)

#plot

p1 <- DimPlot(Myeloid, label = TRUE, shuffle = T, pt.size = 1, label.size = 4, label.box = T) + 
  scale_color_manual(values = color) + 
  scale_fill_manual(values = color) + NoLegend()
p2 <- SpatialDimPlot(Myeloid, label = TRUE, label.size = 3, pt.size.factor = 2) + 
  scale_fill_manual(values = color) + NoLegend()

pdf("./11.ST_EOPC/5.Myeloid/umap_plot.pdf", height = 5, width = 10)
wrap_plots(p1, p2)
dev.off()

saveRDS(Myeloid, "./11.ST_EOPC/5.Myeloid/Myeloid_clustered.Rds")

## cluster markers --------------------------------

library(SCP)

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)

Idents(Myeloid) <- "SCT_snn_res.0.8"
Allmarkers <- FindAllMarkers(Myeloid,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = T)
saveRDS(Allmarkers, "./11.ST_EOPC/5.Myeloid/Markers/cellcluster_markers.Rds")

sigmarkers <- Allmarkers[Allmarkers$p_val_adj < 0.01 & Allmarkers$avg_log2FC > 0,]

ht <- FeatureHeatmap(
  srt = Myeloid, 
  group.by = "SCT_snn_res.0.8", 
  features = sigmarkers$gene, 
  feature_split = sigmarkers$cluster,
  species = "Homo_sapiens", 
  db = c("GO_BP"), 
  anno_terms = TRUE,
  nlabel = 20,
  height = 9, 
  width = 8,
  group_palcolor = color,
  feature_split_palcolor = color,
  topTerm = 10
)

pdf("./11.ST_EOPC/5.Myeloid/Markers/Heatmap_celltype_siggene_topGOBP.pdf", width = 20, height = 20)
ht$plot
dev.off()

## SlingShot ---------------------------------------

library(SCP)

Myeloid <- RunSlingshot(
  srt = Myeloid, 
  group.by = "SCT_snn_res.0.8", 
  reduction = "UMAP")
dev.off()

pdf("./11.ST_EOPC/5.Myeloid/Slingshot/Trajectory.pdf", width = 10, height = 5)
CellDimPlot(Myeloid, 
            group.by = "SCT_snn_res.0.8", 
            reduction = "UMAP",
            dims = c(1, 2), 
            lineages = paste0("Lineage", 1),
            palcolor = color)
dev.off()

## Monocle2 ----------------------------------------

#构造表达及注释数据
Myeloid <- NormalizeData(Myeloid)
exp.matrix <- as(as.matrix(Myeloid@assays$Spatial@data), "sparseMatrix")
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(exp.matrix)
exp_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- Myeloid@meta.data
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
diff_test_res<-differentialGeneTest(exp.monocle,fullModelFormulaStr = "~SCT_snn_res.0.8") 
ordering_genes<-row.names (subset(diff_test_res, qval < 0.05))
saveRDS(ordering_genes, "./11.ST_EOPC/5.Myeloid/Monocle/ordering_genes.Rds")
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
#plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

saveRDS(exp.monocle, "./11.ST_EOPC/5.Myeloid/Monocle/exp.monocle.Rds")

#将不同分组情况的拟时序轨迹图画到一起

plot1<-plot_cell_trajectory(exp.monocle, color_by = "SCT_snn_res.0.8",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
plot2<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size = 0.2, shuffle = T) + 
  scale_color_gradient(low = "#101D44", high = "#9FB2ED")
plot3<-plot_cell_trajectory(exp.monocle, color_by = "SampleType",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
pdf("./5.4.Mono_harmony/Monocle/Nomix_trajectory_plot.pdf",width = 3,height = 9)
CombinePlots(plots = list(plot1,plot2,plot3),legend = NULL,ncol=1)
dev.off()

## Seurat anchoring ----------------------------------------------

Mono_nomix <- readRDS("./5.4.Mono_harmony/Nomix_Myeloid_annoted.Rds")
Mono_nomix <- Mono_nomix %>% subset(celltype != "Macro_FTL")
Mono_nomix$celltype <- as.character(Mono_nomix$celltype)
Mono_nomix <- SCTransform(Mono_nomix, vars.to.regress = "percent.mt")
DefaultAssay(Mono_nomix) <- "RNA"
DefaultAssay(Myeloid) <- "Spatial"

anchors <- FindTransferAnchors(reference = Mono_nomix, query = Myeloid, normalization.method = "LogNormalize")
predictions.assay <- TransferData(anchorset = anchors, refdata = Mono_nomix$celltype, prediction.assay = TRUE, 
                                  weight.reduction = Myeloid[["pca"]], dims = 1:30)
saveRDS(predictions.assay, "./11.ST_EOPC/5.Myeloid/PredictionAssay_Myeloid_celltype_SeuratAnchored.Rds")
Myeloid[["predictions"]] <- predictions.assay
saveRDS(Myeloid, "./11.ST_EOPC/5.Myeloid/Myeloid_celltype_SeuratAnchored.Rds")

#DefaultAssay(Myeloid) <- "predictions"
#pdf("./12.ST_LOPC/3.anchoring/anchor_ST_celltype_Seurat.pdf", width = 12, height = 12)
#SpatialFeaturePlot(Myeloid, rownames(Myeloid), pt.size.factor = 2, ncol = 4, crop = TRUE)
#dev.off()

#Boxplot

cluster <- Myeloid$SCT_snn_res.0.8 %>% as.data.frame()
colnames(cluster) <- "cluster"

mtx <- cbind(t(as.data.frame(predictions.assay@data)), cluster)
mtx <- mtx[,-13]
mtx <- melt(mtx)
mtx$variable <- gsub("-","_",mtx$variable)
mtx$variable <- factor(mtx$variable, levels = c('pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
                                                'Mono_IL1B', 'Mono_FCN1',
                                                'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                                'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A'))

P = ggplot(mtx, aes(x = cluster, y = value, fill = variable)) +
  scale_fill_manual(values = paste0(color, "CC")) +
  geom_boxplot(width = 0.8, outlier.colour = NA, position=position_dodge(0.8), notch = T) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("Predict Score") +
  xlab("") 
pdf("./11.ST_EOPC/5.Myeloid/BoxPlot_SeuratPredict_cluster_split.pdf", width = 7, height = 3)
P
dev.off()

## Myeloid score ----------------------------------------------

Mye_markers <- readRDS("./5.4.Mono_harmony/Markers/Nomix_celltype_markers_nomix.Rds")
Mye_markers <- Mye_markers[which(Mye_markers$avg_log2FC > 0.25 & Mye_markers$p_val_adj < 0.01),]

Mye_gs <- list()
for (i in as.character(unique(Mye_markers$cluster))) {
  Mye_gs[[i]] <- Mye_markers$gene[Mye_markers$cluster == i]
}

expMtx <- as.matrix(Myeloid@assays$Spatial@data)
cells_rankings <- AUCell_buildRankings(expMtx, nCores=10, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(Mye_gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./11.ST_EOPC/5.Myeloid/Myescore_AUC.Rds")

Myeloid[["Mye_AUCell"]] <- CreateAssayObject(counts = cells_AUC@assays@data$AUC)
DefaultAssay(Myeloid) <- "Mye_AUCell"

# boxplot 

AUCell <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()
cluster <- Myeloid$SCT_snn_res.0.8 %>% as.data.frame()
colnames(cluster) <- "cluster"

mtx <- cbind(AUCell, cluster)
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = cluster, y = value, fill = variable)) +
  scale_fill_manual(values = paste0(color, "CC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("AUCell score") +
  xlab("") 
pdf("./11.ST_EOPC/5.Myeloid/BoxPlot_MyeScore_cluster_split.pdf", width = 8, height = 3)
P
dev.off()

## Markers annotate --------------------------------------------------

DefaultAssay(Myeloid) <- "Spatial"
Features <- c("JCHAIN","IGLC2","IGHA1",
              "CLEC9A","WDFY4","CPNE3",
              "CD1C","FCER1A","CD1E",
              "LAMP3","IDO1","BIRC3",
              "IL1B","EREG","CCL20",
              "FCN1","S100A8","S100A9",
              "CXCL11","CXCL10","GBP1",
              "C3","LPAR6","MAF",
              "RGS1","JUN","FOS",
              "FOLR2","SELENOP","MRC1",
              "APOE","CTSD","APOC1",
              "MT1A","MT1G","MT1H",
              "FTL","KLK3","MSMB")

library(scRNAtoolVis)

pdf("./11.ST_EOPC/5.Myeloid/MarkerBubble_cluster.pdf",width = 8, height = 6)
jjDotPlot(Myeloid,
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

## annotate --------------------------------------------------

celltype <- data.frame(cluster = levels(Myeloid@meta.data$SCT_snn_res.0.8), celltype = "Macro_APOE")
celltype$celltype[celltype$cluster %in% c(0)] <- "pDC"

Myeloid@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  Myeloid@meta.data[which(Myeloid@meta.data$SCT_snn_res.0.8 == celltype$cluster[i]),'celltype'] <- celltype$celltype[i]}
table(Myeloid@meta.data$celltype)

Idents(Myeloid) <- "celltype"

saveRDS(Myeloid, "./11.ST_EOPC/5.Myeloid/Myeloid_annotated.Rds")

p1 <- DimPlot(Myeloid, label = TRUE, shuffle = T, pt.size = 1, label.size = 4, label.box = T) + 
  scale_color_manual(values = color[c(3,4)]) + 
  scale_fill_manual(values = color[c(3,4)]) + NoLegend()
p2 <- SpatialDimPlot(Myeloid, label = TRUE, label.size = 3, pt.size.factor = 2) + 
  scale_fill_manual(values = color[c(3,4)]) + NoLegend()

pdf("./11.ST_EOPC/5.Myeloid/celltype_umap_plot.pdf", height = 5, width = 10)
wrap_plots(p1, p2)
dev.off()

## macro score -------------------------------------------

Macrp_gs <- read.table("Macrophage_gs.txt", header = T, sep = "\t")

gs <- list()
for (i in c("M1","M2","MDSC","Angiogenesis","Phagocytosis")) {
  gs[[i]] <- Macrp_gs[Macrp_gs$term == i,]$gene
}

library(AUCell)

expMtx <- GetAssayData(Myeloid, assay = "Spatial", slot = "data")
cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./11.ST_EOPC/5.Myeloid/Macro_gs/Macro_gs_AUC.Rds")

AUCell <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

celltype <- Myeloid$celltype %>% as.character() %>% as.data.frame()
colnames(celltype) <- "celltype"
rownames(celltype) <- Cells(Myeloid)

mtx <- cbind(AUCell, celltype)

mtx <- split(mtx, mtx$celltype)
mean <- c()
for (i in 1:length(mtx)) {
  mean_tmp <- mtx[[i]][,-grep("celltype",colnames(mtx[[i]]))] %>% colMeans()
  mean <- rbind(mean, mean_tmp)
  rownames(mean)[i] <- names(mtx)[i]
}

## pheatmap

pdf("./11.ST_EOPC/5.Myeloid/Macro_gs/Macro_gs_heatmap.pdf", width = 3, height = 1.5)
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
         main = "",
         angle_col = "45",
         border_color = "grey",
         display_numbers = T,
         number_color = "black",
         legend = F
)
dev.off()

pdf("./11.ST_EOPC/5.Myeloid/Macro_gs/Macro_gs_heatmap_rotate.pdf", width = 2, height = 4)
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
pheatmap(t(mean),
         cluster_cols = F,
         cluster_rows = F,
         color = c(colorRampPalette(colors = c("#053E7A","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A52A2A"))(length(bk)/2)),
         show_colnames = T,
         show_rownames = T,
         scale="row",  #矫正
         breaks=bk,
         legend_breaks=seq(-2,2,1),
         fontsize = 8,
         main = "",
         angle_col = "45",
         border_color = "grey",
         display_numbers = T,
         number_color = "black",
         legend = F
)
dev.off()

# Boxplot

mtx <- cbind(AUCell, celltype)
mtx <- melt(mtx)

pdf("./11.ST_EOPC/5.Myeloid/Macro_gs/Boxplot_macro_phenotype.pdf")

mtx$celltype <- factor(mtx$celltype, levels = c("pDC","Macro_APOE"))

ggplot(mtx, aes(x = variable, y = value, fill = celltype)) +
  scale_fill_manual(values = paste0(color, "CC")) +
  geom_boxplot(width = 0.8, outlier.colour = NA, position=position_dodge(0.9)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("AUCell Score") +
  xlab("") +
  scale_fill_manual(values = color[c(3,4)]) +
  rotate_x_text(angle = 45) +
  stat_compare_means(aes(group=celltype),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","")),
                     label = "p.signif")

dev.off()

## FA and cholesterol score -----------------------------------

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

expMtx <- GetAssayData(Myeloid, assay = "Spatial", slot = "data")
GSVA <- gsva(expMtx, gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "11.ST_EOPC/5.Myeloid/Lipid_5gs_GSVA.Rds")

Myeloid[["lipid"]] <- CreateAssayObject(counts = GSVA)
DefaultAssay(Myeloid) <- "lipid"

pdf("11.ST_EOPC/5.Myeloid/lipid/Vlnplot_5gs.pdf")
VlnPlot(Myeloid, features = rownames(Myeloid), cols = color, pt.size = 0)
dev.off()

pdf("11.ST_EOPC/5.Myeloid/lipid/Boxplot_5gs.pdf")

celltype <- Myeloid$celltype %>% as.data.frame()
colnames(celltype) <- "celltype"

mtx <- cbind(t(as.data.frame(Myeloid@assays$lipid@data)), celltype)
mtx <- melt(mtx)
mtx$celltype <- factor(mtx$celltype, levels = c("pDC","Macro_APOE"))

ggplot(mtx, aes(x = variable, y = value, fill = celltype)) +
  scale_fill_manual(values = paste0(color, "CC")) +
  geom_boxplot(width = 0.8, outlier.colour = NA, position=position_dodge(0.9), notch = T) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("GSVA Score") +
  xlab("") +
  scale_fill_manual(values = color[c(3,4)]) +
  rotate_x_text(angle = 45) +
  stat_compare_means(aes(group=celltype),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","")),
                     label = "p.signif")
  
dev.off()

pdf("11.ST_EOPC/5.Myeloid/lipid/Spatialplot_5gs.pdf")
SpatialFeaturePlot(Myeloid, features = rownames(Myeloid), pt.size.factor = 2)
dev.off()

## SPP1 and MDSC ------------------------------------------

P1 <- SpatialFeaturePlot(Myeloid, features=c("APOE","CD163","SPP1","CD84"), ncol=4)
P2 <- VlnPlot(Myeloid, features=c("APOE","CD163","SPP1","CD84"), pt.size = 1, cols = color[c(3,4)], ncol=4)

pdf("11.ST_EOPC/5.Myeloid/APOE_CD163_SPP1_CD84_spatial.pdf", height = 4, width = 46)
P1
dev.off()

pdf("11.ST_EOPC/5.Myeloid/APOE_CD163_SPP1_CD84_exp.pdf", height = 5, width = 16)
P2
dev.off()

wilcox.test(Myeloid@assays$Spatial@data["APOE",colnames(Myeloid)[which(Myeloid$celltype == "Macro_APOE")]],
            Myeloid@assays$Spatial@data["APOE",colnames(Myeloid)[which(Myeloid$celltype == "pDC")]])
wilcox.test(Myeloid@assays$Spatial@data["CD163",colnames(Myeloid)[which(Myeloid$celltype == "Macro_APOE")]],
            Myeloid@assays$Spatial@data["CD163",colnames(Myeloid)[which(Myeloid$celltype == "pDC")]])
wilcox.test(Myeloid@assays$Spatial@data["CD36",colnames(Myeloid)[which(Myeloid$celltype == "Macro_APOE")]],
            Myeloid@assays$Spatial@data["CD36",colnames(Myeloid)[which(Myeloid$celltype == "pDC")]])
wilcox.test(Myeloid@assays$Spatial@data["SPP1",colnames(Myeloid)[which(Myeloid$celltype == "Macro_APOE")]],
            Myeloid@assays$Spatial@data["SPP1",colnames(Myeloid)[which(Myeloid$celltype == "pDC")]])
wilcox.test(Myeloid@assays$Spatial@data["CD84",colnames(Myeloid)[which(Myeloid$celltype == "Macro_APOE")]],
            Myeloid@assays$Spatial@data["CD84",colnames(Myeloid)[which(Myeloid$celltype == "pDC")]])


# 1.3.7.Epi and Myeloid -----------------------------------------

Epi_Mye <- subset(EOPC4_ST, celltype %in% c("Cancer-1","Cancer-2","Cancer-3","Myeloid"))

pdf("11.ST_EOPC/Epi_Mye.pdf", width = 3, height = 3)
SpatialDimPlot(Epi_Mye, label = T, pt.size.factor = 1.5, label.size = 3) + 
  scale_fill_manual(values = color[c(1,5,3,2)]) + NoLegend()
dev.off()

# 2.1.Data processing - LOPC-3-K -------------------------------------------------------

# load data from "filtered_feature_bc_matrix.h5" & "spatial"

LOPC3_ST <- Load10X_Spatial(data.dir = "DATA/PCSC010-K/", 
                            assay = "Spatial", 
                            slice = "slice", 
                            filter.matrix = TRUE)

LOPC3_ST@meta.data$orig.ident <- "LOPC"

# 2.2.ST ---------------------------------------------------------------------------------

# Data visualization 

## Data visualization ----------------------------

plot1 <- VlnPlot(LOPC3_ST, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(LOPC3_ST, features = "nCount_Spatial") + theme(legend.position = "right")

pdf("./12.ST_LOPC/1.QC/nCounts.pdf", height = 5, width = 12)
wrap_plots(plot1, plot2)
dev.off()

plot1 <- VlnPlot(LOPC3_ST, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(LOPC3_ST, features = "nFeature_Spatial") + theme(legend.position = "right")

pdf("./12.ST_LOPC/1.QC/nFeature.pdf", height = 5, width = 12)
wrap_plots(plot1, plot2)
dev.off()

## Normalization -----------------------------

LOPC3_ST[["percent.mt"]] <- PercentageFeatureSet(LOPC3_ST, pattern = "^MT-")
LOPC3_ST <- NormalizeData(LOPC3_ST, verbose = TRUE)
str(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
LOPC3_ST <- CellCycleScoring(LOPC3_ST, g2m.features=g2m.genes, s.features=s.genes)
LOPC3_ST <- SCTransform(LOPC3_ST, assay = "Spatial", vars.to.regress = c("percent.mt"))

## PCA ---------------------------------------

LOPC3_ST <- RunPCA(LOPC3_ST, verbose = FALSE)   %>%
  FindNeighbors(reduction = "pca", dims = 1:30) 

for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 0.9,1,1.5,2.0)) {
  LOPC3_ST <- FindClusters(LOPC3_ST, resolution = res)}

apply(LOPC3_ST@meta.data[, grep("SCT_snn_res", colnames(LOPC3_ST@meta.data))],2,table)

library(clustree)

p2_tree=clustree(LOPC3_ST@meta.data, prefix = "SCT_snn_res.") 
ggsave(plot=p2_tree, filename="./12.ST_LOPC/2.clustering/Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(LOPC3_ST) <- "SCT_snn_res.2"

#tsne & UMAP

LOPC3_ST <- RunUMAP(LOPC3_ST, reduction = "pca", dims = 1:30)
saveRDS(LOPC3_ST, "12.ST_LOPC/2.clustering/LOPC3_ST_clustered.Rds")

#plot

p1 <- DimPlot(LOPC3_ST, label = TRUE, shuffle = T, pt.size = 1, label.size = 4, label.box = T) + NoLegend()
p2 <- SpatialDimPlot(LOPC3_ST, label = TRUE, label.size = 3, pt.size.factor = 2) + NoLegend()

pdf("./12.ST_LOPC/2.clustering/umap_plot.pdf", height = 5, width = 10)
wrap_plots(p1, p2)
dev.off()

pdf("./12.ST_LOPC/2.clustering/spatial_plot_split_by_clusters.pdf", width = 14, height = 14)
p4 <- SpatialDimPlot(LOPC3_ST, 
                     cells.highlight = CellsByIdentities(object = LOPC3_ST, idents = levels(LOPC3_ST$SCT_snn_res.2)),
                     facet.highlight = TRUE, ncol = 5) + NoLegend()
p4
dev.off()

#interactive plot

#SpatialDimPlot(LOPC3_ST, interactive = TRUE)
#SpatialFeaturePlot(LOPC3_ST, features = c("EPCAM", "MLH1", "MSH2", 
#                                    "MSH6", "PMS2"))
#LinkedDimPlot(LOPC3_ST)

#pdf("./2.clustering/cluster0.pdf", height = 5, width = 10)
#LinkedDimPlot(LOPC3_ST)
#dev.off()

## Find spatially variable genes ---------------------------------------

LOPC3_ST <- FindSpatiallyVariableFeatures(LOPC3_ST, assay = "SCT",
                                          features = VariableFeatures(LOPC3_ST)[1:1000], 
                                          selection.method = "markvariogram")
saveRDS(LOPC3_ST, "./12.ST_LOPC/2.clustering/EOPC4_clustered.Rds")

top.features <- head(SpatiallyVariableFeatures(LOPC3_ST, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(LOPC3_ST, features = top.features, ncol = 3, alpha = c(0.1, 1))

SpatialFeaturePlot(LOPC3_ST, features = c("LYZ"), interactive = FALSE)

## cluster markers --------------------------------

library(SCP)

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)

Allmarkers <- FindAllMarkers(LOPC3_ST,
                             min.pct = 0.1, 
                             logfc.threshold = 0,
                             verbose = T)
saveRDS(Allmarkers, "./12.ST_LOPC/2.clustering/cellcluster_markers.Rds")

sigmarkers <- Allmarkers[Allmarkers$p_val_adj < 0.01 & Allmarkers$avg_log2FC > 0,]
topmarkers <- Allmarkers %>% group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

# 2.3.1.anchored (Seurat) ---------------------------------

anchors <- FindTransferAnchors(reference = PCSC_RandomCell, query = LOPC3_ST, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = PCSC_RandomCell$inferCNVCellStat, prediction.assay = TRUE, 
                                  weight.reduction = LOPC3_ST[["pca"]], dims = 1:30)
LOPC3_ST[["predictions"]] <- predictions.assay
saveRDS(LOPC3_ST, "./12.ST_LOPC/3.anchoring/LOPC3_ST_celltype_SeuratAnchored.Rds")

DefaultAssay(LOPC3_ST) <- "predictions"
pdf("./12.ST_LOPC/3.anchoring/anchor_ST_celltype_Seurat.pdf", width = 12, height = 12)
SpatialFeaturePlot(LOPC3_ST, rownames(LOPC3_ST), pt.size.factor = 2, ncol = 4, crop = TRUE)
dev.off()

# 2.3.2.anchored (spacexr) ---------------------------------

library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)

#prepare SC data
PCSC_RandomCell <- readRDS("PCSC_RandomCell_new.Rds")
counts <- as.matrix(PCSC_RandomCell@assays$RNA@counts)
meta_data <- PCSC_RandomCell@meta.data
cell_types <- meta_data$inferCNVCellStat
names(cell_types) <- rownames(meta_data)
cell_types <- as.factor(cell_types)
nUMI <- meta_data$nCount_RNA
names(nUMI) <- rownames(meta_data)
### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)

#prepare ST data
coords <- GetTissueCoordinates(LOPC3_ST, scale = NULL) 
spatialcounts <- as.matrix(LOPC3_ST@assays$Spatial@counts)
nUMI <- LOPC3_ST@meta.data$nCount_Spatial
names(nUMI) <- rownames(LOPC3_ST@meta.data)
### Create SpatialRNA object
puck <- SpatialRNA(coords, spatialcounts, nUMI)

#RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') #doublet_mode = 'doublet'
saveRDS(myRCTD, "./12.ST_LOPC/3.anchoring/spacexr_anchored.Rds")

# 2.3.3. annotation ---------------------------------------------

Features=c("EPCAM","CDH1", 
           "AR", "AMACR", 
           "TP63", "KRT5",
           "PTPRC", 
           "CD3E", "CD3G", "CD2", 
           "MS4A1", "CD79A", 
           "CD14", "FCGR3A", "LYZ", 
           "KIT", "MS4A2", 
           "VIM", 
           "PECAM1", "ENG", "VWF", 
           "FAP", "PDGFRA", "THY1", 
           "COL3A1", "ACTA2", "MYL9")

DefaultAssay(LOPC3_ST) <- "Spatial"
pdf("./12.ST_LOPC/2.clustering/AllmarkerBubble.pdf",width = 10, height = 10)
jjDotPlot(LOPC3_ST,
          gene = Features,
          id = "SCT_snn_res.2",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 6,
          x.text.vjust = 0.5)
dev.off()

for (i in Features) {
  Scatter <- FeaturePlot(LOPC3_ST, features = i, cols = c("grey", "#A52A2A"),
                         pt.size = 0.3, ncol = ) + NoLegend()
  ggsave(plot = Scatter, device = "tiff", filename = paste0("./12.ST_LOPC/2.clustering/Scatter/", i, ".tiff"), width = 2, height = 2)
}


celltype <- data.frame(cluster = levels(LOPC3_ST@meta.data$SCT_snn_res.2), celltype = "Cancer-1")
celltype$celltype[celltype$cluster %in% c(0,17,18)] <- "SMC"
celltype$celltype[celltype$cluster %in% c(1,2,6,20,22)] <- "Fibroblast"
celltype$celltype[celltype$cluster %in% c(3)] <- "Cancer-2"
celltype$celltype[celltype$cluster %in% c(4)] <- "B-infiltr.-Cancer"
celltype$celltype[celltype$cluster %in% c(11)] <- "B-cell"
celltype$celltype[celltype$cluster %in% c(7,8,14,21)] <- "Mixed-cell"
celltype$celltype[celltype$cluster %in% c(9,16)] <- "Cancer-3"
celltype$celltype[celltype$cluster %in% c(12,13)] <- "Myeloid"
celltype$celltype[celltype$cluster %in% c(19)] <- "Cancer-4"

LOPC3_ST@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  LOPC3_ST@meta.data[which(LOPC3_ST@meta.data$SCT_snn_res.2 == celltype$cluster[i]),'celltype'] <- celltype$celltype[i]}
table(LOPC3_ST@meta.data$celltype)

Idents(LOPC3_ST) <- "celltype"
Idents(LOPC3_ST) <- factor(Idents(LOPC3_ST), levels = c("Cancer-1", "Cancer-2", "Cancer-3", "Cancer-4", "Fibroblast", 
                                                        "Myeloid", "SMC", "B-cell", "B-infiltr.-Cancer", "Mixed-cell"))

saveRDS(LOPC3_ST, "./12.ST_LOPC/LOPC3_ST_annotated.Rds")

p1 <- DimPlot(LOPC3_ST, label = TRUE, shuffle = T, pt.size = 1, label.size = 4, label.box = T) + 
  scale_color_manual(values = color[c(1,5,3,14,7,2,4,6,11,12)]) + 
  scale_fill_manual(values = color[c(1,5,3,14,7,2,4,6,11,12)]) + NoLegend()
p2 <- SpatialDimPlot(LOPC3_ST, label = FALSE, pt.size.factor = 2.2) + 
  scale_fill_manual(values = color[c(1,5,3,14,7,2,4,6,11,12)]) + NoLegend()

pdf("./12.ST_LOPC/2.clustering/celltype_umap_plot.pdf", height = 5, width = 10)
wrap_plots(p1, p2)
dev.off()

pdf("./12.ST_LOPC/2.clustering/celltype_Spatial_plot.pdf", height = 4, width = 4)
SpatialDimPlot(LOPC3_ST, label = FALSE, pt.size.factor = 2) + 
  scale_fill_manual(values = color[c(1,5,3,14,7,2,4,6,11,12)])
dev.off()


## Dotplot

LOPC3_ST$celltype <- factor(LOPC3_ST$celltype, levels = 
                              c("Cancer-1","Cancer-2","Cancer-3","Cancer-4",
                                "B-infiltr.-Cancer","B-cell","Myeloid",
                                "Fibroblast","SMC","Mixed-cell"))

pdf("./12.ST_LOPC/2.clustering/Celltype_AllmarkerBubble.pdf",width = 10, height = 8)
jjDotPlot(LOPC3_ST,
          gene = Features,
          id = "celltype",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 6,
          x.text.vjust = 0.5)
dev.off()

# 2.3.4.MP1 score -----------------------------------------

Epi <- subset(LOPC3_ST, celltype %in% c("Cancer-1","Cancer-2","Cancer-3","Cancer-4","B-infiltr.-Cancer"))

pdf("./12.ST_LOPC/4.Epi/celltype_Spatial_plot.pdf", height = 4, width = 4)
SpatialDimPlot(Epi, label = FALSE, pt.size.factor = 2) + 
  scale_fill_manual(values = color[c(1,5,3,14,11)])
dev.off()

celltype = unique(as.character(Epi$celltype))
for (i in celltype) {
  P <- SpatialDimPlot(subset(Epi, celltype == i), label = FALSE, pt.size.factor = 2) + 
    scale_fill_manual(values = color[c(1,5,3,14,11)][grep(i, celltype)]) + NoLegend()
  ggsave(paste0("./12.ST_LOPC/4.Epi/", celltype, "_Spatial_plot.pdf"), P, height = 3, width = 3)
}

cols <- color[c(1,5,3,14,11)]

for (i in 1:length(cols)) {
  
  P <- SpatialDimPlot(Epi, cols = cols[i], cols.highlight = c(cols[i], "grey50"), pt.size.factor = 2.2,
                      cells.highlight = CellsByIdentities(object = Epi, idents = levels(Epi$celltype)),
                      facet.highlight = TRUE, ncol = 3) + NoLegend()
  ggsave(paste0("./12.ST_LOPC/4.Epi/", levels(Epi$celltype)[i], "_Spatial_plot.pdf"), P, height = 9, width = 6)
  
}


library(AUCell)
Markers <- readRDS("./5.1.1.epithelial_NMF/process_k3/MP_markers.Rds")
expMtx <- as.matrix(Epi@assays$Spatial@data)
cells_rankings <- AUCell_buildRankings(expMtx, nCores=10, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./12.ST_LOPC/4.Epi/MPScore_AUC.Rds")

Epi[["MP_AUCell"]] <- CreateAssayObject(counts = cells_AUC@assays@data$AUC)
DefaultAssay(Epi) <- "MP_AUCell"

pdf("./12.ST_LOPC/4.Epi/SpatialPlot_MP1_score.pdf")
SpatialFeaturePlot(Epi, features = "MP1", pt.size.factor = 2.2)
dev.off()

#boxplot

exp <- data.frame(celltype=Epi$celltype, MP1=Epi@assays$MP_AUCell@data["MP1",])
rownames(exp) <- colnames(Epi)
exp$celltype <- factor(exp$celltype, levels = c("Cancer-1","Cancer-2","Cancer-3","Cancer-4","B-infiltr.-Cancer"))
p=ggboxplot(exp, x="celltype", y="MP1",fill= "celltype",
            ylab="MP1 AUCell score",
            xlab="",
            palette = color[c(1,5,3,14,7)]) +
  ylim(c(0,0.5))
p=p + rotate_x_text(45) + theme(axis.title.y = element_text(vjust=2, size=18))
p1=p+stat_compare_means(aes(group=celltype),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","")),
                        label = "p.signif")
p1

pdf("./12.ST_LOPC/4.Epi/BoxPlot_MP1Score_CancerGroup.pdf", width = 3, height = 5)
p1 + NoLegend()
dev.off()

wilcox.test(exp$MP1[exp$celltype == "Cancer-1"], exp$MP1[exp$celltype == "Cancer-2"])
wilcox.test(exp$MP1[exp$celltype == "Cancer-3"], exp$MP1[exp$celltype == "Cancer-1"])
wilcox.test(exp$MP1[exp$celltype == "Cancer-2"], exp$MP1[exp$celltype == "Cancer-3"])
wilcox.test(exp$MP1[exp$celltype == "Cancer-4"], exp$MP1[exp$celltype == "B-infiltr.-Cancer"])

## BMPRs --------------------------------------

## BMPs ----------------------------------------------

BMPRs <- Epi@assays$Spatial@data[c("BMPR1B","BMPR2"),] %>% as.matrix() %>% t() %>% as.data.frame()

celltype <- Epi$celltype %>% as.data.frame()
colnames(celltype) <- "celltype"

mtx <- cbind(BMPRs, celltype) %>% as.data.frame()
mtx$BMPRs <- (mtx$BMPR1B + mtx$BMPR2) / 2
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-1"], mtx$BMPRs[mtx$celltype=="Cancer-2"])  #ns
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-1"], mtx$BMPRs[mtx$celltype=="Cancer-3"])  #***
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-1"], mtx$BMPRs[mtx$celltype=="Cancer-4"])  #**
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-1"], mtx$BMPRs[mtx$celltype=="Cancer-4"])  #**
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-1"], mtx$BMPRs[mtx$celltype=="B-infiltr.-Cancer"])  #**
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-2"], mtx$BMPRs[mtx$celltype=="Cancer-3"])  #**
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-2"], mtx$BMPRs[mtx$celltype=="Cancer-4"])  #**
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-2"], mtx$BMPRs[mtx$celltype=="B-infiltr.-Cancer"])  #**
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-3"], mtx$BMPRs[mtx$celltype=="Cancer-4"])  #**
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-3"], mtx$BMPRs[mtx$celltype=="B-infiltr.-Cancer"])  #**
wilcox.test(mtx$BMPRs[mtx$celltype=="Cancer-4"], mtx$BMPRs[mtx$celltype=="B-infiltr.-Cancer"])  #**

BMPRs_exp <- t(mtx[,-3])
Epi[["BMPRs"]] <- CreateAssayObject(counts = BMPRs_exp)
DefaultAssay(Epi) <- "BMPRs"

Epi$celltype <- factor(Epi$celltype, levels = c("Cancer-1","Cancer-2","Cancer-3","Cancer-4","B-infiltr.-Cancer"))
P = VlnPlot(Epi, "BMPRs", group.by="celltype", cols = color[c(1,5,3,14,7)], pt.size = 0.05) + 
  NoLegend() +
  xlab("") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=8, colour="black"), 
        axis.title.y=element_text(size = 10))

pdf("./12.ST_LOPC/4.Epi/BMPRs_exp.pdf", width = 2.5, height = 3.2)
P
dev.off()

# 2.4.Fibroblast -----------------------------------------------

Fibr <- subset(LOPC3_ST, celltype == "Fibroblast")
Fibr <- SCTransform(Fibr, assay = "Spatial", vars.to.regress = "percent.mt")

## PCA ---------------------------------------

Fibr <- RunPCA(Fibr, verbose = FALSE)   %>%
  FindNeighbors(reduction = "pca", dims = 1:30) 

for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 0.9,1,1.5,2)) {
  Fibr <- FindClusters(Fibr, resolution = res)}

apply(Fibr@meta.data[, grep("SCT_snn_res", colnames(Fibr@meta.data))],2,table)

library(clustree)

p2_tree=clustree(Fibr@meta.data, prefix = "SCT_snn_res.") 
ggsave(plot=p2_tree, filename="./12.ST_LOPC/5.Fibr/Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(Fibr) <- "SCT_snn_res.0.3"

#tsne & UMAP

Fibr <- RunUMAP(Fibr, reduction = "pca", dims = 1:30)
saveRDS(Fibr, "12.ST_LOPC/5.Fibr/Fibr_clustered.Rds")

#plot

p1 <- DimPlot(Fibr, label = TRUE, shuffle = T, pt.size = 1, label.size = 4, label.box = T) + 
  scale_color_manual(values = color) + 
  scale_fill_manual(values = color) + NoLegend()
p2 <- SpatialDimPlot(Fibr, label = TRUE, label.size = 3, pt.size.factor = 2) + 
  scale_fill_manual(values = color) + NoLegend()

pdf("./12.ST_LOPC/5.Fibr/umap_plot.pdf", height = 5, width = 10)
wrap_plots(p1, p2)
dev.off()

## Seurat anchoring ----------------------------------------------

Fibr_nomix <- readRDS("5.5.Fibr_harmony/Nomix_Fibr_annoted.Rds")
Fibr_nomix$celltype <- as.character(Fibr_nomix$celltype)
Fibr_nomix <- Fibr_nomix %>% subset(celltype %in% grep("CAF", unique(Fibr_nomix$celltype), value = T))
#Fibr_nomix <- SCTransform(Fibr_nomix, vars.to.regress = "percent.mt")
#DefaultAssay(Fibr_nomix) <- "RNA"
DefaultAssay(Fibr) <- "Spatial"

anchors <- FindTransferAnchors(reference = Fibr_nomix, query = Fibr, normalization.method = "LogNormalize")
predictions.assay <- TransferData(anchorset = anchors, refdata = Fibr_nomix$celltype, prediction.assay = TRUE, 
                                  weight.reduction = Fibr[["pca"]], dims = 1:30)
Fibr[["predictions"]] <- predictions.assay
saveRDS(Fibr, "./12.ST_LOPC/5.Fibr/Fibr_celltype_SeuratAnchored.Rds")

#DefaultAssay(Myeloid) <- "predictions"
#pdf("./12.ST_LOPC/3.anchoring/anchor_ST_celltype_Seurat.pdf", width = 12, height = 12)
#SpatialFeaturePlot(Myeloid, rownames(Myeloid), pt.size.factor = 2, ncol = 4, crop = TRUE)
#dev.off()

#Boxplot

cluster <- Fibr$SCT_snn_res.0.3 %>% as.data.frame()
colnames(cluster) <- "cluster"

mtx <- cbind(t(as.data.frame(predictions.assay@data)), cluster)
mtx <- mtx[,-3]
mtx <- melt(mtx)

P = ggplot(mtx, aes(x = cluster, y = value, fill = variable)) +
  scale_fill_manual(values = paste0(color, "CC")) +
  geom_boxplot(width = 0.8, outlier.colour = NA, position=position_dodge(0.8), notch = T) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("Predict Score") +
  xlab("") 

## BMPs ----------------------------------------------

BMPs <- Fibr@assays$Spatial@data[c("BMP5","BMP7","BMP4"),] %>% as.matrix() %>% t() %>% as.data.frame()

Fibr$celltype <- "Cancer-2-assoc-CAF"
Fibr$celltype[which(Fibr$SCT_snn_res.0.3 == 1)] <- "Other-CAF"

Idents(Fibr) <- "celltype"

p1 <- DimPlot(Fibr, label = TRUE, shuffle = T, pt.size = 1, label.size = 4, label.box = T) + 
  scale_color_manual(values = color[c(7,12)]) + 
  scale_fill_manual(values = color[c(7,12)]) + NoLegend()
p2 <- SpatialDimPlot(Fibr, label = FALSE, pt.size.factor = 2) + 
  scale_fill_manual(values = color[c(7,12)]) + NoLegend()

pdf("./12.ST_LOPC/5.Fibr/celltype_plot.pdf", height = 5, width = 10)
wrap_plots(p1, p2)
dev.off()

celltype <- Fibr$celltype %>% as.data.frame()
colnames(celltype) <- "celltype"

mtx <- cbind(BMPs, celltype) %>% as.data.frame()
wilcox.test(mtx$BMP5[mtx$celltype=="Cancer-3-assoc-CAF"], mtx$BMP5[mtx$celltype=="Other-CAF"])  #ns
wilcox.test(mtx$BMP7[mtx$celltype=="Cancer-3-assoc-CAF"], mtx$BMP7[mtx$celltype=="Other-CAF"])  #***
wilcox.test(mtx$BMP4[mtx$celltype=="Cancer-3-assoc-CAF"], mtx$BMP4[mtx$celltype=="Other-CAF"])  #**

BMPs_P <- list()
for (i in c("BMP5","BMP7","BMP4")) {
  P = VlnPlot(Fibr, i, group.by="celltype", cols = color[c(7,12)], pt.size = 0.05) + 
    NoLegend() +
    xlab("") +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
          axis.text.y=element_text(size=8, colour="black"), 
          axis.title.y=element_text(size = 10))
  BMPs_P[[i]]<- P
}

library(gridExtra)

pdf("./12.ST_LOPC/5.Fibr/BMPs_exp.pdf", width = 4, height = 3)
grid.arrange(grobs = BMPs_P, ncol=3)
dev.off()
