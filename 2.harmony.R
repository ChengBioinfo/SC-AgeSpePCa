library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library("SingleR")
library(limma)
library(SingleR)
library(monocle)
library(tidyverse)
library(RCurl)
library(ggsci)
library(gridExtra)
library(harmony)
library(scRNAtoolVis)

color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")
           
# Integrating data (harmony) ------------------------------------------------------------------

PCSC_preInteg <- readRDS("./1.data_processing/QC/PCSC_QC.Rds")

folder <- "./2.1.harmony"
if(!dir.exists(folder)){
  dir.create(folder)
}

# 1.Cell cycle scoring ---------------------------------------------------------

folder <- "./2.1.harmony/cellcycle"
if(!dir.exists(folder)){
  dir.create(folder)
}

PCSC_preInteg <- NormalizeData(PCSC_preInteg)
str(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
PCSC_preInteg <- CellCycleScoring(PCSC_preInteg,
                               g2m.features = g2m.genes,
                               s.features = s.genes)
View(PCSC_preInteg@meta.data)

PCSC_preInteg <- FindVariableFeatures(PCSC_preInteg,
                                      selection.method = "vst",
                                      nfeatures = 2000,
                                      verbose = FALSE) %>%
              ScaleData(vars.to.regress = c("percent.mt"))
PCSC_preInteg <- RunPCA(PCSC_preInteg, npcs = 50)

ElbowPlot(PCSC_preInteg, ndims=50)
PCSC_preInteg <- RunUMAP(PCSC_preInteg, dims = 1:22)

pdf("./2.1.harmony/umap_preInteg.pdf", width = 7.5, height = 7)
DimPlot(PCSC_preInteg,
        reduction = "umap",
        group.by= "orig.ident",
        cols = color)
dev.off()

png("./2.1.harmony/umap_preInteg.png", width = 620, height = 480, res = 144)
DimPlot(PCSC_preInteg,
        reduction = "umap",
        group.by= "orig.ident",
        cols = color) +
  theme(plot.title = element_blank())
dev.off()

pdf("./2.1.harmony/cellcycle/cellcycle.pdf",width = 12, height = 4)
DimPlot(PCSC_preInteg,
        reduction = "umap",
        group.by= "Phase",
        split.by = "Phase",
        cols = color)
dev.off()

pdf("./2.1.harmony/cellcycle/cellcycle_SampleSplit.pdf",width = 20, height = 8)
DimPlot(PCSC_preInteg, 
        reduction = "pca", 
        group.by = "orig.ident",
        split.by = "Phase",
        shuffle = T,
        ncol = 5,
        cols = color)
dev.off()

# 2.harmony --------------------------------------------------------

pdf("./2.1.harmony/ElbowPlot_harmony.pdf")
PCSC <- RunHarmony(PCSC_preInteg, "orig.ident", plot_convergence = T)
dev.off()

PCSC <- RunUMAP(PCSC, reduction = "harmony", dims = 1:20)

# 3.clustering -------------------------------------------------------

PCSC <- FindNeighbors(PCSC, reduction = "harmony", 
                      dims = 1:20 )
for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1)) {
  PCSC <- FindClusters(PCSC,  
                       resolution = res)}

apply(PCSC@meta.data[, grep("RNA_snn_res", colnames(PCSC@meta.data))],2,table)

library(clustree)

p2_tree=clustree(PCSC@meta.data, prefix = "RNA_snn_res.") +
  scale_color_npg() + scale_fill_npg()
ggsave(plot=p2_tree, filename="./2.1.harmony/Tree_diff_resolution.pdf",width = 7, height = 12)
p2_tree

Idents(PCSC) <- "RNA_snn_res.0.2"
PCSC@meta.data$orig.ident <- PCSC_preInteg@meta.data$orig.ident

pdf("./2.1.harmony/CellCluster-UMAPPlot_res0.2.pdf",width = 7,height = 7)
DimPlot(PCSC, 
        reduction = "umap", 
        pt.size = 0.1, 
        label = T,
        cols = color) + 
  NoLegend() +
  ggtitle("")
dev.off()

png("./2.1.harmony/CellCluster-UMAPPlot_res0.2.png", width = 480,height = 480, res = 144)
DimPlot(PCSC, 
        reduction = "umap", 
        pt.size = 0.1, 
        label = T,
        cols = color) + 
  NoLegend() +
  theme(plot.title = element_blank())
dev.off()

png("./2.1.harmony/umap_postInteg.png", width = 620, height = 480, res = 144)
DimPlot(PCSC,
        reduction = "umap",
        group.by= "orig.ident",
        cols = color) +
  theme(plot.title = element_blank())
dev.off()

pdf("./2.1.harmony/CellCluster-UMAPPlot_cellcycle_SampleTypeSplit.pdf",width = 13,height = 7)
DimPlot(PCSC, 
        reduction = "umap", 
        pt.size = 0.1, 
        label = T, 
        group.by = "Phase", 
        split.by = "SampleType",
        cols = color) +
  ggtitle("")
dev.off()

pdf("./2.1.harmony/CellCluster-UMAPPlot_cellcycle_SampleType.pdf",width = 7,height = 7)
DimPlot(PCSC, 
        reduction = "umap", 
        pt.size = 0.1, 
        label = T, 
        group.by = "SampleType",
        cols = color) +
  ggtitle("")
dev.off()

pdf("./2.1.harmony/CellCluster-UMAPPlot_cellcycle_Sample.pdf",width = 7,height = 7)
DimPlot(PCSC, 
        reduction = "umap", 
        pt.size = 0.1, 
        label = T, 
        group.by = "orig.ident",
        cols = color) +
  ggtitle("")
dev.off()

# 4. annotation -------------------------------------------------------

# 4.1.SingleR annotation ------------------------------------------------------

refdata <- HumanPrimaryCellAtlasData()        #An external refdata is also acceptable 
#load("D:\\资料\\细胞注释库\\HumanPrimaryCellAtlasData.RData")
testdata <- GetAssayData(PCSC, slot="data")
#testdata <- PCSC@assays$SCT@data
clusters <- PCSC@meta.data$seurat_clusters

#按cluster进行注释------------------------------------------------------------------
cellpred <- SingleR(test = testdata, ref = refdata, 
                    labels = refdata$label.main, 
                    clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype <- data.frame(ClusterID=rownames(cellpred), 
                       celltype=cellpred$labels, stringsAsFactors = F)

PCSC@meta.data$singleR_celltype = "Unknown"
for(i in 1:nrow(celltype)){
  PCSC@meta.data[which(PCSC@meta.data$seurat_clusters == celltype$ClusterID[i]),'singleR_celltype'] <- celltype$celltype[i]
}

plot_anno <- DimPlot(PCSC, group.by="singleR_celltype", pt.size=0.1, label=TRUE, label.size=5, reduction='umap', shuffle = T, cols = color) +
  NoLegend()
ggsave(filename = "./4.annotation/CellCluster-UMAPPlot_res0.2_SingleR.pdf",
       plot = plot_anno)

png("./2.1.harmony/CellCluster-UMAPPlot_res0.2_SingleR.png",width = 480,height = 480, res = 144)
DimPlot(PCSC, 
        reduction = "umap", 
        group.by="singleR_celltype",
        pt.size = 0.1, 
        label = T,
        cols = color,
        repel = T,
        label.size = 3) + 
  NoLegend() +
  theme(plot.title = element_blank())
dev.off()

#按单个细胞进行注释---------------------------------------------------------------------
cellpred1 <- SingleR(test = testdata, ref = refdata, 
                     labels = refdata$label.main, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
PCSC[["SingleR.labels"]]<- as.character(cellpred1$labels)

plot_anno <- DimPlot(PCSC, group.by="SingleR.labels", pt.size=0.5, label=FALSE, label.size=5, reduction='umap', 
                     shuffle = T, cols = color)
ggsave(filename = "./4.annotation/CellCluster-UMAPPlot_SingleR_cell_resolution.pdf",
       plot = plot_anno, width = 10, height = 7)

# 4.2.使用scType进行注释 -----------------------------------------------------------------

library("HGNChelper")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# __secondary----------------------------------------------------------
db_ = "./4.annotation/ScTypeDB_short.xlsx"
tissue = "Prostate_secondary" 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix

es.max = sctype_score(scRNAseqData = PCSC[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = NULL)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(PCSC@meta.data$RNA_snn_res.0.2), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(PCSC@meta.data[PCSC@meta.data$RNA_snn_res.0.2==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(PCSC@meta.data$RNA_snn_res.0.2==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Mixed"
print(sctype_scores[,1:3])

PCSC@meta.data$scType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  PCSC@meta.data$scType[PCSC@meta.data$RNA_snn_res.0.2 == j] = as.character(cl_type$type[1])
}

pdf("./2.1.harmony/CellCluster-UMAPPlot_res0.2_ScType.pdf", )
DimPlot(PCSC, reduction = "umap", label = TRUE, repel = TRUE, 
        pt.size=0.1, group.by = 'scType', shuffle = T, cols = color) + NoLegend()
dev.off()

png("./2.1.harmony/CellCluster-UMAPPlot_res0.2_ScType.png",width = 480,height = 480, res = 144)
DimPlot(PCSC, 
        reduction = "umap", 
        group.by="scType",
        pt.size = 0.1, 
        label = T,
        cols = color,
        repel = T,
        label.size = 3) + 
  NoLegend() +
  theme(plot.title = element_blank())
dev.off()

# 4.3rename clusters ----------------------------

celltype=data.frame(ClusterID=0:13,
                    celltype='Epithelial')   
celltype[celltype$ClusterID %in% c(0,13),2]='T-cell'  
celltype[celltype$ClusterID %in% c(9),2]='B-cell' 
celltype[celltype$ClusterID %in% c(4),2]='Myeloid' 
celltype[celltype$ClusterID %in% c(11),2]='Mast' 
celltype[celltype$ClusterID %in% c(2),2]='Endothelial'
celltype[celltype$ClusterID %in% c(8),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(6),2]='SMC'
head(celltype)

PCSC@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  PCSC@meta.data[which(PCSC@meta.data$RNA_snn_res.0.2 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(PCSC@meta.data$celltype)

Idents(PCSC) <- "celltype"

## pal for celltype
pal_for_celltype <- color[c(1:13)]
names(pal_for_celltype) <- c("Epithelial","Normal","Other_Epi","T-cell","B-cell",
                             "Myeloid","Mast","Endothelial","Fibroblast", "Non1", "Non2", "SMC")

pdf("./2.1.harmony/CellCluster-UMAPPlot_res0.2_celltype.pdf", width = 5, height = 5)
DimPlot(PCSC, 
        label = T, 
        pt.size=0.5, 
        shuffle = T, 
        label.color = "black", 
        label.box = T, 
        cols = pal_for_celltype, 
        label.size = 5) + 
  NoLegend()
dev.off()

pdf("./2.1.harmony/CellCluster-UMAPPlot_res0.2_celltype_legend.pdf", width = 5, height = 5)
DimPlot(PCSC, 
        label = T, 
        pt.size=0.5, 
        shuffle = T, 
        label.color = "black", 
        label.box = T, 
        cols = pal_for_celltype, 
        label.size = 5)
dev.off()

## Marker Scatter plot -----------------------------------------------

Features=c("CDH1", "CD3D", "CD79A", "LYZ", "TPSAB1", "PECAM1", "THY1", "FAP")
for (i in Features) {
  Scatter <- FeaturePlot(PCSC, features = i, cols = c("grey", "#A52A2A"),
                         pt.size = 0.01) + NoLegend() + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          plot.title = element_blank())
  ggsave(plot = Scatter, device = "png", filename = paste0("./2.1.harmony/Scatter/", i, ".png"), width = 2, height = 2)
}

## 针对不同细胞类型画marker的dotplot--------------------------------------------------------------------------

PCSC$celltype <- factor(PCSC$celltype, 
                        levels = c("Epithelial", "T-cell", "B-cell",
                                   "Myeloid", "Mast", "Endothelial",
                                   "Fibroblast", "SMC"))
Idents(PCSC) <- "celltype"

Allmarker=c("EPCAM","CDH1", "PTPRC", 
            "CD3E", "CD3G", "CD3D", "CD2", "MS4A1", "CD79A", "CD14", 
            "FCGR3A", "CD68", "CD163", "LYZ", 
            "KIT", "MS4A2", "TPSAB1", "TPSB2", "VIM", "PECAM1", "ENG", "VWF", 
            "FAP", "PDGFRA", "THY1", "COL3A1", "ACTA2", "MYL9")

pdf("./2.1.harmony/AllmarkerBubble.pdf",width = 10, height = 4.7)
jjDotPlot(PCSC,
          gene = Allmarker,
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

# 5.tSNE/UMAP plots grouped by samples -----------------------------------------------------------

pdf("./2.1.harmony/CellCluster-UMAPPlot_SamGroup.pdf")
DimPlot(object = PCSC, 
        group.by="orig.ident", 
        pt.size=0.5,reduction = "umap",
        shuffle = T, cols = color)
dev.off()

# 6.tSNE/UMAP plots slited by samples -----------------------------------------------------------

pdf(paste0("./2.1.harmony/CellCluster-UMAPPlot_SamSplit.pdf"),width = 20,height = 8)
DimPlot(object = PCSC, 
        split.by="orig.ident", 
        pt.size=0.1,reduction = "umap",
        shuffle = T,
        cols = color,
        ncol = 5,
        label = T) +
  NoLegend()
dev.off()

# 8.cell counts of each sample in clusters ----------------------------------------------------------
table(PCSC@meta.data$orig.ident)
View(PCSC@meta.data)
Sample_cluster <- table(PCSC@meta.data$orig.ident,PCSC@meta.data$RNA_snn_res.0.2)
write.table(table(PCSC@meta.data$orig.ident,PCSC@meta.data$RNA_snn_res.0.2),"./2.1.harmony/cell_cluster.txt",sep="\t",quote = F)

# 9.sample type-------------------------------------------------------------------

folder <- "2.1.harmony/Sampletype"
if(!dir.exists(folder)){
  dir.create(folder)
}

pdf("./2.1.harmony/DimPlot_Sampletype_umap_split.pdf", width = 14, height = 7.5)
DimPlot(object = PCSC, 
        split.by="SampleType", 
        pt.size=0.1,
        reduction = "umap",
        shuffle = T,
        label = T,
        cols = color) +
  NoLegend()
dev.off()

## distribution ---------------

#Roe

tab <- table(PCSC$celltype, PCSC$SampleType)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./2.1.harmony/celltype_Roe_SampleType.pdf", width = 3, height = 1.5)
pheatmap(t(Roe),
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(colors = c("white","#A52A2A"))(100),
         show_colnames = T,
         show_rownames = T,
         scale="none",
         fontsize = 8,
         fontsize_row = 12,
         fontsize_col = 8,
         main = "Ro/e",
         angle_col = "45",
         display_numbers = T,
         legend = F)
dev.off()

pdf("./2.1.harmony/celltype_Roe_SampleType_rev.pdf", width = 2, height = 3)
pheatmap(Roe,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(colors = c("white","#A52A2A"))(100),
         show_colnames = T,
         show_rownames = T,
         scale="none",
         fontsize = 8,
         fontsize_row = 12,
         fontsize_col = 8,
         main = "Ro/e",
         angle_col = "45",
         display_numbers = T,
         legend = F)
dev.off()

#sample - celltype - barplot

Sample_type <- as.data.frame(table(PCSC$orig.ident,PCSC@active.ident))
colnames(Sample_type) <- c("Sample", "CellType", "nCells")
Sample_type$CellType <- factor(Sample_type$CellType, 
                               levels = c("SMC", "Fibroblast", "Endothelial",
                                          "Mast", "Myeloid", "B-cell", "T-cell", "Epithelial"))
pal_for_celltype <- color[c(1:13)]
names(pal_for_celltype) <- c("Epithelial","Normal","Other_Epi","T-cell","B-cell",
                             "Myeloid","Mast","Endothelial","Fibroblast", "Non1", "Non2", "SMC")


pdf("./2.1.harmony/Sampletype/proportion_sample_type_nolegend.pdf", width = 3, height = 4)
Sample_type %>%
  ggplot(aes(x=Sample, y=nCells,fill=CellType)) +
  geom_bar(stat="identity", position = 'fill') +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = pal_for_celltype) +
  NoLegend()
dev.off()

#sampletype - celltype - barplot

Sample_type <- as.data.frame(table(PCSC$SampleType,PCSC@active.ident))
colnames(Sample_type) <- c("SampleType", "CellType", "nCells")
Sample_type$CellType <- factor(Sample_type$CellType, 
                               levels = c("SMC", "Fibroblast", "Endothelial",
                                          "Mast", "Myeloid", "B-cell", "T-cell", "Epithelial"))
pal_for_celltype <- color[c(1:13)]
names(pal_for_celltype) <- c("Epithelial","Normal","Other_Epi","T-cell","B-cell",
                             "Myeloid","Mast","Endothelial","Fibroblast", "Non1", "Non2", "SMC")


pdf("./2.1.harmony/Sampletype/proportion_sampletype_type_nolegend.pdf", width = 3, height = 1.5)
Sample_type %>%
  ggplot(aes(x=SampleType, y=nCells,fill=CellType)) +
  geom_bar(stat="identity", position = 'fill') +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = pal_for_celltype) +
  NoLegend()
dev.off()

#sampletype - celltype - PiePlot
SampleType_type <- as.data.frame(table(PCSC$SampleType,PCSC@active.ident))
colnames(SampleType_type) <- c("SampleType", "CellType", "nCells")
Sample_type$CellType <- factor(Sample_type$CellType, 
                               levels = c("SMC", "Fibroblast", "Endothelial",
                                          "Mast", "Myeloid", "B-cell", "T-cell", "Epithelial"))
pal_for_celltype <- color[c(1:13)]
names(pal_for_celltype) <- c("Epithelial","Normal","Other_Epi","T-cell","B-cell",
                             "Myeloid","Mast","Endothelial","Fibroblast", "Non1", "Non2", "SMC")

pdf("./2.1.harmony/Sampletype/proportion_sampleType_type.pdf", width = 7, height = 5)
SampleType_type %>%
  ggplot(aes(x=SampleType, y=nCells,fill=CellType)) +
  geom_bar(stat="identity", position = 'fill') +
  coord_polar(theta = "y") +
  theme_classic() + 
  theme(axis.text.y = element_text(angle=90, hjust=0.5, vjust=1.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = pal_for_celltype)
dev.off()

#sample - cellgroup - barplot

PCSC_rename <- RenameIdents(
  object = PCSC,
  "SMC" = "Stromal",
  "Fibroblast" = "Stromal",
  "Endothelial" = "Stromal",
  "Mast" = "Immune",
  "Myeloid" = "Immune",
  "B-cell" = "Immune",
  "T-cell" = "Immune",
  "Epithelial" = "Epithelial"
)

Sample_cellgroup <- as.data.frame(table(PCSC_rename@meta.data$orig.ident,PCSC_rename@active.ident))
colnames(Sample_cellgroup) <- c("Sample", "CellGroup", "nCells")
pal_for_cellGroup <- color[1:3]
names(pal_for_cellGroup) <- c("Epithelial","Immune","Stromal")

pdf("./2.1.harmony/Sampletype/proportion_sample_cellGroup.pdf", width = 3, height = 4)
Sample_cellgroup %>%
  ggplot(aes(x=Sample, y=nCells,fill=CellGroup)) +
  geom_bar(stat="identity", position = 'fill') +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = pal_for_cellGroup) +
  NoLegend()
dev.off()

#sampletype - cellgroup - pieplot
PCSC_rename <- RenameIdents(
  object = PCSC,
  "Epithelial" = "Epithelial",
  "T-cell" = "Immune",
  "B-cell" = "Immune",
  "Myeloid" = "Immune",
  "Mast" = "Immune",
  "Endothelial" = "Stromal",
  "SMC" = "Stromal",
  "Fibroblast" = "Stromal"
)
SampleType_cellgroup <- as.data.frame(table(PCSC_rename$SampleType,PCSC_rename@active.ident))
colnames(SampleType_cellgroup) <- c("SampleType", "CellGroup", "nCells")

pdf("./2.1.harmony/Sampletype/proportion_sampleGroup_type.pdf", width = 7, height = 5)
SampleType_cellgroup %>%
  ggplot(aes(x=SampleType, y=nCells,fill=CellGroup)) +
  geom_bar(stat="identity", position = 'fill') +
  coord_polar(theta = "y") +
  theme_classic() +
  theme(axis.text.y = element_text(angle=90, hjust=0.5, vjust=1.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = pal_for_cellGroup)
dev.off()

write.table(table(PCSC@meta.data$orig.ident,PCSC@active.ident), 
            "./2.1.harmony/Sampletype/origidentSplit_celltype_proportion.txt",sep="\t",quote = F)

write.table(table(PCSC@meta.data$SampleType,PCSC@active.ident), 
            "./2.1.harmony/Sampletype/sampleSplit_celltype_proportion.txt",sep="\t",quote = F)

saveRDS(PCSC, "./2.1.harmony/PCSC_annoted.Rds")

#10.DEG of cell types between EO/LO ------------------------------------------

#把mixed细胞都剔除作差异分析

PCSC <- readRDS("./6.0.cellchat_Epi_all/Epi_all_annotated.Rds")
PCSC <- PCSC[-grep("^RP[SL]", rownames(PCSC)),]
PCSC <- PCSC[-grep("^MT-", rownames(PCSC)),]

PCSC$SampleCelltype <- paste(PCSC$inferCNVCellStat, PCSC$SampleType, sep = "_")
Idents(PCSC) <- "SampleCelltype"

celltype <- c("Cancer","T-cell","B-cell", "Myeloid",
              "Mast","Endothelial","Fibroblast", "SMC")

## find markers ----------------------------------------
marker_condition <- data.frame()
for ( ci in celltype ) {
  tmp.marker <- FindMarkers(
    PCSC, 
    logfc.threshold = 0, 
    min.pct = 0.1,
    ident.1 = paste0(ci, "_EOPC"),
    ident.2 = paste0(ci ,"_LOPC")
  )
  
  tmp.marker$gene <- rownames(tmp.marker)
  tmp.marker$condition <- ifelse(tmp.marker$avg_log2FC > 0, paste0(ci,"_EOPC"), paste0(ci,"_LOPC"))
  tmp.marker$cluster <- ci
  
  #tmp.marker=tmp.marker%>%filter(p_val_adj < 0.01) #p_val_adj值不筛选
  tmp.marker <- as.data.frame(tmp.marker)
  tmp.marker <- tmp.marker %>% arrange(desc(avg_log2FC))
  
  marker_condition <- marker_condition %>% rbind(tmp.marker)
}

saveRDS(marker_condition, "./2.1.harmony/celltype_DEG/DEGs_EO_LO_CellType_Split.Rds")
write.table(marker_condition, file = "./2.1.harmony/celltype_DEG/DEGs_EO_LO_CellType_Split.txt",
            quote = F, sep = "\t", row.names = F, col.names = T)

## sig markers ---------------------------------------------

library(stringr)
marker_condition$sig <- "not_sig"
marker_condition$sig[abs(marker_condition$avg_log2FC) > 0.25 & marker_condition$p_val_adj < 0.01] <- 
  str_sub(marker_condition$condition[abs(marker_condition$avg_log2FC) > 0.25 & marker_condition$p_val_adj < 0.01], 1, -6)

#控制顺序
marker_condition$sig <- factor(marker_condition$sig, levels = c("not_sig", celltype))
marker_condition$cluster <- factor(marker_condition$cluster,levels = celltype)
marker_condition <- marker_condition%>%arrange(cluster,sig)
#控制范围
marker_condition$avg_log2FC[marker_condition$avg_log2FC > 3] = 3
marker_condition$avg_log2FC[marker_condition$avg_log2FC < c(-3)] = -3

#配色
color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")

pal_for_celltype <- color[c(1:13)]
names(pal_for_celltype) <- c("Cancer","Normal","Other_Epi","T-cell","B-cell",
                             "Myeloid","Mast","Endothelial","Fibroblast", "Non1", "Non2", "SMC")


#画图
p1 <- marker_condition %>% 
  ggplot(aes(x = cluster, y = avg_log2FC, color = sig)) +
  geom_jitter(width = 0.25, size = 0.5)+
  scale_color_manual(values = c(pal_for_celltype, "not_sig" = "#dee1e6"))+
  scale_y_continuous("EOPC VS LOPC, average log2FC", expand = c(0.02,0))+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.text.y.left = element_text(size = 14, color = "black"),
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_text(size = 16)
  )
#ggsave("DEG_jitter1.pdf",width = 15,height = 12,units = "cm")  

#画图，改进

library(ggrepel)

marker_condition2 <- marker_condition
marker_condition2$logPadj <- -log10(marker_condition2$p_val_adj)
marker_condition2$logPadj <- ifelse(marker_condition2$avg_log2FC > 0,
                                    marker_condition2$logPadj,
                                    -marker_condition2$logPadj)

plot.list <- list()
for (ci in celltype) {
  tmpdf <- marker_condition2 %>% filter(cluster == ci)
  minabs <- abs(min(tmpdf$logPadj))
  maxabs <- max(tmpdf$logPadj)
  thre = 0
  if(minabs < maxabs) {
    tmpdf$logPadj[tmpdf$logPadj > minabs] = minabs
    thre=minabs
  }
  if(minabs > maxabs) {
    tmpdf$logPadj[tmpdf$logPadj < (-maxabs)] = -maxabs
    thre=maxabs
  }
  if(minabs == maxabs & maxabs == Inf) {
    thre = min(
      abs(
        range(
          tmpdf$logPadj[tmpdf$logPadj < Inf & tmpdf$logPadj > -Inf]
        )
      )
    )
    tmpdf$logPadj[tmpdf$logPadj < (-thre)] = -thre
    tmpdf$logPadj[tmpdf$logPadj > thre] = thre
  }
  
  plotdata <- tmpdf
  tmpdf <- tmpdf %>% filter(sig != "not_sig") 
  tmpdf=tmpdf%>%arrange(desc(avg_log2FC))
  tmpdf.a=head(tmpdf%>%filter(avg_log2FC > 0),5)
  tmpdf.a$d=thre*2*0.05+(-thre)-tmpdf.a$logPadj
  tmpdf.b=tail(tmpdf%>%filter(avg_log2FC < 0),5)
  tmpdf.b$d=thre*2*0.95-thre  - tmpdf.b$logPadj
  textdata.down = tmpdf.b
  textdata.up   = tmpdf.a
  
  ###画图
  tmpplot <- plotdata %>%
    ggplot(aes(x = logPadj, y = avg_log2FC))+
    geom_point(aes(color = sig), size=1)+
    geom_hline(yintercept = c(-0.25,0.25),linetype="dashed")+
    geom_text_repel(data = textdata.down,
                    mapping = aes(label=gene),
                    nudge_x = textdata.down$d, size = 3,
                    direction = "y", hjust = 1, segment.size = 0.2)+
    geom_text_repel(data = textdata.up,
                    mapping = aes(label=gene),
                    nudge_x = textdata.up$d, size = 3,
                    direction = "y", hjust = 0, segment.size = 0.2)+
    labs(title = ci) +
    scale_color_manual(values = c(pal_for_celltype, "not_sig"="#dee1e6")) +
    scale_y_continuous("EOPC VS LOPC, average log2FC", expand = c(0.02,0), limits = c(-3,3)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      
      axis.ticks.x.bottom = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.title.x.bottom = element_blank(),
      axis.text.y.left = element_text(size = 10,color = "black"),
      axis.title.y.left = element_text(size = 12),
      
      plot.title = element_text(size = 12,hjust = 0.5)
    )
  
  index <- which(ci == celltype)
  if (index != 1) {
    tmpplot = tmpplot+theme(
      axis.title.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.text.y.left = element_blank()
    )
  }
  if (index == length(sort(unique(as.character(marker_condition2$cluster))))) {
    segment.df <- data.frame(x = c(0 - thre / 5,0 + thre / 5),
                             xend = c(-thre,thre),
                             y = c(-3,-3),
                             yend = c(-3,-3))
    tmpplot <- tmpplot+
      geom_segment(data = segment.df,
                   mapping = aes(x = x, xend = xend, y = y , yend = yend),
                   arrow = arrow(length=unit(0.3, "cm")))
    
  }
  plot.list[[get("index")]]=tmpplot
}

library(patchwork)
p2 <- wrap_plots(plot.list, ncol = 8) & theme(plot.margin = unit(c(0,0,0,0), "cm"))
ggsave("./2.1.harmony/celltype_DEG/DEG_jitter.pdf", plot = p2, width = 7.5, height = 4)

## find markers 0.3 ----------------------------------------
marker_condition <- data.frame()
for ( ci in celltype ) {
  tmp.marker <- FindMarkers(
    PCSC, 
    logfc.threshold = 0, 
    min.pct = 0.3,
    ident.1 = paste0(ci, "_EOPC"),
    ident.2 = paste0(ci ,"_LOPC")
  )
  
  tmp.marker$gene <- rownames(tmp.marker)
  tmp.marker$condition <- ifelse(tmp.marker$avg_log2FC > 0, paste0(ci,"_EOPC"), paste0(ci,"_LOPC"))
  tmp.marker$cluster <- ci
  
  #tmp.marker=tmp.marker%>%filter(p_val_adj < 0.01) #p_val_adj值不筛选
  tmp.marker <- as.data.frame(tmp.marker)
  tmp.marker <- tmp.marker %>% arrange(desc(avg_log2FC))
  
  marker_condition <- marker_condition %>% rbind(tmp.marker)
}

saveRDS(marker_condition, "./2.1.harmony/celltype_DEG/DEGs_EO_LO_CellType_Split.Rds")
write.table(marker_condition, file = "./2.1.harmony/celltype_DEG/DEGs_EO_LO_CellType_Split_proportion0.3.txt",
            quote = F, sep = "\t", row.names = F, col.names = T)

## sig markers ---------------------------------------------

library(stringr)
marker_condition$sig <- "not_sig"
marker_condition$sig[abs(marker_condition$avg_log2FC) > 0.25 & marker_condition$p_val_adj < 0.01] <- 
  str_sub(marker_condition$condition[abs(marker_condition$avg_log2FC) > 0.25 & marker_condition$p_val_adj < 0.01], 1, -6)

#控制顺序
marker_condition$sig <- factor(marker_condition$sig, levels = c("not_sig", celltype))
marker_condition$cluster <- factor(marker_condition$cluster,levels = celltype)
marker_condition <- marker_condition%>%arrange(cluster,sig)
#控制范围
marker_condition$avg_log2FC[marker_condition$avg_log2FC > 3] = 3
marker_condition$avg_log2FC[marker_condition$avg_log2FC < c(-3)] = -3

#配色
color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")

pal_for_celltype <- color[c(1:13)]
names(pal_for_celltype) <- c("Cancer","Normal","Other_Epi","T-cell","B-cell",
                             "Myeloid","Mast","Endothelial","Fibroblast", "Non1", "Non2", "SMC")


#画图
p1 <- marker_condition %>% 
  ggplot(aes(x = cluster, y = avg_log2FC, color = sig)) +
  geom_jitter(width = 0.25, size = 0.5)+
  scale_color_manual(values = c(pal_for_celltype, "not_sig" = "#dee1e6"))+
  scale_y_continuous("EOPC VS LOPC, average log2FC", expand = c(0.02,0))+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.text.y.left = element_text(size = 14, color = "black"),
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_text(size = 16)
  )
#ggsave("DEG_jitter1.pdf",width = 15,height = 12,units = "cm")  

#画图，改进

library(ggrepel)

marker_condition2 <- marker_condition
marker_condition2$logPadj <- -log10(marker_condition2$p_val_adj)
marker_condition2$logPadj <- ifelse(marker_condition2$avg_log2FC > 0,
                                    marker_condition2$logPadj,
                                    -marker_condition2$logPadj)

plot.list <- list()
for (ci in celltype) {
  tmpdf <- marker_condition2 %>% filter(cluster == ci)
  minabs <- abs(min(tmpdf$logPadj))
  maxabs <- max(tmpdf$logPadj)
  thre = 0
  if(minabs < maxabs) {
    tmpdf$logPadj[tmpdf$logPadj > minabs] = minabs
    thre=minabs
  }
  if(minabs > maxabs) {
    tmpdf$logPadj[tmpdf$logPadj < (-maxabs)] = -maxabs
    thre=maxabs
  }
  if(minabs == maxabs & maxabs == Inf) {
    thre = min(
      abs(
        range(
          tmpdf$logPadj[tmpdf$logPadj < Inf & tmpdf$logPadj > -Inf]
        )
      )
    )
    tmpdf$logPadj[tmpdf$logPadj < (-thre)] = -thre
    tmpdf$logPadj[tmpdf$logPadj > thre] = thre
  }
  
  plotdata <- tmpdf
  tmpdf <- tmpdf %>% filter(sig != "not_sig") 
  tmpdf=tmpdf%>%arrange(desc(avg_log2FC))
  tmpdf.a=head(tmpdf%>%filter(avg_log2FC > 0),5)
  tmpdf.a$d=thre*2*0.05+(-thre)-tmpdf.a$logPadj
  tmpdf.b=tail(tmpdf%>%filter(avg_log2FC < 0),5)
  tmpdf.b$d=thre*2*0.95-thre  - tmpdf.b$logPadj
  textdata.down = tmpdf.b
  textdata.up   = tmpdf.a
  
  ###画图
  tmpplot <- plotdata %>%
    ggplot(aes(x = logPadj, y = avg_log2FC))+
    geom_point(aes(color = sig), size=1)+
    geom_hline(yintercept = c(-0.25,0.25),linetype="dashed")+
    geom_text_repel(data = textdata.down,
                    mapping = aes(label=gene),
                    nudge_x = textdata.down$d, size = 3,
                    direction = "y", hjust = 1, segment.size = 0.2)+
    geom_text_repel(data = textdata.up,
                    mapping = aes(label=gene),
                    nudge_x = textdata.up$d, size = 3,
                    direction = "y", hjust = 0, segment.size = 0.2)+
    labs(title = ci) +
    scale_color_manual(values = c(pal_for_celltype, "not_sig"="#dee1e6")) +
    scale_y_continuous("EOPC VS LOPC, average log2FC", expand = c(0.02,0), limits = c(-3,3)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      
      axis.ticks.x.bottom = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.title.x.bottom = element_blank(),
      axis.text.y.left = element_text(size = 10,color = "black"),
      axis.title.y.left = element_text(size = 12),
      
      plot.title = element_text(size = 12,hjust = 0.5)
    )
  
  index <- which(ci == celltype)
  if (index != 1) {
    tmpplot = tmpplot+theme(
      axis.title.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.text.y.left = element_blank()
    )
  }
  if (index == length(sort(unique(as.character(marker_condition2$cluster))))) {
    segment.df <- data.frame(x = c(0 - thre / 5,0 + thre / 5),
                             xend = c(-thre,thre),
                             y = c(-3,-3),
                             yend = c(-3,-3))
    tmpplot <- tmpplot+
      geom_segment(data = segment.df,
                   mapping = aes(x = x, xend = xend, y = y , yend = yend),
                   arrow = arrow(length=unit(0.3, "cm")))
    
  }
  plot.list[[get("index")]]=tmpplot
}

library(patchwork)
p2 <- wrap_plots(plot.list, ncol = 8) & theme(plot.margin = unit(c(0,0,0,0), "cm"))
ggsave("./2.1.harmony/celltype_DEG/DEG_jitter_proportion0.3.pdf", plot = p2, width = 7.5, height = 4)


#11.GSEA of cell types between EO/LO ------------------------------------------

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)

marker_condition <- readRDS("2.1.harmony/celltype_DEG/DEGs_EO_LO_CellType_Split.Rds")
celltype <- c("Cancer","T-cell","B-cell", "Myeloid",
              "Mast","Endothelial","Fibroblast", "SMC")

## cancersea -----------------------------------------------

hsets <- read.gmt("2.1.harmony/cancersea.gmt")

# GSEA

Plist <- list()

for (i in celltype) {
  
  df <- marker_condition[marker_condition$cluster == i,]
  #转换基因ID  如果下载的是Gene Symbols，可以不转
  #df_id<-bitr(df$SYMBOL, #转换的列是df数据框中的SYMBOL列
  #            fromType = "SYMBOL",#需要转换ID类型
  #            toType = "ENTREZID",#转换成的ID类型
  #            OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db
  #df_all<-merge(df,df_id,by="SYMBOL",all=F)
  #### GSEA分析
  
  df <- df[order(-df$avg_log2FC),]  #先按照logFC降序排序
  gene_fc <- df$avg_log2FC  #把foldchange按照从大到小提取出来
  names(gene_fc) <- df$gene
  
  egmt <- GSEA(gene_fc, TERM2GENE = hsets, pAdjustMethod = "fdr", pvalueCutoff = 1)  #GSEA分析
  
  egmt_result_df <- as.data.frame(egmt)  #转换成数据框
  egmt_result_df$Cluster <- i
  saveRDS(egmt_result_df, paste0("2.1.harmony/celltype_enrich/cancersea/GSEA_results/", i, "_GSEA.Rds"))
  
  #### 结果展示
  #气泡图 展示geneset被激活还是抑制
  p <- dotplot(egmt, split=".sign", showCategory = 8, font.size = 16,
               label_format = 25, color = "pvalue", title = i) + 
    facet_grid(~.sign)
  #edit legends
  #+guides(
  #reverse color order (higher value on top)
  #color = guide_colorbar(reverse = TRUE))
  #reverse size order (higher diameter on top) 
  #size = guide_legend(reverse = TRUE))
  pdf(paste0("2.1.harmony/celltype_enrich/cancersea/", i, "_cancersea_pathway_dotplot.pdf"), height=8, width=8)
  print(p)
  dev.off()
  
  Plist[[i]] <- p
}

library(patchwork)
p2 <- wrap_plots(Plist, ncol = 2) & theme(plot.margin = unit(c(0,0,0,0), "cm"))
pdf("2.1.harmony/celltype_enrich/cancersea/cancersea_pathway_dotplot.pdf", height=28, width=16)
print(p2)
dev.off()

### UpSet plot (activated) -------------------------------------------

library(UpSetR)
outFile="intersect_activate.txt"        #输出交集基因文件
outPic="upset_activate.pdf"                  #输出图片
setwd("2.1.harmony/celltype_enrich/cancersea/GSEA_results/")      #设置工作目录

files=dir()                          #获取目录下所有文件
files=grep("Rds$",files,value=T)     #提取.txt结尾的文件
gsList=list()

#获取所有txt文件中的基因信息，保存到geneList
for(i in 1:length(files)){
  inputFile=files[i]
  rt=readRDS(inputFile)               #读取输入文件
  rt=rt[which(rt$NES > 0 & rt$pvalue < 0.05),]
  gsNames=as.vector(rt[,1])               #提取基因名称
  header=unlist(strsplit(inputFile,"_"))
  gsList[[header[1]]]=gsNames
  gsLength=length(gsNames)
  print(paste(header[1],gsLength,sep=" "))
}

#绘制UpSetͼ
upsetData=fromList(gsList)
pdf(file=outPic,onefile = FALSE,width=6,height=5)
upset(upsetData,
      nsets = length(gsList),               #展示多少个数据
      nintersects = NA,                       #展示基因集数目
      order.by = "freq",                      #按照数目排序
      show.numbers = "no",                   #柱状图上方是否显示数值
      number.angles = 20,                     #字体角度
      point.size = 2,                         #点的大小
      matrix.color="red",                     #交集点颜色
      line.size = 0.8,                        #线条粗细
      mainbar.y.label = "Gene Intersections", 
      sets.x.label = "Set Size")
dev.off()

### UpSet plot (suppressed) -------------------------------------------

library(UpSetR)
outFile="intersect_suppressed.txt"        #输出交集基因文件
outPic="upset_suppressed.pdf"                  #输出图片
setwd("2.1.harmony/celltype_enrich/cancersea/GSEA_results/")      #设置工作目录

files=dir()                          #获取目录下所有文件
files=grep("Rds$",files,value=T)     #提取.txt结尾的文件
gsList=list()

#获取所有txt文件中的基因信息，保存到geneList
for(i in 1:length(files)){
  inputFile=files[i]
  rt=readRDS(inputFile)               #读取输入文件
  rt=rt[which(rt$NES < 0 & rt$pvalue < 0.05),]
  gsNames=as.vector(rt[,1])               #提取基因名称
  header=unlist(strsplit(inputFile,"_"))
  gsList[[header[1]]]=gsNames
  gsLength=length(gsNames)
  print(paste(header[1],gsLength,sep=" "))
}

#绘制UpSetͼ
upsetData=fromList(gsList)
pdf(file=outPic,onefile = FALSE,width=6,height=5)
upset(upsetData,
      nsets = length(gsList),               #展示多少个数据
      nintersects = NA,                       #展示基因集数目
      order.by = "freq",                      #按照数目排序
      show.numbers = "no",                   #柱状图上方是否显示数值
      number.angles = 20,                     #字体角度
      point.size = 2,                         #点的大小
      matrix.color="red",                     #交集点颜色
      line.size = 0.8,                        #线条粗细
      mainbar.y.label = "Gene Intersections", 
      sets.x.label = "Set Size")
dev.off()

## KEGG -----------------------------------------------

library(msigdbr)
KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
KEGG$gs_name <- gsub("KEGG_", "", KEGG$gs_name)
KEGG$gs_name <- gsub("_", " ", KEGG$gs_name)
KEGG$gs_name <- str_to_title(KEGG$gs_name)

# GSEA

Plist <- list()

for (i in celltype) {
  
  df <- marker_condition[marker_condition$cluster == i,]
  
  #### GSEA分析
  
  df <- df[order(-df$avg_log2FC),]  #先按照logFC降序排序
  gene_fc <- df$avg_log2FC  #把foldchange按照从大到小提取出来
  names(gene_fc) <- df$gene
  
  egmt <- GSEA(gene_fc, TERM2GENE = KEGG, pAdjustMethod = "fdr", pvalueCutoff = 1)  #GSEA分析
  
  egmt_result_df <- as.data.frame(egmt)  #转换成数据框
  egmt_result_df$Cluster <- i
  saveRDS(egmt_result_df, paste0("2.1.harmony/celltype_enrich/cancersea/GSEA_results/", i, "_GSEA.Rds"))
  
  #### 结果展示
  #气泡图 展示geneset被激活还是抑制
  p <- dotplot(egmt, split=".sign", showCategory = 8, font.size = 16,
               label_format = 25, color = "pvalue", title = i) + 
    facet_grid(~.sign)
  #edit legends
  #+guides(
  #reverse color order (higher value on top)
  #color = guide_colorbar(reverse = TRUE))
  #reverse size order (higher diameter on top) 
  #size = guide_legend(reverse = TRUE))
  pdf(paste0("2.1.harmony/celltype_enrich/cancersea/", i, "_cancersea_pathway_dotplot.pdf"), height=8, width=8)
  print(p)
  dev.off()
  
  Plist[[i]] <- p
}

library(patchwork)
p2 <- wrap_plots(Plist, ncol = 2) & theme(plot.margin = unit(c(0,0,0,0), "cm"))
pdf("2.1.harmony/celltype_enrich/cancersea/cancersea_pathway_dotplot.pdf", height=28, width=16)
print(p2)
dev.off()

### UpSet plot (activated) -------------------------------------------

library(UpSetR)
outFile="intersect_activate.txt"        #输出交集基因文件
outPic="upset_activate.pdf"                  #输出图片
setwd("2.1.harmony/celltype_enrich/cancersea/GSEA_results/")      #设置工作目录

files=dir()                          #获取目录下所有文件
files=grep("Rds$",files,value=T)     #提取.txt结尾的文件
gsList=list()

#获取所有txt文件中的基因信息，保存到geneList
for(i in 1:length(files)){
  inputFile=files[i]
  rt=readRDS(inputFile)               #读取输入文件
  rt=rt[which(rt$NES > 0 & rt$pvalue < 0.05),]
  gsNames=as.vector(rt[,1])               #提取基因名称
  header=unlist(strsplit(inputFile,"_"))
  gsList[[header[1]]]=gsNames
  gsLength=length(gsNames)
  print(paste(header[1],gsLength,sep=" "))
}

#绘制UpSetͼ
upsetData=fromList(gsList)
pdf(file=outPic,onefile = FALSE,width=6,height=5)
upset(upsetData,
      nsets = length(gsList),               #展示多少个数据
      nintersects = NA,                       #展示基因集数目
      order.by = "freq",                      #按照数目排序
      show.numbers = "no",                   #柱状图上方是否显示数值
      number.angles = 20,                     #字体角度
      point.size = 2,                         #点的大小
      matrix.color="red",                     #交集点颜色
      line.size = 0.8,                        #线条粗细
      mainbar.y.label = "Gene Intersections", 
      sets.x.label = "Set Size")
dev.off()

### UpSet plot (activated) -------------------------------------------

library(UpSetR)
outFile="intersect_suppressed.txt"        #输出交集基因文件
outPic="upset_suppressed.pdf"                  #输出图片
setwd("2.1.harmony/celltype_enrich/cancersea/GSEA_results/")      #设置工作目录

files=dir()                          #获取目录下所有文件
files=grep("Rds$",files,value=T)     #提取.txt结尾的文件
gsList=list()

#获取所有txt文件中的基因信息，保存到geneList
for(i in 1:length(files)){
  inputFile=files[i]
  rt=readRDS(inputFile)               #读取输入文件
  rt=rt[which(rt$NES < 0 & rt$pvalue < 0.05),]
  gsNames=as.vector(rt[,1])               #提取基因名称
  header=unlist(strsplit(inputFile,"_"))
  gsList[[header[1]]]=gsNames
  gsLength=length(gsNames)
  print(paste(header[1],gsLength,sep=" "))
}

#绘制UpSetͼ
upsetData=fromList(gsList)
pdf(file=outPic,onefile = FALSE,width=6,height=5)
upset(upsetData,
      nsets = length(gsList),               #展示多少个数据
      nintersects = NA,                       #展示基因集数目
      order.by = "freq",                      #按照数目排序
      show.numbers = "no",                   #柱状图上方是否显示数值
      number.angles = 20,                     #字体角度
      point.size = 2,                         #点的大小
      matrix.color="red",                     #交集点颜色
      line.size = 0.8,                        #线条粗细
      mainbar.y.label = "Gene Intersections", 
      sets.x.label = "Set Size")
dev.off()

#单个通路
pdf(paste0(substr(dir[i],1,6),"single_pathway.pdf"),height=6,width=10)
gseaplot2(egmt, geneSetID = c(3), subplots = 1:3,pvalue_table = TRUE)
dev.off()

#12.scMetabolism of cells ------------------------------------------

library(scMetabolism)
library(ggplot2)
library(rsvd)

#使用去除mixed的数据

PCSC <- readRDS("./6.0.cellchat_Epi_all/Epi_all_annotated.Rds")
PCSC <- PCSC[-grep("^RP[SL]", rownames(PCSC)),]
PCSC <- PCSC[-grep("^MT-", rownames(PCSC)),]

PCSC_Metab <- sc.metabolism.Seurat(obj = PCSC, 
                                   method = "AUCell", 
                                   imputation = F, 
                                   ncores = 2, 
                                   metabolism.type = "KEGG")
saveRDS(PCSC_Metab, "2.1.harmony/scMetabolism/PCSC_Metab_KEGG.Rds")

PCSC_Metab <- sc.metabolism.Seurat(obj = PCSC, 
                                   method = "AUCell", 
                                   imputation = F, 
                                   ncores = 2, 
                                   metabolism.type = "REACTOME")
saveRDS(PCSC_Metab, "2.1.harmony/scMetabolism/PCSC_Metab_REACTOME.Rds")

## REACTOME met celltype EO/LO ---------------------------------------

PCSC <- readRDS("./6.0.cellchat_Epi_all/Epi_all_annotated.Rds")
PCSC_Metab <- readRDS("2.1.harmony/scMetabolism/PCSC_Metab_REACTOME.Rds")
Met_mtx <- PCSC_Metab@assays$METABOLISM$score

celltype <- c("Cancer","T-cell","B-cell", "Myeloid",
              "Mast","Endothelial","Fibroblast", "SMC")

Plist <- list()
for (i in 1:length(celltype)) {
  
  Seu_tmp <- subset(PCSC, inferCNVCellStat == celltype[i])
  type <- Seu_tmp$SampleType
  type <- sort(type)
  
  mtx <- Met_mtx[, names(type)]
  
  ## limma
  library(limma)
  
  # 设置或导入分组
  group <- factor(c(rep("EOPC", table(type)[1]), rep("LOPC", table(type)[2])), levels = c('EOPC', 'LOPC'))
  design <- model.matrix(~0+group)
  colnames(design) = levels(factor(group))
  rownames(design) = colnames(mtx)
  head(design)
  
  # Tunor VS Normal
  compare <- makeContrasts(EOPC - LOPC, levels=design)
  fit <- lmFit(mtx, design)
  fit2 <- contrasts.fit(fit, compare)
  fit3 <- eBayes(fit2)
  Diff <- topTable(fit3, coef=1, number=200)
  head(Diff)
  
  # 排序
  Diff <- Diff %>% arrange(t)
  
  # barplot
  
  dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
  dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
  if (nrow(dat_plot[dat_plot$t < 0,]) < 15 & nrow(dat_plot[dat_plot$t > 0,]) < 15) {
    down <- nrow(dat_plot[dat_plot$t < 0,])
    up <- nrow(dat_plot[dat_plot$t > 0,])
  } else if (nrow(dat_plot[dat_plot$t < 0,]) < 15 & nrow(dat_plot[dat_plot$t > 0,]) >= 15) {
    down <- nrow(dat_plot[dat_plot$t < 0,])
    up <- 15
  } else if (nrow(dat_plot[dat_plot$t < 0,]) >= 15 & nrow(dat_plot[dat_plot$t > 0,]) < 15) {
    down <- 15
    up <- nrow(dat_plot[dat_plot$t > 0,])
  } else {
    up = 15
    down = 15
  }
  dat_plot <- rbind(top_n(dat_plot, down, wt = t), top_n(dat_plot, -up, wt = t))
  dat_plot$threshold = factor(ifelse(dat_plot$t  > 0, 'Up', 'Down'), levels=c('Up','Down','NoSignifi'))
  
  # 变成因子类型
  dat_plot <- dat_plot %>% arrange(t)
  dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
  
  # 绘制
  library(ggplot2)
  library(ggprism)
  
  p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
    geom_col()+
    coord_flip() +
    scale_fill_manual(values = c('Up'= "#A52A2ACC", 'Down'='#053E7ACC')) +
    xlab('') + 
    ylab('t value of AUCell score, MP1_On versus Off') + #注意坐标轴旋转了
    ggtitle(celltype[i]) +
    guides(fill=F)+ # 不显示图例
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  # 添加标签
  
  lim <- max(max(dat_plot$t), abs(min(dat_plot$t)))
  
  p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                         aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black', size = 2.5) +
    geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
              aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black', size = 2.5) +
    ylim(c(-lim,lim))
  
  Plist[[celltype[i]]] <- p_MP1
}

library(gridExtra)
pdf("2.1.harmony/scMetabolism/REACTOME_Met_diff_top15.pdf", width = 16, height = 20)
grid.arrange(grobs = Plist, ncol = 2)
dev.off()

## KEGG met celltype EO/LO ---------------------------------------

PCSC <- readRDS("./6.0.cellchat_Epi_all/Epi_all_annotated.Rds")
PCSC_Metab <- readRDS("2.1.harmony/scMetabolism/PCSC_Metab_KEGG.Rds")
Met_mtx <- PCSC_Metab@assays$METABOLISM$score

celltype <- c("Cancer","T-cell","B-cell", "Myeloid",
              "Mast","Endothelial","Fibroblast", "SMC")

Plist <- list()
for (i in 1:length(celltype)) {
  
  Seu_tmp <- subset(PCSC, inferCNVCellStat == celltype[i])
  type <- Seu_tmp$SampleType
  type <- sort(type)
  
  mtx <- Met_mtx[, names(type)]
  
  ## limma
  library(limma)
  
  # 设置或导入分组
  group <- factor(c(rep("EOPC", table(type)[1]), rep("LOPC", table(type)[2])), levels = c('EOPC', 'LOPC'))
  design <- model.matrix(~0+group)
  colnames(design) = levels(factor(group))
  rownames(design) = colnames(mtx)
  head(design)
  
  # Tunor VS Normal
  compare <- makeContrasts(EOPC - LOPC, levels=design)
  fit <- lmFit(mtx, design)
  fit2 <- contrasts.fit(fit, compare)
  fit3 <- eBayes(fit2)
  Diff <- topTable(fit3, coef=1, number=200)
  head(Diff)
  
  # 排序
  Diff <- Diff %>% arrange(t)
  
  # barplot
  
  dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
  dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
  if (nrow(dat_plot[dat_plot$t < 0,]) < 15 & nrow(dat_plot[dat_plot$t > 0,]) < 15) {
    down <- nrow(dat_plot[dat_plot$t < 0,])
    up <- nrow(dat_plot[dat_plot$t > 0,])
  } else if (nrow(dat_plot[dat_plot$t < 0,]) < 15 & nrow(dat_plot[dat_plot$t > 0,]) >= 15) {
    down <- nrow(dat_plot[dat_plot$t < 0,])
    up <- 15
  } else if (nrow(dat_plot[dat_plot$t < 0,]) >= 15 & nrow(dat_plot[dat_plot$t > 0,]) < 15) {
    down <- 15
    up <- nrow(dat_plot[dat_plot$t > 0,])
  } else {
    up = 15
    down = 15
  }
  dat_plot <- rbind(top_n(dat_plot, down, wt = t), top_n(dat_plot, -up, wt = t))
  dat_plot$threshold = factor(ifelse(dat_plot$t  > 0, 'Up', 'Down'), levels=c('Up','Down','NoSignifi'))
  
  # 变成因子类型
  dat_plot <- dat_plot %>% arrange(t)
  dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
  
  # 绘制
  library(ggplot2)
  library(ggprism)
  
  p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
    geom_col()+
    coord_flip() +
    scale_fill_manual(values = c('Up'= "#A52A2ACC", 'Down'='#053E7ACC')) +
    xlab('') + 
    ylab('t value of AUCell score, MP1_On versus Off') + #注意坐标轴旋转了
    ggtitle(celltype[i]) +
    guides(fill=F)+ # 不显示图例
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  # 添加标签
  
  lim <- max(max(dat_plot$t), abs(min(dat_plot$t)))
  
  p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                         aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black', size = 2.5) +
    geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
              aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black', size = 2.5) +
    ylim(c(-lim,lim))
  
  Plist[[celltype[i]]] <- p_MP1
}

library(gridExtra)
pdf("2.1.harmony/scMetabolism/KEGG_Met_diff_top15.pdf", width = 16, height = 20)
grid.arrange(grobs = Plist, ncol = 2)
dev.off()


#13.scFEA of cells ---------------------------------------------

#使用去除mixed的数据

PCSC <- readRDS("./6.0.cellchat_Epi_all/Epi_all_annotated.Rds")
PCSC <- PCSC[-grep("^RP[SL]", rownames(PCSC)),]
PCSC <- PCSC[-grep("^MT-", rownames(PCSC)),]

exp <- as.matrix(PCSC@assays$RNA@counts)
write.csv(exp, "2.1.harmony/scFEA/input/expr.csv", quote = F, sep = "\t")

#bash ../../SC/2.1.harmony/scFEA/scFEA.bash

scFEA <- read.csv("2.1.harmony/scFEA/output/balance_20230326-225854.csv", row.names =1)
scFEA <- t(scFEA)
PCSC[["FEA"]] <- CreateAssayObject(counts = scFEA)

celltype <- c("Cancer","T-cell","B-cell", "Myeloid",
              "Mast","Endothelial","Fibroblast", "SMC")

Plist <- list()
for (i in 1:length(celltype)) {
  
  Seu_tmp <- subset(PCSC, inferCNVCellStat == celltype[i])
  type <- Seu_tmp$SampleType
  type <- sort(type)
  
  mtx <- scFEA[, names(type)]
  
  ## limma
  library(limma)
  
  # 设置或导入分组
  group <- factor(c(rep("EOPC", table(type)[1]), rep("LOPC", table(type)[2])), levels = c('EOPC', 'LOPC'))
  design <- model.matrix(~0+group)
  colnames(design) = levels(factor(group))
  rownames(design) = colnames(mtx)
  head(design)
  
  # Tunor VS Normal
  compare <- makeContrasts(EOPC - LOPC, levels=design)
  fit <- lmFit(mtx, design)
  fit2 <- contrasts.fit(fit, compare)
  fit3 <- eBayes(fit2)
  Diff <- topTable(fit3, coef=1, number=200)
  head(Diff)
  
  # 排序
  Diff <- Diff %>% arrange(t)
  
  # barplot
  
  dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
  dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
  if (nrow(dat_plot[dat_plot$t < 0,]) < 15 & nrow(dat_plot[dat_plot$t > 0,]) < 15) {
    down <- nrow(dat_plot[dat_plot$t < 0,])
    up <- nrow(dat_plot[dat_plot$t > 0,])
  } else if (nrow(dat_plot[dat_plot$t < 0,]) < 15 & nrow(dat_plot[dat_plot$t > 0,]) >= 15) {
    down <- nrow(dat_plot[dat_plot$t < 0,])
    up <- 15
  } else if (nrow(dat_plot[dat_plot$t < 0,]) >= 15 & nrow(dat_plot[dat_plot$t > 0,]) < 15) {
    down <- 15
    up <- nrow(dat_plot[dat_plot$t > 0,])
  } else {
    up = 15
    down = 15
  }
  dat_plot <- rbind(top_n(dat_plot, down, wt = t), top_n(dat_plot, -up, wt = t))
  dat_plot$threshold = factor(ifelse(dat_plot$t  > 0, 'Up', 'Down'), levels=c('Up','Down','NoSignifi'))
  
  # 变成因子类型
  dat_plot <- dat_plot %>% arrange(t)
  dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
  
  # 绘制
  library(ggplot2)
  library(ggprism)
  
  p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
    geom_col()+
    coord_flip() +
    scale_fill_manual(values = c('Up'= "#A52A2ACC", 'Down'='#053E7ACC')) +
    xlab('') + 
    ylab('t value of scFEA score, EOPC versus LOPC') + #注意坐标轴旋转了
    ggtitle(celltype[i]) +
    guides(fill=F)+ # 不显示图例
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  # 添加标签
  
  lim <- max(max(dat_plot$t), abs(min(dat_plot$t)))
  
  p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                         aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black', size = 2.5) +
    geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
              aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black', size = 2.5) +
    ylim(c(-lim,lim))
  
  Plist[[celltype[i]]] <- p_MP1
}

library(gridExtra)
pdf("2.1.harmony/scFEA/scFEA_diff_top15.pdf", width = 8, height = 16)
grid.arrange(grobs = Plist, ncol = 2)
dev.off()

pdf("2.1.harmony/scFEA/scFEA_diff_top15_noEpi_noBcell.pdf", width = 8, height = 12)
grid.arrange(grobs = Plist[-c(1,3)], ncol = 2)
dev.off()

# 14.Compass analyses -----------------------------------------

PCSC <- readRDS("./6.0.cellchat_Epi_all/Epi_all_annotated.Rds")
PCSC <- PCSC[-grep("^RP[SL]", rownames(PCSC)),]
PCSC <- PCSC[-grep("^MT-", rownames(PCSC)),]

table(PCSC$infe)

mtx <- c()
for (ct in unique(PCSC$inferCNVCellStat)) {
  exp <- subset(PCSC, inferCNVCellStat == ct)@assays$RNA@data %>% as.matrix()
  exp_m <- rowMeans(exp)
  mtx <- cbind(mtx, exp_m)
  colnames(mtx)[length(colnames(mtx))] <- ct
  
  for (st in unique(PCSC$SampleType)) {
    exp <- subset(PCSC, inferCNVCellStat == ct & SampleType == st)@assays$RNA@data %>% as.matrix()
    exp_m <- rowMeans(exp)
    mtx <- cbind(mtx, exp_m)
    colnames(mtx)[length(colnames(mtx))] <- paste(ct, st, sep = "_")
  }
}

Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
DefaultAssay(Epi) <- "RNA"
# MP + SampleType
Epi$SampleType_MP1 <- paste(Epi$SampleType, Epi$MP1_OnOff, sep = "_")
Epi <- Epi[-grep("^RP[SL]", rownames(Epi)),]
Epi <- Epi[-grep("^MT-", rownames(Epi)),]

for (ct in unique(Epi$MP1_OnOff)) {
  exp <- subset(Epi, MP1_OnOff == ct)@assays$RNA@data %>% as.matrix()
  exp_m <- rowMeans(exp)
  mtx <- cbind(mtx, exp_m)
  colnames(mtx)[length(colnames(mtx))] <- paste("Cancer_MP1", ct, sep = "_")
  
  for (st in unique(Epi$SampleType)) {
    exp <- subset(Epi, MP1_OnOff == ct & SampleType == st)@assays$RNA@data %>% as.matrix()
    exp_m <- rowMeans(exp)
    mtx <- cbind(mtx, exp_m)
    colnames(mtx)[length(colnames(mtx))] <- paste("Cancer_MP1", ct, st, sep = "_")
  }
}

write.csv(mtx, "2.1.harmony/compass/celltype_exp_mean.csv")
write.table(mtx, "2.1.harmony/compass/celltype_exp_mean.tsv", quote=FALSE, sep='\t', col.names = NA)

## random cells -----------------------------------------

select_cell <- c()
for (ct in unique(PCSC$inferCNVCellStat)) {
  for (st in unique(PCSC$SampleType)) {
    id <- subset(PCSC, inferCNVCellStat == ct & SampleType == st) %>% Cells()
    if (length(id) / 50 < 50) {id_s <- sample(id, 50)} else
      if (length(id) / 50 > 200) {id_s <- sample(id, 200)} else {id_s <- sample(id, length(id) / 50)}
    select_cell <- c(select_cell, id_s)
  }
}

PCSC_s <- subset(PCSC, cells = select_cell)
exp <- PCSC_s@assays$RNA@data
write.table(exp, "2.1.harmony/compass/select_cells/selected_Cells_exp_mean.tsv", quote=FALSE, sep='\t', col.names = NA)


#compass --data celltype_exp_mean.tsv --num-processes 10 --num-threads 20 --species homo_sapiens
