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
library(tibble)
library(reshape2)

Epi <- readRDS("5.1.epithelial_harmony_TryTwo/Epi_annoted.Rds")
Epi@meta.data <- Epi@meta.data[,c(1:8,19,20,25)]
Tcell <- readRDS("5.2.Tcell_harmony/Tcell_annotated.Rds")
Tcell@meta.data <- Tcell@meta.data[,c(1:8,19,20,25)]
Epi_T <- merge(Epi, Tcell, project = "Epi_T")
Epi_T <- subset(Epi_T, subset = celltype != "Epi_contaminate")
Epi_T <- subset(Epi_T, subset = celltype != "Shared")

saveRDS(Epi_T, "/data_analysis_2/csl/cyf/SC/6.1.cellchat_Epi_Tcell/Epi_T.Rds")

Epi_T@meta.data$celltype[Epi_T@meta.data$celltype == "EOPC_spe"] <- "EOPC-Epi"
Epi_T@meta.data$celltype[Epi_T@meta.data$celltype == "LOPC_spe"] <- "LOPC-Epi"
Epi_T@meta.data$celltype <- gsub("_", "-", Epi_T@meta.data$celltype)


# 1.CellChat ----------------------------------------------------------------------

library(CellChat)
library(patchwork)

## 1.1.EOPC_Epi/LOPC_Epi vs Tcells --------------------------------------

### EOPC ------------------------------------------------

#### prepare ---------------------------------------

EOPC <- subset(Epi_T, SampleType == "EOPC")
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
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
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
write.csv(df.net, "./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/net_lr_EOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/net_pathway_EOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
EOPC <- cellchat
saveRDS(EOPC, "./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/cellchat_obj_EOPC.Rds")

### LOPC ------------------------------------------------

#### prepare ---------------------------------------

LOPC <- subset(Epi_T, SampleType == "LOPC")
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
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
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
write.csv(df.net, "./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/net_lr_LOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/net_pathway_LOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
LOPC <- cellchat
saveRDS(LOPC, "./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/cellchat_obj_LOPC.Rds")

### integrating objects -------------
PCSC.list <- list(EOPC=EOPC, LOPC=LOPC)
cellchat <- mergeCellChat(PCSC.list, add.names = names(PCSC.list), cell.prefix = TRUE)
saveRDS(cellchat, "./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/cellchat_obj_merge.Rds")

### Visualization -----------------------
##overview barplot
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Overview_number_strength.pdf", p, width = 6, height = 4)

#### difference overview--------------
##overview circle
pdf("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Diff_number_strength_net.pdf")
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

##overview heatmap
pdf("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Diff_number_strength_heatmap.pdf")
par(mfrow = c(1,2))
h1 <- netVisual_heatmap(cellchat, measure = "count")
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
dev.off()

#### compare two groups ----------------------
##compare counts
pdf("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Counts_Compare_net.pdf")
par(mfrow = c(1,2))
weight.max <- getMaxWeight(PCSC.list, attribute = c("idents","count"))
for (i in 1:length(PCSC.list)) {
  netVisual_circle(PCSC.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(PCSC.list)[i]))
}
dev.off()

##compare weight
pdf("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Weight_Compare_net.pdf")
par(mfrow = c(1,2))
weight.max <- getMaxWeight(PCSC.list, attribute = c("idents","weight"))
for (i in 1:length(PCSC.list)) {
  netVisual_circle(PCSC.list[[i]]@net$weight, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Weight of interactions - ", names(PCSC.list)[i]))
}
dev.off()

#compare count in particular cell groups
par(mfrow = c(1,2))
s.cell <- c("CD4+ T cells", "CD8+ T cells", "Monocytes")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))
# save as Counts_Compare_select.pdf 10*6.5

#compare weight in particular cell groups
par(mfrow = c(1,2))
s.cell <- c("CD4+ T cells", "CD8+ T cells", "Monocytes")
weight1 <- cco.list[[1]]@net$weight[s.cell, s.cell]
weight2 <- cco.list[[2]]@net$weight[s.cell, s.cell]
weight.max <- max(max(weight1), max(weight2))
netVisual_circle(weight1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Weight of interactions-", names(cco.list)[1]))
netVisual_circle(weight2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Weight of interactions-", names(cco.list)[2]))
# save as Counts_Compare_select.pdf 10*6.5

#### 保守和特异性信号通路的识别与可视化--------------------------------
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Compare_pathway_strengh.pdf", p, width = 8, height = 15)

#### 流行学习识别差异信号通路(失败) ---------------------------------
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
#netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)

saveRDS(cellchat, "./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/cellchat.rds")

#### 细胞信号模式对比 ----------------------

library(ComplexHeatmap)

##overview
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_all.pdf  10*6

##outgoing
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_outgoing.pdf  10*6

##incoming
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_incoming.pdf  10*6

#### 特定信号通路的对比 ------------------------------
pathways.show <- c("TGFB") 
weight.max <- getMaxWeight(PCSC.list, slot.name = c("netP"), attribute = pathways.show) 

##circle
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(PCSC.list)) {
  netVisual_aggregate(PCSC.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(PCSC.list)[i]))
}
# save as Compare_IL16_net.pdf  10*6.5

##heatmap
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# save as Compare_IL16_heatmap.pdf  12*6.5

##chord
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_chord.pdf  10*6.5

#### 配体-受体对比分析 -----------------------------------
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9),  comparison = c(1, 2), angle.x = 45)
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Compare_Epi_to_other_bubble.pdf", p, width = 10, height = 25)

p <- netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), targets.use = c(3,5),  comparison = c(1, 2), angle.x = 45)
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Compare_other_to_Epi_bubble.pdf", p, width = 10, height = 25)

p1 <- netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Increased signaling in EOPC", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Decreased signaling in EOPC", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Compare_Epi_to_other_regulated.pdf", pc, width = 16, height = 20)

p1 <- netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), targets.use = c(3,5), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Increased signaling in EOPC", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), targets.use = c(3,5), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Decreased signaling in EOPC", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Compare_other_to_cancer_regulated.pdf", pc, width = 16, height = 20)

unique(cellchat@netP$EOPC$pathways)
pdf("./6.1.cellchat_Epi_Tcell/cellchat/celltype_compare/Compare_LR_chord.pdf.pdf")
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(PCSC.list)) {
  netVisual_chord_gene(PCSC.list[[i]], sources.use = c(2,3,4,7), targets.use = c(2,3,4,7), signaling = "VEGF", 
                       lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20,
                       title.name = paste0("VEGF - ", names(PCSC.list)[i]))
}
dev.off()

## 1.2.EOPC/LOPC_group vs Tcell -----------------------------------

Epi <- readRDS("5.1.epithelial_harmony/Epi_annoted.Rds")
Epi@meta.data <- Epi@meta.data[,c(1:8,19,20,25)]
Epi$celltype <- "Epi"
Epi_T_group <- merge(Epi, Tcell, project = "Epi_T_group")
Epi_T_group <- subset(Epi_T_group, subset = celltype != "Epi_contaminate")
Epi_T_group@meta.data$celltype <- gsub("_", "-", Epi_T_group@meta.data$celltype)
table(Epi_T_group$celltype)

### EOPC ------------------------------------------------

#### prepare ---------------------------------------

EOPC <- subset(Epi_T_group, SampleType == "EOPC")
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
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
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
write.csv(df.net, "./6.1.cellchat_Epi_Tcell/cellchat/compare/net_lr_EOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "./6.1.cellchat_Epi_Tcell/cellchat/compare/net_pathway_EOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
EOPC <- cellchat
saveRDS(EOPC, "./6.1.cellchat_Epi_Tcell/cellchat/compare/cellchat_obj_EOPC.Rds")

### LOPC ------------------------------------------------

#### prepare ---------------------------------------

LOPC <- subset(Epi_T_group, SampleType == "LOPC")
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
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
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
write.csv(df.net, "./6.1.cellchat_Epi_Tcell/cellchat/compare/net_lr_LOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "./6.1.cellchat_Epi_Tcell/cellchat/compare/net_pathway_LOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
LOPC <- cellchat
saveRDS(LOPC, "./6.1.cellchat_Epi_Tcell/cellchat/compare/cellchat_obj_LOPC.Rds")

### integrating objects -------------
PCSC.list <- list(EOPC=EOPC, LOPC=LOPC)
cellchat <- mergeCellChat(PCSC.list, add.names = names(PCSC.list), cell.prefix = TRUE)
saveRDS(cellchat, "./6.1.cellchat_Epi_Tcell/cellchat/compare/cellchat_obj_merge.Rds")

### Visualization -----------------------
##overview barplot
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/compare/Overview_number_strength.pdf", p, width = 6, height = 4)

#### difference overview--------------
##overview circle
pdf("./6.1.cellchat_Epi_Tcell/cellchat/compare/Diff_number_strength_net.pdf")
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

##overview heatmap
pdf("./6.1.cellchat_Epi_Tcell/cellchat/compare/Diff_number_strength_heatmap.pdf")
par(mfrow = c(1,2))
h1 <- netVisual_heatmap(cellchat, measure = "count")
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
dev.off()

#### compare two groups ----------------------
##compare counts
pdf("./6.1.cellchat_Epi_Tcell/cellchat/compare/Counts_Compare_net.pdf")
par(mfrow = c(1,2))
weight.max <- getMaxWeight(PCSC.list, attribute = c("idents","count"))
for (i in 1:length(PCSC.list)) {
  netVisual_circle(PCSC.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(PCSC.list)[i]))
}
dev.off()

##compare weight
pdf("./6.1.cellchat_Epi_Tcell/cellchat/compare/Weight_Compare_net.pdf")
par(mfrow = c(1,2))
weight.max <- getMaxWeight(PCSC.list, attribute = c("idents","weight"))
for (i in 1:length(PCSC.list)) {
  netVisual_circle(PCSC.list[[i]]@net$weight, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Weight of interactions - ", names(PCSC.list)[i]))
}
dev.off()

#compare count in particular cell groups
par(mfrow = c(1,2))
s.cell <- c("CD4+ T cells", "CD8+ T cells", "Monocytes")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))
# save as Counts_Compare_select.pdf 10*6.5

#compare weight in particular cell groups
par(mfrow = c(1,2))
s.cell <- c("CD4+ T cells", "CD8+ T cells", "Monocytes")
weight1 <- cco.list[[1]]@net$weight[s.cell, s.cell]
weight2 <- cco.list[[2]]@net$weight[s.cell, s.cell]
weight.max <- max(max(weight1), max(weight2))
netVisual_circle(weight1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Weight of interactions-", names(cco.list)[1]))
netVisual_circle(weight2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Weight of interactions-", names(cco.list)[2]))
# save as Counts_Compare_select.pdf 10*6.5

#### 保守和特异性信号通路的识别与可视化--------------------------------
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/compare/Compare_pathway_strengh.pdf", p, width = 8, height = 15)

#### 流行学习识别差异信号通路(失败) ---------------------------------
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
#netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)

saveRDS(cellchat, "./6.1.cellchat_Epi_Tcell/cellchat/compare/cellchat.rds")

#### 细胞信号模式对比 ----------------------

library(ComplexHeatmap)

##overview
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_all.pdf  10*6

##outgoing
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_outgoing.pdf  10*6

##incoming
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_incoming.pdf  10*6

#### 特定信号通路的对比 ------------------------------
pathways.show <- c("VEGF") 
weight.max <- getMaxWeight(PCSC.list, slot.name = c("netP"), attribute = pathways.show) 

##circle
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(PCSC.list)) {
  netVisual_aggregate(PCSC.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(PCSC.list)[i]))
}
# save as Compare_IL16_net.pdf  10*6.5

##heatmap
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# save as Compare_IL16_heatmap.pdf  12*6.5

##chord
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_chord.pdf  10*6.5

#### 配体-受体对比分析 -----------------------------------
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(1,2,5,6,8),  comparison = c(1, 2), angle.x = 45)
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/compare/Compare_Epi_to_other_bubble.pdf", p, width = 10, height = 25)

p <- netVisual_bubble(cellchat, sources.use = c(1,2,4,5,6,7,8), targets.use = c(3),  comparison = c(1, 2), angle.x = 45)
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/compare/Compare_other_to_Epi_bubble.pdf", p, width = 10, height = 25)

p1 <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(1,2,4,5,6,7,8), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Increased signaling in EOPC", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(1,2,4,5,6,7,8), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Decreased signaling in EOPC", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("./6.1.cellchat_Epi_Tcell/cellchat/compare/Compare_Epi_to_other_regulated.pdf", pc, width = 16, height = 20)

unique(cellchat@netP$EOPC$pathways)
pdf("./6.1.cellchat_Epi_Tcell/cellchat/compare/Compare_LR_chord.pdf.pdf")
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(PCSC.list)) {
  netVisual_chord_gene(PCSC.list[[i]], sources.use = c(2,3,4,7), targets.use = c(2,3,4,7), signaling = "VEGF", 
                       lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20,
                       title.name = paste0("VEGF - ", names(PCSC.list)[i]))
}
dev.off()


#2.cellphoneDB (cancer MP1 sep) -------------------------------------

#2.1.prepare data ------------------------------------------------

PCSC <- readRDS("./4.1.inferCNV/InferCancer/PCSC_inferred.Rds")
Epi_T <- subset(PCSC, inferCNVCellStat %in% c("Cancer", "T-cell"))
Epi <- readRDS("./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
meta_epi <- Epi@meta.data
meta_all <- Epi_T@meta.data
meta_all$celltype <- "UnKnown"
meta_all$celltype[rownames(meta_all) %in% rownames(meta_epi)] <- meta_epi$MP1_OnOff
Tcell <- readRDS("./5.2.Tcell_harmony/Tcell_annotated.Rds")
Tcell <- subset(Tcell, celltype != "Epi_contaminate")
meta_T <- Tcell@meta.data
meta_all$celltype[rownames(meta_all) %in% rownames(meta_T)] <- meta_T$celltype

Epi_T@meta.data <- meta_all
Epi_T <- subset(Epi_T, celltype != "UnKnown")
Epi_T$celltype[Epi_T$celltype == "On"] <- "MP1_On"
Epi_T$celltype[Epi_T$celltype == "Off"] <- "MP1_Off"

# all
Epi_T$SampleCelltype <- paste0(Epi_T$orig.ident, "_", Epi_T$celltype)
saveRDS(Epi_T, "6.1.cellchat_Epi_Tcell/Epi_T_for_cpdb.Rds")
expr <- as.matrix(Epi_T@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(Epi_T@meta.data), celltype = Epi_T@meta.data$SampleCelltype)  
write.table(celltype, paste0("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

# individual
for (i in unique(Epi_T$orig.ident)) {
  Obj <- subset(Epi_T, orig.ident == i)
  expr <- as.matrix(Obj@assays$RNA@counts)
  expr <- data.frame(Gene = rownames(expr), expr)
  write.table(expr, paste0("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/input/expr_", i, ".txt"),
              sep='\t', quote=F, row.names = F)  #表达谱
  celltype <- data.frame(Cell = rownames(Obj@meta.data), celltype = Obj@meta.data$celltype)  
  write.table(celltype, paste0("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/input/celltype_", i, ".txt"), 
              sep='\t', quote=F, row.names = F)  #celltype
}

# general
expr <- as.matrix(Epi_T@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(Epi_T@meta.data), celltype = Epi_T@meta.data$celltype)  
write.table(celltype, paste0("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

#2.2.cellphoneDB analysis ----------------------------------------------

cellphonedb method statistical_analysis --counts-data gene_name --output-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all --threads 20 ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/input/celltype_all.txt ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/input/expr_all.txt
cellphonedb plot dot_plot --means-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/means.txt --pvalues-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/pvalues.txt --output-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all
cellphonedb plot heatmap_plot --pvalues-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/pvalues.txt --output-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/input/celltype_all.txt

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
cellphonedb method statistical_analysis --counts-data gene_name --output-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/output_$i --threads 20 ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/input/celltype_$i.txt ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/input/expr_$i.txt;
done

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
cellphonedb plot dot_plot --means-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/output_$i/means.txt --pvalues-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/output_$i/pvalues.txt --output-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/output_$i;
cellphonedb plot heatmap_plot --pvalues-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/output_$i/pvalues.txt --output-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/output_$i ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_individual/input/celltype_$i.txt;
done

cellphonedb method statistical_analysis --counts-data gene_name --output-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all --threads 20 ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/input/celltype_all.txt ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/input/expr_all.txt
cellphonedb plot dot_plot --means-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all/means.txt --pvalues-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all/pvalues.txt --output-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all
cellphonedb plot heatmap_plot --pvalues-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all/pvalues.txt --output-path ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all ./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/input/celltype_all.txt

#2.3.general plotting --------------------------------------------------------

library(ktplots)

color <- RColorBrewer::brewer.pal(11, "Spectral")

pvals <- read.delim("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all/pvalues.txt", check.names = FALSE)
means <- read.delim("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all/means.txt", check.names = FALSE)
decon <- read.delim("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all/deconvoluted.txt", check.names = FALSE)

# chord plot
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)

mynet <- read.delim("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all/count_network.txt", check.names = FALSE)
mynet %>% filter(count>0) -> mynet

mynet_On <- mynet[-which(mynet$SOURCE == "MP1_Off" | mynet$TARGET == "MP1_Off"),]
mynet_Off <- mynet[-which(mynet$SOURCE == "MP1_On" | mynet$TARGET == "MP1_On"),]

# On
net <- graph_from_data_frame(mynet_On)
karate_groups <- cluster_optimal(net) #统计每个端点的和
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局
E(net)$width <- E(net)$count/10  #根据count值设置边的宽 
for (j in 1: length(unique(mynet_On$SOURCE)) ){ #配置发出端的颜色
  E(net)[map(unique(mynet_On$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet_On$SOURCE)[j],x))
  })%>% unlist()]$color <- color[j]
}

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all/MP1_On_chord.pdf")
plot(net, edge.arrow.size=.2,  #箭头大小设置为0.1
     edge.curved=0.2, # 只是调了曲率这个参数
     vertex.color=color,
     vertex.frame.color="#555555", #圆圈颜色
     vertex.label.color="black", #标签颜色
     layout = coords, #网络布局位点
     vertex.label.cex=1.5) #标签大小 
dev.off()

# Off
net <- graph_from_data_frame(mynet_Off)
karate_groups <- cluster_optimal(net) #统计每个端点的和
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局
E(net)$width <- E(net)$count/10  #根据count值设置边的宽 
for (j in 1: length(unique(mynet_Off$SOURCE)) ){ #配置发出端的颜色
  E(net)[map(unique(mynet_Off$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet_Off$SOURCE)[j],x))
  })%>% unlist()]$color <- color[j]
}

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1_general/output_all/MP1_Off_chord.pdf")
plot(net, edge.arrow.size=.2,  #箭头大小设置为0.1
     edge.curved=0.2, # 只是调了曲率这个参数
     vertex.color=color,
     vertex.frame.color="#555555", #圆圈颜色
     vertex.label.color="black", #标签颜色
     layout = coords, #网络布局位点
     vertex.label.cex=1.5) #标签大小
dev.off()

# 2.4.overall compare — all samples ---------------------------------------------------

library(ktplots)

color <- RColorBrewer::brewer.pal(11, "Spectral")

pvals <- read.delim(paste0("6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/pvalues.txt"), check.names = FALSE)
means <- read.delim(paste0("6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/means.txt"), check.names = FALSE)
decon <- read.delim(paste0("6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/deconvoluted.txt"), check.names = FALSE)
mynet <- read.delim("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/count_network.txt", check.names = FALSE)
mynet %>% filter(count>0) -> mynet

## compare counts between MP1 group 

net_S <- mynet[grep("MP1_", mynet$SOURCE),]
net_S <- net_S[-grep("MP1_", net_S$TARGET),]
net_S <- net_S[which(substr(net_S$SOURCE, 1, 5) == substr(net_S$TARGET, 1, 5)),]

net_S$Tcell <- substring(net_S$TARGET, 7)
net_S$Group <- substring(net_S$SOURCE, 7)

data_plot <- net_S[, c(4,5,3)]
data_plot$Group <- factor(data_plot$Group, levels = c("MP1_On", "MP1_Off"))

# boxplot

P = ggplot(data_plot, aes(x = Tcell, y = count, fill = Group)) +
  scale_fill_manual(values = c("#E64B35AA", "#3C5488AA")) +
  geom_boxplot(width = 0.6, outlier.colour = "grey40", position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Cancer & T-cell interactions") +
  ylab("Interaction Pair Counts") +
  xlab("") + 
  stat_compare_means(method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif",
                     label.x.npc = "center")

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/boxplot_overall_counts.pdf", 
    width = 5, height = 4)
P
dev.off()

## A. Cancer - Treg ----------------------------------------------

## 2.4.1.Cancer - Treg pairs heatmap (pval as heatmap, pval as tag), compare in all samples ----------------------------------------------

pvals_Cancer_Treg <- pvals[, grep(".*MP.*Treg", colnames(pvals))]
pvals_Cancer_Treg <- pvals_Cancer_Treg[, which(substr(colnames(pvals_Cancer_Treg), 1, 5) == substr(colnames(pvals_Cancer_Treg), nchar(colnames(pvals_Cancer_Treg))-9, nchar(colnames(pvals_Cancer_Treg))-5))]
pvals_Cancer_Treg <- cbind(pvals[,1:11], pvals_Cancer_Treg)
means_Cancer_Treg <- means[, grep(".*MP.*Treg", colnames(means))]
means_Cancer_Treg <- means_Cancer_Treg[, which(substr(colnames(means_Cancer_Treg), 1, 5) == substr(colnames(means_Cancer_Treg), nchar(colnames(means_Cancer_Treg))-9, nchar(colnames(means_Cancer_Treg))-5))]
means_Cancer_Treg <- cbind(pvals[,1:11], means_Cancer_Treg)

result_On_Treg <- plot_cpdb("MP1_On", "Treg", Epi_T, "celltype", means_Cancer_Treg, pvals_Cancer_Treg,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)
result_Off_Treg <- plot_cpdb("MP1_Off", "Treg", Epi_T, "celltype", means_Cancer_Treg, pvals_Cancer_Treg,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)

result_On_Treg_pvals <- data.frame(cell = result_On_Treg$Var2, pair = result_On_Treg$Var1, pvals = result_On_Treg$pvals)
result_On_Treg_pvals <- dcast(result_On_Treg_pvals, pair ~ cell)
colnames(result_On_Treg_pvals)[2:ncol(result_On_Treg_pvals)] <- paste0(substr(colnames(result_On_Treg_pvals)[2:ncol(result_On_Treg_pvals)],1,5), "_On")

result_Off_Treg_pvals <- data.frame(cell = result_Off_Treg$Var2, pair = result_Off_Treg$Var1, pvals = result_Off_Treg$pvals)
result_Off_Treg_pvals <- dcast(result_Off_Treg_pvals, pair ~ cell)
colnames(result_Off_Treg_pvals)[2:ncol(result_Off_Treg_pvals)] <- paste0(substr(colnames(result_Off_Treg_pvals)[2:ncol(result_Off_Treg_pvals)],1,5), "_Off")

dim(result_On_Treg_pvals)
dim(result_Off_Treg_pvals)

result_pvals <- merge(result_On_Treg_pvals, result_Off_Treg_pvals, by = "pair", all = T)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1_On")) & result_pvals$"EOPC1_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2_On")) & result_pvals$"EOPC2_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3_On")) & result_pvals$"EOPC3_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4_On")) & result_pvals$"EOPC4_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC1_On")) & result_pvals$"LOPC1_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC2_On")) & result_pvals$"LOPC2_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC3_On")) & result_pvals$"LOPC3_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC4_On")) & result_pvals$"LOPC4_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC5_On")) & result_pvals$"LOPC5_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC6_On")) & result_pvals$"LOPC6_On" < 0.05),]
write.table(result_pvals_filter, file="./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Cancer_Treg_pval_On_CoSig.xls",sep="\t",row.names=F,quote=F)

#plotting
rownames(result_pvals_filter) <- NULL
result_pvals_filter <- column_to_rownames(result_pvals_filter, "pair")

result_pvals_sig <- result_pvals_filter
ss <- result_pvals_sig < 0.01
result_pvals_sig[ss] <- '**'
s <- result_pvals_sig >= 0.01 & result_pvals_sig < 0.05
result_pvals_sig[s] <- '*'
ns <- is.na(result_pvals_sig)
result_pvals_sig[ns] <- 'n.s'

result_pvals_filter <- -log10(result_pvals_filter)
ns <- is.na(result_pvals_filter)
result_pvals_filter[ns] <- 0

result_pvals_filter <- as.matrix(result_pvals_filter)
result_pvals_sig <- as.matrix(result_pvals_sig)

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_Cancer_Treg_pval_On_CoSig.pdf", width = 8, height = 4)
pheatmap(result_pvals_filter,
         scale = "none",
         cluster_row = F, 
         cluster_col = F, 
         border=NA,
         display_numbers = result_pvals_sig,
         fontsize_number = 8, 
         number_color = "black",
         cellwidth = 20, 
         cellheight =20,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         legend = F
)
dev.off()

## 2.4.2.Cancer - Treg pairs heatmap (mean as heatmap, pval as tag), compare in all samples -------------------------------

pvals_Cancer_Treg <- pvals[, which(substr(colnames(pvals), 1, 5) == substr(colnames(pvals), nchar(colnames(pvals))-9, nchar(colnames(pvals))-5))]
pvals_Cancer_Treg <- pvals_Cancer_Treg[, -c(1,2)]
pvals_Cancer_Treg <- pvals_Cancer_Treg[, grep(".*MP.*Treg", colnames(pvals_Cancer_Treg))]
pvals_Cancer_Treg <- cbind(pvals[,1:11], pvals_Cancer_Treg)
means_Cancer_Treg <- means[, which(substr(colnames(means), 1, 5) == substr(colnames(means), nchar(colnames(pvals))-9, nchar(colnames(pvals))-5))]
means_Cancer_Treg <- means_Cancer_Treg[, -c(1,2)]
means_Cancer_Treg <- means_Cancer_Treg[, grep(".*MP.*Treg", colnames(means_Cancer_Treg))]
means_Cancer_Treg <- cbind(pvals[,1:11], means_Cancer_Treg)

result_On_Treg <- plot_cpdb("MP1_On", "Treg", Epi_T, "celltype", means_Cancer_Treg, pvals_Cancer_Treg,
                            split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                            keep_significant_only = T, return_table = T)
result_Off_Treg <- plot_cpdb("MP1_Off", "Treg", Epi_T, "celltype", means_Cancer_Treg, pvals_Cancer_Treg,
                             split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                             keep_significant_only = T, return_table = T)

result_On_Treg_pvals <- data.frame(cell = result_On_Treg$Var2, pair = result_On_Treg$Var1, pvals = result_On_Treg$pvals)
result_On_Treg_pvals <- dcast(result_On_Treg_pvals, pair ~ cell)
colnames(result_On_Treg_pvals)[2:ncol(result_On_Treg_pvals)] <- paste0(substr(colnames(result_On_Treg_pvals)[2:ncol(result_On_Treg_pvals)],1,5), "_On")

result_Off_Treg_pvals <- data.frame(cell = result_Off_Treg$Var2, pair = result_Off_Treg$Var1, pvals = result_Off_Treg$pvals)
result_Off_Treg_pvals <- dcast(result_Off_Treg_pvals, pair ~ cell)
colnames(result_Off_Treg_pvals)[2:ncol(result_Off_Treg_pvals)] <- paste0(substr(colnames(result_Off_Treg_pvals)[2:ncol(result_Off_Treg_pvals)],1,5), "_Off")

result_On_Treg_means <- data.frame(cell = result_On_Treg$Var2, pair = result_On_Treg$Var1, means = result_On_Treg$means)
result_On_Treg_means <- dcast(result_On_Treg_means, pair ~ cell)
colnames(result_On_Treg_means)[2:ncol(result_On_Treg_means)] <- paste0(substr(colnames(result_On_Treg_means)[2:ncol(result_On_Treg_means)],1,5), "_On")

result_Off_Treg_means <- data.frame(cell = result_Off_Treg$Var2, pair = result_Off_Treg$Var1, means = result_Off_Treg$means)
result_Off_Treg_means <- dcast(result_Off_Treg_means, pair ~ cell)
colnames(result_Off_Treg_means)[2:ncol(result_Off_Treg_means)] <- paste0(substr(colnames(result_Off_Treg_means)[2:ncol(result_Off_Treg_means)],1,5), "_Off")

dim(result_On_Treg_pvals)
dim(result_Off_Treg_pvals)
dim(result_On_Treg_means)
dim(result_Off_Treg_means)

result_pvals <- merge(result_On_Treg_pvals, result_Off_Treg_pvals, by = "pair", all = T)
result_means <- merge(result_On_Treg_means, result_Off_Treg_means, by = "pair", all = T)

#pvals tag

rownames(result_pvals) <- NULL
result_pvals <- column_to_rownames(result_pvals, "pair")

result_pvals_sig <- result_pvals
ss <- result_pvals_sig < 0.01
result_pvals_sig[ss] <- '**'
s <- result_pvals_sig >= 0.01 & result_pvals_sig < 0.05
result_pvals_sig[s] <- '*'
ns <- is.na(result_pvals_sig)
result_pvals_sig[ns] <- 'n'
result_pvals_sig <- as.matrix(result_pvals_sig)

#diff

rownames(result_means) <- NULL
result_means <- column_to_rownames(result_means, "pair")
ns <- is.na(result_means)
result_means[ns] <- 0

On <- as.matrix(result_means[, grep("_On", colnames(result_means))])
Off <- as.matrix(result_means[, grep("_Off", colnames(result_means))])

result <- c()
for (i in 1:nrow(On)){
  P <- wilcox.test(On[i,],Off[i,])$p.value
  logFC <- log2(mean(On[i,])/mean(Off[i,]))
  res <- data.frame(Gene = rownames(On)[i], MeanT = mean(On[i,]), MeanC = mean(Off[i,]), 
                    MedT = median(On[i,]), MedC = median(Off[i,]), logFC = logFC, Pvalue = P)
  MedT <- median(On[i,])
  MedC <- median(Off[i,])
  diffMed <- MedT-MedC
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    result <- rbind(result,res)
  }
}
result$fdr <- p.adjust(result$Pvalue,method = "fdr")
colnames(result) <- c("Pair", "MeanOn", "MeanOff", "MedOn", "MedOff", "logFC", "Pvalue", "FDR")

write.table(result, "./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Cancer_Treg_mean_diff.xls",sep="\t",row.names=F,quote=F)

result_sig <- result[result$FDR < 0.05,]
result_sig_up <- result_sig[result_sig$logFC > 0,]
result_sig_down <- result_sig[result_sig$logFC < 0,]

#heatmap up

result_means_up <- result_means[result_sig_up$Pair, ]
result_means_up <- as.matrix(result_means_up)

NMS <- function(x){
  return((x-min(x)) / (max(x)-min(x)))
}
result_means_up <- t(apply(result_means_up, 1, NMS))

result_pvals_sig_up <- result_pvals_sig[result_sig_up$Pair, ]

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_Cancer_Treg_mean_diff_sigUp.pdf", height = 7.5)
pheatmap(result_means_up,
         scale = "none",
         cluster_row = T, 
         cluster_col = F,
         display_numbers = result_pvals_sig_up,
         fontsize_number = 7, 
         number_color = "black",
         cellwidth = 10, 
         cellheight =10,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         #width = 5,
         #height = 5
         main = "Differentially interacted Cancer - Treg pairs"
)
dev.off()

#heatmap down

result_means_down <- result_means[result_sig_down$Pair, ]
result_means_down <- as.matrix(result_means_down)

NMS <- function(x){
  return((x-min(x)) / (max(x)-min(x)))
}
result_means_down <- t(apply(result_means_down, 1, NMS))

result_pvals_sig_down <- result_pvals_sig[result_sig_down$Pair, ]

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_Cancer_Treg_mean_diff_sigDown.pdf", height = 7.5)
pheatmap(result_means_down,
         scale = "none",
         cluster_row = T, 
         cluster_col = F,
         display_numbers = result_pvals_sig_down,
         fontsize_number = 7, 
         number_color = "black",
         cellwidth = 10, 
         cellheight =10,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         #width = 5,
         #height = 5
         main = "Differentially interacted Cancer - Treg pairs"
)
dev.off()

## 2.4.3.Treg - Cancer pairs heatmap (pval as heatmap, pval as tag), compare in all samples ----------------------------------------------

pvals_Treg_Cancer <- pvals[, grep(".*Treg.*MP", colnames(pvals))]
pvals_Treg_Cancer <- pvals_Treg_Cancer[, which(substr(colnames(pvals_Treg_Cancer), 1, 5) == substr(colnames(pvals_Treg_Cancer), 12, 16))]
pvals_Treg_Cancer <- cbind(pvals[,1:11], pvals_Treg_Cancer)
means_Treg_Cancer <- means[, grep(".*Treg.*MP", colnames(means))]
means_Treg_Cancer <- means_Treg_Cancer[, which(substr(colnames(means_Treg_Cancer), 1, 5) == substr(colnames(means_Treg_Cancer), 12, 16))]
means_Treg_Cancer <- cbind(pvals[,1:11], means_Treg_Cancer)

result_Treg_On <- plot_cpdb("Treg", "MP1_On", Epi_T, "celltype", means_Treg_Cancer, pvals_Treg_Cancer,
                            split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                            keep_significant_only = T, return_table = T)
result_Treg_Off <- plot_cpdb("Treg", "MP1_Off", Epi_T, "celltype", means_Treg_Cancer, pvals_Treg_Cancer,
                             split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                             keep_significant_only = T, return_table = T)

result_Treg_On_pvals <- data.frame(cell = result_Treg_On$Var2, pair = result_Treg_On$Var1, pvals = result_Treg_On$pvals)
result_Treg_On_pvals <- dcast(result_Treg_On_pvals, pair ~ cell)
colnames(result_Treg_On_pvals)[2:ncol(result_Treg_On_pvals)] <- paste0(substr(colnames(result_Treg_On_pvals)[2:ncol(result_Treg_On_pvals)],1,5), "_On")

result_Treg_Off_pvals <- data.frame(cell = result_Treg_Off$Var2, pair = result_Treg_Off$Var1, pvals = result_Treg_Off$pvals)
result_Treg_Off_pvals <- dcast(result_Treg_Off_pvals, pair ~ cell)
colnames(result_Treg_Off_pvals)[2:ncol(result_Treg_Off_pvals)] <- paste0(substr(colnames(result_Treg_Off_pvals)[2:ncol(result_Treg_Off_pvals)],1,5), "_Off")

dim(result_Treg_On_pvals)
dim(result_Treg_Off_pvals)

result_pvals <- merge(result_Treg_On_pvals, result_Treg_Off_pvals, by = "pair", all = T)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1_On")) & result_pvals$"EOPC1_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2_On")) & result_pvals$"EOPC2_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3_On")) & result_pvals$"EOPC3_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4_On")) & result_pvals$"EOPC4_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC1_On")) & result_pvals$"LOPC1_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC2_On")) & result_pvals$"LOPC2_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC3_On")) & result_pvals$"LOPC3_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC4_On")) & result_pvals$"LOPC4_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC5_On")) & result_pvals$"LOPC5_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC6_On")) & result_pvals$"LOPC6_On" < 0.05),]
write.table(result_pvals_filter, file="./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Treg_Cancer_pval_On_CoSig.xls",sep="\t",row.names=F,quote=F)

#plotting
rownames(result_pvals_filter) <- NULL
result_pvals_filter <- column_to_rownames(result_pvals_filter, "pair")

result_pvals_sig <- result_pvals_filter
ss <- result_pvals_sig < 0.01
result_pvals_sig[ss] <- '**'
s <- result_pvals_sig >= 0.01 & result_pvals_sig < 0.05
result_pvals_sig[s] <- '*'
ns <- is.na(result_pvals_sig)
result_pvals_sig[ns] <- 'n.s'

result_pvals_filter <- -log10(result_pvals_filter)
ns <- is.na(result_pvals_filter)
result_pvals_filter[ns] <- 0

result_pvals_filter <- as.matrix(result_pvals_filter)
result_pvals_sig <- as.matrix(result_pvals_sig)

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_Treg_Cancer_pval_On_CoSig.pdf", width = 8, height = 4)
pheatmap(result_pvals_filter,
         scale = "none",
         cluster_row = F, 
         cluster_col = F,
         display_numbers = result_pvals_sig,
         fontsize_number = 8, 
         number_color = "black",
         cellwidth = 20, 
         cellheight = 20,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         legend = F
)
dev.off()

## 2.4.4.Treg - Cancer pairs heatmap (mean as heatmap, pval as tag), compare in all samples -------------------------------

pvals_Treg_Cancer <- pvals[, grep(".*Treg.*MP", colnames(pvals))]
pvals_Treg_Cancer <- pvals_Treg_Cancer[, which(substr(colnames(pvals_Treg_Cancer), 1, 5) == substr(colnames(pvals_Treg_Cancer), 12, 16))]
pvals_Treg_Cancer <- cbind(pvals[,1:11], pvals_Treg_Cancer)
means_Treg_Cancer <- means[, grep(".*Treg.*MP", colnames(means))]
means_Treg_Cancer <- means_Treg_Cancer[, which(substr(colnames(means_Treg_Cancer), 1, 5) == substr(colnames(means_Treg_Cancer), 12, 16))]
means_Treg_Cancer <- cbind(pvals[,1:11], means_Treg_Cancer)

result_Treg_On <- plot_cpdb("Treg", "MP1_On", Epi_T, "celltype", means_Treg_Cancer, pvals_Treg_Cancer,
                            split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                            keep_significant_only = T, return_table = T)
result_Treg_Off <- plot_cpdb("Treg", "MP1_Off", Epi_T, "celltype", means_Treg_Cancer, pvals_Treg_Cancer,
                             split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                             keep_significant_only = T, return_table = T)

result_Treg_On_pvals <- data.frame(cell = result_Treg_On$Var2, pair = result_Treg_On$Var1, pvals = result_Treg_On$pvals)
result_Treg_On_pvals <- dcast(result_Treg_On_pvals, pair ~ cell)
colnames(result_Treg_On_pvals)[2:ncol(result_Treg_On_pvals)] <- paste0(substr(colnames(result_Treg_On_pvals)[2:ncol(result_Treg_On_pvals)],1,5), "_On")

result_Treg_Off_pvals <- data.frame(cell = result_Treg_Off$Var2, pair = result_Treg_Off$Var1, pvals = result_Treg_Off$pvals)
result_Treg_Off_pvals <- dcast(result_Treg_Off_pvals, pair ~ cell)
colnames(result_Treg_Off_pvals)[2:ncol(result_Treg_Off_pvals)] <- paste0(substr(colnames(result_Treg_Off_pvals)[2:ncol(result_Treg_Off_pvals)],1,5), "_Off")

result_Treg_On_means <- data.frame(cell = result_Treg_On$Var2, pair = result_Treg_On$Var1, means = result_Treg_On$means)
result_Treg_On_means <- dcast(result_Treg_On_means, pair ~ cell)
colnames(result_Treg_On_means)[2:ncol(result_Treg_On_means)] <- paste0(substr(colnames(result_Treg_On_means)[2:ncol(result_Treg_On_means)],1,5), "_On")

result_Treg_Off_means <- data.frame(cell = result_Treg_Off$Var2, pair = result_Treg_Off$Var1, means = result_Treg_Off$means)
result_Treg_Off_means <- dcast(result_Treg_Off_means, pair ~ cell)
colnames(result_Treg_Off_means)[2:ncol(result_Treg_Off_means)] <- paste0(substr(colnames(result_Treg_Off_means)[2:ncol(result_Treg_Off_means)],1,5), "_Off")

dim(result_Treg_On_pvals)
dim(result_Treg_Off_pvals)
dim(result_Treg_On_means)
dim(result_Treg_Off_means)

result_pvals <- merge(result_Treg_On_pvals, result_Treg_Off_pvals, by = "pair", all = T)
result_means <- merge(result_Treg_On_means, result_Treg_Off_means, by = "pair", all = T)

#pvals tag

rownames(result_pvals) <- NULL
result_pvals <- column_to_rownames(result_pvals, "pair")

result_pvals_sig <- result_pvals
ss <- result_pvals_sig < 0.01
result_pvals_sig[ss] <- '**'
s <- result_pvals_sig >= 0.01 & result_pvals_sig < 0.05
result_pvals_sig[s] <- '*'
ns <- is.na(result_pvals_sig)
result_pvals_sig[ns] <- 'n'
result_pvals_sig <- as.matrix(result_pvals_sig)

#diff

rownames(result_means) <- NULL
result_means <- column_to_rownames(result_means, "pair")
ns <- is.na(result_means)
result_means[ns] <- 0

On <- as.matrix(result_means[, grep("_On", colnames(result_means))])
Off <- as.matrix(result_means[, grep("_Off", colnames(result_means))])

result <- c()
for (i in 1:nrow(On)){
  P <- wilcox.test(On[i,],Off[i,])$p.value
  logFC <- log2(mean(On[i,])/mean(Off[i,]))
  res <- data.frame(Gene = rownames(On)[i], MeanT = mean(On[i,]), MeanC = mean(Off[i,]), 
                    MedT = median(On[i,]), MedC = median(Off[i,]), logFC = logFC, Pvalue = P)
  MedT <- median(On[i,])
  MedC <- median(Off[i,])
  diffMed <- MedT-MedC
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    result <- rbind(result,res)
  }
}
result$fdr <- p.adjust(result$Pvalue,method = "fdr")
colnames(result) <- c("Pair", "MeanOn", "MeanOff", "MedOn", "MedOff", "logFC", "Pvalue", "FDR")

write.table(result, "./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Treg_Cancer_mean_diff.xls",sep="\t",row.names=F,quote=F)

result_sig <- result[result$FDR < 0.05,]
result_sig_up <- result_sig[result_sig$logFC > 0,]
result_sig_down <- result_sig[result_sig$logFC < 0,]

#heatmap up

result_means_up <- result_means[result_sig_up$Pair, ]
result_means_up <- as.matrix(result_means_up)

NMS <- function(x){
  return((x-min(x)) / (max(x)-min(x)))
}
result_means_up <- t(apply(result_means_up, 1, NMS))

result_pvals_sig_up <- result_pvals_sig[result_sig_up$Pair, ]

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_Treg_Cancer_mean_diff_sigUp.pdf", height = 9)
pheatmap(result_means_up,
         scale = "none",
         cluster_row = T, 
         cluster_col = F,
         display_numbers = result_pvals_sig_up,
         fontsize_number = 7, 
         number_color = "black",
         cellwidth = 10, 
         cellheight =10,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         #width = 5,
         #height = 5
         main = "Differentially interacted Treg - Cancer pairs"
)
dev.off()

#heatmap down

result_means_down <- result_means[result_sig_down$Pair, ]
result_means_down <- as.matrix(result_means_down)

NMS <- function(x){
  return((x-min(x)) / (max(x)-min(x)))
}
result_means_down <- t(apply(result_means_down, 1, NMS))

result_pvals_sig_down <- result_pvals_sig[result_sig_down$Pair, ]

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_Treg_Cancer_mean_diff_sigDown.pdf", height = 4.5)
pheatmap(result_means_down,
         scale = "none",
         cluster_row = T, 
         cluster_col = F,
         display_numbers = result_pvals_sig_down,
         fontsize_number = 7, 
         number_color = "black",
         cellwidth = 10, 
         cellheight =10,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         #width = 5,
         #height = 5
         main = "Differentially interacted Treg - Cancer pairs"
)
dev.off()

## B. Cancer - CD8 Tn ----------------------------------------------

## 2.4.5.Cancer - CD8Tn pairs heatmap (pval as heatmap, pval as tag), compare in all samples ----------------------------------------------

pvals_Cancer_CD8Tn <- pvals[, grep(".*MP.*CD8_Tn", colnames(pvals))]
pvals_Cancer_CD8Tn <- pvals_Cancer_CD8Tn[, which(substr(colnames(pvals_Cancer_CD8Tn), 1, 5) == substr(colnames(pvals_Cancer_CD8Tn), nchar(colnames(pvals_Cancer_CD8Tn))-11, nchar(colnames(pvals_Cancer_CD8Tn))-7))]
pvals_Cancer_CD8Tn <- cbind(pvals[,1:11], pvals_Cancer_CD8Tn)
means_Cancer_CD8Tn <- means[, grep(".*MP.*CD8_Tn", colnames(means))]
means_Cancer_CD8Tn <- means_Cancer_CD8Tn[, which(substr(colnames(means_Cancer_CD8Tn), 1, 5) == substr(colnames(means_Cancer_CD8Tn), nchar(colnames(means_Cancer_CD8Tn))-11, nchar(colnames(means_Cancer_CD8Tn))-7))]
means_Cancer_CD8Tn <- cbind(pvals[,1:11], means_Cancer_CD8Tn)

result_On_CD8Tn <- plot_cpdb("MP1_On", "CD8_Tn", Epi_T, "celltype", means_Cancer_CD8Tn, pvals_Cancer_CD8Tn,
                            split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                            keep_significant_only = T, return_table = T)
result_Off_CD8Tn <- plot_cpdb("MP1_Off", "CD8_Tn", Epi_T, "celltype", means_Cancer_CD8Tn, pvals_Cancer_CD8Tn,
                             split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                             keep_significant_only = T, return_table = T)

result_On_CD8Tn_pvals <- data.frame(cell = result_On_CD8Tn$Var2, pair = result_On_CD8Tn$Var1, pvals = result_On_CD8Tn$pvals)
result_On_CD8Tn_pvals <- dcast(result_On_CD8Tn_pvals, pair ~ cell)
colnames(result_On_CD8Tn_pvals)[2:ncol(result_On_CD8Tn_pvals)] <- paste0(substr(colnames(result_On_CD8Tn_pvals)[2:ncol(result_On_CD8Tn_pvals)],1,5), "_On")

result_Off_CD8Tn_pvals <- data.frame(cell = result_Off_CD8Tn$Var2, pair = result_Off_CD8Tn$Var1, pvals = result_Off_CD8Tn$pvals)
result_Off_CD8Tn_pvals <- dcast(result_Off_CD8Tn_pvals, pair ~ cell)
colnames(result_Off_CD8Tn_pvals)[2:ncol(result_Off_CD8Tn_pvals)] <- paste0(substr(colnames(result_Off_CD8Tn_pvals)[2:ncol(result_Off_CD8Tn_pvals)],1,5), "_Off")

dim(result_On_CD8Tn_pvals)
dim(result_Off_CD8Tn_pvals)

result_pvals <- merge(result_On_CD8Tn_pvals, result_Off_CD8Tn_pvals, by = "pair", all = T)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1_On")) & result_pvals$"EOPC1_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2_On")) & result_pvals$"EOPC2_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3_On")) & result_pvals$"EOPC3_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4_On")) & result_pvals$"EOPC4_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC1_On")) & result_pvals$"LOPC1_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC2_On")) & result_pvals$"LOPC2_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC3_On")) & result_pvals$"LOPC3_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC4_On")) & result_pvals$"LOPC4_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC5_On")) & result_pvals$"LOPC5_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC6_On")) & result_pvals$"LOPC6_On" < 0.05),]
write.table(result_pvals_filter, file="./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Cancer_CD8Tn_pval_On_CoSig.xls",sep="\t",row.names=F,quote=F)

#plotting
rownames(result_pvals_filter) <- NULL
result_pvals_filter <- column_to_rownames(result_pvals_filter, "pair")

result_pvals_sig <- result_pvals_filter
ss <- result_pvals_sig < 0.01
result_pvals_sig[ss] <- '**'
s <- result_pvals_sig >= 0.01 & result_pvals_sig < 0.05
result_pvals_sig[s] <- '*'
ns <- is.na(result_pvals_sig)
result_pvals_sig[ns] <- 'n.s'

result_pvals_filter <- -log10(result_pvals_filter)
ns <- is.na(result_pvals_filter)
result_pvals_filter[ns] <- 0

result_pvals_filter <- as.matrix(result_pvals_filter)
result_pvals_sig <- as.matrix(result_pvals_sig)

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_Cancer_CD8Tn_pval_On_CoSig.pdf", width = 8, height = 4)
pheatmap(result_pvals_filter,
         scale = "none",
         cluster_row = F, 
         cluster_col = F, 
         border=NA,
         display_numbers = result_pvals_sig,
         fontsize_number = 8, 
         number_color = "black",
         cellwidth = 20, 
         cellheight =20,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         legend = F
)
dev.off()

## 2.4.6.Cancer - CD8Tn pairs heatmap (mean as heatmap, pval as tag), compare in all samples -------------------------------

pvals_Cancer_CD8Tn <- pvals[, grep(".*MP.*CD8_Tn", colnames(pvals))]
pvals_Cancer_CD8Tn <- pvals_Cancer_CD8Tn[, which(substr(colnames(pvals_Cancer_CD8Tn), 1, 5) == substr(colnames(pvals_Cancer_CD8Tn), nchar(colnames(pvals_Cancer_CD8Tn))-11, nchar(colnames(pvals_Cancer_CD8Tn))-7))]
pvals_Cancer_CD8Tn <- cbind(pvals[,1:11], pvals_Cancer_CD8Tn)
means_Cancer_CD8Tn <- means[, grep(".*MP.*CD8_Tn", colnames(means))]
means_Cancer_CD8Tn <- means_Cancer_CD8Tn[, which(substr(colnames(means_Cancer_CD8Tn), 1, 5) == substr(colnames(means_Cancer_CD8Tn), nchar(colnames(means_Cancer_CD8Tn))-11, nchar(colnames(means_Cancer_CD8Tn))-7))]
means_Cancer_CD8Tn <- cbind(pvals[,1:11], means_Cancer_CD8Tn)

result_On_CD8Tn <- plot_cpdb("MP1_On", "CD8_Tn", Epi_T, "celltype", means_Cancer_CD8Tn, pvals_Cancer_CD8Tn,
                             split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                             keep_significant_only = T, return_table = T)
result_Off_CD8Tn <- plot_cpdb("MP1_Off", "CD8_Tn", Epi_T, "celltype", means_Cancer_CD8Tn, pvals_Cancer_CD8Tn,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)

result_On_CD8Tn_pvals <- data.frame(cell = result_On_CD8Tn$Var2, pair = result_On_CD8Tn$Var1, pvals = result_On_CD8Tn$pvals)
result_On_CD8Tn_pvals <- dcast(result_On_CD8Tn_pvals, pair ~ cell)
colnames(result_On_CD8Tn_pvals)[2:ncol(result_On_CD8Tn_pvals)] <- paste0(substr(colnames(result_On_CD8Tn_pvals)[2:ncol(result_On_CD8Tn_pvals)],1,5), "_On")

result_Off_CD8Tn_pvals <- data.frame(cell = result_Off_CD8Tn$Var2, pair = result_Off_CD8Tn$Var1, pvals = result_Off_CD8Tn$pvals)
result_Off_CD8Tn_pvals <- dcast(result_Off_CD8Tn_pvals, pair ~ cell)
colnames(result_Off_CD8Tn_pvals)[2:ncol(result_Off_CD8Tn_pvals)] <- paste0(substr(colnames(result_Off_CD8Tn_pvals)[2:ncol(result_Off_CD8Tn_pvals)],1,5), "_Off")

result_On_CD8Tn_means <- data.frame(cell = result_On_CD8Tn$Var2, pair = result_On_CD8Tn$Var1, means = result_On_CD8Tn$means)
result_On_CD8Tn_means <- dcast(result_On_CD8Tn_means, pair ~ cell)
colnames(result_On_CD8Tn_means)[2:ncol(result_On_CD8Tn_means)] <- paste0(substr(colnames(result_On_CD8Tn_means)[2:ncol(result_On_CD8Tn_means)],1,5), "_On")

result_Off_CD8Tn_means <- data.frame(cell = result_Off_CD8Tn$Var2, pair = result_Off_CD8Tn$Var1, means = result_Off_CD8Tn$means)
result_Off_CD8Tn_means <- dcast(result_Off_CD8Tn_means, pair ~ cell)
colnames(result_Off_CD8Tn_means)[2:ncol(result_Off_CD8Tn_means)] <- paste0(substr(colnames(result_Off_CD8Tn_means)[2:ncol(result_Off_CD8Tn_means)],1,5), "_Off")

dim(result_On_CD8Tn_pvals)
dim(result_Off_CD8Tn_pvals)
dim(result_On_CD8Tn_means)
dim(result_Off_CD8Tn_means)

result_pvals <- merge(result_On_CD8Tn_pvals, result_Off_CD8Tn_pvals, by = "pair", all = T)
result_means <- merge(result_On_CD8Tn_means, result_Off_CD8Tn_means, by = "pair", all = T)

#pvals tag

rownames(result_pvals) <- NULL
result_pvals <- column_to_rownames(result_pvals, "pair")

result_pvals_sig <- result_pvals
ss <- result_pvals_sig < 0.01
result_pvals_sig[ss] <- '**'
s <- result_pvals_sig >= 0.01 & result_pvals_sig < 0.05
result_pvals_sig[s] <- '*'
ns <- is.na(result_pvals_sig)
result_pvals_sig[ns] <- 'n'
result_pvals_sig <- as.matrix(result_pvals_sig)

#diff

rownames(result_means) <- NULL
result_means <- column_to_rownames(result_means, "pair")
ns <- is.na(result_means)
result_means[ns] <- 0

On <- as.matrix(result_means[, grep("_On", colnames(result_means))])
Off <- as.matrix(result_means[, grep("_Off", colnames(result_means))])

result <- c()
for (i in 1:nrow(On)){
  P <- wilcox.test(On[i,],Off[i,])$p.value
  logFC <- log2(mean(On[i,])/mean(Off[i,]))
  res <- data.frame(Gene = rownames(On)[i], MeanT = mean(On[i,]), MeanC = mean(Off[i,]), 
                    MedT = median(On[i,]), MedC = median(Off[i,]), logFC = logFC, Pvalue = P)
  MedT <- median(On[i,])
  MedC <- median(Off[i,])
  diffMed <- MedT-MedC
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    result <- rbind(result,res)
  }
}
result$fdr <- p.adjust(result$Pvalue,method = "fdr")
colnames(result) <- c("Pair", "MeanOn", "MeanOff", "MedOn", "MedOff", "logFC", "Pvalue", "FDR")

write.table(result, "./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Cancer_CD8Tn_mean_diff.xls",sep="\t",row.names=F,quote=F)

result_sig <- result[result$FDR < 0.05,]
result_sig_up <- result_sig[result_sig$logFC > 0,]
result_sig_down <- result_sig[result_sig$logFC < 0,]

#heatmap up

result_means_up <- result_means[result_sig_up$Pair, ]
result_means_up <- as.matrix(result_means_up)

NMS <- function(x){
  return((x-min(x)) / (max(x)-min(x)))
}
result_means_up <- t(apply(result_means_up, 1, NMS))

result_pvals_sig_up <- result_pvals_sig[result_sig_up$Pair, ]

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_Cancer_CD8_Tn_mean_diff_sigUp.pdf", height = 7.5)
pheatmap(result_means_up,
         scale = "none",
         cluster_row = T, 
         cluster_col = F,
         display_numbers = result_pvals_sig_up,
         fontsize_number = 7, 
         number_color = "black",
         cellwidth = 10, 
         cellheight =10,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         #width = 5,
         #height = 5
         main = "Differentially interacted Cancer - CD8+ Tn pairs"
)
dev.off()

#heatmap down

result_means_down <- result_means[result_sig_down$Pair, ]
result_means_down <- as.matrix(result_means_down)

NMS <- function(x){
  return((x-min(x)) / (max(x)-min(x)))
}
result_means_down <- t(apply(result_means_down, 1, NMS))

result_pvals_sig_down <- result_pvals_sig[result_sig_down$Pair, ]

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_CD8_Tn_Treg_mean_diff_sigDown.pdf", height = 7.5)
pheatmap(result_means_down,
         scale = "none",
         cluster_row = T, 
         cluster_col = F,
         display_numbers = result_pvals_sig_down,
         fontsize_number = 7, 
         number_color = "black",
         cellwidth = 10, 
         cellheight =10,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         #width = 5,
         #height = 5
         main = "Differentially interacted Cancer - CD8+ Tn pairs"
)
dev.off()

## 2.4.7.Treg - Cancer pairs heatmap (pval as heatmap, pval as tag), compare in all samples ----------------------------------------------

pvals_CD8Tn_Cancer <- pvals[, grep(".*CD8_Tn.*MP", colnames(pvals))]
pvals_CD8Tn_Cancer <- pvals_CD8Tn_Cancer[, which(substr(colnames(pvals_CD8Tn_Cancer), 1, 5) == substr(colnames(pvals_CD8Tn_Cancer), 14, 18))]
pvals_CD8Tn_Cancer <- cbind(pvals[,1:11], pvals_CD8Tn_Cancer)
means_CD8Tn_Cancer <- means[, grep(".*CD8_Tn.*MP", colnames(means))]
means_CD8Tn_Cancer <- means_CD8Tn_Cancer[, which(substr(colnames(means_CD8Tn_Cancer), 1, 5) == substr(colnames(means_CD8Tn_Cancer), 14, 18))]
means_CD8Tn_Cancer <- cbind(pvals[,1:11], means_CD8Tn_Cancer)

result_CD8Tn_On <- plot_cpdb("CD8_Tn", "MP1_On", Epi_T, "celltype", means_CD8Tn_Cancer, pvals_CD8Tn_Cancer,
                            split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                            keep_significant_only = T, return_table = T)
result_CD8Tn_Off <- plot_cpdb("CD8_Tn", "MP1_Off", Epi_T, "celltype", means_CD8Tn_Cancer, pvals_CD8Tn_Cancer,
                             split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                             keep_significant_only = T, return_table = T)

result_CD8Tn_On_pvals <- data.frame(cell = result_CD8Tn_On$Var2, pair = result_CD8Tn_On$Var1, pvals = result_CD8Tn_On$pvals)
result_CD8Tn_On_pvals <- dcast(result_CD8Tn_On_pvals, pair ~ cell)
colnames(result_CD8Tn_On_pvals)[2:ncol(result_CD8Tn_On_pvals)] <- paste0(substr(colnames(result_CD8Tn_On_pvals)[2:ncol(result_CD8Tn_On_pvals)],1,5), "_On")

result_CD8Tn_Off_pvals <- data.frame(cell = result_CD8Tn_Off$Var2, pair = result_CD8Tn_Off$Var1, pvals = result_CD8Tn_Off$pvals)
result_CD8Tn_Off_pvals <- dcast(result_CD8Tn_Off_pvals, pair ~ cell)
colnames(result_CD8Tn_Off_pvals)[2:ncol(result_CD8Tn_Off_pvals)] <- paste0(substr(colnames(result_CD8Tn_Off_pvals)[2:ncol(result_CD8Tn_Off_pvals)],1,5), "_Off")

dim(result_CD8Tn_On_pvals)
dim(result_CD8Tn_Off_pvals)

result_pvals <- merge(result_CD8Tn_On_pvals, result_CD8Tn_Off_pvals, by = "pair", all = T)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1_On")) & result_pvals$"EOPC1_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2_On")) & result_pvals$"EOPC2_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3_On")) & result_pvals$"EOPC3_On" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4_On")) & result_pvals$"EOPC4_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC1_On")) & result_pvals$"LOPC1_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC2_On")) & result_pvals$"LOPC2_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC3_On")) & result_pvals$"LOPC3_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC4_On")) & result_pvals$"LOPC4_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC5_On")) & result_pvals$"LOPC5_On" < 0.05 &
                                            (!is.na(result_pvals$"LOPC6_On")) & result_pvals$"LOPC6_On" < 0.05),]
write.table(result_pvals_filter, file="./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/CD8Tn_Cancer_pval_On_CoSig.xls",sep="\t",row.names=F,quote=F)

#plotting
rownames(result_pvals_filter) <- NULL
result_pvals_filter <- column_to_rownames(result_pvals_filter, "pair")

result_pvals_sig <- result_pvals_filter
ss <- result_pvals_sig < 0.01
result_pvals_sig[ss] <- '**'
s <- result_pvals_sig >= 0.01 & result_pvals_sig < 0.05
result_pvals_sig[s] <- '*'
ns <- is.na(result_pvals_sig)
result_pvals_sig[ns] <- 'n.s'

result_pvals_filter <- -log10(result_pvals_filter)
ns <- is.na(result_pvals_filter)
result_pvals_filter[ns] <- 0

result_pvals_filter <- as.matrix(result_pvals_filter)
result_pvals_sig <- as.matrix(result_pvals_sig)

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_CD8Tn_Cancer_pval_On_CoSig.pdf", width = 8, height = 4)
pheatmap(result_pvals_filter,
         scale = "none",
         cluster_row = F, 
         cluster_col = F,
         display_numbers = result_pvals_sig,
         fontsize_number = 8, 
         number_color = "black",
         cellwidth = 20, 
         cellheight = 20,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         legend = F
)
dev.off()

## 2.4.8.CD8Tn - Cancer pairs heatmap (mean as heatmap, pval as tag), compare in all samples -------------------------------

pvals_CD8Tn_Cancer <- pvals[, grep(".*CD8_Tn.*MP", colnames(pvals))]
pvals_CD8Tn_Cancer <- pvals_CD8Tn_Cancer[, which(substr(colnames(pvals_CD8Tn_Cancer), 1, 5) == substr(colnames(pvals_CD8Tn_Cancer), 14, 18))]
pvals_CD8Tn_Cancer <- cbind(pvals[,1:11], pvals_CD8Tn_Cancer)
means_CD8Tn_Cancer <- means[, grep(".*CD8_Tn.*MP", colnames(means))]
means_CD8Tn_Cancer <- means_CD8Tn_Cancer[, which(substr(colnames(means_CD8Tn_Cancer), 1, 5) == substr(colnames(means_CD8Tn_Cancer), 14, 18))]
means_CD8Tn_Cancer <- cbind(pvals[,1:11], means_CD8Tn_Cancer)

result_CD8Tn_On <- plot_cpdb("CD8_Tn", "MP1_On", Epi_T, "celltype", means_CD8Tn_Cancer, pvals_CD8Tn_Cancer,
                             split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                             keep_significant_only = T, return_table = T)
result_CD8Tn_Off <- plot_cpdb("CD8_Tn", "MP1_Off", Epi_T, "celltype", means_CD8Tn_Cancer, pvals_CD8Tn_Cancer,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)

result_CD8Tn_On_pvals <- data.frame(cell = result_CD8Tn_On$Var2, pair = result_CD8Tn_On$Var1, pvals = result_CD8Tn_On$pvals)
result_CD8Tn_On_pvals <- dcast(result_CD8Tn_On_pvals, pair ~ cell)
colnames(result_CD8Tn_On_pvals)[2:ncol(result_CD8Tn_On_pvals)] <- paste0(substr(colnames(result_CD8Tn_On_pvals)[2:ncol(result_CD8Tn_On_pvals)],1,5), "_On")

result_CD8Tn_Off_pvals <- data.frame(cell = result_CD8Tn_Off$Var2, pair = result_CD8Tn_Off$Var1, pvals = result_CD8Tn_Off$pvals)
result_CD8Tn_Off_pvals <- dcast(result_CD8Tn_Off_pvals, pair ~ cell)
colnames(result_CD8Tn_Off_pvals)[2:ncol(result_CD8Tn_Off_pvals)] <- paste0(substr(colnames(result_CD8Tn_Off_pvals)[2:ncol(result_CD8Tn_Off_pvals)],1,5), "_Off")

result_CD8Tn_On_means <- data.frame(cell = result_CD8Tn_On$Var2, pair = result_CD8Tn_On$Var1, means = result_CD8Tn_On$means)
result_CD8Tn_On_means <- dcast(result_CD8Tn_On_means, pair ~ cell)
colnames(result_CD8Tn_On_means)[2:ncol(result_CD8Tn_On_means)] <- paste0(substr(colnames(result_CD8Tn_On_means)[2:ncol(result_CD8Tn_On_means)],1,5), "_On")

result_CD8Tn_Off_means <- data.frame(cell = result_CD8Tn_Off$Var2, pair = result_CD8Tn_Off$Var1, means = result_CD8Tn_Off$means)
result_CD8Tn_Off_means <- dcast(result_CD8Tn_Off_means, pair ~ cell)
colnames(result_CD8Tn_Off_means)[2:ncol(result_CD8Tn_Off_means)] <- paste0(substr(colnames(result_CD8Tn_Off_means)[2:ncol(result_CD8Tn_Off_means)],1,5), "_Off")

dim(result_CD8Tn_On_pvals)
dim(result_CD8Tn_Off_pvals)
dim(result_CD8Tn_On_means)
dim(result_CD8Tn_Off_means)

result_pvals <- merge(result_CD8Tn_On_pvals, result_CD8Tn_Off_pvals, by = "pair", all = T)
result_means <- merge(result_CD8Tn_On_means, result_CD8Tn_Off_means, by = "pair", all = T)

#pvals tag

rownames(result_pvals) <- NULL
result_pvals <- column_to_rownames(result_pvals, "pair")

result_pvals_sig <- result_pvals
ss <- result_pvals_sig < 0.01
result_pvals_sig[ss] <- '**'
s <- result_pvals_sig >= 0.01 & result_pvals_sig < 0.05
result_pvals_sig[s] <- '*'
ns <- is.na(result_pvals_sig)
result_pvals_sig[ns] <- 'n'
result_pvals_sig <- as.matrix(result_pvals_sig)

#diff

rownames(result_means) <- NULL
result_means <- column_to_rownames(result_means, "pair")
ns <- is.na(result_means)
result_means[ns] <- 0

On <- as.matrix(result_means[, grep("_On", colnames(result_means))])
Off <- as.matrix(result_means[, grep("_Off", colnames(result_means))])

result <- c()
for (i in 1:nrow(On)){
  P <- wilcox.test(On[i,],Off[i,])$p.value
  logFC <- log2(mean(On[i,])/mean(Off[i,]))
  res <- data.frame(Gene = rownames(On)[i], MeanT = mean(On[i,]), MeanC = mean(Off[i,]), 
                    MedT = median(On[i,]), MedC = median(Off[i,]), logFC = logFC, Pvalue = P)
  MedT <- median(On[i,])
  MedC <- median(Off[i,])
  diffMed <- MedT-MedC
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    result <- rbind(result,res)
  }
}
result$fdr <- p.adjust(result$Pvalue,method = "fdr")
colnames(result) <- c("Pair", "MeanOn", "MeanOff", "MedOn", "MedOff", "logFC", "Pvalue", "FDR")

write.table(result, "./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/CD8Tn_Cancer_mean_diff.xls",sep="\t",row.names=F,quote=F)

result_sig <- result[result$FDR < 0.05,]
result_sig_up <- result_sig[result_sig$logFC > 0,]
result_sig_down <- result_sig[result_sig$logFC < 0,]

#heatmap up

result_means_up <- result_means[result_sig_up$Pair, ]
result_means_up <- as.matrix(result_means_up)

NMS <- function(x){
  return((x-min(x)) / (max(x)-min(x)))
}
result_means_up <- t(apply(result_means_up, 1, NMS))

result_pvals_sig_up <- result_pvals_sig[result_sig_up$Pair, ]

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_CD8Tn_Cancer_mean_diff_sigUp.pdf", height = 9)
pheatmap(result_means_up,
         scale = "none",
         cluster_row = T, 
         cluster_col = F,
         display_numbers = result_pvals_sig_up,
         fontsize_number = 7, 
         number_color = "black",
         cellwidth = 10, 
         cellheight =10,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         #width = 5,
         #height = 5
         main = "Differentially interacted CD8+ Tn - Cancer pairs"
)
dev.off()

#heatmap down

result_means_down <- result_means[result_sig_down$Pair, ]
result_means_down <- as.matrix(result_means_down)

NMS <- function(x){
  return((x-min(x)) / (max(x)-min(x)))
}
result_means_down <- t(apply(result_means_down, 1, NMS))

result_pvals_sig_down <- result_pvals_sig[result_sig_down$Pair, ]

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_AllSamples/Heatmap_CD8Tn_Cancer_mean_diff_sigDown.pdf", height = 4.5)
pheatmap(result_means_down,
         scale = "none",
         cluster_row = T, 
         cluster_col = F,
         display_numbers = result_pvals_sig_down,
         fontsize_number = 7, 
         number_color = "black",
         cellwidth = 10, 
         cellheight =10,
         color=c(colorRampPalette(colors = c("white","#E64B35AA"))(50)),
         border_color="black",
         #width = 5,
         #height = 5
         main = "Differentially interacted CD8+ Tn - Cancer pairs"
)
dev.off()

# 2.5.overall compare — EO vs LO (four groups) ---------------------------------------------------

library(ktplots)

color <- RColorBrewer::brewer.pal(11, "Spectral")

pvals <- read.delim(paste0("6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/pvalues.txt"), check.names = FALSE)
means <- read.delim(paste0("6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/means.txt"), check.names = FALSE)
decon <- read.delim(paste0("6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/deconvoluted.txt"), check.names = FALSE)
mynet <- read.delim("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/count_network.txt", check.names = FALSE)
mynet %>% filter(count>0) -> mynet

## compare counts between MP1 group 

net_S <- mynet[grep("MP1_", mynet$SOURCE),]
net_S <- net_S[-grep("MP1_", net_S$TARGET),]
net_S <- net_S[which(substr(net_S$SOURCE, 1, 5) == substr(net_S$TARGET, 1, 5)),]

#net_S <- net_S[c(which(substr(net_S$SOURCE, 1, 2) == "EO" & substr(net_S$SOURCE, nchar(net_S$SOURCE)-1, nchar(net_S$SOURCE)) == "On"), 
#                 which(substr(net_S$SOURCE, 1, 2) == "LO" & substr(net_S$SOURCE, nchar(net_S$SOURCE)-2, nchar(net_S$SOURCE)) == "Off")), ]

net_S$Tcell <- substring(net_S$TARGET, 7)
net_S$Group <- substring(net_S$SOURCE, 7)
net_S$Group[which(substr(net_S$SOURCE, 1, 2) == "EO")] <- paste0("EOPC_", net_S$Group[which(substr(net_S$SOURCE, 1, 2) == "EO")])
net_S$Group[which(substr(net_S$SOURCE, 1, 2) == "LO")] <- paste0("LOPC_", net_S$Group[which(substr(net_S$SOURCE, 1, 2) == "LO")])

data_plot <- net_S[, c(4,5,3)]

# boxplot

P = ggplot(data_plot, aes(x = Tcell, y = count, fill = Group)) +
  #scale_fill_manual(values = c("#E64B35AA", "#3C5488AA")) +
  geom_boxplot(width = 0.6, outlier.colour = "grey40", position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Cancer & T-cell interactions") +
  ylab("Interaction Pair Counts") +
  xlab("") + 
  stat_compare_means(method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif",
                     label.x.npc = "center")

pdf("./6.1.cellchat_Epi_Tcell/cellphoneDB_MP1/output_all/On_Off_SampleTypeSplit/boxplot_overall_counts.pdf", 
    width = 5, height = 4)
P
dev.off()








