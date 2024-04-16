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
library(ggpubr)

Epi <- readRDS("5.1.epithelial_harmony_TryTwo/Epi_annoted.Rds")
PCSC <- readRDS("4.1.inferCNV/InferCancer/PCSC_inferred.Rds")

metaAll <- data.frame(cbind(rownames(PCSC@meta.data), PCSC@meta.data$inferCNVCellStat))
colnames(metaAll) <- c("ID", "type")
metaCancer <- data.frame(cbind(rownames(Epi@meta.data), Epi@meta.data$celltype))
colnames(metaCancer) <- c("ID", "type")

meta <- left_join(metaAll, metaCancer, by = "ID")
meta$type.x[!is.na(meta$type.y)] <- meta$type.y[!is.na(meta$type.y)]

Epi_all <- PCSC
Epi_all@meta.data$celltype <- meta$type.x

Epi_all@meta.data$celltype[Epi_all@meta.data$celltype == "EOPC_spe"] <- "EOPC-Epi"
Epi_all@meta.data$celltype[Epi_all@meta.data$celltype == "LOPC_spe"] <- "LOPC-Epi"
Epi_all@meta.data$celltype <- gsub("_", "-", Epi_all@meta.data$celltype)

Epi_all <- subset(Epi_all, celltype != "Other-Epi")
Epi_all <- subset(Epi_all, celltype != "Normal")
Epi_all <- subset(Epi_all, celltype != "Shared")

saveRDS(Epi_all, "6.0.cellchat_Epi_all/Epi_all.Rds")

# 1.CellChat ----------------------------------------------------------------------

library(CellChat)
library(patchwork)

## 1.1.EOPC_Epi/LOPC_Epi vs other cells --------------------------------------

### EOPC ------------------------------------------------

#### prepare ---------------------------------------

EOPC <- subset(Epi_all, SampleType == "EOPC")
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
write.csv(df.net, "./6.0.cellchat_Epi_all/cellchat/celltype_compare/net_lr_EOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "./6.0.cellchat_Epi_all/cellchat/celltype_compare/net_pathway_EOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
EOPC <- cellchat
saveRDS(EOPC, "./6.0.cellchat_Epi_all/cellchat/celltype_compare/cellchat_obj_EOPC.Rds")

### LOPC ------------------------------------------------

#### prepare ---------------------------------------

LOPC <- subset(Epi_all, SampleType == "LOPC")
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
write.csv(df.net, "./6.0.cellchat_Epi_all/cellchat/celltype_compare/net_lr_LOPC.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "./6.0.cellchat_Epi_all/cellchat/celltype_compare/net_pathway_LOPC.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
LOPC <- cellchat
saveRDS(LOPC, "./6.0.cellchat_Epi_all/cellchat/celltype_compare/cellchat_obj_LOPC.Rds")

### integrating objects -------------
PCSC.list <- list(EOPC=EOPC, LOPC=LOPC)
cellchat <- mergeCellChat(PCSC.list, add.names = names(PCSC.list), cell.prefix = TRUE)
saveRDS(cellchat, "./6.0.cellchat_Epi_all/cellchat/celltype_compare/cellchat_obj_merge.Rds")

### Visualization -----------------------
##overview barplot
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Overview_number_strength.pdf", p, width = 6, height = 4)

#### difference overview--------------
##overview circle
pdf("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Diff_number_strength_net.pdf")
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

##overview heatmap
pdf("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Diff_number_strength_heatmap.pdf")
par(mfrow = c(1,2))
h1 <- netVisual_heatmap(cellchat, measure = "count")
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
dev.off()

#### compare two groups ----------------------
##compare counts
pdf("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Counts_Compare_net.pdf")
par(mfrow = c(1,2))
weight.max <- getMaxWeight(PCSC.list, attribute = c("idents","count"))
for (i in 1:length(PCSC.list)) {
  netVisual_circle(PCSC.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(PCSC.list)[i]))
}
dev.off()

##compare weight
pdf("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Weight_Compare_net.pdf")
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
ggsave("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Compare_pathway_strengh.pdf", p, width = 8, height = 15)

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

saveRDS(cellchat, "./6.0.cellchat_Epi_all/cellchat/celltype_compare/cellchat.rds")

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
pathways.show <- c("TGFb") 
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
for (i in 1:length(PCSC.list)) {
  netVisual_aggregate(PCSC.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(PCSC.list)[i]))
}
# save as Compare_IL16_chord.pdf  10*6.5

#### 配体-受体对比分析 -----------------------------------
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9),  comparison = c(1, 2), angle.x = 45)
ggsave("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Compare_Epi_to_other_bubble.pdf", p, width = 10, height = 25)

p <- netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), targets.use = c(3,5),  comparison = c(1, 2), angle.x = 45)
ggsave("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Compare_other_to_Epi_bubble.pdf", p, width = 10, height = 25)

p1 <- netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Increased signaling in EOPC", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Decreased signaling in EOPC", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Compare_Epi_to_other_regulated.pdf", pc, width = 16, height = 20)

p1 <- netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), targets.use = c(3,5), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Increased signaling in EOPC", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), targets.use = c(3,5), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Decreased signaling in EOPC", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Compare_other_to_cancer_regulated.pdf", pc, width = 16, height = 20)

unique(cellchat@netP$EOPC$pathways)
pdf("./6.0.cellchat_Epi_all/cellchat/celltype_compare/Compare_LR_chord.pdf.pdf")
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(PCSC.list)) {
  netVisual_chord_gene(PCSC.list[[i]], sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9), signaling = "TGFb", 
                       lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20,
                       title.name = paste0("TGFb - ", names(PCSC.list)[i]))
}
dev.off()

## 1.2.TGF_H_Epi/TGF_L_Epi vs other cells --------------------------------------

## 1.2.1.all ------------------------------------------

### data processing ----------------------------------------------
Epi <- readRDS("5.1.epithelial_harmony_TryTwo/Epi_annoted.Rds")

library(irGSEA)
TGFB_score <- irGSEA.score(object = Epi, assay = "RNA", slot = "data", 
                           seeds = 123, ncores = 20, min.cells = 3, 
                           min.feature = 0, custom = F, geneset = NULL, 
                           msigdb = T, species = "Homo sapiens", category = "H",  
                           subcategory = NULL, geneid = "symbol",
                           method = c("ssgsea"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                           kcdf = 'Gaussian')
DefaultAssay(TGFB_score) <- "ssgsea"
TGFB_value <- as.data.frame(t(TGFB_score@assays$ssgsea@data))
colnames(TGFB_value) <- gsub("-", "_", colnames(TGFB_value))
median <- median(TGFB_value$HALLMARK_TGF_BETA_SIGNALING)
TGFB_value$group[TGFB_value$HALLMARK_TGF_BETA_SIGNALING > median] <- "TGFB_High"
TGFB_value$group[TGFB_value$HALLMARK_TGF_BETA_SIGNALING <= median] <- "TGFB_Low"
TGFB_group <- as.matrix(TGFB_value$group)
rownames(TGFB_group) <- rownames(TGFB_value)
colnames(TGFB_group) <- "TGFB_group"

library(tidyverse)
Epi@meta.data <- merge(Epi@meta.data, TGFB_group, by = 0)
Epi@meta.data <- column_to_rownames(Epi@meta.data, "Row.names")

Epi@meta.data$celltype <- Epi@meta.data$TGFB_group
Epi@meta.data <- Epi@meta.data[,c(1:8,19)]

PCSC <- readRDS("4.1.inferCNV/InferCancer/PCSC_inferred.Rds")
PCSC@meta.data$celltype <- PCSC@meta.data$inferCNVCellStat
PCSC@meta.data <- PCSC@meta.data[,c(1:8,19)]

metaAll <- data.frame(cbind(rownames(PCSC@meta.data), PCSC@meta.data$celltype))
colnames(metaAll) <- c("ID", "type")
metaCancer <- data.frame(cbind(rownames(Epi@meta.data), Epi@meta.data$celltype))
colnames(metaCancer) <- c("ID", "type")

meta <- left_join(metaAll, metaCancer, by = "ID")
meta$type.x[!is.na(meta$type.y)] <- meta$type.y[!is.na(meta$type.y)]

Epi_all <- PCSC
Epi_all@meta.data$celltype <- meta$type.x

Epi_all@meta.data$celltype <- gsub("_", "-", Epi_all@meta.data$celltype)
Epi_all <- subset(Epi_all, celltype != "Other-Epi")
Epi_all <- subset(Epi_all, celltype != "Normal")

saveRDS(Epi_all, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/Epi_all.Rds")

library(CellChat)
library(patchwork)

### prepare ---------------------------------------

data.input <- GetAssayData(Epi_all, slot = "data")
meta = Epi_all@meta.data
unique(meta$celltype)

## create CellChat

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

## database
CellChatDB <- CellChatDB.human
unique(CellChatDB$interaction$pathway_name)
#showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)   # This step is necessary even if using the whole database
future::plan("multicore", workers = 20)     # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

### infer cell communications ------------------------------------------

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/net_lr.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/net_pathway.csv")

saveRDS(cellchat, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/cellchat_obj_all.Rds")

df.net <- subsetCommunication(cellchat)
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) 
#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.


### visualization ------------------------------------------------------

#### overview-------------------------------------
cellchat <- aggregateNet(cellchat)
pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/net_number_strength.pdf")

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

#### overview split by cell type-----------------------
mat <- cellchat@net$count
pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/net_number_individual.pdf")
par(mfrow = c(3,3), xpd=TRUE, mar=c(0,1.4,0,1.4))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

mat <- cellchat@net$weight
pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/net_strength_individual.pdf")
par(mfrow = c(3,3), xpd=TRUE, mar=c(0,1.4,0,1.4))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

#### individual pathway-------------------------------
pathways.show <- c("CXCL") 

##Hierarchy plot
levels(cellchat@idents)
vertex.receiver = c(3,5)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy", vertex.receiver = vertex.receiver)

##Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

##Chord plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

##Heatmap
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

#### 计算配体受体对选定信号通路的贡献值 ---------------------------------
pathways.show <- c("TGFb")

netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE) #提取对TGFb有贡献的所有配体受体 

##Hierarchy plot
LR.show <- pairLR.TGFb[1,] 
vertex.receiver = c(2,3,4,6,8)
netVisual_aggregate(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "hierarchy", vertex.receiver = vertex.receiver)

##Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

##Chord plot
netVisual_aggregate(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

##Heatmap
netVisual_heatmap(cellchat, signaling = pathways.show, pairLR.use = LR.show, color.heatmap = "Reds")

#### 批量保存每个信号通路的互作结果 ------------------------------

# Access all the signaling pathways showing significant communications将所有信号通路找出来
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
#levels(cellchat@idents)
#vertex.receiver = c(1,2,4,6) #不画层次图就不需要这一步

dir.create("./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/all_pathways_com_circle") #创建文件夹保存批量画图结果
setwd("./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"), layout = "circle") #绘制网络图
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
         plot=gg, width = 5, height = 2.5, units = 'in', dpi = 300)
}
setwd("../../../../")

#### 多个配体-受体介导的细胞互作关系可视化 ------------------------

#气泡图（全部配体受体）
levels(cellchat@idents)
#需要指定受体细胞和配体细胞
p = netVisual_bubble(cellchat, sources.use = c(8,9), 
                     targets.use = c(1,2,3,4,5,6,7), remove.isolate = FALSE)
ggsave("./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/interaction_Cancer_to_T.pdf", p, width = 6, height = 12)

levels(cellchat@idents)
#需要指定受体细胞和配体细胞
p = netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), 
                     targets.use = c(3,5), remove.isolate = FALSE)
ggsave("./6.0.cellchat_Epi_all/cellchat/TGFB_group/All/interaction_T_to_Cancer.pdf", p, width = 6, height = 12)

#气泡图（指定通路）
#比如制定CCL和CXCL这两个信号通路
netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6,7,8), targets.use = c(2), 
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)

#参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）

## Plot the signaling gene expression distribution
p = plotGeneExpression(cellchat, signaling = "TGFb")
#ggsave("TGFb_GeneExpression_vln.pdf", p, width = 8, height = 8)
p = plotGeneExpression(cellchat, signaling = "TGFb", type = "dot")
#ggsave("TGFb_GeneExpression_dot.pdf", p, width = 8, height = 6)

##1.2.2.compare of TGF_H/TGF_L to all -----------------------------------------------------------------

Epi_all <- readRDS("./6.0.cellchat_Epi_all/cellchat/TGFB_group/Epi_all.Rds")
TGF_H <- subset(Epi_all, celltype != "TGFB-Low")
TGF_L <- subset(Epi_all, celltype != "TGFB-High")

TGF_H$celltype[TGF_H$celltype == "TGFB-High"] <- "Epi"
TGF_L$celltype[TGF_L$celltype == "TGFB-Low"] <- "Epi"

### TGFB_H ------------------------------------------------

#### prepare ---------------------------------------

data.input <- GetAssayData(TGF_H, slot = "data")
meta = TGF_H@meta.data
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
future::plan("multicore", workers = 20)     # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

#### infer cell communications ------------------------------------------

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/net_lr_TGFB_H.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/net_pathway_TGFB_H.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
EOPC <- cellchat
saveRDS(EOPC, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/cellchat_obj_TGFB_H.Rds")

### TGF_L ------------------------------------------------

#### prepare ---------------------------------------

data.input <- GetAssayData(TGF_L, slot = "data")
meta = TGF_L@meta.data
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
future::plan("multicore", workers = 20)     # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

#### infer cell communications ------------------------------------------

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/net_lr_TGFB_L.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/net_pathway_TGFB_L.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
LOPC <- cellchat
saveRDS(LOPC, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/cellchat_obj_TGFB_L.Rds")

### integrating objects -------------
PCSC.list <- list(TGF_H=TGF_H, TGF_L=TGF_L)
cellchat <- mergeCellChat(PCSC.list, add.names = names(PCSC.list), cell.prefix = TRUE)
saveRDS(cellchat, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/cellchat_obj_merge.Rds")

### Visualization -----------------------
##overview barplot
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Overview_number_strength.pdf", p, width = 6, height = 4)

#### difference overview--------------
##overview circle
pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Diff_number_strength_net.pdf")
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

##overview heatmap
pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Diff_number_strength_heatmap.pdf")
par(mfrow = c(1,2))
h1 <- netVisual_heatmap(cellchat, measure = "count")
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
dev.off()

#### compare two groups ----------------------
##compare counts
pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Counts_Compare_net.pdf")
par(mfrow = c(1,2))
weight.max <- getMaxWeight(PCSC.list, attribute = c("idents","count"))
for (i in 1:length(PCSC.list)) {
  netVisual_circle(PCSC.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(PCSC.list)[i]))
}
dev.off()

##compare weight
pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Weight_Compare_net.pdf")
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
ggsave("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Compare_pathway_strengh.pdf", p, width = 8, height = 15)

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

saveRDS(cellchat, "./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/cellchat.rds")

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
pathways.show <- c("TGFb") 
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
for (i in 1:length(PCSC.list)) {
  netVisual_aggregate(PCSC.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(PCSC.list)[i]))
}
# save as Compare_IL16_chord.pdf  10*6.5

#### 配体-受体对比分析 -----------------------------------
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(1,2,4,5,6,7,8), comparison = c(1, 2), angle.x = 45)
ggsave("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Compare_Epi_to_other_bubble.pdf", p, width = 10, height = 25)

p <- netVisual_bubble(cellchat, sources.use = 3, targets.use = 2, comparison = c(1, 2), angle.x = 45)
ggsave("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Compare_Epi_to_Endo_bubble.pdf", p, width = 10, height = 25)

Epi <- readRDS("5.1.epithelial_harmony_TryTwo/Epi_annoted.Rds")
Idents(Epi) <- "celltype"

pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Epi_Endo_VEGF/VEGF_expression.pdf")
par(mfrow = c(2,1), xpd=TRUE)
VlnPlot(Epi, "VEGFB")
VlnPlot(Epi, "VEGFA")
dev.off()

Endo <- readRDS("./5.3.Endo_harmony/Endo_annoted.Rds")
Idents(Endo) <- "celltype"
pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Epi_Endo_VEGF/VEGFR_expression.pdf")
VlnPlot(Endo, "FLT1")
VlnPlot(Endo, "KDR")
dev.off()

p <- netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), targets.use = c(3,5),  comparison = c(1, 2), angle.x = 45)
ggsave("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Compare_other_to_Epi_bubble.pdf", p, width = 10, height = 25)

p1 <- netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Increased signaling in EOPC", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Decreased signaling in EOPC", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Compare_Epi_to_other_regulated.pdf", pc, width = 16, height = 20)

p1 <- netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), targets.use = c(3,5), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Increased signaling in EOPC", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(1,2,4,6,7,8,9), targets.use = c(3,5), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Decreased signaling in EOPC", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Compare_other_to_cancer_regulated.pdf", pc, width = 16, height = 20)

unique(cellchat@netP$EOPC$pathways)
pdf("./6.0.cellchat_Epi_all/cellchat/TGFB_group/compare/Compare_LR_chord.pdf.pdf")
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(PCSC.list)) {
  netVisual_chord_gene(PCSC.list[[i]], sources.use = c(3,5), targets.use = c(1,2,4,6,7,8,9), signaling = "TGFb", 
                       lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20,
                       title.name = paste0("TGFb - ", names(PCSC.list)[i]))
}
dev.off()


#2.cellphoneDB -------------------------------------------------------------

#2.1.prepare data ------------------------------------------------

PCSC <- readRDS("./4.1.inferCNV/InferCancer/PCSC_inferred.Rds")
Epi_all <- subset(PCSC, inferCNVCellStat != "Normal")
Epi_all <- subset(Epi_all, inferCNVCellStat != "Other_Epi")

for (i in unique(Epi_all$orig.ident)) {
  Obj <- subset(Epi_all, orig.ident == i)
  expr <- as.matrix(Obj@assays$RNA@counts)
  expr <- data.frame(Gene = rownames(expr), expr)
  write.table(expr, paste0("./6.0.cellchat_Epi_all/cellphoneDB/input/expr_", i, ".txt"),
              sep='\t', quote=F, row.names = F)  #表达谱
  celltype <- data.frame(Cell = rownames(Obj@meta.data), celltype = Obj@meta.data$inferCNVCellStat)  
  write.table(celltype, paste0("./6.0.cellchat_Epi_all/cellphoneDB/input/celltype_", i, ".txt"), 
              sep='\t', quote=F, row.names = F)  #celltype
}

Epi_all$celltype <- paste0(Epi_all$orig.ident, "_", Epi_all$inferCNVCellStat)
expr <- as.matrix(Epi_all@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("./6.0.cellchat_Epi_all/cellphoneDB/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(Epi_all@meta.data), celltype = Epi_all@meta.data$celltype)  
write.table(celltype, paste0("./6.0.cellchat_Epi_all/cellphoneDB/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

#2.2.cellphoneDB analysis ----------------------------------------------

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
mkdir output_$i;
done

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
cellphonedb method statistical_analysis --counts-data gene_name --output-path ./6.0.cellchat_Epi_all/cellphoneDB/output_$i --threads 20 ./6.0.cellchat_Epi_all/cellphoneDB/input/celltype_$i.txt ./6.0.cellchat_Epi_all/cellphoneDB/input/expr_$i.txt;
done

cellphonedb method statistical_analysis --counts-data gene_name --output-path ./6.0.cellchat_Epi_all/cellphoneDB/output_all --threads 20 ./6.0.cellchat_Epi_all/cellphoneDB/input/celltype_all.txt ./6.0.cellchat_Epi_all/cellphoneDB/input/expr_all.txt;

for i in EOPC1 EOPC2 EOPC3 EOPC4 LOPC1 LOPC2 LOPC3 LOPC4 LOPC5 LOPC6;do
cellphonedb plot dot_plot --means-path ./6.0.cellchat_Epi_all/cellphoneDB/output_$i/means.txt --pvalues-path ./6.0.cellchat_Epi_all/cellphoneDB/output_$i/pvalues.txt --output-path ./6.0.cellchat_Epi_all/cellphoneDB/output_$i;
cellphonedb plot heatmap_plot --pvalues-path ./6.0.cellchat_Epi_all/cellphoneDB/output_$i/pvalues.txt --output-path ./6.0.cellchat_Epi_all/cellphoneDB/output_$i ./6.0.cellchat_Epi_all/cellphoneDB/input/celltype_$i.txt;
done

cellphonedb plot dot_plot --means-path ./6.0.cellchat_Epi_all/cellphoneDB/output_all/means.txt --pvalues-path ./6.0.cellchat_Epi_all/cellphoneDB/output_all/pvalues.txt --output-path ./6.0.cellchat_Epi_all/cellphoneDB/output_all;
cellphonedb plot heatmap_plot --pvalues-path ./6.0.cellchat_Epi_all/cellphoneDB/output_all/pvalues.txt --output-path ./6.0.cellchat_Epi_all/cellphoneDB/output_all ./6.0.cellchat_Epi_all/cellphoneDB/input/celltype_all.txt;

#2.3.interpet and lotting ---------------------------------------------------

library(ktplots)

color <- brewer.pal(8, "Set1")

for (i in c("EOPC1", "EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  pvals <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_", i, "/pvalues.txt"), check.names = FALSE)
  assign(paste0("pvals_", i), pvals)
  means <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_", i, "/means.txt"), check.names = FALSE)
  assign(paste0("means_", i), means)
  decon <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_", i, "/deconvoluted.txt"), check.names = FALSE)
  assign(paste0("decon_", i), decon)
  mynet <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_", i, "/count_network.txt"), check.names = FALSE)
  assign(paste0("mynet_", i), mynet)
}

pvals <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_all/pvalues.txt"), check.names = FALSE)
means <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_all/means.txt"), check.names = FALSE)
decon <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_all/deconvoluted.txt"), check.names = FALSE)

# chord plot
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)

for (i in c("EOPC1", "EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  mynet <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_", i, "/count_network.txt"), check.names = FALSE)
  mynet %>% filter(count>0) -> mynet
  net <- graph_from_data_frame(mynet)
  karate_groups <- cluster_optimal(net) #统计每个端点的和
  coords <- layout_in_circle(net, order =
                               order(membership(karate_groups)))  # 设置网络布局
  E(net)$width <- E(net)$count/10  #根据count值设置边的宽 
  for (j in 1: length(unique(mynet$SOURCE)) ){ #配置发出端的颜色
    E(net)[map(unique(mynet$SOURCE),function(x) {
      get.edge.ids(net,vp = c(unique(mynet$SOURCE)[j],x))
    })%>% unlist()]$color <- color[j]
  }
  pdf(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_", i, "/chord.pdf"))
  plot(net, edge.arrow.size=.2,  #箭头大小设置为0.1
       edge.curved=0.2, # 只是调了曲率这个参数
       vertex.color=color,
       vertex.frame.color="#555555", #圆圈颜色
       vertex.label.color="black", #标签颜色
       layout = coords, #网络布局位点
       vertex.label.cex=1.5) #标签大小 
  dev.off()
}

# chord plot (Cancer Endo Fibro SMC Myeloid)
for (i in c("EOPC1", "EOPC2", "EOPC3", "EOPC4", "LOPC1", "LOPC2", "LOPC3", "LOPC4", "LOPC5", "LOPC6")) {
  mynet <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_", i, "/count_network.txt"), check.names = FALSE)
  mynet <- mynet[-which(mynet$SOURCE %in% c("T-cell", "B-cell", "Mast")),]
  mynet <- mynet[-which(mynet$TARGET %in% c("T-cell", "B-cell", "Mast")),]
  mynet %>% filter(count>0) -> mynet
  net <- graph_from_data_frame(mynet)
  karate_groups <- cluster_optimal(net) #统计每个端点的和
  coords <- layout_in_circle(net, order =
                               order(membership(karate_groups)))  # 设置网络布局
  E(net)$width <- E(net)$count/10  #根据count值设置边的宽 
  for (j in 1: length(unique(mynet$SOURCE)) ){ #配置发出端的颜色
    E(net)[map(unique(mynet$SOURCE),function(x) {
      get.edge.ids(net,vp = c(unique(mynet$SOURCE)[j],x))
    })%>% unlist()]$color <- color[j]
  }
  pdf(paste0("6.0.cellchat_Epi_all/cellphoneDB/output_", i, "/chord_CEFSM.pdf"))
  plot(net, edge.arrow.size=.2,  #箭头大小设置为0.1
       edge.curved=0.2, # 只是调了曲率这个参数
       vertex.color=color,
       vertex.frame.color="#555555", #圆圈颜色
       vertex.label.color="black", #标签颜色
       layout = coords, #网络布局位点
       vertex.label.cex=1.5) #标签大小 
  dev.off()
}

# Cancer - Endo pairs heatmap (mean as heatmap, pval as tag)

pvals_Cancer_Endo <- pvals[, c(1:11, grep(".*Cancer.*Endothelial", colnames(pvals)))]
means_Cancer_Endo <- means[, c(1:11, grep(".*Cancer.*Endothelial", colnames(means)))]
result <- plot_cpdb("Cancer", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                    split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                    keep_significant_only = F, return_table = T)
result_pvals <- data.frame(cell = result$Var2, pair = result$Var1, means = result$pvals)
result_pvals <- dcast(result_pvals, pair ~ cell)
result_means <- data.frame(cell = result$Var2, pair = result$Var1, means = result$means)
result_means <- dcast(result_means, pair ~ cell)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1_Cancer-EOPC1_Endothelial")) & result_pvals$"EOPC1_Cancer-EOPC1_Endothelial" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2_Cancer-EOPC2_Endothelial")) & result_pvals$"EOPC2_Cancer-EOPC2_Endothelial" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3_Cancer-EOPC3_Endothelial")) & result_pvals$"EOPC3_Cancer-EOPC3_Endothelial" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4_Cancer-EOPC4_Endothelial")) & result_pvals$"EOPC4_Cancer-EOPC4_Endothelial" < 0.05),]
write.table(result_pvals_filter, file="./6.0.cellchat_Epi_all/cellphoneDB/Cancer_Endo_diff/sig.xls",sep="\t",row.names=F,quote=F)

#means
result_means_filter <- result_means[which(result_means$pair %in% result_pvals_filter$pair),]

#plotting
rownames(result_pvals_filter) <- NULL
result_pvals_filter <- column_to_rownames(result_pvals_filter, "pair")
rownames(result_means_filter) <- NULL
result_means_filter <- column_to_rownames(result_means_filter, "pair")

ss <- result_pvals_filter < 0.01
result_pvals_filter[ss] <- '**'
s <- result_pvals_filter >= 0.01 & result_pvals_filter < 0.05
result_pvals_filter[s] <- '*'
ns <- is.na(result_pvals_filter)
result_pvals_filter[ns] <- 'n.s'

result_pvals_filter <- as.matrix(result_pvals_filter)
result_means_filter <- as.matrix(result_means_filter)

pdf("./6.0.cellchat_Epi_all/cellphoneDB/Cancer_Endo_diff/Cancer_Endo_heatmap_EOcoexp_mean.pdf", width = 8, height = 11)
bk <- c(seq(-2,2,by=0.05))
pheatmap(result_means_filter,
         scale = "row",
         cluster_row = T, 
         cluster_col = F, 
         border=NA,
         display_numbers = result_pvals_filter,
         fontsize_number = 8, 
         number_color = "black",
         cellwidth = 20, 
         cellheight =20,
         color=c(colorRampPalette(colors = c("white","#E64B35"))(length(bk))),
         legend_breaks=seq(-2,2,1),
         breaks=bk
)
dev.off()

# Cancer - Endo pairs heatmap (pval as heatmap, mean as tag)

pvals_Cancer_Endo <- pvals[, c(1:11, grep(".*Cancer.*Endothelial", colnames(pvals)))]
means_Cancer_Endo <- means[, c(1:11, grep(".*Cancer.*Endothelial", colnames(means)))]
result <- plot_cpdb("Cancer", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                    split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                    keep_significant_only = F, return_table = T)
result_pvals <- data.frame(cell = result$Var2, pair = result$Var1, means = result$pvals)
result_pvals <- dcast(result_pvals, pair ~ cell)
colnames(result_pvals)[2:ncol(result_pvals)] <- substr(colnames(result_pvals)[2:ncol(result_pvals)], 1, 5)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1")) & result_pvals$"EOPC1" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2")) & result_pvals$"EOPC2" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3")) & result_pvals$"EOPC3" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4")) & result_pvals$"EOPC4" < 0.05),]

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

pdf("./6.0.cellchat_Epi_all/cellphoneDB/Cancer_Endo_diff/Cancer_Endo_heatmap_EOcoexp_pval.pdf", width = 6, height = 9)
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
         color=c(colorRampPalette(colors = c("white","salmon"))(50)),
         border_color="black"
)
dev.off()





#3.cellphoneDB (cancer celltype sep) -------------------------------------

#3.1.prepare data ------------------------------------------------

PCSC <- readRDS("./4.1.inferCNV/InferCancer/PCSC_inferred.Rds")
Epi_all <- subset(PCSC, inferCNVCellStat != "Normal")
Epi_all <- subset(Epi_all, inferCNVCellStat != "Other_Epi")
Epi <- readRDS("./5.1.epithelial_harmony_TryTwo/Epi_annoted.Rds")
meta_epi <- Epi@meta.data
meta_all <- Epi_all@meta.data
meta_all$inferCNVCellStat[rownames(meta_all) %in% rownames(meta_epi)] <- meta_epi$celltype
Epi_all@meta.data <- meta_all

Epi_all$celltype <- paste0(Epi_all$orig.ident, "_", Epi_all$inferCNVCellStat)
expr <- as.matrix(Epi_all@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(Epi_all@meta.data), celltype = Epi_all@meta.data$celltype)  
write.table(celltype, paste0("./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

#3.2.cellphoneDB analysis ----------------------------------------------

cellphonedb method statistical_analysis --counts-data gene_name --output-path ./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/output_all --threads 20 ./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/input/celltype_all.txt ./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/input/expr_all.txt

#3.3.interpet and lotting ---------------------------------------------------

library(ktplots)

color <- brewer.pal(8, "Set1")

pvals <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/output_all/pvalues.txt"), check.names = FALSE)
means <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/output_all/means.txt"), check.names = FALSE)
decon <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/output_all/deconvoluted.txt"), check.names = FALSE)

# Cancer - Endo pairs heatmap (mean as heatmap, pval as tag)

pvals_Cancer_Endo <- pvals[, c(1:11, grep(".*Cancer.*Endothelial", colnames(pvals)))]
means_Cancer_Endo <- means[, c(1:11, grep(".*Cancer.*Endothelial", colnames(means)))]
result <- plot_cpdb("Cancer", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                    split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                    keep_significant_only = F, return_table = T)
result_pvals <- data.frame(cell = result$Var2, pair = result$Var1, means = result$pvals)
result_pvals <- dcast(result_pvals, pair ~ cell)
result_means <- data.frame(cell = result$Var2, pair = result$Var1, means = result$means)
result_means <- dcast(result_means, pair ~ cell)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1_Cancer-EOPC1_Endothelial")) & result_pvals$"EOPC1_Cancer-EOPC1_Endothelial" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2_Cancer-EOPC2_Endothelial")) & result_pvals$"EOPC2_Cancer-EOPC2_Endothelial" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3_Cancer-EOPC3_Endothelial")) & result_pvals$"EOPC3_Cancer-EOPC3_Endothelial" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4_Cancer-EOPC4_Endothelial")) & result_pvals$"EOPC4_Cancer-EOPC4_Endothelial" < 0.05),]
write.table(result_pvals_filter, file="./6.0.cellchat_Epi_all/cellphoneDB/Cancer_Endo_diff/sig.xls",sep="\t",row.names=F,quote=F)

#means
result_means_filter <- result_means[which(result_means$pair %in% result_pvals_filter$pair),]

#plotting
rownames(result_pvals_filter) <- NULL
result_pvals_filter <- column_to_rownames(result_pvals_filter, "pair")
rownames(result_means_filter) <- NULL
result_means_filter <- column_to_rownames(result_means_filter, "pair")

ss <- result_pvals_filter < 0.01
result_pvals_filter[ss] <- '**'
s <- result_pvals_filter >= 0.01 & result_pvals_filter < 0.05
result_pvals_filter[s] <- '*'
ns <- is.na(result_pvals_filter)
result_pvals_filter[ns] <- 'n.s'

result_pvals_filter <- as.matrix(result_pvals_filter)
result_means_filter <- as.matrix(result_means_filter)

pdf("./6.0.cellchat_Epi_all/cellphoneDB/Cancer_Endo_diff/Cancer_Endo_heatmap_EOcoexp_mean.pdf", width = 8, height = 11)
bk <- c(seq(-2,2,by=0.05))
pheatmap(result_means_filter,
         scale = "row",
         cluster_row = T, 
         cluster_col = F, 
         border=NA,
         display_numbers = result_pvals_filter,
         fontsize_number = 8, 
         number_color = "black",
         cellwidth = 20, 
         cellheight =20,
         color=c(colorRampPalette(colors = c("white","#E64B35"))(length(bk))),
         legend_breaks=seq(-2,2,1),
         breaks=bk
)
dev.off()

# Cancer - Endo pairs heatmap (pval as heatmap, mean as tag)

pvals_Cancer_Endo <- pvals[, c(1:11, grep(".*spe.*Endothelial", colnames(pvals)))]
means_Cancer_Endo <- means[, c(1:11, grep(".*spe.*Endothelial", colnames(means)))]

result_EOPC_Endo <- plot_cpdb("EOPC_spe", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)
result_LOPC_Endo <- plot_cpdb("LOPC_spe", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)

result_EOPC_Endo_pvals <- data.frame(cell = result_EOPC_Endo$Var2, pair = result_EOPC_Endo$Var1, pvals = result_EOPC_Endo$pvals)
result_EOPC_Endo_pvals <- dcast(result_EOPC_Endo_pvals, pair ~ cell)
result_EOPC_Endo_pvals <- result_EOPC_Endo_pvals[, c(1:5)]

result_LOPC_Endo_pvals <- data.frame(cell = result_LOPC_Endo$Var2, pair = result_LOPC_Endo$Var1, pvals = result_LOPC_Endo$pvals)
result_LOPC_Endo_pvals <- dcast(result_LOPC_Endo_pvals, pair ~ cell)
result_LOPC_Endo_pvals <- result_LOPC_Endo_pvals[, c(1,6:11)]

dim(result_EOPC_Endo_pvals)
dim(result_LOPC_Endo_pvals)

result_pvals <- merge(result_EOPC_Endo_pvals, result_LOPC_Endo_pvals, by = "pair", all = T)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1")) & result_pvals$"EOPC1" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2")) & result_pvals$"EOPC2" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3")) & result_pvals$"EOPC3" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4")) & result_pvals$"EOPC4" < 0.05),]
write.table(result_pvals_filter, file="./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/Cancer_Endo_diff/sig.xls",sep="\t",row.names=F,quote=F)

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

pdf("./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/Cancer_Endo_diff/Cancer_Endo_heatmap_EOcoexp_pval.pdf", width = 6, height = 10)
pheatmap(result_pvals_filter,
         scale = "none",
         cluster_row = F, 
         cluster_col = F,
         display_numbers = result_pvals_sig,
         fontsize_number = 8, 
         number_color = "black",
         cellwidth = 20, 
         cellheight =20,
         color=c(colorRampPalette(colors = c("white","salmon"))(50)),
         border_color="black"
)
dev.off()


#4.cellphoneDB (cancer MP1 sep) -------------------------------------

#3.1.prepare data ------------------------------------------------

PCSC <- readRDS("./4.1.inferCNV/InferCancer/PCSC_inferred.Rds")
Epi_all <- subset(PCSC, inferCNVCellStat != "Normal")
Epi_all <- subset(Epi_all, inferCNVCellStat != "Other_Epi")
Epi <- readRDS("./5.1.epithelial_harmony_TryTwo/CytoTRACE/Epi_noRib_CytoTrace_result.Rds")
meta_epi <- Epi@meta.data
meta_all <- Epi_all@meta.data
meta_all$inferCNVCellStat[rownames(meta_all) %in% rownames(meta_epi)] <- meta_epi$celltype
Epi_all@meta.data <- meta_all

Epi_all$celltype <- paste0(Epi_all$orig.ident, "_", Epi_all$inferCNVCellStat)
expr <- as.matrix(Epi_all@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(Epi_all@meta.data), celltype = Epi_all@meta.data$celltype)  
write.table(celltype, paste0("./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

#3.2.cellphoneDB analysis ----------------------------------------------

cellphonedb method statistical_analysis --counts-data gene_name --output-path ./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/output_all --threads 20 ./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/input/celltype_all.txt ./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/input/expr_all.txt

#3.3.interpet and lotting ---------------------------------------------------

library(ktplots)

color <- brewer.pal(8, "Set1")

pvals <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/output_all/pvalues.txt"), check.names = FALSE)
means <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/output_all/means.txt"), check.names = FALSE)
decon <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/output_all/deconvoluted.txt"), check.names = FALSE)

# Cancer - Endo pairs heatmap (mean as heatmap, pval as tag)

pvals_Cancer_Endo <- pvals[, c(1:11, grep(".*Cancer.*Endothelial", colnames(pvals)))]
means_Cancer_Endo <- means[, c(1:11, grep(".*Cancer.*Endothelial", colnames(means)))]
result <- plot_cpdb("Cancer", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                    split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                    keep_significant_only = F, return_table = T)
result_pvals <- data.frame(cell = result$Var2, pair = result$Var1, means = result$pvals)
result_pvals <- dcast(result_pvals, pair ~ cell)
result_means <- data.frame(cell = result$Var2, pair = result$Var1, means = result$means)
result_means <- dcast(result_means, pair ~ cell)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1_Cancer-EOPC1_Endothelial")) & result_pvals$"EOPC1_Cancer-EOPC1_Endothelial" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2_Cancer-EOPC2_Endothelial")) & result_pvals$"EOPC2_Cancer-EOPC2_Endothelial" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3_Cancer-EOPC3_Endothelial")) & result_pvals$"EOPC3_Cancer-EOPC3_Endothelial" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4_Cancer-EOPC4_Endothelial")) & result_pvals$"EOPC4_Cancer-EOPC4_Endothelial" < 0.05),]
write.table(result_pvals_filter, file="./6.0.cellchat_Epi_all/cellphoneDB/Cancer_Endo_diff/sig.xls",sep="\t",row.names=F,quote=F)

#means
result_means_filter <- result_means[which(result_means$pair %in% result_pvals_filter$pair),]

#plotting
rownames(result_pvals_filter) <- NULL
result_pvals_filter <- column_to_rownames(result_pvals_filter, "pair")
rownames(result_means_filter) <- NULL
result_means_filter <- column_to_rownames(result_means_filter, "pair")

ss <- result_pvals_filter < 0.01
result_pvals_filter[ss] <- '**'
s <- result_pvals_filter >= 0.01 & result_pvals_filter < 0.05
result_pvals_filter[s] <- '*'
ns <- is.na(result_pvals_filter)
result_pvals_filter[ns] <- 'n.s'

result_pvals_filter <- as.matrix(result_pvals_filter)
result_means_filter <- as.matrix(result_means_filter)

pdf("./6.0.cellchat_Epi_all/cellphoneDB/Cancer_Endo_diff/Cancer_Endo_heatmap_EOcoexp_mean.pdf", width = 8, height = 11)
bk <- c(seq(-2,2,by=0.05))
pheatmap(result_means_filter,
         scale = "row",
         cluster_row = T, 
         cluster_col = F, 
         border=NA,
         display_numbers = result_pvals_filter,
         fontsize_number = 8, 
         number_color = "black",
         cellwidth = 20, 
         cellheight =20,
         color=c(colorRampPalette(colors = c("white","#E64B35"))(length(bk))),
         legend_breaks=seq(-2,2,1),
         breaks=bk
)
dev.off()

# Cancer - Endo pairs heatmap (pval as heatmap, mean as tag)

pvals_Cancer_Endo <- pvals[, c(1:11, grep(".*spe.*Endothelial", colnames(pvals)))]
means_Cancer_Endo <- means[, c(1:11, grep(".*spe.*Endothelial", colnames(means)))]

result_EOPC_Endo <- plot_cpdb("EOPC_spe", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)
result_LOPC_Endo <- plot_cpdb("LOPC_spe", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)

result_EOPC_Endo_pvals <- data.frame(cell = result_EOPC_Endo$Var2, pair = result_EOPC_Endo$Var1, pvals = result_EOPC_Endo$pvals)
result_EOPC_Endo_pvals <- dcast(result_EOPC_Endo_pvals, pair ~ cell)
result_EOPC_Endo_pvals <- result_EOPC_Endo_pvals[, c(1:5)]

result_LOPC_Endo_pvals <- data.frame(cell = result_LOPC_Endo$Var2, pair = result_LOPC_Endo$Var1, pvals = result_LOPC_Endo$pvals)
result_LOPC_Endo_pvals <- dcast(result_LOPC_Endo_pvals, pair ~ cell)
result_LOPC_Endo_pvals <- result_LOPC_Endo_pvals[, c(1,6:11)]

dim(result_EOPC_Endo_pvals)
dim(result_LOPC_Endo_pvals)

result_pvals <- merge(result_EOPC_Endo_pvals, result_LOPC_Endo_pvals, by = "pair", all = T)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1")) & result_pvals$"EOPC1" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2")) & result_pvals$"EOPC2" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3")) & result_pvals$"EOPC3" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4")) & result_pvals$"EOPC4" < 0.05),]
write.table(result_pvals_filter, file="./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/Cancer_Endo_diff/sig.xls",sep="\t",row.names=F,quote=F)

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

pdf("./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/Cancer_Endo_diff/Cancer_Endo_heatmap_EOcoexp_pval.pdf", width = 6, height = 10)
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
         color=c(colorRampPalette(colors = c("white","salmon"))(50)),
         border_color="black"
)
dev.off()

#5.cellphoneDB to all annotated celltypes (cancer MP1 sep) -------------------------------------

#5.1.prepare data ------------------------------------------------

PCSC <- readRDS("./4.1.inferCNV/InferCancer/PCSC_inferred.Rds")
Epi <- readRDS("./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
Tcell <- readRDS("./5.2.Tcell_harmony/Tcell_annotated.Rds")
Mono <- readRDS("./5.4.Mono_harmony/Myeloid_annoted.Rds")
Fibr <- readRDS("./5.5.Fibr_harmony/Fibr_annoted.Rds")
Endo <- readRDS("./5.3.Endo_harmony/Endo_annoted.Rds")

Epi_all <- subset(PCSC, inferCNVCellStat != "Normal")
Epi_all <- subset(Epi_all, inferCNVCellStat != "Other_Epi")
Tcell <- subset(Tcell, celltype != "Epi_mix")
Tcell$celltype <- as.character(Tcell$celltype)
Mono <- subset(Mono, celltype != "Epi_mixed" & celltype != "Endo_mixed")
Mono$celltype <- as.character(Mono$celltype)
Fibr <- subset(Fibr, celltype != "Epi_mix" & celltype != "Endo_mix" & celltype != "Tcell_mix" & celltype != "Mye_mix")
Fibr$celltype <- as.character(Fibr$celltype)
Endo <- subset(Endo, celltype != "Epi-like" & celltype != "Myloid-like" & celltype != "T-cell-like")
Endo$celltype <- as.character(Endo$celltype)
Endo$celltype <- "Endothelial"

meta_epi <- data.frame(ID = rownames(Epi@meta.data), celltype = Epi@meta.data$MP1_OnOff)
meta_Tcell <- data.frame(ID = rownames(Tcell@meta.data), celltype = Tcell@meta.data$celltype)
meta_Mono <- data.frame(ID = rownames(Mono@meta.data), celltype = Mono@meta.data$celltype)
meta_Fibr <- data.frame(ID = rownames(Fibr@meta.data), celltype = Fibr@meta.data$celltype)
meta_Endo <- data.frame(ID = rownames(Endo@meta.data), celltype = Endo@meta.data$celltype)
meta_all <- data.frame(ID = rownames(Epi_all@meta.data), celltype = Epi_all@meta.data$celltype)

meta_all <- merge(meta_all, meta_epi, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

meta_all <- merge(meta_all, meta_Tcell, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

meta_all <- merge(meta_all, meta_Mono, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

meta_all <- merge(meta_all, meta_Fibr, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

meta_all$celltype[meta_all$celltype == "Endothelial"] <- "Endo"
meta_all <- merge(meta_all, meta_Endo, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

Epi_all@meta.data$celltype <- meta_all$celltype
Epi_all <- subset(Epi_all, celltype != "T-cell" & celltype != "Myeloid" & celltype != "Mesenchymal" & celltype != "Endo")

Epi_all@meta.data$celltype[Epi_all@meta.data$celltype == "On"] <- "MP1_On"
Epi_all@meta.data$celltype[Epi_all@meta.data$celltype == "Off"] <- "MP1_Off"
Epi_all@meta.data$celltype[Epi_all@meta.data$celltype == "M2-TAM1"] <- "M2-TAM-1"
Epi_all@meta.data$celltype[Epi_all@meta.data$celltype == "M2-TAM2"] <- "M2-TAM-2"

saveRDS(Epi_all, "6.0.cellchat_Epi_all/Epi_all_annotated.Rds")

## 5.2.EOPC ------------------------------------------------

EO <- subset(Epi_all, SampleType == "EOPC")

expr <- as.matrix(EO@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(EO@meta.data), celltype = EO@meta.data$celltype)  
write.table(celltype, paste0("./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

### 5.2.1.cellphoneDB analysis ----------------------------------------------

cellphonedb method statistical_analysis --counts-data gene_name --output-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output --threads 20 ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/input/celltype_all.txt ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/input/expr_all.txt
cellphonedb plot dot_plot --means-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output/means.txt --pvalues-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output/pvalues.txt --output-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output
cellphonedb plot heatmap_plot --pvalues-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output/pvalues.txt --output-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/input/celltype_all.txt


## 5.3.LOPC ------------------------------------------------

LO <- subset(Epi_all, SampleType == "LOPC")

expr <- as.matrix(LO@assays$RNA@counts)
expr <- data.frame(Gene = rownames(expr), expr)
write.table(expr, paste0("./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/input/expr_all.txt"),
            sep='\t', quote=F, row.names = F)  #表达谱
celltype <- data.frame(Cell = rownames(LO@meta.data), celltype = LO@meta.data$celltype)  
write.table(celltype, paste0("./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/input/celltype_all.txt"), 
            sep='\t', quote=F, row.names = F)  #celltype

### 5.3.1.cellphoneDB analysis ----------------------------------------------

cellphonedb method statistical_analysis --counts-data gene_name --output-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output --threads 20 ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/input/celltype_all.txt ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/input/expr_all.txt
cellphonedb plot dot_plot --means-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output/means.txt --pvalues-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output/pvalues.txt --output-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output
cellphonedb plot heatmap_plot --pvalues-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output/pvalues.txt --output-path ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output ./6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/input/celltype_all.txt

## 5.4.interpet and lotting ---------------------------------------------------

library(ktplots)

color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")
                    
pvals_EO <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output/pvalues.txt"), check.names = FALSE)
means_EO <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output/means.txt"), check.names = FALSE)
decon_EO <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output/deconvoluted.txt"), check.names = FALSE)
means_sig_EO <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/EO/output/significant_means.txt"), check.names = FALSE)

pvals_LO <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output/pvalues.txt"), check.names = FALSE)
means_LO <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output/means.txt"), check.names = FALSE)
decon_LO <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output/deconvoluted.txt"), check.names = FALSE)
means_sig_LO <- read.delim(paste0("6.0.cellchat_Epi_all/cellphoneDB_MP1_EOvsLO/LO/output/significant_means.txt"), check.names = FALSE)

# Cancer - Endo pairs heatmap (mean as heatmap, pval as tag)

## EO sig

means_sig_Cancer_Endo_EO <- means_sig_EO[, c(1:11, grep("MP1.*Endothelial", colnames(means_sig_EO)))]
means_sig_Cancer_Endo_EO <- means_sig_Cancer_Endo_EO[-which(is.na(means_sig_Cancer_Endo_EO$"MP1_Off|Endothelial") &
                                                             is.na(means_sig_Cancer_Endo_EO$"MP1_On|Endothelial")),]
colnames(means_sig_Cancer_Endo_EO)[c(12:13)] <- paste0("EOPC ", colnames(means_sig_Cancer_Endo_EO)[c(12:13)])

## LO sig

means_sig_Cancer_Endo_LO <- means_sig_LO[, c(1:11, grep("MP1.*Endothelial", colnames(means_sig_LO)))]
means_sig_Cancer_Endo_LO <- means_sig_Cancer_Endo_LO[-which(is.na(means_sig_Cancer_Endo_LO$"MP1_Off|Endothelial") &
                                                              is.na(means_sig_Cancer_Endo_LO$"MP1_On|Endothelial")),]
colnames(means_sig_Cancer_Endo_LO)[c(12:13)] <- paste0("LOPC ", colnames(means_sig_Cancer_Endo_LO)[c(12:13)])

## merge sig

name <- colnames(means_sig_Cancer_Endo_LO)[c(1:11)]
means_sig_Cancer_Endo <- full_join(means_sig_Cancer_Endo_EO, means_sig_Cancer_Endo_LO, by = name)

## EO mean/pval

means_Cancer_Endo_EO <- means_EO[, c(1:11, grep("MP1.*Endothelial", colnames(means_EO)))]
pvals_Cancer_Endo_EO <- pvals_EO[, c(1:11, grep("MP1.*Endothelial", colnames(pvals_EO)))]
means_Cancer_Endo_EO <- means_Cancer_Endo_EO[which(means_Cancer_Endo_EO$id_cp_interaction %in% means_sig_Cancer_Endo$id_cp_interaction),]
pvals_Cancer_Endo_EO <- pvals_Cancer_Endo_EO[which(pvals_Cancer_Endo_EO$id_cp_interaction %in% means_sig_Cancer_Endo$id_cp_interaction),]

colnames(means_Cancer_Endo_EO)[c(12:13)] <- paste0("EOPC - ", colnames(means_Cancer_Endo_EO)[c(12:13)])
colnames(pvals_Cancer_Endo_EO)[c(12:13)] <- paste0("EOPC - ", colnames(pvals_Cancer_Endo_EO)[c(12:13)])

## LO mean/pval

means_Cancer_Endo_LO <- means_LO[, c(1:11, grep("MP1.*Endothelial", colnames(means_LO)))]
pvals_Cancer_Endo_LO <- pvals_LO[, c(1:11, grep("MP1.*Endothelial", colnames(pvals_LO)))]
means_Cancer_Endo_LO <- means_Cancer_Endo_LO[which(means_Cancer_Endo_LO$id_cp_interaction %in% means_sig_Cancer_Endo$id_cp_interaction),]
pvals_Cancer_Endo_LO <- pvals_Cancer_Endo_LO[which(pvals_Cancer_Endo_LO$id_cp_interaction %in% means_sig_Cancer_Endo$id_cp_interaction),]

colnames(means_Cancer_Endo_LO)[c(12:13)] <- paste0("LOPC - ", colnames(means_Cancer_Endo_LO)[c(12:13)])
colnames(pvals_Cancer_Endo_LO)[c(12:13)] <- paste0("LOPC - ", colnames(pvals_Cancer_Endo_LO)[c(12:13)])

## Merge mean/pval

means_Cancer_Endo <- merge(means_Cancer_Endo_EO, means_Cancer_Endo_LO, by = name)
pvals_Cancer_Endo <- merge(pvals_Cancer_Endo_EO, pvals_Cancer_Endo_LO, by = name)

means_Cancer_Endo <- means_Cancer_Endo[,c(2,12:ncol(means_Cancer_Endo))]
pvals_Cancer_Endo <- pvals_Cancer_Endo[,c(2,12:ncol(pvals_Cancer_Endo))]

#plotting
rownames(means_Cancer_Endo) <- NULL
means_Cancer_Endo <- column_to_rownames(means_Cancer_Endo, "interacting_pair")
rownames(pvals_Cancer_Endo) <- NULL
pvals_Cancer_Endo <- column_to_rownames(pvals_Cancer_Endo, "interacting_pair")

ss <- pvals_Cancer_Endo < 0.01
pvals_Cancer_Endo[ss] <- '**'
s <- pvals_Cancer_Endo >= 0.01 & pvals_Cancer_Endo < 0.05
pvals_Cancer_Endo[s] <- '*'
ns <- pvals_Cancer_Endo >= 0.05
pvals_Cancer_Endo[ns] <- 'n.s'

## EO_On_spe

pvals_EO_On_spe <- pvals_Cancer_Endo[which(pvals_Cancer_Endo$"EOPC - MP1_On|Endothelial" != "n.s" &
                                             pvals_Cancer_Endo$"LOPC - MP1_On|Endothelial" == "n.s"),]
pvals_LO_On_spe <- pvals_Cancer_Endo[which(pvals_Cancer_Endo$"EOPC - MP1_On|Endothelial" == "n.s" &
                                             pvals_Cancer_Endo$"LOPC - MP1_On|Endothelial" != "n.s"),]
pvals_LO_On_share <- pvals_Cancer_Endo[which(pvals_Cancer_Endo$"EOPC - MP1_On|Endothelial" != "n.s" &
                                             pvals_Cancer_Endo$"LOPC - MP1_On|Endothelial" != "n.s"),]

#****************************************************************************************************************

means_Cancer_Endo <- as.matrix(means_Cancer_Endo)
pvals_Cancer_Endo <- as.matrix(pvals_Cancer_Endo)

pdf("./6.0.cellchat_Epi_all/cellphoneDB/Cancer_Endo_diff/Cancer_Endo_heatmap_EOcoexp_mean.pdf", width = 8, height = 11)
bk <- c(seq(-2,2,by=0.05))
pheatmap(means_Cancer_Endo,
         scale = "row",
         cluster_row = T, 
         cluster_col = F, 
         border=NA,
         display_numbers = pvals_Cancer_Endo,
         fontsize_number = 8, 
         number_color = "black",
         cellwidth = 20, 
         cellheight =20,
         color=c(colorRampPalette(colors = c("white","#E64B35"))(length(bk))),
         legend_breaks=seq(-2,2,1),
         breaks=bk
)
dev.off()

# Cancer - Endo pairs heatmap (pval as heatmap, mean as tag)

pvals_Cancer_Endo <- pvals[, c(1:11, grep(".*spe.*Endothelial", colnames(pvals)))]
means_Cancer_Endo <- means[, c(1:11, grep(".*spe.*Endothelial", colnames(means)))]

result_EOPC_Endo <- plot_cpdb("EOPC_spe", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)
result_LOPC_Endo <- plot_cpdb("LOPC_spe", "Endothelial", Epi_all, "inferCNVCellStat", means_Cancer_Endo, pvals_Cancer_Endo,
                              split.by = "orig.ident", cluster_rows = FALSE, scale = F,
                              keep_significant_only = T, return_table = T)

result_EOPC_Endo_pvals <- data.frame(cell = result_EOPC_Endo$Var2, pair = result_EOPC_Endo$Var1, pvals = result_EOPC_Endo$pvals)
result_EOPC_Endo_pvals <- dcast(result_EOPC_Endo_pvals, pair ~ cell)
result_EOPC_Endo_pvals <- result_EOPC_Endo_pvals[, c(1:5)]

result_LOPC_Endo_pvals <- data.frame(cell = result_LOPC_Endo$Var2, pair = result_LOPC_Endo$Var1, pvals = result_LOPC_Endo$pvals)
result_LOPC_Endo_pvals <- dcast(result_LOPC_Endo_pvals, pair ~ cell)
result_LOPC_Endo_pvals <- result_LOPC_Endo_pvals[, c(1,6:11)]

dim(result_EOPC_Endo_pvals)
dim(result_LOPC_Endo_pvals)

result_pvals <- merge(result_EOPC_Endo_pvals, result_LOPC_Endo_pvals, by = "pair", all = T)

#pvals
result_pvals_filter <- result_pvals[which((!is.na(result_pvals$"EOPC1")) & result_pvals$"EOPC1" < 0.05 &
                                            (!is.na(result_pvals$"EOPC2")) & result_pvals$"EOPC2" < 0.05 &
                                            (!is.na(result_pvals$"EOPC3")) & result_pvals$"EOPC3" < 0.05 &
                                            (!is.na(result_pvals$"EOPC4")) & result_pvals$"EOPC4" < 0.05),]
write.table(result_pvals_filter, file="./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/Cancer_Endo_diff/sig.xls",sep="\t",row.names=F,quote=F)

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

pdf("./6.0.cellchat_Epi_all/cellphoneDB_cancer_subtype/Cancer_Endo_diff/Cancer_Endo_heatmap_EOcoexp_pval.pdf", width = 6, height = 10)
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
         color=c(colorRampPalette(colors = c("white","salmon"))(50)),
         border_color="black"
)
dev.off()

# 6.NicheNet -------------------------------------------------------

Ir_network <- readRDS("./6.0.cellchat_Epi_all/NicheNet/lr_network.rds")
ligand_target_matrix <- readRDS("./6.0.cellchat_Epi_all/NicheNet/ligand_target_matrix.rds")
weighted_networks <- readRDS("./6.0.cellchat_Epi_all/NicheNet/weighted_networks.rds")

PCSC <- readRDS("./4.1.inferCNV/InferCancer/PCSC_inferred.Rds")
Epi <- readRDS("./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
Tcell <- readRDS("./5.2.Tcell_harmony/Nomix_Tcell_annotated.Rds")
Mono <- readRDS("./5.4.Mono_harmony/Nomix_Myeloid_annoted.Rds")
Fibr <- readRDS("./5.5.Fibr_harmony/Nomix_Fibr_annoted.Rds")
Endo <- readRDS("./5.3.Endo_harmony/Endo_annoted.Rds")

Epi_all <- subset(PCSC, inferCNVCellStat != "Normal" & inferCNVCellStat != "Other_Epi")

Endo <- subset(Endo, celltype != "Epi-like" & celltype != "Myloid-like" & celltype != "T-cell-like")
Endo$celltype <- as.character(Endo$celltype)
Endo$celltype <- "Endothelial"

meta_epi <- data.frame(ID = rownames(Epi@meta.data), celltype = as.character(Epi@meta.data$MP1_OnOff))
meta_Tcell <- data.frame(ID = rownames(Tcell@meta.data), celltype = as.character(Tcell@meta.data$celltype))
meta_Mono <- data.frame(ID = rownames(Mono@meta.data), celltype = as.character(Mono@meta.data$celltype))
meta_Fibr <- data.frame(ID = rownames(Fibr@meta.data), celltype = as.character(Fibr@meta.data$celltype))
meta_Endo <- data.frame(ID = rownames(Endo@meta.data), celltype = as.character(Endo@meta.data$celltype))
meta_all <- data.frame(ID = rownames(Epi_all@meta.data), celltype = as.character(Epi_all@meta.data$celltype))

meta_all <- merge(meta_all, meta_epi, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

meta_all <- merge(meta_all, meta_Tcell, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

meta_all <- merge(meta_all, meta_Mono, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

meta_all <- merge(meta_all, meta_Fibr, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

meta_all$celltype[meta_all$celltype == "Endothelial"] <- "Endo"
meta_all <- merge(meta_all, meta_Endo, by = "ID", all = T)
meta_all$celltype <- ifelse(!is.na(meta_all$celltype.y), meta_all$celltype.y, meta_all$celltype.x)
meta_all <- meta_all[, -c(2,3)]

cells <- Cells(Epi_all)
meta <- meta_all[ match(cells, meta_all$ID),]

Epi_all@meta.data$celltype <- meta$celltype
Epi_all <- subset(Epi_all, celltype != "T-cell" & celltype != "Myeloid" & celltype != "SMC" & celltype != "Fibroblast" & celltype != "Endo")

Epi_all@meta.data$celltype[Epi_all@meta.data$celltype == "On"] <- "MP1_On"
Epi_all@meta.data$celltype[Epi_all@meta.data$celltype == "Off"] <- "MP1_Off"

Epi_all$celltype <- factor(Epi_all$celltype, 
                           levels = c("MP1_On","MP1_Off",
                                      "Tcm","CD4_Tn","Treg","CD8_Tn","CD8_Tem","NK","NKT","B-cell",
                                      'pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
                                      'Mono_IL1B', 'Mono_FCN1',
                                      'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                      'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL',
                                      "Mast","Endothelial",
                                      "iCAF_APOD", "iCAF_CCL2",
                                      "myCAF_RGS5", "myCAF_CCL21", "SMC_KLF2", "SMC_HOPX", "Myoblast_MYF5"))

saveRDS(Epi_all, "6.0.cellchat_Epi_all/Epi_all_annotated_new.Rds")
Epi_all <- readRDS("6.0.cellchat_Epi_all/Epi_all_annotated_new.Rds")

## 6.1. predict ligands influence MP1 --------------------------------------

Epi_all$celltype2 <- as.character(Epi_all$celltype)
Epi_all$celltype2[Epi_all$celltype2 %in% c("MP1_On", "MP1_Off")] <- "Epithelial"
TME_cells <- unique(Epi_all$celltype2)[-grep("Epithelial", unique(Epi_all$celltype2))]
TME_cells <- TME_cells[-which(TME_cells == "Myoblast")]

Idents(Epi_all) <- "celltype2"

# overview **************************************

## predict nichenet (TME to Cancer) all ligands---------------------------------

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = Epi_all, 
  receiver = "Epithelial", 
  condition_colname = "SampleType", condition_oi = "EOPC", condition_reference = "LOPC", 
  sender = TME_cells, 
  geneset = "up",
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks,
  filter_top_ligands = F,
  top_n_targets = 2000)

## get MP1 genes and intersect with over-expressed targets -------------------------

MPGeneScoreAvg <- readRDS("./5.1.1.epithelial_NMF/process_k3/MP_Gene_Score_Avg.Rds")
geneset_oi <- MPGeneScoreAvg$MP1$Gene
matrix <- nichenet_output$ligand_target_matrix[,intersect(colnames(nichenet_output$ligand_target_matrix),geneset_oi)]

## get ordering of ligands -----------------------------------------

ligands_order <- nichenet_output$ligand_activities$test_ligand
rownames(matrix) <- gsub("\\.", "-", rownames(matrix))
matrix <- matrix[ligands_order,]
matrix_top <- matrix[c(1:30),]

rev_rownames <- rownames(matrix_top) %>% rev()
lig_target_heatmap <- matrix_top[rev_rownames,] %>%
  make_heatmap_ggplot("Ligands in TME cells",
                      "MP1 genes over-exp. in EOPC-Cancer cells",
                      color = "#7A378B",legend_position = "top", x_axis_position = "bottom",
                      legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "#7A378B", breaks = c(0,0.002,0.004)) +
  theme(axis.text.x = element_text(face = "italic", vjust = 0.5, size = 10, colour = "black"),
        axis.title.x = element_text(vjust = 1, size = 12),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12))
ggsave("./6.0.cellchat_Epi_all/NicheNet/lig_target_heatmap_all_to_Epi.pdf", height = 6, width = 5)

## predict nichenet (TME to Cancer) top ligands ---------------------------------

nichenet_output_top = nichenet_seuratobj_aggregate(
  seurat_obj = Epi_all, 
  receiver = "Epithelial", 
  condition_colname = "SampleType", condition_oi = "EOPC", condition_reference = "LOPC", 
  sender = TME_cells, 
  geneset = "up",
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks,
  top_n_targets = 250,
  top_n_ligands = 30)

### ligand_pearson heatmap ----------------------------------------

ligand_activities <- nichenet_output_top$ligand_activities
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>%
  as.matrix() %>% 
  magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>%
  make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>%
  make.names()
vis_ligand_pearson = ligand_pearson_matrix[c(1:30),] %>% 
  rev() %>% 
  as.matrix(ncol = 1) %>% 
  magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% 
  make_heatmap_ggplot("Prioritized ligands",
                      "Ligand activity", 
                      color = "#CD661D", 
                      legend_position = "top",
                      x_axis_position = "bottom", 
                      legend_title = "Pearson correlation coefficient\n(target gene prediction ability)") +
  theme(legend.text = element_text(size = 9),
        axis.text.x = element_text(angle = 0, size = 10, colour = "black"),
        axis.title.x = element_text(angle = 0, size = 12),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12))
p_ligand_pearson

ggsave("./6.0.cellchat_Epi_all/NicheNet/lig_pearson_heatmap_all_to_Epi.pdf", p_ligand_pearson, height = 6, width = 2)

### ligand_pearson heatmap ----------------------------------------

TME = Epi_all %>% subset(celltype != "MP1_On" & celltype != "MP1_Off")
TME$celltype <- factor(TME$celltype, 
                       levels = c("Tcm","CD4_Tn","Treg","CD8_Tn","CD8_Tem","NK","NKT","B-cell",
                                  'pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
                                  'Mono_IL1B', 'Mono_FCN1',
                                  'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                  'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL',
                                  "Mast","Endothelial",
                                  "CAF_APOD", "CAF_CCL2",
                                  "MyoFib_RGS5", "MyoFib_CCL21", "SMC_KLF2", "SMC_HOPX", "Myoblast_MYF5"))

Idents(TME) <- "celltype"

rotated_dotplot = DotPlot(TME, features = rev_rownames, cols = "RdYlBu", split.by = "SampleType") +
  coord_flip() + theme(legend.text = element_text(size = 10),
                       legend.title = element_text(size = 12),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                       axis.title.x = element_blank())

ggsave("./6.0.cellchat_Epi_all/NicheNet/lig_distribution_all_to_Epi.pdf", rotated_dotplot, height = 6, width = 12)

###  ligand-receptor network ---------------------------------------

exp <- as.matrix(subset(Epi_all, inferCNVCellStat == "Cancer")@assays$RNA@counts)
exp[exp > 0] <- 1
pos_pct <- rowSums(exp) / ncol(exp)
expressed_genes_receiver <- names(pos_pct)[which(pos_pct > 0.1)]

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% rev_rownames & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% rev_rownames & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

## Show a heatmap of the ligand-receptor interactions

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% .[,rev_rownames] %>% t() %>% 
  make_heatmap_ggplot("Ligands in TME cells",
                      "Receptors in Cancer cells", 
                      color = "#7A378B", 
                      x_axis_position = "bottom",
                      legend_title = "Prior interaction potential") +
  theme(axis.text.x = element_text(face = "italic", vjust = 0.5, size = 10, colour = "black"),
        axis.title.x = element_text(vjust = 1, size = 12),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12))
p_ligand_receptor_network

ggsave("6.0.cellchat_Epi_all/NicheNet/lig_lig_receptor_network.pdf", p_ligand_receptor_network, height = 6, width = 5)

### exp: TSLP PTGS2

pdf("6.0.cellchat_Epi_all/NicheNet/exp_TSLP_PTGS2_new.pdf", width = 6, height = 6)
VlnPlot(TME, features = c("TSLP", "PTGS2"), split.by = "SampleType", 
        cols = c("#053E7ACC","#A52A2ACC"), ncol = 1, pt.size = 0.01) +
  theme(axis.text.x = element_text(size = 8),
        axis.title.x = element_blank())
dev.off()

### Bulk data: TSLP/PTGS2 versus MP1 score -----------------------------

#### TCGA ----------------------------

expMtx <- read.table("./0.TCGA/mRNA_Matrix_TPM.txt", row.names = 1, header = T, sep = "\t")
expMtx <- expMtx[c("TSLP", "PTGS2"),]
expMtx_t <- as.data.frame(t(expMtx))

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_TCGA.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx <- AUCell_mtx["MP1",]
AUCell_mtx_t <- data.frame(MP1 = AUCell_mtx)

cor_mtx <- merge(expMtx_t, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
#cor_mtx <- cor_mtx[-grep("\\.11", rownames(cor_mtx)),]

cor.test(cor_mtx[,1], cor_mtx[,3], method = "spearman", exact = F) #sig
cor.test(cor_mtx[,2], cor_mtx[,3], method = "spearman", exact = F) #sig

#TSLP

cor <- cor.test(cor_mtx[,1], cor_mtx[,3], method = "spearman", exact = F)$estimate
pValue <- cor.test(cor_mtx[,1], cor_mtx[,3], method = "spearman", exact = F)$p.value

TSLP <- ggplot(cor_mtx, aes(x=TSLP, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = TSLP, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("TSLP Exp.") + ylab("MP1 AUCell")
TSLP <- ggMarginal(TSLP, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

#TSLP

PTGS2 <- ggplot(cor_mtx, aes(x=PTGS2, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = PTGS2, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("PTGS2 Exp.") + ylab("MP1 AUCell")
PTGS2 <- ggMarginal(PTGS2, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

pdf("6.0.cellchat_Epi_all/NicheNet/TCGA_cor_with_MP1_TSLP_PTGS2.pdf", width = 6, height = 3)
grid.arrange(grobs = list(TSLP, PTGS2), ncol = 2)
dev.off()

#### GSE157547 ----------------------------

expMtx <- read.table("./0.TCGA/GSE157547/MergeExpro_contrib1-GPL5175.txt", header = T, sep = "\t", check.names = F)
expMtx <- expMtx[which(expMtx$Gene %in% c("TSLP", "PTGS2")),]
rownames(expMtx) <- NULL
expMtx <- column_to_rownames(expMtx, "Gene")
expMtx_t <- as.data.frame(t(expMtx))

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE157547.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx <- AUCell_mtx["MP1",]
AUCell_mtx_t <- data.frame(MP1 = AUCell_mtx)

cor_mtx <- merge(expMtx_t, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
#cor_mtx <- cor_mtx[-grep("\\.11", rownames(cor_mtx)),]

cor.test(cor_mtx[,1], cor_mtx[,3], method = "spearman", exact = F) #sig
cor.test(cor_mtx[,2], cor_mtx[,3], method = "spearman", exact = F) #sig

#TSLP

TSLP <- ggplot(cor_mtx, aes(x=TSLP, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = TSLP, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("TSLP Exp.") + ylab("MP1 AUCell")
TSLP <- ggMarginal(TSLP, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

#TSLP

PTGS2 <- ggplot(cor_mtx, aes(x=PTGS2, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = PTGS2, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("PTGS2 Exp.") + ylab("MP1 AUCell")
PTGS2 <- ggMarginal(PTGS2, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

pdf("6.0.cellchat_Epi_all/NicheNet/GSE157547_cor_with_MP1_TSLP_PTGS2.pdf", width = 6, height = 3)
grid.arrange(grobs = list(TSLP, PTGS2), ncol = 2)
dev.off()

#### GSE141551 ----------------------------

exp <- read.table("./0.TCGA/GSE141551/MergeExpro_contrib1-GPL14951.txt", row.names = 1, header = T, sep = "\t")
exp <- as.matrix(exp)
expMtx <- exp[c("TSLP", "PTGS2"),]
expMtx_t <- as.data.frame(t(expMtx))

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE141551.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx <- AUCell_mtx["MP1",]
AUCell_mtx_t <- data.frame(MP1 = AUCell_mtx)

cor_mtx <- merge(expMtx_t, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
#cor_mtx <- cor_mtx[-grep("\\.11", rownames(cor_mtx)),]

cor.test(cor_mtx[,1], cor_mtx[,3], method = "spearman", exact = F) #sig
cor.test(cor_mtx[,2], cor_mtx[,3], method = "spearman", exact = F) #sig

#TSLP

TSLP <- ggplot(cor_mtx, aes(x=TSLP, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = TSLP, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("TSLP Exp.") + ylab("MP1 AUCell")
TSLP <- ggMarginal(TSLP, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

#TSLP

PTGS2 <- ggplot(cor_mtx, aes(x=PTGS2, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = PTGS2, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("PTGS2 Exp.") + ylab("MP1 AUCell")
PTGS2 <- ggMarginal(PTGS2, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

pdf("6.0.cellchat_Epi_all/NicheNet/GSE141551_cor_with_MP1_TSLP_PTGS2.pdf", width = 6, height = 3)
grid.arrange(grobs = list(TSLP, PTGS2), ncol = 2)
dev.off()

#### GSE141551 ----------------------------

exp <- read.table("./0.TCGA/GSE141551/MergeExpro_contrib1-GPL14951.txt", row.names = 1, header = T, sep = "\t")
exp <- as.matrix(exp)
expMtx <- exp[c("TSLP", "PTGS2"),]
expMtx_t <- as.data.frame(t(expMtx))

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE141551.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx <- AUCell_mtx["MP1",]
AUCell_mtx_t <- data.frame(MP1 = AUCell_mtx)

cor_mtx <- merge(expMtx_t, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
#cor_mtx <- cor_mtx[-grep("\\.11", rownames(cor_mtx)),]

cor.test(cor_mtx[,1], cor_mtx[,3], method = "spearman", exact = F) #sig
cor.test(cor_mtx[,2], cor_mtx[,3], method = "spearman", exact = F) #sig

#TSLP

TSLP <- ggplot(cor_mtx, aes(x=TSLP, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = TSLP, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("TSLP Exp.") + ylab("MP1 AUCell")
TSLP <- ggMarginal(TSLP, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

#TSLP

PTGS2 <- ggplot(cor_mtx, aes(x=PTGS2, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = PTGS2, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("PTGS2 Exp.") + ylab("MP1 AUCell")
PTGS2 <- ggMarginal(PTGS2, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

pdf("6.0.cellchat_Epi_all/NicheNet/GSE141551_cor_with_MP1_TSLP_PTGS2.pdf", width = 6, height = 3)
grid.arrange(grobs = list(TSLP, PTGS2), ncol = 2)
dev.off()

#### GSE21034 ----------------------------

exp <- read.table("./0.TCGA/GSE21034/MergeExpro_contrib1-GPL5188.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
expMtx <- exp[c("TSLP", "PTGS2"),]
expMtx_t <- as.data.frame(t(expMtx))

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE21034.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx <- AUCell_mtx["MP1",]
AUCell_mtx_t <- data.frame(MP1 = AUCell_mtx)

cor_mtx <- merge(expMtx_t, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
#cor_mtx <- cor_mtx[-grep("\\.11", rownames(cor_mtx)),]

cor.test(cor_mtx[,1], cor_mtx[,3], method = "spearman", exact = F) #sig
cor.test(cor_mtx[,2], cor_mtx[,3], method = "spearman", exact = F) #sig

#TSLP

TSLP <- ggplot(cor_mtx, aes(x=TSLP, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = TSLP, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("TSLP Exp.") + ylab("MP1 AUCell")
TSLP <- ggMarginal(TSLP, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

#PTGS2

PTGS2 <- ggplot(cor_mtx, aes(x=PTGS2, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = PTGS2, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + 
  xlab("PTGS2 Exp.") + ylab("MP1 AUCell")
PTGS2 <- ggMarginal(PTGS2, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

pdf("6.0.cellchat_Epi_all/NicheNet/GSE21034_cor_with_MP1_TSLP_PTGS2.pdf", width = 6, height = 3)
grid.arrange(grobs = list(TSLP, PTGS2), ncol = 2)
dev.off()

# ***********************************************

Idents(Epi_all) <- "celltype"

MPGeneScoreAvg <- readRDS("./5.1.1.epithelial_NMF/process_k3/MP_Gene_Score_Avg.Rds")
geneset_oi <- MPGeneScoreAvg$MP1$Gene

receiver_IDs <- colnames(subset(Epi_all, inferCNVCellStat == "Cancer"))
exp <- as.matrix(subset(Epi_all, inferCNVCellStat == "Cancer")@assays$RNA@counts)
exp[exp > 0] <- 1
pos_pct <- rowSums(exp) / ncol(exp)
expressed_genes_receiver <- names(pos_pct)[which(pos_pct > 0.1)]

## background genes --------------------------
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

TME_cells <- unique(Epi_all$celltype)[-grep("MP1", unique(Epi_all$celltype))]

for (cell in TME_cells) {
  
  dir <- paste0("./6.0.cellchat_Epi_all/NicheNet/Results_sep/", cell)
  dir.create(dir)
  
}

for (cell in TME_cells) {
  
  dir <- paste0("./6.0.cellchat_Epi_all/NicheNet/Results_sep/", cell)
  
  ## potential ligands ------------------------------
  
  sender_IDs <- colnames(subset(Epi_all, celltype == cell))
  exp <- as.matrix(subset(Epi_all, celltype == cell)@assays$RNA@counts)
  exp[exp > 0] <- 1
  pos_pct <- rowSums(exp) / ncol(exp)
  expressed_genes_sender <- names(pos_pct)[which(pos_pct > 0.1)]
  
  ligands = lr_network %>% pull(from) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
  potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
  
  ## predict ligand activity --------------------------------------
  
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                                background_expressed_genes = background_expressed_genes, 
                                                ligand_target_matrix = ligand_target_matrix, 
                                                potential_ligands = potential_ligands)
  
  saveRDS(ligand_activities, paste0(dir, "/ligand_activities.Rds"))
  
}
  
for (cell in TME_cells[21:26]) {
  
  dir <- paste0("./6.0.cellchat_Epi_all/NicheNet/Results_sep/", cell)
  
  ligand_activities <- readRDS(paste0(dir, "/ligand_activities.Rds"))
  
  ## top 20 ligands ---------------------------------------
  
  best_upstream_ligands = ligand_activities %>% top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
  head(best_upstream_ligands)
  
  p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
    geom_histogram(color="black", fill="#CD661D")  + 
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% 
                                    top_n(20, pearson) %>% 
                                    pull(pearson))), 
               color="#E46B4F", linetype="dashed", linewidth=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  p_hist_lig_activity
  
  ggsave(paste0(dir, "/hist_lig_activity.pdf"), p_hist_lig_activity)
  
  ## Infer target genes of top-ranked ligands and visualize in a heatmap ---------------------------------------
  
  active_ligand_target_links_df = best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
  nrow(active_ligand_target_links_df)
  head(active_ligand_target_links_df)
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                   ligand_target_matrix = ligand_target_matrix, 
                                                                   cutoff = 0.25)
  dim(active_ligand_target_links)
  head(active_ligand_target_links)
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets = active_ligand_target_links_df$target %>% unique()
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% 
    make_heatmap_ggplot(paste0("Prioritized ", cell, "-ligands"),
                        "MP1 genes in cancer cells", 
                        color = "#7A378B",legend_position = "top", x_axis_position = "top",
                        legend_title = "Regulatory potential") + 
    scale_fill_gradient2(low = "whitesmoke",  high = "#7A378B", breaks = c(0,0.005,0.01)) + 
    theme(axis.text.x = element_text(face = "italic"))
  p_ligand_target_network
  
  ggsave(paste0(dir, "/lig_target_network.pdf"), p_ligand_target_network)
  
}

# CAF-2 ----------------------------------------------------

ligand_activities <- readRDS("6.0.cellchat_Epi_all/NicheNet/Results_sep/CAF-2//ligand_activities.Rds")

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

## Ligand-receptor network inference for top-ranked ligands -----------------------------

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
lr_network_top_matrix = lr_network_top_matrix[,best_upstream_ligands]

## Show a heatmap of the ligand-receptor interactions ----------------------

p_ligand_receptor_network = lr_network_top_matrix %>% 
  .[,rev(colnames(.))] %>% t() %>% 
  make_heatmap_ggplot("Prioritized ligands (CAF-2)",
                      "Expressed receptors (Cancer cells)", 
                      color = "#A52A2A", 
                      x_axis_position = "top",
                      legend_title = "Prior interaction potential")
  
p_ligand_receptor_network

ggsave("6.0.cellchat_Epi_all/NicheNet/Results_sep/CAF-2/lig_receptor_network.pdf", p_ligand_receptor_network, height = 5, width = 8)

## Inferring ligand-to-target signaling paths -------------------------------

ligand_tf_matrix = readRDS("6.0.cellchat_Epi_all/NicheNet/ligand_tf_matrix.rds")
sig_network = readRDS("6.0.cellchat_Epi_all/NicheNet/signaling_network.rds")
gr_network = readRDS("6.0.cellchat_Epi_all/NicheNet/gr_network.rds")

ligands_all = c("TSLP", "PTGS2")
targets_all = c("B2M", "FKBP5", "HSP90AA1", "MAF", "NDRG1",
                "KLK3", "TMSB10", "MTDH", "SOD3", "NPY", "MSMB",
                "RPS20", "NDUFB9", "CCK", "SARAF", "DEGS1", "MANF")

#active_ligand_target_links_df = best_upstream_ligands %>% 
#  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
#targets_all = unique(active_ligand_target_links_df$target)

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

# To render the graph: uncomment following line of code
#DiagrammeR::render_graph(graph_min_max, layout = "tree")

#export to cytoscape
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
head(data_source_network) 

output_path = "6.0.cellchat_Epi_all/NicheNet/"
write_output = TRUE 

# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}

weighted_signaling_network <- bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory"))

# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
  data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}

# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

if(write_output){
  bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% write_tsv(paste0(output_path,"annotation_table.txt"))
}

annotation <- bind_rows(specific_annotation_tbl,non_specific_annotation_tbl)
annotation_split <- split(annotation$gene, annotation$annotation)

# sankeyNetwork -----------------------------------


ligand_to_receptor <- weighted_networks$lr_sig %>% 
  filter(from %in% annotation_split$ligand) %>%
  filter(to %in% annotation_split$receptor)
ligand_to_receptor[4,3] <- 0.8

receptor_to_TF <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(to %in% annotation_split$"transcriptional regulator") %>%
  select(!layer)

receptor_to_target <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(to %in% annotation_split$target) %>%
  select(!layer)

TF_to_target <- weighted_signaling_network %>% 
  filter(from %in% receptor_to_TF$to) %>% 
  filter(to %in% annotation_split$target) %>%
  select(!layer)

receptor_to_TF <- receptor_to_TF[-which(receptor_to_TF$to %in% setdiff(receptor_to_TF$to, TF_to_target$from)),]

links <- rbind(ligand_to_receptor,receptor_to_TF,receptor_to_target,TF_to_target)
nodes <- data.frame(gene = unique(c(links$from, links$to)))
order <- nodes$gene
nodes <- merge(nodes, annotation, by="gene")
nodes <- nodes[match(order, nodes$gene),]

Node2index = list()
Node2index[nodes$gene] = 0:(length(nodes$gene)-1)
index <- data.frame(gene = names(Node2index), index = unlist(Node2index))

order_from <- links$from
from2 <- merge(data.frame(gene = links$from), index, by = "gene")
from2 <- from2[match(order_from, from2$gene),]

order_to <- links$to
to2 <- merge(data.frame(gene = links$to), index, by = "gene")
to2 <- to2[match(order_to, to2$gene),]

links <- links %>%
  mutate(source = from2$index) %>%
  mutate(target = to2$index)

links$color <- "grey"
nodes$color <- "#E46B4F"
nodes$color[nodes$annotation == "ligand"] <- "#053E7A"
nodes$color[nodes$annotation == "receptor"] <- "#4A7088"
nodes$color[nodes$annotation == "target"] <- "#A52A2A"
  
write.table(links, "6.0.cellchat_Epi_all/NicheNet/Sankey_links.txt", sep = "\t", quote = F, row.names = F)
write.table(nodes, "6.0.cellchat_Epi_all/NicheNet/Sankey_nodes.txt", sep = "\t", quote = F, row.names = F)

links <- read.table("C:/Users/Administrator/Desktop/Sankey_links.txt",sep = "\t",header=T)
nodes <- read.table("C:/Users/Administrator/Desktop/Sankey_nodes.txt",sep = "\t",header=T)

links$color <- "grey"
nodes$color <- "#E46B4F"
nodes$color[nodes$annotation == "ligand"] <- "#053E7A"
nodes$color[nodes$annotation == "receptor"] <- "#008B8B"
nodes$color[nodes$annotation == "target"] <- "#A52A2A"
color2project = paste(c(unique(nodes$color),unique(links$color)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

library(networkD3)
sankeyNetwork(Links = links, Nodes = nodes, 
              # 指定source、target、value和name：
              Source = "source",
              Target = "target", 
              Value = "weight", 
              NodeID = "gene", # 节点的名字
              NodeGroup = "annotation",
              # 调整配置：
              fontSize = 24, # 节点的字体大小
              nodeWidth = 40, # 节点的宽度
              nodePadding = 10, # 节点之间的距离
              colourScale = JS(my_color),
              width = 1200,
              height = 960,
              sinksRight = F
)




MisLinks = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/link.csv",sep = ",")
MisNodes = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/node.csv",sep = ",")



color2project = paste(c(unique(MisNodes$group_color),unique(MisLinks$group_color2)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

sankeyNetwork(Links = MisLinks, 
              Nodes = MisNodes,
              Source = "source2", 
              Target = "target2",
              Value ="value",
              NodeID = "name", 
              LinkGroup  = "group_color2", 
              colourScale = JS(my_color),
              fontSize = 10
)

## 6.2. predict ligands turn off MP1 --------------------------------------

Epi_all$celltype2 <- as.character(Epi_all$celltype)
Epi_all$celltype2[Epi_all$celltype2 %in% c("MP1_On", "MP1_Off")] <- "Epithelial"
TME_cells <- unique(Epi_all$celltype2)[-grep("Epithelial", unique(Epi_all$celltype2))]
TME_cells <- TME_cells[-which(TME_cells == "Myoblast_MYF5")]

Idents(Epi_all) <- "celltype2"

# overview **************************************

## predict nichenet (TME to Cancer) all ligands---------------------------------

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = Epi_all, 
  receiver = "Epithelial", 
  condition_colname = "SampleType", condition_oi = "LOPC", condition_reference = "EOPC", 
  sender = TME_cells, 
  geneset = "up",
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks,
  filter_top_ligands = F,
  top_n_targets = 500)

saveRDS(nichenet_output, "./6.0.cellchat_Epi_all/NicheNet_Inv/nichenet_Overview_LO_vs_EO.Rds")

## get MP1_Off over-exp genes and intersect with over-expressed targets -------------------------

DefaultAssay(Epi) <- "RNA"
Epi <- Epi[-grep("^RP[SL]", rownames(Epi)),]
Epi <- Epi[-grep("^MT-", rownames(Epi)),]

Idents(Epi) <- "MP1_HL"

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)
HL_markers <- FindMarkers(Epi,
                          ident.1 = "Low",
                          ident.2 = "High",
                          min.pct = 0.3, 
                          logfc.threshold = 0,
                          verbose = T,
                          only.pos = T
)
saveRDS(HL_markers, "./6.0.cellchat_Epi_all/NicheNet_Inv/geneset_oi/Markers_MP1_Low_vs_High_minPct0.3.Rds")
top_HL_markers <- HL_markers[HL_markers$avg_log2FC > 0.25, ]
#geneset_oi <- rownames(top_HL_markers)


Idents(Epi) <- "MP1_OnOff"

plan("multicore", workers = 20)  #multicore
options(future.globals.maxSize = 2000 * 1024^2)
OnOff_markers <- FindMarkers(Epi,
                          ident.1 = "Off",
                          ident.2 = "On",
                          min.pct = 0.3, 
                          logfc.threshold = 0,
                          verbose = T,
                          only.pos = T
)
saveRDS(OnOff_markers, "./6.0.cellchat_Epi_all/NicheNet_Inv/geneset_oi/Markers_MP1_Off_vs_On_minPct0.3.Rds")
OnOff_markers <- readRDS("./6.0.cellchat_Epi_all/NicheNet_Inv/geneset_oi/Markers_MP1_Off_vs_On_minPct0.3.Rds")
top_OnOff_markers <- OnOff_markers[OnOff_markers$avg_log2FC > 0.25, ]
geneset_oi <- rownames(top_OnOff_markers)

matrix <- nichenet_output$ligand_target_matrix[,intersect(colnames(nichenet_output$ligand_target_matrix),geneset_oi)]

## get ordering of ligands -----------------------------------------

ligands_order <- nichenet_output$ligand_activities$test_ligand
matrix <- matrix[match(ligands_order, rownames(matrix)),]
matrix <- matrix[-which(is.na(rownames(matrix))),]
matrix_top <- matrix[c(1:30),]

rev_rownames <- rownames(matrix_top) %>% rev()
lig_target_heatmap <- matrix_top[rev_rownames,] %>%
  make_heatmap_ggplot("Ligands in TME cells",
                      "Over-exp. genes in LOPC-MP1_Off cells",
                      color = "#7A378B",legend_position = "top", x_axis_position = "bottom",
                      legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "#7A378B") +
  theme(axis.text.x = element_text(face = "italic", vjust = 0.5, size = 10, colour = "black"),
        axis.title.x = element_text(vjust = 1, size = 12),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12))
ggsave("./6.0.cellchat_Epi_all/NicheNet_Inv/lig_target_heatmap_all_to_Epi.pdf", height = 6, width = 5)

## predict nichenet (TME to Cancer) top ligands ---------------------------------

nichenet_output_top = nichenet_seuratobj_aggregate(
  seurat_obj = Epi_all, 
  receiver = "Epithelial", 
  condition_colname = "SampleType", condition_oi = "LOPC", condition_reference = "EOPC", 
  sender = TME_cells, 
  geneset = "up",
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks,
  top_n_targets = 250,
  top_n_ligands = 30)

### ligand_pearson heatmap ----------------------------------------

ligand_activities <- nichenet_output_top$ligand_activities
ligand_activities <- ligand_activities[ligand_activities$test_ligand %in% rownames(matrix),]
ligand_pearson_matrix <- ligand_activities %>% select(pearson) %>%
  as.matrix() %>% 
  magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>%
  make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>%
  make.names()
vis_ligand_pearson = ligand_pearson_matrix[c(1:30),] %>% 
  rev() %>% 
  as.matrix(ncol = 1) %>% 
  magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% 
  make_heatmap_ggplot("Prioritized ligands",
                      "Ligand activity", 
                      color = "#CD661D", 
                      legend_position = "top",
                      x_axis_position = "bottom", 
                      legend_title = "Pearson correlation coefficient\n(target gene prediction ability)") +
  theme(legend.text = element_text(size = 9),
        axis.text.x = element_text(angle = 0, size = 10, colour = "black"),
        axis.title.x = element_text(angle = 0, size = 12),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12))
p_ligand_pearson

ggsave("./6.0.cellchat_Epi_all/NicheNet_Inv/lig_pearson_heatmap_all_to_Epi.pdf", p_ligand_pearson, height = 6, width = 2)

### ligand_pearson heatmap ----------------------------------------

TME <- Epi_all %>% subset(celltype != "MP1_On" & celltype != "MP1_Off")
TME$celltype <- factor(TME$celltype, 
                       levels = c("Tcm","CD4_Tn","Treg","CD8_Tn","CD8_Tem","NK","NKT","B-cell",
                                  'pDC_JCHAIN', 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',
                                  'Mono_IL1B', 'Mono_FCN1',
                                  'Macro_CXCL11', 'Macro_C3', 'Macro_RGS1', 
                                  'Macro_FOLR2', 'Macro_APOE', 'Macro_MT1A', 'Macro_FTL',
                                  "Mast","Endothelial",
                                  "CAF_APOD", "CAF_CCL2",
                                  "MyoFib_RGS5", "MyoFib_CCL21", "SMC_KLF2", "SMC_HOPX", "Myoblast_MYF5"))

Idents(TME) <- "celltype"

rotated_dotplot = DotPlot(TME, features = rev_rownames, cols = "RdYlBu", split.by = "SampleType") +
  coord_flip() + theme(legend.text = element_text(size = 8),
                       legend.title = element_text(size = 10),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
                       axis.title.x = element_blank(),
                       axis.text.y = element_text(size = 10),
                       axis.title.y = element_blank())

ggsave("./6.0.cellchat_Epi_all/NicheNet_Inv/lig_distribution_allTME_to_Epi.pdf", rotated_dotplot, height = 6, width = 12)

CAF <- Epi_all %>% subset(celltype %in% c("CAF_APOD","CAF_CCL2"))
CAF$celltype <- factor(CAF$celltype, 
                       levels = c("CAF_APOD", "CAF_CCL2"))
Idents(CAF) <- "celltype"

rotated_dotplot = DotPlot(CAF, features = rev_rownames, cols = "RdYlBu", split.by = "SampleType") +
  coord_flip() + theme(legend.text = element_text(size = 8),
                       legend.title = element_text(size = 10),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
                       axis.title.x = element_blank(),
                       axis.text.y = element_text(size = 10),
                       axis.title.y = element_blank())

ggsave("./6.0.cellchat_Epi_all/NicheNet_Inv/lig_distribution_CAF_to_Epi.pdf", rotated_dotplot, height = 6, width = 4)

BMPs_rotated_dotplot <- DotPlot(TME, features = c("BMP5","BMP7","BMP4"), dot.scale = 4, 
                                cols = "RdYlBu", split.by = "SampleType") +
  coord_flip() + theme(legend.text = element_text(size = 8),
                       legend.title = element_text(size = 10),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
                       axis.title.x = element_blank(),
                       axis.text.y = element_text(size = 10),
                       axis.title.y = element_blank())

ggsave("./6.0.cellchat_Epi_all/NicheNet_Inv/BMPs_distribution_CAF_to_Epi.pdf", BMPs_rotated_dotplot, height = 3, width = 8)

###  ligand-receptor network ---------------------------------------

exp <- as.matrix(subset(Epi_all, inferCNVCellStat == "Cancer")@assays$RNA@counts)
exp[exp > 0] <- 1
pos_pct <- rowSums(exp) / ncol(exp)
expressed_genes_receiver <- names(pos_pct)[which(pos_pct > 0.3)]
expressed_receptors <- expressed_genes_receiver

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% rev_rownames & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% rev_rownames & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

## Show a heatmap of the ligand-receptor interactions

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
vis_ligand_receptor_network <- vis_ligand_receptor_network[,match(rev_rownames, colnames(vis_ligand_receptor_network))]
vis_ligand_receptor_network <- vis_ligand_receptor_network[,-which(is.na(colnames(vis_ligand_receptor_network)))]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% 
  make_heatmap_ggplot("Ligands in TME cells",
                      "Receptors in Cancer cells", 
                      color = "#7A378B", 
                      x_axis_position = "bottom",
                      legend_title = "Prior interaction potential") +
  theme(axis.text.x = element_text(face = "italic", vjust = 0.5, size = 10, colour = "black"),
        axis.title.x = element_text(vjust = 1, size = 12),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12))
p_ligand_receptor_network

ggsave("6.0.cellchat_Epi_all/NicheNet_Inv/lig_lig_receptor_network.pdf", p_ligand_receptor_network, height = 6, width = 5)

## Bulk data: BMPs versus MP1 score -----------------------------

#### TCGA ----------------------------

expMtx <- read.table("./0.TCGA/mRNA_Matrix_TPM.txt", row.names = 1, header = T, sep = "\t")
expMtx <- expMtx[c("BMP5", "BMP4", "BMP7"),]
expMtx_t <- as.data.frame(t(expMtx))

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_TCGA.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx <- AUCell_mtx["MP1",]
AUCell_mtx_t <- data.frame(MP1 = AUCell_mtx)

cor_mtx <- merge(expMtx_t, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
#cor_mtx <- cor_mtx[-grep("\\.11", rownames(cor_mtx)),]

cor.test(cor_mtx[,1], cor_mtx[,4], method = "pearson", exact = F) #sig
cor.test(cor_mtx[,2], cor_mtx[,4], method = "pearson", exact = F) #sig
cor.test(cor_mtx[,3], cor_mtx[,4], method = "pearson", exact = F) #sig

#BMP5

BMP5 <- ggplot(cor_mtx, aes(x=BMP5, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = BMP5, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP5 Exp.") + ylab("MP1 AUCell")

#BMP4

BMP4 <- ggplot(cor_mtx, aes(x=BMP4, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = BMP4, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP4 Exp.") + ylab("MP1 AUCell")

#BMP7

BMP7 <- ggplot(cor_mtx, aes(x=BMP7, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = BMP7, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP7 Exp.") + ylab("MP1 AUCell")

pdf("6.0.cellchat_Epi_all/NicheNet_Inv/TCGA_cor_MP1_with_BMPs.pdf", width = 3, height = 9)
grid.arrange(grobs = list(BMP5, BMP4, BMP7), ncol = 1)
dev.off()

#### GSE157547 ----------------------------

expMtx <- read.table("./0.TCGA/GSE157547/MergeExpro_contrib1-GPL5175.txt", header = T, sep = "\t", check.names = F)
expMtx <- expMtx[which(expMtx$Gene %in% c("BMP5", "BMP4", "BMP7")),]
rownames(expMtx) <- NULL
expMtx <- column_to_rownames(expMtx, "Gene")
expMtx_t <- as.data.frame(t(expMtx))

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE157547.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx <- AUCell_mtx["MP1",]
AUCell_mtx_t <- data.frame(MP1 = AUCell_mtx)

cor_mtx <- merge(expMtx_t, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
#cor_mtx <- cor_mtx[-grep("\\.11", rownames(cor_mtx)),]

cor.test(cor_mtx[,1], cor_mtx[,4], method = "pearson", exact = F) #sig
cor.test(cor_mtx[,2], cor_mtx[,4], method = "pearson", exact = F) #sig
cor.test(cor_mtx[,3], cor_mtx[,4], method = "pearson", exact = F) #sig

#BMP5

BMP5 <- ggplot(cor_mtx, aes(x=BMP5, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = BMP5, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP5 Exp.") + ylab("MP1 AUCell")

#BMP4

BMP4 <- ggplot(cor_mtx, aes(x=BMP4, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = BMP4, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP4 Exp.") + ylab("MP1 AUCell")

#BMP7

BMP7 <- ggplot(cor_mtx, aes(x=BMP7, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = BMP7, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP7 Exp.") + ylab("MP1 AUCell")

pdf("6.0.cellchat_Epi_all/NicheNet_Inv/GSE157547_cor_MP1_with_BMPs.pdf", width = 3, height = 9)
grid.arrange(grobs = list(BMP5, BMP4, BMP7), ncol = 1)
dev.off()

#### GSE141551 ----------------------------

exp <- read.table("./0.TCGA/GSE141551/MergeExpro_contrib1-GPL14951.txt", row.names = 1, header = T, sep = "\t")
exp <- as.matrix(exp)
expMtx <- exp[c("BMP5", "BMP4", "BMP7"),]
expMtx_t <- as.data.frame(t(expMtx))

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE141551.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx <- AUCell_mtx["MP1",]
AUCell_mtx_t <- data.frame(MP1 = AUCell_mtx)

cor_mtx <- merge(expMtx_t, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
#cor_mtx <- cor_mtx[-grep("\\.11", rownames(cor_mtx)),]

cor.test(cor_mtx[,1], cor_mtx[,4], method = "spearman", exact = F) #nosig
cor.test(cor_mtx[,2], cor_mtx[,4], method = "spearman", exact = F) #sig
cor.test(cor_mtx[,3], cor_mtx[,4], method = "spearman", exact = F) #sig

#BMP5

BMP5 <- ggplot(cor_mtx, aes(x=BMP5, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = BMP5, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP5 Exp.") + ylab("MP1 AUCell")

#BMP4

BMP4 <- ggplot(cor_mtx, aes(x=BMP4, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = BMP4, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP4 Exp.") + ylab("MP1 AUCell")

#BMP7

BMP7 <- ggplot(cor_mtx, aes(x=BMP7, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = BMP7, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP7 Exp.") + ylab("MP1 AUCell")

pdf("6.0.cellchat_Epi_all/NicheNet_Inv/GSE141551_cor_MP1_with_BMPs.pdf", width = 3, height = 9)
grid.arrange(grobs = list(BMP5, BMP4, BMP7), ncol = 1)
dev.off()

#### GSE21034 ----------------------------

exp <- read.table("./0.TCGA/GSE21034/MergeExpro_contrib1-GPL5188.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
expMtx <- exp[c("BMP5", "BMP4", "BMP7"),]
expMtx_t <- as.data.frame(t(expMtx))

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE21034.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx <- AUCell_mtx["MP1",]
AUCell_mtx_t <- data.frame(MP1 = AUCell_mtx)

cor_mtx <- merge(expMtx_t, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
#cor_mtx <- cor_mtx[-grep("\\.11", rownames(cor_mtx)),]

cor.test(cor_mtx[,1], cor_mtx[,4], method = "spearman", exact = F) #sig
cor.test(cor_mtx[,2], cor_mtx[,4], method = "spearman", exact = F) #sig
cor.test(cor_mtx[,3], cor_mtx[,4], method = "spearman", exact = F) #sig

#BMP5

BMP5 <- ggplot(cor_mtx, aes(x=BMP5, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = BMP5, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP5 Exp.") + ylab("MP1 AUCell")

#BMP4

BMP4 <- ggplot(cor_mtx, aes(x=BMP4, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = BMP4, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP4 Exp.") + ylab("MP1 AUCell")

#BMP7

BMP7 <- ggplot(cor_mtx, aes(x=BMP7, y=MP1)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = BMP7, y = MP1)) +
  theme(axis.text.x=element_text(colour="black", size = 14),
        axis.text.y=element_text(size=14, colour="black"), 
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid")) + 
  xlab("BMP7 Exp.") + ylab("MP1 AUCell")

pdf("6.0.cellchat_Epi_all/NicheNet_Inv/GSE21034_cor_MP1_with_BMPs.pdf", width = 3, height = 9)
grid.arrange(grobs = list(BMP5, BMP4, BMP7), ncol = 1)
dev.off()

## top ligands and receptors in Cancers ------------------------------------------

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = ligands[match(rev_rownames, ligands)]

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

best_upstream_ligands = potential_ligands
head(best_upstream_ligands)

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% expressed_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
lr_network_top_matrix = lr_network_top_matrix[,match(rev(expressed_ligands), colnames(lr_network_top_matrix))]

# transfer NA ligands to 0
NAs <- rev(expressed_ligands)[match(colnames(lr_network_top_matrix), rev(expressed_ligands)) %>% is.na]
index <- which(is.na(colnames(lr_network_top_matrix)))
colnames(lr_network_top_matrix)[index] <- NAs
na <- is.na(lr_network_top_matrix)
lr_network_top_matrix[na] <- 0

## Show a heatmap of the ligand-receptor interactions ----------------------

p_ligand_receptor_network = lr_network_top_matrix %>% 
  .[,rev(colnames(.))] %>% t() %>% 
  make_heatmap_ggplot("Prioritized ligands (CAF)",
                      "Expressed receptors (Cancer cells)", 
                      color = "#A52A2A", 
                      x_axis_position = "top",
                      legend_title = "Prior interaction potential")

p_ligand_receptor_network

ggsave("6.0.cellchat_Epi_all/NicheNet_Inv/lig_receptor_network.pdf", p_ligand_receptor_network, height = 6, width = 5)

## Inferring ligand-to-target signaling paths -------------------------------

ligand_tf_matrix = readRDS("6.0.cellchat_Epi_all/NicheNet/ligand_tf_matrix.rds")
sig_network = readRDS("6.0.cellchat_Epi_all/NicheNet/signaling_network.rds")
gr_network = readRDS("6.0.cellchat_Epi_all/NicheNet/gr_network.rds")

ligands_all = c("BMP4", "BMP5", "BMP7")
targets_all = c("CLDN4", "CLDN7", "ELF3", "ENC1", "RACK1",
                "TM4SF1", "CKS1B", "ACSL5", "EIF5A", "SOX9", "GATA2",
                "MIF", "CKS2", "STMN1", "ZNF503", "KRT19", "H2AFZ",
                "CXADR","PRSS23","RNF43","TNFSF10","TPI1","NPM1",
                "ACTG1","ISG15","EPCAM","TXN","NACA")

#active_ligand_target_links_df = best_upstream_ligands %>% 
#  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
#targets_all = unique(active_ligand_target_links_df$target)

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

# To render the graph: uncomment following line of code
#DiagrammeR::render_graph(graph_min_max, layout = "tree")

#export to cytoscape
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
head(data_source_network) 

output_path = "6.0.cellchat_Epi_all/NicheNet_Inv/"
write_output = TRUE 

# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}

weighted_signaling_network <- bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory"))

# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
  data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}

# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

if(write_output){
  bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% write_tsv(paste0(output_path,"annotation_table.txt"))
}

annotation <- bind_rows(specific_annotation_tbl,non_specific_annotation_tbl)
annotation_split <- split(annotation$gene, annotation$annotation)

## sankeyNetwork -----------------------------------


ligand_to_receptor <- weighted_networks$lr_sig %>% 
  filter(from %in% annotation_split$ligand) %>%
  filter(to %in% annotation_split$receptor) %>%
  filter(to %in% c("BMPR1B","BMPR2"))

receptor_to_TF <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(from %in% c("BMPR1B","BMPR2")) %>% 
  filter(to %in% annotation_split$"transcriptional regulator") %>%
  select(!layer)

receptor_to_target <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(from %in% c("BMPR1B","BMPR2")) %>% 
  filter(to %in% annotation_split$target) %>%
  select(!layer)

TF_to_target <- weighted_signaling_network %>% 
  filter(from %in% receptor_to_TF$to) %>%
  filter(to %in% annotation_split$target) %>%
  select(!layer)

receptor_to_TF <- receptor_to_TF[-which(receptor_to_TF$to %in% setdiff(receptor_to_TF$to, TF_to_target$from)),]

links <- rbind(ligand_to_receptor,receptor_to_TF,receptor_to_target,TF_to_target)
nodes <- data.frame(gene = unique(c(links$from, links$to)))
order <- nodes$gene
nodes <- merge(nodes, annotation, by="gene")
nodes <- nodes[match(order, nodes$gene),]

Node2index = list()
Node2index[nodes$gene] = 0:(length(nodes$gene)-1)
index <- data.frame(gene = names(Node2index), index = unlist(Node2index))

order_from <- links$from
from2 <- merge(data.frame(gene = links$from), index, by = "gene")
from2 <- from2[match(order_from, from2$gene),]

order_to <- links$to
to2 <- merge(data.frame(gene = links$to), index, by = "gene")
to2 <- to2[match(order_to, to2$gene),]

links <- links %>%
  mutate(source = from2$index) %>%
  mutate(target = to2$index)

write.table(links, "6.0.cellchat_Epi_all/NicheNet_Inv/Sankey_links.txt", sep = "\t", quote = F, row.names = F)
write.table(nodes, "6.0.cellchat_Epi_all/NicheNet_Inv/Sankey_nodes.txt", sep = "\t", quote = F, row.names = F)

links <- read.table("C:/Users/Administrator/Desktop/Sankey_links.txt",sep = "\t",header=T)
nodes <- read.table("C:/Users/Administrator/Desktop/Sankey_nodes.txt",sep = "\t",header=T)

links$color <- "grey"
nodes$color <- "#E46B4F"
nodes$color[nodes$annotation == "ligand"] <- "#053E7A"
nodes$color[nodes$annotation == "receptor"] <- "#008B8B"
nodes$color[nodes$annotation == "target"] <- "#A52A2A"
color2project = paste(c(unique(nodes$color),unique(links$color)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

library(networkD3)
sankeyNetwork(Links = links, Nodes = nodes, 
              # 指定source、target、value和name：
              Source = "source",
              Target = "target", 
              Value = "weight", 
              NodeID = "gene", # 节点的名字
              NodeGroup = "annotation",
              # 调整配置：
              fontSize = 24, # 节点的字体大小
              nodeWidth = 40, # 节点的宽度
              nodePadding = 10, # 节点之间的距离
              colourScale = JS(my_color),
              width = 1200,
              height = 960,
              sinksRight = F
)




MisLinks = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/link.csv",sep = ",")
MisNodes = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/node.csv",sep = ",")



color2project = paste(c(unique(MisNodes$group_color),unique(MisLinks$group_color2)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

sankeyNetwork(Links = MisLinks, 
              Nodes = MisNodes,
              Source = "source2", 
              Target = "target2",
              Value ="value",
              NodeID = "name", 
              LinkGroup  = "group_color2", 
              colourScale = JS(my_color),
              fontSize = 10
)


# ***********************************************

Idents(Epi_all) <- "celltype"

receiver_IDs <- colnames(subset(Epi_all, inferCNVCellStat == "Cancer"))
exp <- as.matrix(subset(Epi_all, inferCNVCellStat == "Cancer")@assays$RNA@counts)
exp[exp > 0] <- 1
pos_pct <- rowSums(exp) / ncol(exp)
expressed_genes_receiver <- names(pos_pct)[which(pos_pct > 0.3)]

## background genes --------------------------
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

TME_cells <- unique(Epi_all$celltype)[-grep("MP1", unique(Epi_all$celltype))]

for (cell in TME_cells) {
  
  dir <- paste0("./6.0.cellchat_Epi_all/NicheNet_Inv/Results_sep/", cell)
  dir.create(dir)
  
}

for (cell in TME_cells) {
  
  dir <- paste0("./6.0.cellchat_Epi_all/NicheNet_Inv/Results_sep/", cell)
  
  ## potential ligands ------------------------------
  
  sender_IDs <- colnames(subset(Epi_all, celltype == cell))
  exp <- as.matrix(subset(Epi_all, celltype == cell)@assays$RNA@counts)
  exp[exp > 0] <- 1
  pos_pct <- rowSums(exp) / ncol(exp)
  expressed_genes_sender <- names(pos_pct)[which(pos_pct > 0.3)]
  
  ligands = lr_network %>% pull(from) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
  potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
  
  ## predict ligand activity --------------------------------------
  
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                                background_expressed_genes = background_expressed_genes, 
                                                ligand_target_matrix = ligand_target_matrix, 
                                                potential_ligands = potential_ligands)
  
  saveRDS(ligand_activities, paste0(dir, "/ligand_activities.Rds"))
  
}

for (cell in TME_cells) {
  
  dir <- paste0("./6.0.cellchat_Epi_all/NicheNet_Inv/Results_sep/", cell)
  
  ligand_activities <- readRDS(paste0(dir, "/ligand_activities.Rds"))
  
  ## top 20 ligands ---------------------------------------
  
  best_upstream_ligands = ligand_activities %>% top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
  head(best_upstream_ligands)
  
  p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
    geom_histogram(color="black", fill="#CD661D")  + 
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% 
                                    top_n(20, pearson) %>% 
                                    pull(pearson))), 
               color="#E46B4F", linetype="dashed", linewidth=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  p_hist_lig_activity
  
  ggsave(paste0(dir, "/hist_lig_activity.pdf"), p_hist_lig_activity)
  
  ## Infer target genes of top-ranked ligands and visualize in a heatmap ---------------------------------------
  
  active_ligand_target_links_df = best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
  nrow(active_ligand_target_links_df)
  head(active_ligand_target_links_df)
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                   ligand_target_matrix = ligand_target_matrix, 
                                                                   cutoff = 0.25)
  dim(active_ligand_target_links)
  head(active_ligand_target_links)
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets = active_ligand_target_links_df$target %>% unique()
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% 
    make_heatmap_ggplot(paste0("Prioritized ", cell, "-ligands"),
                        "Over-exp genes in MP1_Off cancer cells", 
                        color = "#7A378B",legend_position = "top", x_axis_position = "top",
                        legend_title = "Regulatory potential") + 
    scale_fill_gradient2(low = "whitesmoke",  high = "#7A378B", breaks = c(0,0.005,0.01)) + 
    theme(axis.text.x = element_text(face = "italic"))
  p_ligand_target_network
  
  ggsave(paste0(dir, "/lig_target_network.pdf"), p_ligand_target_network)
  
}

# CAF ------------------------------------------------
cell <- c("CAF_APOD", "CAF_CCL2")

dir <- paste0("./6.0.cellchat_Epi_all/NicheNet_Inv/Results_sep/", "CAF")
dir.create(dir)

## potential ligands ------------------------------

sender_IDs <- colnames(subset(Epi_all, celltype %in% cell))
exp <- as.matrix(subset(Epi_all, celltype %in% cell)@assays$RNA@counts)
exp[exp > 0] <- 1
pos_pct <- rowSums(exp) / ncol(exp)
expressed_genes_sender <- names(pos_pct)[which(pos_pct > 0.3)]

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

## predict ligand activity --------------------------------------

ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

saveRDS(ligand_activities, paste0(dir, "/ligand_activities.Rds"))

## top 20 ligands ---------------------------------------

ligand_activities <- readRDS(paste0(dir, "/ligand_activities.Rds"))

best_upstream_ligands = ligand_activities %>% top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="#CD661D")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% 
                                  top_n(20, pearson) %>% 
                                  pull(pearson))), 
             color="#E46B4F", linetype="dashed", linewidth=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

ggsave(paste0(dir, "/hist_lig_activity.pdf"), p_hist_lig_activity)

## Infer target genes of top-ranked ligands and visualize in a heatmap ---------------------------------------

active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, 
                                                                 cutoff = 0.25)
dim(active_ligand_target_links)
head(active_ligand_target_links)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot(paste0("Prioritized ", "CAF", "-ligands"),
                      "Over-exp genes in MP1_Off cancer cells", 
                      color = "#7A378B",legend_position = "top", x_axis_position = "top",
                      legend_title = "Regulatory potential") + 
  scale_fill_gradient2(low = "whitesmoke",  high = "#7A378B", breaks = c(0,0.005,0.01)) + 
  theme(axis.text.x = element_text(face = "italic"))
p_ligand_target_network

ggsave(paste0(dir, "/lig_target_network.pdf"), p_ligand_target_network)

## Inferring ligand-to-target signaling paths -------------------------------

ligand_tf_matrix = readRDS("6.0.cellchat_Epi_all/NicheNet/ligand_tf_matrix.rds")
sig_network = readRDS("6.0.cellchat_Epi_all/NicheNet/signaling_network.rds")
gr_network = readRDS("6.0.cellchat_Epi_all/NicheNet/gr_network.rds")

ligands_all = c("TSLP", "PTGS2")
targets_all = c("B2M", "FKBP5", "HSP90AA1", "MAF", "NDRG1",
                "KLK3", "TMSB10", "MTDH", "SOD3", "NPY", "MSMB",
                "RPS20", "NDUFB9", "CCK", "SARAF", "DEGS1", "MANF")

#active_ligand_target_links_df = best_upstream_ligands %>% 
#  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
#targets_all = unique(active_ligand_target_links_df$target)

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

# To render the graph: uncomment following line of code
#DiagrammeR::render_graph(graph_min_max, layout = "tree")

#export to cytoscape
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
head(data_source_network) 

output_path = "6.0.cellchat_Epi_all/NicheNet/"
write_output = TRUE 

# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}

weighted_signaling_network <- bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory"))

# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
  data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}

# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

if(write_output){
  bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% write_tsv(paste0(output_path,"annotation_table.txt"))
}

annotation <- bind_rows(specific_annotation_tbl,non_specific_annotation_tbl)
annotation_split <- split(annotation$gene, annotation$annotation)

# sankeyNetwork -----------------------------------


ligand_to_receptor <- weighted_networks$lr_sig %>% 
  filter(from %in% annotation_split$ligand) %>%
  filter(to %in% annotation_split$receptor)
ligand_to_receptor[4,3] <- 0.8

receptor_to_TF <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(to %in% annotation_split$"transcriptional regulator") %>%
  select(!layer)

receptor_to_target <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(to %in% annotation_split$target) %>%
  select(!layer)

TF_to_target <- weighted_signaling_network %>% 
  filter(from %in% receptor_to_TF$to) %>% 
  filter(to %in% annotation_split$target) %>%
  select(!layer)

receptor_to_TF <- receptor_to_TF[-which(receptor_to_TF$to %in% setdiff(receptor_to_TF$to, TF_to_target$from)),]

links <- rbind(ligand_to_receptor,receptor_to_TF,receptor_to_target,TF_to_target)
nodes <- data.frame(gene = unique(c(links$from, links$to)))
order <- nodes$gene
nodes <- merge(nodes, annotation, by="gene")
nodes <- nodes[match(order, nodes$gene),]

Node2index = list()
Node2index[nodes$gene] = 0:(length(nodes$gene)-1)
index <- data.frame(gene = names(Node2index), index = unlist(Node2index))

order_from <- links$from
from2 <- merge(data.frame(gene = links$from), index, by = "gene")
from2 <- from2[match(order_from, from2$gene),]

order_to <- links$to
to2 <- merge(data.frame(gene = links$to), index, by = "gene")
to2 <- to2[match(order_to, to2$gene),]

links <- links %>%
  mutate(source = from2$index) %>%
  mutate(target = to2$index)

links$color <- "grey"
nodes$color <- "#E46B4F"
nodes$color[nodes$annotation == "ligand"] <- "#053E7A"
nodes$color[nodes$annotation == "receptor"] <- "#4A7088"
nodes$color[nodes$annotation == "target"] <- "#A52A2A"

write.table(links, "6.0.cellchat_Epi_all/NicheNet/Sankey_links.txt", sep = "\t", quote = F, row.names = F)
write.table(nodes, "6.0.cellchat_Epi_all/NicheNet/Sankey_nodes.txt", sep = "\t", quote = F, row.names = F)

links <- read.table("C:/Users/Administrator/Desktop/Sankey_links.txt",sep = "\t",header=T)
nodes <- read.table("C:/Users/Administrator/Desktop/Sankey_nodes.txt",sep = "\t",header=T)

links$color <- "grey"
nodes$color <- "#E46B4F"
nodes$color[nodes$annotation == "ligand"] <- "#053E7A"
nodes$color[nodes$annotation == "receptor"] <- "#008B8B"
nodes$color[nodes$annotation == "target"] <- "#A52A2A"
color2project = paste(c(unique(nodes$color),unique(links$color)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

library(networkD3)
sankeyNetwork(Links = links, Nodes = nodes, 
              # 指定source、target、value和name：
              Source = "source",
              Target = "target", 
              Value = "weight", 
              NodeID = "gene", # 节点的名字
              NodeGroup = "annotation",
              # 调整配置：
              fontSize = 24, # 节点的字体大小
              nodeWidth = 40, # 节点的宽度
              nodePadding = 10, # 节点之间的距离
              colourScale = JS(my_color),
              width = 1200,
              height = 960,
              sinksRight = F
)




MisLinks = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/link.csv",sep = ",")
MisNodes = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/node.csv",sep = ",")



color2project = paste(c(unique(MisNodes$group_color),unique(MisLinks$group_color2)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

sankeyNetwork(Links = MisLinks, 
              Nodes = MisNodes,
              Source = "source2", 
              Target = "target2",
              Value ="value",
              NodeID = "name", 
              LinkGroup  = "group_color2", 
              colourScale = JS(my_color),
              fontSize = 10
)

# CAF_APOD ----------------------------------------------------

ligand_activities <- readRDS("6.0.cellchat_Epi_all/NicheNet_Inv/Results_sep/CAF_APOD/ligand_activities.Rds")

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

## Ligand-receptor network inference for top-ranked ligands -----------------------------

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
lr_network_top_matrix = lr_network_top_matrix[,best_upstream_ligands]

## Show a heatmap of the ligand-receptor interactions ----------------------

p_ligand_receptor_network = lr_network_top_matrix %>% 
  .[,rev(colnames(.))] %>% t() %>% 
  make_heatmap_ggplot("Prioritized ligands (CAF_APOD)",
                      "Expressed receptors (Cancer cells)", 
                      color = "#A52A2A", 
                      x_axis_position = "top",
                      legend_title = "Prior interaction potential")

p_ligand_receptor_network

ggsave("6.0.cellchat_Epi_all/NicheNet_Inv/Results_sep/CAF_APOD/lig_receptor_network.pdf", p_ligand_receptor_network, height = 5, width = 8)

## Inferring ligand-to-target signaling paths -------------------------------

ligand_tf_matrix = readRDS("6.0.cellchat_Epi_all/NicheNet/ligand_tf_matrix.rds")
sig_network = readRDS("6.0.cellchat_Epi_all/NicheNet/signaling_network.rds")
gr_network = readRDS("6.0.cellchat_Epi_all/NicheNet/gr_network.rds")

ligands_all = c("TSLP", "PTGS2")
targets_all = c("B2M", "FKBP5", "HSP90AA1", "MAF", "NDRG1",
                "KLK3", "TMSB10", "MTDH", "SOD3", "NPY", "MSMB",
                "RPS20", "NDUFB9", "CCK", "SARAF", "DEGS1", "MANF")

#active_ligand_target_links_df = best_upstream_ligands %>% 
#  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
#targets_all = unique(active_ligand_target_links_df$target)

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

# To render the graph: uncomment following line of code
#DiagrammeR::render_graph(graph_min_max, layout = "tree")

#export to cytoscape
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
head(data_source_network) 

output_path = "6.0.cellchat_Epi_all/NicheNet/"
write_output = TRUE 

# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}

weighted_signaling_network <- bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory"))

# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
  data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}

# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

if(write_output){
  bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% write_tsv(paste0(output_path,"annotation_table.txt"))
}

annotation <- bind_rows(specific_annotation_tbl,non_specific_annotation_tbl)
annotation_split <- split(annotation$gene, annotation$annotation)

# sankeyNetwork -----------------------------------


ligand_to_receptor <- weighted_networks$lr_sig %>% 
  filter(from %in% annotation_split$ligand) %>%
  filter(to %in% annotation_split$receptor)
ligand_to_receptor[4,3] <- 0.8

receptor_to_TF <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(to %in% annotation_split$"transcriptional regulator") %>%
  select(!layer)

receptor_to_target <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(to %in% annotation_split$target) %>%
  select(!layer)

TF_to_target <- weighted_signaling_network %>% 
  filter(from %in% receptor_to_TF$to) %>% 
  filter(to %in% annotation_split$target) %>%
  select(!layer)

receptor_to_TF <- receptor_to_TF[-which(receptor_to_TF$to %in% setdiff(receptor_to_TF$to, TF_to_target$from)),]

links <- rbind(ligand_to_receptor,receptor_to_TF,receptor_to_target,TF_to_target)
nodes <- data.frame(gene = unique(c(links$from, links$to)))
order <- nodes$gene
nodes <- merge(nodes, annotation, by="gene")
nodes <- nodes[match(order, nodes$gene),]

Node2index = list()
Node2index[nodes$gene] = 0:(length(nodes$gene)-1)
index <- data.frame(gene = names(Node2index), index = unlist(Node2index))

order_from <- links$from
from2 <- merge(data.frame(gene = links$from), index, by = "gene")
from2 <- from2[match(order_from, from2$gene),]

order_to <- links$to
to2 <- merge(data.frame(gene = links$to), index, by = "gene")
to2 <- to2[match(order_to, to2$gene),]

links <- links %>%
  mutate(source = from2$index) %>%
  mutate(target = to2$index)

links$color <- "grey"
nodes$color <- "#E46B4F"
nodes$color[nodes$annotation == "ligand"] <- "#053E7A"
nodes$color[nodes$annotation == "receptor"] <- "#4A7088"
nodes$color[nodes$annotation == "target"] <- "#A52A2A"

write.table(links, "6.0.cellchat_Epi_all/NicheNet/Sankey_links.txt", sep = "\t", quote = F, row.names = F)
write.table(nodes, "6.0.cellchat_Epi_all/NicheNet/Sankey_nodes.txt", sep = "\t", quote = F, row.names = F)

links <- read.table("C:/Users/Administrator/Desktop/Sankey_links.txt",sep = "\t",header=T)
nodes <- read.table("C:/Users/Administrator/Desktop/Sankey_nodes.txt",sep = "\t",header=T)

links$color <- "grey"
nodes$color <- "#E46B4F"
nodes$color[nodes$annotation == "ligand"] <- "#053E7A"
nodes$color[nodes$annotation == "receptor"] <- "#008B8B"
nodes$color[nodes$annotation == "target"] <- "#A52A2A"
color2project = paste(c(unique(nodes$color),unique(links$color)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

library(networkD3)
sankeyNetwork(Links = links, Nodes = nodes, 
              # 指定source、target、value和name：
              Source = "source",
              Target = "target", 
              Value = "weight", 
              NodeID = "gene", # 节点的名字
              NodeGroup = "annotation",
              # 调整配置：
              fontSize = 24, # 节点的字体大小
              nodeWidth = 40, # 节点的宽度
              nodePadding = 10, # 节点之间的距离
              colourScale = JS(my_color),
              width = 1200,
              height = 960,
              sinksRight = F
)




MisLinks = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/link.csv",sep = ",")
MisNodes = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/node.csv",sep = ",")



color2project = paste(c(unique(MisNodes$group_color),unique(MisLinks$group_color2)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

sankeyNetwork(Links = MisLinks, 
              Nodes = MisNodes,
              Source = "source2", 
              Target = "target2",
              Value ="value",
              NodeID = "name", 
              LinkGroup  = "group_color2", 
              colourScale = JS(my_color),
              fontSize = 10
)

# CAF_CCL2 ----------------------------------------------------

ligand_activities <- readRDS("6.0.cellchat_Epi_all/NicheNet_Inv/Results_sep/CAF_CCL2/ligand_activities.Rds")

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

## Ligand-receptor network inference for top-ranked ligands -----------------------------

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
lr_network_top_matrix = lr_network_top_matrix[,best_upstream_ligands]

## Show a heatmap of the ligand-receptor interactions ----------------------

p_ligand_receptor_network = lr_network_top_matrix %>% 
  .[,rev(colnames(.))] %>% t() %>% 
  make_heatmap_ggplot("Prioritized ligands (CAF_CCL2)",
                      "Expressed receptors (Cancer cells)", 
                      color = "#A52A2A", 
                      x_axis_position = "top",
                      legend_title = "Prior interaction potential")

p_ligand_receptor_network

ggsave("6.0.cellchat_Epi_all/NicheNet_Inv/Results_sep/CAF_CCL2/lig_receptor_network.pdf", p_ligand_receptor_network, height = 5, width = 8)

## Inferring ligand-to-target signaling paths -------------------------------

ligand_tf_matrix = readRDS("6.0.cellchat_Epi_all/NicheNet/ligand_tf_matrix.rds")
sig_network = readRDS("6.0.cellchat_Epi_all/NicheNet/signaling_network.rds")
gr_network = readRDS("6.0.cellchat_Epi_all/NicheNet/gr_network.rds")

ligands_all = c("TSLP", "PTGS2")
targets_all = c("B2M", "FKBP5", "HSP90AA1", "MAF", "NDRG1",
                "KLK3", "TMSB10", "MTDH", "SOD3", "NPY", "MSMB",
                "RPS20", "NDUFB9", "CCK", "SARAF", "DEGS1", "MANF")

#active_ligand_target_links_df = best_upstream_ligands %>% 
#  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
#targets_all = unique(active_ligand_target_links_df$target)

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

# To render the graph: uncomment following line of code
#DiagrammeR::render_graph(graph_min_max, layout = "tree")

#export to cytoscape
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
head(data_source_network) 

output_path = "6.0.cellchat_Epi_all/NicheNet/"
write_output = TRUE 

# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}

weighted_signaling_network <- bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory"))

# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
  data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}

# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

if(write_output){
  bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% write_tsv(paste0(output_path,"annotation_table.txt"))
}

annotation <- bind_rows(specific_annotation_tbl,non_specific_annotation_tbl)
annotation_split <- split(annotation$gene, annotation$annotation)

# sankeyNetwork -----------------------------------


ligand_to_receptor <- weighted_networks$lr_sig %>% 
  filter(from %in% annotation_split$ligand) %>%
  filter(to %in% annotation_split$receptor)
ligand_to_receptor[4,3] <- 0.8

receptor_to_TF <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(to %in% annotation_split$"transcriptional regulator") %>%
  select(!layer)

receptor_to_target <- weighted_signaling_network %>% 
  filter(from %in% ligand_to_receptor$to) %>% 
  filter(to %in% annotation_split$target) %>%
  select(!layer)

TF_to_target <- weighted_signaling_network %>% 
  filter(from %in% receptor_to_TF$to) %>% 
  filter(to %in% annotation_split$target) %>%
  select(!layer)

receptor_to_TF <- receptor_to_TF[-which(receptor_to_TF$to %in% setdiff(receptor_to_TF$to, TF_to_target$from)),]

links <- rbind(ligand_to_receptor,receptor_to_TF,receptor_to_target,TF_to_target)
nodes <- data.frame(gene = unique(c(links$from, links$to)))
order <- nodes$gene
nodes <- merge(nodes, annotation, by="gene")
nodes <- nodes[match(order, nodes$gene),]

Node2index = list()
Node2index[nodes$gene] = 0:(length(nodes$gene)-1)
index <- data.frame(gene = names(Node2index), index = unlist(Node2index))

order_from <- links$from
from2 <- merge(data.frame(gene = links$from), index, by = "gene")
from2 <- from2[match(order_from, from2$gene),]

order_to <- links$to
to2 <- merge(data.frame(gene = links$to), index, by = "gene")
to2 <- to2[match(order_to, to2$gene),]

links <- links %>%
  mutate(source = from2$index) %>%
  mutate(target = to2$index)

links$color <- "grey"
nodes$color <- "#E46B4F"
nodes$color[nodes$annotation == "ligand"] <- "#053E7A"
nodes$color[nodes$annotation == "receptor"] <- "#4A7088"
nodes$color[nodes$annotation == "target"] <- "#A52A2A"

write.table(links, "6.0.cellchat_Epi_all/NicheNet/Sankey_links.txt", sep = "\t", quote = F, row.names = F)
write.table(nodes, "6.0.cellchat_Epi_all/NicheNet/Sankey_nodes.txt", sep = "\t", quote = F, row.names = F)

links <- read.table("C:/Users/Administrator/Desktop/Sankey_links.txt",sep = "\t",header=T)
nodes <- read.table("C:/Users/Administrator/Desktop/Sankey_nodes.txt",sep = "\t",header=T)

links$color <- "grey"
nodes$color <- "#E46B4F"
nodes$color[nodes$annotation == "ligand"] <- "#053E7A"
nodes$color[nodes$annotation == "receptor"] <- "#008B8B"
nodes$color[nodes$annotation == "target"] <- "#A52A2A"
color2project = paste(c(unique(nodes$color),unique(links$color)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

library(networkD3)
sankeyNetwork(Links = links, Nodes = nodes, 
              # 指定source、target、value和name：
              Source = "source",
              Target = "target", 
              Value = "weight", 
              NodeID = "gene", # 节点的名字
              NodeGroup = "annotation",
              # 调整配置：
              fontSize = 24, # 节点的字体大小
              nodeWidth = 40, # 节点的宽度
              nodePadding = 10, # 节点之间的距离
              colourScale = JS(my_color),
              width = 1200,
              height = 960,
              sinksRight = F
)




MisLinks = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/link.csv",sep = ",")
MisNodes = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/sankeyNetwork/node.csv",sep = ",")



color2project = paste(c(unique(MisNodes$group_color),unique(MisLinks$group_color2)),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')

sankeyNetwork(Links = MisLinks, 
              Nodes = MisNodes,
              Source = "source2", 
              Target = "target2",
              Value ="value",
              NodeID = "name", 
              LinkGroup  = "group_color2", 
              colourScale = JS(my_color),
              fontSize = 10
)
