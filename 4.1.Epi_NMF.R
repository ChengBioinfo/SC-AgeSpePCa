color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")

# subset analysis ---------------------------------------------------------------------

PCSC <- readRDS("./4.3.Doublet/PCSC_inferred_doublet.Rds")

# epithelial --------------------------------------------------------------------

folder <- "5.1.1.epithelial_NMF"
if(!dir.exists(folder)){
  dir.create(folder)
}

# A. NMF for all samples ---------------------------------------------

# 1. subset ----------------------------------------------------------------

Epi <- subset(PCSC, inferCNVCellStat == "Cancer")

# 2. NMF ---------------------------------------------------------

Epi <- NormalizeData(Epi)

Epi <- FindVariableFeatures(Epi,
                            selection.method = "vst",
                            nfeatures = 2000,
                            verbose = T) %>%
  ScaleData(do.center = F)
exp <- Epi@assays$RNA@scale.data
res <- nmf(exp, rank = 15, method = "snmf/r", seed = 'nndsvd')
saveRDS(res, "./5.1.1.epithelial_NMF/Epi_nmf_res.Rds")

fs <- extractFeatures(res, 50L)
fs <- lapply(fs, function(x) rownames(res)[x])
fs <- do.call("rbind", fs)
rownames(fs) <- paste0("cluster", 1:15)
fs <- t(fs)
write.csv(fs, "./5.1.1.epithelial_NMF/NMF_TopGenes.csv", row.names = F)
DT::datatable(t(fs))

# 3. harmony ---------------------------------------------------

dim <- 1:15

## 降维
Epi <- RunPCA(Epi, verbose = T)
Epi@reductions$nmf <- Epi@reductions$pca
Epi@reductions$nmf@cell.embeddings <- t(coef(res))    
Epi@reductions$nmf@feature.loadings <- basis(res)

pdf("./5.1.1.epithelial_NMF/harmony_convergence.pdf")
Epi <- RunHarmony(Epi, group.by.vars = "orig.ident", reduction = "nmf", plot_convergence = T)
dev.off()

Epi <- RunUMAP(Epi, reduction = "harmony", dims = 1:15)

# 4.clustering -------------------------------------------------------

Epi <- FindNeighbors(Epi, reduction = "harmony", 
                     dims = 1:15)
for (ress in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1)) {
  Epi <- FindClusters(Epi, resolution = ress)}

apply(Epi@meta.data[, grep("RNA_snn_res", colnames(Epi@meta.data))],2,table)

library(clustree)

p2_tree=clustree(Epi@meta.data, prefix = "RNA_snn_res.") +
  scale_color_npg()
ggsave(plot=p2_tree, filename="./5.1.1.epithelial_NMF/Tree_diff_resolution.pdf",width = 7, height = 12)

Idents(Epi) <- "RNA_snn_res.0.5"

pdf("./5.1.1.epithelial_NMF/CellCluster-UMAPPlot_res0.5.pdf",width = 7,height = 7)
DimPlot(Epi, reduction = "umap", pt.size = 0.1, label = T) + NoLegend()
dev.off()

pdf("./5.1.1.epithelial_NMF/CellCluster-UMAPPlot_SampleType.pdf",width = 7,height = 7)
DimPlot(Epi, reduction = "umap", group.by = "SampleType", pt.size = 0.1, label = T) +
  scale_color_npg()
dev.off()

pdf("./5.1.1.epithelial_NMF/CellCluster-UMAPPlot_Sample.pdf",width = 7,height = 7)
DimPlot(Epi, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F, shuffle = T) + 
  scale_color_npg()
dev.off()

pdf("./5.1.1.epithelial_NMF/cellcycle/CellCluster-UMAPPlot_cellcycle_res0.3.pdf",width = 7,height = 7)
DimPlot(Epi, reduction = "umap", group.by = "Phase", pt.size = 0.1, label = T) + scale_color_npg()
dev.off()

pdf("./5.1.1.epithelial_NMF/cellcycle/CellCluster-UMAPPlot_SampleTypeSplit_cellcycle_res0.3.pdf",width = 13,height = 7)
DimPlot(Epi, reduction = "umap", group.by = "Phase", split.by = "SampleType", pt.size = 0.1, label = T) + scale_color_npg()
dev.off()

pdf("./5.1.1.epithelial_NMF/cellcycle/CellCluster-UMAPPlot_SampleSplit_cellcycle_res0.3.pdf",width = 36,height = 15)
DimPlot(Epi, reduction = "umap", group.by = "Phase", split.by = "orig.ident", pt.size = 0.1, label = T, ncol = 5) + scale_color_npg()
dev.off()

pdf("./5.1.1.epithelial_NMF/cellcycle/CellCluster-UMAPPlot_ClusterSplit_cellcycle_res0.3.pdf",width = 29,height = 22)
DimPlot(Epi, reduction = "umap", group.by = "Phase", split.by = "RNA_snn_res.0.5", pt.size = 0.1, label = T, ncol = 4) + scale_color_npg()
dev.off()

saveRDS(Epi, "./5.1.1.epithelial_NMF/Epi_clustered.Rds")
Epi <- readRDS("./5.1.1.epithelial_NMF/Epi_clustered.Rds")

# B. cNMF for individual sample ---------------------------------------------

# 1. subset ----------------------------------------------------------------

Epi <- subset(PCSC, inferCNVCellStat == "Cancer")
for (i in as.character(unique(Epi@meta.data$orig.ident))) {
  subset <- subset(Epi, subset = orig.ident == i)
  exp <- as.data.frame(subset[["RNA"]]@counts)
  exp <- exp[rowSums(exp) > 0,]
  exp <- exp[!str_detect(rownames(exp), "^MT-"),]
  exp <- data.frame(t(exp))
  write.table(exp, paste0("./5.1.1.epithelial_NMF/count_data/", i, ".count.txt"), quote = F, sep = "\t", row.names = T, col.names = T)
}

# 2.cNMF ---------------------------------------------------

library(reticulate)

### 调用子环境的python
use_condaenv(condaenv = "cnmf_env", required = T)
py_config() #如果显示cnmf_env环境里面的python就OK

# 3:15
source("./5.1.1.epithelial_NMF/1.R")
step1(dir_input = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/count_data", dir_output = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res1", k=3:15, iteration = 200)

source("./5.1.1.epithelial_NMF/2.R")
step2(dir_input = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res1",dir_output = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res2",dir_count = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/count_data",usage_filter = 0.03,top_gene = 30,cor_min = 0,cor_max = 0.6)

# 3:20
source("./5.1.1.epithelial_NMF/1.R")
step1(dir_input = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/count_data", dir_output = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res1.1", k=3:20, iteration = 100)

## k_selection
source("./5.1.1.epithelial_NMF/2.R")
step2(dir_input = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res1.1",dir_output = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res2.1",dir_count = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/count_data",usage_filter = 0.03,top_gene = 30,cor_min = 0,cor_max = 0.6)

## k_selection2
source("./5.1.1.epithelial_NMF/2_k_selection2.R")
step2(dir_input = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res1.1",dir_output = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res2.1.2",dir_count = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/count_data",usage_filter = 0.03,top_gene = 30,cor_min = 0,cor_max = 0.6)

## k_selection3
source("./5.1.1.epithelial_NMF/2_k_selection3.R")
step2(dir_input = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res1.1",dir_output = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res2.1.3",dir_count = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/count_data",usage_filter = 0.03,top_gene = 50,cor_min = 0,cor_max = 0.6)

## k_selection4
source("./5.1.1.epithelial_NMF/2_k_selection4.R")
step2(dir_input = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res1.1",dir_output = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res2.1.4",dir_count = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/count_data",usage_filter = 0.03,top_gene = 50,cor_min = 0,cor_max = 0.6)

## k_selection5
source("./5.1.1.epithelial_NMF/2_k_selection5.R")
step2(dir_input = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res1.1",dir_output = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/res2.1.5",dir_count = "C:/Users/Administrator/Desktop/process/5.1.1.epithelial_NMF/count_data",usage_filter = 0.03,top_gene = 50,cor_min = 0,cor_max = 0.6)

# 3:20, 3000 features
source("./5.1.1.epithelial_NMF/1.R")
step1(dir_input = "./5.1.1.epithelial_NMF/count_data", dir_output = "./5.1.1.epithelial_NMF/res1.2", k=3:20, iteration = 200, feature = 3000)

## k_selection
source("./5.1.1.epithelial_NMF/2.R")
step2(dir_input = "./5.1.1.epithelial_NMF/res1.2",dir_output = "./5.1.1.epithelial_NMF/res2.2",dir_count = "./5.1.1.epithelial_NMF/count_data",usage_filter = 0.03,top_gene = 50,cor_min = 0,cor_max = 0.6)

#3.subsequent analysis --------------------------------------------

## K_selection3 ---------------------------------

topgenes <- read.table("./5.1.1.epithelial_NMF/res2.1.3/program_topngene.txt", row.names = 1, header = T, sep = "\t")

MP1 <- topgenes[, c("EOPC4_1","LOPC3_2","EOPC3_2","LOPC1_1","LOPC2_2","LOPC4_2","EOPC3_1","LOPC1_7","LOPC6_2","EOPC1_13","EOPC2_2","LOPC4_9")]
MP2 <- topgenes[, c("EOPC3_4","EOPC4_4","EOPC1_7","LOPC3_3","LOPC4_7","LOPC5_1","LOPC2_3","LOPC3_8")]
MP3 <- topgenes[, c("EOPC3_3","EOPC4_5","LOPC3_10","EOPC4_3","LOPC3_6")]
MP4 <- topgenes[, c("EOPC1_2","LOPC5_2","LOPC1_4","LOPC6_5")]
MP5 <- topgenes[, c("EOPC1_11","LOPC5_6","EOPC2_7","EOPC4_6","LOPC3_7","LOPC1_2","LOPC4_5")]
MP6 <- topgenes[, c("EOPC2_9","LOPC4_1","LOPC5_5","LOPC1_5","LOPC4_8","EOPC3_7","LOPC1_6")]

MPlist <- list(MP1,MP2,MP3,MP4,MP5,MP6)
MPlist <- lapply(MPlist, function(x) {unique(as.character(as.matrix(x)))})

saveRDS(MPlist, "./5.1.1.epithelial_NMF/process_k3/MPlist.Rds")
MPlist <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPlist.Rds")

### pct_genes -------------------------------------------

mtx <- as.matrix(Epi@assays$RNA@counts)

for (i in 1:length(MPlist)) {
  change <- setdiff(MPlist[[i]], rownames(mtx))
  change <- gsub("\\.", "-", change)
  a <- c(MPlist[[i]], change)
  exp <- mtx[rownames(mtx) %in% a,]
  Markers <- c()
  for (x in 1:nrow(exp)) {
    pos <- ifelse(exp[x,] > 0, "Positive", "Negative")
    pct <- length(pos[pos=="Positive"])/length(pos)
    if (pct > 0.3) {Markers <- c(Markers, rownames(exp)[x])}
  }
  assign(paste0("MP", i, "_Markers"), Markers)
}

Markers <- list(MP1 = MP1_Markers,
                MP2 = MP2_Markers,
                MP3 = MP3_Markers,
                MP4 = MP4_Markers,
                MP5 = MP5_Markers,
                MP6 = MP6_Markers)
saveRDS(Markers, "./5.1.1.epithelial_NMF/process_k3/MP_markers.Rds")
Markers <- readRDS("./5.1.1.epithelial_NMF/process_k3/MP_markers.Rds")

for (i in 1:length(Markers)) {
  
  genes <- as.data.frame(Markers[[i]])
  colnames(genes) <- "Genes"
  write.table(genes, paste0("./5.1.1.epithelial_NMF/process_k3/MP_Markers/Markers", names(Markers)[i], ".txt"),
                            sep = "\t", quote = F, row.names = F)
  
}

### difference of MP score between EO and LO in SC ------------------------------

#### AUCell -----------------------------------------------

expMtx <- as.matrix(Epi@assays$RNA@data)
cells_rankings <- AUCell_buildRankings(expMtx, nCores=10, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC.Rds")

pdf("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_threshould.pdf")
set.seed(123)
par(mfrow=c(2,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
dev.off()

cells_AUC <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC.Rds")
AUCell_mtx <- as.matrix(cells_AUC@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
group <- data.frame(SampleType = Epi$SampleType)
rownames(group) <- rownames(Epi@meta.data)
mtx <- cbind(group, AUCell_mtx_t)

library(reshape2)
data <- melt(mtx)
colnames(data)[2] <- "MP"

library(ggpubr)

P = ggplot(data, aes(x = MP, y = value, fill = SampleType)) +
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
  ggtitle("SC data") +
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

pdf("./5.1.1.epithelial_NMF/process_k3/boxplot_MPScore_AUC_diff.pdf", width=5, height=3)
print(P)
dev.off()

#### irGSEA -------------------------------------------

library(irGSEA)

Epi_irGSEA <- irGSEA.score(object = Epi, assay = "RNA", slot = "data", 
                           seeds = 123, ncores = 20, min.cells = 3, 
                           min.feature = 0, custom = T, geneset = Markers, 
                           method = c("AUCell", "UCell", "singscore", "ssgsea"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                           kcdf = 'Gaussian')

saveRDS(Epi_irGSEA, "./5.1.1.epithelial_NMF/process_k3/Epi_MP_irGSEA.Rds")

# integrate differentially enriched genesets
result.dge <- irGSEA.integrate(object = Epi_irGSEA,
                               group.by = "SampleType",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))
result.dge.unlist <- Reduce(rbind, result.dge[-5])
write.table(result.dge.unlist, "./5.1.1.epithelial_NMF/process_k3/Epi_MP_irGSEA_dge.txt", row.names = F, sep = "\t", quote = F)

### On/Off of MP --------------------------------------------
Epi_MP_AUCell <- CreateSeuratObject(AUCell_mtx,
                                    project = "Epi_MP", 
                                    names.field = 1,
                                    names.delim = "_")
Epi_MP_AUCell@meta.data <- Epi@meta.data
Epi_MP_AUCell@reductions <- Epi@reductions

threshould <- data.frame(MP = c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6"),
                         cutoff = c(0.4, 0.81, 0.51, 0.58, 0.38, 0.54))

for (i in rownames(AUCell_mtx)) {
  Epi_MP_AUCell$proportion <- ifelse(Epi_MP_AUCell@assays$RNA@data[i,] > threshould[which(threshould$MP == i),2], "On", "Off")
  pos <- table(Epi_MP_AUCell$SampleType, Epi_MP_AUCell$proportion)
  chi <- chisq.test(pos)
  p <- signif(chi$p.value,3)
  pos <- melt(pos)
  colnames(pos) <- c("SampleType", "Status", "nCells")
  barplot <- ggplot(pos, aes(x=SampleType, y=nCells,fill=Status)) +
    geom_bar(stat="identity", position = 'fill') +
    theme_bw() +
    scale_fill_manual(values = c("#BBBBBB", "#556B2FCC")) +
    annotate(geom="text", x=1.5, y=1.05, label=paste0("Chi-squre P-value = ", p), size = 2) +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ylab("Pct. On/Off status") +
    xlab("") + 
    ggtitle(i) +
    NoLegend()
  assign(paste0(i, "_OnOff"), barplot)
}

list <- list(MP1_OnOff, MP2_OnOff, MP3_OnOff, MP4_OnOff, MP5_OnOff, MP6_OnOff)

pdf("./5.1.1.epithelial_NMF/process_k3/MP_On_Off_pct.pdf",
    height = 5, width = 6)
grid.arrange(grobs = list, ncol= 3)
dev.off()

### MP enrichment ---------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)

hsets <- read.gmt("./5.1.1.epithelial_NMF/hallmark_cancersea.gmt")
enrich.result <- data.frame()

for (i in 1:length(MPlist)) {
  tmp <- enricher(MPlist[[i]], TERM2GENE = hsets, pAdjustMethod = "fdr")
  if (is.null(tmp)) {
    next
  }
  
  tmp1 <- tmp@result[tmp@result$qvalue < 0.05,]
  tmp1$program <- paste0("MP", i)
  rownames(tmp1)=NULL
  enrich.result <- rbind(enrich.result,tmp1)
}

write.xlsx(enrich.result, "./5.1.1.epithelial_NMF/process_k3/MP_enrichment_HallmarkCancersea.xlsx",row.names = F)

enrich.result <- read.xlsx2("./5.1.1.epithelial_NMF/process_k3/MP_enrichment_HallmarkCancersea.xlsx", sheetIndex = 1)
enrich.result$p.adjust <- as.numeric(enrich.result$p.adjust)
plot <- data.frame(term = enrich.result$Description, logP = -log10(enrich.result$p.adjust), MP = enrich.result$program)

plot$term[which(duplicated(plot$term))] <-
  paste0(plot$term[which(duplicated(plot$term))], "_2")
plot$term[which(duplicated(plot$term))] <-
  paste0(plot$term[which(duplicated(plot$term))], "_3")

plot$term <- factor(plot$term, levels = unique(plot$term))
plot$MP <- factor(plot$MP, levels = unique(plot$MP))

P <- ggplot(data = plot, aes(x = term, y = logP, fill = MP)) +
  geom_bar(stat='identity') +
  ggtitle("MP enrichment (CancerSEA)") +
  xlab("") +
  ylab("-log10(P.adjust)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 10),
        axis.text.y=element_text(size=10, colour="black"), 
        axis.title.y=element_text(size = 12),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#CD7054", "#2F4F4F", "#7A378B", "#008B8B", "#4A7088", "#8B8B00"))

ggsave("MP_enrichment_HallmarkCancersea.pdf", P, width = 10, height = 8)

### intersect of markers and DEGs ----------------------------------------------

DEG_TCGA <- read.table("0.TCGA/Bulk_DEG/TCGA/EOvsLO_all.txt", row.names = 1, header = T, sep = "\t")
diff_TCGA <- DEG_TCGA[DEG_TCGA$Pvalue < 0.01,]
DEG_GSE157547 <- read.table("0.TCGA/Bulk_DEG/GSE157547/EOvsLO_all.txt", row.names = 1, header = T, sep = "\t")
diff_GSE157547 <- DEG_GSE157547[DEG_GSE157547$Pvalue < 0.01,]
DEG_GSE141551 <- read.table("0.TCGA/Bulk_DEG/GSE141551/EOvsLO_all.txt", row.names = 1, header = T, sep = "\t")
diff_GSE141551 <- DEG_GSE141551[DEG_GSE141551$Pvalue < 0.01,]

BulkList <- list(TCGA = DEG_TCGA, GSE157547 = DEG_GSE157547, GSE141551 = DEG_GSE141551)
inte_markers <- data.frame()
name <- c()
for (i in 1:length(Markers)) {
  for (j in 1:length(BulkList)) {
    inte <- BulkList[[j]][intersect(Markers[[i]], rownames(BulkList[[j]])),]
    inte <- inte[order(inte$Pvalue),]
    inte <- inte[inte$Pvalue < 0.05,]
    tmp <- data.frame(gene = rownames(inte))
    tmp <- as.data.frame(t(tmp))
    inte_markers <- rbind.fill(inte_markers, tmp)
    name <- c(name, paste0(names(Markers)[i], "_", names(BulkList)[j]))
  }
}
rownames(inte_markers) <- name
inte_markers <- t(inte_markers)
write.table(inte_markers, "./5.1.1.epithelial_NMF/process_k3/intesect_markers_and_Bulk_DEG.txt",
            row.names = F, sep = "\t")

### difference of MP AUCell score between EO and LO in Bulk ----------------------------------

library(reshape2)

#### TCGA ----------------------------

expMtx <- read.table("../../SC/0.TCGA/mRNA_Matrix_TPM.txt", row.names = 1, header = T, sep = "\t")
expMtx <- as.matrix(expMtx)
cells_rankings <- AUCell_buildRankings(expMtx, nCores=10, plotStats=TRUE)
cells_AUC_TCGA <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC_TCGA, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_TCGA.Rds")

pdf("./5.1.1.epithelial_NMF/process_k3/TCGA_MP_AUC_threshould.pdf")
set.seed(123)
par(mfrow=c(2,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC_TCGA, plotHist=TRUE, assign=TRUE)
dev.off()

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_TCGA.Rds")
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
data_TCGA <- melt(mtx)
colnames(data_TCGA)[2] <- "MP"

mtx_TCGA_logit <- mtx[,c(1,2)]
mtx_TCGA_logit$Type <- ifelse(mtx_TCGA_logit$Type == "EOPC", 1, 0)

P = ggplot(data_TCGA, aes(x = MP, y = value, fill = Type)) +
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
  ggtitle("TCGA") +
  ylab("AUCell score") +
  xlab("") + 
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")

pdf("./5.1.1.epithelial_NMF/process_k3/boxplot_MPScore_AUC_diff_TCGA.pdf", width=5, height=3)
print(P)
dev.off()

#### GSE141551 ---------------------------------

exp <- read.table("./0.TCGA/GSE141551/MergeExpro_contrib1-GPL14951.txt", row.names = 1, header = T, sep = "\t")
exp <- as.matrix(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE141551 <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE141551, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE141551.Rds")

cells_AUC_GSE141551 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE141551.Rds")
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
data_GSE141551 <- melt(mtx)
colnames(data_GSE141551)[2] <- "MP"

mtx_GSE141551_logit <- mtx[,c(1,2)]
mtx_GSE141551_logit$Type <- ifelse(mtx_GSE141551_logit$Type == "EOPC", 1, 0)

#### GSE157547 ----------------------------------------

exp <- read.table("./0.TCGA/GSE157547/MergeExpro_contrib1-GPL5175.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
exp <- 2^exp
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE157547 <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE157547, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE157547.Rds")

cells_AUC_GSE157547 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE157547.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE157547@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE157547/SampleInfo_contrib1-GPL5175.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE157547 <- melt(mtx)
colnames(data_GSE157547)[2] <- "MP"

mtx_GSE157547_logit <- mtx[,c(1,2)]
mtx_GSE157547_logit$Type <- ifelse(mtx_GSE157547_logit$Type == "EOPC", 1, 0)

#### GSE153352 --------------------------------------------------

exp <- read.table("./0.TCGA/GSE153352/MergeExpro_contrib1-GPL5175.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
exp <- 2^exp
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE153352 <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE153352, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE153352.Rds")

cells_AUC_GSE153352 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE153352.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE153352@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE153352/SampleInfo_contrib1-GPL5175.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE153352 <- melt(mtx)
colnames(data_GSE153352)[2] <- "MP"

mtx_GSE153352_logit <- mtx[,c(1,2)]
mtx_GSE153352_logit$Type <- ifelse(mtx_GSE153352_logit$Type == "EOPC", 1, 0)

#### GSE88808 -----------------------------------------

exp <- read.table("./0.TCGA/GSE88808/MergeExpro_contrib1-GPL22571.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE88808 <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE88808, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE88808.Rds")

cells_AUC_GSE88808 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE88808.Rds")
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
data_GSE88808 <- melt(mtx)
colnames(data_GSE88808)[2] <- "MP"

mtx_GSE88808_logit <- mtx[,c(1,2)]
mtx_GSE88808_logit$Type <- ifelse(mtx_GSE88808_logit$Type == "EOPC", 1, 0)

#### GSE62116 ----------------------------------------

exp <- read.table("./0.TCGA/GSE62116/MergeExpro_contrib1-GPL5188.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
exp <- 2^exp
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE62116 <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE62116, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE62116.Rds")

cells_AUC_GSE62116 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE62116.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE62116@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE62116/SampleInfo_contrib1-GPL5188.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE62116 <- melt(mtx)
colnames(data_GSE62116)[2] <- "MP"

mtx_GSE62116_logit <- mtx[,c(1,2)]
mtx_GSE62116_logit$Type <- ifelse(mtx_GSE62116_logit$Type == "EOPC", 1, 0)

#### GSE21034 ----------------------------------------

exp <- read.table("./0.TCGA/GSE21034/MergeExpro_contrib1-GPL5188.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE21034 <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE21034, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE21034.Rds")

cells_AUC_GSE21034 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE21034.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE21034@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE21034/GSE21034_clinical.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"DxAge"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE21034 <- melt(mtx)
colnames(data_GSE21034)[2] <- "MP"

mtx_GSE21034_logit <- mtx[,c(1,2)]
mtx_GSE21034_logit$Type <- ifelse(mtx_GSE21034_logit$Type == "EOPC", 1, 0)

#### GSE183019 ----------------------------------------

exp <- read.table("./0.TCGA/GSE183019/GSE183019_mRNA_processed_TPM.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE183019 <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE183019, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE183019.Rds")

cells_AUC_GSE183019 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE183019.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE183019@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE183019/clinical.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"Age"])
rownames(age) <- gsub("-", ".", clinical$Sample)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE183019 <- melt(mtx)
colnames(data_GSE183019)[2] <- "MP"

mtx_GSE183019_logit <- mtx[,c(1,2)]
mtx_GSE183019_logit$Type <- ifelse(mtx_GSE183019_logit$Type == "EOPC", 1, 0)

#### GSE201284 ----------------------------------------

exp <- read.table("./0.TCGA/GSE201284/GSE201284_mRNA_processed_TPM.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_rankings <- AUCell_buildRankings(exp, nCores=10, plotStats=TRUE)
cells_AUC_GSE201284 <- AUCell_calcAUC(Markers, cells_rankings, nCores = 20)
saveRDS(cells_AUC_GSE201284, "./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE201284.Rds")

cells_AUC_GSE201284 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_GSE201284.Rds")
AUCell_mtx <- as.matrix(cells_AUC_GSE201284@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/GSE201284/clinical.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"Age"])
rownames(age) <- gsub("-", ".", clinical$Sample)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, AUCell_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE201284 <- melt(mtx)
colnames(data_GSE201284)[2] <- "MP"

mtx_GSE201284_logit <- mtx[,c(1,2)]
mtx_GSE201284_logit$Type <- ifelse(mtx_GSE201284_logit$Type == "EOPC", 1, 0)

#### boxplot ----------------------------------------

library(ggpubr)

data_list <- list(TCGA=data_TCGA, GSE183019 = data_GSE183019, GSE201284 = data_GSE201284, 
                  GSE141551=data_GSE141551, GSE153352=data_GSE153352, GSE157547=data_GSE157547, 
                  GSE88808 = data_GSE88808, GSE62116 = data_GSE62116, GSE21034 = data_GSE21034)
Plist <- list()
for (i in 1:length(data_list)) {
  P = ggplot(data_list[[i]], aes(x = MP, y = value, fill = Type)) +
    scale_fill_manual(values = c("#053E7AAA","#A52A2AAA")) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          legend.text=element_text(colour="black", size=8),
          legend.title=element_text(colour="black", size=10),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(names(data_list)[i]) +
    ylab("AUCell score") +
    xlab("") + 
    stat_compare_means(aes(group=Type),
                       method="wilcox.test",
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                       label = "p.signif")
  Plist[[i]] <- P
}

pdf("./5.1.1.epithelial_NMF/process_k3/boxplot_MPScore_AUC_diff_Bulk.pdf", width=16, height=7)
grid.arrange(grobs=Plist, ncol=5)
dev.off()

#### MP1 boxplot --------------------------------------------

MP1_list <- list(GSE183019 = data_GSE183019, GSE201284 = data_GSE201284, 
                 GSE141551=data_GSE141551, GSE153352=data_GSE153352, GSE157547=data_GSE157547, 
                 GSE88808 = data_GSE88808, GSE62116 = data_GSE62116, GSE21034 = data_GSE21034)
for (i in 1:length(MP1_list)) {
  MP1_list[[i]] <- MP1_list[[i]][MP1_list[[i]]$MP == "MP1",]
  MP1_list[[i]]$Type[which(MP1_list[[i]]$Type == "EOPC")] <- paste0("EOPC", " (", table(MP1_list[[i]]$Type)[1], ")")
  MP1_list[[i]]$Type[which(MP1_list[[i]]$Type == "LOPC")] <- paste0("LOPC", " (", table(MP1_list[[i]]$Type)[2], ")")
}

Plist <- list()
for (i in 1:length(MP1_list)) {
  P = ggplot(MP1_list[[i]], aes(x = Type, y = value, fill = Type)) +
    scale_fill_manual(values = c("#053E7AAA","#A52A2AAA")) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(names(MP1_list)[i]) +
    ylab("AUCell score") +
    xlab("") + 
    stat_compare_means(method="wilcox.test",
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "")),
                       label = "p.signif",
                       label.x.npc = "center") +
    NoLegend()
  Plist[[i]] <- P
}

pdf("./5.1.1.epithelial_NMF/process_k3/boxplot_MP1_AUC_diff_Bulk.pdf", width=8, height=7)
grid.arrange(grobs=Plist, ncol=4)
dev.off()

#### logistic regression and meta-analysis ---------------------------------------

MP1_logit_list <- list(TCGA = mtx_TCGA_logit, GSE183019 = mtx_GSE183019_logit, GSE201284 = mtx_GSE201284_logit, 
                       GSE141551=mtx_GSE141551_logit, GSE153352=mtx_GSE153352_logit, GSE157547=mtx_GSE157547_logit, 
                       GSE88808 = mtx_GSE88808_logit, GSE62116 = mtx_GSE62116_logit, GSE21034 = mtx_GSE21034_logit)

result <- c()
for (i in 1:length(MP1_logit_list)) {
  data <- MP1_logit_list[[i]]
  data$MP1 <- scale(data$MP1)
  logit <- glm(Type ~ MP1, data = data, family = binomial(link = "logit"))
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

pdf("./5.1.1.epithelial_NMF/process_k3/bulk_MP1_meta.pdf", width = 10, height = 7)
forest(meta, text.random = "", text.w.random = "", print.subgroup.name = F, lwd = 3, hetlab = "",
       test.subgroup = F, col.study = "#556B2F", leftcols = "studlab", leftlabs = "Dataset",
       print.tau2 = F)
dev.off()


### difference of MP ssGSEA score between EO and LO in Bulk ----------------------------------

library(reshape2)
library(GSVA)

#### TCGA ----------------------------

expMtx <- read.table("./0.TCGA/mRNA_Matrix_TPM.txt", row.names = 1, header = T, sep = "\t")
expMtx <- as.matrix(expMtx)
cells_ssGSEA_TCGA <- gsva(expMtx, Markers, method = "ssgsea", parallel.sz = 20)
saveRDS(cells_ssGSEA_TCGA, "./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_TCGA.Rds")

cells_ssGSEA_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_TCGA.Rds")
ssGSEA_mtx_t <- as.data.frame(t(cells_ssGSEA_TCGA))
clinical <- read.table("./0.TCGA/Clinical.txt", row.names = 1, header = T, sep = "\t")
age <- as.data.frame(clinical[,1])
colnames(age) <- "age"
rownames(age) <- gsub("_", ".", rownames(clinical))
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, ssGSEA_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_TCGA <- melt(mtx)
colnames(data_TCGA)[2] <- "MP"

#### GSE141551 ---------------------------------

exp <- read.table("./0.TCGA/GSE141551/MergeExpro_contrib1-GPL14951.txt", row.names = 1, header = T, sep = "\t")
exp <- as.matrix(exp)
cells_ssGSEA_GSE141551 <- gsva(exp, Markers, method = "ssgsea", parallel.sz = 20)
saveRDS(cells_ssGSEA_GSE141551, "./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE141551.Rds")

cells_ssGSEA_GSE141551 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE141551.Rds")
ssGSEA_mtx_t <- as.data.frame(t(cells_ssGSEA_GSE141551))
clinical <- read.table("./0.TCGA/GSE141551/SampleInfo_contrib1-GPL14951.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(agegroup = clinical[,"X1.AGE_GROUP"])
rownames(age) <- gsub("_", ".", rownames(clinical))
age$ageup <- as.numeric(substr(age$agegroup,4,5))
age$agelow <- as.numeric(substr(age$agegroup,1,2))
age$avg <- apply(age[,c(2,3)], 1, median)
age$Type[age$ageup < 55] <- "EOPC"
age$Type[age$agelow >= 55] <- "LOPC"
mtx <- merge(age, ssGSEA_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-c(1:4)]
data_GSE141551 <- melt(mtx)
colnames(data_GSE141551)[2] <- "MP"

#### GSE157547 ----------------------------------------

exp <- read.table("./0.TCGA/GSE157547/MergeExpro_contrib1-GPL5175.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_ssGSEA_GSE157547 <- gsva(exp, Markers, method = "ssgsea", parallel.sz = 20)
saveRDS(cells_ssGSEA_GSE157547, "./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE157547.Rds")

cells_ssGSEA_GSE157547 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE157547.Rds")
ssGSEA_mtx_t <- as.data.frame(t(cells_ssGSEA_GSE157547))
clinical <- read.table("./0.TCGA/GSE157547/SampleInfo_contrib1-GPL5175.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, ssGSEA_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE157547 <- melt(mtx)
colnames(data_GSE157547)[2] <- "MP"

#### GSE153352 --------------------------------------------------

exp <- read.table("./0.TCGA/GSE153352/MergeExpro_contrib1-GPL5175.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_ssGSEA_GSE153352 <- gsva(exp, Markers, method = "ssgsea", parallel.sz = 20)
saveRDS(cells_ssGSEA_GSE153352, "./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE153352.Rds")

cells_ssGSEA_GSE153352 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE153352.Rds")
ssGSEA_mtx_t <- as.data.frame(t(cells_ssGSEA_GSE153352))
clinical <- read.table("./0.TCGA/GSE153352/SampleInfo_contrib1-GPL5175.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, ssGSEA_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE153352 <- melt(mtx)
colnames(data_GSE153352)[2] <- "MP"

#### GSE88808 -----------------------------------------

exp <- read.table("./0.TCGA/GSE88808/MergeExpro_contrib1-GPL22571.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_ssGSEA_GSE88808 <- gsva(exp, Markers, method = "ssgsea", parallel.sz = 20)
saveRDS(cells_ssGSEA_GSE88808, "./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE88808.Rds")

cells_ssGSEA_GSE88808 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE88808.Rds")
ssGSEA_mtx_t <- as.data.frame(t(cells_ssGSEA_GSE88808))
clinical <- read.table("./0.TCGA/GSE88808/SampleInfo_contrib1-GPL22571.txt", row.names = 1, header = T, sep = "\t")
clinical <- clinical[clinical$"X1.TISSUE" == "tumor",]
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, ssGSEA_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE88808 <- melt(mtx)
colnames(data_GSE88808)[2] <- "MP"

#### GSE62116 ----------------------------------------

exp <- read.table("./0.TCGA/GSE62116/MergeExpro_contrib1-GPL5188.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_ssGSEA_GSE62116 <- gsva(exp, Markers, method = "ssgsea", parallel.sz = 20)
saveRDS(cells_ssGSEA_GSE62116, "./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE62116.Rds")

cells_ssGSEA_GSE62116 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE62116.Rds")
ssGSEA_mtx_t <- as.data.frame(t(cells_ssGSEA_GSE62116))
clinical <- read.table("./0.TCGA/GSE62116/SampleInfo_contrib1-GPL5188.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"X1.AGE"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, ssGSEA_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE62116 <- melt(mtx)
colnames(data_GSE62116)[2] <- "MP"

#### GSE21034 ----------------------------------------

exp <- read.table("./0.TCGA/GSE21034/MergeExpro_contrib1-GPL5188.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_ssGSEA_GSE21034 <- gsva(exp, Markers, method = "ssgsea", parallel.sz = 20)
saveRDS(cells_ssGSEA_GSE21034, "./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE21034.Rds")

cells_ssGSEA_GSE21034 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE21034.Rds")
ssGSEA_mtx_t <- as.data.frame(t(cells_ssGSEA_GSE21034))
clinical <- read.table("./0.TCGA/GSE21034/GSE21034_clinical.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"DxAge"])
rownames(age) <- rownames(clinical)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, ssGSEA_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE21034 <- melt(mtx)
colnames(data_GSE21034)[2] <- "MP"

#### GSE183019 ----------------------------------------

exp <- read.table("./0.TCGA/GSE183019/GSE183019_mRNA_processed_TPM.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_ssGSEA_GSE183019 <- gsva(exp, Markers, method = "ssgsea", parallel.sz = 20)
saveRDS(cells_ssGSEA_GSE183019, "./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE183019.Rds")

cells_ssGSEA_GSE183019 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE183019.Rds")
ssGSEA_mtx_t <- as.data.frame(t(cells_ssGSEA_GSE183019))
clinical <- read.table("./0.TCGA/GSE183019/clinical.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"Age"])
rownames(age) <- gsub("-", ".", clinical$Sample)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, ssGSEA_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE183019 <- melt(mtx)
colnames(data_GSE183019)[2] <- "MP"

#### GSE201284 ----------------------------------------

exp <- read.table("./0.TCGA/GSE201284/GSE201284_mRNA_processed_TPM.txt", header = T, sep = "\t", check.names = F)
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp <- avereps(exp)
cells_ssGSEA_GSE201284 <- gsva(exp, Markers, method = "ssgsea", parallel.sz = 20)
saveRDS(cells_ssGSEA_GSE201284, "./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE201284.Rds")

cells_ssGSEA_GSE201284 <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_ssGSEA_GSE201284.Rds")
ssGSEA_mtx_t <- as.data.frame(t(cells_ssGSEA_GSE201284))
clinical <- read.table("./0.TCGA/GSE201284/clinical.txt", row.names = 1, header = T, sep = "\t")
age <- data.frame(age = clinical[,"Age"])
rownames(age) <- gsub("-", ".", clinical$Sample)
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"
mtx <- merge(age, ssGSEA_mtx_t, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[,-1]
data_GSE201284 <- melt(mtx)
colnames(data_GSE201284)[2] <- "MP"

#### boxplot ----------------------------------------

library(ggpubr)

data_list <- list(TCGA=data_TCGA, GSE183019 = data_GSE183019, GSE201284 = data_GSE201284, 
                  GSE141551=data_GSE141551, GSE153352=data_GSE153352, GSE157547=data_GSE157547, 
                  GSE88808 = data_GSE88808, GSE62116 = data_GSE62116, GSE21034 = data_GSE21034)
Plist <- list()
for (i in 1:length(data_list)) {
  P = ggplot(data_list[[i]], aes(x = MP, y = value, fill = Type)) +
    scale_fill_manual(values = c("#053E7AAA","#A52A2AAA")) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          legend.text=element_text(colour="black", size=8),
          legend.title=element_text(colour="black", size=10),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(names(data_list)[i]) +
    ylab("ssGSEA score") +
    xlab("") + 
    stat_compare_means(aes(group=Type),
                       method="wilcox.test",
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                       label = "p.signif")
  Plist[[i]] <- P
}

pdf("./5.1.1.epithelial_NMF/process_k3/boxplot_MPScore_ssGSEA_diff_Bulk.pdf", width=16, height=7)
grid.arrange(grobs=Plist, ncol=5)
dev.off()

#### MP1 boxplot --------------------------------------------

MP1_list <- list(GSE183019 = data_GSE183019, GSE201284 = data_GSE201284, 
                 GSE141551=data_GSE141551, GSE153352=data_GSE153352, GSE157547=data_GSE157547, 
                 GSE88808 = data_GSE88808, GSE62116 = data_GSE62116, GSE21034 = data_GSE21034)
for (i in 1:length(MP1_list)) {
  MP1_list[[i]] <- MP1_list[[i]][MP1_list[[i]]$MP == "MP1",]
  MP1_list[[i]]$Type[which(MP1_list[[i]]$Type == "EOPC")] <- paste0("EOPC", " (", table(MP1_list[[i]]$Type)[1], ")")
  MP1_list[[i]]$Type[which(MP1_list[[i]]$Type == "LOPC")] <- paste0("LOPC", " (", table(MP1_list[[i]]$Type)[2], ")")
}

Plist <- list()
for (i in 1:length(MP1_list)) {
  P = ggplot(MP1_list[[i]], aes(x = Type, y = value, fill = Type)) +
    scale_fill_manual(values = c("#053E7AAA","#A52A2AAA")) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(names(MP1_list)[i]) +
    ylab("AUCell score") +
    xlab("") + 
    stat_compare_means(method="wilcox.test",
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "")),
                       label = "p.signif",
                       label.x.npc = "center") +
    NoLegend()
  Plist[[i]] <- P
}

pdf("./5.1.1.epithelial_NMF/process_k3/boxplot_MP1_ssGSEA_diff_Bulk.pdf", width=8, height=7)
grid.arrange(grobs=Plist, ncol=4)
dev.off()

### plot MP AUCell and On/Off on umap plot ----------------------------------

Epi_MP_AUCell <- CreateSeuratObject(AUCell_mtx,
                                    project = "Epi_MP", 
                                    names.field = 1,
                                    names.delim = "_")
Epi <- readRDS("./5.1.epithelial_harmony_TryTwo/Epi_clustered.Rds")
Epi_MP_AUCell@meta.data <- Epi@meta.data
Epi_MP_AUCell@reductions <- Epi@reductions

threshould <- data.frame(MP = c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6"),
                         cutoff = c(0.4, 0.81, 0.51, 0.58, 0.38, 0.54))

for (i in rownames(AUCell_mtx)) {
  Epi_MP_AUCell$proportion <- ifelse(Epi_MP_AUCell@assays$RNA@data[i,] > median(Epi_MP_AUCell@assays$RNA@data[i,]), "On", "Off")
  colnames(Epi_MP_AUCell@meta.data)[grep("proportion", colnames(Epi_MP_AUCell@meta.data))] <- paste0(i, "_OnOff")
}

DimPlot(Epi_MP_AUCell, group.by = "MP1_OnOff", split.by = "orig.ident")

Epi <- readRDS("./5.1.epithelial_harmony_TryTwo/Epi_clustered.Rds")
cells_AUC <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC.Rds")
AUCell_mtx <- as.matrix(as.data.frame(cells_AUC@assays@data$AUC))
Epi[["MP_AUCell"]] <- CreateAssayObject(counts = AUCell_mtx)
DefaultAssay(Epi) <- "MP_AUCell"

FeaturePlot(Epi, "MP1")

### MP hub genes ----------------------------------------------

files <- paste0("EOPC", c(1:4), "_program.Zscore.txt")
files <- c(files, paste0("LOPC", c(1:6), "_program.Zscore.txt"))

topgenes <- read.table("./5.1.1.epithelial_NMF/res2.1.3/program_topngene.txt", row.names = 1, header = T, sep = "\t")

MP1 <- c("EOPC4_1","LOPC3_2","EOPC3_2","LOPC1_1","LOPC2_2","LOPC4_2","EOPC3_1","LOPC1_7","LOPC6_2","EOPC1_13","EOPC2_2","LOPC4_9")
MP2 <- c("EOPC3_4","EOPC4_4","EOPC1_7","LOPC3_3","LOPC4_7","LOPC5_1","LOPC2_3","LOPC3_8")
MP3 <- c("EOPC3_3","EOPC4_5","LOPC3_10","EOPC4_3","LOPC3_6")
MP4 <- c("EOPC1_2","LOPC5_2","LOPC1_4","LOPC6_5")
MP5 <- c("EOPC1_11","LOPC5_6","EOPC2_7","EOPC4_6","LOPC3_7","LOPC1_2","LOPC4_5")
MP6 <- c("EOPC2_9","LOPC4_1","LOPC5_5","LOPC1_5","LOPC4_8","EOPC3_7","LOPC1_6")

MPs <- list(MP1 = MP1, MP2 = MP2, MP3 = MP3, MP4 = MP4,MP5 = MP5, MP6 = MP6)

for (i in files) {
  rt <- read.table(paste0("./5.1.1.epithelial_NMF/res2.1.3/", i), header = T, sep = "\t")
  assign(substr(i, 1, 5), rt)
}

ProgramList <- list(EOPC1 = EOPC1, EOPC2 = EOPC2, EOPC3 = EOPC3, EOPC4 = EOPC4, 
                    LOPC1 = LOPC1, LOPC2 = LOPC2, LOPC3 = LOPC3, LOPC4 = LOPC4, 
                    LOPC5 = LOPC5, LOPC6 = LOPC6)

GeneScoreList <- list()
for (i in 1:length(MPs)) {
  list_tmp <- list()
  for (j in 1:length(MPs[[i]])) {
    dt <- ProgramList[[substr(MPs[[i]][j],1,5)]]
    Gene_score <- dt[, c(grep(MPs[[i]][j], colnames(dt)),grep("gene", colnames(dt)))]
    Gene_score <- Gene_score[order(Gene_score[,1], decreasing = T),][1:50,]
    list_tmp[[j]] <- Gene_score
    names(list_tmp)[j] <- MPs[[i]][j]
  }
  GeneScoreList[[i]] <- list_tmp
  names(GeneScoreList)[i] <- names(MPs)[i]
}

saveRDS(GeneScoreList, "./5.1.1.epithelial_NMF/process_k3/MP_GeneScoreList_top50.Rds")

GeneScoreList <- list()
for (i in 1:length(MPs)) {
  list_tmp <- list()
  for (j in 1:length(MPs[[i]])) {
    dt <- ProgramList[[substr(MPs[[i]][j],1,5)]]
    Gene_score <- dt[, c(grep(MPs[[i]][j], colnames(dt)),grep("gene", colnames(dt)))]
    list_tmp[[j]] <- Gene_score
    names(list_tmp)[j] <- MPs[[i]][j]
  }
  GeneScoreList[[i]] <- list_tmp
  names(GeneScoreList)[i] <- names(MPs)[i]
}

saveRDS(GeneScoreList, "./5.1.1.epithelial_NMF/process_k3/MP_GeneScoreList.Rds")

topgenes <- read.table("./5.1.1.epithelial_NMF/res2.1.3/program_topngene.txt", header = T, sep = "\t")
MP1 <- topgenes[, c("EOPC4_1","LOPC3_2","EOPC3_2","LOPC1_1","LOPC2_2","LOPC4_2","EOPC3_1","LOPC1_7","LOPC6_2","EOPC1_13","EOPC2_2","LOPC4_9")]
MP2 <- topgenes[, c("EOPC3_4","EOPC4_4","EOPC1_7","LOPC3_3","LOPC4_7","LOPC5_1","LOPC2_3","LOPC3_8")]
MP3 <- topgenes[, c("EOPC3_3","EOPC4_5","LOPC3_10","EOPC4_3","LOPC3_6")]
MP4 <- topgenes[, c("EOPC1_2","LOPC5_2","LOPC1_4","LOPC6_5")]
MP5 <- topgenes[, c("EOPC1_11","LOPC5_6","EOPC2_7","EOPC4_6","LOPC3_7","LOPC1_2","LOPC4_5")]
MP6 <- topgenes[, c("EOPC2_9","LOPC4_1","LOPC5_5","LOPC1_5","LOPC4_8","EOPC3_7","LOPC1_6")]
MPlist <- list(MP1,MP2,MP3,MP4,MP5,MP6)

Markers <- readRDS("./5.1.1.epithelial_NMF/process_k3/MP_markers.Rds")

MPGeneScoreAvg <- list()
for (i in 1:length(GeneScoreList)) {
  ScoreAvg <- c()
  for (j in 1:length(Markers[[i]])) {
    sum <- 0
    for (n in 1:length(GeneScoreList[[i]])) {
      if (length(which(GeneScoreList[[i]][[n]]$gene == Markers[[i]][j])) != 0) {
        sum <- sum + GeneScoreList[[i]][[n]][which(GeneScoreList[[i]][[n]]$gene == Markers[[i]][j]),1]
      }
    }
    avg <- sum / length(GeneScoreList[[i]])
    tmp <- data.frame(Gene = Markers[[i]][j], ZScoreAvg = avg)
    ScoreAvg <- rbind(ScoreAvg, tmp)
  }
  ScoreAvg <- ScoreAvg[order(ScoreAvg$ZScoreAvg, decreasing = T),]
  MPGeneScoreAvg[[i]] <- ScoreAvg
  names(MPGeneScoreAvg)[i] <- names(GeneScoreList)[i]
}

saveRDS(MPGeneScoreAvg, "./5.1.1.epithelial_NMF/process_k3/MP_Gene_Score_Avg.Rds")

MPGeneScore <- c()
for (i in 1:length(MPGeneScoreAvg)) {
  tmp <- data.frame(MP = names(MPGeneScoreAvg)[i], MPGeneScoreAvg[[i]])
  MPGeneScore <- rbind(MPGeneScore,tmp)
}

### MP1 top 100 genes enrichment ------------------------------------

hsets <- read.gmt("./5.1.1.epithelial_NMF/hallmark_cancersea.gmt")
hsets$term <- gsub("_", " ", hsets$term)
hsets$term <- str_to_title(hsets$term)

library(msigdbr)
msigdbr_show_species()

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
hallmark$gs_name <- gsub("HALLMARK_", "", hallmark$gs_name)
hallmark$gs_name <- gsub("_", " ", hallmark$gs_name)
hallmark$gs_name <- str_to_title(hallmark$gs_name)

BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
BIOCARTA$gs_name <- gsub("BIOCARTA_", "", BIOCARTA$gs_name)
BIOCARTA$gs_name <- gsub("_", " ", BIOCARTA$gs_name)
BIOCARTA$gs_name <- str_to_title(BIOCARTA$gs_name)

PID <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
PID$gs_name <- gsub("PID_", "", PID$gs_name)
PID$gs_name <- gsub("_PATHWAY", "", PID$gs_name)
PID$gs_name <- gsub("_", " ", PID$gs_name)
PID$gs_name <- str_to_title(PID$gs_name)

MP1 <- MPGeneScoreAvg[[1]]$Gene
MP1_ENTRZ <- AnnotationDbi::select(org.Hs.eg.db, keys=MP1, keytype="SYMBOL", columns=c("ENTREZID"))
ENTRAID <- as.character(MP1_ENTRZ$ENTREZID)

enrich_cancersea <- enricher(MP1, TERM2GENE = hsets, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
enrich_kegg <- enrichKEGG(ENTRAID, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 5)
enrich_BIOCARTA <- enricher(MP1, TERM2GENE = BIOCARTA, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
enrich_PID <- enricher(MP1, TERM2GENE = PID, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)

cancersea <- enrich_cancersea@result
cancersea$Type <- "cancersea"
KEGG <- enrich_kegg@result
KEGG$Type <- "KEGG"
KEGG$Description <- str_to_title(KEGG$Description)
BIOCARTA <- enrich_BIOCARTA@result
BIOCARTA$Type <- "BIOCARTA"
PID <- enrich_PID@result
PID$Type <- "PID"

enrich_result <- rbind(cancersea[,c(10,2,5)], KEGG[,c(10,2,5)],
                       BIOCARTA[,c(10,2,5)], PID[,c(10,2,5)])
enrich_result <- enrich_result[which(enrich_result$pvalue < 0.05),]
enrich_result$logP <- -log10(enrich_result$pvalue)
enrich_result$Description[which(duplicated(enrich_result$Description))] <-
  paste0(enrich_result$Description[which(duplicated(enrich_result$Description))], "_2")

enrich_result$Description <- factor(enrich_result$Description, levels = unique(enrich_result$Description))
enrich_result$Type <- factor(enrich_result$Type, levels = unique(enrich_result$Type))

P <- ggplot(data = enrich_result, aes(x = Description, y = logP, fill = Type)) +
  geom_bar(stat='identity') +
  ggtitle("MP1 enrichment") +
  xlab("") +
  ylab("-log10(P-value)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 10),
        axis.text.y=element_text(size=10, colour="black"), 
        axis.title.y=element_text(size = 12),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#CD7054", "#2F4F4F", "#7A378B", "#008B8B"))

pdf("./5.1.1.epithelial_NMF/process_k3/MP1_enrichment_barplot.pdf", width = 12, height = 6)
P
dev.off()

#install.packages("gg.gap")
library(gg.gap)

pdf("./5.1.1.epithelial_NMF/process_k3/MP1_enrichment_barplot_gap.pdf", width = 12, height = 7)
gg.gap(plot=P,
       segments = c(10,20),
       ylim = c(0,25),
       tick_width = c(5,5), 
       rel_heights = c(1,0,0.25))
dev.off()

### cytoTRACE (correlation with MP1) -----------------------------------

Epi <- readRDS("./5.1.epithelial_harmony_TryTwo/Epi_clustered.Rds")
exp <- as.matrix(Epi@assays$RNA@counts)
SampleType <- Epi@meta.data$SampleType
names(SampleType) <- rownames(Epi@meta.data)

Epi_CytoTrace_result <- CytoTRACE(exp, enableFast = TRUE, ncores = 20, subsamplesize = 1000)
saveRDS(Epi_CytoTrace_result, "./5.1.epithelial_harmony_TryTwo/CytoTRACE/Epi_CytoTrace_result.Rds")

exp <- exp[-grep("^RP[SL]", rownames(exp)),]
Epi_CytoTrace_result <- CytoTRACE(exp, enableFast = TRUE, ncores = 20, subsamplesize = 1000)
saveRDS(Epi_CytoTrace_result, "./5.1.epithelial_harmony_TryTwo/CytoTRACE/Epi_noRib_CytoTrace_result.Rds")

Sample <- Epi@meta.data$orig.ident
names(Sample) <- rownames(Epi@meta.data)
Epi_CytoTrace_result <- CytoTRACE(exp, enableFast = TRUE, batch = Sample, ncores = 20, subsamplesize = 1000)
saveRDS(Epi_CytoTrace_result, "./5.1.epithelial_harmony_TryTwo/CytoTRACE/Epi_noRib_debatch_CytoTrace_result.Rds")

emb <- Epi@reductions$umap@cell.embeddings
plotCytoTRACE(Epi_CytoTrace_result, emb = emb, outputDir = "./5.1.1.epithelial_NMF/process_k3/CytoTRACE/noRib_debatch_")
plotCytoGenes(Epi_CytoTrace_result, numOfGenes = 10, outputDir = "./5.1.1.epithelial_NMF/process_k3/CytoTRACE/noRib_debatch_", colors = c("#A52A2A", "#053E7A"))

emb <- Epi@reductions$umap@cell.embeddings

# MP score as assay

cells_AUC <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC.Rds")
AUCell_mtx <- as.matrix(cells_AUC@assays@data$AUC)
Epi[["MP_AUCell"]] <- CreateAssayObject(counts = AUCell_mtx)
DefaultAssay(Epi) <- "MP_AUCell"

threshould <- data.frame(MP = c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6"),
                         cutoff = c(0.4, 0.81, 0.51, 0.58, 0.38, 0.54))

for (i in rownames(AUCell_mtx)) {
  Epi$proportion <- ifelse(Epi@assays$MP_AUCell@data[i,] > threshould[threshould$MP == i,2], "On", "Off")
  colnames(Epi@meta.data)[grep("proportion", colnames(Epi@meta.data))] <- paste0(i, "_OnOff")
}

DimPlot(Epi, group.by="SampleType", cols = c("#A52A2A", "#053E7A"), pt.size=0.01, shuffle=T)

P1 <- DimPlot(Epi, group.by="MP1_OnOff", cols = c("#A52A2A", "#053E7A"), pt.size = 0.5, shuffle = T)

for (i in rownames(AUCell_mtx)) {
  Epi$proportion <- ifelse(Epi@assays$MP_AUCell@data[i,] > median(Epi@assays$MP_AUCell@data[i,]), "High", "Low")
  colnames(Epi@meta.data)[grep("proportion", colnames(Epi@meta.data))] <- paste0(i, "_HL")
}

saveRDS(Epi, "./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")

P2 <- DimPlot(Epi, group.by="MP1_HL", cols = c("#9E0142", "#5E4FA2"), pt.size = 0.5, shuffle = T)

pdf("./5.1.1.epithelial_NMF/process_k3/CytoTRACE/DimPlot_MP1_pct_CytoTRACE.pdf", width = 14, height = 7)
P1 + P2
dev.off()

pdf("./5.1.1.epithelial_NMF/process_k3/CytoTRACE/DimPlot_MP1_OnOff_CytoTRACE.pdf", width = 5, height = 5)
P1 + NoLegend() + theme(plot.title = element_blank())
dev.off()

# CytoTRACE asassay

CytoTRACE_mtx <- as.data.frame(Epi_CytoTrace_result$CytoTRACE)
colnames(CytoTRACE_mtx) <- "CytoTRACE"

MP1_OnOff <- data.frame(MP1_OnOff = Epi@meta.data$MP1_OnOff)
rownames(MP1_OnOff) <- rownames(Epi@meta.data)
MP1_HL <- data.frame(MP1_HL = Epi@meta.data$MP1_HL)
rownames(MP1_HL) <- rownames(Epi@meta.data)

# boxplot_OnOff

mtx <- merge(CytoTRACE_mtx, MP1_OnOff, by = 0) %>% column_to_rownames("Row.names")
data <- melt(mtx)
data$MP1_OnOff <- factor(data$MP1_OnOff, levels = c("On", "Off"))

P = ggplot(data, aes(x = MP1_OnOff, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("CytoTRACE") +
  ylab("Predicted order") +
  xlab("") + 
  stat_compare_means(method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif",
                     label.x.npc = "center") +
  NoLegend()

pdf("./5.1.1.epithelial_NMF/process_k3/CytoTRACE/boxplot_CytoTRACE_MP1_On_Off.pdf", width = 5)
P
dev.off()

# boxplot_HL

mtx <- merge(CytoTRACE_mtx, MP1_HL, by = 0) %>% column_to_rownames("Row.names")
data <- melt(mtx)

P = ggplot(data, aes(x = MP1_HL, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#A52A2AAA", "#053E7AAA")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("CytoTRACE") +
  ylab("Predicted order") +
  xlab("") + 
  stat_compare_means(method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif",
                     label.x.npc = "center") +
  NoLegend()

pdf("./5.1.1.epithelial_NMF/process_k3/CytoTRACE/boxplot_CytoTRACE_MP1_H_L.pdf", width = 5)
P
dev.off()

# correlation of MP1 and CytoTRACE score

MP1 <- Epi@assays$MP_AUCell@data["MP1",]
CytoTRACE <- Epi_CytoTrace_result$CytoTRACE

df <- as.data.frame(t(rbind(MP1, CytoTRACE)))
corP <- cor.test(MP1, CytoTRACE, method = "pearson")
corS <- cor.test(MP1, CytoTRACE, method = "spearman", exact = F)
cor <- corP$estimate
pValue <- corP$p.value
p1 <- ggplot(df, aes(x=MP1, y=CytoTRACE)) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
  theme_bw()+
  stat_cor(method = 'pearson', aes(x = MP1, y = CytoTRACE)) +
  theme(axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) 
p2=ggMarginal(p1, type = "density", xparams = list(fill = "#9E0142AA"),yparams = list(fill = "#5E4FA2AA"))

pdf("./5.1.1.epithelial_NMF/process_k3/CytoTRACE/MP1_CytoTRACE_corPlot_pearson.pdf")
p2
dev.off()

### stemness score based on signature (not favorable) --------------------------------------------------

stem <- read.table("./stemcell.txt", header = T)
stem <- list(stemness = stem$Genes)

expMtx <- as.matrix(Epi@assays$RNA@data)
cells_rankings <- AUCell_buildRankings(expMtx, nCores=10, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(stem, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.1.epithelial_harmony_TryTwo/Stem_cell_score/StemScore_AUC.Rds")

# diff between MP1 On/Off

MP1_OnOff <- data.frame(MP1_OnOff = Epi@meta.data$MP1_OnOff)
rownames(MP1_OnOff) <- rownames(Epi@meta.data)
MP1_HL <- data.frame(MP1_HL = Epi@meta.data$MP1_HL)
rownames(MP1_HL) <- rownames(Epi@meta.data)

stem_cells_AUC <- readRDS("./5.1.epithelial_harmony_TryTwo/Stem_cell_score/StemScore_AUC.Rds")
stem_mtx <- as.data.frame(t(as.data.frame(stem_cells_AUC@assays@data$AUC)))

mtx <- merge(stem_mtx, MP1_HL, by = 0) %>% column_to_rownames("Row.names")
data <- melt(mtx)

P = ggplot(data, aes(x = MP1_HL, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#A52A2AAA", "#053E7AAA")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Stemness (Miranda et al.)") +
  ylab("AUCell score") +
  xlab("") + 
  stat_compare_means(method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif",
                     label.x.npc = "center") +
  NoLegend()

pdf("./5.1.1.epithelial_NMF/process_k3/CytoTRACE/boxplot_CytoTRACE_MP1_On_Off.pdf", width = 5)
P
dev.off()

### stemness score based on OCLR mRNAsi --------------------------------------------------

source("./PCBC/predict.mRNAsi.R")
exp <- as.matrix(Epi@assays$RNA@data)
mRNAsi_score <- predict.mRNAsi(exp)
saveRDS(mRNAsi_score, "./5.1.epithelial_harmony_TryTwo/mRNAsi/mRNAsi_score.Rds")

# diff between MP1 On/Off

MP1_OnOff <- data.frame(MP1_OnOff = Epi@meta.data$MP1_OnOff)
rownames(MP1_OnOff) <- rownames(Epi@meta.data)
MP1_HL <- data.frame(MP1_HL = Epi@meta.data$MP1_HL)
rownames(MP1_HL) <- rownames(Epi@meta.data)

mRNAsi_score <- readRDS("./5.1.epithelial_harmony_TryTwo/mRNAsi/mRNAsi_score.Rds")
mRNAsi_mtx <- (as.data.frame(mRNAsi_score))

mtx <- merge(mRNAsi_mtx, MP1_HL, by = 0) %>% column_to_rownames("Row.names")
data <- melt(mtx)

P = ggplot(data, aes(x = MP1_HL, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#A52A2AAA", "#053E7AAA")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("mRNAsi (Malta et al.)") +
  ylab("Stemness Index") +
  xlab("") + 
  stat_compare_means(method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif",
                     label.x.npc = "center") +
  NoLegend()

pdf("./5.1.1.epithelial_NMF/process_k3/mRNAsi/boxplot_mRNAsi_MP1_H_L.pdf", width = 5)
P
dev.off()

mtx <- merge(mRNAsi_mtx, MP1_OnOff, by = 0) %>% column_to_rownames("Row.names")
data <- melt(mtx)

P = ggplot(data, aes(x = MP1_OnOff, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#A52A2AAA", "#053E7AAA")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("mRNAsi (Malta et al.)") +
  ylab("Stemness Index") +
  xlab("") + 
  stat_compare_means(method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif",
                     label.x.npc = "center") +
  NoLegend()

pdf("./5.1.1.epithelial_NMF/process_k3/mRNAsi/boxplot_mRNAsi_MP1_On_Off.pdf", width = 5)
P
dev.off()

# SampleType - MP1_OnOff

mRNAsi_score <- readRDS("./5.1.epithelial_harmony_TryTwo/mRNAsi/mRNAsi_score.Rds")
mRNAsi_mtx <- (as.data.frame(mRNAsi_score))
Epi$SampleTypeMP1 <- paste(Epi$SampleType,Epi$MP1_OnOff,sep="_")
MP1_OnOff <- Epi$SampleTypeMP1 %>% as.data.frame()
colnames(MP1_OnOff) <- "Group_Onoff"
mtx <- merge(mRNAsi_mtx, MP1_OnOff, by = 0) %>% column_to_rownames("Row.names")
mtx$Group_Onoff <- factor(mtx$Group_Onoff, levels = c("EOPC_On", "LOPC_On", "EOPC_Off", "LOPC_Off"))
mtx <- mtx[order(mtx$Group_Onoff),]
data <- melt(mtx)

P = ggplot(data, aes(x = Group_Onoff, y = value, fill = Group_Onoff)) +
  scale_fill_manual(values = c("#E46B4FCC", "#8B8B00CC", "#B8860BCC", "#483D8BCC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("mRNAsi (Malta et al.)") +
  ylab("Stemness Index") +
  xlab("") +
  NoLegend()

pdf("./5.1.1.epithelial_NMF/process_k3/mRNAsi/boxplot_mRNAsi_SampleType_MP1_On_Off.pdf", width = 5)
P
dev.off()

### correlation with clinical characteristics (TCGA) ---------------------------------------

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_TCGA.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)
AUCell_mtx_t <- as.data.frame(t(AUCell_mtx))
clinical <- read.table("./0.TCGA/clinical_simplify.txt", header = T, row.names = 1, sep = "\t")
rownames(clinical) <- gsub("-", ".", rownames(clinical))
clinical[clinical == ""] <- NA

#M
clinical$clinical_M[grep("M1", clinical$clinical_M)] <- "M1"

#GS
clinical$gleason_score <- factor(clinical$gleason_score, 
                                 levels = c("6","7","8","9","10"))

mtx <- merge(AUCell_mtx_t, clinical, by = 0) %>% column_to_rownames("Row.names")
mtx <- mtx[, c("MP1", "BCR", "clinical_M", "gleason_score", "pathologic_N", "pathologic_T")]

#wilcoxon.test
pairwise.wilcox.test(mtx$MP1, mtx$clinical_M, p.adjust.method = "none")
pairwise.wilcox.test(mtx$MP1, mtx$gleason_score, p.adjust.method = "none")
pairwise.wilcox.test(mtx$MP1, mtx$pathologic_N, p.adjust.method = "none")
pairwise.wilcox.test(mtx$MP1, mtx$pathologic_T, p.adjust.method = "none")
pairwise.wilcox.test(mtx$MP1, mtx$BCR, p.adjust.method = "none")

#barplot

color <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
           "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
           "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")
color <- paste0(color, "DD")

library(reshape2)
library(ggpubr)

Plist <- list()
for (i in 2:ncol(mtx)) {
  data <- mtx[, c(1,i)]
  NAs <- which(is.na(data[,2]))
  if(length(NAs) > 0) {data <- data[-NAs,]} else {data <- data}
  levels <- length(table(data[,2]))
  col <- color[1:levels]
  data <- melt(data)
  colnames(data) <- c("cli", "MP1", "value")
  
  P = ggplot(data, aes(x = cli, y = value, fill = cli)) +
    scale_fill_manual(values = col) +
    geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 60, hjust = 1, colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(colnames(mtx)[i]) +
    ylab("MP1 AUCell") +
    xlab("") + 
    stat_compare_means(method="wilcox.test",
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "")),
                       label = "p.signif",
                       label.x.npc = "center") +
    NoLegend()
  
  Plist[[i-1]] <- P
  names(Plist)[i-1] <- colnames(mtx)[i]
}

for (i in 1:length(Plist)) {
  ggsave(plot = Plist[[i]], 
         filename = paste0("5.1.1.epithelial_NMF/process_k3/MP1_CliBox_", names(Plist)[i], ".pdf"),
         height = 3,
         width = 2 + 0.5*1*length(table(mtx[,(i+1)])))
}

### scMetablism diff EO vs LO ------------------------------------

#### On/Off KEGG --------------------------------

PCSC_Metab <- readRDS("2.1.harmony/scMetabolism/PCSC_Metab_KEGG.Rds")

Epi <- readRDS("./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
type <- Epi$MP1_OnOff
type <- sort(type)

Met_mtx <- PCSC_Metab@assays$METABOLISM$score
Met_mtx <- Met_mtx[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("MP1_Off", table(type)[1]), rep("MP1_On", table(type)[2])), levels = c('MP1_Off', 'MP1_On'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(Met_mtx)
head(design)

# Tunor VS Normal
compare <- makeContrasts(MP1_On - MP1_Off, levels=design)
fit <- lmFit(Met_mtx, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# 排序
Diff <- Diff %>% arrange(t)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 20, wt = t), top_n(dat_plot, -10, wt = t))
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
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black') +
  ylim(c(-85,85))
ggsave("5.1.1.epithelial_NMF/process_k3/scMetabolism/MP1_OnOff_KEGG_Met_AUCell_diff_bar_MP1_top15.pdf", p_MP1, width = 8, height  = 8)

#### High/Low KEGG --------------------------------

PCSC_Metab <- readRDS("2.1.harmony/scMetabolism/PCSC_Metab_KEGG.Rds")

Epi <- readRDS("./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
type <- Epi$MP1_HL
type <- sort(type)

Met_mtx <- PCSC_Metab@assays$METABOLISM$score
Met_mtx <- Met_mtx[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("MP1_High", table(type)[1]), rep("MP1_Low", table(type)[2])), levels = c('MP1_High', 'MP1_Low'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(Met_mtx)
head(design)

# Tunor VS Normal
compare <- makeContrasts(MP1_High - MP1_Low, levels=design)
fit <- lmFit(Met_mtx, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# 排序
Diff <- Diff %>% arrange(t)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 20, wt = t), top_n(dat_plot, -10, wt = t))
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
  ylab('t value of AUCell score, MP1_High versus Low') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black') +
  ylim(c(-80,80))
ggsave("5.1.1.epithelial_NMF/process_k3/scMetabolism/MP1_HL_KEGG_Met_AUCell_diff_bar_MP1_top15.pdf", p_MP1, width = 8, height  = 8)

#### On/Off REACTOME --------------------------------

PCSC_Metab <- readRDS("2.1.harmony/scMetabolism/PCSC_Metab_REACTOME.Rds")

Epi <- readRDS("./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
type <- Epi$MP1_OnOff
type <- sort(type)

Met_mtx <- PCSC_Metab@assays$METABOLISM$score
Met_mtx <- Met_mtx[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("MP1_Off", table(type)[1]), rep("MP1_On", table(type)[2])), levels = c('MP1_Off', 'MP1_On'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(Met_mtx)
head(design)

# Tunor VS Normal
compare <- makeContrasts(MP1_On - MP1_Off, levels=design)
fit <- lmFit(Met_mtx, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# 排序
Diff <- Diff %>% arrange(t)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 20, wt = t), top_n(dat_plot, -10, wt = t))
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
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black') +
  ylim(c(-116,116))
ggsave("5.1.1.epithelial_NMF/process_k3/scMetabolism/MP1_OnOff_REACTOME_Met_AUCell_diff_bar_MP1_top15.pdf", p_MP1, width = 8, height  = 8)

#### High/Low REACTOME --------------------------------

PCSC_Metab <- readRDS("2.1.harmony/scMetabolism/PCSC_Metab_REACTOME.Rds")

Epi <- readRDS("./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
type <- Epi$MP1_HL
type <- sort(type)

Met_mtx <- PCSC_Metab@assays$METABOLISM$score
Met_mtx <- Met_mtx[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("MP1_High", table(type)[1]), rep("MP1_Low", table(type)[2])), levels = c('MP1_High', 'MP1_Low'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(Met_mtx)
head(design)

# Tunor VS Normal
compare <- makeContrasts(MP1_High - MP1_Low, levels=design)
fit <- lmFit(Met_mtx, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# 排序
Diff <- Diff %>% arrange(t)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 20, wt = t), top_n(dat_plot, -10, wt = t))
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
  ylab('t value of AUCell score, MP1_High versus Low') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black') +
  ylim(c(-115,115))
ggsave("5.1.1.epithelial_NMF/process_k3/scMetabolism/MP1_HL_REACTOME_Met_AUCell_diff_bar_MP1_top15.pdf", p_MP1, width = 8, height  = 8)


### scFEA -----------------------------------------------------

#### On/Off ------------------------------------------
scFEA <- read.csv("5.1.epithelial_harmony_TryTwo/scFEA/output/balance_20220830-210644.csv", row.names =1)
scFEA <- t(scFEA)
Epi <- readRDS("./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
Epi[["FEA"]] <- CreateAssayObject(counts = scFEA)


type <- Epi$MP1_OnOff
type <- sort(type)

FEA_mtx <- Epi@assays$FEA@counts
FEA_mtx <- FEA_mtx[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("MP1_Off", table(type)[1]), rep("MP1_On", table(type)[2])), levels = c('MP1_Off', 'MP1_On'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(FEA_mtx)
head(design)

# Tunor VS Normal
compare <- makeContrasts(MP1_On - MP1_Off, levels=design)
fit <- lmFit(FEA_mtx, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# 排序
Diff <- Diff %>% arrange(t)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 15, wt = t), top_n(dat_plot, -15, wt = t))
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
  ylab('t value of scFEA score, MP1_On versus Off') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black') +
  ylim(c(-85,85))
ggsave("5.1.1.epithelial_NMF/process_k3/scFEA/FEA_diff_bar_MP1_top15.pdf", p_MP1, width = 8, height  = 8)

#### High/Low ----------------------------------------------

scFEA <- read.csv("5.1.epithelial_harmony_TryTwo/scFEA/output/balance_20220830-210644.csv", row.names =1)
scFEA <- t(scFEA)
Epi <- readRDS("./5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
Epi[["FEA"]] <- CreateAssayObject(counts = scFEA)

type <- Epi$MP1_HL
type <- sort(type)

FEA_mtx <- Epi@assays$FEA@counts
FEA_mtx <- FEA_mtx[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("MP1_High", table(type)[1]), rep("MP1_Low", table(type)[2])), levels = c('MP1_High', 'MP1_Low'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(FEA_mtx)
head(design)

# Tunor VS Normal
compare <- makeContrasts(MP1_High - MP1_Low, levels=design)
fit <- lmFit(FEA_mtx, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# 排序
Diff <- Diff %>% arrange(t)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 15, wt = t), top_n(dat_plot, -15, wt = t))
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
  ylab('t value of scFEA score, MP1_On versus Off') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black') +
  ylim(c(-100,100))
ggsave("5.1.1.epithelial_NMF/process_k3/scFEA/MP1_HL_FEA_diff_bar_MP1_top15.pdf", p_MP1, width = 8, height  = 8)

### Monocle2 ---------------------------------------

Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
# MP + SampleType
Epi$SampleType_MP1 <- paste(Epi$SampleType, Epi$MP1_OnOff, sep = "_")

#构造表达及注释数据
exp.matrix <- as(as.matrix(Epi@assays$RNA@data), "sparseMatrix")
feature_ann <- data.frame(gene_id=rownames(exp.matrix),
                          gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(exp.matrix)
exp_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- Epi@meta.data
rownames(sample_ann) <- colnames(exp.matrix)
exp_pd <- new("AnnotatedDataFrame", data = sample_ann)

#生成monocle对象
exp.monocle <- newCellDataSet(
  exp.matrix, 
  phenoData = exp_pd,
  featureData = exp_fd, 
  expressionFamily = negbinomial.size())
head(pData(exp.monocle))
head(fData(exp.monocle))

#计算sizefactor
exp.monocle <- estimateSizeFactors(exp.monocle)
exp.monocle <- estimateDispersions(exp.monocle)

#根据seurat cluster计算差异表达基因并挑选用于构建拟时序轨迹的基因
#MP1差异基因
diff_test_res <- differentialGeneTest(exp.monocle,fullModelFormulaStr = "~MP1_OnOff") 
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
saveRDS(ordering_genes, "./5.1.1.epithelial_NMF/process_k3/Monocle/MP1_ordering_genes.Rds")

#MP1模块的MP1差异基因
Markers <- readRDS("./5.1.1.epithelial_NMF/process_k3/MP_markers.Rds")
ordering_genes <- ordering_genes[Markers[[1]]]

exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)

#EO/LO_MP1差异基因
diff_test_res <- differentialGeneTest(exp.monocle,
                                      fullModelFormulaStr = "~SampleType_MP1",
                                      cores = 20) 
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
saveRDS(ordering_genes, "./5.1.1.epithelial_NMF/process_k3/Monocle/SampleType_MP1_ordering_genes.Rds")

#MP1模块的EO/LO_MP1差异基因(结果与上一种无差异，均纳入了所有272个基因)
#Markers <- readRDS("./5.1.1.epithelial_NMF/process_k3/MP_markers.Rds")
#ordering_genes <- ordering_genes[Markers[[1]]]
#exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)

#plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

saveRDS(exp.monocle, "./5.1.1.epithelial_NMF/process_k3/Monocle/MP1_exp.monocle.Rds")
saveRDS(exp.monocle, "./5.1.1.epithelial_NMF/process_k3/Monocle/SampleType_MP1_exp.monocle.Rds")

#将不同分组情况的拟时序轨迹图画到一起

plot1<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size = 0.2, shuffle = T) + 
  scale_color_gradient(low = "#101D44", high = "#9FB2ED")
plot2<-plot_cell_trajectory(exp.monocle, color_by = "MP1_OnOff",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color[c(2,1)])
plot3<-plot_cell_trajectory(exp.monocle, color_by = "SampleType",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
plot4<-plot_cell_trajectory(exp.monocle, color_by = "State",cell_size = 0.2, shuffle = T) + 
  scale_color_manual(values = color)
pdf("./5.1.1.epithelial_NMF/process_k3/Monocle/SampleType_MP1_trajectory_plot.pdf",width = 12,height = 10)
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = NULL,ncol=2)
dev.off()

#### Staes_SampleType_Roe--------------------------------------

tab <- table(exp.monocle$SampleType, exp.monocle$State)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected

pdf("./5.1.1.epithelial_NMF/process_k3/Monocle/State_Roe_SampleType.pdf", width = 2, height = 4)
bk <- c(seq(0.5,1.5,by=0.01))
pheatmap(t(Roe),
         cluster_cols = F,
         cluster_rows = T,
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

### PAM50 score ---------------------------------------------

Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")

source("../PAM50/PAM50_R/subtypePrediction_functions.R")
paramDir <- "../PAM50/PAM50_R/"
trainCentroids <- paste(paramDir,"pam50_centroids.txt",sep="/")
pamout.centroids <- read.table(trainCentroids,sep="\t",header=T,row.names=1)
inputDir <- "./5.1.1.epithelial_NMF/process_k3/PAM50/input"
inputFile <- "Epi_pam50_exp.txt"
short <- "Epi"

#### prepare pam50 input file --------------------------------------

Epi_exp <- GetAssayData(Epi, assay = "RNA", slot = "data")
rownames(Epi_exp)[rownames(Epi_exp) == 'NDC80'] = 'KNTC2'
rownames(Epi_exp)[rownames(Epi_exp) == 'ORC6'] = 'ORC6L'
rownames(Epi_exp)[rownames(Epi_exp) == 'NUF2'] = 'CDCA1'
Epi_50 <- Epi_exp[rownames(Epi_exp) %in% rownames(pamout.centroids),]
Epi_50 <- Epi_50[order(rownames(Epi_50)),]
write.table(Epi_50, paste(inputDir, inputFile, sep = "/"), col.names = NA, 
            row.names = T, sep = "\t", quote = F)

#### calculate PAM50 subtype score --------------------------------

calibrationParameters<- NA
hasClinical<-FALSE
collapseMethod<-"mean"

trainCentroids <- paste(paramDir,"pam50_centroids.txt",sep="/")
trainFile <- paste(paramDir,"220arrays_nonUBCcommon+12normal_50g.txt",sep="/")
proliferationGenes <-c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
stdArray <-T # just for visualization, and only set to F if many missing genes
predFiles <- paste(inputDir,inputFile,sep="/")

# for subtype only model
glthreshold<- -0.15
ghthreshold<-  0.1

# for subtype + proliferation model
gplthreshold<- -0.25
gphthreshold<-  0.1

# only need train data for visualizations
x<-readarray(trainFile,hr=2)
x$xd<-standardize(medianCtr(x$xd))

# load the published centroids for classifcation
pamout.centroids <- read.table(trainCentroids,sep="\t",header=T,row.names=1)
outputDir <- "./5.1.1.epithelial_NMF/process_k3/PAM50/output"

# read in the data file
if(hasClinical){
  xhr=2
}else{
  xhr=1
}
y<-readarray(predFiles,hr=xhr,method=collapseMethod,impute=F)

# normalization
if(is.na(calibrationParameters)){
  y$xd<-medianCtr(y$xd)
}else{
  if(calibrationParameters != -1){
    medians<-readarray(calibrationFile,hr=1)
    print(paste("calibration to:",dimnames(medians$xd)[[2]][calibrationParameters]))
    tm<-overlapSets(medians$xd,y$xd)
    y$xd<-(tm$y-tm$x[,calibrationParameters])
    #y$xd<-(tm$y-tm$x[,calibrationParameters])/tm$x[,15]
  }
}

num.missing<- NA

if(stdArray){
  y$xd<-standardize(y$xd)
}

erScore<-as.vector(t(y$xd["ESR1",]))
her2Score<-as.vector(t(y$xd["ERBB2",]))

# assign the subtype scores and calculate the proliferation score
this.proliferationGenes<-dimnames(y$xd)[[1]] %in% proliferationGenes

prolifScore<-apply(y$xd[this.proliferationGenes,],2,mean,na.rm=T)

##### euclidean distance -------------------------------------------------
out <- sspPredict(pamout.centroids, classes="", y$xd, std=F, 
                  distm="euclidean", centroids=T)
out$distances<- -1*out$distances

call.conf<-c()
for(j in 1:length(out$predictions)){
  call.conf[j]<- 1-cor.test(out$testData[,j],out$centroids[,which(colnames(pamout.centroids)==out$predictions[j])],method="euclidean")$p.value
}
call.conf<-round(call.conf,2)

# calculate the risk scores
genomic <- 0.04210193*out$distances[,1] + 0.12466938*out$distances[,2] + -0.35235561*out$distances[,3] + 0.14213283*out$distances[,4]
genomicWprolif <- -0.0009299747*out$distances[,1] + 0.0692289192*out$distances[,2] + -0.0951505484*out$distances[,3] +  0.0493487685*out$distances[,4] + 0.3385116381*prolifScore
if(hasClinical){
  xT<-as.numeric(as.vector(y$classes$T))
  combined <- 0.0442770*out$distances[,1] + 0.1170297*out$distances[,2] + -0.2608388*out$distances[,3] + 0.1055908*out$distances[,4] + 0.1813751*xT
  combinedWprolif <- -0.009383416*out$distances[,1] +  0.073725503*out$distances[,2] + -0.090436516*out$distances[,3] + 0.053013865*out$distances[,4] + 0.131605960*xT + 0.327259375*prolifScore
}

# threshold the risk score
griskgroups<-genomic
griskgroups[genomic>ghthreshold]<-"high"
griskgroups[genomic>glthreshold & genomic<ghthreshold]<-"med"
griskgroups[genomic<glthreshold]<-"low"
gpriskgroups<-genomicWprolif
gpriskgroups[genomicWprolif>gphthreshold]<-"high"
gpriskgroups[genomicWprolif>gplthreshold & genomicWprolif<gphthreshold]<-"med"
gpriskgroups[genomicWprolif<gplthreshold]<-"low"

genomic<- 100* (genomic + 0.35 ) / 0.85
genomicWprolif<- 100* (genomicWprolif + 0.35 ) / 0.85

# write output files

outtable<-cbind(out$distances, out$predictions, genomic, griskgroups, prolifScore, genomicWprolif, gpriskgroups, erScore, her2Score)
dimnames(outtable)[[2]]<-c("Basal","Her2","LumA","LumB","Normal","Call",
                           "ROR-S (Subtype Only)","ROR-S Group (Subtype Only)","Proliferation Score", 
                           "ROR-P (Subtype + Proliferation)","ROR-P Group (Subtype + Proliferation)",
                           "ER","Her2")

outFile<- paste(outputDir,paste(short,"_euclidean_pam50scores.txt",sep=""),sep="/")
write.table(outtable,outFile,sep="\t",col.names=NA)

##### Boxplot - MP1 OnOff Basal LuminalA LuminalB --------------------------

score <- out$distances[,c(1,3,4)] %>% as.data.frame()
rownames(score) <- rownames(outtable)
colnames(score) <- colnames(outtable)[c(1,3,4)]
MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

mtx <- cbind(score, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
library(reshape2)
data <- melt(mtx)
colnames(data)[2] <- "SubType"

library(ggpubr)
P = ggplot(data, aes(x = SubType, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("PAM50 score") +
  ylab("-Euclidean distance") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf(paste(outputDir,"Epi_MP1_OnOff_PAM50_score_euclidean.pdf",sep = "/"), width = 4, height = 3)
P
dev.off()

##### Boxplot - SampleType - MP1 OnOff Basal LuminalA LuminalB --------------------------

score <- read.table("5.1.1.epithelial_NMF/process_k3/PAM50/output/Epi_euclidean_pam50scores.txt", header = T, row.names = 1)
score <- score[,c(1,3,4)]
Epi$SampleTypeMP1 <- paste(Epi$SampleType,Epi$MP1_OnOff,sep="_")
MP1_OnOff <- Epi$SampleTypeMP1 %>% as.data.frame()
colnames(MP1_OnOff) <- "Group_Onoff"

mtx <- cbind(score, MP1_OnOff)
mtx$Group_Onoff <- factor(mtx$Group_Onoff, levels = c("EOPC_On", "LOPC_On", "EOPC_Off", "LOPC_Off"))
mtx <- mtx[order(mtx$Group_Onoff),]
library(reshape2)
data <- melt(mtx)
colnames(data)[2] <- "SubType"

library(ggpubr)
P = ggplot(data, aes(x = SubType, y = value, fill = Group_Onoff)) +
  scale_fill_manual(values = c("#E46B4FCC", "#8B8B00CC", "#B8860BCC", "#483D8BCC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("PAM50 score") +
  ylab("-Euclidean distance") +
  xlab("")
pdf("5.1.1.epithelial_NMF/process_k3/PAM50/output/Epi_SampleType_MP1_OnOff_PAM50_score_euclidean.pdf", width = 4, height = 3)
P
dev.off()

##### Boxplot - MP1 HL Basal LuminalA LuminalB -----------------------------------

score <- out$distances[,c(1,3,4)] %>% as.data.frame()
rownames(score) <- rownames(outtable)
colnames(score) <- colnames(outtable)[c(1,3,4)]
MP1_HL <- Epi$MP1_HL %>% as.data.frame()
colnames(MP1_HL) <- "MP1_HL"

mtx <- cbind(score, MP1_HL)
mtx$MP1_HL <- factor(mtx$MP1_HL, levels = c("High", "Low"))
mtx <- mtx[order(mtx$MP1_HL),]
library(reshape2)
data <- melt(mtx)
colnames(data)[2] <- "SubType"

library(ggpubr)
P = ggplot(data, aes(x = SubType, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("PAM50 score") +
  ylab("-Euclidean distance") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_HL),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf(paste(outputDir,"Epi_MP1_HL_PAM50_score_euclidean.pdf",sep = "/"), width = 4, height = 3)
P
dev.off()

##### Boxplot - MP1 OnOff HER2 score --------------------------

score <- data.frame(HER2Score = her2Score)
rownames(score) <- rownames(outtable)
MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

mtx <- cbind(score, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(1)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("HER2 score") +
  ylab("-Euclidean distance") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf(paste(outputDir,"Epi_MP1_OnOff_HER2_score_euclidean.pdf",sep = "/"), width = 2.5, height = 3)
P
dev.off()

##### Boxplot - MP1 HL HER2 score -----------------------------------

score <- data.frame(HER2Score = her2Score)
rownames(score) <- rownames(outtable)
MP1_HL <- Epi$MP1_HL %>% as.data.frame()
colnames(MP1_HL) <- "MP1_HL"

mtx <- cbind(score, MP1_HL)
mtx$MP1_HL <- factor(mtx$MP1_HL, levels = c("High", "Low"))
mtx <- mtx[order(mtx$MP1_HL),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("HER2 score") +
  ylab("-Euclidean distance") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_HL),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf(paste(outputDir,"Epi_MP1_HL_HER2_score_euclidean.pdf",sep = "/"), width = 4, height = 3)
P
dev.off()

##### spearman correlation -------------------------------------------------
out <- sspPredict(pamout.centroids, classes="", y$xd, std=F, 
                  distm="spearman", centroids=T)
out$distances<- -1*out$distances

call.conf<-c()
for(j in 1:length(out$predictions)){
  call.conf[j]<- 1-cor.test(out$testData[,j],out$centroids[,which(colnames(pamout.centroids)==out$predictions[j])],method="spearman")$p.value
}
call.conf<-round(call.conf,2)

# calculate the risk scores
genomic <- 0.04210193*out$distances[,1] + 0.12466938*out$distances[,2] + -0.35235561*out$distances[,3] + 0.14213283*out$distances[,4]
genomicWprolif <- -0.0009299747*out$distances[,1] + 0.0692289192*out$distances[,2] + -0.0951505484*out$distances[,3] +  0.0493487685*out$distances[,4] + 0.3385116381*prolifScore
if(hasClinical){
  xT<-as.numeric(as.vector(y$classes$T))
  combined <- 0.0442770*out$distances[,1] + 0.1170297*out$distances[,2] + -0.2608388*out$distances[,3] + 0.1055908*out$distances[,4] + 0.1813751*xT
  combinedWprolif <- -0.009383416*out$distances[,1] +  0.073725503*out$distances[,2] + -0.090436516*out$distances[,3] + 0.053013865*out$distances[,4] + 0.131605960*xT + 0.327259375*prolifScore
}

# threshold the risk score
griskgroups<-genomic
griskgroups[genomic>ghthreshold]<-"high"
griskgroups[genomic>glthreshold & genomic<ghthreshold]<-"med"
griskgroups[genomic<glthreshold]<-"low"
gpriskgroups<-genomicWprolif
gpriskgroups[genomicWprolif>gphthreshold]<-"high"
gpriskgroups[genomicWprolif>gplthreshold & genomicWprolif<gphthreshold]<-"med"
gpriskgroups[genomicWprolif<gplthreshold]<-"low"

genomic<- 100* (genomic + 0.35 ) / 0.85
genomicWprolif<- 100* (genomicWprolif + 0.35 ) / 0.85

# write output files

outtable<-cbind(out$distances, out$predictions, call.conf, genomic, griskgroups, prolifScore, genomicWprolif, gpriskgroups, erScore, her2Score)
dimnames(outtable)[[2]]<-c("Basal","Her2","LumA","LumB","Normal","Call","Confidence",
                           "ROR-S (Subtype Only)","ROR-S Group (Subtype Only)","Proliferation Score", 
                           "ROR-P (Subtype + Proliferation)","ROR-P Group (Subtype + Proliferation)",
                           "ER","Her2")

outFile<- paste(outputDir,paste(short,"_spearman_pam50scores.txt",sep=""),sep="/")
write.table(outtable,outFile,sep="\t",col.names=NA)

##### Boxplot - MP1 OnOff Basal LuminalA LuminalB --------------------------

score <- out$distances[,c(1,3,4)] %>% as.data.frame()
rownames(score) <- rownames(outtable)
colnames(score) <- colnames(outtable)[c(1,3,4)]
MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

mtx <- cbind(score, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
library(reshape2)
data <- melt(mtx)
colnames(data)[2] <- "SubType"

library(ggpubr)
P = ggplot(data, aes(x = SubType, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("PAM50 score") +
  ylab("Pearson Cor. Coef.") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf(paste(outputDir,"Epi_MP1_OnOff_PAM50_score_pearson.pdf",sep = "/"), width = 4, height = 3)
P
dev.off()

##### Boxplot - MP1 HL Basal LuminalA LuminalB -----------------------------------

score <- out$distances[,c(1,3,4)] %>% as.data.frame()
rownames(score) <- rownames(outtable)
colnames(score) <- colnames(outtable)[c(1,3,4)]
MP1_HL <- Epi$MP1_HL %>% as.data.frame()
colnames(MP1_HL) <- "MP1_HL"

mtx <- cbind(score, MP1_HL)
mtx$MP1_HL <- factor(mtx$MP1_HL, levels = c("High", "Low"))
mtx <- mtx[order(mtx$MP1_HL),]
library(reshape2)
data <- melt(mtx)
colnames(data)[2] <- "SubType"

library(ggpubr)
P = ggplot(data, aes(x = SubType, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("PAM50 score") +
  ylab("Pearson Cor. Coef.") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_HL),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf(paste(outputDir,"Epi_MP1_HL_PAM50_score_pearson.pdf",sep = "/"), width = 4, height = 3)
P
dev.off()

### ERBB/PI3K/MAPK pathway ---------------------------------------------

#### ERBB genes ------------------------------------

Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
exp <- Epi@assays$RNA@data[c("EGFR","ERBB2","ERBB3","ERBB4"),] %>% as.data.frame() %>% t()
colnames(exp)[c(1,2)] <- c("EGFR (ERBB1)","ERBB2 (HER2)")

# MP1_OnOff

MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

mtx <- cbind(exp, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("ERBB Genes") +
  ylab("Norm. Exp.") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/ERBB/BoxPlot_Epi_MP1_OnOff_ERBBgenes.pdf", width = 5, height = 4)
P
dev.off()

# MP1_OnOff

SampleType <- Epi$SampleType %>% as.data.frame()
colnames(SampleType) <- "SampleType"

mtx <- cbind(exp, SampleType)
mtx$SampleType <- factor(mtx$SampleType, levels = c("EOPC", "LOPC"))
mtx <- mtx[order(mtx$SampleType),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = SampleType)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("ERBB Genes") +
  ylab("Norm. Exp.") +
  xlab("") + 
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/ERBB/BoxPlot_Epi_MP1_OnOff_ERBBgenes.pdf", width = 5, height = 4)
P
dev.off()

# MP1 HL

MP1_HL <- Epi$MP1_HL %>% as.data.frame()
colnames(MP1_HL) <- "MP1_HL"

mtx <- cbind(exp, MP1_HL)
mtx$MP1_HL <- factor(mtx$MP1_HL, levels = c("High", "Low"))
mtx <- mtx[order(mtx$MP1_HL),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("ERBB Genes") +
  ylab("Norm. Exp.") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_HL),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/ERBB/BoxPlot_Epi_MP1_HL_ERBBgenes_AUCell.pdf", width = 4, height = 4)
P
dev.off()

#### Pathways --------------------------------------------

##### AUCell --------------------------------------

KEGG <- readRDS("kegg_hsa_gmt.Rds")
#PI3K-Akt signaling pathway
#MAPK signaling pathway
#JAK-STAT signaling pathway

gs_selected <- KEGG[KEGG$term %in% c("PI3K-Akt signaling pathway",
                                     "MAPK signaling pathway",
                                     "JAK-STAT signaling pathway"),]

library(AUCell)
gs <- list()
for (i in c("PI3K-Akt signaling pathway",
            "MAPK signaling pathway",
            "JAK-STAT signaling pathway")) {
  gs[[i]] <- gs_selected[which(gs_selected$term == i),]$gene
}

expMtx <- GetAssayData(Epi, assay = "RNA", slot = "data")
cells_rankings <- AUCell_buildRankings(expMtx, nCores=10, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.1.1.epithelial_NMF/process_k3/ERBB/3gs_AUC.Rds")

AUCell <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

#MP1 OnOff
MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

mtx <- cbind(AUCell, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("KEGG pathways") +
  ylab("AUCell score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/ERBB/BoxPlot_Epi_MP1_OnOff_3gs_AUCell.pdf", width = 4, height = 4)
P
dev.off()

#MP1 HL
MP1_HL <- Epi$MP1_HL %>% as.data.frame()
colnames(MP1_HL) <- "MP1_HL"

mtx <- cbind(AUCell, MP1_HL)
mtx$MP1_HL <- factor(mtx$MP1_HL, levels = c("High", "Low"))
mtx <- mtx[order(mtx$MP1_HL),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("KEGG pathways") +
  ylab("AUCell score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_HL),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/ERBB/BoxPlot_Epi_MP1_HL_3gs_AUCell.pdf", width = 4, height = 4)
P
dev.off()

##### GSVA --------------------------------------

Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
KEGG <- readRDS("kegg_hsa_gmt.Rds")
#PI3K-Akt signaling pathway
#MAPK signaling pathway
#JAK-STAT signaling pathway

gs_selected <- KEGG[KEGG$term %in% c("ErbB signaling pathway",
                                     "PI3K-Akt signaling pathway",
                                     "MAPK signaling pathway",
                                     "JAK-STAT signaling pathway"),]
gs <- list()
for (i in c("ErbB signaling pathway",
            "PI3K-Akt signaling pathway",
            "MAPK signaling pathway",
            "JAK-STAT signaling pathway")) {
  gs[[i]] <- gs_selected[which(gs_selected$term == i),]$gene
}

library(GSVA)

expMtx <- GetAssayData(Epi, assay = "RNA", slot = "data")
GSVA <- gsva(expMtx, gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "./5.1.1.epithelial_NMF/process_k3/ERBB/3gs_GSVA.Rds")

#MP1 OnOff
MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

GSVA <- GSVA %>% as.data.frame() %>% t()
mtx <- cbind(GSVA, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("KEGG pathways") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/ERBB/BoxPlot_Epi_MP1_OnOff_3gs_GSVA.pdf", width = 4, height = 4)
P
dev.off()

#MP1 HL
MP1_HL <- Epi$MP1_HL %>% as.data.frame()
colnames(MP1_HL) <- "MP1_HL"

mtx <- cbind(GSVA, MP1_HL)
mtx$MP1_HL <- factor(mtx$MP1_HL, levels = c("High", "Low"))
mtx <- mtx[order(mtx$MP1_HL),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("KEGG pathways") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_HL),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/ERBB/BoxPlot_Epi_MP1_HL_3gs_GSVA.pdf", width = 4, height = 4)
P
dev.off()

### Lipid met GSVA --------------------------------------

Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
KEGG <- readRDS("kegg_hsa_gmt.Rds")
#"Fatty acid biosynthesis",
#"Fatty acid degradation",
#"Biosynthesis of unsaturated fatty acids",
#"Fatty acid metabolism",
#"Cholesterol metabolism"
#"Glycerolipid metabolism"        
#"Glycerophospholipid metabolism"
#"Glycosphingolipid biosynthesis"
#"Sphingolipid metabolism"

gs_selected <- KEGG[KEGG$term %in% c("Fatty acid biosynthesis",
                                     "Fatty acid degradation",
                                     "Biosynthesis of unsaturated fatty acids",
                                     "Fatty acid metabolism",
                                     "Cholesterol metabolism",
                                     "Glycerolipid metabolism",      
                                     "Glycerophospholipid metabolism",
                                     "Glycosphingolipid biosynthesis",
                                     "Sphingolipid metabolism"),]
gs <- list()
for (i in c("Fatty acid biosynthesis",
            "Fatty acid degradation",
            "Biosynthesis of unsaturated fatty acids",
            "Fatty acid metabolism",
            "Cholesterol metabolism",
            "Glycerolipid metabolism",      
            "Glycerophospholipid metabolism",
            "Glycosphingolipid biosynthesis",
            "Sphingolipid metabolism")) {
  gs[[i]] <- gs_selected[which(gs_selected$term == i),]$gene
}

library(GSVA)

expMtx <- GetAssayData(Epi, assay = "RNA", slot = "data")
GSVA <- gsva(expMtx, gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "./5.1.1.epithelial_NMF/process_k3/lipid/9gs_GSVA.Rds")

#MP1 OnOff
MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

GSVA <- GSVA %>% as.data.frame() %>% t()
mtx <- cbind(GSVA, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=8, colour="black"), 
        axis.title.y=element_text(size = 10),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("KEGG pathways") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("./5.1.1.epithelial_NMF/process_k3/lipid/BoxPlot_Epi_MP1_OnOff_3gs_GSVA.pdf", width = 5, height = 4)
P
dev.off()

#MP1 HL
MP1_HL <- Epi$MP1_HL %>% as.data.frame()
colnames(MP1_HL) <- "MP1_HL"

mtx <- cbind(GSVA, MP1_HL)
mtx$MP1_HL <- factor(mtx$MP1_HL, levels = c("High", "Low"))
mtx <- mtx[order(mtx$MP1_HL),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=8, colour="black"), 
        axis.title.y=element_text(size = 10),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("KEGG pathways") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_HL),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("./5.1.1.epithelial_NMF/process_k3/lipid/BoxPlot_Epi_MP1_HL_3gs_GSVA.pdf", width = 5, height = 4)
P
dev.off()

### ARPC_MSPC_NEPC geneset (Han et al. 2022) --------------------------------------------

PC_gs <- read.table("ARPC_MSPC_NEPC_gs.txt", sep="\t", header=T)

library(AUCell)
gs <- list()
for (i in c("ARPC",
            "MSPC",
            "NEPC")) {
  gs[[i]] <- PC_gs[which(PC_gs$term == i),]$gene
}

expMtx <- GetAssayData(Epi, assay = "RNA", slot = "data")
cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.1.1.epithelial_NMF/process_k3/mPCa_gs/mPC_gs_AUC.Rds")

AUCell <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

#MP1 OnOff
MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

mtx <- cbind(AUCell, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Han et al. 2022") +
  ylab("AUCell score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/mPCa_gs/BoxPlot_Epi_MP1_OnOff_mPC_gs_AUCell.pdf", width = 4, height = 4)
P
dev.off()

#MP1 HL
MP1_HL <- Epi$MP1_HL %>% as.data.frame()
colnames(MP1_HL) <- "MP1_HL"

mtx <- cbind(AUCell, MP1_HL)
mtx$MP1_HL <- factor(mtx$MP1_HL, levels = c("High", "Low"))
mtx <- mtx[order(mtx$MP1_HL),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_HL)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Han et al. 2022") +
  ylab("AUCell score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_HL),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/mPCa_gs/BoxPlot_Epi_MP1_HL_mPC_gs_AUCell", width = 4, height = 4)
P
dev.off()

##correlation between MP1 and mPC_score

Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")
MP1 <- Epi@assays$MP_AUCell@data["MP1",] %>% as.data.frame
colnames(MP1) <- "MP1"
mtx <- cbind(AUCell, MP1)

Plist <- list()
for (i in c(1:3)) {

  P <- ggplot(mtx, aes(x=mtx[,i], y=mtx[,4])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 1.5) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = mtx[,i], y = mtx[,4])) +
    theme(axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.y=element_text(size = 14),
          axis.title.x=element_text(size = 14),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(mtx)[i]," AUCell")) + ylab(paste0(colnames(mtx)[4]," AUCell"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/mPCa_gs/CorPlot_MP1_",colnames(mtx)[i],".pdf"), P, width = 2, height = 2)
}

#EO/LO
SampleType <- Epi$SampleType %>% as.data.frame()
colnames(SampleType) <- "SampleType"

mtx <- cbind(AUCell, SampleType)
mtx$SampleType <- factor(mtx$SampleType, levels = c("EOPC", "LOPC"))
mtx <- mtx[order(mtx$SampleType),]
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = SampleType)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Han et al. 2022") +
  ylab("AUCell score") +
  xlab("") + 
  stat_compare_means(aes(group=SampleType),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/mPCa_gs/BoxPlot_Epi_EO_LO_mPC_gs_AUCell", width = 4, height = 4)
P
dev.off()

### EMT pathway --------------------------------------------

Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")

cancersea <- read.gmt("hallmark_cancersea.gmt")

select_gs <- list(EMT = cancersea$gene[cancersea$term == "EMT"],
                  EPITHELIAL_MESENCHYMAL_TRANSITION = cancersea$gene[cancersea$term == "EPITHELIAL_MESENCHYMAL_TRANSITION"])

library(GSVA)

expMtx <- GetAssayData(Epi, assay = "RNA", slot = "data")
GSVA <- gsva(expMtx, select_gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "./5.1.1.epithelial_NMF/process_k3/EMT/GSVA.Rds")

#MP1 OnOff
MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

GSVA <- GSVA %>% as.data.frame() %>% t()
mtx <- cbind(GSVA, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
mtx <- mtx[,-1]
colnames(mtx)[1] <- "HALLMARK_EMT"
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("EMT") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/EMT/BoxPlot_Epi_MP1_OnOff_GSVA.pdf", width = 3, height = 4)
P
dev.off()

### hypoxia pathway --------------------------------------------

Epi <- readRDS("5.1.1.epithelial_NMF/process_k3/Epi_MPAUCell_OnOff.Rds")

cancersea <- read.gmt("hallmark_cancersea.gmt")
BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
BIOCARTA$gs_name <- gsub("BIOCARTA_", "", BIOCARTA$gs_name)

select_gs <- list(Hypoxia = cancersea$gene[cancersea$term == "HYPOXIA"],
                  HIF_Pathway = BIOCARTA$gene[BIOCARTA$gs_name == "HIF_PATHWAY"])

library(GSVA)

expMtx <- GetAssayData(Epi, assay = "RNA", slot = "data")
GSVA <- gsva(expMtx, select_gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "./5.1.1.epithelial_NMF/process_k3/hypoxia/GSVA.Rds")

#MP1 OnOff
MP1_OnOff <- Epi$MP1_OnOff %>% as.data.frame()
colnames(MP1_OnOff) <- "MP1_OnOff"

GSVA <- GSVA %>% as.data.frame() %>% t()
mtx <- cbind(GSVA, MP1_OnOff)
mtx$MP1_OnOff <- factor(mtx$MP1_OnOff, levels = c("On", "Off"))
mtx <- mtx[order(mtx$MP1_OnOff),]
mtx <- mtx[,-2]
colnames(mtx)[1] <- "HALLMARK_Hypoxia"
library(reshape2)
data <- melt(mtx)

library(ggpubr)
P = ggplot(data, aes(x = variable, y = value, fill = MP1_OnOff)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Hypoxia") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=MP1_OnOff),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/hypoxia/BoxPlot_Epi_MP1_OnOff_GSVA.pdf", width = 3, height = 4)
P
dev.off()

### Bulk hypoxia and EMT (MP1 score / EO vs LO) ----------------------------------

cancersea <- read.gmt("hallmark_cancersea.gmt")

select_gs <- list(EMT = cancersea$gene[cancersea$term == "EMT"],
                  Hypoxia = cancersea$gene[cancersea$term == "HYPOXIA"])

#read gene matrix and MP score
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

#GSVA score
GSVA <- gsva(expMtx, select_gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "5.1.1.epithelial_NMF/process_k3/Bulk validation/Hypo_EMT_GSVA.Rds")

#### MP1 and GSVA --------------------------------------------------------

MP1_GSVA <- merge(t(AUCell_mtx), t(GSVA), by = 0) %>% column_to_rownames("Row.names")
MP1_GSVA <- MP1_GSVA[-grep(".11$|.06$", rownames(MP1_GSVA)),]

Plist <- list()
for (i in ((ncol(MP1_GSVA)-1) : (ncol(MP1_GSVA)))) {
  
  P <- ggplot(MP1_GSVA, aes(x=MP1_GSVA[,1], y=MP1_GSVA[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = MP1_GSVA[,1], y = MP1_GSVA[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(MP1_GSVA)[1]," AUCell")) + ylab(paste0(colnames(MP1_GSVA)[i]," GSVA"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/Hypo_EMT/CorPlot_MP1_",colnames(MP1_GSVA)[i],".pdf"), P, width = 2, height = 2)
}

#### EO vs LO and GSVA --------------------------------------------------------

Age_GSVA <- merge(age, t(GSVA), by = 0) %>% column_to_rownames("Row.names")

# correlation

for (i in ((ncol(Age_GSVA)-1) : (ncol(Age_GSVA)))) {
  
  P <- ggplot(Age_GSVA, aes(x=Age_GSVA[,1], y=Age_GSVA[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.76) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = Age_GSVA[,1], y = Age_GSVA[,i]), size = 2.5 ) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_GSVA)[i]," GSVA"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/Hypo_EMT/CorPlot_Age_",colnames(Age_GSVA)[i],".pdf"), P, width = 2, height = 2)
}

# boxplot

data <- melt(Age_GSVA[,c(2,3)])
P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/Bulk validation/Hypo_EMT/BoxPlot_EMT_EOLO.pdf", width = 1.5, height = 4)
P + NoLegend()
dev.off()

data <- melt(Age_GSVA[,c(2,4)])
P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/Bulk validation/Hypo_EMT/BoxPlot_Hypoxia_EOLO.pdf", width = 3, height = 4)
P
dev.off()

### Bulk lipid metabolism (MP1 score / EO vs LO) ----------------------------------

KEGG <- readRDS("kegg_hsa_gmt.Rds")

gs_selected <- KEGG[KEGG$term %in% c("Fatty acid biosynthesis",
                                     "Fatty acid degradation",
                                     "Biosynthesis of unsaturated fatty acids",
                                     "Fatty acid metabolism",
                                     "Cholesterol metabolism",
                                     "Glycerolipid metabolism",      
                                     "Glycerophospholipid metabolism",
                                     "Glycosphingolipid biosynthesis",
                                     "Sphingolipid metabolism"),]
gs <- list()
for (i in c("Fatty acid biosynthesis",
            "Fatty acid degradation",
            "Biosynthesis of unsaturated fatty acids",
            "Fatty acid metabolism",
            "Cholesterol metabolism",
            "Glycerolipid metabolism",      
            "Glycerophospholipid metabolism",
            "Glycosphingolipid biosynthesis",
            "Sphingolipid metabolism")) {
  gs[[i]] <- gs_selected[which(gs_selected$term == i),]$gene
}

names(gs) <- c("FA biosynthesis","FA degradation","Biosynthesis of UFA",
               "FA metabolism","Cho metabolism","GL metabolism",
               "GP metabolism","GSP metabolism","SP metabolism")

#read gene matrix and MP score
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

#GSVA score
GSVA <- gsva(expMtx, gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "5.1.1.epithelial_NMF/process_k3/Bulk validation/lipid_GSVA.Rds")

#### MP1 and GSVA --------------------------------------------------------

MP1_GSVA <- merge(t(AUCell_mtx), t(GSVA), by = 0) %>% column_to_rownames("Row.names")
MP1_GSVA <- MP1_GSVA[-grep(".11$|.06$", rownames(MP1_GSVA)),]

Plist <- list()
for (i in (7 : (ncol(MP1_GSVA)))) {
  
  P <- ggplot(MP1_GSVA, aes(x=MP1_GSVA[,1], y=MP1_GSVA[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = MP1_GSVA[,1], y = MP1_GSVA[,i]), size = 2) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(MP1_GSVA)[1]," AUCell")) + ylab(paste0(colnames(MP1_GSVA)[i]," GSVA"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/lipid/CorPlot_MP1_",colnames(MP1_GSVA)[i],".pdf"), P, width = 1.5, height = 1.5)
}

#### EO vs LO and GSVA --------------------------------------------------------

Age_GSVA <- merge(age, t(GSVA), by = 0) %>% column_to_rownames("Row.names")

# correlation

for (i in (3 : (ncol(Age_GSVA)))) {
  
  P <- ggplot(Age_GSVA, aes(x=Age_GSVA[,1], y=Age_GSVA[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = Age_GSVA[,1], y = Age_GSVA[,i]), size = 2) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_GSVA)[i]," GSVA"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/lipid/CorPlot_Age_",colnames(Age_GSVA)[i],".pdf"), P, width = 1.5, height = 1.5)
}

# boxplot

data <- melt(Age_GSVA[,c(2:ncol(Age_GSVA))])
P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/Bulk validation/lipid/BoxPlot_lipid_EOLO.pdf", width = 4, height = 4)
P + NoLegend()
dev.off()

### Bulk ErbB (MP1 score / EO vs LO) ----------------------------------

#### ErbB genes ---------------------------------------

exp <- expMtx[c("EGFR","ERBB2","ERBB3","ERBB4"),] %>% as.data.frame()
exp <- exp+1
exp <- log(exp)
Age_exp <- merge(age, t(exp), by = 0) %>% column_to_rownames("Row.names")

##### MP1 and ErbB genes -----------------------------------

MP1_exp <- merge(t(AUCell_mtx), t(exp), by = 0) %>% column_to_rownames("Row.names")
MP1_exp <- MP1_exp[-grep(".11$|.06$", rownames(MP1_exp)),]

Plist <- list()
for (i in (7 : (ncol(MP1_exp)))) {
  
  P <- ggplot(MP1_exp, aes(x=MP1_exp[,1], y=MP1_exp[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = MP1_exp[,1], y = MP1_exp[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(MP1_exp)[1]," AUCell")) + ylab(paste0(colnames(MP1_exp)[i]," GSVA"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/ErbB/CorPlot_MP1_",colnames(MP1_exp)[i],".pdf"), P, width = 2, height = 2)
}

##### age and ErbB genes -----------------------------------

# correlation

for (i in (3 : (ncol(Age_exp)))) {
  
  P <- ggplot(Age_exp, aes(x=Age_exp[,1], y=Age_exp[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = Age_exp[,1], y = Age_exp[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_exp)[i]," Exp"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/ErbB/CorPlot_Age_",colnames(Age_exp)[i],".pdf"), P, width = 2, height = 2)
}

# boxplot

data <- melt(Age_exp[,c(2:ncol(Age_exp))])
P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("log(TPM + 1)") +
  xlab("") + 
  stat_compare_means(aes(group=Type),
                     method="t.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/Bulk validation/ErbB/BoxPlot_ErbB_genes_EOLO.pdf", width = 2, height = 4)
P + NoLegend()
dev.off()

#### ErbB pathways ---------------------------------------

KEGG <- readRDS("kegg_hsa_gmt.Rds")
gs_selected <- KEGG[KEGG$term %in% c("PI3K-Akt signaling pathway",
                                     "MAPK signaling pathway",
                                     "JAK-STAT signaling pathway"),]
gs <- list()
for (i in c("PI3K-Akt signaling pathway",
            "MAPK signaling pathway",
            "JAK-STAT signaling pathway")) {
  gs[[i]] <- gs_selected[which(gs_selected$term == i),]$gene
}

#read gene matrix and MP score
expMtx <- read.table("./0.TCGA/mRNA_Matrix_TPM.txt", row.names = 1, header = T, sep = "\t")
expMtx <- as.matrix(expMtx)

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/ErbB/MPScore_AUC_TCGA.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)

clinical <- read.table("./0.TCGA/Clinical.txt", row.names = 1, header = T, sep = "\t")
age <- as.data.frame(clinical[,1])
colnames(age) <- "age"
rownames(age) <- gsub("_", ".", rownames(clinical))
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"

#GSVA score
GSVA <- gsva(expMtx, gs, method = "gsva",
             kcdf="Gaussian",
             verbose=T, 
             parallel.sz = 20)
saveRDS(GSVA, "5.1.1.epithelial_NMF/process_k3/Bulk validation/ErbB_GSVA.Rds")

##### MP1 and GSVA --------------------------------------------------------

MP1_GSVA <- merge(t(AUCell_mtx), t(GSVA), by = 0) %>% column_to_rownames("Row.names")
MP1_GSVA <- MP1_GSVA[-grep(".11$|.06$", rownames(MP1_GSVA)),]

Plist <- list()
for (i in (7 : (ncol(MP1_GSVA)))) {
  
  P <- ggplot(MP1_GSVA, aes(x=MP1_GSVA[,1], y=MP1_GSVA[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = MP1_GSVA[,1], y = MP1_GSVA[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(MP1_GSVA)[1]," AUCell")) + ylab(paste0(colnames(MP1_GSVA)[i]," GSVA"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/ErbB/CorPlot_MP1_",colnames(MP1_GSVA)[i],".pdf"), P, width = 2, height = 2)
}

#### EO vs LO and GSVA --------------------------------------------------------

Age_GSVA <- merge(age, t(GSVA), by = 0) %>% column_to_rownames("Row.names")

# correlation

Plist <- list()
for (i in (3 : (ncol(Age_GSVA)))) {
  
  P <- ggplot(Age_GSVA, aes(x=Age_GSVA[,1], y=Age_GSVA[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = Age_GSVA[,1], y = Age_GSVA[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_GSVA)[i]," GSVA"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/ErbB/CorPlot_Age_",colnames(Age_GSVA)[i],".pdf"), P, width = 2, height = 2)
}

# boxplot

data <- melt(Age_GSVA[,c(2:ncol(Age_GSVA))])
P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("GSVA score") +
  xlab("") + 
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/Bulk validation/ErbB/BoxPlot_lipid_EOLO.pdf", width = 2, height = 6)
P + NoLegend()
dev.off()

### Bulk PAM50 (MP1 score / EO vs LO) ----------------------------------

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

#### calculate PAM50 score ------------------------------
source("../PAM50/PAM50_R/subtypePrediction_functions.R")
paramDir <- "../PAM50/PAM50_R/"
trainCentroids <- paste(paramDir,"pam50_centroids.txt",sep="/")
pamout.centroids <- read.table(trainCentroids,sep="\t",header=T,row.names=1)
inputDir <- "5.1.1.epithelial_NMF/process_k3/Bulk validation/PAM50/input"
inputFile <- "exp.txt"
short <- "TCGA"
outputDir <- "5.1.1.epithelial_NMF/process_k3/Bulk validation/PAM50/output"

#### prepare pam50 input file --------------------------------------

rownames(expMtx)[rownames(expMtx) == 'NDC80'] = 'KNTC2'
rownames(expMtx)[rownames(expMtx) == 'ORC6'] = 'ORC6L'
rownames(expMtx)[rownames(expMtx) == 'NUF2'] = 'CDCA1'
exp_50 <- expMtx[rownames(expMtx) %in% rownames(pamout.centroids),]
exp_50 <- exp_50[order(rownames(exp_50)),]
exp_50 <- log(exp_50 + 1)
write.table(exp_50, "5.1.1.epithelial_NMF/process_k3/Bulk validation/PAM50/input/exp.txt", 
            sep = "\t", row.names = T, col.names = NA, quote = F)
#### calculate PAM50 subtype score --------------------------------

calibrationParameters<- NA
hasClinical<-FALSE
collapseMethod<-"mean"

trainCentroids <- paste(paramDir,"pam50_centroids.txt",sep="/")
trainFile <- paste(paramDir,"220arrays_nonUBCcommon+12normal_50g.txt",sep="/")
proliferationGenes <-c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
stdArray <-T # just for visualization, and only set to F if many missing genes
predFiles <- paste(inputDir,inputFile,sep="/")

# for subtype only model
glthreshold<- -0.15
ghthreshold<-  0.1

# for subtype + proliferation model
gplthreshold<- -0.25
gphthreshold<-  0.1

# only need train data for visualizations
x<-readarray(trainFile,hr=2)
x$xd<-standardize(medianCtr(x$xd))

# load the published centroids for classifcation
pamout.centroids <- read.table(trainCentroids,sep="\t",header=T,row.names=1)

# read in the data file
if(hasClinical){
  xhr=2
}else{
  xhr=1
}
y<-readarray(predFiles,hr=xhr,method=collapseMethod,impute=F)

# normalization
if(is.na(calibrationParameters)){
  y$xd<-medianCtr(y$xd)
}else{
  if(calibrationParameters != -1){
    medians<-readarray(calibrationFile,hr=1)
    print(paste("calibration to:",dimnames(medians$xd)[[2]][calibrationParameters]))
    tm<-overlapSets(medians$xd,y$xd)
    y$xd<-(tm$y-tm$x[,calibrationParameters])
    #y$xd<-(tm$y-tm$x[,calibrationParameters])/tm$x[,15]
  }
}

num.missing<- NA

if(stdArray){
  y$xd<-standardize(y$xd)
}

erScore<-as.vector(t(y$xd["ESR1",]))
her2Score<-as.vector(t(y$xd["ERBB2",]))

# assign the subtype scores and calculate the proliferation score
this.proliferationGenes<-dimnames(y$xd)[[1]] %in% proliferationGenes

prolifScore<-apply(y$xd[this.proliferationGenes,],2,mean,na.rm=T)

##### euclidean distance -------------------------------------------------
out <- sspPredict(pamout.centroids, classes="", y$xd, std=F, 
                  distm="euclidean", centroids=T)
out$distances<- -1*out$distances

# calculate the risk scores
genomic <- 0.04210193*out$distances[,1] + 0.12466938*out$distances[,2] + -0.35235561*out$distances[,3] + 0.14213283*out$distances[,4]
genomicWprolif <- -0.0009299747*out$distances[,1] + 0.0692289192*out$distances[,2] + -0.0951505484*out$distances[,3] +  0.0493487685*out$distances[,4] + 0.3385116381*prolifScore
if(hasClinical){
  xT<-as.numeric(as.vector(y$classes$T))
  combined <- 0.0442770*out$distances[,1] + 0.1170297*out$distances[,2] + -0.2608388*out$distances[,3] + 0.1055908*out$distances[,4] + 0.1813751*xT
  combinedWprolif <- -0.009383416*out$distances[,1] +  0.073725503*out$distances[,2] + -0.090436516*out$distances[,3] + 0.053013865*out$distances[,4] + 0.131605960*xT + 0.327259375*prolifScore
}

# threshold the risk score
griskgroups<-genomic
griskgroups[genomic>ghthreshold]<-"high"
griskgroups[genomic>glthreshold & genomic<ghthreshold]<-"med"
griskgroups[genomic<glthreshold]<-"low"
gpriskgroups<-genomicWprolif
gpriskgroups[genomicWprolif>gphthreshold]<-"high"
gpriskgroups[genomicWprolif>gplthreshold & genomicWprolif<gphthreshold]<-"med"
gpriskgroups[genomicWprolif<gplthreshold]<-"low"

genomic<- 100* (genomic + 0.35 ) / 0.85
genomicWprolif<- 100* (genomicWprolif + 0.35 ) / 0.85

# write output files

outtable<-cbind(out$distances, out$predictions, genomic, griskgroups, prolifScore, genomicWprolif, gpriskgroups, erScore, her2Score)
dimnames(outtable)[[2]]<-c("Basal","Her2","LumA","LumB","Normal","Call",
                           "ROR-S (Subtype Only)","ROR-S Group (Subtype Only)","Proliferation Score", 
                           "ROR-P (Subtype + Proliferation)","ROR-P Group (Subtype + Proliferation)",
                           "ER","Her2")

outFile<- paste(outputDir,paste(short,"_euclidean_pam50scores.txt",sep=""),sep="/")
write.table(outtable,outFile,sep="\t",col.names=NA)

##### MP1 and GSVA --------------------------------------------------------

PAM50 <- read.table("5.1.1.epithelial_NMF/process_k3/Bulk validation/PAM50/output/TCGA_euclidean_pam50scores.txt",
                    header = T, row.names = 1)
PAM50 <- PAM50[,c(1,3,4)]
MP1_PAM50 <- merge(t(AUCell_mtx), PAM50, by = 0) %>% column_to_rownames("Row.names")
MP1_PAM50 <- MP1_PAM50[-grep(".11$|.06$", rownames(MP1_PAM50)),]

for (i in (7 : (ncol(MP1_PAM50)))) {
  
  P <- ggplot(MP1_PAM50, aes(x=MP1_PAM50[,1], y=MP1_PAM50[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = MP1_PAM50[,1], y = MP1_PAM50[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(MP1_PAM50)[1]," AUCell")) + ylab(paste0(colnames(MP1_PAM50)[i]," (-Euc. Dis.)"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/PAM50/CorPlot_MP1_",colnames(MP1_PAM50)[i],".pdf"), P, width = 2, height = 2)
}

##### EO vs LO and GSVA --------------------------------------------------------

Age_PAM50 <- merge(age, PAM50, by = 0) %>% column_to_rownames("Row.names")

# correlation

for (i in (3 : (ncol(Age_PAM50)))) {
  
  P <- ggplot(Age_PAM50, aes(x=Age_PAM50[,1], y=Age_PAM50[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = Age_PAM50[,1], y = Age_PAM50[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_PAM50)[i]," (-Euc. Dis.)"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/PAM50/CorPlot_Age_",colnames(Age_PAM50)[i],".pdf"), P, width = 2, height = 2)
}

for (i in (3 : (ncol(Age_PAM50)))) {
  
  P <- ggplot(Age_PAM50, aes(x=Age_PAM50[,1], y=Age_PAM50[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'pearson', aes(x = Age_PAM50[,1], y = Age_PAM50[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_PAM50)[i]," (-Euc. Dis.)"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/PAM50/CorPlot_pearson_Age_",colnames(Age_PAM50)[i],".pdf"), P, width = 2, height = 2)
}

# boxplot

data <- melt(Age_PAM50[,c(2:ncol(Age_PAM50))])
P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("-Euclidean distance") +
  xlab("") + 
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/Bulk validation/PAM50/BoxPlot_pam50_EOLO.pdf", width = 2, height = 4)
P + NoLegend()
dev.off()

### Bulk mRNAsi (MP1 score / EO vs LO) ----------------------------------

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

#### calculate mRNAsi score ------------------------------
### stemness score based on OCLR mRNAsi

source("./PCBC/predict.mRNAsi.R")
mRNAsi_score <- predict.mRNAsi(expMtx)
saveRDS(mRNAsi_score, "./5.1.1.epithelial_NMF/process_k3/Bulk validation/mRNAsi/mRNAsi_score.Rds")

##### MP1 and mRNAsi --------------------------------------------------------

mRNAsi <- as.data.frame(mRNAsi_score)
MP1_mRNAsi <- merge(t(AUCell_mtx), mRNAsi, by = 0) %>% column_to_rownames("Row.names")
MP1_mRNAsi <- MP1_mRNAsi[-grep(".11$|.06$", rownames(MP1_mRNAsi)),]

P <- ggplot(MP1_mRNAsi, aes(x=MP1_mRNAsi[,1], y=MP1_mRNAsi[,7])) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = MP1_mRNAsi[,1], y = MP1_mRNAsi[,7]), size = 2.5) +
  theme(axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        axis.title.x=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
  xlab("MP1 AUCell") + ylab("mRNAsi")

ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/mRNAsi/CorPlot_MP1_mRNAsi.pdf"), P, width = 2, height = 2)

##### EO vs LO and mRNAsi --------------------------------------------------------

Age_mRNAsi <- merge(age, mRNAsi, by = 0) %>% column_to_rownames("Row.names")

# correlation

i=3

P <- ggplot(Age_mRNAsi, aes(x=Age_mRNAsi[,1], y=Age_mRNAsi[,i])) +
  geom_point(size = 0.5, color = "grey")+ 
  geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
  theme_bw()+
  stat_cor(method = 'spearman', aes(x = Age_mRNAsi[,1], y = Age_mRNAsi[,i]), size = 2.5) +
  theme(axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        axis.title.x=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
  xlab("Age") + ylab("mRNAsi")

ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/mRNAsi/CorPlot_Age_mRNAsi.pdf"), P, width = 2, height = 2)

# boxplot

data <- melt(Age_mRNAsi[,c(2:ncol(Age_mRNAsi))])
P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=6, colour="black"), 
        axis.title.y=element_text(size = 8),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("") +
  ylab("mRNAsi") +
  xlab("") + 
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")),
                     label = "p.signif")
pdf("5.1.1.epithelial_NMF/process_k3/Bulk validation/mRNAsi/BoxPlot_mRANsi_EOLO.pdf", width = 2, height = 2)
P + NoLegend()
dev.off()

### Bulk mPC score (MP1 score / EO vs LO) ----------------------------------

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

#### calculate mPC score ------------------------------

PC_gs <- read.table("ARPC_MSPC_NEPC_gs.txt", sep="\t", header=T)

library(AUCell)
gs <- list()
for (i in c("ARPC",
            "MSPC",
            "NEPC")) {
  gs[[i]] <- PC_gs[which(PC_gs$term == i),]$gene
}

expMtx <- log(expMtx + 1)
cells_rankings <- AUCell_buildRankings(expMtx, nCores=20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gs, cells_rankings, nCores = 20)
saveRDS(cells_AUC, "./5.1.1.epithelial_NMF/process_k3/Bulk validation/mPC/mPC_gs_AUC.Rds")

mPC <- cells_AUC@assays@data$AUC %>% t() %>% as.data.frame()

##### MP1 and mPC --------------------------------------------------------

MP1_mPC <- merge(t(AUCell_mtx), mPC, by = 0) %>% column_to_rownames("Row.names")
MP1_mPC <- MP1_mPC[-grep(".11$|.06$", rownames(MP1_mPC)),]

for (i in (7 : (ncol(MP1_mPC)))) {
  
  P <- ggplot(MP1_mPC, aes(x=MP1_mPC[,1], y=MP1_mPC[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = MP1_mPC[,1], y = MP1_mPC[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab(paste0(colnames(MP1_mPC)[1]," AUCell")) + ylab(paste0(colnames(MP1_mPC)[i]," AUCell"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/mPC/CorPlot_MP1_",colnames(MP1_mPC)[i],".pdf"), P, width = 2, height = 2)
}

##### EO vs LO and mPC --------------------------------------------------------

Age_mPC <- merge(age, mPC, by = 0) %>% column_to_rownames("Row.names")

# correlation

for (i in (3 : (ncol(Age_mPC)))) {
  
  P <- ggplot(Age_mPC, aes(x=Age_mPC[,1], y=Age_mPC[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'spearman', aes(x = Age_mPC[,1], y = Age_mPC[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_mPC)[i]," AUCell"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/mPC/CorPlot_Age_",colnames(Age_mPC)[i],".pdf"), P, width = 2, height = 2)
}

for (i in (3 : (ncol(Age_mPC)))) {
  
  P <- ggplot(Age_mPC, aes(x=Age_mPC[,1], y=Age_mPC[,i])) +
    geom_point(size = 0.5, color = "grey")+ 
    geom_smooth(method="lm", formula = y ~ x, color = "#46732E", size = 0.75) + 
    theme_bw()+
    stat_cor(method = 'pearson', aes(x = Age_mPC[,1], y = Age_mPC[,i]), size = 2.5) +
    theme(axis.text.x=element_text(colour="black", size = 6),
          axis.text.y=element_text(size=6, colour="black"), 
          axis.title.y=element_text(size = 8),
          axis.title.x=element_text(size = 8),
          panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
    xlab("Age") + ylab(paste0(colnames(Age_mPC)[i]," AUCell"))
  
  ggsave(paste0("5.1.1.epithelial_NMF/process_k3/Bulk validation/mPC/CorPlot_pearson_Age_",colnames(Age_mPC)[i],".pdf"), P, width = 2, height = 2)
}

# boxplot

data <- melt(Age_mPC[,c(2:ncol(Age_mPC))])
P = ggplot(data, aes(x = variable, y = value, fill = Type)) +
  scale_fill_manual(values = c("#053E7ACC","#A52A2ACC")) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, colour="black", size = 8),
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
pdf("5.1.1.epithelial_NMF/process_k3/Bulk validation/mPC/BoxPlot_mPC_EOLO.pdf", width = 2, height = 6)
P + NoLegend()
dev.off()

### Bulk scFEA (MP1 EO/LO) --------------------------------------

#### read gene matrix and MP score --------------------
expMtx <- read.table("./0.TCGA/mRNA_Matrix_Count.txt", row.names = 1, header = T, sep = "\t")
expMtx <- expMtx[-grep("^RP[SL]", rownames(expMtx)),]
expMtx <- expMtx[-grep("^MT-", rownames(expMtx)),]
expMtx <- as.matrix(expMtx)
write.csv(expMtx, "5.1.1.epithelial_NMF/process_k3/Bulk validation/scFEA/input/TCGA_expr.csv", quote = F, sep = "\t")

cells_AUC_TCGA <- readRDS("./5.1.1.epithelial_NMF/process_k3/MPScore_AUC_TCGA.Rds")
AUCell_mtx <- as.matrix(cells_AUC_TCGA@assays@data$AUC)

clinical <- read.table("./0.TCGA/Clinical.txt", row.names = 1, header = T, sep = "\t")
age <- as.data.frame(clinical[,1])
colnames(age) <- "age"
rownames(age) <- gsub("_", ".", rownames(clinical))
age$Type[age$age <= 55] <- "EOPC"
age$Type[age$age > 55] <- "LOPC"

#### calculate scFEA score ------------------------------------

python src/scFEA.py --input_dir input --res_dir output --sc_imputation True --test_file TCGA_expr.csv

scFEA <- read.csv("5.1.1.epithelial_NMF/process_k3/Bulk validation/scFEA/output/balance_20230526-154102.csv", row.names =1)
scFEA <- t(scFEA)

##### MP1 H/L ---------------------------------------------

AUCell_mtx <- AUCell_mtx[,-grep(".11$|.06$", colnames(AUCell_mtx))]
MP1 <- AUCell_mtx[1,] %>% as.data.frame()
colnames(MP1) <- "MP1"
med <- median(MP1$MP1)
MP1$MP1_HL <- "MP1_L"
MP1$MP1_HL[MP1$MP1 > med] <- "MP1_H"
type <- MP1[,2]
names(type) <- rownames(MP1)
type <- sort(type)

scFEA <- scFEA[, names(type)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("MP1_H", table(type)[1]), rep("MP1_L", table(type)[2])), levels = c('MP1_H', 'MP1_L'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(scFEA)
head(design)

# Tunor VS Normal
compare <- makeContrasts(MP1_H - MP1_L, levels=design)
fit <- lmFit(scFEA, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# 排序
Diff <- Diff %>% arrange(t)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
dat_plot <- rbind(top_n(dat_plot, 15, wt = t), top_n(dat_plot, -15, wt = t))
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
  ylab('t value of scFEA score, MP1_On versus Off') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black')
ggsave("5.1.1.epithelial_NMF/process_k3/Bulk validation/scFEA/FEA_diff_bar_MP1_top15.pdf", p_MP1, width = 8, height  = 8)

##### MP1 H/L ---------------------------------------------

type <- age[,2]
names(type) <- rownames(age)
type <- sort(type)

scFEA <- scFEA[, match(names(type),colnames(scFEA))]
scFEA <- scFEA[, !is.na(colnames(scFEA))]

type <- type[names(type) %in% colnames(scFEA)]

## limma
library(limma)

# 设置或导入分组
group <- factor(c(rep("EOPC", table(type)[1]), rep("LOPC", table(type)[2])), levels = c('EOPC', 'LOPC'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(scFEA)
head(design)

# Tunor VS Normal
compare <- makeContrasts(EOPC - LOPC, levels=design)
fit <- lmFit(scFEA, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

# 排序
Diff <- Diff %>% arrange(t)

# barplot

dat_plot <- data.frame(id = row.names(Diff), t = Diff$t, adj.P.Val = Diff$adj.P.Val)
dat_plot <- dat_plot[which(dat_plot$adj.P.Val < 0.05),]
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
  guides(fill=F)+ # 不显示图例
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

# 添加标签

p_MP1 <- p + geom_text(data = dat_plot[1:nrow(dat_plot[dat_plot$t < 0,]),], 
                       aes(x = id, y = 0.1, label = id), hjust = 0, color = 'black') +
  geom_text(data = dat_plot[(nrow(dat_plot[dat_plot$t < 0,])+1):nrow(dat_plot),],
            aes(x = id, y = -0.1, label = id), hjust = 1, color = 'black')
ggsave("5.1.1.epithelial_NMF/process_k3/Bulk validation/scFEA/FEA_diff_bar_EO_LO.pdf", p_MP1, width = 8, height  = 8)
