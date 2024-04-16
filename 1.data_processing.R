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
library(ggsci)
library(gridExtra)

# Data processing------------------------------------------------------------------

# 1.load data ------------------------------------------------------------------

folder <- "1.data processing"
if(!dir.exists(folder)){
  dir.create(folder)
}

for (file in c("PCSC002", "PCSC003", "PCSC005", "PCSC007", "PCSC008", "PCSC009", "PCSC010", "PCSC012", "PCSC013", "PCSC015")) {
  data <- Read10X(paste0(file))
  assign(file, data)
}

# 2.averege for duplicates ----------------------------------------------------------------

length(unique(rownames(PCSC002)))-length(rownames(PCSC002))

# 3.add labels for samples and merge data ----------------------------------------------------------------

colnames(PCSC002) <- substr(colnames(PCSC002), 1, nchar(colnames(PCSC002))-2)
colnames(PCSC002) <- paste0("EOPC1_", colnames(PCSC002), "_002")

colnames(PCSC003) <- substr(colnames(PCSC003), 1, nchar(colnames(PCSC003))-2)
colnames(PCSC003) <- paste0("LOPC1_", colnames(PCSC003), "_003")

colnames(PCSC005) <- substr(colnames(PCSC005), 1, nchar(colnames(PCSC005))-2)
colnames(PCSC005) <- paste0("EOPC2_", colnames(PCSC005), "_005")

colnames(PCSC007) <- paste0("EOPC3_", colnames(PCSC007), "_007")

colnames(PCSC008) <- paste0("EOPC4_", colnames(PCSC008), "_008")

colnames(PCSC009) <- paste0("LOPC2_", colnames(PCSC009), "_009")

colnames(PCSC010) <- paste0("LOPC3_", colnames(PCSC010), "_010")

colnames(PCSC012) <- paste0("LOPC4_", colnames(PCSC012), "_012")

colnames(PCSC013) <- paste0("LOPC5_", colnames(PCSC013), "_013")

colnames(PCSC015) <- paste0("LOPC6_", colnames(PCSC015), "_015")

merge <- cbind(PCSC002, PCSC003, PCSC005, PCSC007, PCSC008, PCSC009, PCSC010, PCSC012, PCSC013, PCSC015)
dim(merge)
length(grep("_008", colnames(merge)))

# 4.create Seurat objects based on merge data for QC ----------------------------------------------------------

PCSC_preCCA <- CreateSeuratObject(
               merge,
               project = "PCSC_preCCA", 
               min.cells = 10,
               min.features = 200,
               names.field = 1,
               names.delim = "_")

dim(PCSC_preCCA)
View(PCSC_preCCA@meta.data)
saveRDS(PCSC_preCCA, file = "./1.data_processing/PCSC_preCCA.Rds")
rm(data, file, PCSC002, PCSC003, PCSC005, PCSC007, PCSC008, PCSC009, PCSC010, PCSC012, PCSC013, PCSC015, merge)

PCSC_preCCA <- readRDS("./1.data_processing/PCSC_preCCA.Rds")

# 5.QC ----------------------------------------------------------------------------

folder <- "1.data processing/QC"
if(!dir.exists(folder)){
  dir.create(folder)
}

# 5.1.cell counts --------------------------------------------------------------------

before_sep <- as.data.frame(table(PCSC_preCCA@meta.data$orig.ident))

Ncell_before <- PCSC_preCCA@meta.data %>%
  ggplot(aes(x=orig.ident, fill=orig.ident)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_fill_npg() +
  ggtitle("NCells") +
  theme(legend.position = "none")

# 5.2.nCount_RNA ----------------------------------------------------------------

count_den <- PCSC_preCCA@meta.data %>%
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  scale_fill_npg() +
  geom_vline(xintercept = 1000) +
  theme(legend.position = "none")

# 5.3.nFeature_RNA----------------------------------------------------------

feature_den <- PCSC_preCCA@meta.data %>%
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ylab("Cell density") +
  scale_x_log10() +
  scale_fill_npg() +
  geom_vline(xintercept = 500) +
  theme(legend.position = "none")

# 5.4. mito percentage, nCount vs. nFeature ----------------------------------------------------

PCSC_preCCA[["percent.mt"]] <- PercentageFeatureSet(PCSC_preCCA,pattern = "^MT-")

mt_den <- PCSC_preCCA@meta.data %>%
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ylab("Cell density")  +
  scale_fill_npg() +
  geom_vline(xintercept = 15) +
  theme(legend.position = "none")

pdf("./1.data_processing/QC/QC_nFeature_nCount.pdf")
PCSC_preCCA@meta.data %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 300) +
  geom_hline(yintercept = 500) +
  scale_fill_npg() +
  facet_wrap(~orig.ident)
dev.off()

# 5.5.number of genes detected per UMI ----------------------------------------------------------------

PCSC_preCCA$log10GenesPerUMI <- log10(PCSC_preCCA$nFeature_RNA)/log10(PCSC_preCCA$nCount_RNA)

perUMI_den <- PCSC_preCCA@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ylab("Cell density") +
  scale_fill_npg() +
  geom_vline(xintercept = 0.8) +
  theme(legend.position = "none")

pdf("./1.data_processing/QC/QC_density.pdf",width = 11,height = 4.5)  
grid.arrange(count_den,feature_den,mt_den,perUMI_den,nrow=1,ncol=4)
dev.off() 

pdf("./1.data_processing/QC/density_legend.pdf")
PCSC_preCCA@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ylab("Cell density") +
  scale_fill_npg() +
  geom_vline(xintercept = 0.8)
dev.off()

# 5.6.VlnPlot ------------------------------------------------------------
P1 <- PCSC_preCCA@meta.data %>%
  ggplot(aes(x=orig.ident, y=nFeature_RNA, fill=orig.ident)) +
  geom_violin() +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature") +
  theme(legend.position = "none") 

P2 <- PCSC_preCCA@meta.data %>%
  ggplot(aes(x=orig.ident, y=nCount_RNA, fill=orig.ident)) +
  geom_violin() +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount") +
  theme(legend.position = "none") 

P3 <- PCSC_preCCA@meta.data %>%
  ggplot(aes(x=orig.ident, y=percent.mt, fill=orig.ident)) +
  geom_violin() +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("percent.mt") +
  theme(legend.position = "none") 

P4 <- PCSC_preCCA@meta.data %>%
  ggplot(aes(x=orig.ident, y=log10GenesPerUMI, fill=orig.ident)) +
  geom_violin() +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("log10GenesPerUMI") +
  theme(legend.position = "none")

pdf("./1.data_processing/QC/QC_VlnPlot.pdf",width = 11,height = 4.5)  
grid.arrange(P1,P2,P3,P4,nrow=1,ncol=4)
dev.off() 

pdf("./1.data_processing/QC/Vln_legend.pdf")
PCSC_preCCA@meta.data %>%
  ggplot(aes(x=orig.ident, y=log10GenesPerUMI, fill=orig.ident)) +
  geom_violin() +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("log10GenesPerUMI")
dev.off()

# 5.6.distribution of nFeature, nCount, mito percentage and log10GenesPerUMI ---------------------------------------------------

# per 5 percentage
gene.freq <- do.call("cbind", tapply(PCSC_preCCA@meta.data$nFeature_RNA,PCSC_preCCA@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(PCSC_preCCA@meta.data$nCount_RNA,PCSC_preCCA@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(PCSC_preCCA@meta.data$percent.mt,PCSC_preCCA@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
GPU.freq <- do.call("cbind", tapply(PCSC_preCCA@meta.data$log10GenesPerUMI,PCSC_preCCA@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq,GPU.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"),
                            paste(colnames(GPU.freq),"Gene","per","UMI",sep = "_"))
write.table(freq.combine,file = "./1.data_processing/QC/QC_freq.txt",quote = F,sep = "\t")
rm(gene.freq,rna.freq,mt.freq)
View(freq.combine)

# mito percentage vs. nCount & nGene vs. nCount
 
plot1 <- FeatureScatter(PCSC_preCCA, feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T, pt.size = 0.5) + 
  scale_color_npg() +
  NoLegend()
plot2 <- FeatureScatter(PCSC_preCCA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle = T, pt.size = 0.5) + 
  scale_color_npg() +
  NoLegend()
pdf("./1.data_processing/QC/QC-FeatureScatter.pdf", width = 6, height = 4.5)
grid.arrange(plot1,plot2,nrow=1,ncol=4)
dev.off()
rm(plot1,plot2)

pdf("./1.data_processing/QC/FeatureScatter_legend.pdf")
plot1 <- FeatureScatter(PCSC_preCCA, feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T, pt.size = 0.5) + 
  scale_color_npg()
dev.off()

# 6.filter cells and genes

# 6.1.filter cells ----------------------------------------------------------------

before <- length(rownames(PCSC_preCCA@meta.data))

PCSC_preCCA <- subset(PCSC_preCCA, 
                      subset = 
                        nFeature_RNA > 500 & 
                        nCount_RNA > 1000 &
                        percent.mt < 15 &
                        log10GenesPerUMI > 0.80)

after <- length(rownames(PCSC_preCCA@meta.data))

after_sep <- as.data.frame(table(PCSC_preCCA@meta.data$orig.ident))
change_sep <- rbind(before_sep, after_sep)
change_sep$stat <- rep(c("A", "B"), each = 10)
colnames(change_sep) <- c("Sample", "nCells", "stat")

pdf("./1.data_processing/QC/Cell_change.pdf")
change_sep %>%
  ggplot(aes(x=Sample, y=nCells,fill=stat)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_npg()
dev.off()

saveRDS(PCSC_preCCA, "./1.data_processing/QC/PCSC_QC.Rds")

change <- matrix(c(before,after))
rownames(change) <- c("before", "after")
write.table(change, "./1.data_processing/QC/CellChange.txt", quote = F, sep = "\t", col.names = F)

rm(change, before, after, freq.combine)

Ncell_after <- PCSC_preCCA@meta.data %>%
  ggplot(aes(x=orig.ident, fill=orig.ident)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_fill_npg() +
  ggtitle("NCells") +
  theme(legend.position = "none")

# 6.2.filter genes

counts <- GetAssayData(object = PCSC_preCCA, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
PCSC_preCCA <- CreateSeuratObject(filtered_counts, meta.data = PCSC_preCCA@meta.data)
saveRDS(PCSC_preCCA, "./1.data_processing/QC/PCSC_QC.Rds")

# 7.create objects for each sample -----------------------------------------------------------------

#for (i in c("002", "003", "005", "007", "008", "009")) {
#  object <- subset(PCSC_preCCA, orig.ident == i)
#  assign(paste0("Seurat", i), object)
#}

#saveRDS(Seurat002, "./1.data processing/Seurat002.Rds")
#saveRDS(Seurat003, "./1.data processing/Seurat003.Rds")
#saveRDS(Seurat005, "./1.data processing/Seurat005.Rds")
#saveRDS(Seurat007, "./1.data processing/Seurat007.Rds")
#saveRDS(Seurat008, "./1.data processing/Seurat008.Rds")
#saveRDS(Seurat009, "./1.data processing/Seurat009.Rds")
