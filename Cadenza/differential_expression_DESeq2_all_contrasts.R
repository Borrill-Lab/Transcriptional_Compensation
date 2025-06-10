#Delfi Dorussen
#Aim to calculate fold change values in gene expression for comparisons between each of the
#TILLING lines
#This should be run after import_kallisto_quantifications.R

BiocManager::install("DESeq2")
library("DESeq2")
library("dplyr")
library("tidyr")
library(ggplot2)
library(rio)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Output")
abundances <- read.table("TA_lines_tpm.tsv")
head(abundances)

#Remove low confidence (LC) genes
abundances <- cbind(rownames(abundances), data.frame(abundances, row.names=NULL))
colnames(abundances)[1] <- "gene"
head(abundances)

abundances <- abundances[!grepl("LC",abundances$gene),]

#Importing count data and subsetting in the same way as the tpm data
counts <- read.table("TA_lines_count.tsv")
head(counts)

counts <- counts[rownames(counts) %in% abundances$gene,]

#Providing metadata to the count matrices
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project")
samples <- read.csv("sample_index.csv")
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Output")
rownames(samples) <- samples$SampleName
samples <- samples[,c("Genotype","Type")]
samples$Genotype <- factor(samples$Genotype)
samples$Type <- factor(samples$Type)

#Check if metadata is in the correct order - should give output 'TRUE'
all(rownames(samples) == colnames(counts))

##DESeq2 Analysis
dds1 <- DESeqDataSetFromMatrix(countData=round(counts), colData=samples, design=~Genotype)
dds1$Type <- relevel(dds1$Genotype,"CWT") #Sets WT as reference level
dds1 <- DESeq(dds1)

#WT vs other genotypes
res_WT_vs_C0604 <- results(dds1, contrast=c("Genotype","CWT","C0604"))
res_WT_vs_C0895 <- results(dds1, contrast=c("Genotype","CWT","C0895"))
res_WT_vs_C1015 <- results(dds1, contrast=c("Genotype","CWT","C1015"))
res_WT_vs_C1704 <- results(dds1, contrast=c("Genotype","CWT","C1704"))

#Other genotypes vs WT
res_C0604_vs_WT <- results(dds1, contrast=c("Genotype","C0604","CWT"))
res_C0895_vs_WT <- results(dds1, contrast=c("Genotype","C0895","CWT"))
res_C1015_vs_WT <- results(dds1, contrast=c("Genotype","C1015","CWT"))
res_C1704_vs_WT <- results(dds1, contrast=c("Genotype","C1704","CWT"))

#C0604 vs other genotypes
res_C0604_vs_c0895 <- results(dds1, contrast=c("Genotype","C0604","C0895"))
res_C0604_vs_C1015 <- results(dds1, contrast=c("Genotype","C0604","C1015"))
res_C0604_vs_C1704 <- results(dds1, contrast=c("Genotype","C0604","C1704"))

#C0895 vs other genotypes
res_C0895_vs_C0604 <- results(dds1, contrast=c("Genotype","C0895","C0604"))
res_C0895_vs_C1015 <- results(dds1, contrast=c("Genotype","C0895","C1015"))
res_C0895_vs_C1704 <- results(dds1, contrast=c("Genotype","C0895","C1704"))

#C1015 vs other genotypes
res_C1015_vs_C0604 <- results(dds1, contrast=c("Genotype","C1015","C0604"))
res_C1015_vs_C0895 <- results(dds1, contrast=c("Genotype","C1015","C0895"))
res_C1015_vs_C1704 <- results(dds1, contrast=c("Genotype","C1015","C1704"))

#C1704 vs other genotypes
res_C1704_vs_C0604 <- results(dds1, contrast=c("Genotype","C1704","C0604"))
res_C1704_vs_C0895 <- results(dds1, contrast=c("Genotype","C1704","C0895"))
res_C1704_vs_C1015 <- results(dds1, contrast=c("Genotype","C1704","C1015"))

contrast_list <- list(res_WT_vs_C0604, res_WT_vs_C0895, res_WT_vs_C1015, res_WT_vs_C1704,
               res_C0604_vs_WT, res_C0895_vs_WT, res_C1015_vs_WT, res_C1704_vs_WT,
               res_C0604_vs_c0895, res_C0604_vs_C1015, res_C0604_vs_C1704,
               res_C0895_vs_C0604, res_C0895_vs_C1015, res_C0895_vs_C1704,
               res_C1015_vs_C0604, res_C1015_vs_C0895, res_C1015_vs_C1704,
               res_C1704_vs_C0604, res_C1704_vs_C0895, res_C1704_vs_C1015)

contrast_list <- lapply(contrast_list, as.data.frame) 

rownames <- rownames(res_WT_vs_C0604)
for(i in seq_along(contrast_list)){
  contrast_list[[i]]$Gene <- rownames
}

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Output/DESeq2_All_Contrasts/")
export_list(contrast_list, c("res_WT_vs_C0604.csv", "res_WT_vs_C0895.csv", "res_WT_vs_C1015.csv", "res_WT_vs_C1704.csv",
                              "res_C0604_vs_WT.csv", "res_C0895_vs_WT.csv", "res_C1015_vs_WT.csv", "res_C1704_vs_WT.csv",
                              "res_C0604_vs_c0895.csv", "res_C0604_vs_C1015.csv", "res_C0604_vs_C1704.csv",
                              "res_C0895_vs_C0604.csv", "res_C0895_vs_C1015.csv", "res_C0895_vs_C1704.csv",
                              "res_C1015_vs_C0604.csv", "res_C1015_vs_C0895.csv", "res_C1015_vs_C1704.csv",
                              "res_C1704_vs_C0604.csv", "res_C1704_vs_C0895.csv", "res_C1704_vs_C1015.csv"))
