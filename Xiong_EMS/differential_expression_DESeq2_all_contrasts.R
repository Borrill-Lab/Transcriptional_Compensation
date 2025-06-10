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

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/Output")
abundances <- read.table("TA_lines_tpm.tsv")
head(abundances)

#Remove low confidence (LC) genes - could add additional filtering here
abundances <- cbind(rownames(abundances), data.frame(abundances, row.names=NULL))
colnames(abundances)[1] <- "gene"
head(abundances)

abundances <- abundances[!grepl("LC",abundances$gene),]

#Importing count data and subsetting in the same way as the tpm data
counts <- read.table("TA_lines_count.tsv")
head(counts)

counts <- counts[rownames(counts) %in% abundances$gene,]

#Providing metadata to the count matrices
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/")
samples <- read.csv("sample_index.csv")
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/Output")
rownames(samples) <- samples$SampleName
samples <- samples[,c("Genotype","Type")]
samples$Genotype <- factor(samples$Genotype)
samples$Type <- factor(samples$Type)

#Check if metadata is in the correct order - these commands should give output 'TRUE'
all(rownames(samples) == colnames(counts))
rownames(samples) <- colnames(counts)

##DESeq2 Analysis
dds1 <- DESeqDataSetFromMatrix(countData=round(counts), colData=samples, design=~Genotype) #check if it is ok to round counts here!
dds1$Type <- relevel(dds1$Genotype,"WT") #Sets WT as reference level
dds1 <- DESeq(dds1)

#WT vs other genotypes
res_WT_vs_dm3 <- results(dds1, contrast=c("Genotype","WT","dm3"))
res_WT_vs_dm4 <- results(dds1, contrast=c("Genotype","WT","dm4"))

#Other genotypes vs WT
res_dm3_vs_WT <- results(dds1, contrast=c("Genotype","dm3","WT"))
res_dm4_vs_WT <- results(dds1, contrast=c("Genotype","dm4","WT"))

#dm3 vs other genotypes
res_dm3_vs_dm4 <- results(dds1, contrast=c("Genotype","dm3","dm4"))

#dm4 vs other genotypes
res_dm4_vs_dm3 <- results(dds1, contrast=c("Genotype","dm4","dm3"))

contrast_list <- list(res_WT_vs_dm3, res_WT_vs_dm4,
               res_dm3_vs_WT, res_dm4_vs_WT,
               res_dm3_vs_dm4, res_dm4_vs_dm3)

contrast_list <- lapply(contrast_list, as.data.frame) 

rownames <- rownames(res_WT_vs_dm3)
for(i in seq_along(contrast_list)){
  contrast_list[[i]]$Gene <- rownames
}

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/Output/DESeq2_All_Contrasts/")
export_list(contrast_list, c("res_WT_vs_dm3.csv", "res_WT_vs_dm4.csv",
                              "res_dm3_vs_WT.csv", "res_dm4_vs_WT.csv",
                              "res_dm3_vs_dm4.csv", "res_dm4_vs_dm3.csv"))
