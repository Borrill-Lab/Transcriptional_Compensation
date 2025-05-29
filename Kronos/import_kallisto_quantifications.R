#Delfi Dorussen
#Aim to import the kallisto RNAseq quantifiction files and summarise the data at the gene level (rather than transcript)
#Adapted from https://mbp-tech-talks.github.io/2018-2019/intro-differential-expression/
#and https://github.com/Borrill-Lab/NAM_RNAi_Senescence/blob/main/scripts/02_tximport_summarise_counts_tpm_per_gene.R

install.packages("BiocManager")
library(BiocManager)
BiocManager::install("tximportData", force=TRUE)
BiocManager::install("tximport")
BiocManager::install("rhdf5")
install.packages("readr")
library(tximportData)
library(readr)
library(tximport)
library(rhdf5)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/")

mydir <- getwd()

#Specifying file locations of the kallisto output files
samples <- read.csv("sample_index.csv",header=TRUE)
files <- file.path(mydir,samples$SampleName,"abundance.tsv")
names(files) <- samples$SampleName
all(file.exists(files)) #should be TRUE if all samples are present

#Creating gene to transcript table
url <- "https://raw.githubusercontent.com/Borrill-Lab/WheatFlagLeafSenescence/master/data/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt"
destfile <- file.path(mydir,"transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt")
download.file(url,destfile)
tx2gene <- read.table("transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt")

#Import Kallisto transcript quantifications
txi <- tximport(files, type="kallisto", tx2gene=tx2gene)
names(txi)  #should give abundance, counts, length, countsFromAbundance
head(txi$counts)

#Save count and tpm files
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/Output")
write.table(txi$counts, file="TA_lines_count.tsv")
write.table(txi$abundance, file="TA_lines_tpm.tsv")

#Calculate average gene length across samples
head(txi$length)
gene_lengths <- as.data.frame(rowMeans(txi$length))
#head(gene_lengths)
colnames(gene_lengths) <- c("length")
write.csv(gene_lengths, file="TA_gene_lengths.csv")
