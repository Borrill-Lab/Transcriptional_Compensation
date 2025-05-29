#Delfi Dorussen
#Aim to calculate fold change values in gene expression for comparisons between each of the
#TILLING lines
#This should be run after import_kallisto_quantifications.R

library(BiocManager)
BiocManager::install("DESeq2")
library("DESeq2")
library("dplyr")
library("tidyr")
library(ggplot2)
library(rio)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/Output/")
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
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/")
samples <- read.csv("sample_index.csv")
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/Output")
rownames(samples) <- samples$SampleName
samples <- samples[,c("Genotype","Type","Subgenome")]
samples$Genotype <- factor(samples$Genotype)
samples$Type <- factor(samples$Type)

#Check if metadata is in the correct order - should give output 'TRUE'
all(rownames(samples) == colnames(counts))

#Comparing first batch of samples to KWT1-3
counts_batch1 <- counts %>%
  select(K427_2:K427_15, K2619_1:KWT_3)

samples_batch1 <- samples[c(1:3, 7:21),]
samples_batch1$Genotype <- factor(samples_batch1$Genotype)
samples_batch1$Type <- factor(samples_batch1$Type)

all(rownames(samples_batch1) == colnames(counts_batch1))

##DESeq2 Analysis
dds1 <- DESeqDataSetFromMatrix(countData=round(counts_batch1), colData=samples_batch1, design=~Genotype)
dds1$Type <- relevel(dds1$Genotype,"KWT") #Sets WT as reference level
dds1 <- DESeq(dds1)

#Other genotypes vs WT
res_K427_vs_WT <- results(dds1, contrast=c("Genotype","K427","KWT"))
res_K2619_vs_WT <- results(dds1, contrast=c("Genotype","K2619","KWT"))
res_K2864_vs_WT <- results(dds1, contrast=c("Genotype","K2864","KWT"))
res_K3239_vs_WT <- results(dds1, contrast=c("Genotype","K3239","KWT"))
res_K4533_vs_WT <- results(dds1, contrast=c("Genotype","K4533","KWT"))

#Comparing second batch of samples to KWT4-6
counts_batch2 <- counts %>%
  select(K774_2:K774_25, KWT_4:KWT_6)

samples_batch2 <- samples[c(4:6, 22:24),]
samples_batch2$Genotype <- factor(samples_batch2$Genotype)
samples_batch2$Type <- factor(samples_batch2$Type)

all(rownames(samples_batch2) == colnames(counts_batch2))

##DESeq2 Analysis
dds2 <- DESeqDataSetFromMatrix(countData=round(counts_batch2), colData=samples_batch2, design=~Genotype) #check if it is ok to round counts here!
dds2$Type <- relevel(dds2$Genotype,"KWT") #Sets WT as reference level
dds2 <- DESeq(dds2)

res_K774_vs_WT <- results(dds2, contrast=c("Genotype","K774","KWT"))


contrast_list <- list(res_K427_vs_WT, res_K2619_vs_WT, res_K2864_vs_WT, res_K3239_vs_WT,
                      res_K4533_vs_WT, res_K774_vs_WT)

contrast_list <- lapply(contrast_list, as.data.frame) 

rownames <- rownames(res_K427_vs_WT)
for(i in seq_along(contrast_list)){
  contrast_list[[i]]$Gene <- rownames
}

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/Output/DESeq2_All_Contrasts/")
export_list(contrast_list, c("res_K427_vs_WT.csv", "res_K2619_vs_WT.csv", "res_K2864_vs_WT.csv", "res_K3239_vs_WT.csv",
                              "res_K4533_vs_WT.csv", "res_K774_vs_WT.csv"))



#Looking at specific genes of interest
counts_dds1 <- counts(dds1, normalized=TRUE)
counts_dds1_PHS1 <- counts_dds1[row.names(counts_dds1) %in% c("TraesCS5A02G395200", "TraesCS5B02G400000"),]
counts_dds1_PHS1 <- as.data.frame(counts_dds1_PHS1)
counts_dds1_PHS1$gene <- rownames(counts_dds1_PHS1)
counts_dds1_PHS1 <- counts_dds1_PHS1 %>%
  pivot_longer(cols=c(1:18), names_to="Genotype") %>%
  separate_wider_delim(col="Genotype", names=c("Genotype","Rep"), delim="_")

PHS1_summary <- counts_dds1_PHS1 %>%
  group_by(gene, Genotype) %>%
  summarise(mean=mean(value), SE=sd(value)/sqrt(3))

PHS1_summary_min <- subset(PHS1_summary, Genotype=="K4533"|Genotype=="K2864"|Genotype=="KWT")
counts_dds1_PHS1_min <- subset(counts_dds1_PHS1, Genotype=="K4533"|Genotype=="K2864"|Genotype=="KWT")

colours <- c("#ABABAB", "#58a643", "#4a2976")

PHS1_summary_min$Genotype <- factor(PHS1_summary_min$Genotype, levels=c("KWT","K4533","K2864"))

PHS1_bar <- ggplot(PHS1_summary_min, aes(x=Genotype, y=mean, fill=Genotype)) +
  geom_bar(stat="identity", alpha=0.8, colour="black") +
  geom_jitter(data=counts_dds1_PHS1_min, aes(x=Genotype, y=value), width=0.1, size=2) +
  scale_fill_manual(values=colours) +
  facet_grid(rows=NULL, cols=vars(gene)) +
  geom_errorbar(aes(ymin=mean-(1.96*SE), ymax=mean+(1.96*SE)), width=.2) +
  theme_bw() +
  theme(legend.position="bottom") +
  ylab("Normalised Counts")
ggsave("PHS1_normalised_counts.pdf", width=5, height=5)


counts_dds1_CSP41a <- counts_dds1[row.names(counts_dds1) %in% c("TraesCS6A02G025700", "TraesCS6B02G036400"),]
counts_dds1_CSP41a <- as.data.frame(counts_dds1_CSP41a)
counts_dds1_CSP41a$gene <- rownames(counts_dds1_CSP41a)
counts_dds1_CSP41a <- counts_dds1_CSP41a %>%
  pivot_longer(cols=c(1:18), names_to="Genotype") %>%
  separate_wider_delim(col="Genotype", names=c("Genotype","Rep"), delim="_")

CSP41a_summary <- counts_dds1_CSP41a %>%
  group_by(gene, Genotype) %>%
  summarise(mean=mean(value), SE=sd(value)/sqrt(3))

CSP41a_summary_min <- subset(CSP41a_summary, Genotype=="K3239"|Genotype=="K2619"|Genotype=="KWT")
counts_dds1_CSP41a_min <- subset(counts_dds1_CSP41a, Genotype=="K3239"|Genotype=="K2619"|Genotype=="KWT")

CSP41a_summary_min$Genotype <- factor(CSP41a_summary_min$Genotype, levels=c("KWT","K3239","K2619"))

CSP41a_bar <- ggplot(CSP41a_summary_min, aes(x=Genotype, y=mean, fill=Genotype)) +
  geom_bar(stat="identity", alpha=0.8, colour="black") +
  geom_jitter(data=counts_dds1_CSP41a_min, aes(x=Genotype, y=value), width=0.1, size=2) +
  scale_fill_manual(values=colours) +
  facet_grid(rows=NULL, cols=vars(gene)) +
  geom_errorbar(aes(ymin=mean-(1.96*SE), ymax=mean+(1.96*SE)), width=.2) +
  theme_bw() +
  theme(legend.position="bottom") +
  ylab("Normalised Counts")
ggsave("CSP41a_normalised_counts.pdf", width=5, height=5)



