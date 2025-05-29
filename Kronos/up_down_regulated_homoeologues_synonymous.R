#Delfi Dorussen
#Aim to combine variant information and differential expression information to identify
#potential homoeologous transcriptional compensation

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/Output/")

#Load in syn_variants_diads_triads and subset abundances to only keep triads/diads in this list
syn_diads_triads <- read.csv("../syn_variants_diads_triads.csv")
head(syn_diads_triads)
syn_diads_triads <- select(syn_diads_triads, X:B)

#Pivot syn_triads_diads to have each homoeologue as a separate row
syn_diads_triads <- pivot_longer(syn_diads_triads, 24:25, names_to="Subgenome", values_to="Homoeologues")
syn_diads_triads <- subset(syn_diads_triads, is.na(Homoeologues)==FALSE)

#Specify whether the homoeologue has the PTC variant
syn_diads_triads$Mutant <- NA
for(i in 1:nrow(syn_diads_triads)){
  print(i)
  if(syn_diads_triads$GENE[i]==syn_diads_triads$Homoeologues[i]){
    syn_diads_triads$Mutant[i] <- "YES"
  } else{
    syn_diads_triads$Mutant[i] <- "NO"
  }
}

#Subset variants to keep syns that are in a TILLING line, but not KWT
syn_diads_triads_unique_K2619 <- subset(syn_diads_triads, K2619=="1/1"&KWT=="0/0")
syn_diads_triads_unique_K2864 <- subset(syn_diads_triads, K2864=="1/1"&KWT=="0/0")
syn_diads_triads_unique_K3239 <- subset(syn_diads_triads, K3239=="1/1"&KWT=="0/0")
syn_diads_triads_unique_K427 <- subset(syn_diads_triads, K427=="1/1"&KWT=="0/0")
syn_diads_triads_unique_K4533 <- subset(syn_diads_triads, K4533=="1/1"&KWT=="0/0")
syn_diads_triads_unique_K774 <- subset(syn_diads_triads, K774=="1/1"&KWT=="0/0")

#Load in contrasts to WT
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/Output/DESeq2_All_Contrasts/")
res_K2619_vs_CWT <- read.csv("res_K2619_vs_WT.csv")
res_K2864_vs_CWT <- read.csv("res_K2864_vs_WT.csv")
res_K3239_vs_CWT <- read.csv("res_K3239_vs_WT.csv")
res_K427_vs_CWT <- read.csv("res_K427_vs_WT.csv")
res_K4533_vs_CWT <- read.csv("res_K4533_vs_WT.csv")
res_K774_vs_CWT <- read.csv("res_K774_vs_WT.csv")

#Merge FC data with variant data
syn_diads_triads_unique_K2619 <- merge(syn_diads_triads_unique_K2619, res_K2619_vs_CWT, by.x="Homoeologues", by.y="Gene")
syn_diads_triads_unique_K2864 <- merge(syn_diads_triads_unique_K2864, res_K2864_vs_CWT, by.x="Homoeologues", by.y="Gene")
syn_diads_triads_unique_K3239 <- merge(syn_diads_triads_unique_K3239, res_K3239_vs_CWT, by.x="Homoeologues", by.y="Gene")
syn_diads_triads_unique_K427 <- merge(syn_diads_triads_unique_K427, res_K427_vs_CWT, by.x="Homoeologues", by.y="Gene")
syn_diads_triads_unique_K4533 <- merge(syn_diads_triads_unique_K4533, res_K4533_vs_CWT, by.x="Homoeologues", by.y="Gene")
syn_diads_triads_unique_K774 <- merge(syn_diads_triads_unique_K774, res_K774_vs_CWT, by.x="Homoeologues", by.y="Gene")

#Getting number of homoeolog groups affected by a PTC in each TILLING line
syn_diads_triads_unique_K2619$group_num <- factor(syn_diads_triads_unique_K2619$group_num)
nlevels(syn_diads_triads_unique_K2619$group_num) #591

syn_diads_triads_unique_K2864$group_num <- factor(syn_diads_triads_unique_K2864$group_num)
nlevels(syn_diads_triads_unique_K2864$group_num) #693

syn_diads_triads_unique_K3239$group_num <- factor(syn_diads_triads_unique_K3239$group_num)
nlevels(syn_diads_triads_unique_K3239$group_num) #213

syn_diads_triads_unique_K427$group_num <- factor(syn_diads_triads_unique_K427$group_num)
nlevels(syn_diads_triads_unique_K427$group_num) #239

syn_diads_triads_unique_K4533$group_num <- factor(syn_diads_triads_unique_K4533$group_num)
nlevels(syn_diads_triads_unique_K4533$group_num) #290

syn_diads_triads_unique_K774$group_num <- factor(syn_diads_triads_unique_K774$group_num)
nlevels(syn_diads_triads_unique_K774$group_num) #576

#Add column to specify which comparison is being made
syn_diads_triads_unique_K2619$Comparison <- "K2619vsWT"
syn_diads_triads_unique_K2864$Comparison <- "K2864vsWT"
syn_diads_triads_unique_K3239$Comparison <- "K3239vsWT"
syn_diads_triads_unique_K427$Comparison <- "K427vsWT"
syn_diads_triads_unique_K4533$Comparison <- "K4533vsWT"
syn_diads_triads_unique_K774$Comparison <- "K774vsWT"

#Combine into one dataframe
syn_diads_triads_unique <- rbind(syn_diads_triads_unique_K2619, syn_diads_triads_unique_K2864,
                                 syn_diads_triads_unique_K3239, syn_diads_triads_unique_K427,
                                 syn_diads_triads_unique_K4533, syn_diads_triads_unique_K774)
syn_diads_triads_unique$group_num <- factor(syn_diads_triads_unique$group_num)
nlevels(syn_diads_triads_unique$group_num) #1783 - so quite a few homoeolog groups must be shared amongst the TILLING lines

#Indicate whether the gene is up or downregulated, or no significant change
syn_diads_triads_unique$Expression <- ifelse(syn_diads_triads_unique$log2FoldChange>0, "UP", "DOWN")
syn_diads_triads_unique$Expression <- ifelse(syn_diads_triads_unique$padj>0.05, "NO CHANGE", syn_diads_triads_unique$Expression)
syn_diads_triads_unique <- subset(syn_diads_triads_unique, is.na(Expression)==FALSE)

#Make plots
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/Output")

summary_syn <- syn_diads_triads_unique %>%
  group_by(Comparison, Mutant, Expression) %>%
  summarise(count=n())

summary_syn$Mutant <- factor(summary_syn$Mutant, levels=c("YES","NO"))

barplot <- ggplot(summary_syn, aes(x=Mutant, y=count, fill=Expression)) +
  geom_bar(position="fill", stat="identity", alpha=0.8, colour="black") +
  theme_bw() +
  scale_fill_manual(values=c("#F0E442","#ABABAB","#009E73"), labels=c("Down", "No Change", "Up")) +
  scale_x_discrete(labels=c("Genes with Syn","Homoeologues")) +
  xlab(NULL) +
  ylab("Proportion of Genes") +
  labs(fill="Expression Relative to WT") +
  facet_grid(cols=vars(Comparison)) +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1))
ggsave("up_down_regulated_homoeologues_stacked_padj0.05_syn.pdf", width=8.5, height=5, units="in", dpi=300)

#Count number of dyads with upregulation
up_homoeologs <- subset(syn_diads_triads_unique, Mutant=="NO"&Expression=="UP")
up_homoeologs$group_num <- factor(up_homoeologs$group_num)
nlevels(up_homoeologs$group_num) #28
