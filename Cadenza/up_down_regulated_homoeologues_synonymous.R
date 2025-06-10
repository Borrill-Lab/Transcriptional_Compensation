#Delfi Dorussen
#Aim to combine variant information and differential expression information to identify
#potential homoeologous transcriptional compensation

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Output")

#Load in syn_variants_diads_triads and subset abundances to only keep triads/diads in this list
syn_diads_triads <- read.csv("../synonymous_variants_diads_triads.csv")
head(syn_diads_triads)

#Pivot syn_triads_diads to have each homoeologue as a separate row
syn_diads_triads <- pivot_longer(syn_diads_triads, 23:25, names_to="Subgenome_Hom", values_to="Homoeologues")
syn_diads_triads <- subset(syn_diads_triads, is.na(Homoeologues)==FALSE)

#Specify whether the homoeologue has the syn variant
syn_diads_triads$Mutant <- NA
for(i in 1:nrow(syn_diads_triads)){
  print(i)
  if(syn_diads_triads$GENE[i]==syn_diads_triads$Homoeologues[i]){
    syn_diads_triads$Mutant[i] <- "YES"
  } else{
    syn_diads_triads$Mutant[i] <- "NO"
  }
}


#Subset variants to keep syns that are only in one line
syn_diads_triads_unique_C0604 <- subset(syn_diads_triads, C0604=="1/1"&CWT=="0/0")
syn_diads_triads_unique_C0895 <- subset(syn_diads_triads, C0895=="1/1"&CWT=="0/0")
syn_diads_triads_unique_C1015 <- subset(syn_diads_triads, C1015=="1/1"&CWT=="0/0")
syn_diads_triads_unique_C1704 <- subset(syn_diads_triads, C1704=="1/1"&CWT=="0/0")

#Load in contrasts to WT
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Output/DESeq2_All_Contrasts/")
res_C0604_vs_CWT <- read.csv("res_C0604_vs_WT.csv")
res_C0895_vs_CWT <- read.csv("res_C0895_vs_WT.csv")
res_C1015_vs_CWT <- read.csv("res_C1015_vs_WT.csv")
res_C1704_vs_CWT <- read.csv("res_C1704_vs_WT.csv")

#Merge FC data with variant data
syn_diads_triads_unique_C0604 <- merge(syn_diads_triads_unique_C0604, res_C0604_vs_CWT, by.x="Homoeologues", by.y="Gene")
syn_diads_triads_unique_C0895 <- merge(syn_diads_triads_unique_C0895, res_C0895_vs_CWT, by.x="Homoeologues", by.y="Gene")
syn_diads_triads_unique_C1015 <- merge(syn_diads_triads_unique_C1015, res_C1015_vs_CWT, by.x="Homoeologues", by.y="Gene")
syn_diads_triads_unique_C1704 <- merge(syn_diads_triads_unique_C1704, res_C1704_vs_CWT, by.x="Homoeologues", by.y="Gene")

#Getting number of homoeolog groups affected by a syn in each TILLING line
syn_diads_triads_unique_C0604$group_num <- factor(syn_diads_triads_unique_C0604$group_num)
nlevels(syn_diads_triads_unique_C0604$group_num) #643

syn_diads_triads_unique_C0895$group_num <- factor(syn_diads_triads_unique_C0895$group_num)
nlevels(syn_diads_triads_unique_C0895$group_num) #598

syn_diads_triads_unique_C1015$group_num <- factor(syn_diads_triads_unique_C1015$group_num)
nlevels(syn_diads_triads_unique_C1015$group_num) #1093

syn_diads_triads_unique_C1704$group_num <- factor(syn_diads_triads_unique_C1704$group_num)
nlevels(syn_diads_triads_unique_C1704$group_num) #627

#Remove duplicates due to having more than one homoeolog with a variant from the table 
syn_C0604_YES <- subset(syn_diads_triads_unique_C0604, Mutant=="YES")
syn_C0604_NO <- subset(syn_diads_triads_unique_C0604, Mutant=="NO")

syn_C0604_NO_dups_removed <- subset(syn_C0604_NO, !Homoeologues%in%syn_C0604_YES$Homoeologues)
syn_C0604_NO_dups_removed <- syn_C0604_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_C0604_dups_removed <- rbind(syn_C0604_YES, syn_C0604_NO_dups_removed)
nlevels(as.factor(syn_diads_triads_unique_C0604_dups_removed$group_num)) #643

syn_C0895_YES <- subset(syn_diads_triads_unique_C0895, Mutant=="YES")
syn_C0895_NO <- subset(syn_diads_triads_unique_C0895, Mutant=="NO")

syn_C0895_NO_dups_removed <- subset(syn_C0895_NO, !Homoeologues%in%syn_C0895_YES$Homoeologues)
syn_C0895_NO_dups_removed <- syn_C0895_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_C0895_dups_removed <- rbind(syn_C0895_YES, syn_C0895_NO_dups_removed)
nlevels(as.factor(syn_diads_triads_unique_C0895_dups_removed$group_num)) #598

syn_C1015_YES <- subset(syn_diads_triads_unique_C1015, Mutant=="YES")
syn_C1015_NO <- subset(syn_diads_triads_unique_C1015, Mutant=="NO")

syn_C1015_NO_dups_removed <- subset(syn_C1015_NO, !Homoeologues%in%syn_C1015_YES$Homoeologues)
syn_C1015_NO_dups_removed <- syn_C1015_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_C1015_dups_removed <- rbind(syn_C1015_YES, syn_C1015_NO_dups_removed)
nlevels(as.factor(syn_diads_triads_unique_C1015_dups_removed$group_num)) #1093

syn_C1704_YES <- subset(syn_diads_triads_unique_C1704, Mutant=="YES")
syn_C1704_NO <- subset(syn_diads_triads_unique_C1704, Mutant=="NO")

syn_C1704_NO_dups_removed <- subset(syn_C1704_NO, !Homoeologues%in%syn_C1704_YES$Homoeologues)
syn_C1704_NO_dups_removed <- syn_C1704_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_C1704_dups_removed <- rbind(syn_C1704_YES, syn_C1704_NO_dups_removed)
nlevels(as.factor(syn_diads_triads_unique_C1704_dups_removed$group_num)) #627


#Add column to specify which comparison is being made
syn_diads_triads_unique_C0604_dups_removed$Comparison <- "C0604vsWT"
syn_diads_triads_unique_C0895_dups_removed$Comparison <- "C0895vsWT"
syn_diads_triads_unique_C1015_dups_removed$Comparison <- "C1015vsWT"
syn_diads_triads_unique_C1704_dups_removed$Comparison <- "C1704vsWT"

#Combine into one dataframe
syn_diads_triads_unique <- rbind(syn_diads_triads_unique_C0604_dups_removed, syn_diads_triads_unique_C0895_dups_removed,
                                 syn_diads_triads_unique_C1015_dups_removed, syn_diads_triads_unique_C1704_dups_removed)
syn_diads_triads_unique$group_num <- factor(syn_diads_triads_unique$group_num)
nlevels(syn_diads_triads_unique$group_num) #2630

#Indicate whether the gene is up or downregulated, or no significant change
syn_diads_triads_unique$Expression <- ifelse(syn_diads_triads_unique$log2FoldChange>0, "UP", "DOWN")
syn_diads_triads_unique$Expression <- ifelse(syn_diads_triads_unique$padj>0.05, "NO CHANGE", syn_diads_triads_unique$Expression)
syn_diads_triads_unique <- subset(syn_diads_triads_unique, is.na(Expression)==FALSE)

#Calculate number of groups where at least one homoeolog is upregulated
syn_upregulated <- subset(syn_diads_triads_unique, Expression=="UP"&Mutant=="NO")
syn_upregulated$group_num <- as.numeric(syn_upregulated$group_num)
syn_upregulated$group_num <- as.factor(syn_upregulated$group_num)
nlevels(syn_upregulated$group_num)

#Make plots
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Output")

summary_syn <- syn_diads_triads_unique %>%
  group_by(Comparison, Mutant, Expression) %>%
  summarise(count=n())

summary_syn$Mutant <- factor(summary_syn$Mutant, levels=c("YES","NO"))

barplot <- ggplot(summary_syn, aes(x=Mutant, y=count, fill=Expression)) +
  geom_bar(position="fill", stat="identity", alpha=0.8, colour="black") +
  theme_bw() +
  scale_fill_manual(values=c("#F0E442","#ABABAB","#009E73"), labels=c("Down", "No Change", "Up")) +
  scale_x_discrete(labels=c("Genes with syn","Homoeologues")) +
  xlab(NULL) +
  ylab("Proportion of Genes") +
  labs(fill="Expression Relative to WT") +
  facet_grid(cols=vars(Comparison)) +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1))
ggsave("up_down_regulated_homoeologues_stacked_padj0.05_syn_dups_removed.pdf", width=6, height=5, units="in", dpi=300)
