#Delfi Dorussen
#Aim to combine variant information and differential expression information to identify
#potential homoeologous transcriptional compensation

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/Output")

#Load in syn_variants_diads_triads and subset abundances to only keep triads/diads in this list
syn_diads_triads <- read.csv("../synonymous_variants_diads_triads.csv")
head(syn_diads_triads)

#Pivot syn_triads_diads to have each homoeologue as a separate row
syn_diads_triads <- pivot_longer(syn_diads_triads, 20:22, names_to="Subgenome_Hom", values_to="Homoeologues")
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
syn_diads_triads_unique_dm3 <- subset(syn_diads_triads, dm3=="1/1"&WT=="0/0")
syn_diads_triads_unique_dm4 <- subset(syn_diads_triads, dm4=="1/1"&WT=="0/0")

#Load in contrasts to WT
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/Output/DESeq2_All_Contrasts/")
res_dm3_vs_WT <- read.csv("res_dm3_vs_WT.csv")
res_dm4_vs_WT <- read.csv("res_dm4_vs_WT.csv")

#Merge FC data with variant data
syn_diads_triads_unique_dm3 <- merge(syn_diads_triads_unique_dm3, res_dm3_vs_WT, by.x="Homoeologues", by.y="Gene")
syn_diads_triads_unique_dm4 <- merge(syn_diads_triads_unique_dm4, res_dm4_vs_WT, by.x="Homoeologues", by.y="Gene")

#Getting number of homoeolog groups affected by a syn in each TILLING line
syn_diads_triads_unique_dm3$group_num <- factor(syn_diads_triads_unique_dm3$group_num)
nlevels(syn_diads_triads_unique_dm3$group_num) #194

syn_diads_triads_unique_dm4$group_num <- factor(syn_diads_triads_unique_dm4$group_num)
nlevels(syn_diads_triads_unique_dm4$group_num) #363

#Remove duplicates due to having more than one homoeolog with a variant from the table 
syn_dm3_NO <- subset(syn_diads_triads_unique_dm3, Mutant=="NO")
syn_dm3_YES <- subset(syn_diads_triads_unique_dm3, Mutant=="YES")
syn_dm3_NO_dups_removed <- subset(syn_dm3_NO, !Homoeologues %in% syn_dm3_YES$Homoeologues)
syn_dm3_NO_dups_removed <- syn_dm3_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_dm3_YES <- subset(syn_dm3_YES, group_num %in% syn_dm3_NO_dups_removed$group_num)
syn_dm3_YES_dups_removed <- syn_dm3_YES %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_dm3_dups_removed <- rbind(syn_dm3_YES_dups_removed, syn_dm3_NO_dups_removed)
nlevels(as.factor(syn_diads_triads_unique_dm3_dups_removed$group_num)) #194


syn_dm4_NO <- subset(syn_diads_triads_unique_dm4, Mutant=="NO")
syn_dm4_YES <- subset(syn_diads_triads_unique_dm4, Mutant=="YES")
syn_dm4_NO_dups_removed <- subset(syn_dm4_NO, !Homoeologues %in% syn_dm4_YES$Homoeologues)
syn_dm4_NO_dups_removed <- syn_dm4_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_dm4_YES <- subset(syn_dm4_YES, group_num %in% syn_dm4_NO_dups_removed$group_num)
syn_dm4_YES_dups_removed <- syn_dm4_YES %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_dm4_dups_removed <- rbind(syn_dm4_YES_dups_removed, syn_dm4_NO_dups_removed)
nlevels(as.factor(syn_diads_triads_unique_dm4_dups_removed$group_num)) #363

#Add column to specify which comparison is being made
syn_diads_triads_unique_dm3_dups_removed$Comparison <- "dm3vsWT"
syn_diads_triads_unique_dm4_dups_removed$Comparison <- "dm4vsWT"

#Combine into one dataframe
syn_diads_triads_unique <- rbind(syn_diads_triads_unique_dm3_dups_removed, syn_diads_triads_unique_dm4_dups_removed)
syn_diads_triads_unique$group_num <- factor(syn_diads_triads_unique$group_num)
nlevels(syn_diads_triads_unique$group_num) #536

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
