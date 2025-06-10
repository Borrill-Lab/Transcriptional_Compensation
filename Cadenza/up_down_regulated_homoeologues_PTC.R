#Delfi Dorussen
#Aim to combine variant information and differential expression information to identify
#potential homoeologous transcriptional compensation

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Output")

#Load in PTC_variants_diads_triads and subset abundances to only keep triads/diads in this list
PTC_diads_triads <- read.csv("../PTC_variants_diads_triads_2.csv")
head(PTC_diads_triads)

#Pivot PTC_triads_diads to have each homoeologue as a separate row
PTC_diads_triads <- pivot_longer(PTC_diads_triads, 22:24, names_to="Subgenome", values_to="Homoeologues")
PTC_diads_triads <- subset(PTC_diads_triads, is.na(Homoeologues)==FALSE)

#Specify whether the homoeologue has the PTC variant
PTC_diads_triads$Mutant <- NA
for(i in 1:nrow(PTC_diads_triads)){
  print(i)
  if(PTC_diads_triads$GENE[i]==PTC_diads_triads$Homoeologues[i]){
    PTC_diads_triads$Mutant[i] <- "YES"
  } else{
    PTC_diads_triads$Mutant[i] <- "NO"
  }
}

#Subset variants to keep PTCs that are only in one line
PTC_diads_triads_unique_C0604 <- subset(PTC_diads_triads, C0604=="1/1"&CWT=="0/0")
PTC_diads_triads_unique_C0895 <- subset(PTC_diads_triads, C0895=="1/1"&CWT=="0/0")
PTC_diads_triads_unique_C1015 <- subset(PTC_diads_triads, C1015=="1/1"&CWT=="0/0")
PTC_diads_triads_unique_C1704 <- subset(PTC_diads_triads, C1704=="1/1"&CWT=="0/0")

#Load in contrasts to WT
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Output/DESeq2_All_Contrasts/")
res_C0604_vs_CWT <- read.csv("res_C0604_vs_WT.csv")
res_C0895_vs_CWT <- read.csv("res_C0895_vs_WT.csv")
res_C1015_vs_CWT <- read.csv("res_C1015_vs_WT.csv")
res_C1704_vs_CWT <- read.csv("res_C1704_vs_WT.csv")

#Merge FC data with variant data
PTC_diads_triads_unique_C0604 <- merge(PTC_diads_triads_unique_C0604, res_C0604_vs_CWT, by.x="Homoeologues", by.y="Gene")
PTC_diads_triads_unique_C0895 <- merge(PTC_diads_triads_unique_C0895, res_C0895_vs_CWT, by.x="Homoeologues", by.y="Gene")
PTC_diads_triads_unique_C1015 <- merge(PTC_diads_triads_unique_C1015, res_C1015_vs_CWT, by.x="Homoeologues", by.y="Gene")
PTC_diads_triads_unique_C1704 <- merge(PTC_diads_triads_unique_C1704, res_C1704_vs_CWT, by.x="Homoeologues", by.y="Gene")

#Getting number of homoeolog groups affected by a PTC in each TILLING line
PTC_diads_triads_unique_C0604$group_num <- factor(PTC_diads_triads_unique_C0604$group_num)
nlevels(PTC_diads_triads_unique_C0604$group_num) #34

PTC_diads_triads_unique_C0895$group_num <- factor(PTC_diads_triads_unique_C0895$group_num)
nlevels(PTC_diads_triads_unique_C0895$group_num) #48 - 2 groups (723 and 15749) have multiple PTCs in the same homoeolog so need to remove dups

PTC_C0895_NO <- subset(PTC_diads_triads_unique_C0895, Mutant=="NO")
PTC_C0895_YES <- subset(PTC_diads_triads_unique_C0895, Mutant=="YES")
PTC_C0895_NO_dups_removed <- subset(PTC_C0895_NO, !Homoeologues %in% PTC_C0895_YES$Homoeologues)
PTC_C0895_NO_dups_removed <- PTC_C0895_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

PTC_C0895_YES <- subset(PTC_C0895_YES, group_num %in% PTC_C0895_NO_dups_removed$group_num)
PTC_C0895_YES_dups_removed <- PTC_C0895_YES %>%
  distinct(Homoeologues, .keep_all=TRUE)

PTC_diads_triads_unique_C0895_dups_removed <- rbind(PTC_C0895_YES_dups_removed, PTC_C0895_NO_dups_removed)

PTC_diads_triads_unique_C1015$group_num <- factor(PTC_diads_triads_unique_C1015$group_num)
nlevels(PTC_diads_triads_unique_C1015$group_num) #38 - here there is one group that has 2 homoeologs with PTC (20359) so need to remove dups

PTC_C1015_YES <- subset(PTC_diads_triads_unique_C1015, Mutant=="YES")
PTC_C1015_NO <- subset(PTC_diads_triads_unique_C1015, Mutant=="NO")

PTC_C1015_NO_dups_removed <- subset(PTC_C1015_NO, !Homoeologues%in%PTC_C1015_YES$Homoeologues)
PTC_C1015_NO_dups_removed <- PTC_C1015_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

PTC_diads_triads_unique_C1015_dups_removed <- rbind(PTC_C1015_YES, PTC_C1015_NO_dups_removed)
nlevels(as.factor(PTC_diads_triads_unique_C1015_dups_removed$group_num)) #38

PTC_diads_triads_unique_C1704$group_num <- factor(PTC_diads_triads_unique_C1704$group_num)
nlevels(PTC_diads_triads_unique_C1704$group_num) #40

#Add column to specify which comparison is being made
PTC_diads_triads_unique_C0604$Comparison <- "C0604vsWT"
PTC_diads_triads_unique_C0895_dups_removed$Comparison <- "C0895vsWT"
PTC_diads_triads_unique_C1015_dups_removed$Comparison <- "C1015vsWT"
PTC_diads_triads_unique_C1704$Comparison <- "C1704vsWT"

#Combine into one dataframe
PTC_diads_triads_unique <- rbind(PTC_diads_triads_unique_C0604, PTC_diads_triads_unique_C0895_dups_removed,
                                 PTC_diads_triads_unique_C1015_dups_removed, PTC_diads_triads_unique_C1704)
PTC_diads_triads_unique$group_num <- as.factor(PTC_diads_triads_unique$group_num)
nlevels(PTC_diads_triads_unique$group_num) #158 - so 2 homoeolog groups must be shared amongst the TILLING lines

#Indicate whether the gene is up or downregulated, or no significant change
PTC_diads_triads_unique$Expression <- ifelse(PTC_diads_triads_unique$log2FoldChange>0, "UP", "DOWN")
PTC_diads_triads_unique$Expression <- ifelse(PTC_diads_triads_unique$padj>0.05, "NO CHANGE", PTC_diads_triads_unique$Expression)
PTC_diads_triads_unique <- subset(PTC_diads_triads_unique, is.na(Expression)==FALSE)


#Make plots
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Output")

summary_PTC <- PTC_diads_triads_unique %>%
  group_by(Comparison, Mutant, Expression) %>%
  summarise(count=n())

summary_PTC$Mutant <- factor(summary_PTC$Mutant, levels=c("YES","NO"))

barplot <- ggplot(summary_PTC, aes(x=Mutant, y=count, fill=Expression)) +
  geom_bar(position="fill", stat="identity", alpha=0.8, colour="black") +
  theme_bw() +
  scale_fill_manual(values=c("#F0E442","#ABABAB","#009E73"), labels=c("Down", "No Change", "Up")) +
  scale_x_discrete(labels=c("Genes with PTC","Homoeologues")) +
  xlab(NULL) +
  ylab("Proportion of Genes") +
  labs(fill="Expression Relative to WT") +
  facet_grid(cols=vars(Comparison)) +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1))
ggsave("up_down_regulated_homoeologues_stacked_padj0.05_dups_removed.pdf", width=6, height=5, units="in", dpi=300)
