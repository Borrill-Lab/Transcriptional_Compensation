#Delfi Dorussen
#Aim to combine variant information and differential expression information to identify
#potential homoeologous transcriptional compensation

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/Output/")

#Load in PTC_variants_diads_triads and subset abundances to only keep triads/diads in this list
PTC_diads_triads <- read.csv("../PTC_variants_diads_triads_2.csv")
head(PTC_diads_triads)

#Pivot PTC_triads_diads to have each homoeologue as a separate row
PTC_diads_triads <- pivot_longer(PTC_diads_triads, 20:22, names_to="Subgenome", values_to="Homoeologues")
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
PTC_diads_triads_unique_dm3 <- subset(PTC_diads_triads, dm3=="1/1"&WT=="0/0")
PTC_diads_triads_unique_dm4 <- subset(PTC_diads_triads, dm4=="1/1"&WT=="0/0")

#Load in contrasts to WT
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/Output/DESeq2_All_Contrasts/")
res_dm3_vs_CWT <- read.csv("res_dm3_vs_WT.csv")
res_dm4_vs_CWT <- read.csv("res_dm4_vs_WT.csv")

#Merge FC data with variant data
PTC_diads_triads_unique_dm3 <- merge(PTC_diads_triads_unique_dm3, res_dm3_vs_CWT, by.x="Homoeologues", by.y="Gene")
PTC_diads_triads_unique_dm4 <- merge(PTC_diads_triads_unique_dm4, res_dm4_vs_CWT, by.x="Homoeologues", by.y="Gene")

#Getting number of homoeolog groups affected by a PTC in each TILLING line
PTC_diads_triads_unique_dm3$group_num <- factor(PTC_diads_triads_unique_dm3$group_num)
nlevels(PTC_diads_triads_unique_dm3$group_num) #7

PTC_diads_triads_unique_dm4$group_num <- factor(PTC_diads_triads_unique_dm4$group_num)
nlevels(PTC_diads_triads_unique_dm4$group_num) #23

#Add column to specify which comparison is being made
PTC_diads_triads_unique_dm3$Comparison <- "dm3vsWT"
PTC_diads_triads_unique_dm4$Comparison <- "dm4vsWT"

#Combine into one dataframe
PTC_diads_triads_unique <- rbind(PTC_diads_triads_unique_dm3, PTC_diads_triads_unique_dm4)
PTC_diads_triads_unique$group_num <- factor(PTC_diads_triads_unique$group_num)
nlevels(PTC_diads_triads_unique$group_num) #30 homoeolog groups 

#Indicate whether the gene is up or downregulated, or no significant change
PTC_diads_triads_unique$Expression <- ifelse(PTC_diads_triads_unique$log2FoldChange>0, "UP", "DOWN")
PTC_diads_triads_unique$Expression <- ifelse(PTC_diads_triads_unique$padj>0.05, "NO CHANGE", PTC_diads_triads_unique$Expression)
PTC_diads_triads_unique <- subset(PTC_diads_triads_unique, is.na(Expression)==FALSE)

#Make plots
summary_PTC <- PTC_diads_triads_unique %>%
  group_by(Comparison, Mutant, Expression) %>%
  summarise(count=n())

summary_PTC$Mutant <- factor(summary_PTC$Mutant, levels=c("YES","NO"))

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/Output/")
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
ggsave("up_down_regulated_homoeologues_stacked_padj0.05.pdf", width=6, height=5, units="in", dpi=300)
