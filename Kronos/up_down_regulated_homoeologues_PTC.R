#Delfi Dorussen
#Aim to combine variant information and differential expression information to identify
#potential homoeologous transcriptional compensation

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/Output/")

#Load in PTC_variants_diads_triads and subset abundances to only keep triads/diads in this list
PTC_diads_triads <- read.csv("../PTC_variants_diads_triads.csv")
head(PTC_diads_triads)
PTC_diads_triads <- select(PTC_diads_triads, X:B)

#Pivot PTC_triads_diads to have each homoeologue as a separate row
PTC_diads_triads <- pivot_longer(PTC_diads_triads, 24:25, names_to="Subgenome", values_to="Homoeologues")
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

#Subset variants to keep PTCs that are in a TILLING line, but not KWT
PTC_diads_triads_unique_K2619 <- subset(PTC_diads_triads, K2619=="1/1"&KWT=="0/0")
PTC_diads_triads_unique_K2864 <- subset(PTC_diads_triads, K2864=="1/1"&KWT=="0/0")
PTC_diads_triads_unique_K3239 <- subset(PTC_diads_triads, K3239=="1/1"&KWT=="0/0")
PTC_diads_triads_unique_K427 <- subset(PTC_diads_triads, K427=="1/1"&KWT=="0/0")
PTC_diads_triads_unique_K4533 <- subset(PTC_diads_triads, K4533=="1/1"&KWT=="0/0")
PTC_diads_triads_unique_K774 <- subset(PTC_diads_triads, K774=="1/1"&KWT=="0/0")

#Load in contrasts to WT
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/Output/DESeq2_All_Contrasts/")
res_K2619_vs_CWT <- read.csv("res_K2619_vs_WT.csv")
res_K2864_vs_CWT <- read.csv("res_K2864_vs_WT.csv")
res_K3239_vs_CWT <- read.csv("res_K3239_vs_WT.csv")
res_K427_vs_CWT <- read.csv("res_K427_vs_WT.csv")
res_K4533_vs_CWT <- read.csv("res_K4533_vs_WT.csv")
res_K774_vs_CWT <- read.csv("res_K774_vs_WT.csv")

#Merge FC data with variant data
PTC_diads_triads_unique_K2619 <- merge(PTC_diads_triads_unique_K2619, res_K2619_vs_CWT, by.x="Homoeologues", by.y="Gene")
PTC_diads_triads_unique_K2864 <- merge(PTC_diads_triads_unique_K2864, res_K2864_vs_CWT, by.x="Homoeologues", by.y="Gene")
PTC_diads_triads_unique_K3239 <- merge(PTC_diads_triads_unique_K3239, res_K3239_vs_CWT, by.x="Homoeologues", by.y="Gene")
PTC_diads_triads_unique_K427 <- merge(PTC_diads_triads_unique_K427, res_K427_vs_CWT, by.x="Homoeologues", by.y="Gene")
PTC_diads_triads_unique_K4533 <- merge(PTC_diads_triads_unique_K4533, res_K4533_vs_CWT, by.x="Homoeologues", by.y="Gene")
PTC_diads_triads_unique_K774 <- merge(PTC_diads_triads_unique_K774, res_K774_vs_CWT, by.x="Homoeologues", by.y="Gene")

#Getting number of homoeolog groups affected by a PTC in each TILLING line
PTC_diads_triads_unique_K2619$group_num <- factor(PTC_diads_triads_unique_K2619$group_num)
nlevels(PTC_diads_triads_unique_K2619$group_num) #17

PTC_diads_triads_unique_K2864$group_num <- factor(PTC_diads_triads_unique_K2864$group_num)
nlevels(PTC_diads_triads_unique_K2864$group_num) #37

PTC_diads_triads_unique_K3239$group_num <- factor(PTC_diads_triads_unique_K3239$group_num)
nlevels(PTC_diads_triads_unique_K3239$group_num) #9

PTC_diads_triads_unique_K427$group_num <- factor(PTC_diads_triads_unique_K427$group_num)
nlevels(PTC_diads_triads_unique_K427$group_num) #12

PTC_diads_triads_unique_K4533$group_num <- factor(PTC_diads_triads_unique_K4533$group_num)
nlevels(PTC_diads_triads_unique_K4533$group_num) #6

PTC_diads_triads_unique_K774$group_num <- factor(PTC_diads_triads_unique_K774$group_num)
nlevels(PTC_diads_triads_unique_K774$group_num) #31

#Add column to specify which comparison is being made
PTC_diads_triads_unique_K2619$Comparison <- "K2619vsWT"
PTC_diads_triads_unique_K2864$Comparison <- "K2864vsWT"
PTC_diads_triads_unique_K3239$Comparison <- "K3239vsWT"
PTC_diads_triads_unique_K427$Comparison <- "K427vsWT"
PTC_diads_triads_unique_K4533$Comparison <- "K4533vsWT"
PTC_diads_triads_unique_K774$Comparison <- "K774vsWT"

#Combine into one dataframe
PTC_diads_triads_unique <- rbind(PTC_diads_triads_unique_K2619, PTC_diads_triads_unique_K2864,
                                 PTC_diads_triads_unique_K3239, PTC_diads_triads_unique_K427,
                                 PTC_diads_triads_unique_K4533, PTC_diads_triads_unique_K774)
PTC_diads_triads_unique$group_num <- factor(PTC_diads_triads_unique$group_num)
nlevels(PTC_diads_triads_unique$group_num) #101 - so a couple of homoeolog groups must be shared amongst the TILLING lines (expect 112)

#Indicate whether the gene is up or downregulated, or no significant change
PTC_diads_triads_unique$Expression <- ifelse(PTC_diads_triads_unique$log2FoldChange>0, "UP", "DOWN")
PTC_diads_triads_unique$Expression <- ifelse(PTC_diads_triads_unique$padj>0.05, "NO CHANGE", PTC_diads_triads_unique$Expression)
PTC_diads_triads_unique <- subset(PTC_diads_triads_unique, is.na(Expression)==FALSE)

#Calculate number of groups where at least one homoeolog is upregulated
PTC_upregulated <- subset(PTC_diads_triads_unique, Expression=="UP"&Mutant=="NO")
PTC_upregulated$group_num <- as.numeric(PTC_upregulated$group_num)
PTC_upregulated$group_num <- as.factor(PTC_upregulated$group_num)
nlevels(PTC_upregulated$group_num) #3

#Make plots
setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/kallisto_noD/Output")

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
ggsave("up_down_regulated_homoeologues_stacked_padj0.05.pdf", width=8.5, height=5, units="in", dpi=300)
