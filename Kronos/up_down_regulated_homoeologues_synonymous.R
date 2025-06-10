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

#Remove duplicates due to having more than one homoeolog with a variant from the table
#And remove any homoeolog groups which don't have one homoeolog as WT
syn_K2619_NO <- subset(syn_diads_triads_unique_K2619, Mutant=="NO")
syn_K2619_YES <- subset(syn_diads_triads_unique_K2619, Mutant=="YES")
syn_K2619_NO_dups_removed <- subset(syn_K2619_NO, !Homoeologues %in% syn_K2619_YES$Homoeologues)
syn_K2619_NO_dups_removed <- syn_K2619_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_K2619_YES <- subset(syn_K2619_YES, group_num %in% syn_K2619_NO_dups_removed$group_num)
syn_K2619_YES_dups_removed <- syn_K2619_YES %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_K2619_dups_removed <- rbind(syn_K2619_YES_dups_removed, syn_K2619_NO_dups_removed)
syn_diads_triads_unique_K2619_dups_removed$group_num <- as.character(syn_diads_triads_unique_K2619_dups_removed$group_num)
nlevels(as.factor(syn_diads_triads_unique_K2619_dups_removed$group_num)) #589

syn_K2864_NO <- subset(syn_diads_triads_unique_K2864, Mutant=="NO")
syn_K2864_YES <- subset(syn_diads_triads_unique_K2864, Mutant=="YES")
syn_K2864_NO_dups_removed <- subset(syn_K2864_NO, !Homoeologues %in% syn_K2864_YES$Homoeologues)
syn_K2864_NO_dups_removed <- syn_K2864_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_K2864_YES <- subset(syn_K2864_YES, group_num %in% syn_K2864_NO_dups_removed$group_num)
syn_K2864_YES_dups_removed <- syn_K2864_YES %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_K2864_dups_removed <- rbind(syn_K2864_YES_dups_removed, syn_K2864_NO_dups_removed)
syn_diads_triads_unique_K2864_dups_removed$group_num <- as.character(syn_diads_triads_unique_K2864_dups_removed$group_num)
nlevels(as.factor(syn_diads_triads_unique_K2864_dups_removed$group_num)) #685

syn_K3239_NO <- subset(syn_diads_triads_unique_K3239, Mutant=="NO")
syn_K3239_YES <- subset(syn_diads_triads_unique_K3239, Mutant=="YES")
syn_K3239_NO_dups_removed <- subset(syn_K3239_NO, !Homoeologues %in% syn_K3239_YES$Homoeologues)
syn_K3239_NO_dups_removed <- syn_K3239_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_K3239_YES <- subset(syn_K3239_YES, group_num %in% syn_K3239_NO_dups_removed$group_num)
syn_K3239_YES_dups_removed <- syn_K3239_YES %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_K3239_dups_removed <- rbind(syn_K3239_YES_dups_removed, syn_K3239_NO_dups_removed)
syn_diads_triads_unique_K3239_dups_removed$group_num <- as.character(syn_diads_triads_unique_K3239_dups_removed$group_num)
nlevels(as.factor(syn_diads_triads_unique_K3239_dups_removed$group_num)) #210

syn_K427_NO <- subset(syn_diads_triads_unique_K427, Mutant=="NO")
syn_K427_YES <- subset(syn_diads_triads_unique_K427, Mutant=="YES")
syn_K427_NO_dups_removed <- subset(syn_K427_NO, !Homoeologues %in% syn_K427_YES$Homoeologues)
syn_K427_NO_dups_removed <- syn_K427_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_K427_YES <- subset(syn_K427_YES, group_num %in% syn_K427_NO_dups_removed$group_num)
syn_K427_YES_dups_removed <- syn_K427_YES %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_K427_dups_removed <- rbind(syn_K427_YES_dups_removed, syn_K427_NO_dups_removed)
syn_diads_triads_unique_K427_dups_removed$group_num <- as.character(syn_diads_triads_unique_K427_dups_removed$group_num)
nlevels(as.factor(syn_diads_triads_unique_K427_dups_removed$group_num)) #238

syn_K4533_NO <- subset(syn_diads_triads_unique_K4533, Mutant=="NO")
syn_K4533_YES <- subset(syn_diads_triads_unique_K4533, Mutant=="YES")
syn_K4533_NO_dups_removed <- subset(syn_K4533_NO, !Homoeologues %in% syn_K4533_YES$Homoeologues)
syn_K4533_NO_dups_removed <- syn_K4533_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_K4533_YES <- subset(syn_K4533_YES, group_num %in% syn_K4533_NO_dups_removed$group_num)
syn_K4533_YES_dups_removed <- syn_K4533_YES %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_K4533_dups_removed <- rbind(syn_K4533_YES_dups_removed, syn_K4533_NO_dups_removed)
syn_diads_triads_unique_K4533_dups_removed$group_num <- as.character(syn_diads_triads_unique_K4533_dups_removed$group_num)
nlevels(as.factor(syn_diads_triads_unique_K4533_dups_removed$group_num)) #288

syn_K774_NO <- subset(syn_diads_triads_unique_K774, Mutant=="NO")
syn_K774_YES <- subset(syn_diads_triads_unique_K774, Mutant=="YES")
syn_K774_NO_dups_removed <- subset(syn_K774_NO, !Homoeologues %in% syn_K774_YES$Homoeologues)
syn_K774_NO_dups_removed <- syn_K774_NO_dups_removed %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_K774_YES <- subset(syn_K774_YES, group_num %in% syn_K774_NO_dups_removed$group_num)
syn_K774_YES_dups_removed <- syn_K774_YES %>%
  distinct(Homoeologues, .keep_all=TRUE)

syn_diads_triads_unique_K774_dups_removed <- rbind(syn_K774_YES_dups_removed, syn_K774_NO_dups_removed)
syn_diads_triads_unique_K774_dups_removed$group_num <- as.character(syn_diads_triads_unique_K774_dups_removed$group_num)
nlevels(as.factor(syn_diads_triads_unique_K774_dups_removed$group_num)) #564

#Add column to specify which comparison is being made
syn_diads_triads_unique_K2619_dups_removed$Comparison <- "K2619vsWT"
syn_diads_triads_unique_K2864_dups_removed$Comparison <- "K2864vsWT"
syn_diads_triads_unique_K3239_dups_removed$Comparison <- "K3239vsWT"
syn_diads_triads_unique_K427_dups_removed$Comparison <- "K427vsWT"
syn_diads_triads_unique_K4533_dups_removed$Comparison <- "K4533vsWT"
syn_diads_triads_unique_K774_dups_removed$Comparison <- "K774vsWT"

#Combine into one dataframe
syn_diads_triads_unique <- rbind(syn_diads_triads_unique_K2619_dups_removed, syn_diads_triads_unique_K2864_dups_removed,
                                 syn_diads_triads_unique_K3239_dups_removed, syn_diads_triads_unique_K427_dups_removed,
                                 syn_diads_triads_unique_K4533_dups_removed, syn_diads_triads_unique_K774_dups_removed)
syn_diads_triads_unique$group_num <- factor(syn_diads_triads_unique$group_num)
nlevels(syn_diads_triads_unique$group_num) #1772

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
ggsave("up_down_regulated_homoeologues_stacked_padj0.05_syn_dups_removed.pdf", width=8.5, height=5, units="in", dpi=300)

#Count number of dyads with upregulation
up_homoeologs <- subset(syn_diads_triads_unique, Mutant=="NO"&Expression=="UP")
up_homoeologs$group_num <- as.character(up_homoeologs$group_num)

up_homoeologs$group_num <- factor(up_homoeologs$group_num)
nlevels(up_homoeologs$group_num) #26
