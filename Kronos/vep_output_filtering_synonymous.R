#Delfi Dorussen
#Aim to further filter the VEP output to keep only synonymous variants that are homozygous mutant
#in at least one of the TILLING lines, affect canonical transcripts, and are in a diad/triad


library(tidyr)
library(dplyr)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Kronos_TILLING_lines/")

syn_variants <- read.table("variant_effect_output_synonymous.txt")
colnames(syn_variants) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","K2619","K2864","K3239","K427","K4533","K774","KWT")
syn_variants <- separate_wider_delim(syn_variants, INFO, delim="CSQ=", names=c("INFO","VARIANTS"))
syn_variants <- separate_wider_delim(syn_variants, VARIANTS, delim=",", too_few="align_start", names=c("Variant1","Variant2","Variant3","Variant4","Variant5",
                                                                            "Variant6","Variant7","Variant8","Variant9","Variant10","Variant11","Variant12",
                                                                            "Variant13","variant14","Variant15","Variant16","Variant17","Variant18","Variant19",
                                                                            "Variant20","Variant21","Variant22","Variant23","Variant24","Variant25","Variant26"))
syn_variants <- pivot_longer(syn_variants, Variant1:Variant26, values_to="Variant")
syn_variants <- subset(syn_variants, is.na(Variant)==FALSE)

syn_variants <- subset(syn_variants, grepl("synonymous_variant", Variant)==TRUE)

#Only need to keep the first 3 characters in these columns (0/0, 0/1, 1/1 etc.)
syn_variants$K2619 <- substr(syn_variants$K2619, 0, 3)
syn_variants$K2864 <- substr(syn_variants$K2864, 0, 3)
syn_variants$K3239 <- substr(syn_variants$K3239, 0, 3)
syn_variants$K427 <- substr(syn_variants$K427, 0, 3)
syn_variants$K4533 <- substr(syn_variants$K4533, 0, 3)
syn_variants$K774 <- substr(syn_variants$K774, 0, 3)
syn_variants$KWT <- substr(syn_variants$KWT, 0, 3)

#Subset such that at least one of the lines has the non-PTC allele, and at least one has the PTC allele (not het)
syn_variants_1 <- subset(syn_variants, K2619=="0/0"|K2864=="0/0"|K3239=="0/0"|K427=="0/0"|K4533=="0/0"|K774=="0/0"|KWT=="0/0")
syn_variants_1 <- subset(syn_variants_1, K2619=="1/1"|K2864=="1/1"|K3239=="1/1"|K427=="1/1"|K4533=="1/1"|K774=="1/1"|KWT=="1/1")
#This leaves 11,739 variants

#Filter to keep only variants in canonical transcripts
syn_variants_1_canonical <- subset(syn_variants_1, grepl("YES",Variant))
#8116 variants are kept

#Read in list of genes in triads
genes_in_triads <- scan("genes_in_triads.txt", what="character")

#Loop to add column specifying whether the variant is in a gene that is in a triad
syn_variants_1_canonical$TRIAD <- NA
syn_variants_1_canonical$GENE <- NA
i <- 399
for(i in 1:40786){
  #print(i)
  triad_gene <- genes_in_triads[i]
  syn_variants_1_canonical$TRIAD[grepl(triad_gene, syn_variants_1_canonical$Variant)==TRUE] <- "YES"
  syn_variants_1_canonical$GENE[grepl(triad_gene, syn_variants_1_canonical$Variant)==TRUE] <- triad_gene
}

#Read in list of genes in diads
genes_in_diads <- scan("genes_in_diads.txt", what="character")

#Loop to add column specifying whether the variant is in a gene that is in a diad
syn_variants_1_canonical$DIAD <- NA
i <- 399
for(i in 1:5290){
  #print(i)
  diad_gene <- genes_in_diads[i]
  syn_variants_1_canonical$DIAD[grepl(diad_gene, syn_variants_1_canonical$Variant)==TRUE] <- "YES"
  syn_variants_1_canonical$GENE[grepl(diad_gene, syn_variants_1_canonical$Variant)==TRUE] <- diad_gene
}

#New dataframe only with variants in diads or triads
syn_variants_diads_triads <- subset(syn_variants_1_canonical, TRIAD=="YES"|DIAD=="YES")
#Keeps 4453 variants

#Add the homoeolog group to each affected gene
homoeologues <- read.csv("v1.1_genes_TF_homoeolog_info.csv")
homoeologues_groups <- homoeologues[,c(1,3)]
syn_variants_diads_triads <- merge(syn_variants_diads_triads, homoeologues_groups, by.x="GENE", by.y="v1.1_ID")

#Add columns with the homoeologues of the affected genes
homoeologues_groups <- subset(homoeologues_groups, is.na(group_num)==FALSE)
homoeologues_groups$subgenome <- NA
for(i in 1:79871){
  #print(i)
  gene <- homoeologues_groups[i,]
  if(grepl("A0", gene$v1.1_ID)==TRUE){
    homoeologues_groups[i,3] <- "A"
  } else if(grepl("B0", gene$v1.1_ID)==TRUE){
    homoeologues_groups[i,3] <- "B"
  }
}

homoeologues_1 <- homoeologues_groups %>% pivot_wider(values_from=v1.1_ID, names_from=subgenome)
syn_variants_diads_triads <- merge(syn_variants_diads_triads, homoeologues_1, by.x="group_num", by.y="group_num")

write.csv(syn_variants_diads_triads, file="syn_variants_diads_triads.csv")
