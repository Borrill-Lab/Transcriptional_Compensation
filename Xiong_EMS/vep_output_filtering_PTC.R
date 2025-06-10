#Delfi Dorussen
#Aim to further filter the VEP output to keep only PTC variants that are homozygous mutant
#in at least one of the TILLING lines, affect canonical transcripts, and are in a diad/triad

library(tidyr)
library(dplyr)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/Xiong_data_analysis/")

PTC_variants <- read.table("variant_effect_output_PTC.txt")
colnames(PTC_variants) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","dm3","dm4","WT")
PTC_variants <- separate_wider_delim(PTC_variants, INFO, delim="CSQ=", names=c("INFO","VARIANTS"))
PTC_variants <- separate_wider_delim(PTC_variants, VARIANTS, delim=",", too_few="align_start", names=c("Variant1","Variant2","Variant3","Variant4","Variant5",
                                                                            "Variant6","Variant7","Variant8","Variant9","Variant10"))
PTC_variants <- pivot_longer(PTC_variants, Variant1:Variant10, values_to="Variant")
PTC_variants <- subset(PTC_variants, is.na(Variant)==FALSE)

PTC_variants <- subset(PTC_variants, grepl("stop_gained", Variant)==TRUE)

#Only need to keep the first 3 characters in these columns (0/0, 0/1, 1/1 etc.)
PTC_variants$dm3 <- substr(PTC_variants$dm3, 0, 3)
PTC_variants$dm4 <- substr(PTC_variants$dm4, 0, 3)
PTC_variants$WT <- substr(PTC_variants$WT, 0, 3)

#Subset such that at least one of the lines has the non-PTC allele, and at least one has the PTC allele (not het)
PTC_variants_1 <- subset(PTC_variants, dm3=="0/0"|dm4=="0/0"|WT=="0/0")
PTC_variants_1 <- subset(PTC_variants_1, dm3=="1/1"|dm4=="1/1"|WT=="1/1")
#This leaves 63 variants

#Filter to keep only variants in canonical transcripts
PTC_variants_1_canonical <- subset(PTC_variants_1, grepl("YES",Variant))
#42 variants are kept

#Read in list of genes in triads
genes_in_triads <- scan("genes_in_triads.txt", what="character")

#Loop to add column specifying whether the variant is in a gene that is in a triad
PTC_variants_1_canonical$TRIAD <- NA
PTC_variants_1_canonical$GENE <- NA
i <- 399
for(i in 1:61179){
  print(i)
  triad_gene <- genes_in_triads[i]
  PTC_variants_1_canonical$TRIAD[grepl(triad_gene, PTC_variants_1_canonical$Variant)==TRUE] <- "YES"
  PTC_variants_1_canonical$GENE[grepl(triad_gene, PTC_variants_1_canonical$Variant)==TRUE] <- triad_gene
}

#Read in list of genes in diads
genes_in_diads <- scan("genes_in_diads.txt", what="character")

#Loop to add column specifying whether the variant is in a gene that is in a diad
PTC_variants_1_canonical$DIAD <- NA
i <- 399
for(i in 1:18692){
  print(i)
  diad_gene <- genes_in_diads[i]
  PTC_variants_1_canonical$DIAD[grepl(diad_gene, PTC_variants_1_canonical$Variant)==TRUE] <- "YES"
  PTC_variants_1_canonical$GENE[grepl(diad_gene, PTC_variants_1_canonical$Variant)==TRUE] <- diad_gene
}

#New dataframe only with variants in diads or triads
PTC_variants_diads_triads <- subset(PTC_variants_1_canonical, TRIAD=="YES"|DIAD=="YES")
#Keeps 37 variants

#Add the homoeolog group to each affected gene
homoeologues <- read.csv("v1.1_genes_TF_homoeolog_info.csv")
homoeologues_groups <- homoeologues[,c(1,3)]
PTC_variants_diads_triads <- merge(PTC_variants_diads_triads, homoeologues_groups, by.x="GENE", by.y="v1.1_ID")
PTC_variants_diads_triads <- distinct(PTC_variants_diads_triads)

#Add columns with the homoeologues of the affected genes
homoeologues_groups <- subset(homoeologues_groups, is.na(group_num)==FALSE)
homoeologues_groups$subgenome <- NA
for(i in 1:79871){
  print(i)
  gene <- homoeologues_groups[i,]
  if(grepl("A0", gene$v1.1_ID)==TRUE){
    homoeologues_groups[i,3] <- "A"
  } else if(grepl("B0", gene$v1.1_ID)==TRUE){
    homoeologues_groups[i,3] <- "B"
  } else if(grepl("D0", gene$v1.1_ID)==TRUE){
    homoeologues_groups[i,3] <- "D"
  }
}

homoeologues_1 <- homoeologues_groups %>% pivot_wider(values_from=v1.1_ID, names_from=subgenome)
PTC_variants_diads_triads <- merge(PTC_variants_diads_triads, homoeologues_1, by.x="group_num", by.y="group_num")

write.csv(PTC_variants_diads_triads, file="PTC_variants_diads_triads_2.csv")
