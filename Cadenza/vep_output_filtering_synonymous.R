#Delfi Dorussen
#Aim to further filter the VEP output to keep only synonymous variants that are homozygous mutant
#in at least one of the TILLING lines, affect canonical transcripts, and are in a diad/triad

library(tidyr)
library(dplyr)

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project")

synonymous_variants <- read.table("variant_effect_output_synonymous.txt")
colnames(synonymous_variants) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","C0604","C0895","C1015","C1704","CWT")
synonymous_variants <- separate_wider_delim(synonymous_variants, INFO, delim="CSQ=", names=c("INFO","VARIANTS"))
synonymous_variants <- separate_wider_delim(synonymous_variants, VARIANTS, delim=",", too_few="align_start", names=c("Variant1","Variant2","Variant3","Variant4","Variant5",
                                                                            "Variant6","Variant7","Variant8","Variant9","Variant10","Variant11","Variant12","Variant13",
                                                                            "Variant14","Variant15","Variant16","Variant17","Variant18","Variant19","Variant20","Variant21"))
synonymous_variants <- pivot_longer(synonymous_variants, Variant1:Variant21, values_to="Variant")
synonymous_variants <- subset(synonymous_variants, is.na(Variant)==FALSE)

synonymous_variants <- subset(synonymous_variants, grepl("synonymous_variant", Variant)==TRUE)

#Only need to keep the first 3 characters in these columns (0/0, 0/1, 1/1 etc.)
synonymous_variants$C0604 <- substr(synonymous_variants$C0604, 0, 3)
synonymous_variants$C0895 <- substr(synonymous_variants$C0895, 0, 3)
synonymous_variants$C1015 <- substr(synonymous_variants$C1015, 0, 3)
synonymous_variants$C1704 <- substr(synonymous_variants$C1704, 0, 3)
synonymous_variants$CWT <- substr(synonymous_variants$CWT, 0, 3)

#Subset such that at least one of the lines has the non-PTC allele, and at least one has the PTC allele (not het)
synonymous_variants_1 <- subset(synonymous_variants, C0604=="0/0"|C0895=="0/0"|C1015=="0/0"|C1704=="0/0"|CWT=="0/0")
synonymous_variants_1 <- subset(synonymous_variants_1, C0604=="1/1"|C0895=="1/1"|C1015=="1/1"|C1704=="1/1"|CWT=="1/1")
#This leaves 13,861 variants

#Filter to keep only variants in canonical transcripts
synonymous_variants_1_canonical <- subset(synonymous_variants_1, grepl("YES",Variant))
#10,323 variants are kept

#Read in list of genes in triads
genes_in_triads <- scan("genes_in_triads.txt", what="character")

#Loop to add column specifying whether the variant is in a gene that is in a triad
synonymous_variants_1_canonical$TRIAD <- NA
synonymous_variants_1_canonical$GENE <- NA
i <- 399
for(i in 1:61179){
  print(i)
  triad_gene <- genes_in_triads[i]
  synonymous_variants_1_canonical$TRIAD[grepl(triad_gene, synonymous_variants_1_canonical$Variant)==TRUE] <- "YES"
  synonymous_variants_1_canonical$GENE[grepl(triad_gene, synonymous_variants_1_canonical$Variant)==TRUE] <- triad_gene
}

#Read in list of genes in diads
genes_in_diads <- scan("genes_in_diads.txt", what="character")

#Loop to add column specifying whether the variant is in a gene that is in a diad
synonymous_variants_1_canonical$DIAD <- NA
i <- 399
for(i in 1:18692){
  print(i)
  diad_gene <- genes_in_diads[i]
  synonymous_variants_1_canonical$DIAD[grepl(diad_gene, synonymous_variants_1_canonical$Variant)==TRUE] <- "YES"
  synonymous_variants_1_canonical$GENE[grepl(diad_gene, synonymous_variants_1_canonical$Variant)==TRUE] <- diad_gene
}

#New dataframe only with variants in diads or triads
synonymous_variants_diads_triads <- subset(synonymous_variants_1_canonical, TRIAD=="YES"|DIAD=="YES")
#Keeps 8474 variants

#Add the homoeolog group to each affected gene
homoeologues <- read.csv("v1.1_genes_TF_homoeolog_info.csv")
homoeologues_groups <- homoeologues[,c(1,3)]
synonymous_variants_diads_triads <- merge(synonymous_variants_diads_triads, homoeologues_groups, by.x="GENE", by.y="v1.1_ID")

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
synonymous_variants_diads_triads <- merge(synonymous_variants_diads_triads, homoeologues_1, by.x="group_num", by.y="group_num")
synonymous_variants_diads_triads <- synonymous_variants_diads_triads[!duplicated(synonymous_variants_diads_triads[,c("GENE","C0604","C0895","C1015", "C1704", "CWT")]),]
#Keep 5466 variants

write.csv(synonymous_variants_diads_triads, file="synonymous_variants_diads_triads.csv")
