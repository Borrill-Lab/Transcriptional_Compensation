#Delfi Dorussen
#Aim to create lists of genes that are either in diads or triads

setwd("U:/Year1/Objective 3/Transcriptional Adaptation Project/")
mydir <- getwd()

###Homoeologue groups
url <- "https://raw.githubusercontent.com/Borrill-Lab/TF_Triads/main/data_files/v1.1_genes_TF_homoeolog_info.csv"
destfile <- file.path(mydir, "v1.1_genes_TF_homoeolog_info.csv")
download.file(url, destfile)
homoeologues <- read.csv("v1.1_genes_TF_homoeolog_info.csv")
head(homoeologues)
tail(homoeologues)

homoeologues_triads <- subset(homoeologues, type=="1:1:1") #Only keep genes in 1:1:1 triads
homoeologues_triads_genes <- homoeologues_triads$v1.1_ID
write.table(homoeologues_triads_genes, "genes_in_triads.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

homoeologues_diads <- subset(homoeologues, type=="1:1:0" | type=="1:0:1" | type=="0:1:1")
homoeologues_diads_genes <- homoeologues_diads$v1.1_ID
write.table(homoeologues_diads_genes, "genes_in_diads.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
