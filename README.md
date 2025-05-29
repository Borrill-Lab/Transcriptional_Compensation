# Transcriptional_Compensation

## Data
Raw read files for RNA-seq of the Cadenza and Kronos EMS-mutagenised lines can be found in the ENA under project PRJEB89501.

## Analysis of RNA-seq Data
### Preparing Reference Transcriptome (Kronos only)
- **Filter Chinese Spring reference transcriptome to remove D-subgenome transcripts**: filter_reference_fasta.sh, requires the Chinese Spring reference transcriptome available [here]()
- **Prepare kallisto index**: kallisto_index.sh

### Differential Expression Analysis
- **Pseudoalignment and quantification with kallisto**: kallisto.pl
- **Import quantifications and summarise at the gene level**: import_kallisto_quantifications.R
- **Identify differentially expressed genes in the mutagenised lines**: differential_expression_DESeq2_all_contrasts.R

### Variant Calling
