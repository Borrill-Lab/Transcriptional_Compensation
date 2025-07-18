# Transcriptional_Compensation

## Data
Raw read files for RNA-seq of the Cadenza and Kronos EMS-mutagenised lines can be found in the ENA under project PRJEB89501.

The data from Xiong et al. (2020) is described in this paper: [Transcriptome sequencing reveals hotspot mutation regions and dwarfing mechanisms in wheat mutants induced by γ-ray irradiation and EMS](https://doi.org/10.1093/jrr/rrz075)

## Analysis of RNA-seq Data
The same pipeline was used for the Cadenza lines (scripts [here](https://github.com/Borrill-Lab/Transcriptional_Compensation/tree/main/Cadenza)), Kronos lines (scripts [here](https://github.com/Borrill-Lab/Transcriptional_Compensation/tree/main/Kronos)), and the EMS-mutagenised lines from Xiong et al. (2020) (scripts [here](https://github.com/Borrill-Lab/Transcriptional_Compensation/tree/main/Xiong_EMS)). 

### Preparing Reference Transcriptome (Kronos only)
- **Filter Chinese Spring reference transcriptome to remove D-subgenome transcripts**: filter_reference_fasta.sh, requires the Chinese Spring reference transcriptome available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/)
- **Prepare kallisto index**: kallisto_index.sh

### Differential Expression Analysis
- **Pseudoalignment and quantification with kallisto**: kallisto.pl (requires the Chinese Spring reference transcriptome available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/)
- **Import quantifications and summarise at the gene level**: import_kallisto_quantifications.R
- **Identify differentially expressed genes in the mutagenised lines**: differential_expression_DESeq2_all_contrasts.R

### Variant Calling
- **Trim reads**: trim.pl
- **Map reads to Chinese Spring reference genome**: hisat.pl, requires the Chinese Spring reference genome available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/)
- **Call variants**: call_variants_by_chrom.pl
- **Predict variant effects**: vep.sh
- **Filter variants to keep only PTC variants or synonymous variants**: vep_filter_PTC.sh *or* vep_filter_synonymous.sh
- **Create lists of genes in diads or triads**: genes_in_diads_triads.R
- **Further filtering to keep variants that are homozygous mutant in at least one of the EMS-mutagenised lines, affect canonical transcripts, and are in a diad/triad**: vep_output_filtering_PTC.R *or* vep_output_filtering_synonymous.R

### Combine Differential Expression and Variant Calling
- **Identify potential homoeologous transcriptional compensation associated with PTC or synonymous mutations and plot**: up_down_regulated_homoeologues_PTC.R *or* up_down_regulated_homoeologues_synonymous.R
