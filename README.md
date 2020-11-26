# Data and codes of the paper "Whole Genome Duplication Promotes Immediate Differentiation for Speciation of Budding Yeast"
## DEG analysis
### Trans_id
Genome annotation from NCBI (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/885/GCF_001417885.1_Kmar_1.0/GCF_001417885.1_Kmar_1.0_genomic.gff.gz) was used to annotate gene symbols.  
The corresponding R script is 1_Trans_id.R.
### Grouping and normalizing
We grouped 21 sequencing samples into 7 groups according to their genetic characteristics.Subsequently, we normalized each group of data.  
The corresponding R script is 2_Grouping_normalizing.R.
### DEG
Then we identified the differentially expressed genes (DEGs) for all group pairs utilizing the DEseq2 package.  
The corresponding R script is DEseq_vst/3_DEG_analysis.R.

## GO enrichment analysis
Genes identified using DEseq2 were analyzed for functional enrichment using the Gene Ontology (GO) database, and the identification of enriched pathways in K.marxianus gene sets was performed with the Bioconductor R package clusterProfiler (Yu et al., 2012). The Bioconductor R package topGO (Alexa and Rahnenfuhrer, 2016) was used to determine GO terms significantly enriched in different gene sets.   
The corresponding R scripts are DEseq_vst/DEseq_DEG.csv/4_Annotation.R and DEseq_vst/DEseq_DEG.csv/5_KM_GO.R.
## Comparison with S. cerevisiae
We compared the RNA sequencing data of the K. marxianus strains with published data for 39 samples of S. cerevisiae (GSE129483). To avoid ambiguous comparisons, only 1345 identified unique genes shared by both K. marxianus and S. cerevisiae were included in our study after quality control and filtering. We employed principle component analysis (PCA) on a Euclidean distance matrix of paired samples.  
The corresponding Data,R scripts and figures are all in the folder k_f_sce
## Analysis of alternations of biological processes
We identified the genes of K. marxianus with significant transcriptional changes due to WGD.After that, we performed an enrichment analysis of alternative genes for gene ontology. The following 6 biological processes are highlighted：DNA replication (GO: 0006260; 142 genes), transmembrane transport (GO: 0022857; 230 genes), cytoskeleton organization (GO: 0007010; 218 genes), oxoacid metabolism (GO: 0043436; 338 genes), and mitotic cell cycle (GO: 0000278; 317 genes). In each biological process，we compared the value of |log2-fold change (FC)| between these three groups：the haploid to MAT a/a and α/α diploid strains，and the diploid MAT a /α to tetploid MAT aa/αα strains. Distributions of the transcriptional alternation between these groups were visualized through the violin plots using the R package ggplot2.  
The corresponding R scripts is DEseq_vst/DEseq_DEG.csv/6_violin_plot.R, and the figures are in the folder DEseq_vst/DEseq_DEG.csv/GO_violin.
