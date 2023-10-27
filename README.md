# Kreysing_Franze

## Mechanical regulation of electrical maturation in hippocampal neurons

## Abstract

In this project, we cultured rat hippocampal neurons in mechanically different environments and studied electrical maturation. We compared WT neurons with two different Piezo1 knockdown conditions. These two knockdown conditions were generated in two independent CRISPR-Cas9 KD assays. Each assay is based on four different CRISPR-Cas9 guides targeting the Piezo1 gene. RNA sequencing was used to analyse the pathway leading to the stiffness dependent maturation behaviour.

## SampleTable
SampleTable_summary.csv [[CSV](SampleTable_summary.csv)]

## Raw fastqs and processed counts data Availability

All the associated data presented in this paper are available from the corresponding author upon reasonable request.
Raw RNA-sequencing data and processed counts data are accessible through the EMBL-EBI ArrayExpress with accession number E-MTAB-13503. (https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13503) Mapping data summary and some differential analysis results are shown in Supplementary_Table_26102023.xlsx[[xlsx](Figures_Tables/Supplementary_Table_26102023.xlsx)]

## Tool Arguments for Alignments
Tool Name      |   Arguments          |
---------------|----------------------|
TrimGalore | 	–illumina; –gzip; –fastqc; –fastqc_args ‘–nogroup –extract’|
STAR 	|–runThreadN 4; –outSAMtype BAM SortedByCoordinate; –readFilesCommand zcat|
HTSeq |	-a 10; -m union; -s no; -t exon|

## Reference files links

**1) Reference fasta file**

ftp://ftp.ensembl.org/pub/release-101/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz

**2)Gene model reference (GTF file)**

ftp://ftp.ensembl.org/pub/release-101/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.101.gtf.gz

## Code availability

Custom codes for RNA-sequencing analysis used in this paper are available on on Github (https://github.com/CTR-BFX/Kreysing_Franze, DOI:xxxxxxx).

## RNASeq analysis Methods

Hippocampal tissue was collected from E17-18 embryos from four pregnant rats. These four biological replicates are referred to as R1-4. For each replicate, the hippocampal tissue was pooled and neurons were isolated. Neurons were divided into 3 groups for each replicate: H1, H2, and HWT. H1 and H2 refer to Piezo1 KD(Knockdown) cells whereas HWT represents a wildtype control. H1 and H2 were electroporated with a set of 4 independent Piezo1 guides and Cas9 proteins, HWT was electroporated with CRISPR RNA and Cas9 protein but without a guide sequence as described in above. These groups of cells were plated on soft and stiff hydrogels (100 Pa and 10 kPa) in 35 mm Petri dishes as described above. <br>

After 7 days, the hydrogels were scraped off the Petri dishes and all material was processed with the RNeasy Plus Micro Kit (Qiagen #74034). RNA reads from replicate R4 were excluded in our analysis, since they originated from a much bigger litter (21 embryos) compared to the other replicates (13, 15, 15 embryos, respectively) which might have impacted the developmental stage of the embryos. <br>

The RNA samples were sequenced by Cambridge Genomic Services(CGS) via their low input RNA assay. Raw fastq files were submitted to EMBL-EBI ArrayExpress with accession number E-MTAB-13503(https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13503). Single end, 75 bp length RNA sequence quality control were performed using fastqc (version 0.11.9)<sup>1</sup>, TrimGalore (version 0.6.6) <sup>2</sup>, then aligned to the Rattus Norvegicus genome (Rnor_6.0) using STAR(version 2.6.1d)<sup>3</sup>, and gene counts were generated using HTSeq (subread version 2.0.1)<sup>4</sup>. The above analysis was performed by CGS. The summary of the mapping statistics and number of genes identified for each library is given in Supplementary Table 1(STable1). <br>


The differential expression (DE) analysis was mainly performed using R software (version 4.2.2)<sup>5</sup> *edgeR* (version 3.40.2)<sup>6</sup> pipeline. Firstly, *ComBat_seq*(sva version 3.50.0)<sup>7</sup> to remove any batch effects and then filtered the genes out where there are less than 4 samples with normalized counts greater than or equal to 5. Secondly, calculate the normalization factors and correct for GC content and gene length bias for all genes which pass filtering using *cqn*(version 1.48.0)<sup>8</sup>. Finally, We used the **glm-LR** test *edgeR* approach to make pairwise comparisons between groups. We then adjust for multiple testing via the FDR (Benjamini-Hochberg) approach used by *edgeR*. <br>

Fully DEGs normalised counts and paired significant DEGs' lists are in Supplementary Table 2-11 and also be available in GitHub https://github.com/xz289/Kreysing_Franze(DOI:xxx). Individual gene counts plot used *DESeq2* (version 1.38.3)<sup>9</sup> median of ratios normalisation (counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene) method,  and then a log2(normalised counts + 1) transformation was performed.

[1] S Andrews et al. Fastqc: A quality control tool for high throughput sequence data. Reference Source, 2010.

[2] F Krueger. Trim galore. A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation Bisufite- Seq) libraries, 2013.

[3] Alexander Dobin, Carrie A Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R Gingeras. Star: ultrafast universal rna-seq aligner. Bioinformatics, 29(1):15–21, 2013.

[4] Simon Anders, Paul Theodor Pyl, and Wolfgang Huber. Htseq—a python framework to work with high-throughput sequencing data. Bioinformatics, 31(2):166–169, 2015.

[5] R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

[6] Mark D Robinson, Davis J McCarthy, and Gordon K Smyth. edger: a bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1):139–140, 2010.

[7]	Zhang, Y., Parmigiani, G. & Johnson, W. E. ComBat-seq: batch effect adjustment for RNA-seq count data. NAR Genomics Bioinforma. 2, lqaa078 (2020).

[8] Hansen KD, Irizarry RA, Wu Z (2012). “Removing technical variability in RNA-seq data using conditional quantile normalization.” Biostatistics, 13(2), 204–216.

[9] Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.

## Figures and Tables with Legend

1) **paired comparison PCA and Volcano Plot**


Design          |   PCA(Download)    |PCA(Image) |  Volcano(Download) | Volcano(Image)|
--------        |  --------------------|------|----------|------|
H1_10 vs H1_100 | [[PDF](Figures_Tables/CTR_kf284_0004-PCAplot_H1_10vsH1_100_Top500_rm4R.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-PCAplot_H1_10vsH1_100_Top500_rm4R.png" width=300px> |[[PDF](Figures_Tables/CTR_kf284_0004-edgeR_H1_10vsH1_100_rm4R_volcano_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-edgeR_H1_10vsH1_100_rm4R_volcano_plot_new.png" width=300px> |
H2_10 vs H2_100 | [[PDF](Figures_Tables/CTR_kf284_0004-PCAplot_H2_10vsH2_100_Top500_rm4R.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-PCAplot_H2_10vsH2_100_Top500_rm4R.png" width=300px> |[[PDF](Figures_Tables/CTR_kf284_0004-edgeR_H2_10vsH2_100_rm4R_volcano_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-edgeR_H2_10vsH2_100_rm4R_volcano_plot_new.png" width=300px> |
WT_10 vs WT_100 | [[PDF](Figures_Tables/CTR_kf284_0004-PCAplot_WT_10vsWT_100_Top500_rm4R.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-PCAplot_WT_10vsWT_100_Top500_rm4R.png" width=300px> |[[PDF](Figures_Tables/CTR_kf284_0004-edgeR_WT_10vsWT_100_rm4R_volcano_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-edgeR_WT_10vsWT_100_rm4R_volcano_plot_new.png" width=300px> |
H2_10 vs H1_10 | [[PDF](Figures_Tables/CTR_kf284_0004-PCAplot_H2_10vsH1_10_Top500_rm4R.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-PCAplot_H2_10vsH1_10_Top500_rm4R.png" width=300px> |[[PDF](Figures_Tables/CTR_kf284_0004-edgeR_H2_10vsH1_10_rm4R_volcano_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-edgeR_H2_10vsH1_10_rm4R_volcano_plot_new.png" width=300px> |
H2_100 vs H1_100 | [[PDF](Figures_Tables/CTR_kf284_0004-PCAplot_H2_100vsH1_100_Top500_rm4R.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-PCAplot_H2_100vsH1_100_Top500_rm4R.png" width=300px> |[[PDF](Figures_Tables/CTR_kf284_0004-edgeR_H2_100vsH1_100_rm4R_volcano_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-edgeR_H2_100vsH1_100_rm4R_volcano_plot_new.png" width=300px> |
H1_10 vs WT_10 | [[PDF](Figures_Tables/CTR_kf284_0004-PCAplot_H1_10vsWT_10_Top500_rm4R.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-PCAplot_H1_10vsWT_10_Top500_rm4R.png" width=300px> |[[PDF](Figures_Tables/CTR_kf284_0004-edgeR_H1_10vsWT_10_rm4R_volcano_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-edgeR_H1_10vsWT_10_rm4R_volcano_plot_new.png" width=300px> |
H2_10 vs WT_10 | [[PDF](Figures_Tables/CTR_kf284_0004-PCAplot_H2_10vsWT_10_Top500_rm4R.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-PCAplot_H2_10vsWT_10_Top500_rm4R.png" width=300px> |[[PDF](Figures_Tables/CTR_kf284_0004-edgeR_H2_10vsWT_10_rm4R_volcano_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-edgeR_H2_10vsWT_10_rm4R_volcano_plot_new.png" width=300px> |
H1_100 vs WT_100 | [[PDF](Figures_Tables/CTR_kf284_0004-PCAplot_H1_100vsWT_100_Top500_rm4R.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-PCAplot_H1_100vsWT_100_Top500_rm4R.png" width=300px> |[[PDF](Figures_Tables/CTR_kf284_0004-edgeR_H1_100vsWT_100_rm4R_volcano_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-edgeR_H1_100vsWT_100_rm4R_volcano_plot_new.png" width=300px> |
H2_100 vs WT_100 | [[PDF](Figures_Tables/CTR_kf284_0004-PCAplot_H2_100vsWT_100_Top500_rm4R.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-PCAplot_H2_100vsWT_100_Top500_rm4R.png" width=300px> |[[PDF](Figures_Tables/CTR_kf284_0004-edgeR_H2_100vsWT_100_rm4R_volcano_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-edgeR_H2_100vsWT_100_rm4R_volcano_plot_new.png" width=300px> |

2) **Selected individual genes log2(Normalised Counts) plot**

GeneName         |   Download    | Image |
-----------------|---------------|-------|
Ttr|[[PDF](Figures_Tables/CTR_kf284_0004-Ttr_Count_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-Ttr_Count_plot_new.png" width=300px> |
Thbs1|[[PDF](Figures_Tables/CTR_kf284_0004-Thbs1_Count_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-Thbs1_Count_plot_new.png" width=300px> |
Pcsk9|[[PDF](Figures_Tables/CTR_kf284_0004-Pcsk9_Count_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-Pcsk9_Count_plot_new.png" width=300px> |
Optc|[[PDF](Figures_Tables/CTR_kf284_0004-Optc_Count_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-Optc_Count_plot_new.png" width=300px> |
Gfap|[[PDF](Figures_Tables/CTR_kf284_0004-Gfap_Count_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-Gfap_Count_plot_new.png" width=300px> |
Adam3a|[[PDF](Figures_Tables/CTR_kf284_0004-Adam3a_Count_plot_new.pdf)]|<IMG SRC="Figures_Tables/CTR_kf284_0004-Adam3a_Count_plot_new.png" width=300px> |
