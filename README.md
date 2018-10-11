Scripts used in the Holo-seq data analysing.
Note that all the files should be provided with full path if not in the same directory as the working directory.


###########################################

#RNAseq.sh

#Description: This is a shell script used to map fastq data to genome(mm9 or hg19) and get gene read count for mRNA and total RNA sequencing libraries.

#Usage: 

	RNASeq.sh $1
$1 is the basename of fastq file of the sample used call small RNA.

#Note:

make sure that GTF file, genome index file, CRG-Flow software and gene_exon_loci.txt file are properly prepared, details see the script and FileFormatDemo.txt.

The working directory needs to be the same with the fastq file's.


###########################################

#SmallRNASeq.sh

#Description: This is a shell script used to call miRNA, snoRNA and tsRNA from small RNA sequencing libraries.

#Usage: 

	SmallRNASeq.sh $1
$1 is the basename of fastq file of the sample used call small RNA.

#Note:

make sure that GTF file, genome index file, small RNA index files, miRNA_preMiRNA_loci.txt file and CountMatureMirna.R file are 
properly prepared, details see the script and FileFormatDemo.txt.

bbduk.sh is a sub script of BBmap.

The working directory needs to be the same with the fastq file's.

###########################################

#FindNegPairs.R

#Description: This is an R script used to find the miRNA and mRNA target pairs expressed in significant negtive correlation with miRNA-
mRNA dual sequencing data.

#Usage: 

	R CMD BATCH "--args $ReadCountMatrix_mRNA.csv $ReadCountMatrix_miRNA.csv $miRNA_mRNA_StrongProvedPairs.txt"
$ReadCountMatrix_mRNA.csv and $ReadCountMatrix_miRNA.csv files are read count matrix files of mRNA and miRNA, see FileFormatDemo.txt.

$miRNA_mRNA_StrongProvedPairs.txt file contains the miRNA and target mRNA pairs that have been proved by multiple experiments.


###########################################

#HierarchicalClustering_heatmap.R

#Description: This is an R script used to do hierarchical clustering and draw heatmaps

#Usage: 

	R CMD BATCH "--args $ReadCountMatirx.csv $RNA $cut1 $cut2"

$ReadCountMatirx.csv is a read count matrix file of mRNA or miRNA, see FileFormatDemo.txt.

$RNA can be "miRNA" or "mRNA", depends on your data.

$cut1 and cut2 are parameters used to select genes/miRNAs that have no less than $cut1 mapped reads in no less than $cut2 samples.


###########################################

#FindSE.R

#Description: This is an R script used to find the super enhancer closest to a gene/miRNA

#Usage: 

	R CMD BATCH "--args $RNA_loci.csv $SE_loci.csv"
$RNA_loci.csv file contains the loci of a gene or miRNA on genome, see FileFormatDemo.txt.

$SE_loci.csv file contains the loci of a super enhancer, see FileFormatDemo.txt.

#Note:

the output gene/miRNA and super enhancer pairs need to be further corrected by fingdings of published studies(see method in paper).


###########################################

#CoverageCalculate.sh

#Description: This is a shell script to calculate the coverage of genes belonging to different exon length range.

#Usage: 

	CoverageCalculate.sh $1
$1 is the basename of fastq file of the sample used to calculate the coverage.

#Note: 

make sure that GTF file, genome index files, refexon_bed.txt file and CoverageCalculate.R file are properly prepared, details see the script and FileFormatDemo.txt.

The working directory needs to be the same with the fastq file's.

###########################################

#FindTopVariatedMiRNAandGenes.R

#Description: This is an R script used to find the top variated miRNAs or genes based on their expression with miRNA-mRNA dual sequencing data.

#Usage: 

	R CMD BATCH "--args $ReadCountMatrix_miRNA.csv $ReadCountMatrix_mRNA.csv $cut1 $cut2 $cut3 $cut4"
$ReadCountMatrix_miRNA.csv is a read count matrix file of miRNA, see FileFormatDemo.txt.

$ReadCountMatrix_mRNA.csv is a read count matrix file of mRNA, see FileFormatDemo.txt.

$cut1 and $cut2 are parameters used to select miRNAs that have no less than $cut1 mapped reads in no less than $cut2 samples.

$cut3 and $cut4 are parameters used to select mRNAs that have no less than $cut3 mapped reads in no less than $cut3 samples.


###########################################

#tSNEclustering.R

#Description: This is an R script used to do t-SNE clustering with PCA result as input.

#Usage: 

	R CMD BATCH "--args $ReadCountMatirx.csv $cut1 $cut2 $perplexity $groupsize1 $groupsize2 $groupsize3"
$ReadCountMatirx.csv is a read count matrix file of mRNA or miRNA, see FileFormatDemo.txt.

$cut1 and $cut2 are parameters used to select genes/miRNAs that have no less than $cut1 mapped reads in no less than $cut2 samples.

$perplexity is the perplexity parameter value used to do the t-SNE analysis.

$groupsize1, $groupsize2, $groupsize3 are the sample size of each group, presuming that these samples belong to three group.



###################
software version:
TopHat v2.0.11, bedtools v2.19.1, bowtie2 v2.2.2, samtools v0.1.19-44428cd, BBMap v35.85
###################

