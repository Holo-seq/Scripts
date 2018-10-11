#!/bin/bash
#this script is used to calculate the coverage of genes belonging to different exon length range
###
###<GTF file> A set of gene model annotations and/or known transcripts of mm9 or hg19
###<genome_index_base> <tsRNA_index_base> <miRNA_index_base> <snoRNA_index_base> The basename of the genome index, tsRNA reference index, miRNA reference index, snoRNA reference index of mm9 or hg19 to be searched
###<CoverageCalculate.R> The R scripts to calculate the coverage for genes of different exon length range, please note that it needs different 
###<refexon_loci.txt> is the bed file with every exon loci of every gene in each line, see FileFormatDemo.txt

#!/bin/bash
tophat -o ./$1 -G <GTF file>  -p 10 --no-coverage-search <genome_index_base> $1.fastq &>$1_tophat.log
cd $1

samtools view accepted_hits.bam | grep "NH:i:1$" accepted_hits.sam > accepted_hits_unique.sam
samtools view -H accepted_hits.bam > accepted_hits.header
cat accepted_hits.header accepted_hits_unique.sam | samtools view -S -b accepted_hits_unique_withheader.sam > accepted_hits_unique.bam
bedtools bamTobed -tag NH -i accepted_hits_unique.bam > $1.bed

R CMD BATCH "--args $1 <refexon_loci.txt>" <CoverageCalculate.R> CoverageCalculate.Rout



