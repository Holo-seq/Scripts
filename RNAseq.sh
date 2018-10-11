#!/bin/bash
###this script is used to map raw Hiseq RNA-seq data to mm9 or hg19 via tophat
###<GTF file> A set of gene model annotations and/or known transcripts of mm9 or hg19
###<genome_index_base> The basename of the mm9 or hg19 genome index to be searched
###<CRG-flow> A C++ software used to count unique and multiple mapped reads to each gene
###<gene_exon_loci> a bed file of which each line contains the information of the start and end loci of each gene, as well as the startsite and length of every exon belong to these genes(mm9 and hg19)
###$1 is the base name of input fastq file, $1.tagcount is the output read count file, in which the first column is gene No, the second column is the number of unique mapped reads to each gene, the third column is the multiple mapped reads to each gene. 


tophat -o ./$1 -G <GTF file>  -p 10 --no-coverage-search <genome_index_base> $1.fastq &>$1_tophat.log
samtools view $1/accepted_hits.bam | <CRG-Flow> <gene_exon_loci.txt> $1.tagcount >$1.count

