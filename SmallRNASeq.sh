#!/bin/bash
#this script is used to call tsRNA, miRNA and snoRNA from small RNA libraries
###
###<GTF file> A set of gene model annotations and/or known transcripts of mm9 or hg19
###<genome_index_base> <tsRNA_index_base> <miRNA_index_base> <snoRNA_index_base> The basename of the genome index, tsRNA reference index, miRNA reference index, snoRNA reference index of mm9 or hg19 to be searched
###<miRNA_preMiRNA_loci> The file contains miRNAs and their start&end loci in their corresponding pre-miRNAs
###<CountMatureMirna.R> The R script to count mature miRNA
###$1 is the base name of input fastq file, *.count files are small RNA count files


tophat -o ./$1 -G <GTF file>  -p 10 --no-coverage-search <genome_index_base> $1.fastq &>$1_tophat.log
samtools bam2fq $1/unmapped.bam > unmapped.fq

bbduk.sh in=unmapped.fq literal=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC out=adapter.unmapped.fq outm=unmapped.adapter.fq k=10 overwrite=t qtrim=f

##tsRNA
awk '{if((NR%4==2)||(NR%4==0)) print substr($1,1,40)}' unmapped.adapter.fq > unmapped.trimmed_40bp.fq
bowtie2 --local -N 1 -L 30 -p 2 -x <tsRNA_index_base> -U unmapped.trimmed_40bp.fq -S $1_tsRNA40bp.sam --no-unal --un unmapped.40bp.fq
awk '{if((NR%4==2)||(NR%4==0)) print substr($1,1,30)}' unmapped.40bp.fq > unmapped.trimmed.30bp.fq
bowtie2 --local -N 1 -L 18 -p 2 -x <tsRNA_index_base> -U unmapped.trimmed.30bp.fq -S $1_tsRNA30bp.sam --no-unal 
cat $1_tsRNA40bp.sam $1_tsRNA30bp.sam > $1_tsRNA.sam
awk '{if ($4>0) print $3}' $1_tsRNA.sam | sort | uniq -c > $1_tsRNA.count

##miRNA
awk '{if((NR%4==2)||(NR%4==0)) print substr($1,1,30)}' unmapped_adapter.fq > unmapped.trimmed.fq
bowtie2 --local -N 1 -L 16 -p 2 -x <pre-miRNA_index_base> -U unmapped.trimmed.fq -S $1_preMiRNA.sam
awk '{if(match($6,"17M")||match($6,"18M")||match($6,"19M")||match($6,"20M")||match($6,"21M")||match($6,"22M")||match($6,"23M")||match($6,"24M")||match($6,"25M"))print }' $1_preMiRNA.sam > $1_preMiRNA.filtered.sam
awk '{print $3,$4}' $1_preMiRNA.filtered.sam > $1_extMiRNA.filtered.txt
R CMD BATCH "--args $1_extMiRNA.filtered.txt <miRNA_preMiRNA_loci.txt> $1_matureMiRNA.count " <CountMatureMirna.R> CountMatureMirna.Rout


##snoRNA
bowtie2 -x <snoRNA_index_base> -N 1 -p 2 $1.fq -S $1_snoRNA.sam --no-unal
awk '{if ($4>0) print $3}' $1_snoRNA.sam | sort | uniq -c > $1_snoRNA.count

sed -i '/bowtie/d' *.count
sed -i "s/^ *//g" *.count
sed -i "s/ /\t/g" *.count
rm unmapped*

