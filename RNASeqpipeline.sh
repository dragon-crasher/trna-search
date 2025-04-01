#!/bin/bash

echo "current shell: $SHELL"
SECONDS=0

#change working directory
cd /mnt/d/bioinformatics/RNAseq_pipeline/

#STEP 1: Run fastqc

#fastqc data/demo.fastq -o data/

#run trimmomatic to trim reads with poor quality
#java -jar /mnt/d/bioinformatics/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 data/demo.fastq data/demo_trimmed.fastq TRAILING:10 -phred33
#echo "Trimmommatic finished running!"

#fastqc data/demo_trimmed.fastq -o data/

#STEP 2: Run HISAT2

#mkdir HISAT2
#get the genome indices
#wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

#run alignment 
hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U data/demo_trimmed.fastq | samtools sort -o HISAT2/demo_trimmed.bam
#STEP 3: Run featureCounts - Quantification

# get gtf
#wget ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
featureCounts -s 2 -a mnt/d/bioinformatics/RNAseq_pipeline/hg38/Homo_sapiens.GRCh38.113.chr_patch_hapl_scaff.gtf.gz -o quants/demo_featurecounts.txt HISAT2/demo_trimmed.bam

duration=$SECONDS
echo "$(($duration/60)) minutes and $(($duration % 60)) seconds elapsed."
