#!/bin/bash

trimmomatic=$(find $HOME -name trimmomatic-0.39.jar)
adapters=$(find $HOME -type d -name "Trimmomatic-0.39")
adapters="${adapters}/adapters"
PICARD=$(find $HOME -name picard.jar)

read_length=42

read1="raw-data/SRR10044669_GSM4043375_GATA1_CD34_rep1_CD34_differentiated_Homo_sapiens_OTHER_1.fastq.gz"
read2="${read1%_1.fastq.gz}_2.fastq.gz"

read1o="trim/${read1##*/}"
read2o="trim/${read2##*/}"

output_read1_paired="${read1o%_1.fastq.gz}_paired_1.fastq.gz"
output_read1_unpaired="${read1o%_1.fastq.gz}_unpaired_1.fastq.gz"


output_read2_paired="${read2o%_2.fastq.gz}_paired_2.fastq.gz"
output_read2_unpaired="${read2o%_2.fastq.gz}_unpaired_2.fastq.gz"


mkdir -p trim

#first trimming
java -jar "$trimmomatic" PE -threads 40 -phred33 "$read1" "$read2" "$output_read1_paired" "$output_read1_unpaired" "$output_read2_paired" "$output_read2_unpaired" ILLUMINACLIP:"$adapters/TruSeq3-PE.fa":2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

#second trimming:
#kseq_test has to be manually compiled

mkdir -p trim2
output_read1_paired2="trim2/${output_read1_paired##*/}"
output_read2_paired2="trim2/${output_read2_paired##*/}"

kseq_test "$output_read1_paired" "$read_length" "$output_read1_paired2"
kseq_test "$output_read2_paired" "$read_length" "$output_read2_paired2"

#alignment

mkdir -p bam
bowtie2_output="bam/SRR10044669_GSM4043375_GATA1_CD34_rep1_CD34_differentiated_Homo_sapiens_OTHER.bam"
bowtie2 -p 42 --dovetail --phred33 -x "/home/niek/Documents/references/bowtie2-index/hg19/bowtie2-hg19" -1 "$output_read1_paired2" -2 "$output_read2_paired2" 2> align.log | samtools view -bS - > "$bowtie2_output"

mkdir -p {sorted,dup.marked,dedup}

#sorting bam files
java -jar "$PICARD" SortSam INPUT="$bowtie2_output" OUTPUT="sorted/$bowtie2_output" SORT_ORDER=coordinate

#marking duplicates
java -jar "$PICARD" MarkDuplicates INPUT="sorted/$bowtie2_output" OUTPUT="dup.marked/$bowtie2_output" METRICS_FILE=metrics.txt

#removing duplicates
java -jar "$PICARD" MarkDuplicates INPUT="sorted/$bowtie2_output" OUTPUT="dedup/$bowtie2_output" METRICS_FILE=metrics."$bowtie2_output".txt REMOVE_DUPLICATES=true

mkdir -p {sorted.120bp dup.marked.120bp dedup.120bp}



