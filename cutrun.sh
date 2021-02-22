#!/bin/bash

### Author: Niek Wit (University of Cambridge) 2021 ###

SCRIPT_DIR=$(find $HOME -type d -name "CUT-RUN")
threads=""
working_dir=$(pwd)
align_mm=0
rename_config=""
genome=""
PICARD=$(find $HOME -name picard.jar)
dedup=""

usage() {                                    
	echo "Usage: $0 [-g <genome build>] [-r OPTIONAL:renames NGS files] [-m <INT> mismatches allowed for alignment (standard is zero) OPTIONAL] [-t <INT> number of CPU threads to be used]"
	exit 2
}

while getopts 't:g:rdm:?h' c
do
	case $c in
		t) 
			threads=$OPTARG;;
		g) 
			genome="$OPTARG";;
		r)	
			rename_config="rename.config";;
		d)	
			dedup="yes";;
		m)  
			align_mm=$OPTARG;;	
		h|?) 	
			usage;;
	esac
done

if [[ $align_mm == 0 ]] || [[ $align_mm == 1 ]];
	then
		:
	else
		echo "ERROR: -m parameter must be 0 or 1"
		usage
		exit 1
fi

if [[ -z "$threads" ]];
	then
		threads=$(nproc --all)
fi

#renames files if -r flag and rename template have been given
if [[ ! -z "$rename_config" ]];
	then
		input=$rename_config
		while IFS= read -r line
		do
			original_file=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into original file name
			new_file=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into new file name
			mv "raw-data/${original_file}" "raw-data/${new_file}"
		done < "$input"
fi

#checks if file extension is .fq or .fastq
input_format=$(ls raw-data/*.gz | head -1)
if [[ $input_format == *"fastq"* ]]
	then
  		input_extension="fastq.gz"
elif [[ $input_format == *"fq"* ]]
	then
  		input_extension="fq.gz"		
fi

#Fastq quality control
fastqc_folder="fastqc/"
if [[ ! -d  "$fastqc_folder" ]]; 
	then
		echo "Performing FastQC"
		mkdir fastqc
		fastqc -q --threads "$threads" -o fastqc/ raw-data/*$input_extension
		echo "Performing MultiQC"
		multiqc -o fastqc/ fastqc/ . 2>> multiqc.log
fi

#trimming samples
trim_folder="trim/"
if [[ ! -d  "$trim_folder" ]];
	then
		echo "Performing read trimming"
		mkdir -p trim
		for read1 in raw-data/*"_1.$input_extension"
		do 
			echo $read1 >> trim.log
			read2="${read1%_1."$input_extension"}_2."$input_extension""
			trim_galore -j 4 -o ./trim --paired $read1 $read2 2>> trim.log
		done
fi

#aligning samples (bowtie2 settings:Skene et al 2018 Nature Protocols)
index_path=$(cat "$SCRIPT_DIR/config.yml" | shyaml get-value $genome.index_path)
#black_list_path=$(cat "$SCRIPT_DIR/config.yml" | shyaml get-value $genome.black_list_path)

bam_folder="bam/"
if [[ ! -d  "$bam_folder" ]];
	then
		mkdir -p bam
		touch align.log
		for read1 in trim/*"_1_val_1.fq.gz"
		do 
			read2="${read1%_1_val_1.fq.gz}_2_val_2.fq.gz"
			file_name="${read1##*/}"
			base_name="${file_name%_1_val_1.fq.gz}"
			echo "$base_name" >> align.log
			align_output="bam/$base_name.bam"
			echo "Aligning $base_name"
			bowtie2 -p "$threads" --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x "$index_path" -1 "$read1" -2 "$read2" 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$threads" - | samtools sort -@ "$threads" - > "$align_output"
		done
fi                                                                                                                                      

#deduplication
if [[ "$dedup" == "yes" ]];
	then
		mkdir -p deduplication
		for bam in bam/*.bam
		do
			echo "Removing duplicates $bam"
			dedup_output="${file%.bam}-dedup.bam"
			dedup_output="deduplication/${dedup_output##*/}"
			java -jar $PICARD MarkDuplicates INPUT=$file OUTPUT=$dedup_output REMOVE_DUPLICATES=TRUE METRICS_FILE=$dedup_output-metric.txt 2>> deduplication.log
		done
fi


