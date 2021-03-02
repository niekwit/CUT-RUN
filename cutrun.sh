#!/bin/bash

### Author: Niek Wit (University of Cambridge) 2021 

threads=""
rename_config=""
genome=""
dedup=""
cutoff=0
align_mm=0
big_wig=""
reads_length=""
peak=""
ngsplot=""

PICARD=$(find $HOME -name picard.jar)
SCRIPT_DIR=$(find $HOME -type d -name "CUT-RUN")
WORK_DIR=$(pwd)/

function usage {
echo "Usage: $0 [-g <genome build>] [-r OPTIONAL:renames NGS files] [-m <INT> mismatches allowed for alignment (standard is zero) OPTIONAL] [-t <INT> number of CPU threads to be used]"
exit 2
}

while getopts 't:g:rdm:f:blpn?h' c
do
	case $c in
		t)
			threads=$OPTARG;;
		g)
			genome="$OPTARG";;
		r)
			rename_config="rename.config";;
		d)
			dedup="TRUE";;
		m)
			align_mm=$OPTARG;;
		f)
			cutoff=$OPTARG;;
		b)
			big_wig="TRUE";;
		l)
			reads_length="TRUE";;
		p)
			peak="TRUE";;
		n)
			ngsplot="TRUE";;
		h|?)
			usage;;
	esac
done

###checking input variables
if [[ $align_mm == 0 ]] || [[ $align_mm == 1 ]];
	then
		:
	else
		echo "ERROR: -m parameter must be 0 or 1"
		usage
		exit 1
fi

if [[ ! $cutoff =~ ^-?[0-9]+$ ]];
	then
		echo "ERROR: filter cutoff should be an integer."
		exit 1
fi

if [[ -z "$threads" ]] 
	then
		threads=$(nproc --all)
elif [[ ! -z "$threads" ]] && [[ ! $threads =~ ^-?[0-9]+$ ]]
	then
		echo "ERROR: number of threads should be an integer."
		exit 1
fi

###renames files if -r flag, and rename template is present
if [[ ! -z "$rename_config" ]] && [[ -e "rename.config" ]];
	then
		input=$rename_config
		while IFS= read -r line
		do
			original_file=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into original file name
			new_file=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into new file name
			mv "raw-data/${original_file}" "raw-data/${new_file}"
		done < "$input"
fi

###checks if file extension is .fq or .fastq
input_format=$(ls raw-data/*.gz | head -1)
if [[ $input_format == *"fastq"* ]]
	then
  		input_extension="fastq.gz"
elif [[ $input_format == *"fq"* ]]
	then
  		input_extension="fq.gz"		
fi

###Fastq quality control
fastqc_folder="fastqc/"
if [[ ! -d "$fastqc_folder" ]]; 
	then
		echo "Performing FASTQC"
		mkdir fastqc
		fastqc -q --threads "$threads" -o fastqc/ raw-data/*$input_extension
		echo "Performing MultiQC"
		multiqc -o fastqc/ fastqc/ . 2>> multiqc.log
else
	echo "FASTQC/MultiQC already performed"
fi

###trimming samples
trim_folder="trim/"
if [[ ! -d "$trim_folder" ]];
	then
		echo "Performing read trimming"
		mkdir -p trim
		for read1 in raw-data/*"_1.$input_extension"
		do 
			echo $read1 >> trim.log
			read2="${read1%_1."$input_extension"}_2."$input_extension""
			trim_galore -j 4 -o ./trim --paired $read1 $read2 2>> trim.log
		done
else
	echo "Trimming already performed"	
fi

###aligning samples 
bam_folder="bam/"

if [[ ! -d "$bam_folder" ]];
	then
		mkdir -p bam
		#touch align.log
		index_path=$(cat "$SCRIPT_DIR/genome.yml" | shyaml get-value genome.$genome.0 | awk '{print$2}')
		#black_list_path=$(cat "$SCRIPT_DIR/genome.yml" | shyaml get-value genome.$genome.1 | awk '{print$2}')
		for read1 in trim/*"_1_val_1.fq.gz"
		do 
			read2="${read1%_1_val_1.fq.gz}_2_val_2.fq.gz"
			file_name="${read1##*/}"
			base_name="${file_name%_1_val_1.fq.gz}"
			echo "$base_name" >> align.log
			align_output="bam/$base_name.bam"
			echo "Aligning $base_name"
			bowtie2 -p "$threads" --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x "$index_path" -1 "$read1" -2 "$read2" 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$threads" - | samtools sort -@ "$threads" - > "$align_output"
			samtools index -b -@$threads "$align_output"
		done
else
	echo "Alignment already performed"
fi

#--dovetail --phred33
#--local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700

###deduplication
dedup_folder="deduplication/"
if [[ "$dedup" == "TRUE" ]] && [[ ! -d "$dedup_folder" ]];
	then
		mkdir -p deduplication
		touch deduplication.log
		for bam in bam/*.bam
		do
			echo "Removing duplicates $bam"
			dedup_output="${bam%.bam}-dedup.bam"
			dedup_output=deduplication/${dedup_output##*/}
			java -jar $PICARD MarkDuplicates INPUT="$bam" OUTPUT="$dedup_output" REMOVE_DUPLICATES=TRUE  METRICS_FILE=$dedup_output-metric.txt 2>> deduplication.log
			samtools index -b -@$threads "$dedup_output"
		done
else
    echo "Deduplication already performed/not selected"
fi

###generate histogram of reads lengths
#GSE125988

function reads_length {
#get frequency of read lengths
base_file="${file%.bam}"
base_file="${base_file##*/}"
out_file="$outdir/read-length-$base_file.txt"
echo "read_length" >> "$out_file"
samtools view -@$threads "$file" | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort >> "$out_file" 
#plot histogram of frequency of read lengths
Rscript $SCRIPT_DIR/read_length.R "$WORK_DIR" "$out_file" "$base_file" "$outdir"
}

reads_length_folder=read-lengths
reads_length_folder_dedup="$reads_length_folder"-dedup
if [[ "$reads_length" == "TRUE" ]] && [[ ! -d "$reads_length_folder" ]];
	then
		echo "Generating histograms of read lengths"
		mkdir "$reads_length_folder"
		outdir="$reads_length_folder"
		for file in bam/*.bam
		do 
			reads_length
		done
elif [[ "$reads_length" == "TRUE" ]] && [[ ! -d "$reads_length_folder_dedup" ]];
	then
		echo "Generating histograms of read lengths"
		mkdir "$reads_length_folder_dedup"
		outdir="$reads_length_folder_dedup"
		for file in $dedup_folder/*.bam
		do 
			reads_length
		done
else
		echo "Read length histograms already generated"
fi


###create bigWig files
#load bamCoverage settings:
if [[ "$big_wig" == "TRUE" ]];
	then
		binsize=$(cat settings.yaml | shyaml get-value bigWig.binSize)
		normalizeusing=$(cat settings.yaml | shyaml get-value bigWig.normalizeUsing)
		extendreads=$(cat settings.yaml | shyaml get-value bigWig.extendReads)
		effectivegenomesize=$(cat settings.yaml | shyaml get-value bigWig.effectiveGenomeSize)
		#checking bamCoverage settings:
		if [[ ! $binsize =~ ^-?[0-9]+$ ]];
			then
				echo "ERROR: binSize should be an integer."
				exit 1
		fi

		if [[ $normalizeusing != "RPKM" ]] && [[ $normalizeusing != "CPM" ]] && [[ $normalizeusing != "BPM" ]] && [[ $normalizeusing != "RPGC" ]] && [[ $normalizeusing != "None" ]];
			then
				echo "ERROR: invalid normalisation method chosen."
				echo "Available methods: RPKM, CPM, BPM, RPGC or None"
				exit 1
		fi

		if [[ ! $extendreads =~ ^-?[0-9]+$ ]];
			then
				echo "ERROR: extendReads should be an integer."
				exit 1
		fi

		if [[ ! $effectivegenomesize =~ ^-?[0-9]+$ ]];
			then
				echo "ERROR: effectiveGenomeSize should be an integer."
				echo "See https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html"
				exit 1
		fi

		#creates bigWig files
		bigwig_dir=bigwig_bin$binsize
		bigwig_dir_dedup=bigwig_dedup_bin$binsize
		
		function big_wig {
		#bigwig_output="${file%-dedupl-sort-bl.bam}-norm.bw"
		bigwig_output=${bigwig_output##*/}
		bamCoverage -p $threads --binSize "$binsize" --normalizeUsing "$normalizeusing" --extendReads "$extendreads" --effectiveGenomeSize "$effectivegenomesize" -b $file -o $bigwig_dir/$bigwig_output 2>> bigwig.log
		}  

		if [[ -d "$bam_folder" ]] && [[ ! -d "$bigwig_dir" ]];
			then
				echo "Creating BigWig files"
				mkdir -p "$bigwig_dir"
				for file in bam/*bam
				do 
					bigwig_output="${file%.bam}-norm.bw"
					big_wig
				done
		elif [[ -d "$dedup_folder" ]] && [[ ! -d "$bigwig_dir_dedup" ]];
			then
				echo "Creating BigWig files"
				mkdir -p "$bigwig_dir_dedup"
				bigwig_dir="$bigwig_dir_dedup"
				for file in "$dedup_folder"/*bam
				do 
					bigwig_output="${file%.bam}-dedup-norm.bw"
					big_wig
				done
		fi
fi


###filter out reads below cutoff
filter_folder="filter/"

function filter_reads {
	file_name="${bam##*/}"
	base_name=${file_name%.bam}
	samtools view -h -@ "$threads" "$bam" | awk -v x=$cutoff 'length($10) < x || $1 ~ /^@/' | samtools view -bS -@ "$threads" - > "filter/$base_name.${cutoff}bp.bam"
} #function to filter reads out below given cutoff

if [[ "$cutoff" == 0 ]];
	then
		echo "No size filtering of reads selected"
elif [[ ! -d "$filter_folder" ]] && [[ -d "$dedup_folder" ]];
	then
		mkdir filter
		echo "Filtering reads <${cutoff} bp"
		for bam in deduplication/*.bam
		do
			filter_reads
		done
elif [[ ! -d "$filter_folder" ]] && [[ ! -d "$dedup_folder" ]];
	then
		mkdir filter
		echo "Filtering reads <${cutoff} bp"
		for bam in bam/*.bam
		do
			filter_reads
		done
else
	echo "Size filtering already performed"
fi


###create metagene plots
ngsplot_folder="ngsplot"
ngsplot_folder_dedup="$ngsplot_folder-dedup"
feature=$(cat settings.yaml | shyaml get-value ngsplot.feature)
window=$(cat settings.yaml | shyaml get-value ngsplot.window)

function ngs_plot {
echo "Generating metagene plots and heatmaps"
for file in $bam_folder/*.bam
do
	case $file in 
		*input*) continue;;
		*) 
			ngs_output=${file%.bam}
			ngs_output=${ngs_output##*/}
			out_dir="$ngsplot_folder/$ngs_output"
			mkdir -p "$out_dir"
			ngs.plot.r -G "$genome" -R "$feature" -C $file -O "$out_dir/${ngs_output}_$feature" -T $ngs_output -L "$window" 
	esac
done
}

if [[ "$ngsplot" == "TRUE" ]];
	then
		#check input variables
		if [[ $feature != "tss" ]] && [[ $feature != "tes" ]] && [[ $feature != "genebody" ]] && [[ $feature != "exon" ]] && [[ $feature != "cgi" ]] && [[ $feature != "enhancer" ]] && [[ $feature != "dhs" ]] && [[ $feature != "bed" ]];
		then
			echo "ERROR: wrong genomic feature selected for ngs.plot.R"
			echo "Available options: tss, tes, genebody, exon, cgi, enhancer, dhs or bed"
			exit 1
		fi 
		if [[ ! $window =~ ^-?[0-9]+$ ]];
		then
			echo "ERROR: window size should be an integer for ngs.plot.R"
			exit 1
		fi
		
		#plot
		if [[ -d "$bam_folder" ]] && [[ ! -d "$ngsplot_folder" ]];
			then
				ngs_plot
		elif [[ -d "$dedup_folder" ]] && [[ ! -d "$ngsplot_folder_dedup" ]];
			then
				bam_folder=$dedup_folder
				ngsplot_folder=$ngsplot_folder_dedup
				ngs_plot
		fi
fi


###peak calling
if [[ "$peak" == "TRUE" ]];
	then
		echo "Calling and annotating peaks"
fi



