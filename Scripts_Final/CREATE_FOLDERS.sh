# Create Directories and samples.txt file for analysis in 'working' directory

create_working_environment(){
	create_folders
	create_samples.txt
}

create_folders(){
	echo "Creating sub-folders ..."
	mkdir -p working
	mkdir -p results/{alignment,bam_aln,bam_index,mapping,paired_unaligned_fastq,sam_aln,samtools_stats,unaligned_reads,unmapped_bam,unmapped_fastq, kraken2_out, metaphlan_bowtie, metaphlan}
	mkdir -p data/{aedes_ref,fastq}
	echo "DONE creating folders and sub-folders"
}

create_samples.txt(){
	echo "Creating fastq_samples.txt file..."
	ls ../data/*gz | cut -f3 -d "/" | cut -d "." -f1 | uniq > fastq_samples.txt
	echo "DONE"
}

