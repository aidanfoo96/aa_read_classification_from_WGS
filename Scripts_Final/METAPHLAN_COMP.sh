# Script for Taxaonmic Composition of Unmapped Reads 

for sample in $(fastq_samples.txt)
do 
	metaphlan ../results/paired_unaligned_fastq/"$sample"_unmapped_r1.fastq,..results/paired_unaligned_fastq/"$sample"_unmapped_r2.fastq\
	--bowtie2out ..results/metaphlan/"$sample".bowtie2.bz2 --nproc 5 --input_type fastq -o ../results/metaphlan/"$sample"_prof.txt
done

merge_metaphlan_tables.py ../results/metaphlan/*_prof.txt > merged_abundance_table.txt

metaphlan ../results/paired_unaligned_fastq/SRR1523067_unmapped_r1.fastq,../results/paired_unaligned_fastq/SRR1523067_unmapped_r2.fastq --bowtie2out ../results/metaphlan/SRR1523067.bowtie2.bz2 --nproc 5 --input_type fastq -o SRR1523067_prof.txt