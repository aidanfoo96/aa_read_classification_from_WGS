# get results of group of samples and put into one text file

for sample in $(cat fastq_samples.txt)
do 
	samtools stats ../results/bam_aln/"$sample".bam > ../results/samtools_stats/"$sample".txt
done

# Parse to R

