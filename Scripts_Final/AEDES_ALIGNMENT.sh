## Align fastq files to aedes_reference 

do_alignment()
for sample in $(cat fastq_sample.txt)
do 
	hisat2 -x ../data/aedes_ref/aegypti_ref\
	-1 ../data/fastq/"$sample"_1.fastq.gz\
	-2 ../data/fastq/"$sample"_2.fastq.gz\
	-p 4 -S ../results/alignment/"$sample".sam -t
done

echo "Done Alignment...moving on to samtools processing"

## Convert to bam, sort, index
for sample in $(cat fastq_sample.txt)
do 
	samtools view -b ../results/alignment/"$sample".sam > ../results/bam_aln/"$sample"-raw.bam 
	samtools sort ../results/bam_aln/"$sample"-raw.bam > ../results/bam_index/"$sample".bam 
	samtools index ../results/bam_index/"$sample".bam
done

echo "Finished conversion, sorting and indexing"

