## Extract unaligned reads!

for sample in $(cat fastq_samples.txt)
do
	samtools view -u -f 4 -F 264 ../results/bam_aln/"$sample"-raw.bam > ../results/bam_aln/"$sample"_tmps1.bam
	samtools view -u -f 8 -F 260 ../results/bam_aln/"$sample"-raw.bam > ../results/bam_aln/"$sample"_tmps2.bam
	samtools view -u -f 12 -F 256 ../results/bam_aln/"$sample"-raw.bam > ..results/bam_aln/"$sample"_tmps3.bam
	samtools merge -u ../results/merged_bam/"$sample"_merged.bam ../results/bam_aln/"$sample"_tmps1.bam \
	../results/bam_aln/"$sample"_tmps2.bam ..results/bam_aln/"$sample"_tmps3.bam 
	samtools sort -n ../results/merged_bam/"$sample"_merged.bam > ../results/merged_bam/"$sample"_sorted.bam
done 

echo "Done unaligned read extraction...now converting to paired fastq reads"

for sample in $(cat fastq_samples.txt)
do
	bamToFastq -i ../results/merged_bam/"$sample"_sorted.bam \
	-fq ../results/paired_unaligned_fastq/"$sample"_r1.fastq \
	-fq2 ../results/paried_unaligned_fastq/"$sample"r2.fastq
done

