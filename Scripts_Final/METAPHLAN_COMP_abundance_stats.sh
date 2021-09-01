### Metaphlan analysis with rel_ab_w_read_stats and --unknown_estimation

for sample in $(cat fastq_samples.txt)
do
	metaphlan ../results/metaphlan_bowtie/"$sample".bowtie2 \
	--nproc 5 --input_type bowtie2out -t rel_ab_w_read_stats \
	-o ../results/metaphlan/"$sample"_abundance_stats.txt \
	echo "Done for sample "$sample""
done

