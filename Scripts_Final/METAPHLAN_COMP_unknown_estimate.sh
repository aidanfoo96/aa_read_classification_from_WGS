for sample in $(cat fastq_samples.txt)
do
	metaphlan ../results/metaphlan_bowtie/"$sample".bowtie2 \
	--nproc 5 --input_type bowtie2out -t --unknown_estimation \
	-o ../results/metaphlan/"$sample"_unknown_estimate.txt \
	echo "Done for sample "$sample""
done
