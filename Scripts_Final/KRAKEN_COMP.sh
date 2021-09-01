## Perform kraken2 taxonomic classification accross all samples

for sample in $(cat fastq_samples.txt)
do
	kraken2 --db /home/db/Kraken2_db/ --threads 5 --output ../results/kraken2_out/"$sample"_kraken_run \
	--report ../results/kraken2_out/"$sample"_report.txt --use-mpa-style \
	--paired ../results/paired_unaligned_fastq/"$sample"_r1.fastq ../results/paired_unaligned_fastq/"$sample"_r2.fastq
done