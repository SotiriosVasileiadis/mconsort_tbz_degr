# fetch the datasets with the fastq-dump command and save them in the newly created 1_sequence datasets folder
mkdir 1_sequence_datasets
cd 1_sequence_datasets

for i in `cat ../list`
do
	fastq-dump --split-files --gzip ${i}
done

cd ..


# quality trim the data with the fastp algorithm as follows
mkdir 2_qc_trimmed_data

for i in `cat ../list`
do
	fastp -i 1_sequence_datasets/${i}_1.fastq.gz -I 1_sequence_datasets/${i}_2.fastq.gz -o 2_qc_trimmed_data/${i}_R1_paired.fastq.gz -O 2_qc_trimmed_data/${i}_R2_paired.fastq.gz --unpaired1 2_qc_trimmed_data/${i}_R1_unpaired.fastq.gz --unpaired2 2_qc_trimmed_data/${i}_R2_unpaired.fastq.gz --failed_out 2_qc_trimmed_data/${i}_failed_out.fastq.gz --detect_adapter_for_pe -3 --cut_tail_window_size 4 --cut_tail_mean_quality 18 -w 16 -h 2_qc_trimmed_data/${i}_QA_report.html -j 2_qc_trimmed_data/${i}_QA_report.json
done

