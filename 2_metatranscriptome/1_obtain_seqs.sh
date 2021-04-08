# fetch the datasets with the fastq-dump command and save them in the newly created 1_sequence datasets folder
mkdir 1_sequence_datasets
cd 1_sequence_datasets

for i in `cat ../list`
do
	fastq-dump --split-files --gzip ${i}
done

cd ..


