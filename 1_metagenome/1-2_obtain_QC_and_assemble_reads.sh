# fetch all short reads metegenome data
mkdir 1_sequence_datasets
cd 1_sequence_datasets
fastq-dump --split-files --gzip SRR7135606 
fastq-dump --split-files --gzip SRR7135607
fastq-dump --split-files --gzip SRR7135608
fastq-dump --split-files --gzip SRR7135609
fastq-dump --split-files --gzip SRR7135610
fastq-dump --split-files --gzip SRR7135611
fastq-dump --split-files --gzip SRR7135612
fastq-dump --split-files --gzip SRR14070200

# run nonpareil on the data (unfortunately it is not possible to run the command in compressed files)
for i in SRR7135606_1.fastq.gz SRR7135607_1.fastq.gz SRR7135608_1.fastq.gz SRR7135612_1.fastq.gz
do
	zcat ${i} > ${i%.*}
	nonpareil -s ${i%.*} -T kmer -f fastq -b nonpareil_${i%.*}
	rm ${i%.*}
done

Rscript nonpareil.r

cd ..

# combine all short reads and long reads with corresponding files
cat 1_sequence_datasets/SRR7135606_1.fastq.gz 1_sequence_datasets/SRR7135607_1.fastq.gz 1_sequence_datasets/SRR7135608_1.fastq.gz 1_sequence_datasets/SRR7135612_1.fastq.gz > 1_sequence_datasets/short1.fastq.gz
cat 1_sequence_datasets/SRR7135606_2.fastq.gz 1_sequence_datasets/SRR7135607_2.fastq.gz 1_sequence_datasets/SRR7135608_1.fastq.gz 1_sequence_datasets/SRR7135612_2.fastq.gz > 1_sequence_datasets/short2.fastq.gz
zcat 1_sequence_datasets/SRR7135609_1.fastq.gz 1_sequence_datasets/SRR7135610_1.fastq.gz 1_sequence_datasets/SRR7135611_1.fastq.gz 1_sequence_datasets/SRR14070200_1.fastq.gz > 1_sequence_datasets/long.fastq

# The assembly version at the contig names causes runtime issues and needs removing (with the piped perl oneliner)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/513/095/GCA_006513095.1_ASM651309v1/GCA_006513095.1_ASM651309v1_genomic.fna.gz ; gunzip -c *.fna.gz | perl -pe "s/\.1 / /" > 1_sequence_datasets/prev_assemb.fna


mkdir 2_assembly
cd 2_assembly
mkdir 1_trimmed_data

# quality control the sequences with fastp
fastp \
	-i ../1_sequence_datasets/short1.fastq.gz \
	-I ../1_sequence_datasets/short2.fastq.gz \
	-o 1_trimmed_data/R1_paired.fastq.gz \
	-O 1_trimmed_data/R2_paired.fastq.gz \
	--unpaired1 1_trimmed_data/R1_unpaired.fastq.gz \
	--unpaired2 1_trimmed_data/R2_unpaired.fastq.gz \
	--failed_out 1_trimmed_data/failed_out.fastq.gz \
	--detect_adapter_for_pe \
	-3 \
	--cut_tail_window_size 4 \
	--cut_tail_mean_quality 18 \
	-w 16 \
	-h 1_trimmed_data/QA_report.html \
	-j 1_trimmed_data/QA_report.json


# run OPERA-MS with the previous assembly as reference.
OPERA-MS.pl \
	--contig-file ../1_sequence_datasets/prev_assemb.fna \
	--short-read1 1_trimmed_data/R1_paired.fastq.gz \
	--short-read2 1_trimmed_data/R2_paired.fastq.gz \
	--long-read 1_sequence_datasets/long.fastq \
	--num-processors 56 \
	--no-gap-filling \ --polishing \
	--out-dir 2_opera_assembly_out 2> 2_opera_assemb_log.err
mv 2_opera_assemb_log.err 2_opera_assembly_out/



