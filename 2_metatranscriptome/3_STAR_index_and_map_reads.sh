# prepare the file structure 
mkdir 3_star_seq_mapping
cd 3_star_seq_mapping

# prepare the reference sequence for mapping
mkdir ref_forRNAseq
cp ../../1_metagenome/4_annotation/ANNOTATION_DIR/genome.fna ./ref_forRNAseq/
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ./ref_forRNAseq/ --genomeFastaFiles ./ref_forRNAseq/genome.fna --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 12

# start the transciript mapping
for i in `cat ../list`
do
	# WARNING!!! do not provide the gtf or the (meta-)genome fasta (again) here because it creates duplicate contig headers and all reads map as duplicates
	mkdir ${i}; cd ${i}
	STAR --runThreadN 32 \
	--genomeDir ../ref_forRNAseq/ \
	--readFilesIn ../../2_trimmed_data/${i}_R1_paired.fastq.gz ../../2_trimmed_data/${i}_R2_paired.fastq.gz \
	--readFilesCommand zcat

	samtools view -Sb Aligned.out.sam > Aligned.out.bam
	samtools sort -o Aligned.out.sorted.bam Aligned.out.bam
	samtools index Aligned.out.sorted.bam
	# remove the un-sorted/indexed sam and bam file
	rm Aligned.out.sam Aligned.out.bam
	cd ../
done

