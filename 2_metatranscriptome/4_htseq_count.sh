
#mkdir 4_combined_stringtie_prokka_lt5kfeat
mkdir -p 4_htseq_counts/strandedR
mkdir -p 4_htseq_counts/unstranded
mkdir -p 4_htseq_counts/stranded
# remove the sequence part from gff3... htseq-count does not like it
# save the number of lines till that point
FEATUREEND=`egrep -n -m 1 '##FASTA' ../1_metagenome/4_annotation/ANNOTATION_DIR/genome.gff | perl -pe "s/\:.+//g"`
# keep the feature lines in a new gff version
head -n ${FEATUREEND} ../1_metagenome/4_annotation/ANNOTATION_DIR/genome.gff > 4_htseq_counts/genome_noseqs.gff
# although the sequencing was stranded I have run all possible strandedness versions for validation (the correct one, as seen also in the results, was the reverse strandedR output)
for i in `cat list`
do
	htseq-count -m union -f bam --nonunique all -s reverse -t CDS -i ID ./3_star_seq_mapping/${i}/Aligned.out.sorted.bam ./4_htseq_counts/genome_noseqs.gff > ./4_htseq_counts/strandedR/counts_${i}.txt
	htseq-count -m union -f bam --nonunique all -s no -t CDS -i ID ./3_star_seq_mapping/${i}/Aligned.out.sorted.bam ./4_htseq_counts/genome_noseqs.gff > ./4_htseq_counts/unstranded/counts_${i}.txt
	htseq-count -m union -f bam --nonunique all -s yes -t CDS -i ID ./3_star_seq_mapping/${i}/Aligned.out.sorted.bam ./4_htseq_counts/genome_noseqs.gff > ./4_htseq_counts/stranded/counts_${i}.txt
done



# I repeated the above counting including all feature types including CDS... for doing so I had to create a new gff with amended feature types since htseq-count can count a single feature at the time (I named the features as MYFEATURE)
cat 4_htseq_counts/genome_noseqs.gff | perl -pe "s/\tCDS\t/\tMYFEATURE\t/g" | perl -pe "s/\ttRNA\t/\tMYFEATURE\t/g" | perl -pe "s/\ttmRNA\t/\tMYFEATURE\t/g" | perl -pe "s/\trRNA\t/\tMYFEATURE\t/g" | perl -pe "s/\trepeat_region\t/\tMYFEATURE\t/g" > 4_htseq_counts/genome_noseqs_amended_features.gff
# create the folders for saving these files
mkdir -p 4_htseq_counts/strandedR_all_features
mkdir -p 4_htseq_counts/unstranded_all_features
mkdir -p 4_htseq_counts/stranded_all_features
# although the sequencing was stranded I have run all possible strandedness versions for validation (the correct one, as seen also in the results, was the reverse strandedR output)
for i in `cat list`
do
	htseq-count -m union -f bam --nonunique all -s reverse -t MYFEATURE -i ID ./3_star_seq_mapping/${i}/Aligned.out.sorted.bam ./4_htseq_counts/genome_noseqs_amended_features.gff > ./4_htseq_counts/strandedR_all_features/counts_${i}.txt
	htseq-count -m union -f bam --nonunique all -s no -t MYFEATURE -i ID ./3_star_seq_mapping/${i}/Aligned.out.sorted.bam ./4_htseq_counts/genome_noseqs_amended_features.gff > ./4_htseq_counts/unstranded_all_features/counts_${i}.txt
	htseq-count -m union -f bam --nonunique all -s yes -t MYFEATURE -i ID ./3_star_seq_mapping/${i}/Aligned.out.sorted.bam ./4_htseq_counts/genome_noseqs_amended_features.gff > ./4_htseq_counts/stranded_all_features/counts_${i}.txt
done
