# test for the best hits between the binned contigs (use the annotated ones) and the plasmid sequence of Inoue et al 2004 containing all the car related genes of Sphingomonas sp. KA1 plasmid pCAR3
# prepare the genome database
cp ../4_annotation/ANNOTATION_DIR/genome.fna ../z_accessory_files_dbs_and_executbles/
makeblastdb -in ../z_accessory_files_dbs_and_executbles/genome.fna -dbtype nucl

# run the blast to find the highly homologous hits and sort out the most important (use UGENE to visualize)
blastn -query ../z_accessory_files_dbs_and_executbles/AB270530.1.fasta -db ../z_accessory_files_dbs_and_executbles/genome.fna -outfmt 6 -evalue 1e-20 -max_target_seqs 1 > blast.txt

# the BLAST output suggests that only contig groopm_bin_50___opera_contig_243 contains associated hits
# I have therefore extracted the contig in fasta and genbank format

# perform the blast to be used in the genoplotR
makeblastdb -in ../z_accessory_files_dbs_and_executbles/AB270530.1.fasta -dbtype nucl
blastn -query ../z_accessory_files_dbs_and_executbles/groopm_bin_50___opera_contig_243.fna -db ../z_accessory_files_dbs_and_executbles/AB270530.1.fasta -outfmt 6 -evalue 1e-20 > groopm_bin_50___opera_contig_243.blout.txt
