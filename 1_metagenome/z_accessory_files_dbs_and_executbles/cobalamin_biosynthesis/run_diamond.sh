# download the latest Uniprot SwissProt database
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
# run diamond
diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta.dmnd
diamond blastp -d uniprot_sprot.fasta.dmnd -q ../../4_annotation/ANNOTATION_DIR/protein.faa -o uniprot_sprot.fasta.dmnd.m8_1e-10_qcov50 -f 6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore -k 1 -p 32 --evalue 1e-10 --query-cover 50