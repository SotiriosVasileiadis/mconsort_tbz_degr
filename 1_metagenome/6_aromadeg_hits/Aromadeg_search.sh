# copy the BLAST preformated aromadeg db from the dbs folder and the predicted proetein sequences from the annotations file
cp -r ../z_accessory_files_dbs_and_executbles/aromadeg_search/* ./
cp ../4_annotation/ANNOTATION_DIR/protein.faa ./
# run BLAST for obtaining the best hit
blastp -query protein.faa -db aromadb/AromaDeg_all_unaligned.fasta -max_target_seqs 1 -evalue 1e-20 -outfmt 6 -num_threads 32 > blastouttbl.txt
# do the same for reviewgin the alignments and obtaining other close hits as well
blastp -query protein.faa -db aromadb/AromaDeg_all_unaligned.fasta -max_target_seqs 10 -evalue 1e-20 -num_threads 32 > blastout.txt 
