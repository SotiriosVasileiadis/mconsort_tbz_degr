mkdir strandedness_test
# this run was performed to identify the strandedness of features for which the start-stop signatures were missing and were identified according to blast results only and not structural gene information
# for this run I have changed the NovosphingSphing which was annotated according to NCBI BLAST back to the Agrobacterium contig name, since the merged bam file I originally created had this annotation which meant that these counts would be lost if I did not change this back and the changing procedure would take 3-4 hours.
# run the stranded version
htseq-count -m intersection-nonempty -f bam -s yes -t gene -i ID ../merged_bams/all_Aligned.out.sorted.insert_lt_601_gt_neg601.bam ../../metawatt/0_append_annot_final_files/final_bins_annotated_sel_l500_cv5_final_no_novsph.gff > strandedness_test/counts_merged_stranded.txt
# run the unstranded_version
htseq-count -m intersection-nonempty -f bam -s no -t gene -i ID ../merged_bams/all_Aligned.out.sorted.insert_lt_601_gt_neg601.bam ../../metawatt/0_append_annot_final_files/final_bins_annotated_sel_l500_cv5_final_no_novsph.gff > strandedness_test/counts_merged_unstranded.txt
# after that compare the count differences in R and check the feature orientation and verify with tablet or IGV the orientation of the mapping of the RNA seq data