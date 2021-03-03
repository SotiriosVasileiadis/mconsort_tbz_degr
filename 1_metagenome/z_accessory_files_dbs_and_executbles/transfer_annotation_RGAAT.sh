# convert to the ratt-desired embl format transferring also the functional annotation
#seqret -sequence Fusarium_solani_FsK.gbk -sformat genbank -outseq embl/FsK.embl -osformat embl -feature Y
#seqret -sequence Fusarium_solani_FsK.gbk -sformat genbank -outseq embl/FsK.gff -osformat gff -feature Y
#seqret -sequence Fusarium_solani_FsK.gbk -sformat genbank -outseq FsK.fasta -osformat fasta

#cp ../pilon_polishing/pilon_corrected/fs.dmo.cns.pilon.fasta fs.dmo.cns.pilon.fasta
# run RGAAT for the annotation transferring
RGAAT.pl -g final_bins_annotated_sel_l500_cv5_final.fsa -a gff final_bins_annotated_sel_l500_cv5_final.gff -n contigs.polished.fasta -o rgaat_out

# map the old contigs to the new ones... this is just a test
#minimap2 -ax map-pd -t 56 contigs.polished.fasta final_bins_annotated_sel_l500_cv5_final.fsa > mapped_old_to_new_assemb.sam
# option "--cs" is recommended as paftools.js may need it
minimap2 -cx asm5 -t 56 --cs final_bins_annotated_sel_l500_cv5_final.fsa contigs.polished.fasta > mapped_old_to_new_assemb_cs.sam

# prep the annotation files
#cd embl
#awk 'BEGIN{RS="\n?ID   sc"} (NR-1){print "ID   sc" $0 > ("output-file_" NR ".embl")}' FsK.embl
#mv Fsk.embl ../
#cd ..

# run the transfer with flo you need to place also the relevant fasta and the gff files in the working directory and make the necessary changes on the flo_opts.yaml file concerning these files
# the result files will be placed in the run directory which needs to be removed avery time the software is run
mkdir flo
cd flo
cp ~/software/annotation/flo/flo/opts_example.yaml flo_opts.yaml 
rake -f ~/software/annotation/flo/flo/Rakefile

# I have cleaned the gff from duplicated attributes (genome tools don't like that)
GFFcleaner FsK.gff2
# I have tidied up the gff3 file using the genome tools
gt gff3 -sort -tidy -retainids FsK_clean.gff > FsK_clean_tidy.gff
# then I run the flo executable
rake -f ~/software/annotation/flo/flo/Rakefile
### the folders "run" and "run_clean_tidy" contain the same output files
















######## the rest did not work #######


# I could (haven't done it cause it provides weird behaviour... removes CDSs along) run remove the gene or mRNA annotation which interferes with the transcripts (the gene feature characterization is used instead of mRNA)... prepare the gene-less gff, prepare the options file and run again
# you can find the tip at https://github.com/wurmlab/flo
# I have replaced gaps with underscores since the software introduces tabs in the case of gaps
perl -pe "s/ /_/g" FsK.gff2 > FsK_nogap.gff2
perl -pi -e "s/_/ /g if /^##sequence-region/" FsK_nogap.gff2
grep -v "biological_region" FsK_nogap.gff2 > FsK_nogap_nbr.gff2
~/software/annotation/flo/flo/gff_remove_feats.rb mRNA FsK_nogap_nbr.gff2 > FsK_nogap_nbr_nogene.gff2
gt gff3 -sort -tidy -retainids FsK_nogap_nbr_nogene.gff2 > FsK_nogap_nbr_nogene_tidy.gff2


# run the ratt
#cp ../pilon_polishing/pilon_corrected/fs.dmo.cns.pilon.fasta fs.dmo.cns.pilon.fasta
# convert the genbank format to gff3 and that to the embl format necessary for running the ratt suite (direct convertsion with seqret did not work)
##keep in mind that the molecule type is not provided in te gff file and this creates problems with the RATTwithGFF notrun below because this information has to be added interactively
##EMBLmyGFF3 FsK.gff FsK.fasta -o embl/FsK.embl -p PRJ123456 -i "FUN" -s "Fusarium solani" -t linear -r 1
cd embl
awk 'BEGIN{RS="\n?ID   sc"} (NR-1){print "ID   sc" $0 > ("output-file_" NR ".embl")}' FsK.embl

awk 'BEGIN{RS="\n?>"} (NR-1){print ">" $0 > ("query_" NR ".fasta")}' ../fs.dmo.cns.pilon.fasta
mkdir run_ratt
cd run_ratt
start.ratt.sh ../embl ../fs.dmo.cns.pilon.fasta Fsk_nanopore_polished_ratted Assembly
#start.ratt.sh embl fs.dmo.cns.pilon.fasta Fsk_nanopore_polished_ratted Assembly
#RATTwithGFF.py FsK.gff FsK.fasta fs.dmo.cns.pilon.fasta refToQuery Assembly
