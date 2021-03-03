mkdir 4_annotation
cd 4_annotation
cp ../3_binning/3_binning/dastool/contigs_binned.fasta ./

dfast -g contigs_binned.fasta \
--database ../z_accessory_files_dbs_and_executbles/AromadegDB/AromadegDfastDB.txt \
--references ../z_accessory_files_dbs_and_executbles/final_bins_annotated_sel_l500_cv5_final.gbf \
--config ../z_accessory_files_dbs_and_executbles/config_wth_locus_DFAST.py \
--cpu 32 \
--force

cd ..

## search the SEED database according to https://github.com/transcript/subsystemshttps://github.com/transcript/subsystems 
cd 4_annotation
mkdir SEED_annotation
cd SEED_annotation
diamond blastp --max-target-seqs 1 -d /home/sotosv/dbs/figfams/dbbuild/db_files/notabs.subsystems.complex.merged.reduced.db.dmnd -q ../ANNOTATION_DIR/protein.faa -o myresults.m8
# the Figs extraction script expects the figs at column 3... for this reason I have run the following two commands before extracting the annotations
cut -f1 myresults.m8 > myresults.m8_ids
paste myresults.m8_ids myresults.m8 > myresults.m8_wth_ids

## extract the annotations with
python2.7 /home/sotosv/dbs/figfams/dbbuild/subsystems/fig_swapper.py /home/sotosv/dbs/figfams/dbbuild/db_files/subsystems2peg myresults.m8_wth_ids /home/sotosv/dbs/figfams/dbbuild/db_files/subsystems2role /home/sotosv/dbs/figfams/dbbuild/db_files/subsys.txt
cd ../../

## run the plasmidome annotation with plasflow as well... for installation and running instructions see https://github.com/smaegol/PlasFlow
cd 4_annotation
mkdir plasflow
cd plasflow
# activate the environment
source activate plasflow # you can use also "conda activate plaflow"
cp ../ANNOTATION_DIR/genome.fna ./
filter_sequences_by_length.pl -input genome.fna -output filtered_genome.fna -thresh 1000
PlasFlow.py --input filtered_genome.fna --output filtered_genome_plas_pred.txt
source activate plasflow
cd ../../


## annotation was performed also with the emapper by submitting the sequences to the up-to-date http://eggnog-mapper.embl.de
# this could be performed locally with the following command (used in the actual annotation process as described in the info.txt document)
# python2 "$EMAPPERPATH"/emapper.py --cpu "6" -i "/data/shared/emapper_jobs/user_data/MM_meyjcu7j/query_seqs.fa" --output "query_seqs.fa" --output_dir "/data/shared/emapper_jobs/user_data/MM_meyjcu7j" -m "diamond" -d "none" --tax_scope "auto"        --go_evidence "non-electronic" --target_orthologs "all" --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query-cover 20 --subject-cover 0 --override --temp_dir "/data/shared/emapper_jobs/user_data/MM_meyjcu7j"
# however, due to the possible issues with python 2 which is nowadays being abandoned, the online approach was preferred