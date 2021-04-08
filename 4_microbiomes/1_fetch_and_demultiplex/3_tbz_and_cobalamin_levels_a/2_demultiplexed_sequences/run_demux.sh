# obtain the computer processors
MY_PROCS=20
# call the demultiplexing pipeline
./DemuxOwnBCsys_absPATH.sh demux_out /mnt/21B044BE29552220/sotos/TBZ/mconsort_tbz_degr/4_microbiomes/3_tbz_and_cobalamin_levels_a/2_demultiplexed_sequences/R1.fastq /mnt/21B044BE29552220/sotos/TBZ/mconsort_tbz_degr/4_microbiomes/3_tbz_and_cobalamin_levels_a/2_demultiplexed_sequences/R2.fastq bcfile.txt ${MY_PROCS}