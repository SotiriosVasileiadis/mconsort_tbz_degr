# obtain the computer processors
MY_PROCS=`nproc`
# call the demultiplexing pipeline
./DemuxOwnBCsys_absPATH.sh demux_out /mnt/21B044BE29552220/sotos/TBZ/mconsort_tbz_degr/3_microbiomes/1_tbz_vs_suc/2_demultiplexed_sequences/R1.fastq /mnt/21B044BE29552220/sotos/TBZ/mconsort_tbz_degr/3_microbiomes/1_tbz_vs_suc/2_demultiplexed_sequences/R2.fastq bcfile.txt ${MY_PROCS}