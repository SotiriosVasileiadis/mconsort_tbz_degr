# obtain the computer processors
MY_PROCS=`nproc`
# call the demultiplexing pipeline
./DemuxOwnBCsys_absPATH.sh demux_out absolutepath/R1.fastq absolutepath/R2.fastq bcfile.txt ${MY_PROCS}