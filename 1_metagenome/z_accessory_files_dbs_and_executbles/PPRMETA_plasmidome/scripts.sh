# used the PPR-meta from https://doi.org/10.1093/gigascience/giz066
# the softare is looking for the executables in the current folder... So you need to copy them here before running
cp -r ~/software/annotation/plasmidome/PPR-Meta/* ./
./PPR_Meta contigs.polished.fasta contigs.polished.plasmidome.csv
