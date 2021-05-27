
mkdir 3_binning
cd 3_binning

mkdir 1_qc_controlled_reads

# quality control the short sequences (names provided in the "list" file) with fastp
for i in `cat list`
do
	fastp \
		-i ../1_sequence_datasets/${i}_1.fastq.gz \
		-I ../1_sequence_datasets/${i}_2.fastq.gz \
		-o 1_qc_controlled_reads/${i}_R1_paired.fastq.gz \
		-O 1_qc_controlled_reads/${i}_R2_paired.fastq.gz \
		--unpaired1 1_qc_controlled_reads/${i}_R1_unpaired.fastq.gz \
		--unpaired2 1_qc_controlled_reads/${i}_R2_unpaired.fastq.gz \
		--failed_out 1_qc_controlled_reads/${i}_failed_out.fastq.gz \
		--detect_adapter_for_pe \
		-3 \
		--cut_tail_window_size 4 \
		--cut_tail_mean_quality 18 \
		-w 16 \
		-h 1_qc_controlled_reads/${i}_QA_report.html \
		-j 1_qc_controlled_reads/${i}_QA_report.json
done

## mapr all reads against the OPERA-MS assembly
mkdir 2_read_mapping

# First map the short reads
cd 2_read_mapping
cp ../../2_assembly/2_opera_assembly_out/contigs.fasta ./
bowtie2-build contigs.fasta contigs.fasta.bt2i

for i in `cat ../list`
do
	bowtie2 -x contigs.fasta.bt2i -X 700 -1 ../1_qc_controlled_reads/${i}_R1_paired.fastq.gz -2 ../1_qc_controlled_reads/${i}_R2_paired.fastq.gz -S ${i}_out.sam -p 32
	samtools view -bS ${i}_out.sam > ${i}_out.bam
	samtools sort ${i}_out.bam -o ${i}_out.sorted.bam
	samtools index ${i}_out.sorted.bam
done

# This script was used for aligning the long reads to the reference and calculating the per contig sequence depth as required from solidbin
# nanopore reads
minimap2 -ax map-ont contigs.fasta ../../1_sequence_datasets/porechopped_nanofilt_fastq.fastq > aln_nanopo_out.sam
samtools view -bS aln_nanopo_out.sam > aln_nanopo_out.bam
samtools sort -o aln_nanopo_out.sorted.bam aln_nanopo_out.bam
samtools index aln_nanopo_out.sorted.bam

# pacbio reads (the three datasets derived from the same sample were merged)
cat ../../1_sequence_datasets/SRR71356{09,10,11}_1.fastq.gz > ../../1_sequence_datasets/pacbio.fastq.gz
minimap2 -ax map-pb contigs.fasta ../../1_sequence_datasets/pacbio.fastq.gz > aln_pacbio_out.sam
samtools view -bS aln_pacbio_out.sam > aln_pacbio_out.bam
samtools sort -o aln_pacbio_out.sorted.bam aln_pacbio_out.bam
samtools index aln_pacbio_out.sorted.bam

cd ..


## perform tha actual binning
mkdir 3_binning
cd 3_binning

## start with GroopM
# I had to create a docker VM with plenty of issues untile prep, which does not accept symbolic links and I therefore had to copy the fasta and bam mapping files and their indexed files
mkdir groopm
cd groopm

cp ../../../2_assembly/2_opera_assembly_out/contigs.fasta ./
cp ../../2_read_mapping/*.sorted.bam* ./

# generate the db.pm database parsing all necessary data
GROOPMDATA=$PWD/
# the --shm-size flag can be set to lower values than 100g in case of space lack... I used this to make sure that the container would not stop running in case of large files
docker run --rm -it --shm-size 100g -v $GROOPMDATA:/data ubuntu16.04/groopm:vs2 groopm parse -t 32 data/db.gm data/contigs.fasta data/MD_3_S11_L001_R_out.sorted.bam data/SRR7135606_out.sorted.bam data/SRR7135607_out.sorted.bam data/SRR7135608_out.sorted.bam data/SRR7135612_out.sorted.bam data/aln_nanopo_out.sorted.bam data/aln_pacbio_out.sorted.bam
# create the core bins
docker run --rm -it --shm-size 100g -v $GROOPMDATA:/data ubuntu16.04/groopm:vs2 groopm core data/db.gm
# I did not evenutally manually corrected the bins due to Tkinder problems with the image...
#docker run --rm -it --shm-size 100g -v $GROOPMDATA:/data ubuntu16.04/groopm:vs2 groopm refine data/db.gm
# recruit unbinned contigs
docker run --rm -it --shm-size 100g -v $GROOPMDATA:/data ubuntu16.04/groopm:vs2 groopm recruit data/db.gm
# extract the contig bins
docker run --rm -it --shm-size 100g -v $GROOPMDATA:/data ubuntu16.04/groopm:vs2 groopm extract data/db.gm data/contigs.fasta -p data/groopm

# get the membership to be used with dastool
for i in groopm_bin*.fna
do
	#myfile=`echo bin.${i}.fa`
	grep ">" ${i} > tmp_cntg_memb
	perl -pi -e "s/>//g" tmp_cntg_memb
	perl -pi -e "s/^/${i}\t/g" tmp_cntg_memb
	perl -pi -e "s/^^/\n/g" tmp_cntg_memb
	cat tmp_cntg_memb >> tmp_cntg_memb2
done
sed '/^[[:space:]]*$/d' tmp_cntg_memb2 | awk '{ print $2 "\t" $1}' | perl -pe "s/\.fna//g" > groopm_cntg_memb.txt
rm tmp_cntg_mem*

cd ..

# then go on with metabat2
mkdir metabat2
cd metabat2
runMetaBat.sh --saveCls --unbinned ../../2_read_mapping/contigs.fasta ../../2_read_mapping/SRR7135606_out.sorted.bam ../../2_read_mapping/SRR7135607_out.sorted.bam ../../2_read_mapping/SRR7135608_out.sorted.bam ../../2_read_mapping/SRR7135612_out.sorted.bam ../../2_read_mapping/aln_nanopo_out.sorted.bam ../../2_read_mapping/aln_pacbio_out.sorted.bam ../../2_read_mapping/MD_3_S11_L001_R_out.sorted.bam
# get the membership to be used with dastool
cd contigs.fasta.metabat-bins--unbinned-20210122_223128/
for i in $(seq 1 28)
do
	#myfile=`echo bin.${i}.fa`
	grep ">" bin.${i}.fa > tmp_cntg_memb
	perl -pi -e "s/>//g" tmp_cntg_memb
	perl -pi -e "s/^/metabat_bin${i}\t/g" tmp_cntg_memb
	perl -pi -e "s/^^/\n/g" tmp_cntg_memb
	cat tmp_cntg_memb >> tmp_cntg_memb2
done
sed '/^[[:space:]]*$/d' tmp_cntg_memb2 | awk '{ print $2 "\t" $1}' > ../metabat2_cntg_memb.txt
rm tmp_cntg_mem*
cd ..

cd ..


#### get the coverage requested by MaxBin2 and Solidbin
samtools faidx ../2_read_mapping/contigs.fasta
awk -v OFS='\t' {'print $1,$2'} ../2_read_mapping/contigs.fasta.fai > ../2_read_mapping/contigs.fasta.length.txt

mkdir z_cov_file_prep
cd z_cov_file_prep/

for i in `cat ../../list_all_bam`
do
	bedtools genomecov -ibam ../../2_read_mapping/${i}.sorted.bam > ${i}_cov.txt
done

mkdir cov_files

for i in *_cov.txt
do
   echo $i
   stub=${i%_cov.txt}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > cov_files/${stub}_cov.csv
done

# I used Collate.pl for combining the coverage files which is included in the Solidbin suit
perl ../../../z_accessory_files_dbs_and_executbles/Collate.pl cov_files > coverage.tsv

perl -pe "s/,/\t/g;" coverage.tsv > coverage_new.tsv

perl -pi -e "s/cov_files\///g" coverage_new.tsv

cd ..

## run the maxbin2 binner
mkdir maxbin2
cd maxbin2
cp cp -r ../z_cov_file_prep/cov_files ./
ls cov_files/ > my_cov_list
perl -pi -e "s/^/cov_files\//g" my_cov_list
cp ../../2_read_mapping/contigs.fasta ./
run_MaxBin.pl -contig contigs.fasta -out maxbinout -abund_list my_cov_list

# prepare the membership file to be used with dastool in the end
for i in maxbinout*.fasta
do
	myhead=${i%.*}
	grep ">" ${i} > tmp_cntg_memb
	perl -pi -e "s/>//g" tmp_cntg_memb
	perl -pi -e "s/^/${myhead}\t/g" tmp_cntg_memb
	perl -pi -e "s/^^/\n/g" tmp_cntg_memb
	cat tmp_cntg_memb >> tmp_cntg_memb2
done
sed '/^[[:space:]]*$/d' tmp_cntg_memb2 | awk '{ print $2 "\t" $1}' > maxbin2_cntg_memb.txt
rm tmp_cntg_mem*

cd ..

## and then the solidbin binner
mkdir solidbin
cd solidbin
mkdir input
mkdir output
cd input

cp ../../../2_read_mapping/contigs.fasta ./
# I used gen_kmer.py for calculateing 4-mers which I found in Solidbin
python3 ../../../../z_accessory_files_dbs_and_executbles/gen_kmer.py contigs.fasta 500 4
cp ../../z_cov_file_prep/coverage_new.tsv ./


cd ..

## run the solidbin
# copy some necessary files which, for some reason, are seeked in an auxiliary folder next to the input dir
cp -r /home/sotosv/software/binning/solidbin/SolidBin/auxiliary ./
# run solidbin using the installed version at software (for some reason this worked the previous time)
# at some point it requires checkm which can be found at the Solidbin miniconda environment (this solidbin version did not work with the rest software versions but it has checkm installed... alternatievly the right checkm installation should aslo work)
python ~/software/binning/solidbin/SolidBin/SolidBin.py  --contig_file input/contigs.fasta --coverage_profiles input/coverage_new.tsv --composition_profiles input/kmer_4_f500.csv --output output/solidbin_out.txt --log output/log.txt

cd output/ori_result
for i in ori_result*.bin
do
	myhead=${i%.*}
	grep ">" ${i} > tmp_cntg_memb
	perl -pi -e "s/>//g" tmp_cntg_memb
	perl -pi -e "s/^/${myhead}\t/g" tmp_cntg_memb
	perl -pi -e "s/^^/\n/g" tmp_cntg_memb
	perl -pi -e "s/^ori_result/solidbin/g" tmp_cntg_memb
	cat tmp_cntg_memb >> tmp_cntg_memb2
done
sed '/^[[:space:]]*$/d' tmp_cntg_memb2 | awk '{ print $2 "\t" $1}' > ../../solidbin_cntg_memb.txt
rm tmp_cntg_mem*
cd ../../



cd ..


### finally run the dastool refiner
mkdir dastool

cd dastool

cp ../../2_read_mapping/contigs.fasta ./
cp ../*/*_cntg_memb.txt ./
DAS_Tool -i maxbin2_cntg_memb.txt,metabat2_cntg_memb.txt,groopm_cntg_memb.txt,solidbin_cntg_memb.txt -c contigs.fasta -o dastoolout -l maxbin2,metabat,groopm,solidbin -t 32 --write_bins 1 --write_unbinned 1

# prepare an overal contig file with the bin annotation
cd dastoolout_DASTool_bins
for i in *.fa
do
	perl -pe "s/^>/>${i%.*}___/g" ${i} >> ../contigs_binned.fasta
done

# The loci of interest of the previous assembly version that previously showed high expression under the TBZ treatment were compared with this binning approach and all of them belonged in the groopm_bin_50 bin except for the two contig containing the two genes coding for AndAd and AndAa.
# andAc that previously had the highest expression of the locus was classified in this bin. I therefore moved the contig containing the other two "and" locus genes here
cd ../../../..

# obtain the contig membership information and place it in the coverage file
cd 3_binning/3_binning/dastool/dastoolout_DASTool_bins
grep "^>" *.fa | perl -pe "s/ .+//" | perl -pe "s/\.fa:>/\t/g" > ../contig_membership.txt
