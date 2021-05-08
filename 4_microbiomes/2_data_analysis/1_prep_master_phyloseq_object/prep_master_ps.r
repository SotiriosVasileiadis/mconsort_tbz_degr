#### package installation ----
# install the necessary packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2")

#### package loading ----
# load the package and check the package version
library(dada2); packageVersion("dada2")

mythreads <- 32

# copy the demultiplexed sequences in a single folder -

system("cp ../../1_fetch_and_demultiplex/*/2_demultiplexed_sequences/demux_out/analysis_ready/*.fastq ./1_demultiplexed_files/")

# input file path variable
path <- "1_demultiplexed_files/" 
# list files to verify
list.files(path)
# set the variables containing all the forward and the reverse paths to the files of interest with the list.files command
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- gsub("_lib_.+","",basename(fnFs))

# Plot the per-base qualities
pdf(file = "seq_quality_scores.pdf", height = 4, width = 8, onefile = T)
for(i in 1:length(fnFs)){
  print(i)
  plot <- plotQualityProfile(c(fnFs[i],fnRs[i]))
  print(plot)
}
dev.off()

#### sequence quality filtering and control, error modelling, and dereplication ----
# set the file paths where the quality controlled sequences are going to be saved
filtFs <- file.path("filtered", paste(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste(sample.names, "_R_filt.fastq.gz"))


# filter the sequences and save then in the folders provided above and get their statistics in a table
# use matchIDs=F because otherwise there will be issues with the ID replacements provided by SRA
# also truncate the sequences to something like 220 bp which is compatible with the V4 amplicon length and also prevents issues with the MiSeq 2x300bp reads mixing with the HiSeq 2x250bp
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft = 11,
                     compress=TRUE, multithread=mythreads, matchIDs=F, truncLen=c(220,220), minLen = 150)
head(out)

# learns error rates using a machine learning algorithm
errF <- learnErrors(filtFs, multithread=mythreads)
plotErrors(errF, nominalQ=TRUE)
errR <- learnErrors(filtRs, multithread=mythreads)
plotErrors(errR, nominalQ=TRUE)


# dereplication of each one of the red pairs to unique sequences (collapsing of the identical sequences for each pair per sample)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# sample composition inference after read correction
dadaFs <- dada(derepFs, err=errF, multithread=mythreads)
dadaRs <- dada(derepRs, err=errR, multithread=mythreads)


# merge read pairs retaining the per sequence sample information
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#### construct the sequence table, remove the chimeras, and create a summary ----
# construct sequence table
seqtab <- makeSequenceTable(mergers)

# chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mythreads, verbose=TRUE)
dim(seqtab.nochim)

# record the portion of good sequences out of the total prior the chimera removal
sum(seqtab.nochim)/sum(seqtab)

# track reads through the pipeline
getN <- function(x) {sum(getUniques(x))}
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# save the read quality control data (if you want to download the table at your computer, you need to go to the "More" option and select export)
write.table(track, file = "readQC.txt", col.names = NA, sep = "\t", quote = FALSE)

#### taxonomically classify the sequences ----
# obtain the taxonomy information for each sequence at all possible levels
taxa <- assignTaxonomy(seqtab.nochim, minBoot = 80, "../../dbs/silva_nr99_v138.1_train_set.fa.gz", multithread=mythreads)



library(phangorn)
library(DECIPHER)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)

# prepare the Neighbor joining tree for the phyloseq object
treeNJ <- NJ(dm) # Note, tip order != sequence order



# load the phyloseq package which is quite useful for phylogenetic marker diversity studies
library(phyloseq); packageVersion("phyloseq")

## load the experimental design file
samdf <- read.table("design",header = T, row.names = 1, sep = "\t")

#### construct the main phyloseq object ----
ps_orig <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE) 
               ,sample_data(samdf) 
               ,tax_table(taxa)
               ,phy_tree(treeNJ)
               )
# replace th sequences as taxon names with something shorter and more meaninful 
# first install Biostrings which is required for saving the sequences as fasta
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
library("Biostrings")
sequences <- Biostrings::DNAStringSet(taxa_names(ps_orig))
names(sequences) <- taxa_names(ps_orig)
ps <- merge_phyloseq(ps_orig, sequences)

# replace the taxon names (the sequences of each ASV with something easier to read)
library(stringr)
taxa_names(ps) <- paste("ASV",str_pad(1:length(taxa_names(ps)),3, pad = "0"),sep = "")


## remove undezired taxa (unknown, Eukaryotes, mitochondria and chloroplasts)
# Show available ranks in the dataset
rank_names(ps)

# Create table, number of features for each phyla
table(tax_table(ps)[, "Kingdom"], exclude = NULL)

# The following ensures that features with ambiguous phylum annotation are removed (this is optional in case you have good domain (Kingdom referred here) confidence. Here I also remove non target taxa. Further down there is also subsetting and removal of low in abundance phyla.
# remove the uncharacterized or empty taxa
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized") & !Kingdom %in% c("Eukaryota","Archaea") & !is.na(Kingdom))
# use the following command to check which unique values still exist (change Family to Order etc)
#unique(tax_table(ps0)[,"Family"])[order(unique(tax_table(ps0)[,"Family"]))]
# remove Chloroplasts
ps0 <- subset_taxa(ps0, !Order %in% c("Chloroplast"))
# remove Mitochondria
ps0 <- subset_taxa(ps0, !Family %in% c("Mitochondria"))

# select taxa by the deliverable metadata factor 

## prepare the taxonomy that gives at best genus level
ps0_best_tax <- ps0
for(i in 2:ncol(tax_table(ps0_best_tax))){
  for(j in 1:nrow(tax_table(ps0_best_tax))){
    if(is.na(tax_table(ps0_best_tax)[j,i])){
      tax_table(ps0_best_tax)[j,i] <- tax_table(ps0_best_tax)[j,i-1]
    } 
  }
}


# add one more annotation column parted by the ASV and genus tags
mytxtbl_ammended <- data.frame(tax_table(ps0_best_tax), stringsAsFactors = F)
mytxtbl_ammended$ASV_genus <- paste(row.names(mytxtbl_ammended), mytxtbl_ammended$Genus)

# and rename the most dominant genus level taxa for consistency with the MAG annotations
repl_names <- c("Novosphingobium","Sinobacteraceae","Novosphingobium","Sphingomonas","Bradyrhizobium","Sinobacteraceae","Sphingomonas","Bradyrhizobium","Hydrogenophaga","Hydrogenophaga","Filimonas","Thiobacillus","Filimonas","Thiobacillus","Hyphomicrobium","Flavobacterium","Sphingopyxis","Shinella","Hyphomicrobium","Sphingopyxis","Pedobacter","Microbacterium","Shinella","Hydrogenophaga","Flavobacterium","Sphingopyxis","Pedobacter","Hydrogenophaga","Hydrogenophaga","Microbacterium")

mytxtbl_ammended_fin <- mytxtbl_ammended
for(i in 1:length(repl_names)){
  mytxtbl_ammended_fin$Genus[i] <- repl_names[i]
  mytxtbl_ammended_fin$ASV_genus[i] <- paste(row.names(mytxtbl_ammended[i,]),repl_names[i])
}

mytxtbltmp <- tax_table(mytxtbl_ammended_fin)
row.names(mytxtbltmp) <- row.names(mytxtbl_ammended_fin)
colnames(mytxtbltmp) <- colnames(mytxtbl_ammended_fin)
ps0_best_tax2 <- ps0_best_tax
tax_table(ps0_best_tax2) <- mytxtbltmp
ps0_best_tax3 <- ps0_best_tax2
taxa_names(ps0_best_tax3) <- as.data.frame(tax_table(ps0_best_tax2), stringsAsFactors = F)$ASV_genus
#### write/read the final master phyloseq object ----
## write the phyloseq object on the disk with the following command
saveRDS(ps0_best_tax3, "my_master_ps.rds")

# prepare the percent participation phyloseq object
ps0_best_tax3_ra  = transform_sample_counts(ps0_best_tax3, function(x) 100 * x / sum(x))


cairo_pdf("barplot_top50.pdf", height = 8, width = 14)
par(mar = c(17,4,4,4))
barplot(colMeans(data.frame(check.names = F, otu_table(ps0_best_tax3_ra)))[1:50], las = 2, ylab = "RA (%)")
dev.off()
