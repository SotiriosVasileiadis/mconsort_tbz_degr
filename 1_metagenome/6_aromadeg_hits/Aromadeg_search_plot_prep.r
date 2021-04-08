# import and prepare the annotation, database and blast tables for merging into a master table
mat <- read.table("../4_annotation/ANNOTATION_DIR/genome.tsv", header = T, quote = "", check.names = F, sep = "\t", comment.char = "")
mat$binID <- gsub("___.+","",mat$seqid)
blastout <- read.table("blastouttbl.txt", header = T, quote = "", check.names = F, sep = "\t", comment.char = "")
blastout$ID <- gsub("\\|.+","",blastout$qaccver)
blastout$sacc <- gsub("\\|.+","",blastout$saccver)
aromaclass <- read.table("aroma_class.txt", header = T, quote = "", check.names = F, sep = "\t", comment.char = "")
bintaxID <- read.table("../z_accessory_files_dbs_and_executbles/binid_taxonomy.txt", header = T)

# merge the tables to a master table
matblast <- merge(mat, blastout, by.x = "ID", by.y = "ID", all = T)
matblastclass_pre <- merge(matblast, aromaclass, by.x = "sacc", by.y = "accNo", all.x = T)
matblastclass <- merge(matblastclass_pre, bintaxID, by.x = "binID", by.y = "binid", all = T)

# do some housekeeping of the new matrix
mataromafoc <- matblastclass[complete.cases(matblastclass$saccver),]
mataromafoc <- mataromafoc[order(mataromafoc$BinTaxID,mataromafoc$supergrp,mataromafoc$grp),]
my_suprgrps <- levels(mataromafoc$supergrp)

# prepare the matrix for the bubbleplot
mataromafoc$genecnts <- rep(1,nrow(mataromafoc))
mataromafoc$AromaDeg <- paste(mataromafoc$supergrp,mataromafoc$grp,mataromafoc$descr, sep = "---")

# prepare the plot
library(ggplot2)
g1 <- ggplot(subset(mataromafoc), aes(x=BinTaxID,y=AromaDeg)) +
  geom_point(aes(size=genecnts)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "bin", y = paste("Superfamily --- group --- best hit",sep = "")) +
  ggtitle(paste("Metagenome"))
g1$labels$size <- "gene counts"


cairo_pdf(paste("bubble_plot_COGs_TBZ_vs_SUC_gene_trans.pdf", sep = ""), width = 10, height = 30)
g1
dev.off()
