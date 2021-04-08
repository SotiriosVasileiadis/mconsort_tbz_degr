# read and prep the necessary matrices
cov_mat <- read.table("coverage_new.tsv", sep = "\t", header = T)
bin_table <- read.table("../dastool/contig_membership.txt")
colnames(bin_table) <- c("binID","contigID")

# merge the tables
merger <- merge(cov_mat, bin_table, by.x = "contig", by.y = "contigID")
merger_fin <- merger[,grep("SRR|contig|binID",colnames(merger))]

# calculate the mean coverage per bin for each of the short read SRR datasets
aggregated <- aggregate(merger_fin[,grep("SRR",colnames(merger_fin))], by = list(merger_fin$binID), FUN = mean)
aggregated_fin <- data.frame(aggregated[,grep("SRR",colnames(aggregated))],row.names = aggregated$Group.1)

# obtain the relative abundances of the bins and save the table
library(vegan)
aggregated_RA <- round(100*decostand(aggregated_fin, MARGIN = 2, method = "total"),3)
write.table(aggregated_RA, file = "bin_RA.txt", sep = "\t", quote = F, col.names = NA)
