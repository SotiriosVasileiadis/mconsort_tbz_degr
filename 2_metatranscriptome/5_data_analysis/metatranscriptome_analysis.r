#### load and prepare the data ----
## prep the design objects
# set the working directory to the source file location
design_tab <- read.table("input_files/design", header = T, sep = "\t", stringsAsFactors = F)
design_tab$time <- as.factor(design_tab$time)
design_tab_sel_pre <- design_tab
myint <- factor(interaction(design_tab_sel_pre[,2:3]), levels = c("TBZ.57","TBZ.73","TBZ.109","SUC.57","SUC.73","SUC.109"))
design_tab_sel_pre$fileName <- paste(design_tab_sel_pre$name,".txt",sep = "")
design_tab_sel_pre$sampleName <- paste(design_tab_sel_pre$treatment, design_tab_sel_pre$time,"_",rep(c("a","b","c"),6), sep = "")
design_tab_sel_pre$treatment_time <- factor(paste(design_tab_sel_pre$treatment, "_", design_tab_sel_pre$time, sep = ""), levels = c("TBZ_57","TBZ_73","TBZ_109","SUC_57","SUC_73","SUC_109"))
design_tab_sel_pre <- design_tab_sel_pre[,c(7:6,1:5)]
design_tab_sel_pre$treatment <- factor(design_tab_sel_pre$treatment,levels = c("TBZ","SUC"))
design_tab_sel_pre$fileName <- paste("counts_",design_tab_sel_pre$fileName, sep = "")
design_tab_sel <- design_tab_sel_pre


## load the table with all gene annotations and sequences
# obtain the dfast generated annotation (it was generated out of the gff file with the script "prep_genome_tsv.R and was manually curated due to the newline generation because of special characters... this could also be mended through the R script but was eventually corrected through a text editor for better control of potential issues)
my_annot_tbl_pre <- read.table("../../1_metagenome/4_annotation/ANNOTATION_DIR/genome.tsv", header = T, quote = "", comment.char = "", sep = "\t", stringsAsFactors = F)
# replace the useless X column with the row index number
my_annot_tbl_pre$X <- row.names(my_annot_tbl_pre)
colnames(my_annot_tbl_pre)[which(colnames(my_annot_tbl_pre) == "X")] <- "idx"
# then, the SEED annotation (prepare it by adding a header and also by indexing for merging)
my_seed_annot <- read.table("../../1_metagenome/4_annotation/SEED_annotation/myresults.m8_wth_ids.converted", header = F, quote = "", comment.char = "", sep = "\t", stringsAsFactors = F, fill = T)
colnames(my_seed_annot) <- c("qseqid","qseqid.1","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","SEED4","SEED3","SEED2","SEED1")
my_seed_annot$idxseed <- gsub("\\|.+","",my_seed_annot$qseqid)

# then, the eggnog annotation (required some manual curation as well: the header had to be uncommented, while some field names had to be added)
my_eggnog_annot <- read.table("../../1_metagenome/4_annotation/eggnog_mapper_annotation/query_seqs.fa.emapper.annotations_franalysis", header = T, quote = "", comment.char = "#", sep = "\t", stringsAsFactors = F)
# prepare it for merging as well
my_eggnog_annot$idxeggnog <- gsub("\\|.+","",my_eggnog_annot$query_name)

# finally, the plasflow annotation
my_plasmid_annot <- read.table("../../1_metagenome/4_annotation/plasflow/filtered_genome_plas_pred.txt", header = T, quote = "", comment.char = "#", sep = "\t", stringsAsFactors = F)

# prep the final annotation table
my_annot_tbl_pre1 <- merge(my_annot_tbl_pre, my_plasmid_annot, by.x = "seqid", by.y = "contig_name", all = T)
my_annot_tbl_pre2 <- merge(my_annot_tbl_pre1, my_seed_annot, by.x = "ID", by.y = "idxseed", all = T)
my_annot_tbl_pre3 <- merge(my_annot_tbl_pre2, my_eggnog_annot, by.x = "ID", by.y = "idxeggnog", all = T)

# save the new
my_annot_tbl <- my_annot_tbl_pre3

# add the bin information and the row.names
my_annot_tbl$binid <- gsub("___.+","",my_annot_tbl$seqid)
row.names(my_annot_tbl) <- my_annot_tbl$ID
# create also the gc content column
library(stringr)
GCcounts <- str_count(my_annot_tbl$myfna, "g|c")
Allcounts <- str_count(my_annot_tbl$myfna, "a|t|g|c")
my_annot_tbl$gc_content <- GCcounts/Allcounts

# back the final table up and save it in the annotation dir
my_annot_tbl_bckp <- my_annot_tbl
write.table(my_annot_tbl, file = "../../1_metagenome/4_annotation/annotation_table_master_dfast_plasflow_seed_eggnog.tsv", quote = F, sep = "\t", col.names = NA)

## load the htseq files
# set the directory of the count files
directory <- "../4_htseq_counts/strandedR"
library("DESeq2")
# load the files (beware that the sampleName column must be first and followed by the fileName column)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = design_tab_sel[,c(2,7,4,1,3,5,6)],
                                       directory = directory,
                                       design= ~ treatment)
# incorporate the transcript count information at the annotation table
my_cnt_annot_tab <- merge(my_annot_tbl, ddsHTSeq@assays@data@listData$counts, by = "row.names", all.x = T)
my_cnt_annot_tab$transcripts <- rep(1,nrow(my_cnt_annot_tab))
# add the gene length information
my_cnt_annot_tab$gene_length <- my_cnt_annot_tab$end - my_cnt_annot_tab$start

# save the file containing the annotations and transcript counts
write.table(my_cnt_annot_tab, file = "annotation_table_master_dfast_plasflow_seed_eggnog_extended.tsv", quote = F, sep = "\t", col.names = NA)


# prepare the final design matrix based on the interactions of the grouping variables of treatment and time
design.grps <- interaction(design_tab_sel$treatment, design_tab_sel$time)
design.grps <- factor(design.grps, levels = c("TBZ.57", "SUC.57", "TBZ.73", "SUC.73", "TBZ.109", "SUC.109"))





#### run edgeR ----
library(edgeR)
## prep the dgelist object ----
library(DEFormats)
my_dgelist <- as.DGEList(ddsHTSeq)


## As a rule of thumb, we require that a gene have a count of at least 10–15 in at least some libraries before it is considered to be expressed in the study. We could explicitly select for genes that have at least a couple of counts of 10 or more, but it is slightly better to base the filtering on count-per-million (CPM) values so as to avoid favoring genes that are expressed in larger libraries over those expressed in smaller libraries. For the current analysis, we keep genes that have CPM values above 4 (=> 16 copies per library of about 4 million reads) in at least two libraries

my_dgelist_good <- rowSums(cpm(my_dgelist) > 4) >= 2

# is close to but better than
#keep <- rowSums(my_dgelist$counts) > 50

## Here the cutoff of 0.5 for the CPM has been chosen because it is roughly equal to 10/L where L is the minimum library size in millions. The library sizes here are 20–25 million. We used a round value of 0.5 just for simplicity; the exact value is not important because the downstream differential expression analysis is not sensitive to the small changes in this parameter. The requirement of ≥2 libraries is because each group contains two replicates. This ensures that a gene will be retained if it is expressed in both the libraries belonging to any of the six groups.
## The above filtering rule attempts to keep the maximum number of interesting genes in the analysis, but other sensible filtering criteria are also possible. For example keep <- rowSums(y\$counts) > 50 is a very simple criterion that would keep genes with a total read count of more than 50. This would give similar downstream results for this dataset to the filtering actually used. Whatever the filtering rule, it should be independent of the information in the targets file. It should not make any reference to which RNA libraries belong to which group, because doing so would bias the subsequent differential expression analysis.


my_dgelist_fin <- my_dgelist[my_dgelist_good, ]

my_dgelist_fin <- calcNormFactors(my_dgelist_fin)

## The normalization factors of all the libraries multiply to unity. A normalization factor below one indicates that a small number of high count genes are monopolizing the sequencing, causing the counts for other genes to be lower than would be usual given the library size. As a result, the effective library size will be scaled down for that sample. Here we see that the luminal-lactating samples have low normalization factors. This is a sign that these samples contain a number of very highly upregulated genes.
## Note. In general, we find TMM normalization to be satisfactory for almost all well-designed mRNA gene expression experiments. Single-cell RNA-seq is an exception, for which specialized normalization methods are needed (Lun, Bach, and Marioni 2016). Another, less common, type of study requiring special treatment is that with global differential expression, with more than half of the genome differentially expressed between experimental conditions in the same direction (D. Wu et al. 2013). Global differential expression should generally be avoided in well designed experiments. When it can’t be avoided, then some normalization reference such as spike-ins needs to be built into the experiment for reliable normalization to be done (Risso et al. 2014).

## MDS ----
pch <- c(21,22,23,21,22,23) #16,17,18)
colors <- rep(c("grey20", "grey20","grey20","grey70", "grey70","grey70"), 2)

library(vegan) # for the ordiellipse function
cairo_pdf("0_lfc_based_all_samples_MDS.pdf")
myMDSplot <- plotMDS(my_dgelist_fin, bg=colors[1:6], col = "black", pch=pch[1:6], bty = "n", ylim = c(-5,5), xlim = c(-6,5))
ordiellipse(myMDSplot$cmdscale.out, groups = design_tab_sel$treatment, kind = "ehull", lty = 2)
legend("topleft", legend=levels(design_tab_sel$treatment_time), pch=pch, pt.bg=colors, col = "black", ncol=1, bty = "n")
dev.off()



#### canonical analysis ----
# prepare the plot of the canonical analysis of the normalized data in order to assess the importance of the treatment and time

# prep the vegan compatible plot
my_CanNorm_mat <- t(cpm(my_dgelist_fin))

# select the appropriate method according to the detrended correspondence analysis
dcares<-decorana(decostand(my_CanNorm_mat,"hellinger"))
dcares
# Call:
#   decorana(veg = decostand(my_CanNorm_mat, "hellinger")) 
# 
# Detrended correspondence analysis with 26 segments.
# Rescaling of axes with 4 iterations.
# 
# DCA1    DCA2    DCA3     DCA4
# Eigenvalues     0.2139 0.08059 0.04954 0.061087
# Decorana values 0.2159 0.06328 0.01875 0.008098
# Axis lengths    2.2145 1.05910 0.96369 0.892494
# results suggest linear gradients and therefore rda performance

myrda<-rda(decostand(my_CanNorm_mat,"hellinger") ~ treatment + time, design_tab_sel)

anova=anova(myrda)
capture.output(anova, file="0_RDA_anova.txt")
anova_axis=anova(myrda, by="axis",perm.max=9999)
capture.output(anova_axis, file="0_RDA_anova_axis.txt")
anova_terms=anova(myrda, by="terms",permu=9999)
capture.output(anova_terms, file="0_RDA_anova_terms.txt")
k=anova
k$Pr[1]

# get the explained variance

tot_expl<-round(100*sum(anova_terms$Variance[1:2])/sum(anova_terms$Variance),1)
treat_expl<-round(100*sum(anova_terms$Variance[1])/sum(anova_terms$Variance[1:2]),1)
time_expl<-round(100*sum(anova_terms$Variance[2])/sum(anova_terms$Variance[1:2]),1)

PERMANOVA <- adonis2(decostand(my_CanNorm_mat,"hellinger") ~ treatment + time, data = design_tab_sel, permutations = 9999)
capture.output(PERMANOVA, file = "0_RDA_PERMANOVA.txt")

RDA1=100*myrda$CCA$eig[1]/sum(myrda$CCA$eig);
RDA2=100*myrda$CCA$eig[2]/sum(myrda$CCA$eig);
RDA1=format(round(RDA1,1), nsmall = 1);
RDA2=format(round(RDA2,1), nsmall = 1);

sc.1 <- scores(myrda,display="sites",scaling=1)
df<-data.frame(sc.1,design_tab_sel[,4:5])
df$time <- as.numeric(as.character(df$time))

library(ggplot2)
library(PerformanceAnalytics)
library(RColorBrewer)
rbriewer_set<-brewer.pal(9, "Paired")

cairo_pdf(file="0_RDA.pdf",height=6,width=8)
theme_set(theme_bw())
p<-ggplot (df, aes (x = RDA1, y = RDA2))
print(p   + 
        labs(title=paste("PERMANOVA: comm. struct. ~ treat + time   (",round(100*sum(PERMANOVA$R2[1:2]), 1),"% of total var.,P ",round(max(PERMANOVA$`Pr(>F)`, na.rm = T),3),") \n part. of explained variance: treat/time = ",round(100*(PERMANOVA$R2[1]/sum(PERMANOVA$R2[1:2])), 1),"%, P ",round(PERMANOVA$`Pr(>F)`[1],4)," / ",round(100*(PERMANOVA$R2[2]/sum(PERMANOVA$R2[1:2])), 1),"%, P ",round(PERMANOVA$`Pr(>F)`[2],3),sep=""), x=paste("RDA1 (",RDA1," %)",sep=""), y=paste("RDA2 (",RDA2," %)",sep="")) + 
        geom_point(aes(fill = treatment, shape = treatment, size = time)) + 
        scale_size_continuous("time",breaks=c(57,73,109)) + 
        # scale_colour_manual ("treatment", values = brewer.pal(9,"Paired")) + 
        scale_fill_manual ("treatment", values = c("red","steelblue3")) + 
        scale_colour_manual (values = "black") + 
        stat_ellipse(aes(fill=as.factor(df$treatment)),type = "t", linetype = 2, colour="grey40") +
        scale_shape_manual ("treatment",values= c(21:22)))
dev.off()




### fold change plot of each sample vs the rest average ----
cairo_pdf("0_lfc_based_all_samples_MD_vs_average.pdf", height = 14, width = 8)
par(mfrow = c(6,3))
for(i in 1:ncol(my_dgelist_fin$counts)){
  plotMD(my_dgelist_fin, column=i, bty = "n", main = colnames(my_dgelist_fin$counts)[i])
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()


treatment_time <- my_dgelist_fin$samples$treatment_time
## prepare the design/model matrix ----
mydesign <- model.matrix(~0+treatment_time)



## estimate the dispersion and rest gene expression variance parameters ----
my_dgelist_fin <- estimateDisp(my_dgelist_fin, mydesign, robust=TRUE)


cairo_pdf("0_all_samples_BCV.pdf")
plotBCV(my_dgelist_fin, bty = "n")
dev.off()

# quasi likelihood accounting for gene specific variability
fit <- glmQLFit(my_dgelist_fin, mydesign, robust=TRUE)


cairo_pdf("0_all_samples_QLDisp.pdf")
plotQLDisp(fit, bty = "n")
dev.off()


#### Test for differential expression in a pairwise manner ----
# create all combinations of possible tests according to the mydesign values
mycombs <- t(combn(colnames(mydesign),2))
mycombs_names <- paste(mycombs[,1],mycombs[,2], sep = "_")


# prep tha matrix to be amended with the rnaseq comparison results
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]

# add also the cpm for each gene
my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,cpm(my_dgelist_fin), by = "row.names", all.x =T)
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]

# prep also the interaction/sample combination
design.grps_descr <- data.frame(samples = design_tab_sel$sample, grps = design.grps, stringsAsFactors = F)

# run the analysis
for(i in 1:length(mycombs_names)){
  contr0 <- makeContrasts(contrasts = c(mycombs[i,2],mycombs[i,1]), levels=mydesign)
  mycomb <- data.frame(Contrasts = contr0[,1] - contr0[,2])
  
  
  # We use QL F-tests instead of the more usual likelihood ratio tests (LRT) as they give stricter error rate control by accounting for the uncertainty in dispersion estimation: 
  res <- glmQLFTest(fit, contrast=mycomb)
  # get the fold change values etc and add them to the main table
  my_stats_tbl <- res$table
  my_stats_tbl$q.value <- p.adjust(my_stats_tbl$PValue,"fdr")
  colnames(my_stats_tbl) <- paste(gsub("design.grps","",paste(mycombs[i,1],"_vs_",mycombs[i,2], sep  ="")), colnames(my_stats_tbl), sep = "_")
  
  
  mynumerator <- gsub("_",".",gsub("treatment_time","",mycombs[i,1]))
  mydenominator <- gsub("_",".",gsub("treatment_time","",mycombs[i,2]))
  mynumersamples <- design.grps_descr[which(design.grps_descr$grps==mynumerator),1]
  mydenomsamples <- design.grps_descr[which(design.grps_descr$grps==mydenominator ),1]
  mysamples <- design.grps_descr[which(design.grps_descr$grps==mynumerator | design.grps_descr$grps==mydenominator ),1]
  
  # get the log2FC table and add the comparisons in the header
  mylogFC_tbl <- log2((cpm(my_dgelist_fin)[,mysamples] + 0.01)/rowMeans(cpm(my_dgelist_fin)[,mydenomsamples] + 0.01))
  colnames(mylogFC_tbl) <- paste(gsub("design.grps","",paste(mycombs[i,1],"_vs_",mycombs[i,2],"_log2FC", sep  ="")), colnames(mylogFC_tbl), sep = "_")
  
  # merge the tables (mask this step in case you need only to plot things)
  my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,mylogFC_tbl, by = "row.names", all =T)
  row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
  my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]
  my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,my_stats_tbl, by = "row.names", all =T)
  row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
  my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]
  
  
  # prepare the rought plots for the DE genes
  is.de <- decideTestsDGE(res, p.value = 0.01 , lfc = 1)
  
  cairo_pdf(paste("0_pairwise_DE_ME_",gsub("treatment_time","",paste(mycombs[i,2]," vs ",mycombs[i,1], sep  ="")),"_q0.01_lfc1.pdf", sep = ""), height = 5, width = 7)
  plotMD(res, status=is.de, values=c(1,-1), col=c("red","darkblue"), legend=F, bty = "n", hl.cex = 0.3, pt.cex = 1, main = gsub("treatment_time","",paste(mycombs[i,2]," vs ",mycombs[i,1], sep  ="")), ylab = "log2FC", xlab = "mean log CPM", ylim = c(-15,15))
  abline(h = 0, lty = 2, col = "purple")
  text(12, -10, labels = paste(gsub("treatment_time","",mycombs[i,1]), " (", summary(is.de)[1]," genes)", sep = ""), col = "blue", cex = 1.2)
  text(12, 0, labels = paste("(", summary(is.de)[2]," genes)", sep = ""), col = "black", cex = 1.2)
  text(12, 10, labels = paste(gsub("treatment_time","",mycombs[i,2]), " (", summary(is.de)[3]," genes)", sep = ""), col = "red", cex = 1.2)
  dev.off()
}



# save the table with the sample/date contrasts

write.table(my_cnt_annot_tab_wth_sts, file = "00_my_final_tab_with_seqs_and_sampletime_contrasts.tsv", sep = "\t", quote = F, col.names = NA)




















#### Time series analysis for the TBZ ----

mydesigntime <- mydesign
# just change the colnames in the myde
colnames(mydesigntime) <- gsub("_","\\.",gsub("treatment_time","",colnames(mydesigntime)))

con <- makeContrasts(
  TBZ.57.73 = TBZ.57 - TBZ.73,
  TBZ.57.109 = TBZ.57 - TBZ.109,
  TBZ.73.109 = TBZ.73 - TBZ.109, levels=mydesigntime)


res <- glmQLFTest(fit, contrast=con)
# get the fold change values etc and add them to the main table
my_stats_tbl <- res$table
my_stats_tbl$q.value <- p.adjust(my_stats_tbl$PValue,"fdr")


#prepare for the log2FC matrix use the initial samples (day 57) for comparins 
mydenomsamples <- design.grps_descr[which(design.grps_descr$grps== "TBZ.57" ),1]
mysamples <- design.grps_descr[grep("TBZ",design.grps_descr$grps),1]

# get the log2FC table and add the comparisons in the header
mylogFC_tbl <- log2((cpm(my_dgelist_fin)[,mysamples] + 0.01)/rowMeans(cpm(my_dgelist_fin)[,mydenomsamples] + 0.01))
colnames(mylogFC_tbl) <- paste(colnames(mylogFC_tbl), "_vs_TBZ.57_mean_log2FC", sep = "")

## prepare a heatmap ----
my_annotab_for_heatmap <- my_cnt_annot_tab_wth_sts[,grep("SEED",colnames(my_cnt_annot_tab_wth_sts))]
my_annotab_for_heatmap2 <- as.data.frame(paste(my_annotab_for_heatmap$SEED1,my_annotab_for_heatmap$SEED2,my_annotab_for_heatmap$SEED3,my_annotab_for_heatmap$SEED4,sep = "  /  "))
row.names(my_annotab_for_heatmap2) <- row.names(my_annotab_for_heatmap)
heattab <- merge(mylogFC_tbl,my_stats_tbl, by = "row.names")
row.names(heattab) <- heattab$Row.names
heattab2 <- merge(heattab, my_annotab_for_heatmap2, by = "row.names")
heattab_names <- paste(heattab2[,1],heattab2[,ncol(heattab2)])
heattab3 <- heattab2[,3:(ncol(heattab2)-1)]
row.names(heattab3) <- heattab_names
heattab3_ord <- heattab3[order(heattab3$q.value),]

# select the highly DE and expressed genes with cutoffs: q.value of 0.001; abs log2fold change 2; logCPM 4
heattab3_ord_sel <- heattab3_ord[heattab3_ord$q.value < 0.001 & abs(heattab3_ord$logFC.TBZ.57.73) > 2,]
# heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logFC.TBZ.57.109) > 2,]
heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logCPM) > 4,]
heattab3_ord_sel_fin <- heattab3_ord_sel[,grep("_vs_",colnames(heattab3_ord_sel))]


# prepare the heatmap
library(pheatmap)
# set the breaks and colours
mybreaks <- seq(-max(abs(as.matrix(heattab3_ord_sel_fin))),max(abs(as.matrix(heattab3_ord_sel_fin))), length.out = 120)
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[6:1])(n = 119) # the nature colors
# open the graphics device, plot and close the device
cairo_pdf(file = paste("0_Heatmap_TBZ_time_q0.001_lFC2_logCPM2.pdf", sep = ""),height = 100, width = 25)
pheatmap(as.matrix(heattab3_ord_sel_fin), color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
                   #, cluster_rows = F, cluster_cols = F
                   , treeheight_row = 60
                   , treeheight_col = 60
)
dev.off()


# merge the tables
my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,mylogFC_tbl, by = "row.names", all =T)
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]
colnames(my_stats_tbl) <- paste("TBZsamples_vs_TBZ.57_",colnames(my_stats_tbl),sep = "")
my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,my_stats_tbl, by = "row.names", all =T)
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]








#### Time series analysis for the SUC ----

mydesigntime <- mydesign
# just change the colnames in the myde
colnames(mydesigntime) <- gsub("_","\\.",gsub("treatment_time","",colnames(mydesigntime)))

con <- makeContrasts(
  SUC.57.73 = SUC.57 - SUC.73,
  SUC.57.109 = SUC.57 - SUC.109,
  SUC.73.109 = SUC.73 - SUC.109, levels=mydesigntime)


res <- glmQLFTest(fit, contrast=con)
# get the fold change values etc and add them to the main table
my_stats_tbl <- res$table
my_stats_tbl$q.value <- p.adjust(my_stats_tbl$PValue,"fdr")


#prepare for the log2FC matrix use the initial samples (day 57) for comparisons 
mydenomsamples <- design.grps_descr[which(design.grps_descr$grps== "SUC.57" ),1]
mysamples <- design.grps_descr[grep("SUC",design.grps_descr$grps),1]

# get the log2FC table and add the comparisons in the header
mylogFC_tbl <- log2((cpm(my_dgelist_fin)[,mysamples] + 0.01)/rowMeans(cpm(my_dgelist_fin)[,mydenomsamples] + 0.01))
colnames(mylogFC_tbl) <- paste(colnames(mylogFC_tbl), "_vs_SUC.57_mean_log2FC", sep = "")

## prepare a heatmap ----
my_annotab_for_heatmap <- my_cnt_annot_tab_wth_sts[,grep("SEED",colnames(my_cnt_annot_tab_wth_sts))]
my_annotab_for_heatmap2 <- as.data.frame(paste(my_annotab_for_heatmap$SEED1,my_annotab_for_heatmap$SEED2,my_annotab_for_heatmap$SEED3,my_annotab_for_heatmap$SEED4,sep = "  /  "))
row.names(my_annotab_for_heatmap2) <- row.names(my_annotab_for_heatmap)
heattab <- merge(mylogFC_tbl,my_stats_tbl, by = "row.names")
row.names(heattab) <- heattab$Row.names
heattab2 <- merge(heattab, my_annotab_for_heatmap2, by = "row.names")
heattab_names <- paste(heattab2[,1],heattab2[,ncol(heattab2)])
heattab3 <- heattab2[,3:(ncol(heattab2)-1)]
row.names(heattab3) <- heattab_names
heattab3_ord <- heattab3[order(heattab3$q.value),]

# select the highly DE and expressed genes with cutoffs: q.value of 0.001; abs log2fold change 2; logCPM 4
heattab3_ord_sel <- heattab3_ord[heattab3_ord$q.value < 0.001 & abs(heattab3_ord$logFC.SUC.57.73) > 2,]
# heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logFC.SUC.57.109) > 2,]
heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logCPM) > 4,]
heattab3_ord_sel_fin <- heattab3_ord_sel[,grep("_vs_",colnames(heattab3_ord_sel))]



# prepare the heatmap
library(pheatmap)
# set the breaks and colours
mybreaks <- seq(-max(abs(as.matrix(heattab3_ord_sel_fin))),max(abs(as.matrix(heattab3_ord_sel_fin))), length.out = 120)
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[6:1])(n = 119) # the nature colors
# open the graphics device, plot and close the device
cairo_pdf(file = paste("0_Heatmap_SUC_time_q0.001_lFC2_logCPM4.pdf", sep = ""),height = 100, width = 25)
pheatmap(as.matrix(heattab3_ord_sel_fin), color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
         #, cluster_rows = F, cluster_cols = F
         , treeheight_row = 60
         , treeheight_col = 60
)
dev.off()


# merge the tables
my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,mylogFC_tbl, by = "row.names", all =T)
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]
colnames(my_stats_tbl) <- paste("SUCsamples_vs_SUC.57_",colnames(my_stats_tbl),sep = "")
my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,my_stats_tbl, by = "row.names", all =T)
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]





# prepare the MD plots for the DE genes
contr0 <- makeContrasts(contrasts = c("SUC.57","SUC.109"), levels=mydesigntime)
con <- data.frame(Contrasts = contr0[,1] - contr0[,2])


res <- glmQLFTest(fit, contrast=con)

cairo_pdf(paste("0_SUC_57vs109samples_DE_MA_q.0.01_lfc1.pdf", sep = ""), height = 5, width = 6)
is.de <- decideTestsDGE(res, p.value = 0.01 , lfc = 1)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend=NULL, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "SUC.109 vs SUC.57", ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15))
abline(h = 0, lty = 2, col = "purple")
text(12, -10, labels = paste("SUC.109", " (", summary(is.de)[1]," genes)", sep = ""), col = "blue", cex = 1.2)
text(12, 0, labels = paste("(", summary(is.de)[2]," genes)", sep = ""), col = "black", cex = 1.2)
text(12, 10, labels = paste("SUC.57", " (", summary(is.de)[3]," genes)", sep = ""), col = "red", cex = 1.2)
dev.off()


# prepare the SUC 57 vs 73 plot
contr0 <- makeContrasts(contrasts = c("SUC.57","SUC.73"), levels=mydesigntime)
con <- data.frame(Contrasts = contr0[,1] - contr0[,2])
res <- glmQLFTest(fit, contrast=con)

cairo_pdf(paste("0_SUC_57vs73samples_DE_MA_q.0.01_lfc1.pdf", sep = ""), height = 5, width = 6)
is.de <- decideTestsDGE(res, p.value = 0.01 , lfc = 1)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend=NULL, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "SUC.73 vs SUC.57",ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15))
abline(h = 0, lty = 2, col = "purple")
text(12, -10, labels = paste("SUC.73", " (", summary(is.de)[1]," genes)", sep = ""), col = "blue", cex = 1.2)
text(12, 0, labels = paste("(", summary(is.de)[2]," genes)", sep = ""), col = "black", cex = 1.2)
text(12, 10, labels = paste("SUC.57", " (", summary(is.de)[3]," genes)", sep = ""), col = "red", cex = 1.2)
dev.off()


# prepare the SUC 73 vs 109 plot
contr0 <- makeContrasts(contrasts = c("SUC.73","SUC.109"), levels=mydesigntime)
con <- data.frame(Contrasts = contr0[,1] - contr0[,2])
res <- glmQLFTest(fit, contrast=con)

cairo_pdf(paste("0_SUC_73vs109samples_DE_MA_q.0.01_lfc1.pdf", sep = ""), height = 5, width = 6)
is.de <- decideTestsDGE(res, p.value = 0.01 , lfc = 1)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend=NULL, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "SUC.109 vs SUC.73", ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15))
abline(h = 0, lty = 2, col = "purple")
text(12, -10, labels = paste("SUC.109", " (", summary(is.de)[1]," genes)", sep = ""), col = "blue", cex = 1.2)
text(12, 0, labels = paste("(", summary(is.de)[2]," genes)", sep = ""), col = "black", cex = 1.2)
text(12, 10, labels = paste("SUC.73", " (", summary(is.de)[3]," genes)", sep = ""), col = "red", cex = 1.2)
dev.off()


### prepare the above plots for TBZ as well
# prepare the MA plots for the DE genes
contr0 <- makeContrasts(contrasts = c("TBZ.57","TBZ.109"), levels=mydesigntime)
con <- data.frame(Contrasts = contr0[,1] - contr0[,2])
res <- glmQLFTest(fit, contrast=con)

cairo_pdf(paste("0_TBZ_57vs109samples_DE_MA_q.0.01_lfc1.pdf", sep = ""), height = 5, width = 6)
is.de <- decideTestsDGE(res, p.value = 0.01 , lfc = 1)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend=NULL, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "TBZ.109 vs TBZ.57", ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15))
abline(h = 0, lty = 2, col = "purple")
text(12, -10, labels = paste("TBZ.109", " (", summary(is.de)[1]," genes)", sep = ""), col = "blue", cex = 1.2)
text(12, 0, labels = paste("(", summary(is.de)[2]," genes)", sep = ""), col = "black", cex = 1.2)
text(12, 10, labels = paste("TBZ.57", " (", summary(is.de)[3]," genes)", sep = ""), col = "red", cex = 1.2)
dev.off()


# prepare the TBZ 57 vs 73 plot
contr0 <- makeContrasts(contrasts = c("TBZ.57","TBZ.73"), levels=mydesigntime)
con <- data.frame(Contrasts = contr0[,1] - contr0[,2])
res <- glmQLFTest(fit, contrast=con)

cairo_pdf(paste("0_TBZ_57vs73samples_DE_MA_q.0.01_lfc1.pdf", sep = ""), height = 5, width = 6)
is.de <- decideTestsDGE(res, p.value = 0.01 , lfc = 1)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend=NULL, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "TBZ.73 vs TBZ.57",ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15))
abline(h = 0, lty = 2, col = "purple")
text(12, -10, labels = paste("TBZ.73", " (", summary(is.de)[1]," genes)", sep = ""), col = "blue", cex = 1.2)
text(12, 0, labels = paste("(", summary(is.de)[2]," genes)", sep = ""), col = "black", cex = 1.2)
text(12, 10, labels = paste("TBZ.57", " (", summary(is.de)[3]," genes)", sep = ""), col = "red", cex = 1.2)
dev.off()


# prepare the TBZ 73 vs 109 plot
contr0 <- makeContrasts(contrasts = c("TBZ.73","TBZ.109"), levels=mydesigntime)
con <- data.frame(Contrasts = contr0[,1] - contr0[,2])
res <- glmQLFTest(fit, contrast=con)

cairo_pdf(paste("0_TBZ_73vs109samples_DE_MA_q.0.01_lfc1.pdf", sep = ""), height = 5, width = 6)
is.de <- decideTestsDGE(res, p.value = 0.01 , lfc = 1)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend=NULL, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "TBZ.109 vs TBZ.73", ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15))
abline(h = 0, lty = 2, col = "purple")
text(12, -10, labels = paste("TBZ.109", " (", summary(is.de)[1]," genes)", sep = ""), col = "blue", cex = 1.2)
text(12, 0, labels = paste("(", summary(is.de)[2]," genes)", sep = ""), col = "black", cex = 1.2)
text(12, 10, labels = paste("TBZ.73", " (", summary(is.de)[3]," genes)", sep = ""), col = "red", cex = 1.2)
dev.off()






#### TBZ vs SUC ----

contr0 <- model.matrix(~0+design_tab_sel$treatment)
colnames(contr0) <- c("SUC","TBZ")
con <- makeContrasts(TBZvsSUC = TBZ - SUC, levels = contr0)

fit <- glmQLFit(my_dgelist_fin, contr0, robust=TRUE)
res <- glmQLFTest(fit, contrast = con)
# get the fold change values etc and add them to the main table
my_stats_tbl <- res$table
my_stats_tbl$q.value <- p.adjust(my_stats_tbl$PValue,"fdr")


#prepare for the log2FC matrix use the initial samples (day 57) for comparisons 
mydenomsamples <- design_tab_sel[which(design_tab_sel$treatment== "SUC" ),"sampleName"]
mysamples <- design_tab_sel[grep("SUC|TBZ",design_tab_sel$treatment),"sampleName"]

# get the log2FC table and add the comparisons in the header
mylogFC_tbl <- log2((cpm(my_dgelist_fin)[,mysamples] + 0.01)/rowMeans(cpm(my_dgelist_fin)[,mydenomsamples] + 0.01))
colnames(mylogFC_tbl) <- paste(colnames(mylogFC_tbl), "_vs_SUC_mean_log2FC", sep = "")

## prepare a heatmap ----
my_annotab_for_heatmap <- my_cnt_annot_tab_wth_sts[,grep("SEED",colnames(my_cnt_annot_tab_wth_sts))]
my_annotab_for_heatmap2 <- as.data.frame(paste(my_annotab_for_heatmap$SEED1,my_annotab_for_heatmap$SEED2,my_annotab_for_heatmap$SEED3,my_annotab_for_heatmap$SEED4,sep = "  /  "))
row.names(my_annotab_for_heatmap2) <- row.names(my_cnt_annot_tab_wth_sts)
heattab <- merge(mylogFC_tbl,my_stats_tbl, by = "row.names")
row.names(heattab) <- heattab$Row.names
heattab2 <- merge(heattab, my_annotab_for_heatmap2, by = "row.names")
heattab_names <- paste(heattab2[,1],heattab2[,ncol(heattab2)])
heattab3 <- heattab2[,3:(ncol(heattab2)-1)]
row.names(heattab3) <- heattab_names
heattab3_ord <- heattab3[order(heattab3$q.value),]

# select the highly DE and expressed genes with cutoffs: q.value of 0.001; abs log2fold change 2; logCPM 4
heattab3_ord_sel <- heattab3_ord[heattab3_ord$q.value < 0.001 & abs(heattab3_ord$logFC) > 2,]
# heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logFC) > 2,]
heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logCPM) > 4,]
heattab3_ord_sel_fin <- heattab3_ord_sel[,1:(ncol(heattab3_ord_sel)-5)]


# prepare the heatmap
library(pheatmap)
# set the breaks and colours
mybreaks <- seq(-max(abs(as.matrix(heattab3_ord_sel_fin))),max(abs(as.matrix(heattab3_ord_sel_fin))), length.out = 120)
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[6:1])(n = 119) # the nature colors
# open the graphics device, plot and close the device
cairo_pdf(file = paste("0_Heatmap_TBZ_vs_SUC_q0.001_lFC2_logCPM4.pdf", sep = ""),height = 250, width = 25)
pheatmap(as.matrix(heattab3_ord_sel_fin), color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
         #, cluster_rows = F, cluster_cols = F
         , treeheight_row = 60
         , treeheight_col = 60
)
dev.off()


# select the highly DE and expressed genes 
heattab3_ord_sel <- heattab3_ord[heattab3_ord$q.value < 0.00001 & abs(heattab3_ord$logFC) > 4,]
# heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logFC) > 4,]
heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logCPM) > 5,]
heattab3_ord_sel_fin <- heattab3_ord_sel[,grep("_vs_",colnames(heattab3_ord_sel))]

# prepare the heatmap
library(pheatmap)
# set the breaks and colours
mybreaks <- seq(-max(abs(as.matrix(heattab3_ord_sel_fin))),max(abs(as.matrix(heattab3_ord_sel_fin))), length.out = 120)
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[6:1])(n = 119) # the nature colors
# open the graphics device, plot and close the device
cairo_pdf(file = paste("0_Heatmap_TBZ_vs_SUC_q0.00001_lFC4_logCPM5.pdf", sep = ""),height = 250, width = 25)
pheatmap(as.matrix(heattab3_ord_sel_fin), color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
         #, cluster_rows = F, cluster_cols = F
         , treeheight_row = 60
         , treeheight_col = 60
)
dev.off()



# merge the tables
my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,mylogFC_tbl, by = "row.names", all =T)
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]
colnames(my_stats_tbl) <- paste("TBZ_vs_SUC_",colnames(my_stats_tbl),sep = "")
my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,my_stats_tbl, by = "row.names", all =T)
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]


# save the table
write.table(my_cnt_annot_tab_wth_sts, file = "00_my_final_tab_with_seqs_and_sampletime_contrasts_final.tsv", sep = "\t", quote = F, col.names = NA)







# prepare the MA plots for the DE genes
contr0 <- model.matrix(~0+design_tab_sel$treatment)
colnames(contr0) <- c("TBZ","SUC")
con <- makeContrasts(TBZvsSUC = TBZ - SUC, levels = contr0)

fit <- glmQLFit(my_dgelist_fin, contr0, robust=TRUE)
res <- glmQLFTest(fit, contrast = con)

cairo_pdf(paste("0_TBZvsSUCsamples_DE_MA_q.0.01_lfc1.pdf", sep = ""), height = 5, width = 6)
is.de <- decideTestsDGE(res, p.value = 0.01 , lfc = 1)

plotMD(res, status=is.de, values=c(1,-1), col=c("red","steelblue3"), legend=F, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "TBZ vs SUC", ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15))
#myselegns <- read.table("genes_of_interest_carb_cat.txt", stringsAsFactors = F)[,1]

abline(h = 0, lty = 2, col = "purple")
text(12, -10, labels = paste("SUC", " (", summary(is.de)[1]," genes)", sep = ""), col = "steelblue3", cex = 1.2)
text(12, 0, labels = paste("(", summary(is.de)[2]," genes)", sep = ""), col = "black", cex = 1.2)
text(12, 10, labels = paste("TBZ", " (", summary(is.de)[3]," genes)", sep = ""), col = "red", cex = 1.2)

dev.off()


#### Prepare the same plot for selected genes ----
is.de <- decideTestsDGE(res, p.value = 0.01 , lfc = 1)
cairo_pdf(paste("000_TBZvsSUCsamples_DE_MA_q.0.01_lfc1_for_presentation1.pdf", sep = ""), height = 6, width = 10)
plotMD(res, status=is.de, values=c(1,0,-1), col=c("red","grey70","steelblue3"), legend=F, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "TBZ vs SUC", ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15), xlim = c(0,15))
text(x=15,y=10,length(which(is.de[,1] == 1)), col = "red")
text(x=15,y=0,length(which(is.de[,1] == 0)), col = "black")
text(x=15,y=-10,length(which(is.de[,1] == -1)), col = "steelblue3")
abline(h = 0, lty = 2, col = "purple", lwd = 2)
dev.off()

cairo_pdf(paste("000_TBZvsSUCsamples_DE_MA_q.0.01_lfc1_for_presentation2.pdf", sep = ""), height = 6, width = 10)
plotMD(res, status=is.de, values=c(1,0,-1), col=c("red","grey70","steelblue3"), legend=F, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "TBZ vs SUC", ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15), xlim = c(0,15))
myselegns <- read.table("../../1_metagenome/z_accessory_files_dbs_and_executbles/genes_of_interest_carb_cat.txt", stringsAsFactors = F, row.names = 1, sep = "\t")
abline(h = 0, lty = 2, col = "purple", lwd = 2)
points(res$table[myselegns[1:8,1],2:1], cex = 2, lwd = 2)
text(res$table[myselegns[1:8,1],2:1], labels = row.names(myselegns)[1:8], pos = 1)
text(x=15,y=10,length(which(is.de[,1] == 1)), col = "red")
text(x=15,y=0,length(which(is.de[,1] == 0)), col = "black")
text(x=15,y=-10,length(which(is.de[,1] == -1)), col = "steelblue3")
dev.off()

cairo_pdf(paste("000_TBZvsSUCsamples_DE_MA_q.0.01_lfc1_for_presentation3.pdf", sep = ""), height = 6, width = 10)
plotMD(res, status=is.de, values=c(1,0,-1), col=c("red","grey70","steelblue3"), legend=F, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "TBZ vs SUC", ylab = expression(paste("log"[2],"(FC)"),sep = ""), xlab = expression(paste("mean log"[2],"(CPM)"),sep = ""), ylim = c(-15,15), xlim = c(0,15))
myselegns <- read.table("../../1_metagenome/z_accessory_files_dbs_and_executbles/genes_of_interest_carb_cat.txt", stringsAsFactors = F, row.names = 1, sep = "\t")
abline(h = 0, lty = 2, col = "purple", lwd = 2)
points(res$table[myselegns[1:8,1],2:1], cex = 2, lwd = 2)
text(res$table[myselegns[1:8,1],2:1], labels = row.names(myselegns)[1:8], pos = 1)
points(res$table[myselegns[9:17,1],2:1], pch = 2, cex = 2, lwd = 2)
text(res$table[myselegns[9:17,1],2:1], labels = row.names(myselegns)[9:17], pos = 3)
points(res$table[myselegns[18:21,1],2:1], pch = 5, cex = 2, lwd = 2)
text(res$table[myselegns[18:21,1],2:1], labels = row.names(myselegns)[18:21], pos = 3)
text(x=15,y=10,length(which(is.de[,1] == 1)), col = "red")
text(x=15,y=0,length(which(is.de[,1] == 0)), col = "black")
text(x=15,y=-10,length(which(is.de[,1] == -1)), col = "steelblue3")
dev.off()


# prep the boxplots
mydsnrelded <- design_tab_sel
mybrplttbl <- my_cnt_annot_tab_wth_sts[myselegns[,1],grep("\\.y",colnames(my_cnt_annot_tab_wth_sts))]


cairo_pdf(paste("000_TBZvsSUCsamples_DE_MA_q.0.001_lfc1_for_presentation_boxplotslog.pdf", sep = ""), height = 24, width = 3.5)
par(mfrow = c(11,2), mar = c(2,6,2,2))
for(k in 1:nrow(mybrplttbl)){
  mypltdf <- data.frame(t(mybrplttbl[k,]),fact = mydsnrelded$treatment[which(mydsnrelded$sample%in%gsub("\\.y","",colnames(mybrplttbl)))])
  colnames(mypltdf)[1] <- "var"
  mypltdf$fact <- factor(mypltdf$fact, levels = c("TBZ", "SUC"))
  boxplot(log2(mypltdf$var+0.01)~mypltdf$fact, col = c("red", "steelblue3"),frame=F, ylab = expression(paste("log"[2],"(CPM)",sep = "")), ylim = c(-10,10), main = row.names(myselegns)[which(myselegns[,1] %in% row.names(mybrplttbl)[k])], cex.main = 2, xlab = "")
}
dev.off()






























#### bubble-plots ----
## add the BinTax_id column to the large tsv file
mybinidtax <- read.table("../../1_metagenome/z_accessory_files_dbs_and_executbles/binid_taxonomy.txt", sep = "\t", header = T, stringsAsFactors = F)
library(plyr) # for the mapvalues command renaming levels of factor or character strings
my_cnt_annot_tab_wth_sts$BinTaxID <- mapvalues(my_cnt_annot_tab_wth_sts$binid, from = mybinidtax$binid, to = mybinidtax$BinTaxID)

## prepare the dataset with all the information including the gene incidence column
mydatabubble <- data.frame(row.names = paste(row.names(my_cnt_annot_tab_wth_sts),my_cnt_annot_tab_wth_sts$SEED4, sep = " "), BinTaxID = my_cnt_annot_tab_wth_sts$BinTaxID, SEED1 = my_cnt_annot_tab_wth_sts$SEED1, SEED2 = my_cnt_annot_tab_wth_sts$SEED2, SEED3 = my_cnt_annot_tab_wth_sts$SEED3, SEED4 = my_cnt_annot_tab_wth_sts$SEED4, feature_cnts = my_cnt_annot_tab_wth_sts$transcripts)
## add the meanCPM of the various tests (take care of the NAs)
tmp_cpm_vals <- my_cnt_annot_tab_wth_sts[,grep("\\.y",colnames(my_cnt_annot_tab_wth_sts))]

# add the TBZ vs SUC contrast related info
mydatabubble$TBZperSUC.TBZ.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZperSUC.TBZ.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]),grep("TBZ", design_tab_sel$treatment)])
mydatabubble$TBZperSUC.SUC.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZperSUC.SUC.mean.CPM[complete.cases(tmp_cpm_vals[,1])]  <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]),grep("SUC", design_tab_sel$treatment)])
mydatabubble$TBZperSUC.q.value <- my_cnt_annot_tab_wth_sts$TBZ_vs_SUC_q.value


# add the TBZ.73 and 109 vs TBZ.57 contrasts related info (the 73 and 109 are to be plotted in relation to the 57)
mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "57")])
mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "73")])
mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "109")])
mydatabubble$TBZsamples_vs_TBZ.57_q.value <- my_cnt_annot_tab_wth_sts$TBZsamples_vs_TBZ.57_q.value


# add the SUC.73 and 109 vs SUC.57 contrasts related info (the 73 and 109 are to be plotted in relation to the 57)
mydatabubble$SUC.73.109perSUC.57.SUC57.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$SUC.73.109perSUC.57.SUC57.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "SUC" & design_tab_sel$time == "57")])
mydatabubble$SUC.73.109perSUC.57.SUC73.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$SUC.73.109perSUC.57.SUC73.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "SUC" & design_tab_sel$time == "73")])
mydatabubble$SUC.73.109perSUC.57.SUC109.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$SUC.73.109perSUC.57.SUC109.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "SUC" & design_tab_sel$time == "109")])
mydatabubble$SUCsamples_vs_SUC.57_q.value <- my_cnt_annot_tab_wth_sts$SUCsamples_vs_SUC.57_q.value



#### bubbule-plots in loops ----

# set all mydatabubble categories CPM values where the q.value is above the threshold to zero
mydatabubble_bckp <- mydatabubble
# I performed this step to avoid re-running the above each time the thresholds change
mydatabubble <- mydatabubble_bckp
q_threshold <- 0.01

mydatabubble$TBZperSUC.TBZ.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0
mydatabubble$TBZperSUC.SUC.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0


mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0
mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0
mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0

mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM[mydatabubble$SUCsamples_vs_SUC.57_q.value > q_threshold] <- 0
mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM[mydatabubble$SUCsamples_vs_SUC.57_q.value > q_threshold] <- 0


# prepare the metagenomics associated plot for the bins
library(ggplot2)
library(ggtree)

for(i in 1:4){
  mychoice <- paste("SEED",i, sep = "")
  
  
  # prep the collective table with the aggregated features and CPMs and calculate the log2FCs and add the up/downreg colours
  mydatabubble_for_metagenome <- aggregate(mydatabubble[complete.cases(mydatabubble[,mychoice]),c(6:8,10:12,14:16,18:19)], by = list(mydatabubble[complete.cases(mydatabubble[,mychoice]),mychoice],mydatabubble[complete.cases(mydatabubble[,mychoice]),"BinTaxID"]), sum)
  
  # TBZ vs SUC
  mydatabubble_for_metagenome$TBZperSUC.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZperSUC.TBZ.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZperSUC.SUC.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  mydatabubble_for_metagenome$TBZperSUC.log2FC[which(mydatabubble_for_metagenome$TBZperSUC.log2FC == 0)] <- NA # I did this to avoid the white colour bubbles
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZperSUC.log2FC)
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls[which(mydatabubble_for_metagenome$TBZperSUC.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls[which(mydatabubble_for_metagenome$TBZperSUC.log2FC < 0)] <- "green"
  # add the absolute values
  mydatabubble_for_metagenome$TBZperSUC.log2FC <- abs(mydatabubble_for_metagenome$TBZperSUC.log2FC)
  

  # TBZ 73 vs TBZ 57
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC == 0)] <- NA
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC)
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC < 0)] <- "green"
  
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC <- abs(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC)

    
  
  # TBZ 109 vs TBZ 57
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  
  if(length(which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC == 0)) > 0){
    mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC == 0)] <- NA
  }
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC)
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC < 0)] <- "green" 
  
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC <- abs(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC)

    
  
  
  # SUC 73 vs SUC 57
  mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$SUC.73.109perSUC.57.SUC73.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$SUC.73.109perSUC.57.SUC57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FC[which(mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FC == 0)] <- NA
  mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FCcls <- as.character(mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FC)
  mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FCcls <- NA
  mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FCcls[which(mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FCcls[which(mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FC < 0)] <- "green"
  
  mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FC <- abs(mydatabubble_for_metagenome$SUC73perSUC57.SUC73.log2FC)
  
  
  
  # SUC 109 vs SUC 57
  mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$SUC.73.109perSUC.57.SUC109.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$SUC.73.109perSUC.57.SUC57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  
  if(length(which(mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FC == 0)) > 0){
    mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FC[which(mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FC == 0)] <- NA
  }
  mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FCcls <- as.character(mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FC)
  mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FCcls <- NA
  mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FCcls[which(mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FCcls[which(mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FC < 0)] <- "green" 
  
  mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FC <- abs(mydatabubble_for_metagenome$SUC109perSUC57.SUC109.log2FC)
  
  


  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[,-grep("mean.CPM",colnames(mydatabubble_for_metagenome))]
  
  colnames(mydatabubble_for_metagenome)[1:3] <- c("functional group","bin", "counts")
  # order the matrix according to the binid
  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[order(as.numeric(gsub("[A-Z].+|_.+","",gsub("Bin_","",gsub("unbinned","1000000",gsub("unbinned_Hydrogenophaga","999999",mydatabubble_for_metagenome$bin))))),gsub("Bin_","",mydatabubble_for_metagenome$bin)),]
  mydatabubble_for_metagenome$bin <- factor(mydatabubble_for_metagenome$bin, levels = unique(mydatabubble_for_metagenome$bin))
  
  
  
  if(mychoice == "SEED4"){
    myheight <- 150
    mywidth <- 120
  } else if(mychoice == "SEED3") {
    myheight <- 85
    mywidth <- 70
  } else if (mychoice == "SEED2") {
    myheight <- 20
    mywidth <- 65
  } else if (mychoice == "SEED1") {
    myheight <- 7
    mywidth <- 55
  }
  
  cairo_pdf(paste("bubble_plot_",mychoice,"_all_q",q_threshold,".pdf", sep = ""), width = mywidth, height = myheight)
  
  g1 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`counts`)) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("Metagenome"))
  
  g2 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    #geom_point(aes(size=`TBZperSUC.log2FC`, colour = `TBZperSUC.log2FCcls`)) + 
    geom_point(aes(size=`TBZperSUC.log2FC`, colour = I(TBZperSUC.log2FCcls))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    # aes(fill = TBZperSUC.log2FCcls) +
    # scale_fill_manual(values = TBZperSUC.log2FCclsRGB) +
    ggtitle(paste("TBZ/SUC gene expression"))
  
  g3 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`TBZ73perTBZ57.TBZ73.log2FC`, colour = I(TBZ73perTBZ57.TBZ73.log2FCcls))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("TBZ73/TBZ57 gene expression"))
  
  g4 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`TBZ109perTBZ57.TBZ109.log2FC`, colour = I(TBZ109perTBZ57.TBZ109.log2FCcls))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("TBZ109/TBZ57 gene expression"))
  

  g5 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`SUC73perSUC57.SUC73.log2FC`, colour = I(SUC73perSUC57.SUC73.log2FCcls))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("SUC73/SUC57 gene expression"))
  
  g6 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`SUC109perSUC57.SUC109.log2FC`, colour = I(SUC109perSUC57.SUC109.log2FCcls))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("SUC109/SUC57 gene expression"))
  
  
  
    # g5 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    # geom_point(aes(size=`SUC57perSUC109.log2FC`, colour = I(SUC57perSUC109.log2FCcls))) + 
    # #theme_bw() + 
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    # labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    # ggtitle(paste("SUC57/SUC109 gene expression"))
  
  
  multiplot(g1, g2, g3, g4, g5, g6, ncol=6)
  
  
  
  
  dev.off()
  
}




















































#### bubble-plots for the enriched SEED category----
## add the BinTax_id column to the large tsv file
mybinidtax <- read.table("../../1_metagenome/z_accessory_files_dbs_and_executbles/binid_taxonomy.txt", sep = "\t", header = T, stringsAsFactors = F)
library(plyr) # for the mapvalues command renaming levels of factor or character strings
my_cnt_annot_tab_wth_sts$BinTaxID <- mapvalues(my_cnt_annot_tab_wth_sts$binid, from = mybinidtax$binid, to = mybinidtax$BinTaxID)

## prepare the dataset with all the information including the gene incidence column
mydatabubble <- data.frame(row.names = paste(row.names(my_cnt_annot_tab_wth_sts),my_cnt_annot_tab_wth_sts$SEED4, sep = " "), BinTaxID = my_cnt_annot_tab_wth_sts$BinTaxID, SEED1 = my_cnt_annot_tab_wth_sts$SEED1, SEED2 = my_cnt_annot_tab_wth_sts$SEED2, SEED3 = my_cnt_annot_tab_wth_sts$SEED3, SEED4 = my_cnt_annot_tab_wth_sts$SEED4, feature_cnts = my_cnt_annot_tab_wth_sts$transcripts)

## add the meanCPM of the various tests (take care of the NAs)
tmp_cpm_vals <- my_cnt_annot_tab_wth_sts[,grep("\\.y",colnames(my_cnt_annot_tab_wth_sts))]

# add the TBZ vs SUC contrast related info
mydatabubble$TBZperSUC.TBZ.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZperSUC.TBZ.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]),grep("TBZ", design_tab_sel$treatment)])
mydatabubble$TBZperSUC.SUC.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZperSUC.SUC.mean.CPM[complete.cases(tmp_cpm_vals[,1])]  <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]),grep("SUC", design_tab_sel$treatment)])
mydatabubble$TBZperSUC.q.value <- my_cnt_annot_tab_wth_sts$TBZ_vs_SUC_q.value


# add the TBZ.73 and 109 vs TBZ.57 contrasts related info (the 73 and 109 are to be plotted in relation to the 57)
mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "57")])
mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "73")])
mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "109")])
mydatabubble$TBZsamples_vs_TBZ.57_q.value <- my_cnt_annot_tab_wth_sts$TBZsamples_vs_TBZ.57_q.value


# add the SUC.57 vs SUC.109 contrast related info 
mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "SUC" & design_tab_sel$time == "57")])
mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "SUC" & design_tab_sel$time == "109")])
mydatabubble$SUCsamples_vs_SUC.57_q.value <- my_cnt_annot_tab_wth_sts$SUCsamples_vs_SUC.57_q.value



#### bubble-plots in loops ----

# set all mydatabubble categories CPM values where the q.value is above the threshold to zero
mydatabubble_bckp <- mydatabubble
# I performed this step to avoid re-running the above each time the thresholds change
mydatabubble <- mydatabubble_bckp
q_threshold <- 0.05

mydatabubble$TBZperSUC.TBZ.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0
mydatabubble$TBZperSUC.SUC.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0


mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0
mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0
mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0

mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM[mydatabubble$SUCsamples_vs_SUC.57_q.value > q_threshold] <- 0
mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM[mydatabubble$SUCsamples_vs_SUC.57_q.value > q_threshold] <- 0


# prepare the metagenomics associated plot for the bins
library(ggplot2)
library(ggtree)

for(i in 1:4){
  mychoice <- paste("SEED",i, sep = "")
  
  
  # prep the collective table with the aggregated features and CPMs and calculate the log2FCs and add the up/downreg colours
  mydatabubble_for_metagenome <- aggregate(mydatabubble[complete.cases(mydatabubble[,mychoice]),c(6:8,10:12,14:15)], by = list(mydatabubble[complete.cases(mydatabubble[,mychoice]),mychoice],mydatabubble[complete.cases(mydatabubble[,mychoice]),"BinTaxID"]), sum)
  
  # TBZ vs SUC
  mydatabubble_for_metagenome$TBZperSUC.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZperSUC.TBZ.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZperSUC.SUC.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  mydatabubble_for_metagenome$TBZperSUC.log2FC[which(mydatabubble_for_metagenome$TBZperSUC.log2FC == 0)] <- NA # I did this to avoid the white colour bubbles
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZperSUC.log2FC)
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls[which(mydatabubble_for_metagenome$TBZperSUC.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls[which(mydatabubble_for_metagenome$TBZperSUC.log2FC < 0)] <- "green"
  # add the absolute values
  mydatabubble_for_metagenome$TBZperSUC.log2FC <- abs(mydatabubble_for_metagenome$TBZperSUC.log2FC)
  
  
  # TBZ 73 vs TBZ 57
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC == 0)] <- NA
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC)
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC < 0)] <- "green"
  
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC <- abs(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC)
  
  
  
  # TBZ 109 vs TBZ 57
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  
  if(length(which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC == 0)) > 0){
    mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC == 0)] <- NA
  }
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC)
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC < 0)] <- "green" 
  
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC <- abs(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC)
  
  
  
  
  # SUC SUC 57 vs 109
  mydatabubble_for_metagenome$SUC57perSUC109.log2FC[complete.cases(mydatabubble_for_metagenome)] <- as.numeric(log2((mydatabubble_for_metagenome$SUC.109perSUC.57.SUC57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$SUC.109perSUC.57.SUC109.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)))
  mydatabubble_for_metagenome$SUC57perSUC109.log2FC <- as.numeric(mydatabubble_for_metagenome$SUC57perSUC109.log2FC)
  
  if(length(which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC == 0))>0){
    mydatabubble_for_metagenome$SUC57perSUC109.log2FC[which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC == 0)] <- NA
    
  }
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls <- as.character(mydatabubble_for_metagenome$SUC57perSUC109.log2FC)
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls <- NA
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls[which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls[which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC < 0)] <- "green"
  
  mydatabubble_for_metagenome$SUC57perSUC109.log2FC <- abs(mydatabubble_for_metagenome$SUC57perSUC109.log2FC)  
  
  
  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[,-grep("mean.CPM",colnames(mydatabubble_for_metagenome))]
  
  colnames(mydatabubble_for_metagenome)[1:3] <- c("functional group","bin", "counts")
  # order the matrix according to the binid
  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[order(as.numeric(gsub("[A-Z].+|_.+","",gsub("Bin_","",gsub("unbinned","1000000",gsub("unbinned_Hydrogenophaga","999999",mydatabubble_for_metagenome$bin))))),gsub("Bin_","",mydatabubble_for_metagenome$bin)),]
  mydatabubble_for_metagenome$bin <- factor(mydatabubble_for_metagenome$bin, levels = unique(mydatabubble_for_metagenome$bin))
  
  
  
  if(mychoice == "SEED4"){
    myheight <- 150
    mywidth <- 100
  } else if(mychoice == "SEED3") {
    myheight <- 85
    mywidth <- 60
  } else if (mychoice == "SEED2") {
    myheight <- 20
    mywidth <- 55
  } else if (mychoice == "SEED1") {
    myheight <- 7
    mywidth <- 45
  }
  
  cairo_pdf(paste("bubble_plot_enriched_",mychoice,"_all_q",q_threshold,".pdf", sep = ""), width = mywidth, height = myheight)
  
  g1 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`counts`)) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("Metagenome"))
  
  g2 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    #geom_point(aes(size=`TBZperSUC.log2FC`, colour = `TBZperSUC.log2FCcls`)) + 
    geom_point(aes(size=`TBZperSUC.log2FC`, colour = I(TBZperSUC.log2FCcls))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    # aes(fill = TBZperSUC.log2FCcls) +
    # scale_fill_manual(values = TBZperSUC.log2FCclsRGB) +
    ggtitle(paste("TBZ/SUC gene expression"))
  
  g3 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`TBZ73perTBZ57.TBZ73.log2FC`, colour = I(TBZ73perTBZ57.TBZ73.log2FCcls))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("TBZ73/TBZ57 gene expression"))
  
  g4 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`TBZ109perTBZ57.TBZ109.log2FC`, colour = I(TBZ109perTBZ57.TBZ109.log2FCcls))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("TBZ109/TBZ57 gene expression"))
  
  g5 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`SUC57perSUC109.log2FC`, colour = I(SUC57perSUC109.log2FCcls))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("SUC57/SUC109 gene expression"))
  
  
  multiplot(g1, g2, g3, g4, g5, ncol=5)
  
  
  
  
  dev.off()
  
}

















#### SEED bubble plots with keyword search ----
## prepare the dataset with all the information including the gene incidence column

mydatabubble <- data.frame(row.names = paste(my_cnt_annot_tab_wth_sts$ID,my_cnt_annot_tab_wth_sts$SEED4, sep = " "), BinTaxID = my_cnt_annot_tab_wth_sts$BinTaxID, SEED1 = my_cnt_annot_tab_wth_sts$SEED1, SEED2 = my_cnt_annot_tab_wth_sts$SEED2, SEED3 = my_cnt_annot_tab_wth_sts$SEED3, SEED4 = my_cnt_annot_tab_wth_sts$SEED4, feature_cnts = my_cnt_annot_tab_wth_sts$transcripts)
## add the meanCPM of the various tests (take care of the NAs)
tmp_cpm_vals <- my_cnt_annot_tab_wth_sts[,grep("\\.y",colnames(my_cnt_annot_tab_wth_sts))]

# add the TBZ vs SUC contrast related info
mydatabubble$TBZperSUC.TBZ.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZperSUC.TBZ.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]),grep("TBZ", design_tab_sel$treatment)])
mydatabubble$TBZperSUC.SUC.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZperSUC.SUC.mean.CPM[complete.cases(tmp_cpm_vals[,1])]  <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]),grep("SUC", design_tab_sel$treatment)])
mydatabubble$TBZperSUC.q.value <- my_cnt_annot_tab_wth_sts$TBZ_vs_SUC_q.value


# add the TBZ.73 and 109 vs TBZ.57 contrasts related info (the 73 and 109 are to be plotted in relation to the 57)
mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "57")])
mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "73")])
mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "109")])
mydatabubble$TBZsamples_vs_TBZ.57_q.value <- my_cnt_annot_tab_wth_sts$TBZsamples_vs_TBZ.57_q.value


# add the SUC.57 vs SUC.109 contrast related info 
mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "SUC" & design_tab_sel$time == "57")])
mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM <- tmp_cpm_vals[,1]
mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "SUC" & design_tab_sel$time == "109")])
mydatabubble$SUCsamples_vs_SUC.57_q.value <- my_cnt_annot_tab_wth_sts$SUCsamples_vs_SUC.57_q.value



#### bubble-plots in loops ----

# set all mydatabubble categories CPM values where the q.value is above the threshold to zero
mydatabubble_bckp <- mydatabubble
# I performed this step to avoid re-running the above each time the thresholds change
mydatabubble <- mydatabubble_bckp
q_threshold <- 0.05

mydatabubble$TBZperSUC.TBZ.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0
mydatabubble$TBZperSUC.SUC.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0


mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0
mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0
mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0

mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM[mydatabubble$SUCsamples_vs_SUC.57_q.value > q_threshold] <- 0
mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM[mydatabubble$SUCsamples_vs_SUC.57_q.value > q_threshold] <- 0


# prepare the metagenomics associated plot for the bins
library(ggplot2)
library(ggtree)

#for(i in 1:4){
i <- 4

  mychoice <- paste("SEED",i, sep = "")
  
  
  # prep the collective table with the aggregated features and CPMs and calculate the log2FCs and add the up/downreg colours
  mydatabubble_for_metagenome <- aggregate(mydatabubble[complete.cases(mydatabubble[,mychoice]),c(6:8,10:12,14:15)], by = list(mydatabubble[complete.cases(mydatabubble[,mychoice]),mychoice],mydatabubble[complete.cases(mydatabubble[,mychoice]),"BinTaxID"]), sum)
  
  # TBZ vs SUC
  mydatabubble_for_metagenome$TBZperSUC.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZperSUC.TBZ.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZperSUC.SUC.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  mydatabubble_for_metagenome$TBZperSUC.log2FC[which(mydatabubble_for_metagenome$TBZperSUC.log2FC == 0)] <- NA # I did this to avoid the white colour bubbles
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZperSUC.log2FC)
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls[which(mydatabubble_for_metagenome$TBZperSUC.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls[which(mydatabubble_for_metagenome$TBZperSUC.log2FC < 0)] <- "green"
  # add the absolute values
  mydatabubble_for_metagenome$TBZperSUC.log2FC <- abs(mydatabubble_for_metagenome$TBZperSUC.log2FC)
  
  
  # TBZ 73 vs TBZ 57
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC == 0)] <- NA
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC)
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC < 0)] <- "green"
  
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC <- abs(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC)
  
  
  
  # TBZ 109 vs TBZ 57
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  
  if(length(which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC == 0)) > 0){
    mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC == 0)] <- NA
  }
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC)
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC < 0)] <- "green" 
  
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC <- abs(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC)
  
  
  
  
  # SUC SUC 57 vs 109
  mydatabubble_for_metagenome$SUC57perSUC109.log2FC[complete.cases(mydatabubble_for_metagenome)] <- as.numeric(log2((mydatabubble_for_metagenome$SUC.109perSUC.57.SUC57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$SUC.109perSUC.57.SUC109.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)))
  mydatabubble_for_metagenome$SUC57perSUC109.log2FC <- as.numeric(mydatabubble_for_metagenome$SUC57perSUC109.log2FC)
  
  if(length(which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC == 0))>0){
    mydatabubble_for_metagenome$SUC57perSUC109.log2FC[which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC == 0)] <- NA
    
  }
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls <- as.character(mydatabubble_for_metagenome$SUC57perSUC109.log2FC)
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls <- NA
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls[which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls[which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC < 0)] <- "green"
  
  mydatabubble_for_metagenome$SUC57perSUC109.log2FC <- abs(mydatabubble_for_metagenome$SUC57perSUC109.log2FC)  
  
  
  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[,-grep("mean.CPM",colnames(mydatabubble_for_metagenome))]
  
  colnames(mydatabubble_for_metagenome)[1:3] <- c("functional group","bin", "counts")
  
  # order the matrix according to the binid
  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[order(as.numeric(gsub("[A-Z].+|_.+","",gsub("Bin_","",gsub("unbinned","1000000",gsub("unbinned_Hydrogenophaga","999999",mydatabubble_for_metagenome$bin))))),gsub("Bin_","",mydatabubble_for_metagenome$bin)),]
  mydatabubble_for_metagenome$bin <- factor(mydatabubble_for_metagenome$bin, levels = unique(mydatabubble_for_metagenome$bin))
  
  ### perform the keyword selection
  # create a checkpoint
  mydatabubble_for_metagenome_2 <- mydatabubble_for_metagenome
  
  
  # set the keywords
  my_grps = c("arom","signal","auxotr","mobile","transport")
  for(my_grp in my_grps){
    if(my_grp == "arom"){
      assign(paste("my_keywds.",my_grp,sep = ""),"degrad|arom|toluen|catecho|benzo|polycy|napht|phenanth|pyren|imidazol|gallat|fluoren|gentisate|nicotinat|nitroaro|phenol|protocatech|quinat|triazine|shikimat|urate|vanillin|xylene|dehydrodico|bipheny|quinoli|phenone|cresol|cumate|cymene|assimi|chlorin|fluor|ketoadipate|cleavag")
    } else if(my_grp == "signal"){
      assign(paste("my_keywds.",my_grp,sep = ""),"sign|homoser|lactone|biofilm|recepto|transduc|activat|finger")
    } else if(my_grp == "auxotr"){
      assign(paste("my_keywds.",my_grp,sep = ""),"methion|aminoacid|cobalam|aminoacid|biofilm|auxotr|transduc|vitamin|B12|B6|glucosinol|thiamin|ascorb|flavi|pantothen|quinol|quinon|coenzyme|biotin")
    } else if(my_grp == "mobile"){
      assign(paste("my_keywds.",my_grp,sep = ""),"mobile|plasmid|vir|horiz|secretion|conjug|pilus|fertil|phage|transpos|insertion|integras")
    } else if(my_grp == "transport"){
      assign(paste("my_keywds.",my_grp,sep = ""),"transpor|pump|efflux|influx|export|antiport")
    }
    
    
    
    
    mydatabubble_for_metagenome <- mydatabubble_for_metagenome_2[grep(get(paste("my_keywds.",my_grp,sep = "")),mydatabubble_for_metagenome_2$`functional group`,ignore.case=TRUE),]
    
    #if(mychoice == "SEED4"){
    # myheight <- 120
    myheight <- 10
    mywidth <- 60
    # } else if(mychoice == "SEED3") {
    #   myheight <- 70
    #   mywidth <- 65
    # } else if (mychoice == "SEED2") {
    #   myheight <- 20
    #   mywidth <- 65
    # } else if (mychoice == "SEED1") {
    #   myheight <- 30
    #   mywidth <- 65
    # }
    
    cairo_pdf(paste("bubble_plot_",my_grp,"_keywd_",mychoice,"_all_q",q_threshold,".pdf", sep = ""), width = mywidth, height = myheight)
    
    g1 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
      geom_point(aes(size=`counts`)) + 
      #theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
      ggtitle(paste("Metagenome"))
    
    g2 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
      #geom_point(aes(size=`TBZperSUC.log2FC`, colour = `TBZperSUC.log2FCcls`)) + 
      geom_point(aes(size=`TBZperSUC.log2FC`, colour = I(TBZperSUC.log2FCcls))) + 
      #theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
      # aes(fill = TBZperSUC.log2FCcls) +
      # scale_fill_manual(values = TBZperSUC.log2FCclsRGB) +
      ggtitle(paste("TBZ/SUC gene expression"))
    
    g3 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
      geom_point(aes(size=`TBZ73perTBZ57.TBZ73.log2FC`, colour = I(TBZ73perTBZ57.TBZ73.log2FCcls))) + 
      #theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
      ggtitle(paste("TBZ73/TBZ57 gene expression"))
    
    g4 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
      geom_point(aes(size=`TBZ109perTBZ57.TBZ109.log2FC`, colour = I(TBZ109perTBZ57.TBZ109.log2FCcls))) + 
      #theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
      ggtitle(paste("TBZ109/TBZ57 gene expression"))
    
    g5 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
      geom_point(aes(size=`SUC57perSUC109.log2FC`, colour = I(SUC57perSUC109.log2FCcls))) + 
      #theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
      ggtitle(paste("SUC57/SUC109 gene expression"))
    
    
    multiplot(g1, g2, g3, g4, g5, ncol=5)

    dev.off()
  }
  
  
#}



















  
  
  
  
  
  #### bubble plots with keyword search only TBZ vs SUC----
  ## prepare the dataset with all the information including the gene incidence column
  mydatabubble <- data.frame(row.names = paste(my_cnt_annot_tab_wth_sts$ID,my_cnt_annot_tab_wth_sts$SEED4, sep = " "), BinTaxID = my_cnt_annot_tab_wth_sts$BinTaxID, SEED1 = my_cnt_annot_tab_wth_sts$SEED1, SEED2 = my_cnt_annot_tab_wth_sts$SEED2, SEED3 = my_cnt_annot_tab_wth_sts$SEED3, SEED4 = my_cnt_annot_tab_wth_sts$SEED4, feature_cnts = my_cnt_annot_tab_wth_sts$transcripts)
  ## add the meanCPM of the various tests (take care of the NAs)
  tmp_cpm_vals <- my_cnt_annot_tab_wth_sts[,grep("\\.y",colnames(my_cnt_annot_tab_wth_sts))]
  
  # add the TBZ vs SUC contrast related info
  mydatabubble$TBZperSUC.TBZ.mean.CPM <- tmp_cpm_vals[,1]
  mydatabubble$TBZperSUC.TBZ.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]),grep("TBZ", design_tab_sel$treatment)])
  mydatabubble$TBZperSUC.SUC.mean.CPM <- tmp_cpm_vals[,1]
  mydatabubble$TBZperSUC.SUC.mean.CPM[complete.cases(tmp_cpm_vals[,1])]  <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]),grep("SUC", design_tab_sel$treatment)])
  mydatabubble$TBZperSUC.q.value <- my_cnt_annot_tab_wth_sts$TBZ_vs_SUC_q.value
  
  
  # add the TBZ.73 and 109 vs TBZ.57 contrasts related info (the 73 and 109 are to be plotted in relation to the 57)
  mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM <- tmp_cpm_vals[,1]
  mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "57")])
  mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM <- tmp_cpm_vals[,1]
  mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "73")])
  mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM <- tmp_cpm_vals[,1]
  mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "TBZ" & design_tab_sel$time == "109")])
  mydatabubble$TBZsamples_vs_TBZ.57_q.value <- my_cnt_annot_tab_wth_sts$TBZsamples_vs_TBZ.57_q.value
  
  
  # add the SUC.57 vs SUC.109 contrast related info 
  mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM <- tmp_cpm_vals[,1]
  mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "SUC" & design_tab_sel$time == "57")])
  mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM <- tmp_cpm_vals[,1]
  mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM[complete.cases(tmp_cpm_vals[,1])] <- rowMeans(tmp_cpm_vals[complete.cases(tmp_cpm_vals[,1]), which(design_tab_sel$treatment == "SUC" & design_tab_sel$time == "109")])
  mydatabubble$SUCsamples_vs_SUC.57_q.value <- my_cnt_annot_tab_wth_sts$SUCsamples_vs_SUC.57_q.value
  
  
  
  #### bubble-plots in loops ----
  
  # set all mydatabubble categories CPM values where the q.value is above the threshold to zero
  mydatabubble_bckp <- mydatabubble
  # I performed this step to avoid re-running the above each time the thresholds change
  mydatabubble <- mydatabubble_bckp
  q_threshold <- 0.01
  
  mydatabubble$TBZperSUC.TBZ.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0
  mydatabubble$TBZperSUC.SUC.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0
  
  
  mydatabubble$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0
  mydatabubble$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0
  mydatabubble$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[mydatabubble$TBZsamples_vs_TBZ.57_q.value > q_threshold] <- 0
  
  mydatabubble$SUC.109perSUC.57.SUC57.mean.CPM[mydatabubble$SUCsamples_vs_SUC.57_q.value > q_threshold] <- 0
  mydatabubble$SUC.109perSUC.57.SUC109.mean.CPM[mydatabubble$SUCsamples_vs_SUC.57_q.value > q_threshold] <- 0
  
  
  # prepare the metagenomics associated plot for the bins
  library(ggplot2)
  library(ggtree)
  
  #for(i in 1:4){
  i <- 4
  
  mychoice <- paste("SEED",i, sep = "")
  
  
  # prep the collective table with the aggregated features and CPMs and calculate the log2FCs and add the up/downreg colours
  mydatabubble_for_metagenome <- aggregate(mydatabubble[complete.cases(mydatabubble[,mychoice]),c(6:8,10:12,14:15)], by = list(mydatabubble[complete.cases(mydatabubble[,mychoice]),mychoice],mydatabubble[complete.cases(mydatabubble[,mychoice]),"BinTaxID"]), sum)
  
  # TBZ vs SUC
  mydatabubble_for_metagenome$TBZperSUC.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZperSUC.TBZ.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZperSUC.SUC.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  mydatabubble_for_metagenome$TBZperSUC.log2FC[which(mydatabubble_for_metagenome$TBZperSUC.log2FC == 0)] <- NA # I did this for preventing the white colour bubbles
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZperSUC.log2FC)
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls[which(mydatabubble_for_metagenome$TBZperSUC.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZperSUC.log2FCcls[which(mydatabubble_for_metagenome$TBZperSUC.log2FC < 0)] <- "green"
  # add the absolute values
  mydatabubble_for_metagenome$TBZperSUC.log2FC <- abs(mydatabubble_for_metagenome$TBZperSUC.log2FC)
  
  
  # TBZ 73 vs TBZ 57
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ73.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC == 0)] <- NA
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC)
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FCcls[which(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC < 0)] <- "green"
  
  mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC <- abs(mydatabubble_for_metagenome$TBZ73perTBZ57.TBZ73.log2FC)
  
  
  
  # TBZ 109 vs TBZ 57
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC[complete.cases(mydatabubble_for_metagenome)] <- log2((mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ109.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$TBZ.73.109perTBZ.57.TBZ57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001))
  
  if(length(which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC == 0)) > 0){
    mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC == 0)] <- NA
  }
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls <- as.character(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC)
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls <- NA
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FCcls[which(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC < 0)] <- "green" 
  
  mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC <- abs(mydatabubble_for_metagenome$TBZ109perTBZ57.TBZ109.log2FC)
  
  
  
  
  # SUC SUC 57 vs 109
  mydatabubble_for_metagenome$SUC57perSUC109.log2FC[complete.cases(mydatabubble_for_metagenome)] <- as.numeric(log2((mydatabubble_for_metagenome$SUC.109perSUC.57.SUC57.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)/(mydatabubble_for_metagenome$SUC.109perSUC.57.SUC109.mean.CPM[complete.cases(mydatabubble_for_metagenome)] + 0.001)))
  mydatabubble_for_metagenome$SUC57perSUC109.log2FC <- as.numeric(mydatabubble_for_metagenome$SUC57perSUC109.log2FC)
  
  if(length(which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC == 0))>0){
    mydatabubble_for_metagenome$SUC57perSUC109.log2FC[which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC == 0)] <- NA
    
  }
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls <- as.character(mydatabubble_for_metagenome$SUC57perSUC109.log2FC)
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls <- NA
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls[which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC > 0)] <- "red"
  mydatabubble_for_metagenome$SUC57perSUC109.log2FCcls[which(mydatabubble_for_metagenome$SUC57perSUC109.log2FC < 0)] <- "green"
  
  mydatabubble_for_metagenome$SUC57perSUC109.log2FC <- abs(mydatabubble_for_metagenome$SUC57perSUC109.log2FC)  
  
  
  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[,-grep("mean.CPM",colnames(mydatabubble_for_metagenome))]
  
  colnames(mydatabubble_for_metagenome)[1:3] <- c("functional group","bin", "counts")
  
  # order the matrix according to the binid
  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[order(as.numeric(as.factor(mydatabubble_for_metagenome$bin))),]
  mydatabubble_for_metagenome$bin <- factor(mydatabubble_for_metagenome$bin, levels = unique(mydatabubble_for_metagenome$bin))
  
  ### perform the keyword selection
  # create a checkpoint
  mydatabubble_for_metagenome_2 <- mydatabubble_for_metagenome
  
  
  my_grps = c("arom","signal","auxotr","mobile","transport")
  for(my_grp in my_grps){
    if(my_grp == "arom"){
      assign(paste("my_keywds.",my_grp,sep = ""),"degrad|arom|toluen|catecho|benzo|polycy|napht|phenanth|pyren|imidazol|gallat|fluoren|gentisate|nicotinat|nitroaro|phenol|protocatech|quinat|triazine|shikimat|urate|vanillin|xylene|dehydrodico|bipheny|quinoli|phenone|cresol|cumate|cymene|assimi|chlorin|fluor|ketoadipate|cleavag")
      myheight <- 6
      mywidth <- 8
    } else if(my_grp == "signal"){
      assign(paste("my_keywds.",my_grp,sep = ""),"sign|homoser|lactone|biofilm|recepto|transduc|activat|finger")
      myheight <- 4
      mywidth <- 4
    } else if(my_grp == "auxotr"){
      assign(paste("my_keywds.",my_grp,sep = ""),"methion|aminoacid|cobalam|aminoacid|biofilm|auxotr|transduc|vitamin|B12|B6|glucosinol|thiamin|ascorb|flavi|pantothen|quinol|quinon|coenzyme|biotin")
      myheight <- 6
      mywidth <- 7
    } else if(my_grp == "mobile"){
      assign(paste("my_keywds.",my_grp,sep = ""),"mobile|plasmid|vir|horiz|secretion|conjug|pilus|fertil|phage|transpos|insertion|integras")
    } else if(my_grp == "transport"){
      assign(paste("my_keywds.",my_grp,sep = ""),"transpor|pump|efflux|influx|export|antiport")
      myheight <- 6
      mywidth <- 11
    }
    
    
    
    
    mydatabubble_for_metagenome <- mydatabubble_for_metagenome_2[grep(get(paste("my_keywds.",my_grp,sep = "")),mydatabubble_for_metagenome_2$`functional group`,ignore.case=TRUE),]
    
    #if(mychoice == "SEED4"){
    # } else if(mychoice == "SEED3") {
    #   myheight <- 70
    #   mywidth <- 65
    # } else if (mychoice == "SEED2") {
    #   myheight <- 20
    #   mywidth <- 65
    # } else if (mychoice == "SEED1") {
    #   myheight <- 30
    #   mywidth <- 65
    # }
    
    cairo_pdf(paste("bubble_plot_",my_grp,"_keywd_",mychoice,"_all_q_only_TBZ_vs_SUC",q_threshold,".pdf", sep = ""), width = mywidth, height = myheight)
    
    # g1 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    #   geom_point(aes(size=`counts`)) + 
    #   #theme_bw() + 
    #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    #   labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    #   ggtitle(paste("Metagenome"))
    # 
    mydatabubble_for_metagenome_for_plot <- mydatabubble_for_metagenome[complete.cases(mydatabubble_for_metagenome$TBZperSUC.log2FC),]
    g2 <- ggplot(subset(mydatabubble_for_metagenome_for_plot), aes(x=bin,y=`functional group`)) +
      #geom_point(aes(size=`TBZperSUC.log2FC`, colour = `TBZperSUC.log2FCcls`)) + 
      geom_point(aes(size=`TBZperSUC.log2FC`, colour = I(TBZperSUC.log2FCcls))) + 
      #theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
      # aes(fill = TBZperSUC.log2FCcls) +
      # scale_fill_manual(values = TBZperSUC.log2FCclsRGB) +
      ggtitle(paste("TBZ/SUC gene expression"))
    
    # g3 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    #   geom_point(aes(size=`TBZ73perTBZ57.TBZ73.log2FC`, colour = I(TBZ73perTBZ57.TBZ73.log2FCcls))) + 
    #   #theme_bw() + 
    #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    #   labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    #   ggtitle(paste("TBZ73/TBZ57 gene expression"))
    # 
    # g4 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    #   geom_point(aes(size=`TBZ109perTBZ57.TBZ109.log2FC`, colour = I(TBZ109perTBZ57.TBZ109.log2FCcls))) + 
    #   #theme_bw() + 
    #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    #   labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    #   ggtitle(paste("TBZ109/TBZ57 gene expression"))
    # 
    # g5 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    #   geom_point(aes(size=`SUC57perSUC109.log2FC`, colour = I(SUC57perSUC109.log2FCcls))) + 
    #   #theme_bw() + 
    #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    #   labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    #   ggtitle(paste("SUC57/SUC109 gene expression"))
    
    
    #multiplot(g1, g2, g3, g4, g5, ncol=5)
    print(g2)
    
    
    
    dev.off()
    
    #}
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
### perform the network analysis ----
  
  library(igraph)
  library(vegan)
  
  

  #set the relative expression cutoff values for the features to be considered
  cutoffs<-c(0.0001,0.00001,0.001)
  # set the BH cutoff
  myBH <- 0.01
  # set the minimum gene_numbers_for_plotting_bins
  mymingenes <- 100

  
  for(cutoff in cutoffs){
    
    if(cutoff == 0.0001){
      mycorrcutoffs <- c(0.5,0.7)
    } else {
      mycorrcutoffs <- c(0.5,0.7)
    }
    for(mycorrcutoff in mycorrcutoffs){
      
      system(paste("mkdir -p network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,sep=""))
      #for(corr_test in c("spearman","pearson")){
      for(corr_test in c("spearman")){
        system(paste("mkdir -p network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,sep=""))
        
        
        
        ## Prepare the tables (check if the orientation is correct)
        
        final_tab <- t(my_cnt_annot_tab_wth_sts[,grep("\\.y",colnames(my_cnt_annot_tab_wth_sts))])[,complete.cases(my_cnt_annot_tab_wth_sts[,grep("\\.y",colnames(my_cnt_annot_tab_wth_sts))])]
        
        
        library(vegan)
        final_tab_ra<-decostand(final_tab,"total")
        #select features with equal or higher relative abundance to the cutoff in at least 2 samples and prepare a table where all rare features are considered as one
        final_tab_abund<-final_tab_ra[,apply(final_tab_ra,2,function(x) any(length(which(x>=cutoff))>=2))]
        final_tab_rare<-rowSums(final_tab_ra[,apply(final_tab_ra,2,function(x) any(length(which(x>=cutoff))<2))])
        final_mat<-decostand(data.frame(final_tab_abund,rare_features=final_tab_rare),"total")
        final_tab_with_rare<-data.frame(final_tab[,colnames(final_tab_abund)],rare_features=rowSums(final_tab[,colnames(final_tab_ra[,apply(final_tab_ra,2,function(x) any(length(which(x>=cutoff))<2))])]))
        
        
        
        #find mean site similarity
        dist.jaccard<-vegdist(final_mat,method="jaccard")
        mean.dist.jaccard<-mean(dist.jaccard)
        stdev.dist.jaccard<-sd(dist.jaccard)
        mean.sim.jaccard<- 1 - mean.dist.jaccard
        
        ## Perform the Spearman correlation test using the raw values as suggested in Berry and Widder (2014)
        library(Hmisc)
        correl<-rcorr(as.matrix(final_tab_with_rare), type = corr_test)
        r_vals<-correl$r
        p_vals<-correl$P
        diag(p_vals)<-0
        
        # adjust the p values using the BH method
        p_vals_arr <- array(p_vals)
        p_vals_arr_adj <- p.adjust(p_vals_arr, method = "BH")
        p_vals_arr_adj_mat <- matrix(p_vals_arr_adj, nrow = nrow(p_vals), ncol = ncol(p_vals))
        colnames(p_vals_arr_adj_mat) <- row.names(p_vals_arr_adj_mat) <- row.names(p_vals)
        
        # adjust the r values in order to remove (weak negative and weak positive correlations and) correlations where the adjusted p values are higher than 0.05
        r_vals_padj<-r_vals
        r_vals_padj <- ifelse(abs(r_vals_padj) < mycorrcutoff,0,r_vals)
        Sys.time()
        r_vals_padj_mat <- ifelse(p_vals_arr_adj_mat > myBH,0,r_vals_padj)
        Sys.time()
        
        library(igraph)
        g<-graph.adjacency(r_vals_padj_mat, weighted=TRUE, diag=FALSE
                           , mode="undirected"
                           #, mode="directed"
        )
        
        ##Add the bin and seed annotation as labels of the network vertices
        featurename_replacements4 <- my_cnt_annot_tab_wth_sts[colnames(final_mat),grep("BinTaxID|SEED", colnames(my_cnt_annot_tab_wth_sts))]
        ## create the colors
        # select the bin taxonomy info (this choice can be moved above as well)
        sel_tax_lev<-"BinTaxID"
        # create the color vector out of my beloved color collection
        
        V(g)$color1<-as.numeric(as.factor(featurename_replacements4[,sel_tax_lev]))        
        # add the taxon names as a vertice attribute which will be carried further on in subsequent analysis
        V(g)$color.labels<-as.character(featurename_replacements4[,sel_tax_lev])
        
        library(PerformanceAnalytics) # for the rainbowXequal
        library(RColorBrewer) # for the brewer.pal
        V(g)$color<-c(rainbow10equal,brewer.pal(8, name="Dark2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Set3"))[as.numeric(V(g)$color1)]
        
        
        #### get the keystoneness information
        keyg<-simplify(g, remove.loops=T)
        keygdegree<-degree(keyg)
        # get the indirect connections degree
        key2degree<-keygdegree
        for(connectedfeature in names(keygdegree)){
          key2degree[connectedfeature]<-sum(degree(keyg,neighbors(keyg, connectedfeature)))
        }
        
        # copy graph and keep edges with positive only weights for the betweeness and closeness centrality
        keyg_2 <- delete.edges(keyg, which(E(keyg)$weight <=0))
        betcentr<-betweenness(keyg_2)
        clocentr<-closeness(keyg_2,mode="all")
        tranc<-transitivity(keyg,type="local")
        # get the subpopulations according to density and add it to the final table
        gg <- keyg_2
        gg.com <- fastgreedy.community(keyg_2)
        mbrshp <- V(gg)$color.memb <- gg.com$membership
        
        # do the same consireding also negative weights
        gg2 <- keyg
        #gg.com <- fastgreedy.community(keyg)
        mbrshp2 <- V(gg)$color.memb <- gg.com$membership
        
        # compile the table
        df_net<-data.frame(degree=keygdegree,`indirect degree`=key2degree,`betweenness centrality`=betcentr,`closeness centrality`=clocentr,transitivity=tranc,`subcommunity membership`=mbrshp)
        
        ### add the network information to the master-table and save the table
        #my_cnt_annot_tab_wth_sts_ntwrk <- merge(my_cnt_annot_tab_wth_sts_ntwrk,df_net, all.x = T)
        #next
        
        # check the subnetworks where each of the bins is participating
        
        mybinscol <- unique(featurename_replacements4$BinTaxID)[!is.na(unique(featurename_replacements4$BinTaxID))]
        #### per bin heatmaps within network ----
        for(binid in mybinscol){
          mymemberships <- unique(df_net[row.names(featurename_replacements4[featurename_replacements4$BinTaxID == binid,]),6])
          mymemberships <- mymemberships[complete.cases(mymemberships)]
          myassociatedgenes <- row.names(df_net)[which(df_net$subcommunity.membership %in% mymemberships)]
          
          
          final_df_net<-merge(featurename_replacements4,df_net,by="row.names",all.y =T)
          
          #get the sphingomonas specific genes and memberships and extract the table from the master table
          sphing_specific <- final_df_net[which(final_df_net$Row.names%in%myassociatedgenes),]
          sphing_spec_master_tbl <- merge(my_cnt_annot_tab_wth_sts,sphing_specific, by.x = "ID",by.y = "Row.names", all.y = T)
          write.table(sphing_spec_master_tbl,paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/0_",binid,"_sig_expr_gene_membership_BH_", myBH,"_mingenes_",mymingenes,".txt",sep=""),row.names=F,quote=F,sep="\t")
          
          if(length(grep(binid,myassociatedgenes)) <= mymingenes){
            next # go to next bin in case that less than mingenes genes of the main bin participate in the bin with significant gene expression
          }
          
          ## prepare also heatmaps depicting the correlations between the selected bin and the rest membership bins ----
          
          if(length(unique(sphing_specific$BinTaxID)) > 1){
            mybinidgenes <- sphing_specific$Row.names[grep(binid,sphing_specific$BinTaxID)]
            myrestbins <- unique(sphing_specific$BinTaxID)[-which(unique(sphing_specific$BinTaxID) == binid)]
            myrestsbinidgenes <- sphing_specific$Row.names[grep(paste(myrestbins,collapse = "|"),sphing_specific$BinTaxID)]
            myselbin_tab <- r_vals_padj[mybinidgenes,myrestsbinidgenes]
            #set the 0 values to NA cause otherwise they will participate in the means
            myselbin_tab[which(myselbin_tab==0)] <- NA
            
            # get the grouping variables for the row names and aggregate 
            myselbin_tabrowgrps <- sphing_spec_master_tbl$SEED3.y[which(sphing_spec_master_tbl$ID%in%row.names(myselbin_tab))]
            
            myselbin_tab_agg1 <- aggregate(myselbin_tab, list(myselbin_tabrowgrps), mean)
            row.names(myselbin_tab_agg1) <- myselbin_tab_agg1$Group.1
            myselbin_tab_agg1 <- myselbin_tab_agg1[,2:ncol(myselbin_tab_agg1)]
            
            # do the same for the column names (this is slightly more complicated due to the potential multiple bins)
            myselbin_tab_agg1T <- t(myselbin_tab_agg1)
            
            mySEED3annot <- as.character(sphing_spec_master_tbl$SEED3.y[which(sphing_spec_master_tbl$ID%in%colnames(myselbin_tab))])
            mySEED3annot[is.na(mySEED3annot)] <- "other"
            
            myselbin_tabrowgrps2 <- paste(as.character(sphing_spec_master_tbl$BinTaxID.x[which(sphing_spec_master_tbl$ID%in%colnames(myselbin_tab))]),mySEED3annot, sep = " ")
            
            myselbin_tab_agg1T2 <- aggregate(myselbin_tab_agg1T, list(myselbin_tabrowgrps2), mean)
            row.names(myselbin_tab_agg1T2) <- myselbin_tab_agg1T2$Group.1
            myselbin_tab_agg1T2 <- myselbin_tab_agg1T2[,2:ncol(myselbin_tab_agg1T2)]
            
            # the final mat
            myselbin_tab_agg_fin <- t(myselbin_tab_agg1T2)
            myselbin_tab_agg_fin2 <- myselbin_tab_agg_fin
            #set the NAs back to zeros cause they are screwing the script
            myselbin_tab_agg_fin2[is.na(myselbin_tab_agg_fin)] <- 0
            
            # prepare the heatmap
            
            library(gplots)
            tryCatch({ 
              cairo_pdf(file = paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/0_",binid,"_sig_expr_gene_membership_heatmap_BH_", myBH,"_mingenes_",mymingenes,".pdf",sep=""),height = dim(myselbin_tab_agg_fin2)[1]/5, width = dim(myselbin_tab_agg_fin)[2]/2.5)
              # set the breaks and colours
              mybreaks <- seq(-1,1, length.out = 120)
              mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[6:1])(n = 119) # the nature colors
              # open the graphics device, plot and close the device
              pheatmap(as.matrix(myselbin_tab_agg_fin2), color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
                       , cluster_cols = F
                       , treeheight_row = 60
                       , treeheight_col = 60
              )
              dev.off()
            }, error=function(e) if(is.null(dev.list()) == F) {dev.off()})
            ## do the same as above for the positively correlated functional bins/categories with at least 10% of the main bin functional categories
            myselbin_tab_agg_fin3 <- myselbin_tab_agg_fin2
            myselbin_tab_agg_fin3 <- ifelse(myselbin_tab_agg_fin2 > 0 , 1, 0)
            #create a column indexing array
            myselbin_good_cols <- colSums(myselbin_tab_agg_fin3)/nrow(myselbin_tab_agg_fin3)
            myselbin_tab_agg_fin4 <- myselbin_tab_agg_fin2[,which(myselbin_good_cols >= 0.1)]
            
            tryCatch({ 
              
              cairo_pdf(file = paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/0_",binid,"_sig_expr_gene_membership_heatmap_best10perc_positive_BH_", myBH,"_mingenes_",mymingenes,"2.pdf",sep=""),height = dim(myselbin_tab_agg_fin4)[1]/5, width = dim(myselbin_tab_agg_fin4)[2]/2.5)
              # set the breaks and colours
              mybreaks <- seq(0,1, length.out = 120)
              mycolors <- colorRampPalette(c("white",RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[3:1]))(n = 119) # the nature colors
              # open the graphics device, plot and close the device
              pheatmap(as.matrix(myselbin_tab_agg_fin4), color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
                       #, cluster_rows = F
                       , cluster_cols = F
                       , treeheight_row = 60
                       , treeheight_col = 60
              )
              dev.off()
              
              
              
              
            }, error=function(e) if(is.null(dev.list()) == F) {dev.off()})           
            ## do the same as above for the negatively correlated functional bins/categories with at least 10% of the main bin functional categories
            myselbin_tab_agg_fin3 <- myselbin_tab_agg_fin2
            myselbin_tab_agg_fin3 <- ifelse(myselbin_tab_agg_fin2 < 0 , 1, 0)
            #create a column indexing array
            myselbin_good_cols <- colSums(myselbin_tab_agg_fin3)/nrow(myselbin_tab_agg_fin3)
            myselbin_tab_agg_fin4 <- myselbin_tab_agg_fin2[,which(myselbin_good_cols >= 0.1)]
            tryCatch({ 
              
              cairo_pdf(file = paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/0_",binid,"_sig_expr_gene_membership_heatmap_best10perc_negative_BH_", myBH,"_mingenes_",mymingenes,".pdf",sep=""),height = dim(myselbin_tab_agg_fin4)[1]/5, width = dim(myselbin_tab_agg_fin4)[2])
              # set the breaks and colours
              mybreaks <- seq(-1,1, length.out = 120)
              mycolors <- colorRampPalette(c(RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[6:1]))(n = 119) # the nature colors
              # open the graphics device, plot and close the device
              pheatmap(as.matrix(myselbin_tab_agg_fin4), color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
                       #, cluster_rows = F
                       , cluster_cols = F
                       , treeheight_row = 60
                       , treeheight_col = 60
              )
              dev.off()
              
              
            }, error=function(e) if(is.null(dev.list()) == F) {dev.off()})
            
          }
          
          
        }
        
        
        
        # now save the table
        write.table(final_df_net,paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/Keystonness_BH_", myBH,"_mingenes_",mymingenes,".txt",sep=""),row.names=F,quote=F,sep="\t")
        
        
        
        
        
        
        
        
        ################################
        ##### perform the analysis of the keystoneness values
        ################################
        
        tryCatch({ # use the tryCatch function in cause in some cases (at high cutoff values) the following part failed
          
          
          cairo_pdf(file=paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/barplot_panels_keystone_indices_BH_", myBH,"_mingenes_",mymingenes,".pdf",sep=""),height=15,width=10)
          par(mfrow=c(3,2))
          
          for(keyst_index in c("degree","indirect.degree","betweenness.centrality","closeness.centrality","transitivity")){
            counts_tab2<-final_df_net[,keyst_index]
            
            ## prepare the groups table according to the produced DCA
            design_for_mnse_renamed<-final_df_net[,sel_tax_lev] # in order to match the names of the p_adj table
            counts_tab1<-data.frame(treatment=design_for_mnse_renamed,counts=counts_tab2)
            counts_tab<-counts_tab1
            counts_tab[which(is.nan(counts_tab[,2])),2]<-0
            
            
            # perform the anova tests and prepare the significance group letters
            library(multcompView)
            lett<-array(dim=c(ncol(counts_tab)-1,length(levels(factor(counts_tab$treatment)))))
            row.names(lett)<-colnames(counts_tab)[1:(ncol(counts_tab)-1)]
            colnames(lett)<-levels(factor(counts_tab$treatment))
            
            # replace the dashes with underscores to avoid issues with the aov and multcompletters
            counts_tab[,"treatment"] <- factor(counts_tab[,"treatment"])
            
            for(myname_i in 2){
              form<-paste(colnames(counts_tab)[myname_i],"~ treatment")
              if(summary(aov(as.formula(form),counts_tab))[[1]]$Pr[1]<=0.05){
                tuk_res<-TukeyHSD(aov(as.formula(form),counts_tab))
                lett<-multcompLetters2(as.formula(form),tuk_res$treatment[,"p adj"],counts_tab)$Letters[levels(counts_tab$treatment)]
              } else {
                lett<-rep("a",length(levels(counts_tab$treatment)))
              }
              
            }
            
            
            
            lett_plot=lett
            if (summary(lett[]=="a")[[3]]==0)
            {
              lett_plot[]=""
            }
            
            
            # prepare the plot
            library(gplots)
            shared=counts_tab[,2] # the sample names do not come out in the same order as in the design... this is why I had to reorder them.
            design_for_mnse=as.factor(design_for_mnse_renamed)  
            
            means_sh=matrix(nrow=length(levels(design_for_mnse)),ncol=1)
            row.names(means_sh)=levels(design_for_mnse)
            colnames(means_sh)="counts"
            sterr_sh=matrix(nrow=length(levels(design_for_mnse)),ncol=1)
            row.names(sterr_sh)=levels(design_for_mnse)
            colnames(sterr_sh)="counts"
            
            for (i_mat in 1:length(levels(design_for_mnse)))
            {
              for (j_mat in 1)
              {
                means_sh[i_mat,j_mat]=mean(shared[which(design_for_mnse==levels(design_for_mnse)[i_mat])])
                sterr_sh[i_mat,j_mat]=sd(shared[which(design_for_mnse==levels(design_for_mnse)[i_mat])])/sqrt(length(design_for_mnse[which(design_for_mnse==levels(design_for_mnse)[i_mat])]))
              }
            }
            
            library(gdata)
            mat=t(interleave(means_sh,sterr_sh))
            mynames<-as.character(interleave(levels(design_for_mnse),levels(design_for_mnse)))
            colnames(mat)=paste(mynames,c(""," se"),sep="")
            
            
            mat_fin=mat
            
            # prepare the plot
            
            plot_mat=as.data.frame(t(mat_fin[,as.integer(seq(1,ncol(mat_fin),by=2))]))
            plot_err=as.data.frame(t(mat_fin[,as.integer(seq(2,ncol(mat_fin),by=2))]))
            plot_err[which(is.na(plot_err))]<-0
            write.table(plot_mat,paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/barplot_BH_", myBH,"_mingenes_",mymingenes,".txt",sep=""),quote=F,col.names=NA,sep="\t")
            write.table(plot_err,paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/barplot_err_BH_", myBH,"_mingenes_",mymingenes,".txt",sep=""),quote=F,col.names=NA,sep="\t")
            
            plot_mat<-as.matrix(read.table(paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/barplot_BH_", myBH,"_mingenes_",mymingenes,".txt",sep=""),header=T,check.names=F,quote="",row.names=1,sep="\t"))
            
            plot_err<-as.matrix(read.table(paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/barplot_err_BH_", myBH,"_mingenes_",mymingenes,".txt",sep=""),header=T,check.names=F,quote="",row.names=1,sep="\t"))
            
            library(PerformanceAnalytics) # for the rainbowXequal
            library(RColorBrewer) # for the brewer.pal
            par(mar=c(5,20,2,4),xpd=T)
            library(gplots)
            barplot2(plot_mat[nrow(plot_mat):1,],
                     horiz=T,
                     beside=T,
                     las=1,
                     cex.names=1.3,
                     col=c(rainbow10equal,brewer.pal(8, name="Dark2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Set3"))[1:ncol(plot_mat)],
                     plot.ci = TRUE,
                     ci.u=plot_mat[nrow(plot_mat):1,]+plot_err[nrow(plot_err):1,],
                     ci.l=plot_mat[nrow(plot_mat):1,]-plot_err[nrow(plot_err):1,],
                     xlab=gsub("\\."," ",keyst_index),
                     cex.lab=2,
                     font.lab=4,
                     xaxt="n",
                     cex.axis=0.81,
                     font.axis=4
            )
            
            max_axis=max(plot_mat+plot_err)
            axis(1
                 ,cex.axis=1.5)
            
            xvals=1*(plot_mat+plot_err)
            labs=lett_plot
            valsy=seq(1,1.2*(length(plot_mat)),by=1.2)-0.25
            yvals=valsy 
            # the orders of the values were inversed on purpose
            text(x=xvals,y=yvals,labels=labs,cex=1.2,font=2,pos=4)
            
            
          }
          dev.off()
          
        }, error=function(e) if(is.null(dev.list()) == F) {dev.off()})
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        # get the histogram of the communities according to participation within the graph
        cairo_pdf(paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/Membership_hist_BH_", myBH,"_mingenes_",mymingenes,".pdf",sep=""),height=4,width=5)
        doubletons_and_above<-V(gg)$color.memb[which(V(gg)$color.memb%in%which(table(V(gg)$color.memb) > 2))]
        hist(doubletons_and_above,breaks=0:(max(doubletons_and_above))+0.5, col="grey80",xlab="Sub-community #",ylab="gene abund.", main="gene Membership",font.lab=2)
        dev.off()
        
        relabund=t(final_mat)
        
        # first community (delete all the rest vertices)
        for(membership in 1:max(doubletons_and_above)){
          gg.one <- delete.vertices(gg, V(gg)[V(gg)$color.memb != membership ])
          #and simplify the graph for plotting by removing the loops
          gg.one <- simplify(gg.one, remove.loops=T)
          
          edge_width=E(gg.one)$weight
          #create a log transformed range of the values to be used for the plot vertice sizes
          vertsize1<-log10(apply((relabund[row.names(relabund)%in%V(gg.one)$name,])+0.000001,1,mean))
          vertsize2<-abs(min(vertsize1))+vertsize1
          vertsize3<-vertsize2/max(vertsize2)
          vertsize4<-vertsize3*20
          
          # calculate layout in advance for better debugging
          mylayout <- layout_with_fr(gg.one)
          #mylayout <- layout_with_dh(gg.one)
          tryCatch({ # use the tryCatch function in cause in some cases (at high cutoff values) the following part failed
            cairo_pdf(file=paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/Membership_",membership,"_network_BH_", myBH,"_mingenes_",mymingenes,".pdf",sep=""),width=10,height=7)
            #I moved the legend to a different screen to avoid overlaps
            split.screen(rbind(c(0.1,0.3,0.1, 0.98), c(0.5, 0.95, 0.1, 0.98)))
            screen(1)
            par(mar=c(1,0.5,1,0),xpd=T)
            
            plot.new()
            legend("topleft",gsub("Bacteria","unclassified Bacteria",gsub(" Chloroplast","",gsub("_"," ",levels(as.factor(V(g)$color.labels))[which(levels(as.factor(V(g)$color.labels))%in%levels(as.factor(V(gg.one)$color.labels)))]))), pch=21, pt.bg=c(rainbow10equal,brewer.pal(8, name="Dark2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Set3"))[1:length(levels(as.factor(V(g)$color.labels)))][which(levels(as.factor(V(g)$color.labels))%in%levels(as.factor(V(gg.one)$color.labels)))], pt.cex=2, cex=.8, bty="n", ncol=1)
            text(-0.05,0.05,labels=paste("sub-community ",membership,sep=""),font=2,cex=1.2,pos=4)
            text(-0.05,0,labels=paste("richness: ",length(unique(V(gg.one)$color.labels)),sep=""),font=2,cex=1.2,pos=4)
            
            
            #plot
            library(PerformanceAnalytics)
            
            screen(2)
            par(mar=c(1,1,1,1),xpd=T)
            plot.igraph(gg.one
                        # the layout was pre-calculated (see mylayout above and also"?layout_"... yes!!! use the underscore)
                        ,layout=mylayout
                        ,edge.width=(2*E(gg.one)$weight)
                        ,vertex.size=vertsize4
                        ,vertex.color=V(gg.one)$color
                        #,vertex.label=gsub("feature_","",names(vertsize4))
                        #,vertex.label.dist=4
                        #,vertex.label=1:length(V(gg.one))
                        ,vertex.label.color="black"
                        ,vertex.label.cex=0.7
                        ,edge.curved=F
                        #edge.label=round(E(gg.one)$weight-1,digits=2),
                        #edge.label.color="grey30",
                        #edge.label.cex=1.2,
                        ,edge.label.font=2
                        ,edge.color="grey70"
            )
            
            dev.off()
            
          }, error=function(e) if(is.null(dev.list()) == F) {dev.off()})
          
          
          
          
          
          tryCatch({ # use the tryCatch function in cause in some cases (at high cutoff values) the following part failed
            #### same plot as above without the feature names
            cairo_pdf(file=paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/Membership_",membership,"_network_no_feature_names_BH_", myBH,"_mingenes_",mymingenes,".pdf",sep=""),width=10,height=7)
            #i moved the legend to a different screen to avoid overlaps
            split.screen(rbind(c(0.1,0.3,0.1, 0.98), c(0.5, 0.95, 0.1, 0.98)))
            screen(1)
            par(mar=c(1,0,1,0),xpd=T)
            
            plot.new()
            legend("topleft",gsub("Bacteria","unclassified Bacteria",gsub(" Chloroplast","",gsub("_"," ",levels(as.factor(V(g)$color.labels))[which(levels(as.factor(V(g)$color.labels))%in%levels(as.factor(V(gg.one)$color.labels)))]))), pch=21, pt.bg=c(rainbow10equal,brewer.pal(8, name="Dark2"))[1:length(levels(as.factor(V(g)$color.labels)))][which(levels(as.factor(V(g)$color.labels))%in%levels(as.factor(V(gg.one)$color.labels)))], pt.cex=2, cex=.8, bty="n", ncol=1)
            text(-0.05,0.05,labels=paste("sub-community ",membership,sep=""),font=2,cex=1.2,pos=4)
            text(-0.05,0,labels=paste("richness: ",length(unique(V(gg.one)$color.labels)),sep=""),font=2,cex=1.2,pos=4)
            
            
            #plot
            library(PerformanceAnalytics)
            
            screen(2)
            par(mar=c(1,1,1,1),xpd=T)
            plot.igraph(gg.one
                        # the layout was pre-calculated (see mylayout above and also"?layout_"... yes!!! use the underscore)
                        ,layout=mylayout
                        ,edge.width=(2*E(gg.one)$weight)
                        ,vertex.size=vertsize4
                        ,vertex.color=V(gg.one)$color
                        ,vertex.label=""
                        #,vertex.label.dist=4
                        #,vertex.label=1:length(V(gg.one))
                        ,vertex.label.color="black"
                        ,vertex.label.cex=0.7
                        ,edge.curved=F
                        #edge.label=round(E(gg.one)$weight-1,digits=2),
                        #edge.label.color="grey30",
                        #edge.label.cex=1.2,
                        ,edge.label.font=2
                        ,edge.color="grey70"
            )
            
            dev.off()
            
          }, error=function(e) if(is.null(dev.list()) == F) {dev.off()})
          
          
          ### write the graph
          write.graph(gg.one, file=paste("network/corrcutoff_", mycorrcutoff ,"_relabund_cutoff_",cutoff,"/",corr_test,"/Membership_",membership,"_network_BH_", myBH,"_mingenes_",mymingenes,".graphml",sep=""), format="graphml")
        }
      }
    }
  }
  
    
# save the universal table including the network memberships
#write.table(my_cnt_annot_tab_wth_sts_ntwrk, "000_universal_table.txt", sep = "\t", col.names = NA)
  
  
#### prepare the cobalamin biosynthetic pathway genes expression plot ----
  # this section requires running in advance the bash scripts found in the 1_metagenome/z_accessory_files_dbs_and_executbles/cobalamin_biosynthesis folder that searches the manually curated UniProt SWISS-prot database
  myb12biohits <- read.table("../../1_metagenome/z_accessory_files_dbs_and_executbles/cobalamin_biosynthesis/uniprot_sprot.fasta.dmnd.m8_1e-10_qcov50", sep = "\t", stringsAsFactors = F, quote = "", comment.char = "")
  colnames(myb12biohits) <- c("qseqid_B12","stitle_B12","pident_B12","length_B12","mismatch_B12","gapopen_B12","qstart_B12","qend_B12","sstart_B12","send_B12","evalue_B12","bitscore_B12")
  
  # do some cleanup
  myb12biohits$ID <- gsub("\\|.+","",myb12biohits$qseqid_B12)
  myb12biohits$gene_hits_B12 <- gsub(".+ cobS .+","",gsub(".+ cobO .+","cobO",gsub(".+ BtuF .+","BtuF",gsub(".+ CobB\\/CobQ-like.+","CobB",gsub(" PE=.+","",gsub(".+ GN=","",gsub(".+ GN= ","",myb12biohits$stitle))))))) 
  myb12biohits$gene_hits_B12[grep("Uncharacterized protein",myb12biohits$gene_hits_B12)] <- NA
  
  myb12biohits <- myb12biohits[complete.cases(myb12biohits$gene_hits_B12),]
  
  # prep the directory for storring the results
  system("mkdir cobalamin_path_expr")
  compl.tbl <- merge(my_cnt_annot_tab_wth_sts,myb12biohits, by = "ID", all = T)
  # merge the old gene annotation and replace with the gene_hits_B12 only if there is no other annotation in place
  compl.tbl$gene <- gsub("_[0-9]+$","",compl.tbl$gene)
  compl.tbl$my_new_gene <- compl.tbl$gene
  for(myrownum in 1:nrow(compl.tbl)){
    if(is.na(compl.tbl$my_new_gene[myrownum])){
      compl.tbl$my_new_gene[myrownum] <- compl.tbl$gene_hits_B12[myrownum]
    }
  }
  design <- design_tab_sel
  # set the names of the genes of interest
  mygenes <- c("cobA","cobIJ","cobI","cobJ","cobG","cobM","cobF","cobK","cobL","cobH","cobBC","cobN","cobS","cobT","cobO","cobQ","cobD","cobP","cobU","cobV","btuC","btuF","btuD","btuB")
  
  ## bubble plots
  qcutoff <- 0.01
  myannfrb120 <- compl.tbl[grep(paste(mygenes, sep = "", collapse = "|"), compl.tbl$my_new_gene),c(which(colnames(compl.tbl)%in%c("ID","my_new_gene","TBZ_vs_SUC_logFC","TBZ_vs_SUC_q.value","BinTaxID")),grep("\\.y",colnames(compl.tbl)))]
  myannfrb120$genecps <- rep(1, nrow(myannfrb120))
  myannfrb120$CPM.TBZ <- rowMeans(myannfrb120[,grep(paste(design$sampleName[which(design$treatment == "TBZ")], collapse = "|"),colnames(myannfrb120))])
  myannfrb120$CPM.SUC <- rowMeans(myannfrb120[,grep(paste(design$sampleName[which(design$treatment == "SUC")], collapse = "|"),colnames(myannfrb120))])
  myannfrb12b <- aggregate(myannfrb120[,"genecps"], by = list(myannfrb120$my_new_gene, myannfrb120$BinTaxID), sum, na.rm=T)
  row.names(myannfrb12b) <- paste(myannfrb12b$Group.1,myannfrb12b$Group.2)
  myannfrb12c <- aggregate(myannfrb120[,grep("TBZ_vs_SUC",colnames(myannfrb120))], by = list(myannfrb120$my_new_gene, myannfrb120$BinTaxID), mean, na.rm=T)
  row.names(myannfrb12c) <- paste(myannfrb12c$Group.1,myannfrb12c$Group.2)
  myannfrb12d <- merge(myannfrb12b, myannfrb12c, by = "row.names")
  myannfrb12d$Group.1.x <- factor(myannfrb12d$Group.1.x, levels = mygenes[length(mygenes):1])
  colnames(myannfrb12d)[2:4] <- c("gene","bin","genecnts")
  myannfrb12d$clr <- "black"
  myannfrb12d$clr[which(myannfrb12d$TBZ_vs_SUC_q.value > qcutoff)] <- NA
  myannfrb12d$clr[which(myannfrb12d$TBZ_vs_SUC_logFC<0)] <- "red"
  myannfrb12d$clr[which(myannfrb12d$TBZ_vs_SUC_logFC>0)] <- "blue"
  
  myannfrb12d$TBZ_vs_SUC_logFC <- abs(myannfrb12d$TBZ_vs_SUC_logFC)
  
  library(ggplot2)
  library(ggtree)
  
  cairo_pdf(paste("cobalamin_path_expr/bubble_plotnopar_1e-10_qcov50.pdf", sep = ""), width = 10, height = 6)
  theme_set(theme_bw())
  g1 <- ggplot(subset(myannfrb12d), aes(x=bin,y=gene)) +
    geom_point(aes(size=`genecnts`, colour = I(clr))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("gene ",sep = "")) +
    ggtitle(paste("Gene counts"))
  
  g1$labels$size <- "gene cnts"
  g1$data$clr[grep("groopm_bin_50_Sphingomonas",g1$data$bin)] <- "black"
  g1$data$clr[-grep("groopm_bin_50_Sphingomonas",g1$data$bin)] <- "grey70"
  
  g2 <- ggplot(subset(myannfrb12d), aes(x=bin,y=gene)) +
    geom_point(aes(size=`TBZ_vs_SUC_logFC`, colour = I(clr))) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("gene ",sep = "")) +
    ggtitle(paste("TBZ/SUC gene expression"))
  
  g2$labels$size <- "log2FC"
  g3 <- g2
  for(mycolggplot in 1:length(g2$data$clr)){
    if(g2$data$clr[mycolggplot] == "blue" & g2$data$bin[mycolggplot] != "groopm_bin_50_Sphingomonas" & g2$data$bin[mycolggplot] != "solidbin_5229_sub_Hydrogenophaga"){
      g3$data$clr[mycolggplot] <- rgb(0, 0, 255, max = 255, alpha = 50)
    } else if(g2$data$clr[mycolggplot] == "blue" & g2$data$bin[mycolggplot] %in% c("groopm_bin_50_Sphingomonas","solidbin_5229_sub_Hydrogenophaga")){
      g3$data$clr[mycolggplot] <- rgb(0, 0, 255, max = 255, alpha = 255)
    } else if(g2$data$clr[mycolggplot] == "red" & g2$data$bin[mycolggplot] != "groopm_bin_50_Sphingomonas" & g2$data$bin[mycolggplot] != "solidbin_5229_sub_Hydrogenophaga"){
      g3$data$clr[mycolggplot] <- rgb(255, 0, 0, max = 255, alpha = 50)
    } else if(g2$data$clr[mycolggplot] == "red" & g2$data$bin[mycolggplot] %in% c("groopm_bin_50_Sphingomonas","solidbin_5229_sub_Hydrogenophaga")){
      g3$data$clr[mycolggplot] <- rgb(255, 0, 0, max = 255, alpha = 255)
    }
  }

  multiplot(g1, g3, ncol=2)
  
  dev.off()
  
  



    
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### prepare the catechol degradation pathway genes expression plot ----
  # this section exploits the files and scripts run in advance found in the 1_metagenome/z_accessory_files_dbs_and_executbles/cobalamin_biosynthesis folder that searches the manually curated UniProt SWISS-prot database... essentially I exploit the initial material of the previous section. Due to the fact that the need for this came late, I did not have the time to screen and change all necessary file-names to something more meaningful for catechol degradation.
  myb12biohits <- read.table("../../1_metagenome/z_accessory_files_dbs_and_executbles/cobalamin_biosynthesis/uniprot_sprot.fasta.dmnd.m8_1e-10_qcov50", sep = "\t", stringsAsFactors = F, quote = "", comment.char = "")
  
  colnames(myb12biohits) <- c("qseqid_B12","stitle_B12","pident_B12","length_B12","mismatch_B12","gapopen_B12","qstart_B12","qend_B12","sstart_B12","send_B12","evalue_B12","bitscore_B12")
  
  # do some cleanup (I kept the original cobalamin associated cleanup in order to avoid extra work... the cat specific follows)
  myb12biohits$ID <- gsub("\\|.+","",myb12biohits$qseqid_B12)
  myb12biohits$gene_hits_B12 <- gsub(".+ cobS .+","",gsub(".+ cobO .+","cobO",gsub(".+ BtuF .+","BtuF",gsub(".+ CobB\\/CobQ-like.+","CobB",gsub(" PE=.+","",gsub(".+ GN=","",gsub(".+ GN= ","",myb12biohits$stitle))))))) 
  myb12biohits$gene_hits_B12[grep("Uncharacterized protein",myb12biohits$gene_hits_B12)] <- NA
  
  myb12biohits$gene_hits_B12 <- gsub("pcaR","catR",gsub("cat[1,3]","catR",gsub("cat1 operon transcr.+","catR",gsub(".+Probable ","",gsub("catRI+","catR",gsub("pcaJ","catJ",gsub("pcaF","catF",gsub("pcaI","catI",gsub("pcaC","catC",gsub("pcaB","catB",gsub("catC[0-9]","catC",gsub("catB[0-9]","catB",gsub("catA[0-9]","catA",myb12biohits$gene_hits_B12))))))))))))) 
  
  
  myb12biohits <- myb12biohits[complete.cases(myb12biohits$gene_hits_B12),]
  
  # prep the directory for storring the results
  system("mkdir catechol_path_expr")
  compl.tbl <- merge(my_cnt_annot_tab_wth_sts,myb12biohits, by = "ID", all = T)
  # merge the old gene annotation and replace with the gene_hits_B12 only if there is no other annotation in place
  compl.tbl$gene <- gsub("_[0-9]+$","",compl.tbl$gene)
  compl.tbl$my_new_gene <- compl.tbl$gene
  for(myrownum in 1:nrow(compl.tbl)){
    if(is.na(compl.tbl$my_new_gene[myrownum])){
      compl.tbl$my_new_gene[myrownum] <- compl.tbl$gene_hits_B12[myrownum]
    }
  }
  design <- design_tab_sel
  # set the names of the genes of interest
  mygenes <- c("catD","catA","catC","catB","catR","catI","catJ","catF","catE")
  
  ## bubble plots
  qcutoff <- 0.01
  compl.tbl$my_new_gene <- gsub("catRI+","catR",compl.tbl$my_new_gene)
  myannfrb120 <- compl.tbl[grep(paste(mygenes, sep = "", collapse = "|"), compl.tbl$my_new_gene),c(which(colnames(compl.tbl)%in%c("ID","my_new_gene","TBZ_vs_SUC_logFC","TBZ_vs_SUC_q.value","BinTaxID")),grep("\\.y",colnames(compl.tbl)))]
  myannfrb120$genecps <- rep(1, nrow(myannfrb120))
  myannfrb120$CPM.TBZ <- rowMeans(myannfrb120[,grep(paste(design$sampleName[which(design$treatment == "TBZ")], collapse = "|"),colnames(myannfrb120))])
  myannfrb120$CPM.SUC <- rowMeans(myannfrb120[,grep(paste(design$sampleName[which(design$treatment == "SUC")], collapse = "|"),colnames(myannfrb120))])
  myannfrb12b <- aggregate(myannfrb120[,"genecps"], by = list(myannfrb120$my_new_gene, myannfrb120$BinTaxID), sum, na.rm=T)
  row.names(myannfrb12b) <- paste(myannfrb12b$Group.1,myannfrb12b$Group.2)
  myannfrb12c <- aggregate(myannfrb120[,grep("TBZ_vs_SUC",colnames(myannfrb120))], by = list(myannfrb120$my_new_gene, myannfrb120$BinTaxID), mean, na.rm=T)
  row.names(myannfrb12c) <- paste(myannfrb12c$Group.1,myannfrb12c$Group.2)
  
  myannfrb12d <- merge(myannfrb12b, myannfrb12c, by = "row.names")
  myannfrb12d$Group.1.x <- factor(myannfrb12d$Group.1.x, levels = mygenes[length(mygenes):1])
  colnames(myannfrb12d)[2:4] <- c("gene","bin","genecnts")
  myannfrb12d$clr <- "black"
  myannfrb12d$clr[which(myannfrb12d$TBZ_vs_SUC_logFC < 0)] <- "red"
  myannfrb12d$clr[which(myannfrb12d$TBZ_vs_SUC_logFC > 0)] <- "blue"
  myannfrb12d$clr[which(myannfrb12d$TBZ_vs_SUC_q.value > qcutoff)] <- NA
  
  myannfrb12d$TBZ_vs_SUC_logFC <- abs(myannfrb12d$TBZ_vs_SUC_logFC)
  
  library(ggplot2)
  library(ggtree)
  
  cairo_pdf(paste("catechol_path_expr/bubble_plotnopar_1e-10_qcov50_FDRcutof_",qcutoff,".pdf", sep = ""), width = 10, height = 5)
  theme_set(theme_bw())
  g1 <- ggplot(subset(myannfrb12d), aes(x=bin,y=gene)) +
    geom_point(aes(size=`genecnts`, colour = I(clr))) + 
    #theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("gene ",sep = "")) +
    ggtitle(paste("Gene counts"))
  
  g1$labels$size <- "gene cnts"
  g1$data$clr[grep("groopm_bin_50_Sphingomonas",g1$data$bin)] <- "black"
  g1$data$clr[-grep("groopm_bin_50_Sphingomonas",g1$data$bin)] <- "grey70"
  
  g2 <- ggplot(subset(myannfrb12d), aes(x=bin,y=gene)) +
    geom_point(aes(size=`TBZ_vs_SUC_logFC`, colour = I(clr))) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("gene ",sep = "")) +
    ggtitle(paste("TBZ/SUC gene expression"))
  
  g2$labels$size <- "log2FC"
  g3 <- g2
  for(mycolggplot in 1:length(g2$data$clr)){
    if(is.na(g2$data$clr[mycolggplot])){
      next
    } else if(g2$data$clr[mycolggplot] == "blue" & g2$data$bin[mycolggplot] != "groopm_bin_50_Sphingomonas"){
      g3$data$clr[mycolggplot] <- rgb(0, 0, 255, max = 255, alpha = 50)
    } else if(g2$data$clr[mycolggplot] == "blue" & g2$data$bin[mycolggplot] %in% c("groopm_bin_50_Sphingomonas")){
      g3$data$clr[mycolggplot] <- rgb(0, 0, 255, max = 255, alpha = 255)
    } else if(g2$data$clr[mycolggplot] == "red" & g2$data$bin[mycolggplot] != "groopm_bin_50_Sphingomonas"){
      g3$data$clr[mycolggplot] <- rgb(255, 0, 0, max = 255, alpha = 50)
    } else if(g2$data$clr[mycolggplot] == "red" & g2$data$bin[mycolggplot] %in% c("groopm_bin_50_Sphingomonas")){
      g3$data$clr[mycolggplot] <- rgb(255, 0, 0, max = 255, alpha = 255)
    }
  }
  
  multiplot(g1, g3, ncol=2)
  
  dev.off()
  
  