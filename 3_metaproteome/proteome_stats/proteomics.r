#### load and prepare the data ----
## prep the design objects
design_tab <- read.table("../design", header = T, sep = "\t", stringsAsFactors = F)
design_tab$time <- as.factor(design_tab$time)
design_tab_sel <- design_tab[which(design_tab$read_amount == "ok"),]
myint <- interaction(design_tab_sel[,2:3])

## load the table with all gene annotations and sequences
my_annot_tbl_pre1 <- read.table("../../1_metagenome/4_annotation/ANNOTATION_DIR/genome.tsv", header = T, quote = "", comment.char = "", sep = "\t", stringsAsFactors = F)

row.names(my_annot_tbl_pre1) <- my_annot_tbl_pre1$ID
# create also the gc content column
library(stringr)
GCcounts <- str_count(my_annot_tbl_pre1$myfna, "G|C")
Allcounts <- str_count(my_annot_tbl_pre1$myfna, "A|T|G|C")
my_annot_tbl_pre1$gc_content <- GCcounts/Allcounts

# then, add the SEED annotation (prepare it by adding a header and also by indexing for merging)
my_seed_annot <- read.table("../../1_metagenome/4_annotation/SEED_annotation/myresults.m8_wth_ids.converted", header = F, quote = "", comment.char = "", sep = "\t", stringsAsFactors = F, fill = T)
colnames(my_seed_annot) <- c("qseqid","qseqid.1","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","SEED4","SEED3","SEED2","SEED1")
my_seed_annot$idxseed <- gsub("\\|.+","",my_seed_annot$qseqid)
#merge the tables
my_annot_tbl_pre2 <- merge(my_annot_tbl_pre1, my_seed_annot, by.x = "ID", by.y = "idxseed", all = T)

# then, the eggnog annotation (required some manual curation as well: the header had to be uncommented, while some field names had to be added)
my_eggnog_annot <- read.table("../../1_metagenome/4_annotation/eggnog_mapper_annotation/query_seqs.fa.emapper.annotations_franalysis", header = T, quote = "", comment.char = "#", sep = "\t", stringsAsFactors = F)
# prepare it for merging as well
my_eggnog_annot$idxeggnog <- gsub("\\|.+","",my_eggnog_annot$query_name)
#merge the tables
my_annot_tbl_pre3 <- merge(my_annot_tbl_pre2, my_eggnog_annot, by.x = "ID", by.y = "idxeggnog", all = T)




## load the protein relative_abundance data
my_prot_counts_pre <- read.table("../prot_relabun_mat", header = T, row.names = 1)
# use the following to merge the proteomics results (based on the previous assembly/annotation version) with the post nanopore assembly annotations... this leaves out 139 annotations of the first version
my_cnt_annot_tab <- merge(my_annot_tbl_pre3, my_prot_counts_pre, by.x = "old_idx", by.y = "row.names", all.x = T)
# order according to the new annotation index and rename the rows with the DFAST annotation IDs
my_cnt_annot_tab <- my_cnt_annot_tab[order(my_cnt_annot_tab$new_idx),]
row.names(my_cnt_annot_tab) <- my_cnt_annot_tab$ID
# also prepare the my_prot_counts for other uses
my_prot_counts <- my_cnt_annot_tab[complete.cases(my_cnt_annot_tab[,grep("SR",colnames(my_cnt_annot_tab))]),grep("SR",colnames(my_cnt_annot_tab))]

# remove NA rows according to the protein counts matrix
my_cnt_annot_tab_fin <- my_cnt_annot_tab[complete.cases(my_cnt_annot_tab[,grep("SR",colnames(my_cnt_annot_tab))]),]

# prepare the final counts table
cnts.tbl <- my_cnt_annot_tab_fin[,grep("SR",colnames(my_cnt_annot_tab_fin))]
# prepare the final design matrix based on the interactions of the grouping variables of treatment and time
design.grps <- interaction(design_tab_sel[,2:3])
design.grps <- factor(design.grps, levels = c("TBZ.57", "SUC.57", "TBZ.73", "SUC.73", "TBZ.109", "SUC.109"))

# given the pepetide survey method of the orbitrap-LC-MS (each event has equal chances of being analyzed), this equals to having genes of equal length during expression data analysis... for this reason I set tis parameter to 1 for all analyzed genes
my_cnt_annot_tab_sel_fin$gene_length <- 1
gene_lngts <- my_cnt_annot_tab_sel_fin$gene_length

#### run edgeR analysis ----
library(edgeR)
## prep the DGElist object
my_dgelist <- DGEList(cnts.tbl, group = design.grps, genes = gene_lngts)


##'As a rule of thumb, we require that a gene have a protein count of at least 10–15 in at least some libraries before it is considered to be expressed in the study. We could explicitly select for genes that have at least a couple of counts of 10 or more, but it is slightly better to base the filtering on count-per-million (CPM) values so as to avoid favoring genes that are expressed in larger libraries over those expressed in smaller libraries. For the current analysis, we keep genes that have CPM values above 5 (=> 5 copies per library of about 1 million reads) in at least two libraries

my_dgelist_good <- rowSums(cpm(my_dgelist) > 10) >= 5


my_dgelist_fin <- my_dgelist[my_dgelist_good, , keep.lib.sizes=FALSE]
my_dgelist_fin <- calcNormFactors(my_dgelist_fin)


#### canonical analysis ----
# prepare the plot of the canonical analysis of the normalized data in order to assess the importance of the treatment and time

library(vegan)
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
# DCA1    DCA2     DCA3     DCA4
# Eigenvalues     0.1870 0.01666 0.013457 0.005608
# Decorana values 0.1871 0.01433 0.007001 0.003091
# Axis lengths    1.0860 0.53979 0.481408 0.345629
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
df<-data.frame(sc.1,design_tab_sel[,2:3])
df$treatment <- factor(df$treatment,levels = c("TBZ","SUC"))
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
cairo_pdf("0_lfc_based_all_samples_MD_vs_average.pdf", height = 10, width = 10)
par(mfrow = c(5,4))
for(i in 1:ncol(my_dgelist_fin$counts)){
  plotMD(my_dgelist_fin, column=i, bty = "n", main = design.grps[i], values = c(-2,2), hl.col = c("blue","red"))
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()



## prepare the design/model matrix ----
mydesign <- model.matrix(~0+design.grps)



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


# prep the matrix to be amended with the rnaseq comparison results
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]

# add also the cpm for each gene
my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts,cpm(my_dgelist_fin), by = "row.names", all.x =T)
row.names(my_cnt_annot_tab_wth_sts) <- my_cnt_annot_tab_wth_sts$Row.names
my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[,2:ncol(my_cnt_annot_tab_wth_sts)]

# prep also the interaction/sample combination
design.grps_descr <- data.frame(samples = design_tab_sel$sample, grps = design.grps, stringsAsFactors = F)

# run the thing
cairo_pdf(paste("0_all_samples_DE_ME3.pdf", sep = ""), height = 12, width = 12)
par(mfrow = c(4,4))
for(i in 1:length(mycombs_names)){
  contr0 <- makeContrasts(contrasts = c(mycombs[i,2],mycombs[i,1]), levels=mydesign)
  mycomb <- data.frame(Contrasts = contr0[,1] - contr0[,2])

  # We will use QL F-tests instead of the more usual likelihood ratio tests (LRT) as they give stricter error rate control by accounting for the uncertainty in dispersion estimation: 
  res <- glmQLFTest(fit, contrast=mycomb)
  # get the fold change values etc and add them to the main table
  my_stats_tbl <- res$table
  my_stats_tbl$q.value <- p.adjust(my_stats_tbl$PValue,"fdr")
  colnames(my_stats_tbl) <- paste(gsub("design.grps","",paste(mycombs[i,1],"_vs_",mycombs[i,2], sep  ="")), colnames(my_stats_tbl), sep = "_")
  
  
  mynumerator <- gsub("design.grps","",mycombs[i,1])
  mydenominator <- gsub("design.grps","",mycombs[i,2])
  mynumersamples <- design.grps_descr[which(design.grps_descr$grps==mynumerator),1]
  mydenomsamples <- design.grps_descr[which(design.grps_descr$grps==mydenominator ),1]
  mysamples <- design.grps_descr[which(design.grps_descr$grps==mynumerator | design.grps_descr$grps==mydenominator ),1]
  
  # get the log2FC table and add the comparisons in the header
  mylogFC_tbl <- log2((cpm(my_dgelist_fin)[,mysamples] + 0.01)/rowMeans(cpm(my_dgelist_fin)[,mydenomsamples] + 0.01))
  colnames(mylogFC_tbl) <- paste(gsub("design.grps","",paste(mycombs[i,1],"_vs_",mycombs[i,2],"_log2FC", sep  ="")), colnames(mylogFC_tbl), sep = "_")
  
  myplottbl1 <- res$table[,c(1,4)]
  myplottbl2 <- merge(rowMeans(my_prot_counts), myplottbl1, by="row.names", all = T)
  myplottbl3 <- myplottbl2[,-which("Row.names"%in%colnames(myplottbl2))]
  row.names(myplottbl3) <- myplottbl2[,which("Row.names"%in%colnames(myplottbl2))]
  myplottbl4 <- myplottbl3[complete.cases(myplottbl3),]
  myplottbl4$col <- "grey70"
  # correct the PValue name to the adjusted pvalue name and the x axis into protein RA
  colnames(myplottbl4)[3] <- "FDR"
  colnames(myplottbl4)[1] <- "protein mean RA"
  # set the coloring conditions
  fdrcutoff <- 0.001
  lfccutoff <- 2
  myplottbl4$col[which(myplottbl4$logFC<= -lfccutoff & myplottbl4$FDR <= fdrcutoff)] <- "blue"
  myplottbl4$col[which(myplottbl4$logFC>= lfccutoff & myplottbl4$FDR <= fdrcutoff)] <- "red"
  plot(log10(myplottbl4$`protein mean RA`/10e6), myplottbl4$logFC, bty = "n", pch=16, cex = 0.3, main = gsub("design.grps","",paste(mycombs[i,2]," vs ",mycombs[i,1], sep  ="")), ylab = bquote("log"[2]*"FC"), xlab = bquote("log"[10]*"("~.(colnames(myplottbl4)[1])~")"), col = myplottbl4$col, ylim = c(-10,10))
  #abline(h = 0, lty = 2, col = "purple")
  abline(h = c(-lfccutoff,0,lfccutoff), lty = 2, col = "grey50")
  myxplotval <-mean(c(summary(log10(myplottbl4$`protein mean RA`/10e6))[5],summary(log10(myplottbl4$`protein mean RA`/10e6))[6]))
  text(myxplotval, -6, labels = paste(gsub("design.grps","",mycombs[i,1]), " (", length(myplottbl4$col[which(myplottbl4$logFC<= -lfccutoff & myplottbl4$FDR <= fdrcutoff)])," proteins)", sep = ""), col = "blue", cex = 1.2)
  text(myxplotval, 0, labels = paste("(", nrow(myplottbl4) - length(which(myplottbl4$logFC<= -lfccutoff & myplottbl4$FDR <= fdrcutoff)) - length(which(myplottbl4$logFC>= lfccutoff & myplottbl4$FDR <= fdrcutoff))," proteins)", sep = ""), col = "grey30", cex = 1.2)
  text(myxplotval, 6, labels = paste(gsub("design.grps","",mycombs[i,2]), " (", length(myplottbl4$col[which(myplottbl4$logFC >= lfccutoff & myplottbl4$FDR <= fdrcutoff)])," proteins)", sep = ""), col = "red", cex = 1.2)
  
  #dev.off()
}
dev.off()


# save the table with the sample/date contrasts

write.table(my_cnt_annot_tab_wth_sts, file = "00_my_final_tab_with_seqs_and_sampletime_contrasts.tsv", sep = "\t", quote = F, col.names = NA)

















#### TBZ vs SUC ----

contr0 <- model.matrix(~0+design_tab_sel$treatment)
colnames(contr0) <- c("SUC","TBZ")
con <- makeContrasts(TBZvsSUC = TBZ - SUC, levels = contr0)

fit1 <- glmQLFit(my_dgelist_fin, contr0, robust=TRUE)
res <- glmQLFTest(fit1, contrast = con)
# get the fold change values etc and add them to the main table
my_stats_tbl <- res$table
my_stats_tbl$q.value <- p.adjust(my_stats_tbl$PValue,"fdr")


#prepare for the log2FC matrix use the initial samples (day 57) for comparins 
mydenomsamples <- design_tab_sel[which(design_tab_sel$treatment== "SUC" ),1]
mysamples <- design_tab_sel[grep("SUC|TBZ",design_tab_sel$treatment),1]

# get the log2FC table and add the comparisons in the header
mylogFC_tbl <- log2((cpm(my_dgelist_fin)[,mysamples] + 0.01)/rowMeans(cpm(my_dgelist_fin)[,mydenomsamples] + 0.01))
colnames(mylogFC_tbl) <- paste(colnames(mylogFC_tbl), "_vs_SUC_mean_log2FC", sep = "")

## prepare a heatmap ----
my_annotab_for_heatmap <- my_cnt_annot_tab_wth_sts[,grep("SEED",colnames(my_cnt_annot_tab_wth_sts))]
row.names(my_annotab_for_heatmap) <- my_cnt_annot_tab_wth_sts$ID
my_annotab_for_heatmap2 <- as.data.frame(paste(my_annotab_for_heatmap$SEED1,my_annotab_for_heatmap$SEED2,my_annotab_for_heatmap$SEED3,my_annotab_for_heatmap$SEED4,sep = "  /  "))
row.names(my_annotab_for_heatmap2) <- row.names(my_annotab_for_heatmap)
heattab <- merge(mylogFC_tbl,my_stats_tbl, by = "row.names")
row.names(heattab) <- heattab$Row.names
heattab2 <- merge(heattab, my_annotab_for_heatmap2, by = "row.names")
heattab_names <- paste(heattab2[,1],heattab2[,ncol(heattab2)])
heattab3 <- heattab2[,3:(ncol(heattab2)-1)]
row.names(heattab3) <- heattab_names
heattab3_ord <- heattab3[order(heattab3$q.value),]

# select the highly DE and expressed genes with cutoffs: q.value of 0.001; abs log2fold change 2; logCPM 2
heattab3_ord_sel <- heattab3_ord[heattab3_ord$q.value < 0.001 & abs(heattab3_ord$logFC) > 1,]
heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logFC) > 2,]
heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logCPM) > 2,]
heattab3_ord_sel_fin <- heattab3_ord_sel[,1:(ncol(heattab3_ord_sel)-5)]


mycolorval <- max(abs(min(heattab3_ord_sel_fin)),abs(max(heattab3_ord_sel_fin)))

library(gplots)
cairo_pdf(file = paste("0_Heatmap_TBZ_vs_SUC_q0.001_lFC2_logCPM2b.pdf", sep = ""),height = 100, width = 15)
heatmap.2(as.matrix(heattab3_ord_sel_fin), main=""
          , col=colorpanel(n=(length(seq(-mycolorval,mycolorval,by=0.1))-1),low="green",mid="white",high="orange")
          , breaks=seq(-mycolorval,mycolorval,by=0.1)
          , trace="none", margins=c(5,70), keysize=0.4,denscol="black"
          , lmat=rbind( c(0,3), c(2,1), c(0,4) ), lhei=c(1, 20, 5), lwid=c(1,9),symkey=T,symbreaks=T, Colv = F, dendrogram = "row")
dev.off()

cairo_pdf(file = paste("0_Heatmap_TBZ_vs_SUC_q0.001_lFC2_logCPM2_bb.pdf", sep = ""),height = 100, width = 15)
heatmap.2(as.matrix(heattab3_ord_sel_fin), main=""
          , col=colorpanel(n=(length(seq(-mycolorval,mycolorval,by=0.1))-1),low="green",mid="white",high="orange")
          , breaks=seq(-mycolorval,mycolorval,by=0.1)
          , trace="none", margins=c(5,70), keysize=0.4,denscol="black"
          , lmat=rbind( c(0,3), c(2,1), c(0,4) ), lhei=c(1, 20, 5), lwid=c(1,9),symkey=T,symbreaks=T
          #, Colv = F
          , dendrogram = "row"
          )
dev.off()


# select the highly DE and expressed genes with cutoffs: q.value of 0.001; abs log2fold change 2; logCPM 2
heattab3_ord_sel <- heattab3_ord[heattab3_ord$q.value < 0.0008 & abs(heattab3_ord$logFC) > 1,]
heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logFC) > 4,]
heattab3_ord_sel <- heattab3_ord_sel[abs(heattab3_ord_sel$logCPM) > 3,]
heattab3_ord_sel_fin <- heattab3_ord_sel[,1:(ncol(heattab3_ord_sel)-5)]


mycolorval <- max(abs(min(heattab3_ord_sel_fin)),abs(max(heattab3_ord_sel_fin)))

cairo_pdf(file = paste("0_Heatmap_TBZ_vs_SUC_q0.0008_lFC4_logCPM3b.pdf", sep = ""),height = 40, width = 15)
heatmap.2(as.matrix(heattab3_ord_sel_fin), main=""
          , col=colorpanel(n=(length(seq(-mycolorval,mycolorval,by=0.1))-1),low="green",mid="white",high="orange")
          , breaks=seq(-mycolorval,mycolorval,by=0.1)
          , trace="none", margins=c(5,70), keysize=0.4,denscol="black"
          , lmat=rbind( c(0,3), c(2,1), c(0,4) ), lhei=c(1, 20, 5), lwid=c(1,9),symkey=T,symbreaks=T, Colv = F, dendrogram = "row")
dev.off()

cairo_pdf(file = paste("0_Heatmap_TBZ_vs_SUC_q0.0008_lFC4_logCPM3_bb.pdf", sep = ""),height = 40, width = 15)
heatmap.2(as.matrix(heattab3_ord_sel_fin), main=""
          , col=colorpanel(n=(length(seq(-mycolorval,mycolorval,by=0.1))-1),low="green",mid="white",high="orange")
          , breaks=seq(-mycolorval,mycolorval,by=0.1)
          , trace="none", margins=c(5,70), keysize=0.4,denscol="black"
          , lmat=rbind( c(0,3), c(2,1), c(0,4) ), lhei=c(1, 20, 5), lwid=c(1,9),symkey=T,symbreaks=T
          , dendrogram = "row"
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







# prepare the ME plots for the DE genes
cairo_pdf(paste("0_TBZvsSUCsamples_DE_ME_q.0.001_lfc1.pdf", sep = ""), height = 5, width = 6)
is.de <- decideTestsDGE(res, p.value = 0.001 , lfc = 1)

plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend=F, bty = "n", hl.cex = 0.3, pt.cex = 1, main = "TBZ vs SUC", ylab = "log2FC", xlab = "mean log CPM", ylim = c(-15,15))
abline(h = 0, lty = 2, col = "purple")
abline(h = c(-1,1), lty = 2, col = "grey40")
text(12, -10, labels = paste("SUC", " (", summary(is.de)[1]," proteins)", sep = ""), col = "blue", cex = 1.2)
text(12, 0, labels = paste("(", summary(is.de)[2]," proteins)", sep = ""), col = "black", cex = 1.2)
text(12, 10, labels = paste("TBZ", " (", summary(is.de)[3]," proteins)", sep = ""), col = "red", cex = 1.2)

dev.off()


#### overall DE plot TBZ vs SUC ----
cairo_pdf(paste("0_TBZvsSUCsamples_DE_ME_q.0.001_lfc2.pdf", sep = ""), height = 5, width = 6)

myplottbl1 <- res$table[,c(1,4)]
myplottbl2 <- merge(rowMeans(my_prot_counts), myplottbl1, by="row.names", all = T)
myplottbl3 <- myplottbl2[,-which("Row.names"%in%colnames(myplottbl2))]
row.names(myplottbl3) <- myplottbl2[,which("Row.names"%in%colnames(myplottbl2))]
myplottbl4 <- myplottbl3[complete.cases(myplottbl3),]
myplottbl4$col <- "black"
# correct the PValue name to the adjusted pvalue name and the x axis into protein RA
colnames(myplottbl4)[3] <- "FDR"
colnames(myplottbl4)[1] <- "protein mean RA"
# set the coloring conditions
fdrcutoff <- 0.001
lfccuroff <- 2
is.de <- decideTestsDGE(res, p.value = fdrcutoff , lfc = lfccuroff)
myplottbl4$col[which(myplottbl4$FDR <= fdrcutoff & myplottbl4$logFC >= lfccuroff)] <- "red"
myplottbl4$col[which(myplottbl4$FDR <= fdrcutoff & myplottbl4$logFC <= -lfccuroff)] <- "blue"
plot(log10(myplottbl4$`protein mean RA`/10e6), myplottbl4$logFC, bty = "n", pch=16, cex = 0.5, main = "TBZ vs SUC", ylab = bquote("log"[2]*"FC"), xlab = bquote("log"[10]*"("~.(colnames(myplottbl4)[1])~")"), col = myplottbl4$col)
abline(h = c(-lfccuroff,0,lfccuroff), lty = 2, col = "grey50")
myxplotval <-mean(c(summary(log10(myplottbl4$`protein mean RA`/10e6))[5],summary(log10(myplottbl4$`protein mean RA`/10e6))[6]))
text(myxplotval, mean(summary(myplottbl4$logFC)[1:2]), labels = paste("SUC", " (", summary(is.de)[1]," proteins)", sep = ""), col = "blue", cex = 1.2)
text(myxplotval, 0, labels = paste("(", summary(is.de)[2]," proteins)", sep = ""), col = "black", cex = 1.2)
text(myxplotval, -mean(summary(myplottbl4$logFC)[1:2]), labels = paste("TBZ", " (", summary(is.de)[3]," proteins)", sep = ""), col = "red", cex = 1.2)

dev.off()







#### bubble-plots ----
## add the BinTax_id column to the large tsv file

my_cnt_annot_tab_wth_sts$binid <- gsub("___.+","",my_cnt_annot_tab_wth_sts$seqid)
mybintaxonomies <- read.table("../../1_metagenome/z_accessory_files_dbs_and_executbles/binid_taxonomy.txt", header = T, sep = "\t", stringsAsFactors = F)

my_cnt_annot_tab_wth_sts$BinTaxID <- my_cnt_annot_tab_wth_sts$binid
for(mybinidno in 1:length(my_cnt_annot_tab_wth_sts$binid)){
  my_cnt_annot_tab_wth_sts$BinTaxID[mybinidno] <- mybintaxonomies[grep(my_cnt_annot_tab_wth_sts$binid[mybinidno],mybintaxonomies$binid),"BinTaxID"]
}

my_cnt_annot_tab_wth_sts <- my_cnt_annot_tab_wth_sts[complete.cases(my_cnt_annot_tab_wth_sts$ID),]
my_cnt_annot_tab_wth_sts <- merge(my_cnt_annot_tab_wth_sts, cpm(my_dgelist_fin), by.x="ID", by.y="row.names", all = T)
## prepare the dataset with all the information including the gene incidence column
mydatabubble <- data.frame(row.names = paste(my_cnt_annot_tab_wth_sts$ID,my_cnt_annot_tab_wth_sts$SEED4, sep = " "), BinTaxID = my_cnt_annot_tab_wth_sts$BinTaxID, SEED1 = my_cnt_annot_tab_wth_sts$SEED1, SEED2 = my_cnt_annot_tab_wth_sts$SEED2, SEED3 = my_cnt_annot_tab_wth_sts$SEED3, SEED4 = my_cnt_annot_tab_wth_sts$SEED4, feature_cnts = rep(1,length(my_cnt_annot_tab_wth_sts$SEED4)))
## add the meanCPM of the various tests (take care of the NAs)
tmp_cpm_vals <- my_cnt_annot_tab_wth_sts[,colnames(cpm(my_dgelist_fin))]



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
mydatabubble <- mydatabubble_bckp[complete.cases(mydatabubble_bckp[,7:8]),]
q_threshold <- 0.001

mydatabubble$TBZperSUC.TBZ.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0
mydatabubble$TBZperSUC.SUC.mean.CPM[mydatabubble$TBZperSUC.q.value > q_threshold] <- 0

mydatabubble <- mydatabubble[,1:9]

# prepare the metagenomics associated plot for the bins
library(ggplot2)
library(ggtree)

for(i in 1:4){
  mychoice <- paste("SEED",i, sep = "")
  
  
  # prep the collective table with the aggregated features and CPMs and calculate the log2FCs and add the up/downreg colours
  mydatabubble_for_metagenome <- aggregate(mydatabubble[complete.cases(mydatabubble[,mychoice]),c(6:8)], by = list(mydatabubble[complete.cases(mydatabubble[,mychoice]),mychoice],mydatabubble[complete.cases(mydatabubble[,mychoice]),"BinTaxID"]), sum)
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

  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[,-grep("mean.CPM",colnames(mydatabubble_for_metagenome))]
  
  colnames(mydatabubble_for_metagenome)[1:3] <- c("functional group","bin", "counts")
  # order the matrix according to the binid
  mydatabubble_for_metagenome <- mydatabubble_for_metagenome[order(as.numeric(gsub("[A-Z].+|_.+","",gsub("Bin_","",gsub("unbinned","1000000",gsub("unbinned_Hydrogenophaga","999999",mydatabubble_for_metagenome$bin))))),gsub("Bin_","",mydatabubble_for_metagenome$bin)),]
  mydatabubble_for_metagenome$bin <- factor(mydatabubble_for_metagenome$bin, levels = unique(mydatabubble_for_metagenome$bin))
  
  
  
  if(mychoice == "SEED4"){
    myheight <- 65
    mywidth <- 25
  } else if(mychoice == "SEED3") {
    myheight <- 50
    mywidth <- 24
  } else if (mychoice == "SEED2") {
    myheight <- 15
    mywidth <- 20
  } else if (mychoice == "SEED1") {
    myheight <- 7
    mywidth <- 17
  }
  
  cairo_pdf(paste("bubble_plot_",mychoice,"_all_q",q_threshold,".pdf", sep = ""), width = mywidth, height = myheight)
  
  g1 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +
    geom_point(aes(size=`counts`)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("Metagenome"))
  
  g2 <- ggplot(subset(mydatabubble_for_metagenome), aes(x=bin,y=`functional group`)) +

    geom_point(aes(size=`TBZperSUC.log2FC`, colour = I(TBZperSUC.log2FCcls))) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(x = "bin", y = paste("SEED level ",mychoice,sep = ""),gsub("[A-Z,a-z]","",mychoice)," functional groups", sep ="") +
    ggtitle(paste("TBZ/SUC gene expression"))
  
  multiplot(g1, g2, ncol=2)
  
  
  
  
  dev.off()
  
}


