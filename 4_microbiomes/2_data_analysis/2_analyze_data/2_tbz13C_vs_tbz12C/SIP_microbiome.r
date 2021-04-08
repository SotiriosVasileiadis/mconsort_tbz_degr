library(phyloseq)
## read the master file and select the experimental data
ps_fin_ASV_pre <- readRDS("../../1_prep_master_phyloseq_object/my_master_ps.rds")
ps_fin_ASV_pre2 <- subset_samples(ps_fin_ASV_pre, exp == "tbz13C_vs_tbz12C")
ps_fin_ASV <- prune_taxa(taxa_sums(ps_fin_ASV_pre2) > 0, ps_fin_ASV_pre2)

### prepare the relative abundance table
ps_fin_ASV_ra  = transform_sample_counts(ps_fin_ASV, function(x) x / sum(x))

# plot mean relative abundances prior taxon prunning
cairo_pdf("Barplot_of_rel_abund.pdf", height = 5, width = 25)
par(mar = c(14,5,4,4))
barplot(100*colMeans(as.data.frame(otu_table(ps_fin_ASV_ra))), las =2, cex.names = 0.5, ylab = "% relative abundance")
dev.off()

cairo_pdf("Barplot_of_rel_abund_ylog.pdf", height = 5, width = 25)
par(mar = c(14,5,4,4))
barplot(100*colMeans(as.data.frame(otu_table(ps_fin_ASV_ra))), log = "y", las =2, cex.names = 0.5, ylab = "% relative abundance")
par(xpd = F)
abline(h = 0.1)
par(xpd = T)
text(0,0.1, labels = "0.1", pos = 2,offset = 2)
dev.off()



#### prepare the plot ----
## prepare the tables to be used in plotting

counts_bac16S <- data.frame(otu_table(ps_fin_ASV))


library(vegan)
mymatnames_forplot <- c("bacterial 16S rRNA gene")


#reduce the number of the plotable OTUs by setting a relative participation cutoff of 0.001, 0.01 and 0.05
cutoff <- 0.05
  

### tables
ps_fin_ASV_ra_glom <- tax_glom(ps_fin_ASV_ra, taxrank="Genus")
final_tab<-data.frame(otu_table(ps_fin_ASV_ra_glom), check.names = F)
colnames(final_tab) <- gsub("ASV[0-9]+ ","",colnames(final_tab))

#create a rare ASVs collection for ASVs with abundance lower than the cutoff in at least a number of samples defined here:
mynumsamples <- 3

library(vegan)
final_tab_ra<-decostand(final_tab,"total")
final_tab_abund<-final_tab_ra[,apply(final_tab_ra,2,function(x) any(length(which(x>=cutoff))>=mynumsamples))]
final_tab_rare<-rowSums(final_tab_ra[,apply(final_tab_ra,2,function(x) any(length(which(x>=cutoff))<mynumsamples))])
final_mat<- 100*data.frame(final_tab_abund,rare=final_tab_rare, check.names = F)
colnames(final_mat)<-c(colnames(final_tab_abund),paste("≤ ",100*cutoff," % rel. ab.",sep=""))


######### run the plotting command
library(PerformanceAnalytics) # for the rainbowXequal
library(RColorBrewer) # for the brewer.pal

# prepare the design file
mydesign<-data.frame(isotope = gsub("[0-9]h.+","",gsub("[0-9]+_","",sample_names(ps_fin_ASV_ra_glom))), time =  gsub("f.+","",gsub("^.+C[0-9]","",sample_names(ps_fin_ASV_ra_glom))), fraction =  gsub("^.+f","f",sample_names(ps_fin_ASV_ra_glom)), row.names = sample_names(ps_fin_ASV_ra_glom))


#prepare the matrices
mymat1 <- final_mat[row.names(mydesign[which(mydesign$isotope == "12C"),]),]
mymat2 <- final_mat[row.names(mydesign[which(mydesign$isotope == "13C" & mydesign$fraction == "f9"),]),]
mymat3 <- final_mat[row.names(mydesign[which(mydesign$isotope == "13C" & mydesign$fraction == "f5"),]),]
# impute values by duplicating the sample 21_13C1h117f5 for replacing the missing 22 sample
mymat3 <- rbind(mymat3[1:6,],mymat3[6:nrow(mymat3),])
row.names(mymat3)[6:7] <- c("21_13C1h117f5","21_13C2h117f5") 
# prep the design matrix for table 3 according to the previous 
mydesign3 <- mydesign[which(mydesign$isotope == "13C" & mydesign$fraction == "f5"),]
mydesign3 <- rbind(mydesign3[1:6,],mydesign3[6:nrow(mydesign3),])
row.names(mydesign3)[6:7] <- c("21_13C1h117f5","21_13C2h117f5") 

# prepare also the mean and the sd value mats for the stacked barplot with error bars
mymat1means_a <- aggregate(mymat1,list(as.factor(mydesign$time[which(mydesign$isotope == "12C")])), mean)
mymat1means_b <- mymat1means_a[,2:ncol(mymat1means_a)]
row.names(mymat1means_b) <- mymat1means_a[,1]
mymat1means <- t(mymat1means_b)
mymat1means <- mymat1means[,order(as.numeric(gsub("h","",colnames(as.matrix(mymat1means)))))]

mymat1sds_a <- aggregate(mymat1,list(as.factor(mydesign$time[which(mydesign$isotope == "12C")])), sd)
mymat1sds_b <- mymat1sds_a[,2:ncol(mymat1sds_a)]
row.names(mymat1sds_b) <- mymat1sds_a[,1]
mymat1sds <- t(mymat1sds_b)
mymat1sds <- mymat1sds[,order(as.numeric(gsub("h","",colnames(as.matrix(mymat1sds)))))]

mymat2means_a <- aggregate(mymat2,list(as.factor(mydesign$time[which(mydesign$isotope == "13C" & mydesign$fraction == "f9")])), mean)
mymat2means_b <- mymat2means_a[,2:ncol(mymat2means_a)]
row.names(mymat2means_b) <- mymat2means_a[,1]
mymat2means <- t(mymat2means_b)
mymat2means <- mymat2means[,order(as.numeric(gsub("h","",colnames(as.matrix(mymat2means)))))]

mymat2sds_a <- aggregate(mymat2,list(as.factor(mydesign$time[which(mydesign$isotope == "13C" & mydesign$fraction == "f9")])), sd)
mymat2sds_b <- mymat2sds_a[,2:ncol(mymat2sds_a)]
row.names(mymat2sds_b) <- mymat2sds_a[,1]
mymat2sds <- t(mymat2sds_b)
mymat2sds <- mymat2sds[,order(as.numeric(gsub("h","",colnames(as.matrix(mymat2sds)))))]

mymat3means_a <- aggregate(mymat3,list(as.factor(mydesign3$time[which(mydesign3$isotope == "13C" & mydesign3$fraction == "f5")])), mean)
mymat3means_b <- mymat3means_a[,2:ncol(mymat3means_a)]
row.names(mymat3means_b) <- mymat3means_a[,1]
mymat3means <- t(mymat3means_b)
mymat3means <- mymat3means[,order(as.numeric(gsub("h","",colnames(as.matrix(mymat3means)))))]

mymat3sds_a <- aggregate(mymat3,list(as.factor(mydesign3$time[which(mydesign3$isotope == "13C" & mydesign3$fraction == "f5")])), sd)
mymat3sds_b <- mymat3sds_a[,2:ncol(mymat3sds_a)]
row.names(mymat3sds_b) <- mymat3sds_a[,1]
mymat3sds <- t(mymat3sds_b)
mymat3sds <- mymat3sds[,order(as.numeric(gsub("h","",colnames(as.matrix(mymat3sds)))))]

#### prepare the stacked barpots ----

cairo_pdf(file=paste("stacked_barplot_with_error_cutoff_",cutoff,".pdf",sep=""),width=8,height=3.5)
layout(rbind(c(1,2,3),c(4,4,4)),heights = c(3,2))
par(mar = c(0,3,4,0))


my_brplot <- barplot(as.matrix(mymat1means), horiz=F, cex.axis=1, col=c(rainbow10equal,brewer.pal(8, name="Dark2"))
                     #, xaxt="n"
                     #, xlab = "tpi (h)"
                     , main = "12C-TBZ", cex.main = 1.5
                     , axes=TRUE, axisnames=TRUE, las=1, font=4, axis.lty=4, cex.lab = 1.5)



for(coln in 1: ncol(mymat1means)){
  for(rown in 1:nrow(mymat1means)){
    myxval <- my_brplot[coln]
    arrows(myxval, sum(mymat1means[1:rown,coln]),myxval,sum(mymat1means[1:rown,coln]) + mymat1sds[rown,coln], length = 0.1, angle = 90, col = "grey70", lwd=1)
    arrows(myxval, sum(mymat1means[1:rown,coln]),myxval,sum(mymat1means[1:rown,coln]) - mymat1sds[rown,coln], length = 0.1, angle = 90, col = "grey70", lwd=1)
    
  }
  
}


my_brplot <- barplot(as.matrix(mymat2means), horiz=F, cex.axis=1, col=c(rainbow10equal,brewer.pal(8, name="Dark2"))
                     #, xaxt="n"
                     , xlab = "tpi (h)"
                     , main = "13C-TBZ light fraction", cex.main = 1.5
                     , axes=TRUE, axisnames=TRUE, las=1, font=4, axis.lty=4, cex.lab = 1.5)



for(coln in 1: ncol(mymat2means)){
  for(rown in 1:nrow(mymat2means)){
    myxval <- my_brplot[coln]
    arrows(myxval, sum(mymat2means[1:rown,coln]),myxval,sum(mymat2means[1:rown,coln]) + mymat2sds[rown,coln], length = 0.1, angle = 90, col = "grey70", lwd=1)
    arrows(myxval, sum(mymat2means[1:rown,coln]),myxval,sum(mymat2means[1:rown,coln]) - mymat2sds[rown,coln], length = 0.1, angle = 90, col = "grey70", lwd=1)
    
  }
  
}

my_brplot <- barplot(as.matrix(mymat3means), horiz=F, cex.axis=1, col=c(rainbow10equal,brewer.pal(8, name="Dark2"))
                     #, xaxt="n"
                     , xlab = "tpi (h)"
                     , main = "13C-TBZ heavy fraction", cex.main = 1.5
                     , axes=TRUE, axisnames=TRUE, las=1, font=4, axis.lty=4, cex.lab = 1.5)



for(coln in 1: ncol(mymat3means)){
  for(rown in 1:nrow(mymat3means)){
    myxval <- my_brplot[coln]
    arrows(myxval, sum(mymat3means[1:rown,coln]),myxval,sum(mymat3means[1:rown,coln]) + mymat3sds[rown,coln], length = 0.1, angle = 90, col = "grey70", lwd=1)
    arrows(myxval, sum(mymat3means[1:rown,coln]),myxval,sum(mymat3means[1:rown,coln]) - mymat3sds[rown,coln], length = 0.1, angle = 90, col = "grey70", lwd=1)
    
  }
  
}


plot.new()
par(mar = c(0,0,0,0), xpd = T)
graphics::legend(-0.08,1, legend=row.names(mymat1means), col=c(rainbow10equal,brewer.pal(8, name="Dark2")), ncol=4, text.width=0.28, bty="n", pch=15, pt.cex=2.2, cex=1.3, y.intersp = 1.2, title = "tpi (h)\n"
                 #, title.cex = 1.3
)


dev.off()



