library(phyloseq)
## read the master file and select the experimental data
ps_fin_ASV_pre <- readRDS("../../1_prep_master_phyloseq_object/my_master_ps.rds")
ps_fin_ASV_pre2 <- subset_samples(ps_fin_ASV_pre, exp == "cob_test")
ps_fin_ASV <- prune_taxa(taxa_sums(ps_fin_ASV_pre2) > 0, ps_fin_ASV_pre2)

#### run ahclust to review the results
# check on the sample abundances
par(mar = c(8,4,4,4))
mybplot <- barplot(sample_sums(ps_fin_ASV), las = 2)
par(xpd = T)
text(mybplot+0.6, y = sample_sums(ps_fin_ASV)+3000, labels = (sample_sums(ps_fin_ASV)), las = 2, srt = 90, pos = 3)
# remove the sample with the very low abundance and also remove C1_47_NB12_3 that seems to be an outlier in nMDS
ps_fin_pre2 <- prune_samples(sample_sums(ps_fin_ASV) > 500,ps_fin_ASV)
ps_fin_pre2 <- prune_samples(sample_names(ps_fin_pre2)!="C1_47_NB12_3",ps_fin_pre2)

# prepare the phyloseq object with all remaining taxa with above zero counts.
ps_fin_pre3 <- prune_taxa(taxa_sums(ps_fin_pre2) > 0, ps_fin_pre2)

## identify the low abundance taxa and prune them
ps_finra_fr_prunning  = transform_sample_counts(ps_fin_pre3, function(x) x / sum(x))
# plot mean relative abundances prior taxon prunning
cairo_pdf("Barplot_of_rel_abund_prior_pruning.pdf", height = 5, width = 8)
par(mar = c(4,5,4,4))
barplot(100*colMeans(as.data.frame(otu_table(ps_finra_fr_prunning))), las =2, cex.names = 0.15, ylab = "% relative abundance")
dev.off()

cairo_pdf("Barplot_of_rel_abund_prior_pruning_ylog.pdf", height = 5, width = 8)
par(mar = c(4,5,4,4))
barplot(100*colMeans(as.data.frame(otu_table(ps_finra_fr_prunning))), log = "y", las =2, cex.names = 0.15, ylab = "% relative abundance")
par(xpd = F)
abline(h = 0.1)
par(xpd = T)
text(0,0.1, labels = "0.1", pos = 2,offset = 2)
dev.off()

### select the top 500 ASVs for analysis 
mynames_sel_for_prun <- names((as.data.frame(otu_table(ps_finra_fr_prunning)))[which(colMeans(as.data.frame(otu_table(ps_finra_fr_prunning))) >= 0.001)])

ps_finra_pruned_0.001 <- prune_taxa(mynames_sel_for_prun, ps_finra_fr_prunning) 
# hierarhical clustering plot of the complete dataset
library(vegan)
cairo_pdf("Hierarchical_clustering_UPGMA_bray.pdf", height = 5, width = 25)
plot(hclust(vegdist(data.frame(ps_finra_pruned_0.001@otu_table, check.names = F)), method = "average"), xlab = "UPGMA / Bray-Curtis", sub = "")
dev.off()

# prepare the between cycle and treatment barplot of the in use microorganisms
cairo_pdf("Barplots_per_cobalamin_content.pdf", height = 16, width = 16, onefile = T)
par(mfrow = c(4,4))
mycobalamins <- c(0,0.01,1,100)
myCcycles <- c("C1","C2","C3","C4") 
for(myCycle in myCcycles){
  for(mycobalamin in mycobalamins){
    matfrbp <- colMeans(as.data.frame(otu_table(subset_samples(ps_finra_pruned_0.001, cobalamin == mycobalamin & culture_cycle == myCycle))))
    names(matfrbp) <- data.frame(tax_table(ps_finra_pruned_0.001))$ASV_genus
    par(mar = c(15,4,4,4))
    barplot(100*matfrbp, las = 2, ylab = "% relative abundance", ylim = c(0,40), main = bquote(.(myCycle)*","~.(mycobalamin)~"ng µl"^-1~"B12"))
  }
}
dev.off()

cairo_pdf("Barplot_0.001_cutoff.pdf", height = 6, width = 8)
matfrbp <- colMeans(as.data.frame(otu_table(ps_finra_pruned_0.001)))
names(matfrbp) <- data.frame(tax_table(ps_finra_pruned_0.001))$ASV_genus
par(mar = c(15,8,2,2))
mybp <- barplot(100*matfrbp, las = 2, ylab = "% relative abundance", ylim = c(0,40), names.arg="")
text(cex=1, x=mybp+0.5, y=-1.25, names(matfrbp), xpd=TRUE, srt=45, pos = 2)
dev.off()

# prepare the between cycle and treatment barplot of the in use microorganisms for selected time-points
mytimes <- c(47,120,178,241)
ps_finra_pruned_0.001_60to100perc_degr <- subset_samples(ps_finra_pruned_0.001, hp_1st_i %in% mytimes)
cairo_pdf("Barplots_per_cobalamin_content_60_100_TBZ_degradation.pdf", height = 16, width = 16)
par(mfrow = c(4,4))
mycobalamins <- c(0,0.01,1,100)
myCcycles <- c("C1","C2","C3","C4") 
for(myCycle in myCcycles){
  for(mycobalamin in mycobalamins){
    matfrbp <- colMeans(as.data.frame(otu_table(subset_samples(ps_finra_pruned_0.001_60to100perc_degr, cobalamin == mycobalamin & culture_cycle == myCycle))))
    names(matfrbp) <- data.frame(tax_table(ps_finra_pruned_0.001_60to100perc_degr))$ASV_genus
    par(mar = c(15,4,4,4))
    barplot(100*matfrbp, las = 2, ylab = "% relative abundance", ylim = c(0,40), main = bquote(.(myCycle)*","~.(mycobalamin)~"ng µl"^-1~"B12"))
  }
}
dev.off()


cairo_pdf("Phylogenetic_tree.pdf", height = 10, width = 20)
plot_tree(ps_finra_pruned_0.001, color = "culture_cycle", label.tips = "ASV_genus", ladderize = "left", justify = "left" , size = "Abundance")
dev.off()

    
## check for positive correlations between the dominant Sphigomonad and the rest ASVs
library(Hmisc)
rcorr(as.matrix(data.frame(otu_table(subset_samples(ps_finra_pruned_0.001_60to100perc_degr,cobalamin == 0)))),type = "spearman")[[1]][,1]



### NMDS ----

library(plyr)
ps_finra_pruned_0.001_60to100perc_degr_fr_multivar <- ps_finra_pruned_0.001_60to100perc_degr
taxa_names(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar) <- data.frame(tax_table(ps_finra_pruned_0.001_60to100perc_degr))$ASV_genus
mynmds <- ordinate(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar, method="NMDS", distance="bray")

mynmdssit <- data.frame(mynmds$points)
mynmdsspe <- data.frame(mynmds$species)

mymicrmatre <- ps_finra_pruned_0.001_60to100perc_degr_fr_multivar@otu_table
mymicrmatredomOTUs <- tail(names(sort(colMeans(mymicrmatre), decreasing = T)),n=1)
minradomOTUs <- round(100*colMeans(mymicrmatre)[mymicrmatredomOTUs],1)
maxradomOTUs <- round(100*colMeans(mymicrmatre)[1],1)

mynmdssit$cycl <- factor(sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)$culture_cycle)

mynmdssit$dose <- sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)$cobalamin

meannms1 <- plyr::mapvalues(mynmdssit$cycl, from=levels(mynmdssit$cycl), to=1:length(levels(mynmdssit$cycl)))
mynmdssitmeanssdpre <- aggregate(. ~ cycl, mynmdssit[,1:3], function(x) c(mean = mean(x), sd = sd(x)))
mynmdssitmeanssd <- data.frame(mynmdssitmeanssdpre$cycl, mynmdssitmeanssdpre$MDS1, mynmdssitmeanssdpre$MDS2)

# prepare the OTU names

mytaxnms_tbl <- data.frame(tax_table(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar), stringsAsFactors = F)
for(i in 2:ncol(mytaxnms_tbl)){
  for(j in 1:nrow(mytaxnms_tbl)){
    if(is.na(mytaxnms_tbl[j,i])){
      mytaxnms_tbl[j,i] <- mytaxnms_tbl[j,i-1]
    } 
  }
}

mynmdslabs <- mytaxnms_tbl$ASV_genus
names(mynmdslabs) <- row.names(mytaxnms_tbl)

### add also the permanova test result
sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)$cobalamin <- factor(sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)$cobalamin)

library(vegan)
# PERMNOVA
mypermanova <- adonis2(data.frame(otu_table(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)) ~ culture_cycle+cobalamin, data.frame(sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)), permutations = 999, method = "bray",
                       sqrt.dist = FALSE, add = FALSE, by = "terms",
                       parallel = 4)


cairo_pdf(paste("NMDS.pdf",sep=""), height = 6, width = 7)
myplotcols <- RColorBrewer::brewer.pal(n = 4, name = 'RdBu')
# prep the plot
plot(mynmdssit[,1:2], bg = myplotcols[meannms1], frame = F, cex = 0, pch = 21, xlim = c(min(mynmdssit[,1]),max(1.3*mynmdssit[,1])), main = bquote(atop("comm ~ culture cycle + B12 treatment: R"^2~.(round(100*sum(mypermanova$R2[1:2]),1))*"%, P "*.(mypermanova$`Pr(>F)`[1])*phantom(),"cycle: R"^2~.(round(100*mypermanova$R2[1],1))*"%, P "*.(mypermanova$`Pr(>F)`[1])*" / treatment: R"^2~.(round(100*mypermanova$R2[2],1))*"%, P "*.(mypermanova$`Pr(>F)`[2]))), xlab = "nMDS1", ylab = "nMDS2")
par(xpd = T)
vegan::ordiellipse(as.matrix(mynmdssit[,1:2]), factor(mynmdssit$cycl), kind = "ehull", lty = 2, lwd=1)
points(mynmdssit[,1:2], bg = myplotcols[meannms1], cex = 1.5, pch = 20+as.numeric(factor(mynmdssit$dose)))
par(adj = 0)
title(sub = paste("stress ", round(mynmds$stress,2), sep = ""), cex.sub = 1.2)
par(adj = 1)
# title(sub = paste(" (mean RA of presented ASVs: ",minradomOTUs,"-",maxradomOTUs,"%)", sep = ""))
par(adj = .5)
library(TeachingDemos)
par(xpd = T)
# arrows(0,0,0.5*mynmdsspe[,1] , 0.5*mynmdsspe[,2], angle = 25, length = 0.15, col = rgb(153,153,153, max = 255, alpha = 100))
# # #plotrix::thigmophobe.labels(0.7*mynmdsspe[mymicrmatredomOTUs,1], 0.7*mynmdsspe[mymicrmatredomOTUs,2], labels = mynmdslabs[mymicrmatredomOTUs], cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
# text(0.5*mynmdsspe[,1], 0.5*mynmdsspe[,2], labels = mynmdslabs, cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # the color is equivalent to "grey60", but transparent
graphics::legend("topright",bty = "n", legend = c(levels(mynmdssit$cycl),levels(as.factor(as.character(mynmdssit$dose)))), pch = c(rep(22,length(levels(mynmdssit$cycl))),21,22,23,24), pt.bg = c(myplotcols[1:length(levels(mynmdssit$cycl))],rep(rgb(0, 0, 0, max = 255, alpha = 0),4)), pt.cex = 1.5)

dev.off()













### canonical analysis ----

library(plyr)
ps_finra_pruned_0.001_60to100perc_degr_fr_multivar <- ps_finra_pruned_0.001_60to100perc_degr
taxa_names(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar) <- data.frame(tax_table(ps_finra_pruned_0.001_60to100perc_degr))$ASV_genus
# check which method is more appropriate among RDA (linear gradients) and CCA (unimodal gradients) usings DCA
mydca <- ordinate(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar, method="DCA")
# Call:
#   decorana(veg = veganifyOTU(physeq)) 
# 
# Detrended correspondence analysis with 26 segments.
# Rescaling of axes with 4 iterations.
# 
# DCA1     DCA2     DCA3     DCA4
# Eigenvalues     0.04545 0.013695 0.009219 0.006607
# Decorana values 0.04783 0.009517 0.003095 0.002248
# Axis lengths    0.70392 0.364570 0.445069 0.301841
sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)$cobalamin <- factor(sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)$cobalamin)
myrda <- ordinate(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar, method="RDA", formula~culture_cycle+cobalamin, distance="bray")

library(vegan)
# PERMNOVA
mypermanova <- adonis2(data.frame(otu_table(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)) ~ culture_cycle+cobalamin, data.frame(sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)), permutations = 999, method = "bray",
                       sqrt.dist = FALSE, add = FALSE, by = "terms",
                       parallel = 4)

myrdasit <- data.frame(myrda$CCA$wa[,1:2])
myrdaspe <- data.frame(myrda$CCA$v[,1:2])

mymicrmatre <- ps_finra_pruned_0.001_60to100perc_degr_fr_multivar@otu_table
mymicrmatredomOTUs <- tail(names(sort(colMeans(mymicrmatre), decreasing = T)),n=1)
minradomOTUs <- round(100*colMeans(mymicrmatre)[mymicrmatredomOTUs],1)
maxradomOTUs <- round(100*colMeans(mymicrmatre)[1],1)

myrdasit$cycl <- factor(sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)$culture_cycle)

myrdasit$dose <- factor(sample_data(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar)$cobalamin)


meannms1 <- plyr::mapvalues(myrdasit$cycl, from=levels(myrdasit$cycl), to=1:length(levels(myrdasit$cycl)))
myrdasitmeanssdpre <- aggregate(. ~ cycl, myrdasit[,1:3], function(x) c(mean = mean(x), sd = sd(x)))
myrdasitmeanssd <- data.frame(myrdasitmeanssdpre$cycl, myrdasitmeanssdpre$RDA1[,2], myrdasitmeanssdpre$RDA2[,2])

# prepare the OTU names

mytaxnms_tbl <- data.frame(tax_table(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar), stringsAsFactors = F)
for(i in 2:ncol(mytaxnms_tbl)){
  for(j in 1:nrow(mytaxnms_tbl)){
    if(is.na(mytaxnms_tbl[j,i])){
      mytaxnms_tbl[j,i] <- mytaxnms_tbl[j,i-1]
    } 
  }
}

myrdalabs <- mytaxnms_tbl$ASV_genus
names(myrdalabs) <- row.names(mytaxnms_tbl)

#create a colour list to plot the 16 different plant season combinations

cairo_pdf(paste("RDA_SParrows.pdf",sep=""), height = 6, width = 7)
myplotcols <- RColorBrewer::brewer.pal(n = 4, name = 'RdBu')
# prep the axis variances
RDA1 <- round(100*myrda$CCA$eig[1]/sum(myrda$CCA$eig),1)
RDA2 <- round(100*myrda$CCA$eig[2]/sum(myrda$CCA$eig),1)
# prep the plot
plot(myrdasit[,1:2], bg = myplotcols[meannms1], frame = F, cex = 0, pch = 21, xlim = c(min(myrdasit[,1]),max(1.3*myrdasit[,1])), xlab = paste("RDA1(",RDA1,"%)", sep = ""), ylab = paste("RDA2(",RDA2,"%)", sep = ""), main = bquote(atop("comm ~ culture cycle + B12 treatment: R"^2~.(round(100*sum(mypermanova$R2[1:2]),1))*"%, P "*.(mypermanova$`Pr(>F)`[1])*phantom(),"cycle: R"^2~.(round(100*mypermanova$R2[1],1))*"%, P "*.(mypermanova$`Pr(>F)`[1])*" / treatment: R"^2~.(round(100*mypermanova$R2[2],1))*"%, P "*.(mypermanova$`Pr(>F)`[2]))))
par(xpd = T)
vegan::ordiellipse(as.matrix(myrdasit[,1:2]), factor(myrdasit$cycl), kind = "ehull", lty = 2, lwd=1)
points(myrdasit[,1:2], bg = myplotcols[meannms1], cex = 1.5, pch = 20+as.numeric(factor(myrdasit$dose)))
par(adj = 0)
title(sub = paste("stress ", round(myrda$stress,2), sep = ""), cex.sub = 1.2)
par(adj = 1)
title(sub = paste(" (mean RA of presented ASVs: ",minradomOTUs,"-",maxradomOTUs,"%)", sep = ""))
par(adj = .5)

graphics::legend("topright",bty = "n", legend = c(levels(myrdasit$cycl),levels(as.factor(as.character(myrdasit$dose)))), pch = c(rep(22,length(levels(myrdasit$cycl))),21,22,23,24), pt.bg = c(myplotcols[1:length(levels(myrdasit$cycl))],rep(rgb(0, 0, 0, max = 255, alpha = 0),4)), pt.cex = 1.5)

dev.off()














#### alpha div ----

mytimes <- c(47,120,178,241)


ps_fin_pre3_4_cycles_pre <- subset_samples(ps_fin_pre3, hp_1st_i %in% mytimes)
ps_fin_pre3_4_cycles <- prune_taxa(taxa_sums(ps_fin_pre3_4_cycles_pre)>0,ps_fin_pre3_4_cycles_pre)

adiv <- plot_richness(ps_fin_pre3_4_cycles, measures=c("Shannon","InvSimpson","Fisher","Observed","ACE"))
## calculate Good's coverage estimate as well
library(entropart)
good <- MetaCommunity(t(ps_fin_pre3_4_cycles@otu_table))$SampleCoverage.communities
good_tbl <- data.frame(samples = names(good), variable = rep("coverage",length(good)), value=good)

alpha_long <- rbind(adiv$data[,c("samples", "variable", "value")], good_tbl)
# convert long to wide
library(reshape2)
alpha_wide <- dcast(alpha_long, samples ~ variable, value.var="value")
row.names(alpha_wide) <- alpha_wide$samples
alpha_wide <- alpha_wide[,c("coverage","Shannon","InvSimpson","Fisher","Observed","ACE")]
colnames(alpha_wide) <- c("coverage","Shannon","Inv. Simpson","Fisher's α", "Observed", "ACE")

alpha_wide_fact <- merge(alpha_wide, data.frame(sample_data(ps_fin_pre3_4_cycles)), by = "row.names")
row.names(alpha_wide_fact) <- alpha_wide_fact$Row.names
alpha_wide_fact_final_collective <- alpha_wide_fact[-which(colnames(alpha_wide_fact)%in%"Row.names")]
alpha_wide_fact_final_collective$cobalamin <- factor(alpha_wide_fact_final_collective$cobalamin)



cairo_pdf(paste("Alphadiv_plots.pdf",sep=""), height = 7, width = 8)
par(mfrow = c(4,5))

for(mytime in mytimes){
  ps_fin_pre3_sel_pre <- subset_samples(ps_fin_pre3, hp_1st_i == mytime)
  ps_fin_pre3_sel <- prune_taxa(taxa_sums(ps_fin_pre3_sel_pre)>0,ps_fin_pre3_sel_pre)
  
  
  ### alpha diversity analysis ----
  
  adiv <- plot_richness(ps_fin_pre3_sel, measures=c("Shannon","InvSimpson","Fisher","Observed","ACE"))
  ## calculate Good's coverage estimate as well
  library(entropart)
  good <- MetaCommunity(t(ps_fin_pre3_sel@otu_table))$SampleCoverage.communities
  good_tbl <- data.frame(samples = names(good), variable = rep("coverage",length(good)), value=good)
  
  alpha_long <- rbind(adiv$data[,c("samples", "variable", "value")], good_tbl)
  # convert long to wide
  library(reshape2)
  alpha_wide <- dcast(alpha_long, samples ~ variable, value.var="value")
  row.names(alpha_wide) <- alpha_wide$samples
  alpha_wide <- alpha_wide[,c("coverage","Shannon","InvSimpson","Fisher","Observed","ACE")]
  colnames(alpha_wide) <- c("coverage","Shannon","Inv. Simpson","Fisher's α", "Observed", "ACE")
  
  alpha_wide_fact <- merge(alpha_wide, data.frame(sample_data(ps_fin_pre3_sel)), by = "row.names")
  row.names(alpha_wide_fact) <- alpha_wide_fact$Row.names
  alpha_wide_fact <- alpha_wide_fact[-which(colnames(alpha_wide_fact)%in%"Row.names")]
  alpha_wide_fact$cobalamin <- factor(alpha_wide_fact$cobalamin)
  
  
  write.table(alpha_wide, paste("T_",mytime,"_alpha_indices.txt",sep=""), quote = F, col.names = NA)
  
  library("agricolae")
  
  
  
  
  
  #### perform anova or equivalent for the alpha diversity indices ----
  
  
  library("agricolae")
  
  mystatsout <- list()
  mytestvars <- colnames(alpha_wide)[2:6]
  
  for(mytestvar in mytestvars){
    
    # create the aov matrix
    myaovmatrix <- alpha_wide_fact
    
    
    
    # run a shapiro wilk test to select parametric or non parametric analysis
    shap_out <- shapiro.test(myaovmatrix[,mytestvar])
    mystatsout[[paste(mytestvar, sep = " // ")]][["shap"]] <- shap_out
    mytestfact <- "cobalamin"
    # run the parametric or non-parametric analysis according to the shapiro.test results
    if(shap_out$p.value < 0.05){
      # non-parametric
      mykrusk <- kruskal(myaovmatrix[,mytestvar], myaovmatrix[,mytestfact], group = T, p.adj = "BH")
      
      mystatsout[[paste(mytestvar, sep = " // ")]][["krusk"]] <- mykrusk
      # prepare the barplot
      mytestvarord <- levels(myaovmatrix[,mytestfact])
      par(mar = c(4,4,4,2))
      barerrplt <- bar.err(mykrusk$means[mytestvarord[length(mytestvarord):1],], variation="SD", xlim=c(0,max(alpha_wide_fact_final_collective[,mytestvar], na.rm = T)),horiz=TRUE, bar=T, col="grey60", space=1, main= paste(mytestvar), las=1)

      # delete the significance group letters in the case that the anova was not significant
      par(xpd = T)
      if(mykrusk$statistics$p.chisq <= 0.05){
        text(x = mykrusk$means[mytestvarord[length(mytestvarord):1],1] + mykrusk$means[mytestvarord[length(mytestvarord):1],3], barerrplt$x,labels = mykrusk$groups[mytestvarord[length(mytestvarord):1],2], pos = 4, font = 3)
        par(xpd = F)
      }
    } else{
      # perform the parametric
      
      # select the alphadiv matrix and design rows
      myform <- as.formula(paste("`",mytestvar,"` ~ ",mytestfact, sep = ""))
      mymod <- aov(myform, data = myaovmatrix)
      mysumaov <- summary(mymod)
      
      mystatsout[[paste(mytestvar, sep = " // ")]][["ANOVA"]] <- mysumaov
      
      # order the matrices etc
      mytestvarord <- levels(myaovmatrix[,mytestfact])
      
      # run the Tukey test
      myHSDtest <- HSD.test(mymod, mytestfact, group=T)
      
      mystatsout[[paste(mytestvar, sep = " // ")]][["HSD test"]] <- myHSDtest
      
      # prepare the barplot
      par(mar = c(4,4,4,2))
      barerrplt <- bar.err(myHSDtest$means[mytestvarord[length(mytestvarord):1],], variation="SD", xlim=c(0,max(alpha_wide_fact_final_collective[,mytestvar], na.rm = T)),horiz=TRUE, bar=T, col="grey60", space=1, main= paste(mytestvar), las=1)

      # delete the significance group letters in the case that the anova was not significant
      par(xpd = T)
      if(mysumaov[[1]]$`Pr(>F)`[1] <= 0.05){
        text(x = myHSDtest$means[mytestvarord[length(mytestvarord):1],1] + myHSDtest$means[mytestvarord[length(mytestvarord):1],2], barerrplt$x,labels = myHSDtest$groups[mytestvarord[length(mytestvarord):1],2], pos = 4)
        par(xpd = F)
      }
    }
    
  }
  
  capture.output(mystatsout,file = paste("Alphadiv_T_",mytime,"_stats.txt",sep=""))
  
}
dev.off()










##### DA analysis Krusk ----



library(gtools)
library(agricolae)

## choose single timepoint samples only
## define the number of OTUs to be used in the analysis (in this case the 19 most abundant being clearly far more dominant than the rest according to the "Barplot_of_rel_abund_prior_pruning_ylog.pdf" were used)

mynumOTUs <- 19
# mytimes <- levels(as.factor(mydesign$time_factor))
# 
# 
# for(mytime in mytimes){
#   mydesign_sel <- mydesign[grep(mytime,mydesign$time_factor),]
#   counts_final_sel <- counts_final[row.names(mydesign_sel),]
# }
mytimes <- c(47,120,178,241)


for(mytimefrDA in mytimes){

  ps_robust_sel_time <- subset_samples(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar, hp_1st_i %in% mytimefrDA)
  ps_robust_sel_time <- prune_taxa(taxa_sums(ps_robust_sel_time) > 0,ps_robust_sel_time)
  mydesign_sel <- data.frame(sample_data(ps_robust_sel_time))
  ps_robust_sel_time <- prune_taxa(names(taxa_sums(ps_robust_sel_time)[order(taxa_sums(ps_robust_sel_time), decreasing = T)][1:mynumOTUs]),ps_robust_sel_time)
  
  counts_final_sel <- decostand(data.frame(otu_table(ps_robust_sel_time), check.names = F)[row.names(mydesign_sel),], method = "total")
  counts_final1 <- counts_final_sel
  mydesign1 <- mydesign_sel[row.names(counts_final1),]
  
  myotus <- colnames(counts_final1)
  
  mystatsout <- list()
  for(myotu in myotus){
    if(sum(counts_final1[,myotu]) == 0){
      mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
    } else {
      mykrusk <- kruskal(counts_final1[,myotu], mydesign1[row.names(counts_final1),]$cobalamin, group = T)
      mystatsout[[myotu]][["krusk"]] <- mykrusk
    }
  }
  
  mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
  names(mystatsoutkruskpvals) <- myotus
  mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "none")
  # prepare the barplot
  
  
  # prepare the ASV names
  mytax <- data.frame(tax_table(ps_robust_sel_time), stringsAsFactors = F)
  for(myrowtx in 1:nrow(mytax)){
    for(mycoltx in 2:ncol(mytax)){
      if(is.na(mytax[myrowtx,mycoltx])){
        mytax[myrowtx,mycoltx] <- mytax[myrowtx,mycoltx-1]
      }
    }
  }
  mytxplot <- data.frame(OTU = row.names(mytax), txplt = row.names(mytax))
  
  
  
  adjpvals = "pvals.adj"
  if(length(which(get(paste("mystatsoutkrusk",adjpvals, sep = "")) <= 1)) > 0){
    myotu_sigs <- names(which(get(paste("mystatsoutkrusk",adjpvals, sep = "")) <= 1))
    
    library(gtools) # for the stars.pval command
    
    # prepare the barplot with the habitat and time index
    cairo_pdf(paste("DA_kruskal_",adjpvals,"_",mytimefrDA,".pdf", sep = ""), height = 10, width = 16, onefile = T)
    par(mfrow = c(5,4))
    for(myotu_sig in myotu_sigs){
      mytestvarord <- levels(data.frame(sample_data(ps_robust_sel_time))$cobalamin)
      par(mar = c(4,10,4,8))
      barerrplt <- bar.err(100*mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],], variation="SD", xlim=c(0,110*max(mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],1] + mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],3])),horiz=TRUE, bar=T, col="grey60", space=0.5, main= paste(gsub("_"," ",mytxplot$txplt[which(mytxplot$OTU == myotu_sig)]), "\n  FDR -", round(get(paste("mystatsoutkrusk",adjpvals, sep = ""))[myotu_sig],3), stars.pval(get(paste("mystatsoutkrusk",adjpvals, sep = ""))[myotu_sig])), las=1, xlab = "%")

      # delete the significance group letters in the case that the anova was not significant
      par(xpd = T)
      if(mystatsout[[myotu_sig]]$krusk$statistics$p.chisq <= 0.05){
        text(x = 100*mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],1] + 100*mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],3], barerrplt$x,labels = mystatsout[[myotu_sig]]$krusk$groups[mytestvarord[4:1],2], pos = 4, font = 3)
        par(xpd = F)
      }
    }
    dev.off()
    
    
    
    ### prepare the heatmap for the X most abundant microorganisms (includes bacteria) ---- 
    # legend: Heatmap of the 0 to 1 transformed ranges per OTU throughout all samples of the relative abundances of the OTUs within each samples. The relative abundances of the 50 most dominant OTUs were scaled from 0 to 1 between samples.  
    library(pheatmap)
    counts_final3 <- decostand(counts_final1[,1:mynumOTUs], method = "range", MARGIN = 2)
    counts_final4 <- counts_final3
    counts_final5 <- counts_final4
    
    counts_final5t <- t(counts_final5)
    counts_final6 <- merge(mytxplot,counts_final5t, by.x = "OTU", by.y = "row.names", all = T)
    mystarspheat <- array()
    for(myotu_sig in mytxplot$OTU){
      mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkrusk",adjpvals, sep = ""))[myotu_sig])
    }
    mystarspheat <- mystarspheat[complete.cases(mystarspheat)]
    counts_final6_pre <- merge(data.frame(star=mystarspheat, row.names = names(mystarspheat)),counts_final6, by.x = "row.names",by.y = "OTU", all = T)
    counts_final6_pre$nms_star <- paste(counts_final6_pre$txplt,counts_final6_pre$star)
    
    counts_final6_fin <- counts_final6_pre[,4:(ncol(counts_final6_pre)-1)]
    row.names(counts_final6_fin) <- counts_final6_pre$nms_star
    # prepare the table with collapsed sample replicates
    counts_final6_fin_agg_reps_pre_t <- t(counts_final6_fin)
    counts_final6_fin_agg_reps_pre <- aggregate(counts_final6_fin_agg_reps_pre_t, by = list(mydesign1$cobalamin), mean)
    row.names(counts_final6_fin_agg_reps_pre) <- counts_final6_fin_agg_reps_pre$Group.1
    counts_final6_fin_agg_reps <- t(counts_final6_fin_agg_reps_pre[,-grep("Group.1",colnames(counts_final6_fin_agg_reps_pre))])
    counts_final6_fin_agg_reps <- counts_final6_fin_agg_reps[,mytestvarord]
    
    myrelabund <- data.frame("mean RA %" = rowMeans(counts_final6_fin), check.names = F)
    counts_final6_fin_plt_pre <- merge(counts_final6_fin,myrelabund, by = "row.names", all.x = T)
    counts_final6_fin_plt <- counts_final6_fin_plt_pre
    counts_final6_fin_plt <- counts_final6_fin_plt_pre[,grep("C[0-9]_",colnames(counts_final6_fin_plt_pre))]
    row.names(counts_final6_fin_plt) <- counts_final6_fin_plt_pre$Row.names
    colnames(counts_final6_fin_plt) <- paste(colnames(counts_final6_fin_plt),mydesign1$treat)

    
    # do the same for the aggregated table
    counts_final6_fin_agg_reps_plt_pre <- merge(counts_final6_fin_agg_reps,myrelabund, by = "row.names", all.x = T)
    counts_final6_fin_agg_reps_plt <- counts_final6_fin_agg_reps_plt_pre[,grep(paste(levels(mydesign1$treat),collapse = "|"),colnames(counts_final6_fin_agg_reps_plt_pre))]
    row.names(counts_final6_fin_agg_reps_plt) <- paste(counts_final6_fin_agg_reps_plt_pre$Row.names)

    #prep the plot
    mycolors <- colorRampPalette(c("grey95","lightsteelblue4","darkred","darkorange2"))(n = 119)
    
    if(mytimefrDA == "T0"){
      mycuttreecols <- 2
      mycuttreerows <- 2
    } else {
      mycuttreecols <- 2
      mycuttreerows <- 2
    }
    
    mat_fr_plot <- counts_final6_fin_plt
    myrelabund_fr_plot <- data.frame("mean RA %" = myrelabund[row.names(mat_fr_plot),], row.names = row.names(mat_fr_plot), check.names = F)
    # plot the time-specific samples
    cairo_pdf(paste("DA_pheatmap_",mytimefrDA,".pdf", sep = ""), width = 8, height = 5)
    pheatmap(mat_fr_plot, color = mycolors, cluster_rows = T, cutree_cols = mycuttreecols, cutree_rows = mycuttreerows, annotation_row = myrelabund, border_color = "grey30",cellwidth = 10, cellheight = 10)
    dev.off()
    
    
    # plot the aggregated samples for T04
    if(mytimefrDA == "T0"){
      mycuttreecols <- 2
      mycuttreerows <- 2
    } else {
      mycuttreecols <- 2
      mycuttreerows <- 2
    }
    
    mat_fr_plot <- counts_final6_fin_agg_reps_plt
    mat_fr_plot <- mat_fr_plot[,-which(colnames(mat_fr_plot)%in%c("Row.names","mean RA %"))]
    myrelabund_fr_plot <- data.frame("mean RA %" = myrelabund[row.names(mat_fr_plot),], row.names = row.names(mat_fr_plot), check.names = F)
    cairo_pdf(paste("DA_pheatmap_total_reps_agg_",mytimefrDA,".pdf", sep = ""), width = 6, height = 5)
    pheatmap(decostand(mat_fr_plot, "range", MARGIN = 1), color = mycolors, cluster_rows = T, cutree_cols = mycuttreecols, cutree_rows = mycuttreecols, annotation_row = myrelabund_fr_plot, border_color = "grey30",cellwidth = 10, cellheight = 10)
    dev.off()
    
    
    
    
  }
}























#### perform the DA without separating the samples according to culture cycles ----

ps_robust_sel_time <- subset_samples(ps_finra_pruned_0.001_60to100perc_degr_fr_multivar, hp_1st_i %in% c(40,120,241))
ps_robust_sel_time <- prune_taxa(taxa_sums(ps_robust_sel_time) > 0,ps_robust_sel_time)
mydesign_sel <- data.frame(sample_data(ps_robust_sel_time))
ps_robust_sel_time <- prune_taxa(names(taxa_sums(ps_robust_sel_time)[order(taxa_sums(ps_robust_sel_time), decreasing = T)][1:mynumOTUs]),ps_robust_sel_time)

counts_final_sel <- decostand(data.frame(otu_table(ps_robust_sel_time), check.names = F)[row.names(mydesign_sel),], method = "total")
counts_final1 <- counts_final_sel
mydesign1 <- mydesign_sel[row.names(counts_final1),]

myotus <- colnames(counts_final1)

mystatsout <- list()
for(myotu in myotus){
  if(sum(counts_final1[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(counts_final1[,myotu], mydesign1[row.names(counts_final1),]$cobalamin, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "none")
# prepare the barplot


# prepare the ASV names
mytax <- data.frame(tax_table(ps_robust_sel_time), stringsAsFactors = F)
for(myrowtx in 1:nrow(mytax)){
  for(mycoltx in 2:ncol(mytax)){
    if(is.na(mytax[myrowtx,mycoltx])){
      mytax[myrowtx,mycoltx] <- mytax[myrowtx,mycoltx-1]
    }
  }
}
mytxplot <- data.frame(OTU = row.names(mytax), txplt = row.names(mytax))



adjpvals = "pvals.adj"
if(length(which(get(paste("mystatsoutkrusk",adjpvals, sep = "")) <= 1)) > 0){
  myotu_sigs <- names(which(get(paste("mystatsoutkrusk",adjpvals, sep = "")) <= 1))
  
  library(gtools) # for the stars.pval command
  
  # prepare the barplot with the habitat and time index
  cairo_pdf(paste("DA_kruskal_",adjpvals,"_alltimes.pdf", sep = ""), height = 10, width = 16, onefile = T)
  par(mfrow = c(5,4))
  for(myotu_sig in myotu_sigs){
    mytestvarord <- levels(data.frame(sample_data(ps_robust_sel_time))$cobalamin)
    par(mar = c(4,10,4,8))
    barerrplt <- bar.err(100*mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],], variation="SD", xlim=c(0,110*max(mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],1] + mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],3])),horiz=TRUE, bar=T, col="grey60", space=0.5, main= paste(gsub("_"," ",mytxplot$txplt[which(mytxplot$OTU == myotu_sig)]), "\n  FDR -", round(get(paste("mystatsoutkrusk",adjpvals, sep = ""))[myotu_sig],3), stars.pval(get(paste("mystatsoutkrusk",adjpvals, sep = ""))[myotu_sig])), las=1, xlab = "%")

    # delete the significance group letters in the case that the anova was not significant
    par(xpd = T)
    if(mystatsout[[myotu_sig]]$krusk$statistics$p.chisq <= 0.05){
      text(x = 100*mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],1] + 100*mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],3], barerrplt$x,labels = mystatsout[[myotu_sig]]$krusk$groups[mytestvarord[4:1],2], pos = 4, font = 3)
      par(xpd = F)
    }
  }
  dev.off()
  
  
  
  ### prepare the heatmap for the X most abundant microorganisms (includes bacteria) ---- 
  # legend: Heatmap of the 0 to 1 transformed ranges per OTU throughout all samples of the relative abundances of the OTUs within each samples. The relative abundances of the 50 most dominant OTUs were scaled from 0 to 1 between samples.  
  library(pheatmap)
  counts_final3 <- decostand(counts_final1[,1:mynumOTUs], method = "range", MARGIN = 2)
  counts_final4 <- counts_final3
  counts_final5 <- counts_final4
  counts_final5t <- t(counts_final5)
  counts_final6 <- merge(mytxplot,counts_final5t, by.x = "OTU", by.y = "row.names", all = T)
  mystarspheat <- array()
  for(myotu_sig in mytxplot$OTU){
    mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkrusk",adjpvals, sep = ""))[myotu_sig])
  }
  mystarspheat <- mystarspheat[complete.cases(mystarspheat)]
  counts_final6_pre <- merge(data.frame(star=mystarspheat, row.names = names(mystarspheat)),counts_final6, by.x = "row.names",by.y = "OTU", all = T)
  counts_final6_pre$nms_star <- paste(counts_final6_pre$txplt,counts_final6_pre$star)
  
  counts_final6_fin <- counts_final6_pre[,4:(ncol(counts_final6_pre)-1)]
  row.names(counts_final6_fin) <- counts_final6_pre$nms_star
  # prepare the table with collapsed sample replicates
  counts_final6_fin_agg_reps_pre_t <- t(counts_final6_fin)
  counts_final6_fin_agg_reps_pre <- aggregate(counts_final6_fin_agg_reps_pre_t, by = list(mydesign1$cobalamin), mean)
  row.names(counts_final6_fin_agg_reps_pre) <- counts_final6_fin_agg_reps_pre$Group.1
  counts_final6_fin_agg_reps <- t(counts_final6_fin_agg_reps_pre[,-grep("Group.1",colnames(counts_final6_fin_agg_reps_pre))])
  counts_final6_fin_agg_reps <- counts_final6_fin_agg_reps[,mytestvarord]
  
  myrelabund <- data.frame("mean RA %" = rowMeans(counts_final6_fin), check.names = F)
  counts_final6_fin_plt_pre <- merge(counts_final6_fin,myrelabund, by = "row.names", all.x = T)
  counts_final6_fin_plt <- counts_final6_fin_plt_pre
  counts_final6_fin_plt <- counts_final6_fin_plt_pre[,grep("C[0-9]_",colnames(counts_final6_fin_plt_pre))]
  row.names(counts_final6_fin_plt) <- counts_final6_fin_plt_pre$Row.names
  colnames(counts_final6_fin_plt) <- paste(colnames(counts_final6_fin_plt),mydesign1$treat)

  
  # do the same for the aggregated table
  counts_final6_fin_agg_reps_plt_pre <- merge(counts_final6_fin_agg_reps,myrelabund, by = "row.names", all.x = T)
  counts_final6_fin_agg_reps_plt <- counts_final6_fin_agg_reps_plt_pre[,grep(paste(levels(mydesign1$treat),collapse = "|"),colnames(counts_final6_fin_agg_reps_plt_pre))]
  row.names(counts_final6_fin_agg_reps_plt) <- paste(counts_final6_fin_agg_reps_plt_pre$Row.names)
  mycolors <- colorRampPalette(c("grey95","lightsteelblue4","darkred","darkorange2"))(n = 119)
  
  if(mytimefrDA == "T0"){
    mycuttreecols <- 2
    mycuttreerows <- 2
  } else {
    mycuttreecols <- 2
    mycuttreerows <- 2
  }
  
  mat_fr_plot <- counts_final6_fin_plt
  myrelabund_fr_plot <- data.frame("mean RA %" = myrelabund[row.names(mat_fr_plot),], row.names = row.names(mat_fr_plot), check.names = F)
  # plot the time-specific samples
  cairo_pdf(paste("DA_pheatmap_alltimes.pdf", sep = ""), width = 8, height = 5)
  pheatmap(mat_fr_plot, color = mycolors, cluster_rows = T, cutree_cols = mycuttreecols, cutree_rows = mycuttreerows, annotation_row = myrelabund, border_color = "grey30",cellwidth = 10, cellheight = 10)
  dev.off()
  
  
  # plot the aggregated samples for T04
  if(mytimefrDA == "T0"){
    mycuttreecols <- 2
    mycuttreerows <- 2
  } else {
    mycuttreecols <- 2
    mycuttreerows <- 2
  }
  
  mat_fr_plot <- counts_final6_fin_agg_reps_plt
  mat_fr_plot <- mat_fr_plot[,-which(colnames(mat_fr_plot)%in%c("Row.names","mean RA %"))]
  myrelabund_fr_plot <- data.frame("mean RA %" = myrelabund[row.names(mat_fr_plot),], row.names = row.names(mat_fr_plot), check.names = F)
  cairo_pdf(paste("DA_pheatmap_total_reps_agg_alltimes.pdf", sep = ""), width = 6, height = 5)
  pheatmap(decostand(mat_fr_plot, "range", MARGIN = 1), color = mycolors, cluster_rows = T, cutree_cols = mycuttreecols, cutree_rows = mycuttreecols, annotation_row = myrelabund_fr_plot, border_color = "grey30",cellwidth = 10, cellheight = 10)
  dev.off()
}

























#### perform the DA without separating the samples according to culture cycles using samples of all degradation stages----




ps_robust_sel_time <- ps_finra_pruned_0.001
ps_robust_sel_time <- prune_taxa(taxa_sums(ps_robust_sel_time) > 0,ps_robust_sel_time)
sample_data(ps_robust_sel_time)$cobalamin <- factor(sample_data(ps_robust_sel_time)$cobalamin)
mydesign_sel <- data.frame(sample_data(ps_robust_sel_time))
ps_robust_sel_time <- prune_taxa(names(taxa_sums(ps_robust_sel_time)[order(taxa_sums(ps_robust_sel_time), decreasing = T)][1:mynumOTUs]),ps_robust_sel_time)

counts_final_sel <- decostand(data.frame(otu_table(ps_robust_sel_time), check.names = F)[row.names(mydesign_sel),], method = "total")
counts_final1 <- counts_final_sel
mydesign1 <- mydesign_sel[row.names(counts_final1),]

myotus <- colnames(counts_final1)

mystatsout <- list()
for(myotu in myotus){
  if(sum(counts_final1[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(counts_final1[,myotu], mydesign1[row.names(counts_final1),]$cobalamin, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "none")
# prepare the barplot


# prepare the ASV names
mytax <- data.frame(tax_table(ps_robust_sel_time), stringsAsFactors = F)
for(myrowtx in 1:nrow(mytax)){
  for(mycoltx in 2:ncol(mytax)){
    if(is.na(mytax[myrowtx,mycoltx])){
      mytax[myrowtx,mycoltx] <- mytax[myrowtx,mycoltx-1]
    }
  }
}
mytxplot <- data.frame(OTU = mytax$ASV_genus, txplt = mytax$ASV_genus)



adjpvals = "pvals.adj"
if(length(which(get(paste("mystatsoutkrusk",adjpvals, sep = "")) <= 1)) > 0){
  myotu_sigs <- names(which(get(paste("mystatsoutkrusk",adjpvals, sep = "")) <= 1))
  
  library(gtools) # for the stars.pval command
  
  # prepare the barplot with the habitat and time index
  cairo_pdf(paste("DA_kruskal_",adjpvals,"_alltimes_all_degr_stages.pdf", sep = ""), height = 10, width = 16, onefile = T)
  par(mfrow = c(5,4))
  for(myotu_sig in myotu_sigs){
    mytestvarord <- levels(data.frame(sample_data(ps_robust_sel_time))$cobalamin)
    par(mar = c(4,10,4,8))
    barerrplt <- bar.err(100*mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],], variation="SD", xlim=c(0,110*max(mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],1] + mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],3], na.rm = T)),horiz=TRUE, bar=T, col="grey60", space=0.5, main= paste(gsub("_"," ",mytxplot$txplt[which(mytxplot$OTU == myotu_sig)]), "\n  FDR -", round(get(paste("mystatsoutkrusk",adjpvals, sep = ""))[myotu_sig],3), stars.pval(get(paste("mystatsoutkrusk",adjpvals, sep = ""))[myotu_sig])), las=1, xlab = "%")

    # delete the significance group letters in the case that the anova was not significant
    par(xpd = T)
    if(mystatsout[[myotu_sig]]$krusk$statistics$p.chisq <= 0.05){
      text(x = 100*mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],1] + 100*mystatsout[[myotu_sig]]$krusk$means[mytestvarord[4:1],3], barerrplt$x,labels = mystatsout[[myotu_sig]]$krusk$groups[mytestvarord[4:1],2], pos = 4, font = 3)
      par(xpd = F)
    }
  }
  dev.off()
  
  
  
  ### prepare the heatmap for the X most abundant microorganisms (includes bacteria) ---- 
  # legend: Heatmap of the 0 to 1 transformed ranges per OTU throughout all samples of the relative abundances of the OTUs within each samples. The relative abundances of the 50 most dominant OTUs were scaled from 0 to 1 between samples.  
  library(pheatmap)
  counts_final3 <- decostand(counts_final1[,1:mynumOTUs], method = "range", MARGIN = 2)
  counts_final4 <- counts_final3
  counts_final5 <- counts_final4
  counts_final5t <- t(counts_final5)
  counts_final6 <- merge(mytxplot,counts_final5t, by.x = "OTU", by.y = "row.names", all = T)
  mystarspheat <- array()
  for(myotu_sig in mytxplot$OTU){
    mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkrusk",adjpvals, sep = ""))[myotu_sig])
  }
  mystarspheat <- mystarspheat[complete.cases(mystarspheat)]
  counts_final6_pre <- merge(data.frame(star=mystarspheat, row.names = names(mystarspheat)),counts_final6, by.x = "row.names",by.y = "OTU", all = T)
  counts_final6_pre$nms_star <- paste(counts_final6_pre$txplt,counts_final6_pre$star)
  
  counts_final6_fin <- counts_final6_pre[,4:(ncol(counts_final6_pre)-1)]
  row.names(counts_final6_fin) <- counts_final6_pre$nms_star
  # prepare the table with collapsed sample replicates
  counts_final6_fin_agg_reps_pre_t <- t(counts_final6_fin)
  counts_final6_fin_agg_reps_pre <- aggregate(counts_final6_fin_agg_reps_pre_t, by = list(mydesign1$cobalamin), mean)
  row.names(counts_final6_fin_agg_reps_pre) <- counts_final6_fin_agg_reps_pre$Group.1
  counts_final6_fin_agg_reps <- t(counts_final6_fin_agg_reps_pre[,-grep("Group.1",colnames(counts_final6_fin_agg_reps_pre))])
  counts_final6_fin_agg_reps <- counts_final6_fin_agg_reps[,mytestvarord]
  
  myrelabund <- data.frame("mean RA %" = rowMeans(counts_final6_fin), check.names = F)
  counts_final6_fin_plt_pre <- merge(counts_final6_fin,myrelabund, by = "row.names", all.x = T)
  counts_final6_fin_plt <- counts_final6_fin_plt_pre
  counts_final6_fin_plt <- counts_final6_fin_plt_pre[,grep("C[0-9]_",colnames(counts_final6_fin_plt_pre))]
  row.names(counts_final6_fin_plt) <- counts_final6_fin_plt_pre$Row.names
  colnames(counts_final6_fin_plt) <- paste(colnames(counts_final6_fin_plt),mydesign1$treat)

  # do the same for the aggregated table
  counts_final6_fin_agg_reps_plt_pre <- merge(counts_final6_fin_agg_reps,myrelabund, by = "row.names", all.x = T)
  counts_final6_fin_agg_reps_plt <- counts_final6_fin_agg_reps_plt_pre[,grep(paste(levels(mydesign1$treat),collapse = "|"),colnames(counts_final6_fin_agg_reps_plt_pre))]
  row.names(counts_final6_fin_agg_reps_plt) <- paste(counts_final6_fin_agg_reps_plt_pre$Row.names)

  
  #prep the plot
  mycolors <- colorRampPalette(c("grey95","lightsteelblue4","darkred","darkorange2"))(n = 119)
  
  if(mytimefrDA == "T0"){
    mycuttreecols <- 2
    mycuttreerows <- 2
  } else {
    mycuttreecols <- 2
    mycuttreerows <- 2
  }
  
  mat_fr_plot <- counts_final6_fin_plt
  myrelabund_fr_plot <- data.frame("mean RA %" = myrelabund[row.names(mat_fr_plot),], row.names = row.names(mat_fr_plot), check.names = F)
  # plot the time-specific samples
  cairo_pdf(paste("DA_pheatmap_alltimes_all_degr_stages.pdf", sep = ""), width = 8, height = 5)
  pheatmap(mat_fr_plot, color = mycolors, cluster_rows = T, cutree_cols = mycuttreecols, cutree_rows = mycuttreerows, annotation_row = myrelabund, border_color = "grey30",cellwidth = 10, cellheight = 10)
  dev.off()
  
  
  # plot the aggregated samples for T04
  if(mytimefrDA == "T0"){
    mycuttreecols <- 2
    mycuttreerows <- 2
  } else {
    mycuttreecols <- 2
    mycuttreerows <- 2
  }
  
  mat_fr_plot <- counts_final6_fin_agg_reps_plt
  mat_fr_plot <- mat_fr_plot[,-which(colnames(mat_fr_plot)%in%c("Row.names","mean RA %"))]
  myrelabund_fr_plot <- data.frame("mean RA %" = myrelabund[row.names(mat_fr_plot),], row.names = row.names(mat_fr_plot), check.names = F)
  cairo_pdf(paste("DA_pheatmap_total_reps_agg_alltimes_all_degr_stages.pdf", sep = ""), width = 7, height = 5)
  pheatmap(decostand(mat_fr_plot, "range", MARGIN = 1), color = mycolors, cluster_rows = T, cutree_cols = mycuttreecols, cutree_rows = mycuttreecols, annotation_row = myrelabund_fr_plot, border_color = "grey30",cellwidth = 10, cellheight = 10)
  dev.off()
  
  
  
  
}
























#### prep the correlations ----
ps_robust_sel_time <- ps_finra_pruned_0.001
ps_robust_sel_time <- prune_taxa(taxa_sums(ps_robust_sel_time) > 0,ps_robust_sel_time)
sample_data(ps_robust_sel_time)$cobalamin <- factor(sample_data(ps_robust_sel_time)$cobalamin)
mydesign_sel <- data.frame(sample_data(ps_robust_sel_time))
ps_robust_sel_time <- prune_taxa(names(taxa_sums(ps_robust_sel_time)[order(taxa_sums(ps_robust_sel_time), decreasing = T)][1:mynumOTUs]),ps_robust_sel_time)



cairo_pdf("Correlation_tests.pdf",height = 4, width = 8)
par(mfrow = c(1,2))
ps_noB12 <- subset_samples(ps_robust_sel_time, cobalamin == 0)
sphing <- 100*data.frame(otu_table(ps_noB12))$ASV001.Sphingomonas + 100*data.frame(otu_table(ps_noB12))$ASV002.Sphingomonas
hydrog <- 100*data.frame(otu_table(ps_noB12))$ASV009.Hydrogenophaga + 100*data.frame(otu_table(ps_noB12))$ASV010.Hydrogenophaga
mymod <- lm(sphing~hydrog)
myrcorr <- rcorr(sphing, hydrog)
plot(hydrog, sphing, bty="n", pch = 16, ylab = "Sphingomonas RA (%)", xlab = "Hydrogenophaga RA (%)", main = paste("LM: y = ",round(mymod$coefficients["hydrog"],2),"x + ",round(mymod$coefficients["(Intercept)"],1),", r ",round(myrcorr$r[1,2],2),", P ",round(myrcorr$P[1,2],2),"\nno B12 samples", sep = ""))
abline(mymod)


ps_rest <- subset_samples(ps_robust_sel_time, cobalamin != 0)
sphing <- 100*data.frame(otu_table(ps_rest))$ASV001.Sphingomonas + 100*data.frame(otu_table(ps_noB12))$ASV002.Sphingomonas
hydrog <- 100*data.frame(otu_table(ps_rest))$ASV009.Hydrogenophaga + 100*data.frame(otu_table(ps_noB12))$ASV010.Hydrogenophaga
mymod <- lm(sphing~hydrog)
myrcorr <- rcorr(sphing, hydrog)
plot(hydrog, sphing, bty="n", pch = 16, ylab = "Sphingomonas RA (%)", xlab = "Hydrogenophaga RA (%)", main = paste("LM: y = ",round(mymod$coefficients["hydrog"],2),"x + ",round(mymod$coefficients["(Intercept)"],1),", r ",round(myrcorr$r[1,2],2),", P ",round(myrcorr$P[1,2],2),"\nall samples", sep = ""))
abline(mymod)

dev.off()








#### kinetics ----

## prep the kinetics data table according to the required format

mywidekindta <- read.table("my_degr_data.txt", header = T, sep = "\t", check.names = F)
library(reshape2)
mylongkindta <- melt(mywidekindta, id.vars=c("post_inoc", "culture"))
mylongkindta$rep <- gsub("sample","",gsub("_.+","",mylongkindta$variable))
mylongkindta$variable <- gsub("sample[0-9]","treat",mylongkindta$variable)
mylongwidekindta <- dcast(mylongkindta, culture + post_inoc + rep ~ variable, value.var="value")



# prep the table with the DT50s
myDT50tbl <- data.frame(cycle=base::rep(c("C1","C2","C3","C4"),each =4), treat = paste("treat_",base::rep(c("nb12","0.01","1","100"), 4), sep = ""), value = rep(NA, 4, each =4))

library(mkin)

mycultures <- c("C1","C2","C3","C4")


cairo_pdf(file=paste("TBZ_kinetics.pdf",sep = ""),height=8,width=12)
par(xpd=T,bty="n",mar=c(5,5,2,2),mfcol=c(4,4))

for(myculture in mycultures){
  mylongwidekindtacult <- mylongwidekindta[which(mylongwidekindta$culture == myculture),c(2,7,4:6)]

    # use reshape to convert the tables from wide to long and import them in the following analysis
    library("reshape2")

    for(mytreat in colnames(mylongwidekindtacult)[2:5]){
      mynamenum <- which(colnames(mylongwidekindtacult)%in%mytreat)
      C <- mylongwidekindtacult[,c(1,mynamenum)]
      
      colnames(C) <- c("t","parent")
      
      C$t <- as.numeric(as.character(C$t))
      
      # Import the data
      
      # the data should follow the format provided in the example below (tab delimited text)
      # t	parent
      # 0	1.90
      # 0	2.01
      # 0	1.90
      # 4	0.25
      # 4	0.22
      # 4	0.19
      # 7	0.09
      # 7	0.11
      # 7	0.12
      # 15	0.11
      # 15	0.14
      # 15	0.17
      # 21	0.16
      # 21	0.15
      # 21	0.14
      # 28	0.16
      # 28	0.14
      # 28	0.14
      
      # the following line is not necessary for this script
      #C=read.table("data.txt", header=T,dec=",")
      
      # prepare the format for the mkin package
      
      C_mkin=data.frame(
        name = rep("parent",nrow(C)),
        time = C[,1],
        value = C[,2] + 0.01 # had to add a small amount in order for for some zeroed values to be plotted
      )
      
      
      # use mkin for obtaining the start values for the constants (steps: 1) define the models; 2) run iterations; 3) get the model summary)
      
      

      HS=mkinmod(parent=list(type="HS"))
      HS.fit=mkinfit(HS,C_mkin,plot=T)
      summary(HS.fit)
      capture.output(summary(HS.fit),file=paste("kinetics_",myculture,"_",mytreat,"_summary_HS.txt",sep = ""))
      
      
      # set parameteres (most important the pecticide name, the residual range and the legend location)
      pesticide_name=mytreat
      plot(C_mkin[,2:3], ylab=bquote(atop(.(gsub("treat_"," ",pesticide_name))~"ng ml"^-1*phantom(),"% initial conc")),xlab="time (h)", xlim = c(0,100), ylim = c(0,100), main = myculture, cex.lab = 1.5, cex.main = 2)
      
      mysumm <- summary(HS.fit)
      mylinesfrplt <- stats::aggregate(predicted ~ time, HS.fit$data[,c(1,4)], mean)
      lines(mylinesfrplt,lty=1, lwd = 2)
      arrows(mysumm$distimes$DT50,0,mysumm$distimes$DT50,50, length = 0, lty = 2)
      text(mysumm$distimes$DT50, 50, labels = bquote("DT"[50]~.(round(mysumm$distimes$DT50, 1))~h), pos = 4, cex = 1.5)
      
      # fill the table for the ANOVA test
      myDT50tbl[which(myDT50tbl$cycle == myculture & myDT50tbl$treat == mytreat),3] <- mysumm$distimes$DT50
      

      
    }
    
    

}

dev.off()






##### prep the barplot for the DT50s ----

write.table(myDT50tbl, file = "my_DT50s.txt", sep = "\t", quote =F)
myDT50tbl





myDT50tbl_noC1 <- myDT50tbl[myDT50tbl$cycle != "C1",]

library("agricolae")

mystatsout <- list()
mytestvars <- colnames(myDT50tbl_noC1)[3]
cairo_pdf(paste("my_DT50s.pdf",sep=""), height = 3, width = 3)
for(mytestvar in mytestvars){
  
  # create the aov matrix
  myaovmatrix <- myDT50tbl_noC1
  myaovmatrix$treat <- factor(myaovmatrix$treat, levels = c("treat_nb12","treat_0.01","treat_1","treat_100"))
  
  
  # run a shapiro wilk test to select parametric or non parametric analysis
  shap_out <- shapiro.test(myaovmatrix[,mytestvar])
  mystatsout[[paste(mytestvar, sep = " // ")]][["shap"]] <- shap_out
  mytestfact <- "treat"
  # run the parametric or non-parametric analysis according to the shapiro.test results
  if(shap_out$p.value < 0.05){
    # non-parametric
    mykrusk <- kruskal(myaovmatrix[,mytestvar], myaovmatrix[,mytestfact], group = T, p.adj = "BH")
    
    mystatsout[[paste(mytestvar, sep = " // ")]][["krusk"]] <- mykrusk
    # prepare the barplot
    mytestvarord <- levels(myaovmatrix[,mytestfact])
    par(mar = c(4,4,4,2))
    myerrbrpltmat <- mykrusk$means[mytestvarord[length(mytestvarord):1],]
    row.names(myerrbrpltmat) <- gsub("treat_","",row.names(myerrbrpltmat))
    barerrplt <- bar.err(myerrbrpltmat, variation="SD", xlim=c(0,1.2*max(myaovmatrix[,mytestvar], na.rm = T)),horiz=TRUE, bar=T, col="grey60", space=1, main= bquote("EC"[50]*"s"), las=1, xlab = "time (h)")

    # delete the significance group letters in the case that the anova was not significant
    if(mykrusk$statistics$p.chisq <= 0.05){
      par(xpd = T)
      text(x = mykrusk$means[mytestvarord[length(mytestvarord):1],1] + mykrusk$means[mytestvarord[length(mytestvarord):1],3], barerrplt$x,labels = mykrusk$groups[mytestvarord[length(mytestvarord):1],2], pos = 4, font = 3)
      par(xpd = F)
    }
  } else{
    # perform the parametric
    
    # select the alphadiv matrix and design rows
    myform <- as.formula(paste("`",mytestvar,"` ~ ",mytestfact, sep = ""))
    mymod <- aov(myform, data = myaovmatrix)
    mysumaov <- summary(mymod)
    
    mystatsout[[paste(mytestvar, sep = " // ")]][["ANOVA"]] <- mysumaov
    
    # order the matrices etc
    mytestvarord <- levels(myaovmatrix[,mytestfact])
    
    # run the Tukey test
    myHSDtest <- HSD.test(mymod, mytestfact, group=T)
    
    mystatsout[[paste(mytestvar, sep = " // ")]][["HSD test"]] <- myHSDtest
    
    # prepare the barplot
    par(mar = c(4,4,4,2))
    myerrbrpltmat <- myHSDtest$means[mytestvarord[length(mytestvarord):1],]
    row.names(myerrbrpltmat) <- gsub("treat_","",row.names(myerrbrpltmat))
    barerrplt <- bar.err(myerrbrpltmat, variation="SD", xlim=c(0,1.2*max(myaovmatrix[,mytestvar], na.rm = T)),horiz=TRUE, bar=T, col="grey60", space=1, main=  bquote("EC"[50]*"s"), las=1, xlab = "time (h)")

    # delete the significance group letters in the case that the anova was not significant
    if(mysumaov[[1]]$`Pr(>F)`[1] <= 0.05){
      par(xpd = T)
      text(x = myHSDtest$means[mytestvarord[length(mytestvarord):1],1] + myHSDtest$means[mytestvarord[length(mytestvarord):1],2], barerrplt$x,labels = myHSDtest$groups[mytestvarord[length(mytestvarord):1],2], pos = 4)
      par(xpd = F)
    }
  }
  
}

capture.output(mystatsout,file = paste("my_DT50s_stats.txt",sep=""))

dev.off()

