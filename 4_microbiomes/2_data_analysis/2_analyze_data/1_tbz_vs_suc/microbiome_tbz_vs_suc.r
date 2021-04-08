library(phyloseq)
## read the master file and select the experimental data
ps_fin_ASV_pre <- readRDS("../../1_prep_master_phyloseq_object/my_master_ps.rds")
ps_fin_ASV_pre2 <- subset_samples(ps_fin_ASV_pre, exp == "tbz_vs_suc")
ps_fin_ASV <- prune_taxa(taxa_sums(ps_fin_ASV_pre2) > 0, ps_fin_ASV_pre2)

### prepare the relative abundance table and remove low abundance ASVs
ps_fin_ASV_ra_fr_prunning  = transform_sample_counts(ps_fin_ASV, function(x) x / sum(x))

# plot mean relative abundances prior taxon prunning
cairo_pdf("Barplot_of_rel_abund_prior_pruning.pdf", height = 5, width = 18)
par(mar = c(4,5,4,4))
barplot(100*colMeans(as.data.frame(otu_table(ps_fin_ASV_ra_fr_prunning))), las =2, cex.names = 0.15, ylab = "% relative abundance")
dev.off()

cairo_pdf("Barplot_of_rel_abund_prior_pruning_ylog.pdf", height = 5, width = 18)
par(mar = c(4,5,4,4))
barplot(100*colMeans(as.data.frame(otu_table(ps_fin_ASV_ra_fr_prunning))), log = "y", las =2, cex.names = 0.15, ylab = "% relative abundance")
par(xpd = F)
abline(h = 0.1)
par(xpd = T)
text(0,0.1, labels = "0.1", pos = 2,offset = 2)
dev.off()

# calculate the average mean relative abundance of the top 30 ASVs throughout the dataset
sum(colMeans(as.data.frame(otu_table(ps_fin_ASV_ra_fr_prunning)))[1:30])
#the top 20 cover the 0.9701311 of the total reads while the top 30 cover the 0.8951424
# since these ASVs cover most of the abundance of the total number of ASVs and, in their majority they cover the identifie MAGs, you can focus on them
top30ASV_names <- colnames(as.data.frame(otu_table(ps_fin_ASV_ra_fr_prunning)))[1:30]

ps_fin_ASV_prunned <- prune_taxa(top30ASV_names,ps_fin_ASV_ra_fr_prunning)
ps_fin_ASV_prunned_notRA <- prune_taxa(top30ASV_names,ps_fin_ASV)




### perform the differential abundance test
library(metacoder)

metacoder_obj <- parse_phyloseq(ps_fin_ASV_prunned)
metacoder_obj$data$tax_abund <- calc_taxon_abund(metacoder_obj, "otu_table",
                                       cols = metacoder_obj$data$sample_data$sample_id)
metacoder_obj$data$diff_table <- compare_groups(metacoder_obj, data = "tax_abund", cols = metacoder_obj$data$sample_data$sample_id, groups = metacoder_obj$data$sample_data$treatment)


cairo_pdf("metacoder_tree_top30.pdf",height = 6, width = 7)
heat_tree(metacoder_obj, 
          node_label = taxon_names,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = log2_median_ratio, # A column from `metacoder_obj$data$diff_table`
          node_color_interval = c(-8, 8), # The range of `log2_median_ratio` to display
          node_color_range = c("steelblue3", "grey90", "red"), # The color palette used
          node_size_axis_label = "ASV count",
          node_color_axis_label = "Log 2 ratio TBZ/SUC",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

dev.off()


#### prepare an MA plot
myMAtable <- merge(metacoder_obj$data$diff_table, metacoder_obj$data$otu_table, by="taxon_id", all.y = T)
myMAtable$meanRA <- round(100*rowMeans(myMAtable[,grep("^[0-9]+",colnames(myMAtable))]),1)
myMAtable$log2_median_ratio[which(myMAtable$log2_median_ratio == Inf)] <- 6
cairo_pdf("MAplot_top30.pdf")
plot(myMAtable$meanRA,myMAtable$log2_median_ratio, type = "n", bty = "n", ylim = c(-8,8), ylab = "log2(TBZ/SUC)", xlab = "meanRA (%)")
cutoff <- 2
abline(h = 0)
abline(h = c(-cutoff,cutoff), lty = 2)
points(myMAtable$meanRA[which(myMAtable$log2_median_ratio >= cutoff)],myMAtable$log2_median_ratio[which(myMAtable$log2_median_ratio >= cutoff)], pch = 21, bg = "red")
par(xpd = T)
text(myMAtable$meanRA[which(myMAtable$log2_median_ratio >= cutoff)],myMAtable$log2_median_ratio[which(myMAtable$log2_median_ratio >= cutoff)], labels = myMAtable$otu_id[which(myMAtable$log2_median_ratio >= cutoff)], pos = 3)
par(xpd = F)
points(myMAtable$meanRA[which(myMAtable$log2_median_ratio > 0 & myMAtable$log2_median_ratio < cutoff)],myMAtable$log2_median_ratio[which(myMAtable$log2_median_ratio > 0 & myMAtable$log2_median_ratio < cutoff)], pch = 21, bg = rgb(255, 0, 0, 60, maxColorValue = 255))
points(myMAtable$meanRA[which(myMAtable$log2_median_ratio < 0 & myMAtable$log2_median_ratio > -cutoff)],myMAtable$log2_median_ratio[which(myMAtable$log2_median_ratio < 0 & myMAtable$log2_median_ratio > -cutoff)], pch = 21, bg = rgb(0, 0, 255, 60, maxColorValue = 255))
points(myMAtable$meanRA[which(myMAtable$log2_median_ratio <= -cutoff)],myMAtable$log2_median_ratio[which(myMAtable$log2_median_ratio <= -cutoff)], pch = 21, bg = "blue")
par(xpd = T)
text(myMAtable$meanRA[which(myMAtable$log2_median_ratio <= -cutoff)],myMAtable$log2_median_ratio[which(myMAtable$log2_median_ratio <= -cutoff)], labels = myMAtable$otu_id[which(myMAtable$log2_median_ratio <= -cutoff)], pos = 1)
par(xpd = F)
dev.off()







#### multivariate tests ----
# load ade4 for the s.class plotting option
library(ade4)
# ad the vegan package for several data analysis and processing functions
library(vegan)
# check if linear or unimodal responses to the environmental gradients exist concerning the ASVs and use appropriate methods according to the 1st DCA axis length as described in CANOCO manual (Leps and Smileur, 2003) 
dca1 <- decorana(data.frame(otu_table(ps_fin_ASV_prunned), check.names = F))
dca2 <- decorana(decostand(data.frame(otu_table(ps_fin_ASV_prunned), check.names = F),"hellinger"))
mylist <- list(dca1,dca2)
capture.output(mylist, file="DCA_results.txt")
# the results show less than 2 SDs for both transformed and untransformed data suggesting the use of linear methods like PCA and RDA as most appropriate n either case.  

#### PCA and fitting of environmental variables ----
cairo_pdf("PCA_envfit.pdf", height=5, width=5.5)
mymat <- data.frame(otu_table(ps_fin_ASV_prunned), check.names = F)
design_for_plot <- data.frame(sample_data(ps_fin_ASV_prunned), check.names = F)
mymat_helli <- decostand(mymat,"hellinger")
mypca <- rda(mymat_helli)
varX <- 100 * round(mypca$CA$eig[1]/sum(mypca$CA$eig),3)
varY <- 100 * round(mypca$CA$eig[2]/sum(mypca$CA$eig),3)
sites <- scores(mypca, display = "sites")
mycols <- c("red","blue")
myrgbcols <- matrix(ncol = 3, nrow = 26)
for(mycolorrow in 1:(nrow(myrgbcols)/2)){
  myrgbcols[mycolorrow,1:3] <- c(col2rgb(mycols)[,1])
}
for(mycolorrow in (nrow(myrgbcols)/2+1):(nrow(myrgbcols))){
  myrgbcols[mycolorrow,1:3] <- c(col2rgb(mycols)[,2])
}
s.class(sites
        , fac = as.factor(design_for_plot$treatment):as.factor(design_for_plot$time)
        , grid=F
        , cell=2
        , label = ""
        , cstar=1
        , xlim = c(-1.3*max(abs(sites)),1.3*max(abs(sites)))
        , ylim = c(-1.3*max(abs(sites)),1.3*max(abs(sites)))
        , sub = paste("Variance: x ",varX,"%, y ",varY,"% (colour intensity  associated with tpi - h)",sep=""), csub = 1, possub = "bottomleft"
        , col = rgb(myrgbcols, alpha = seq(from = 60,to = 255,length.out = 13), maxColorValue = 255)
)

ef <- envfit(mypca, design_for_plot[,1:5], perm=1000)
capture.output(ef$vectors, file=paste("PCA_envfit.txt",sep=""))

par(xpd=T)
graphics::legend("topleft",title = paste("treatment (R2 ",round(100*ef$factors$r,0),"%, P ",round(ef$factors$pvals, 3),")", sep = ""), legend = levels(factor(design_for_plot$treatment, levels = c("TBZ","SUC"))), pch = 20, col = c("red","blue"), cex=0.8, y.intersp = 0.6, pt.cex = 2, bty="n")


dev.off()

