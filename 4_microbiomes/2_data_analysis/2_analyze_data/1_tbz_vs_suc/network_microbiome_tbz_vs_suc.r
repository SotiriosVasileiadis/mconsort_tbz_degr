library(phyloseq)
## read the master file and select the experimental data
ps_fin_ASV_pre <- readRDS("../../1_prep_master_phyloseq_object/my_master_ps.rds")
ps_fin_ASV_pre2 <- subset_samples(ps_fin_ASV_pre, exp == "tbz_vs_suc")
ps_fin_ASV <- prune_taxa(taxa_sums(ps_fin_ASV_pre2) > 0, ps_fin_ASV_pre2)

# merge the Sphingomonads that seem to belong to the same strain
ps_fin_ASV <- merge_taxa(ps_fin_ASV, eqtaxa = c("ASV004 Sphingomonas","ASV007 Sphingomonas"), archetype=1)
ps_fin_ASV <- merge_taxa(ps_fin_ASV, eqtaxa = c("ASV005 Bradyrhizobium","ASV008 Bradyrhizobium"), archetype=1)
ps_fin_ASV <- merge_taxa(ps_fin_ASV, eqtaxa = c("ASV017 Sphingopyxis","ASV020 Sphingopyxis"), archetype=1)

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
sum(colMeans(as.data.frame(otu_table(ps_fin_ASV_ra_fr_prunning)))[1:20])
#the top 20 cover the 0.9629547 of the total reads while the top 30 cover the 0.9799492
# since these ASVs cover most of the abundance of the total number of ASVs and, in their majority they cover the identifie MAGs, you can focus on them
top20ASV_names <- colnames(as.data.frame(otu_table(ps_fin_ASV_ra_fr_prunning)))[1:20]

ps_fin_ASV_prunned <- prune_taxa(top20ASV_names,ps_fin_ASV_ra_fr_prunning)
ps_fin_ASV_prunned_notRA <- prune_taxa(top20ASV_names,ps_fin_ASV)




#### Analyze and plot the data ----
### Several nice tips for preparing the network plot at: https://www.biostars.org/p/285296/
library(phyloseq) # for handling phyloseq objects
library(Hmisc) # library for the rcorr function
library(igraph) # for preparing and plotting networks

mytreatments <- c("all","TBZ","SUC")
mypvaladjs <- c("none","BH")
myrcutoffs <- c(0,0.5)
mypvalcutoffs <- c(0.05, 0.01)


for(mytreatment in mytreatments){
  # create the directory to save the files based on the treatment
  system(paste("mkdir ",mytreatment, sep = ""))
  
  for(mypvaladj in mypvaladjs){
    for(myrcutoff in myrcutoffs){
      for(mypvalcutoff in mypvalcutoffs){
        
        
        
        if(mytreatment == "all"){
          mypsobj <- ps_fin_ASV_prunned
        } else {
          mypsobj <- ps_fin_ASV_prunned
          sample_data(mypsobj)$treatment <- factor(sample_data(mypsobj)$treatment) # set to factor
          mypsobj <- subset_samples(mypsobj,treatment == mytreatment) # subset according to the treatment
          mypsobj <- prune_taxa(taxa_sums(mypsobj) > 0, mypsobj) # remove orphan taxa after the sample removal
        }
        
        # make a final check to test if the number of samples left is above the minimum required for the correlations
        myNsamp <- length(sample_sums(mypsobj))
        if(myNsamp < 5){
          next
        }
        
        # obtain the Spearman correlation values
        correl <- rcorr(as.matrix(otu_table(mypsobj)), type = "spearman")
        
        # save the r and p values in objects for downstream manipulation according to the set cutoffs
        r_vals<-correl$r
        p_vals<-correl$P
        # set the diagonal of the p_values object to one (to avoid self loops.... although this I remove also with the simplify command)
        diag(p_vals)<-1
        
        # adjust the p-values for multiple hypothesis testing
        p_vals_arr <- array(p_vals) # convert to an array
        p_vals_arr_adj <- p.adjust(p_vals_arr, method = mypvaladj) # adjust the p-values according to the selected method
        p_vals_arr_adj_mat <- matrix(p_vals_arr_adj, nrow = nrow(p_vals), ncol = ncol(p_vals)) # convert back to a matrix
        colnames(p_vals_arr_adj_mat) <- row.names(p_vals_arr_adj_mat) <- row.names(p_vals) # set the row and column names
        
        # set the r_vals according to decided adjusted p-value and r/rho cutoffs
        r_vals_padj <- r_vals # assign the new r_values
        r_vals_padj_mat <- ifelse(p_vals_arr_adj_mat > mypvalcutoff,0,r_vals_padj) # adjust by setting to zero non-significant values
        
        # in case where all r values are set to zero you need to abort the current network testing/plotting and proceed with the next test
        if(all(r_vals_padj_mat == 0)){
          next
        }
        
        # if all went well run the graph algorithm
        ## Several nice tips for preparing the network plot at: https://www.biostars.org/p/285296/
        g<-graph.adjacency(r_vals_padj_mat
                           , weighted=TRUE # for the edge calculation using continuous values
                           , diag=FALSE # we have also set the r_values to zero further up... so not necessary, but better safe than sorry
                           , mode="undirected" # undirected is the mode of choice
        )
        
        # modify plotting parameters
        E(g)[which(E(g)$weight<0)]$color <- "darkred" # Colour negative correlation edges as red
        E(g)[which(E(g)$weight>0)]$color <- "darkgreen" # Colour positive correlation edges as green
        E(g)$weight <- abs(E(g)$weight) # Convert edge weights to absolute values
        V(g)$shape <- "sphere" # Change shape of graph vertices
        V(g)$color <- "skyblue" # Change colour of graph vertices
        V(g)$vertex.frame.color <- "white" # Change colour of vertex frames
        vSizes <- 100 * colMeans(data.frame(otu_table(mypsobj)), na.rm = T) # Calculate the size of the vertices to be proportional to the relative abundance of each taxon represented by each vertex
        edgeweights <- E(g)$weight * 16 # Amplify or decrease the width of the edges
        
        # Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
        mst <- mst(g, algorithm="prim")
        
        # Identify sub-communities
        mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE) # calculate the betweenness centrality for edges in order to use for the cluster calculations 
        mst.clustering <- make_clusters(mst, membership=mst.communities$membership) # calculate the memberships
        myvertcols <- tax_table(mypsobj)[,6]
        myvertcols[which(tax_table(mypsobj)[,6] == "Sphingomonas")] <- "skyblue"
        myvertcols[which(tax_table(mypsobj)[,6] != "Sphingomonas")] <- "skyblue"
        V(mst)$color <- myvertcols
        
        # do a last sanity check
        if(length(V(mst)$color) == 0){
          next
        }
        
        
        # prepare the plot
        # start a graphics device there
        pdf(paste(mytreatment,"/corPadj_",mypvaladj,"_rCutoff_",myrcutoff,"_PValCutoff_",mypvalcutoff,".pdf", sep = ""), height = 8, width = 8)
        plot(
          mst
          , mark.groups=communities(mst.clustering)
          , layout=layout.fruchterman.reingold
          , edge.curved=TRUE
          , vertex.size=vSizes,
          , vertex.label.dist=ifelse(vSizes/5 > 2.7, 2.7, vSizes/5)
          , vertex.label.color="black"
          , asp=FALSE
          , vertex.label.cex=1.5
          , edge.width=edgeweights
          , edge.arrow.mode=0
          , main=paste("corPadj ",mypvaladj,", rCutoff ",myrcutoff,", PValCutoff ",mypvalcutoff,", n samp. ",myNsamp, sep = "")
        )
        

        # close the graphics device
        dev.off()
        
      }
    }
  }
}
