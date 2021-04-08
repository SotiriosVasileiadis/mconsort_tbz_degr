library(genoPlotR)

gpm50ctg243<-read_dna_seg_from_file("../../z_accessory_files_dbs_and_executbles/groopm_bin_50___opera_contig_243.gbk")
# curate the annotation according to database hits and product description
gpm50ctg243$gene[1] <- "andR part." # correct the annotation providing the nearly identical match according to blast against the BAF03250.1 GenPept accession
gpm50ctg243$product[1] <- "regulatory"
gpm50ctg243[c(11,12,13,15,22,23,24,28,32,33,34,35,38,40,42),"gene"] <- c("hypoth.","IS element","IS element","TonB-dep. rec.","hypoth.","hypoth.","hypoth.","hypoth.","hypoth.","hypoth.","hypoth.","membrane p.","haloacid dehalog.","hypoth.","TnpA repr.")
KA1pCAR3 <- read_dna_seg_from_file("../../z_accessory_files_dbs_and_executbles/AB270530.1.gbk")

# the comparison files below were prepared with the second segment as database
mycmp <- read_comparison_from_blast("../blast.txt",color_scheme="red_blue")
mycmp$col <- "darkred"
mycmp <- reverse(mycmp, side = 1)
# remove the negative values from the KA1pCAR3 associated coordinates since you are inversing the segment order
mycmp[,1:2] <- abs(mycmp[,1:2])
# prep the final files as lists (required by the plotting function)
cmplst <- list(mycmp)
dna_segs <- list(KA1pCAR3,gpm50ctg243)

# replace the bars with arrows
for(i in 1:length(dna_segs)){
  dna_segs[[i]]$gene_type <- "arrows"
}

# get the annotations by first setting the midpoint positions of the cds's and then defining the text
mid_pos<-list()
annot<-list()

my_colours<-c("purple","steelblue","lightblue","red2","springgreen3","white")

# perform some curation for the colouring
dna_segs[[1]]$product[35] <- dna_segs[[1]]$gene[35]

for(i in 1:length(dna_segs)){
  dna_segs[[i]]
  mid_pos[[i]]<-middle(dna_segs[[i]])
  annot[[i]]<-annotation(x1=mid_pos[[i]],text=dna_segs[[i]]$gene,rot=45)
  #change also the fill colors
  dna_segs[[i]]$col <- "black"
  dna_segs[[i]]$fill<-my_colours[6]
  dna_segs[[i]]$fill[grep("[C,c]ar",dna_segs[[i]]$gene)]<-my_colours[1]
  dna_segs[[i]]$fill[grep("cat|pca",dna_segs[[i]]$gene)]<-my_colours[2]
  dna_segs[[i]]$fill[grep("czc|cnr|ORF58|ORF59",dna_segs[[i]]$gene)]<-my_colours[3]
  dna_segs[[i]]$fill[grep("kinase|regulat|invertase|RNA|GGDEF|DNA-binding|receptor|regulator|ORF35",dna_segs[[i]]$product)]<-my_colours[4]
  dna_segs[[i]]$fill[unique(c(grep("[T,t]ranspos|[I,i]ntegrase|Tra[A-Z]|replication|plasmid|Vir[A-Z]|resolvase|Mobile|conjugal|integration|partit|recombi|mobA|pinE|transposase",dna_segs[[i]]$product),grep("[T,t]ranspos|[I,i]ntegrase|Tra[A-Z]|replication|IS5376|plasmid|Vir[A-Z]|resolvase|Mobile|conjugal|integration|partit|recombi|mobA|pinE|transposase|parA|Mob\\.|TnpA repr\\.|IS element",dna_segs[[i]]$gene)))]<-my_colours[5]
}

cairo_pdf(paste("Fig_comparison_locus_ctg243_vs_AB270530.pdf", sep = ""),height=2.5,width=20, onefile = T)
par(mar=c(4,4,8,4))
plot_gene_map(dna_segs, 
              dna_seg_labels=c("KA1 pCAR3","gpm 50 contig 243"),
              limit_to_longest_dna_seg=T,
              comparisons=cmplst, 
              dna_seg_line=rep("FALSE", length(dna_segs)),
              dna_seg_scale=rep(TRUE, length(dna_segs)), 
              scale=TRUE,
              n_scale_ticks=7,
              override_color_schemes=F,
              tree = NULL,
              tree_width = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              annotations=annot,
              annotation_height =3,
              annotation_cex=0.7,
              dna_seg_label_cex=1.5,
              xlims = list(c(18747,70000),NULL)
)


# prepare the colour key
my_colours<-c("purple","steelblue","lightblue","red2","springgreen3","white")
mylegend <- c("Carbazole degradation", "Catechol degradation", "Co/Zn/Cd efflux", "Regulatory", "Mobile", "Other/unknown")
plot.new()
par(mar=c(0,0,0,0), xpd = F)
legend(0,10,legend=mylegend[6:1], pch = 22, pt.bg=my_colours[c(6:1)],bty="n", pt.cex=2, cex=1, y.intersp = 7)
dev.off()





### prepare also an unannotated one
cairo_pdf(paste("Fig_comparison_locus_ctg243_vs_AB270530_stripped_from_annotations.pdf", sep = ""),height=1.3,width=15)
par(mar=c(4,4,8,4))
plot_gene_map(dna_segs, 
              dna_seg_labels=c("KA1 pCAR3","gpm 50 contig 243"),
              limit_to_longest_dna_seg=T,
              comparisons=cmplst, 
              dna_seg_line=rep("FALSE", length(dna_segs)),
              dna_seg_scale=rep(TRUE, length(dna_segs)), 
              scale=TRUE,
              n_scale_ticks=7,
              override_color_schemes=F,
              tree = NULL,
              tree_width = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              #annotations=annot,
              annotation_height =3,
              annotation_cex=0.7,
              dna_seg_label_cex=1.5,
              xlims = list(c(18747,70000),NULL)
)

dev.off()





