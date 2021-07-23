# load the files
mycoverage <- read.table("coverage_new.tsv", header = T, row.names = 1, sep = "\t")
mymembership <- read.table("../dastool/dastoolout_DASTool_scaffolds2bin.txt", row.names = 1, sep = "\t")
# merge the coverage and membership tables
myfin_pre <- merge(mymembership, mycoverage, by = "row.names", all = T)
colnames(myfin_pre)[2] <- "bin"

myfin <- myfin_pre[,2:ncol(myfin_pre)]
row.names(myfin) <- myfin_pre$Row.names
myfin$bin <- factor(myfin$bin)

myfinfin <- myfin[,-grep("MD",colnames(myfin))]


# prepare the plot for all bins
cairo_pdf("complete_dataset_line_plot.pdf", height = 5, width = 6)
par(mar = c(5,4,4,10))
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = 'Dark2')[-4])(n = length(levels(myfinfin$bin))) # the nature colors
plot(x = 2:ncol(myfinfin)-1,y = myfinfin[1,2:ncol(myfinfin)] +1, type = "l", col = mycolors[which(levels(myfinfin[1,"bin"])%in%myfinfin[1,1])], frame.plot = F, ylim = c(1,10000), log = "y", xlab = "sample", ylab = "mean coverage")
par(xpd = T)
for(i in 2:nrow(myfinfin)){
  lines(x = 2:ncol(myfinfin)-1,y = myfinfin[i,2:ncol(myfinfin)]+1, col = mycolors[which(levels(myfinfin[1,"bin"])%in%myfinfin[i,1])], frame.plot = F)
}
legend(x = 6.2, y = 10000, legend = levels(myfinfin$bin), pch = "-", pt.cex = 2, col = mycolors[1:length(levels(myfinfin$bin))], bty = "n", ncol = 1)
dev.off()


# do the same with ggplot2... I chose palette on purpose to prevent plotting all bins but plot as many as there are colours... you can change that on will
library(reshape2)
myfinfinfin <- myfinfin
myfinfinfin$contig <- row.names(myfinfinfin)
myfinfin_long <- melt(myfinfinfin, id = c("contig","bin"))
myfinfin_long$value <- myfinfin_long$value + 1

library(ggplot2)
p <- ggplot(data = myfinfin_long) +
  geom_line(aes(x = variable, y = value, colour = bin, group = contig)) +
  scale_y_continuous(trans='log10') +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + 
  ylab("mean coverage")

cairo_pdf("complete_dataset_line_plot_ggplot_not_everything.pdf", height = 6, width = 6)
print(p)
dev.off()

# plot only the Novosphingobium and Sphingomonas contigs
myfinfinfin_long_nov_sph <- myfinfin_long[which(myfinfin_long$bin%in%c("groopm_bin_50","groopm_bin_30")),]

# set the colors
library(plyr)
myfinfinfin_long_nov_sph$CLRs <- mapvalues(myfinfinfin_long_nov_sph$bin,from = c("groopm_bin_30","groopm_bin_50"), to = c(rgb(0, 0, 138, max = 255, alpha = 125), rgb(138, 0, 0, max = 255, alpha = 30)))

myfinfinfin_long_nov_sph$CLRs <- as.character(myfinfinfin_long_nov_sph$CLRs)
# set the colours of the contigs of interest
myfinfinfin_long_nov_sph$CLRs[which(myfinfinfin_long_nov_sph$contig == "opera_contig_243")] <- rgb(138, 0, 0, max = 255, alpha = 255) 
myfinfinfin_long_nov_sph$CLRs[which(myfinfinfin_long_nov_sph$contig%in% c("opera_contig_224","opera_contig_395"))] <- rgb(128,128,128, max = 255, alpha = 255) 
myfinfinfin_long_nov_sph$CLRs[which(myfinfinfin_long_nov_sph$contig%in% "opera_contig_2")] <- rgb(255,69,0, max = 255, alpha = 255) 


p <- ggplot(data = myfinfinfin_long_nov_sph) +
  geom_line(aes(x = variable, y = value, colour = CLRs, group = contig)) +
  scale_y_continuous(trans='log10') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + 
  ylab("mean coverage") +
  scale_color_identity()

cairo_pdf("Nov_Sph_dataset_line_plot_ggplot.pdf", height = 6, width = 6)
print(p)
dev.off()



#### the same as above for percent coverages ----
myfinfin <- myfin[,-grep("MD",colnames(myfin))]

myfinfin_perc_pre <- 100*vegan::decostand(myfinfin[,2:ncol(myfinfin)], MARGIN = 2, method = "total")
myfinfin_perc <- data.frame(bin = myfinfin[,1],myfinfin_perc_pre)

# plot as many as possible bins with ggplot2... I chose palette on purpose to prevent plotting all bins but plot as many as there are colours... you can change that on will
cairo_pdf("perc_complete_dataset_line_plot.pdf", height = 5, width = 6)
par(mar = c(5,4,4,10))
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = 'Dark2')[-4])(n = length(levels(myfinfin_perc$bin))) # the nature colors
plot(x = 2:ncol(myfinfin_perc)-1,y = myfinfin_perc[1,2:ncol(myfinfin_perc)] +1, type = "l", col = mycolors[which(levels(myfinfin_perc[1,"bin"])%in%myfinfin_perc[1,1])], frame.plot = F, ylim = c(1,50), log = "y", xlab = "sample", ylab = "mean coverage (%)")
par(xpd = T)
for(i in 2:nrow(myfinfin_perc)){
  lines(x = 2:ncol(myfinfin_perc)-1,y = myfinfin_perc[i,2:ncol(myfinfin_perc)]+1, col = mycolors[which(levels(myfinfin_perc[1,"bin"])%in%myfinfin_perc[i,1])], frame.plot = F)
}
legend(x = 6.2, y = 50, legend = levels(myfinfin_perc$bin), pch = "-", pt.cex = 2, col = mycolors[1:length(levels(myfinfin_perc$bin))], bty = "n", ncol = 1)
dev.off()



library(reshape2)
myfinfin_percfin <- myfinfin_perc
myfinfin_percfin$contig <- row.names(myfinfin_percfin)
myfinfin_perc_long <- melt(myfinfin_percfin, id = c("contig","bin"))
myfinfin_perc_long$value <- myfinfin_perc_long$value + 0.0001

library(ggplot2)
p <- ggplot(data = myfinfin_perc_long) +
  geom_line(aes(x = variable, y = value, colour = bin, group = contig)) +
  scale_y_continuous(trans='log10') +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + 
  ylab("mean relative coverage (%)")

cairo_pdf("perc_complete_dataset_line_plot_ggplot_not_everything.pdf", height = 6, width = 6)
print(p)
dev.off()

# plot only the Novosphingobium and Sphingomonas contigs

# set the colors
myfinfin_percfin_long_nov_sph <- myfinfin_perc_long[which(myfinfin_perc_long$bin%in%c("groopm_bin_50","groopm_bin_30")),]
library(plyr)
myfinfin_percfin_long_nov_sph$CLRs <- mapvalues(myfinfin_percfin_long_nov_sph$bin,from = c("groopm_bin_30","groopm_bin_50"), to = c(rgb(0, 0, 138, max = 255, alpha = 125), rgb(138, 0, 0, max = 255, alpha = 30)))

myfinfin_percfin_long_nov_sph$CLRs <- as.character(myfinfin_percfin_long_nov_sph$CLRs)
# set the colours of the contigs of interest
myfinfin_percfin_long_nov_sph$CLRs[which(myfinfin_percfin_long_nov_sph$contig == "opera_contig_243")] <- rgb(138, 0, 0, max = 255, alpha = 255) 
myfinfin_percfin_long_nov_sph$CLRs[which(myfinfin_percfin_long_nov_sph$contig%in% c("opera_contig_224","opera_contig_395"))] <- rgb(128,128,128, max = 255, alpha = 255) 

myfinfin_percfin_long_nov_sph$CLRs[which(myfinfin_percfin_long_nov_sph$contig%in% "opera_contig_2")] <- rgb(255,69,0, max = 255, alpha = 255) 

p <- ggplot(data = myfinfin_percfin_long_nov_sph) +
  geom_line(aes(x = variable, y = value, colour = CLRs, group = contig)) +
  scale_y_continuous(trans='log10') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + 
  ylab("mean relative coverage (%)") +
  scale_color_identity()

cairo_pdf("perc_Nov_Sph_dataset_line_plot_ggplot.pdf", height = 6, width = 6)
print(p)
dev.off()
