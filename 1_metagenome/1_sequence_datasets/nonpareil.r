library(Nonpareil)

cairo_pdf(paste("nonpareil_results.pdf", sep = ""),height = 10, width = 3.5, onefile = T)
par(mfrow = c(4,1))
mysites <- c("SRR7135606","SRR7135607","SRR7135608","SRR7135612")
# prep the table
mycovtbl <- array(dim=c(length(mysites),2))
row.names(mycovtbl) <- mysites
colnames(mycovtbl) <- c("NP_coverage","Nd")

for(mysample in mysites){
  mat <- read.table(paste("nonpareil_",mysample,"_1.fastq.npo", sep = ""))
  if(which(mysites%in%mysample) == 1){
    myylab <- "Estimated mean coverage"
  } else {
    myylab <- ""
  }
  np <- Nonpareil.curve(paste("nonpareil_",mysample,"_1.fastq.npo", sep = ""), col = "darkblue", main = paste(gsub("_.+","",mysample), sep = ""), ylab = myylab)
  lines(rep(np@x.adj[length(np@x.adj)],2),c(0,np@y.cov[length(np@y.cov)]), col = "red")
  par(xpd = T)
  text(np@x.adj[length(np@x.adj)],np@y.cov[length(np@y.cov)], labels = paste("est. cov.", round(np@y.cov[length(np@y.cov)],3), "/ Nd", round(np@diversity,2)), cex = 1.4, pos = 3)
  text(np@x.adj[length(np@x.adj)],np@y.cov[length(np@y.cov)], labels = paste(formatC(mat$V1[length(mat$V1)], format = "e", digits = 1), "read pairs /",formatC(np@x.adj[length(np@x.adj)], format = "e", digits = 1),"bp"), cex = 1.2, pos = 1)
  par(xpd = F)
  mycovtbl[mysample,] <- c(np@y.cov[length(np@y.cov)],np@diversity)
}

dev.off()


write.table(mycovtbl,"nonpareil_results.txt", col.names = NA, quote = F, sep = "\t")


