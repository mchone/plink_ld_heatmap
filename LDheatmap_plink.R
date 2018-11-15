#!/public/home/hwu/mcbin/R325/bin/R
args<-commandArgs(T)
library(grid)
library(LDheatmap)
library(genetics)

ag <- unlist(strsplit(args[1], "\\."))
outfile <- paste0(ag[1],"_pruned.tif")

#snp <- data.frame((read.csv(args[1],header = T, stringsAsFactors = F, sep = "\t")))
#snp <- read.table(args[1],header = F, stringsAsFactors = F, sep = "\t")
dist <- read.csv(args[2],header = F, stringsAsFactors = F, sep = "\t")
dis <- as.vector(t(dist))
dis <- as.numeric(dis)

#print(length(colnames(snp)))
#print(dist[1,])

#num<-ncol(snp)        # SNPdata is a data.frame which looks like the CEUSNP above.
#colnames(snp)<-dist[1,] # Name the SNPdata columns with the 2nd column of map file(produced by vcftools).

#class(snp[,1]) # You will get "factor".

#for(i in 1:num){
#  snp[,i]<-as.genotype(snp[,i]) # convert the columns of the original data frame into  genotype objects
#}

#class(snp[,1])  #Here you can see the class becomes "genotype factor".Then you can put the SNPdata into LDheatmap().

tiff(file= outfile,width=72, height=85, units='mm',res=1600,compression='lzw')
#setEPS()
#postscript("rb.eps")
snp <- as.matrix(read.table(args[1],header=T))
MyHeatmap <- LDheatmap(snp, dis, LDmeasure = "r",
                       title = "Pairwise LD in r^2",
                       add.map = TRUE,
                       color = colorRampPalette(c("red","yellow","white"))(250),
                       name = "myLDgrob",
                       add.key = TRUE,
                       #text = TRUE,
                       flip = TRUE
)
grid.edit(gPath("myLDgrob", "geneMap", "diagonal"), gp=gpar(lwd =0.1, col ="black"))
grid.edit(gPath("myLDgrob", "geneMap", "segments"), gp=gpar(lwd =0.1, col ="black"))

dev.off()
