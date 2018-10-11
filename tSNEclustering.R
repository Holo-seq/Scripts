#Rscript --argument $ReadCountMatirx.csv $cut1 $cut2 $perplexity $groupsize1 $groupsize2 $groupsize3 
#say we have three sample groups in this analysis

library(edgeR)
library(Rtsne)
args <- commandArgs(trailingOnly= TRUE)
names(args) <- c("ReadCountMatirx.csv","cut1","cut2","perplexity","groupsize1","groupsize2","groupsize3")

groupsize1 <- args["groupsize1"]
groupsize2 <- args["groupsize2"]
groupsize3 <- args["groupsize3"]
cut1 <- args["cut1"]
cut2 <- args["cut2"]
perp <- args["perplexity"]
wholeframe_tagcount_info <- read.csv(args["ReadCountMatirx.csv"],header=TRUE)

sampleid <- colnames(wholeframe_tagcount_info)
sampleid <- sampleid[,4:length(wholeframe_tagcount_info[1,])]

wholeframe_tagcount_info_picked <- wholeframe_tagcount_info[rowSums(wholeframe_tagcount_info[,4:length(wholeframe_tagcount_info[1,])]>=cut1)>=cut2,] ###pick expressed genes for rpkm and further calculation according to your demand
wholeframe_rpkm_picked <- rpkm(wholeframe_tagcount_info_picked[,4:length(wholeframe_tagcount_info[1,])],gene.length = as.numeric(wholeframe_tagcount_info_picked[,3]),log = FALSE)
wholeframe_rpkm_picked_pca <- prcomp(t(wholeframe_rpkm_picked),scale. = FALSE)
plot(wholeframe_rpkm_picked_pca)
sdev_sum <- sum(wholeframe_rpkm_picked_pca$sdev)
sdev_cut <- which.min(wholeframe_rpkm_picked_pca$sdev>sdev_sum*0.85)

par(pin=c(4,3.5))
for_tsne <- wholeframe_rpkm_picked_pca$x[,1:sdev_cut]
groupsize <- args[5:length(args)]
tsne_out <- Rtsne(for_tsne,dims = 2,perplexity = perp,pca = FALSE) ##perplexity could be adjusted to proper value
color_vec = sample(colours(),length(groupsize))
png("tsneClust.png",width = 5*300,height = 5*300,res = 300,pointsize = 8) 
plot(tsne_out$Y,main = "t-SNE clustering plot",xlab = "t-sne1",ylab = "t-sne2")
for(i in 1:length(groupsize)){
	points(tsne_out$Y[(sum(groupsize[1:i])-groupsize[i]+1):sum(groupsize[1:i]),1:2],col=color_vec[i],pch=16,cex=1) #group 1, A is the index of the last sample of group 1
}
#legend(LegendLoc,c(group 1,group 2,group 3),pch = c(16,16,16),col=color_vec,cex = 1) #LegendLoc is the location of legend
dev.off()

