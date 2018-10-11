#Rscript --argument $ReadCountMatirx.csv $RNA $cut1 $cut2

library(edgeR)
library(gplots)
args <- commandArgs(trailingOnly= TRUE)
names(args) <- c("ReadCountMatirx.csv","RNA","cut1","cut2")

RNA <- args["RNA"]
cut1 <- args["cut1"]
cut2 <- args["cut2"]
wholeframe_tagcount_info <- read.csv(args["ReadCountMatirx.csv"], header=TRUE)
sampleid <- colnames(wholeframe_tagcount_info)
sampleid <- sampleid[4:length(wholeframe_tagcount_info[1,])]

wholeframe_tagcount_info_picked <- wholeframe_tagcount_info[rowSums(wholeframe_tagcount_info[,4:length(wholeframe_tagcount_info[1,])]>=cut1)>=cut2,] ###pick expressed genes for rpkm calculation according to your demand
if(RNA == "mRNA"){
	wholeframe_exp_picked <- rpkm(wholeframe_tagcount_info_picked[,4:length(wholeframe_tagcount_info[1,])],gene.length = as.numeric(wholeframe_tagcount_info_picked[,3]),log = FALSE)
}
if(RNA == "miRNA"){
	count_miRNA_matrix <- wholeframe_tagcount_info[,4:length(wholeframe_tagcount_info[1,])]
	wholeframe_exp_picked <- matrix(0,ncol=length(sampleid),nrow=wholeframe_tagcount_info_picked[,1])
	for(i in 1:length(sampleid)){
		count_vector <- count_miRNA_matrix[,i]
		count_sum <- sum(count_vector)
		wholeframe_exp_picked[,i] <- (count_vector*1000000)/count_sum
	}
}

max_vector <- apply(wholeframe_exp_picked,1,max)
wholeframe_exp_picked_norm <- wholeframe_exp_picked/max_vector
geneid_picked <- as.character(wholeframe_tagcount_info_picked[,2])
png("HierClusHeatmap.png",width = 5*300,height = 5*300,res = 300,pointsize = 8) 
heatmap.2(wholeframe_exp_picked_norm,col = colorpanel(250,low = "green",high = "red",mid = "black"),scale = "none",key = TRUE,keysize = 1,symkey = FALSE,density.info = "none",trace = "none",srtCol = 90,labRow=geneid_picked,labCol = sampleid,margins = c(10,15))
dev.off()


