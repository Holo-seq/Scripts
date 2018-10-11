#Rscript --argument $RNA_loci.csv $SE_loci.csv

library(edgeR)
args <- commandArgs(trailingOnly= TRUE)
names(args) <- c("RNA_loci.csv","SE_loci.csv")

RNA_loci.csv <- args["RNA_loci.csv"]
SE_loci.csv <- args["SE_loci.csv"]

RNA_loci <- read.csv(RNA_loci.csv,header=TRUE)
SE_loci <- read.csv(SE_loci.csv,header=TRUE)


RNA_loci <- read.csv(<RNA_loci.csv>,header=TRUE)
colnames(RNA_loci) <- c("RNA_id","RNA_chr","RNA_start")
RNA_id <- as.character(RNA_loci[,1])
enhancer_loci <- read.csv(<SE_loci.csv>,header=TRUE)
enhancer_loci <- enhancer_loci[,1:4]
colnames(enhancer_loci) <- c("eh_id","eh_chr","eh_start","eh_end")
enhancer_id <- as.character(enhancer_loci[,1])
formode <- data.frame(RNA_loci[1,],enhancer_loci[1,])
RNA_eh_pair_min_frame <- formode[1,]
colnames(RNA_eh_pair_min_frame) <- c("RNA_id","RNA_chr","RNA_start","eh_id","eh_chr","eh_start","eh_end")
distance_vector <- c()
for (i in 1:length(RNA_id)){
	chr <- as.character(RNA_loci[i,2])
	enhancer_temp <- enhancer_loci[which(as.character(enhancer_loci[,2])==chr),]
	index_min <- which(abs((enhancer_temp[,3]+(enhancer_temp[,4]-enhancer_temp[,3])/2)- RNA_loci[i,3])== min(abs(enhancer_temp[,3]+(enhancer_temp[,4]-enhancer_temp[,3])/2- RNA_loci[i,3])))
	combined_row <- c(t(RNA_loci[i,]),t(enhancer_temp[index_min,]))
	RNA_eh_pair_min_frame <- rbind(RNA_eh_pair_min_frame,combined_row)
}
RNA_eh_pair_min_frame <- RNA_eh_pair_min_frame[-1,]
distance <- as.numeric(RNA_eh_pair_min_frame[,3])-((as.numeric(RNA_eh_pair_min_frame[,6])+(as.numeric(RNA_eh_pair_min_frame[,7])-as.numeric(RNA_eh_pair_min_frame[,6]))/2))
RNA_eh_pair_min_frame$distance <- distance
index <- which(abs(distance) < 100000)
RNA_eh_pair_within100kb_frame <- RNA_eh_pair_min_frame[index,]
write.csv(RNA_eh_pair_within100kb_frame,file="RNA_eh_pair_within100kb_frame.csv")

#############
#############
