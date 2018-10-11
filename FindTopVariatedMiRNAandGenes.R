#Rscript --argument $ReadCountMatrix_miRNA.csv $ReadCountMatrix_mRNA.csv $cut1 $cut2 $cut3 $cut4

library(edgeR)
library(gplots)
args <- commandArgs(trailingOnly= TRUE)
names(args) <- c("ReadCountMatrix_miRNA.csv","ReadCountMatrix_mRNA.csv","cut1","cut2","cut3","cut4")

ReadCountMatrix_miRNA.csv <- args["ReadCountMatrix_miRNA.csv"]
ReadCountMatrix_mRNA.csv <- args["ReadCountMatrix_mRNA.csv"]
cut1 <- args["cut1"]
cut2 <- args["cut2"]
cut3 <- args["cut3"]
cut4 <- args["cut4"]

wholeframe_count_mRNA_info <- read.csv(ReadCountMatrix_mRNA.csv,header=TRUE)
wholeframe_count_miRNA_info <- read.csv(ReadCountMatrix_miRNA.csv,header=TRUE)


#################find top variated miRNA 
miRNA_id <- as.character(wholeframe_count_miRNA_info[,2])
vector1 <- colnames(wholeframe_count_miRNA_info)
vector1 <- vector1[4:length(wholeframe_count_miRNA_info[1,])]
#######
#######
count_miRNA_matrix <- as.matrix(wholeframe_count_miRNA_info[,4:length(wholeframe_count_miRNA_info[1,])])
miRNA_sum <- apply(count_miRNA_matrix,2,sum)
miRNA_rpm_matrix <- matrix(ncol=length(vector1),nrow=length(miRNA_id))
colnames(miRNA_rpm_matrix) <- vector1
for(i in 1:length(vector1)){
  miRNA_rpm_matrix[,i] <- (count_miRNA_matrix[,i]/miRNA_sum[i])*1000000
}
miRNA_rpm_frame <- data.frame(miRNA_id,miRNA_rpm_matrix)
#################################################################################
#################################################################################
miRNA_rpm_frame_exp <- miRNA_rpm_frame[rowSums(count_miRNA_matrix[,1:length(vector1)]>=cut1)>=cut2,]		####pick expressed miRNA according to your demand
wholeframe_logrpm <- log2(miRNA_rpm_frame_exp[,2:(length(vector1)+1)]+1)
miRNAid_exp <- as.character(miRNA_rpm_frame_exp[,1])
wholeframe_logrpm_id <- data.frame(miRNAid_exp,wholeframe_logrpm)
mark <- c(1:length(miRNAid_exp))
wholeframe_logrpm_marked <- data.frame(mark,wholeframe_logrpm)
variance <- apply(wholeframe_logrpm_marked[,2:(length(vector1)+1)],1,var)
mean <- apply(wholeframe_logrpm_marked[,2:(length(vector1)+1)],1,mean)
for_zscore_frame <- data.frame(mark,miRNAid_exp,mean,variance)
z_scoring <- mean/variance
for_zscore_frame <- data.frame(for_zscore_frame,z_scoring)
for_zscore_frame_beneath1 <- for_zscore_frame[which(z_scoring<1),]
for_zscore_frame_picked_VM <- for_zscore_frame_beneath1[which(for_zscore_frame_beneath1[,3]>=1),]

for_zscore_frame_picked_VM_ordered <- for_zscore_frame_picked_VM[order(for_zscore_frame_picked_VM[,5],decreasing = FALSE),]
index_VM_ordered <- for_zscore_frame_picked_VM_ordered[,1]
wholeframe_logrpm_id_VM_ordered <- wholeframe_logrpm_id[index_VM_ordered,]


mir26a5p <- wholeframe_logrpm_id_VM_ordered[which(as.character(wholeframe_logrpm_id_VM_ordered[,1])=="hsa-miR-26a-5p"),]
forheatmap_miRNA <- rbind(mir26a5p,wholeframe_logrpm_id_VM_ordered[1:99,])

##########################################################################################################################################
##########################################################################################################################################

#################find top variated mRNA
vector2 <- colnames(wholeframe_count_mRNA_info)
vector2 <- vector2[4:length(wholeframe_count_mRNA_info[1,])]
gene_id <- as.character(wholeframe_count_mRNA_info[,2])
gene_len <- as.numeric(wholeframe_count_mRNA_info[,3])

#######
#######
count_mRNA_matrix <- as.matrix(wholeframe_count_mRNA_info[,4:length(wholeframe_count_mRNA_info[1,])])
mRNA_rpkm_matrix <- rpkm(count_mRNA_matrix,gene.length=gene_len,log=FALSE)
mRNA_rpkm_frame <- data.frame(gene_id,mRNA_rpkm_matrix)
#################################################################################
#################################################################################
mRNA_rpkm_frame_exp <- mRNA_rpkm_frame[rowSums(count_mRNA_matrix[,1:length(vector2)]>=cut3)>=cut4,]		####pick expressed mRNA according to your demand
wholeframe_logrpkm <- log2(mRNA_rpkm_frame_exp[,2:(length(vector2)+1)]+1)
geneid_exp <- as.character(mRNA_rpkm_frame_exp[,1])
wholeframe_logrpkm_id <- data.frame(geneid_exp,wholeframe_logrpkm)
mark <- c(1:length(geneid_exp))
wholeframe_logrpkm_marked <- data.frame(mark,wholeframe_logrpkm)
variance <- apply(wholeframe_logrpkm_marked[,2:(length(vector2)+1)],1,var)
mean <- apply(wholeframe_logrpkm_marked[,2:(length(vector2)+1)],1,mean)
for_zscore_frame <- data.frame(mark,geneid_exp,mean,variance)
z_scoring <- mean/variance
for_zscore_frame <- data.frame(for_zscore_frame,z_scoring)
for_zscore_frame_beneath1 <- for_zscore_frame[which(z_scoring<1),]
for_zscore_frame_picked_VM <- for_zscore_frame_beneath1[which(for_zscore_frame_beneath1[,3]>=1),]

for_zscore_frame_picked_VM_ordered <- for_zscore_frame_picked_VM[order(for_zscore_frame_picked_VM[,5],decreasing = FALSE),]
index_VM_ordered <- for_zscore_frame_picked_VM_ordered[,1]
wholeframe_logrpkm_id_VM_ordered <- wholeframe_logrpkm_id[index_VM_ordered,]

forheatmap_mRNA <- wholeframe_logrpkm_id_VM_ordered[1:500,]

colname1 <- c("id",vector1)

colnames(forheatmap_miRNA) <- colname1
colnames(forheatmap_mRNA) <- colname1
frame_mRNAmiRNAtop600 <- rbind(forheatmap_miRNA[,2:(length(vector1)+1)],forheatmap_mRNA[,2:(length(vector1)+1)])
id_miRNAmRNA <- c(as.character(forheatmap_miRNA[,1]),as.character(forheatmap_mRNA[,1]))
png("HierClusHeatmap_top600VariatedmiRNAmRNA.png",width = 5*300,height = 5*300,res = 300,pointsize = 8) 
heatmap.2(as.matrix(frame_mRNAmiRNAtop600),col = colorpanel(250,low = "purple",high = "yellow",mid = "black"),
                     scale = "row",key = TRUE,keysize = 1,symkey = FALSE,density.info = "none",trace = "none",
                     srtCol = 90,labCol = NULL,margins = c(10,15),labRow = id_miRNAmRNA)
dev.off()
