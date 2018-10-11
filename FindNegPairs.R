#Rscript --argument $ReadCountMatrix_mRNA.csv $ReadCountMatrix_miRNA.csv $miRNA_mRNA_StrongProvedPairs.txt 

library(edgeR)
args <- commandArgs(trailingOnly= TRUE)
names(args) <- c("ReadCountMatrix_mRNA.csv","ReadCountMatrix_miRNA.csv","miRNA_mRNA_StrongProvedPairs.txt")

ReadCountMatrix_mRNA.csv <- args["ReadCountMatrix_mRNA.csv"]
ReadCountMatrix_miRNA.csv <- args["ReadCountMatrix_miRNA.csv"]
miRNA_mRNA_StrongProvedPairs.txt <- args["miRNA_mRNA_StrongProvedPairs.txt"]

wholeframe_count_mRNA_info <- read.csv(ReadCountMatrix_mRNA.csv,header=TRUE)
wholeframe_count_miRNA_info <- read.csv(ReadCountMatrix_miRNA.csv,header=TRUE)
miRNA_id <- as.character(wholeframe_count_miRNA_info[,2])
gene_id <- as.character(wholeframe_count_mRNA_info[,2])
gene_len <- as.numeric(wholeframe_count_mRNA_info[,3])
vector1 <- colnames(wholeframe_count_miRNA_info)
vector1 <- vector1[4:length(wholeframe_count_miRNA_info[1,])]
vector2 <- colnames(wholeframe_count_mRNA_info)
vector2 <- vector2[4:length(wholeframe_count_mRNA_info[1,])]

count_miRNA_matrix <- as.matrix(wholeframe_count_miRNA_info[,4:length(wholeframe_count_miRNA_info[1,])])
count_mRNA_matrix <- as.matrix(wholeframe_count_mRNA_info[,4:length(wholeframe_count_mRNA_info[1,])])
miRNA_rpm_matrix <- matrix(nrow=length(miRNA_id),ncol=length(vector1))
mRNA_rpkm_matrix <- rpkm(count_mRNA_matrix,gene.length=gene_len,log=FALSE)
colnames(count_miRNA_matrix) <- vector1
colnames(count_mRNA_matrix) <- vector2
colnames(miRNA_rpm_matrix) <- vector1
colnames(mRNA_rpkm_matrix) <- vector2
for(i in 1:length(vector1)){
  count_vector <- count_miRNA_matrix[,i]
  count_sum <- sum(count_vector)
  miRNA_rpm_matrix[,i] <- (count_vector*1000000)/count_sum
}

frame_miRNA_rpm <- data.frame(miRNA_id,miRNA_rpm_matrix)
frame_mRNA_rpkm <- data.frame(gene_id,mRNA_rpkm_matrix)


index_miRNA_read_exp <- which(rowSums(count_miRNA_matrix[,1:length(vector1)]=="0")<=round(length(vector1)*2/3)) #1/3 of sample size
index_mRNA_read_exp <- which(rowSums(count_mRNA_matrix[,1:length(vector2)]=="0")<=round(length(vector1)*2/3)) #1/3


frame_miRNA_rpm_exp <- frame_miRNA_rpm[index_miRNA_read_exp,]
frame_mRNA_rpkm_exp <- frame_mRNA_rpkm[index_mRNA_read_exp,]

pairs_table <- read.table(miRNA_mRNA_StrongProvedPairs.txt,sep='\t',header=TRUE)
index_miRNA_pair_picked <- c()
index_mRNA_pair_picked <- c() 
for (i in 1:length(pairs_table[,1])){
  miRNA_pairid <- as.character(pairs_table[i,1])
  mRNA_pairid <- as.character(pairs_table[i,2])
  index_miRNA_pair <- which(as.character(frame_miRNA_rpm_exp[,1]) == miRNA_pairid)
  index_mRNA_pair <- which(as.character(frame_mRNA_rpkm_exp[,1]) == mRNA_pairid)
  if(length(index_miRNA_pair)!=0 & length(index_mRNA_pair)!=0){
    if(length(index_mRNA_pair)>1){
      index_miRNA_pair <- rep(index_miRNA_pair,length(index_mRNA_pair))
    }
    index_miRNA_pair_picked <- c(index_miRNA_pair_picked,index_miRNA_pair)
    index_mRNA_pair_picked <- c(index_mRNA_pair_picked,index_mRNA_pair)
  }
}

miRNA_exp_frame_mRNA_rpkm_exp_picked <- data.frame(frame_miRNA_rpm_exp[index_miRNA_pair_picked,],frame_mRNA_rpkm_exp[index_mRNA_pair_picked,])


#####################

cor_rpkm <- c()
cor_rpkm_p <- c()

for (i in 1:length(miRNA_exp_mRNA_exp_norm_frame_picked[,1])){
  add2 <- cor.test(as.numeric(miRNA_exp_frame_mRNA_rpkm_exp_picked[i,2:(length(vector1)+1)]),as.numeric(miRNA_exp_frame_mRNA_rpkm_exp_picked[i,(length(vector1)+3):((length(vector1)+1)*2)]),method="spearman")
  cor_rpkm <- c(cor_rpkm,add2$estimate)
  cor_rpkm_p <- c(cor_rpkm_p,add2$p.value)
}
frame_pairs_cor_p <- data.frame(miRNA_exp_frame_mRNA_rpkm_exp_picked,cor_rpkm,cor_rpkm_p)
index1 <- which(cor_rpkm<0)
index2 <- which(cor_rpkm_p<0.05)
index_SigNegPairs <- intersect(index1,index2)
frame_pars_cor_p_SigNeg <- frame_pars_cor_p[index_SigNegPairs,]
write.csv(frame_pars_cor_p_SigNeg,file="frame_pars_cor_p_SigNeg.csv")
