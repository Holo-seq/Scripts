#!/usr/bin/R

args <- commandArgs(TRUE)

sample_id <- args[1]
sample_id <- as.character(sample_id)
infile_bed <- paste0(sample_id,".bed")
refexon_bed <- args[2]
refexon_bed <- as.character(refexon_bed)

refFlat_exons<- read.table(refexon_bed,sep='\t',header=FALSE)
sample_bed<- read.table(infile_bed,sep='\t',header=FALSE)
chr_vector<- refFlat_exons[,1]
chr_vector_unique <- unique(chr_vector)
chr_random_index<- grep("random",chr_vector_unique)
chr_vector_unique <- as.character(chr_vector_unique[-chr_random_index])
gene_list_name <- list()
gene_list_total<- list()
gene_list_rev_total<- list()
genelen_list<- list()
total<-0


for (i in 1:length(chr_vector_unique)) {
	chr<- chr_vector_unique[i]
	index_chr_exon<- which(refFlat_exons[,1]==chr)
	index_chr_bed<- which(sample_bed[,1]==chr)
	frame1_exon<- refFlat_exons[index_chr_exon,]
	frame1_bed<- sample_bed[index_chr_bed,]
	exoncount<- length(frame1_exon[,1])
	bedlength<- length(frame1_bed[,1])
	count_vector<- rep(0,300000000)
	if (bedlength==0){
		frame1_bed <- sample_bed[1:2,]
		frame1_bed[1:2,2]<- c(1)
		frame1_bed[1:2,3]<- c(2)
		bedlength<- length(frame1_bed[,1])
	}
	for (j in 1:bedlength) {
		start<- frame1_bed[j,2]
		end<- frame1_bed[j,3]
		count_vector[start:end]<- count_vector[start:end]+1
	}
	count_vector_matrix<- as.matrix(count_vector)

	combined<- matrix(nrow=exoncount,ncol=300000)

	geneid0 <- refFlat_exons[index_chr_exon[1],4]
	gene_strand0 <- refFlat_exons[index_chr_exon[1],5]
	genename0 <- refFlat_exons[index_chr_exon[1],6]
	index_chr_bed <- which(sample_bed[,1]==chr)
	column <- matrix(0)
	num <- 0
	genelen <-0
	gene_list <- list()
	gene_list_rev <- list()
	for (k in 1:exoncount) {
		geneid <- frame1_exon[k,4]
		strand <- refFlat_exons[index_chr_exon[1],5]
		genename <- as.character(frame1_exon[k,6])
		start<- frame1_exon[k,9]
		end<- frame1_exon[k,10]
		cut<- count_vector_matrix[start:end,1]
		cut<- as.matrix(cut)
		if(geneid==geneid0){
			column<- rbind(column,cut)
		}
		if(geneid!=geneid0){
			num<-num+1
			total<- total+1
			cat(total,'\t',file="log.txt",append=TRUE)
			gene_list_name[[total]]<- genename0
			genelen<- length(column[,1])
			combined[num,1:genelen] <- t(column)
			gene_list[[num]]<-c(genename0,t(column))
			gene_list_total[[total]]<-c(genename0,t(column))
			gene_list_rev[[num]]<- c(genename0,rev(t(column)))
			if(strand=="+"){
				gene_list_rev_total[[total]]<- c(genename0,rev(t(column)))
			}
			if(strand=="-"){
				gene_list_rev_total[[total]]<- c(genename0,t(column))
			}
			genelen_list[[total]]<-genelen
			column <- cut
			geneid0<- geneid
			genename0<- genename
			#save.image()
		}
	}
	num<- num+1
	total<- total+1
	cat(total,'\t',file="log.txt",append=TRUE)
	genename0 <- genename
	gene_list_name[[total]]<- genename0
	genelen<- length(column[,1])
	combined[num,1:genelen]<- t(column)
	gene_list[[num]]<-c(genename0,t(column))
	gene_list_total[[total]]<-c(genename0,t(column))
	gene_list_rev[[num]]<- c(genename0,rev(t(column)))
	gene_list_rev_total[[total]]<- c(genename0,rev(t(column)))
	genelen_list[[total]] <- genelen
	rm(combined)
	#save.image()
}
save.image()
###################################################
genelen_vector <- unlist(genelen_list)
gene_name_vector <- unlist(gene_list_name)
genelen_vector0 <- genelen_vector
gene_list_rev_total0 <- gene_list_rev_total

genenum_whole <- length(genelen_vector)
counted<- 0
for(i in 1:genenum_whole){
	add <- sum(as.numeric(gene_list_rev_total[[i]][-1]))
	counted<- counted+add
}
sum_forNormalize<- counted


##########################################################################
genenum_picked_vector <- c()
perc_weighted <- rep(0,100)
length_kb <- c(0,1,2,5)
for(i in 1:length(length_kb)){
	filename=paste0(sample_id,"_perc_",length_kb[i]*2-2,"-",length_kb[i]*2,"kb.csv")
	index_picked <- which(genelen_vector>=(length_kb[i]*2000-2000) & genelen_vector<(length_kb[i]*2000))
	gene_list_picked <- gene_list_rev_total[index_picked]
	genenum_picked <- length(index_picked)
	gene_matrix_picked <- matrix(nrow=genenum_picked,ncol=100)

	genenum_picked_vector <- c(genenum_picked_vector,genenum_picked)

	gene_list_picked_norm <- list()
	##############################

	for(i in 1:genenum_picked){
		gene_list_picked_norm[[i]] <- as.numeric(gene_list_picked[[i]][-1])/sum_forNormalize
	}
	#################################
	gene_list_picked_perc<- list()
	for(i in 1:genenum_picked){
		gene<- gene_list_picked_norm[[i]]
		genelen<- length(gene)
		gene_perc<- c()
		sum0<-0
		binlen=0.01*genelen
		for(j in 1:100){
			if(j < 100){
				sum_int<- sum(gene[1:floor(j*binlen)])
				sum_left <- (j*binlen-floor(j*binlen))*gene[(j*binlen)+1]
				sum1<- sum_int+sum_left

				sum=sum1-sum0
				gene_perc<- c(gene_perc,sum)
				sum0 <- sum1
			}
			if(j == 100){
				sum=sum(gene)-sum0
				gene_perc<- c(gene_perc,sum)
			}
		}
		gene_perc<- gene_perc/binlen    #normalize
		gene_list_picked_perc[[i]]<- gene_perc
	}

	##############################################
	gene_perc_base<- rep(0,100)
	for (i in 1:100){
		c<-0
		for(j in 1:genenum_picked){
			c <- c+gene_list_picked_perc[[j]][i]
		}
		gene_perc_base[i]<- c
	}

	gene_perc_base_smooth <- rep(0,100)
	for(i in 3:98){
		gene_perc_base_smooth[i]<- (gene_perc_base[i-2]+gene_perc_base[i-1]+gene_perc_base[i]+gene_perc_base[i+1]+gene_perc_base[i+2])/5
	}
	gene_perc_base_smooth[1] <- sum(gene_perc_base[1:3])/3
	gene_perc_base_smooth[2] <- sum(gene_perc_base[1:4])/4
	gene_perc_base_smooth[99] <- sum(gene_perc_base[97:100])/4
	gene_perc_base_smooth[100] <- sum(gene_perc_base[98:100])/3


	sum_base_smooth <- sum(gene_perc_base_smooth)
	sum_base <- sum(gene_perc_base)
	norm_gene_perc_base_smooth <- gene_perc_base_smooth/sum_base_smooth
	norm_gene_perc_base <- gene_perc_base/sum_base
	base <- c(1:100)
	frame_base_picked_norm <- data.frame(norm_gene_perc_base,norm_gene_perc_base_smooth,base)
	write.csv(frame_base_picked_norm,file=filename)
	
	perc_weighted <- perc_weighted+norm_gene_perc_base_smooth*(genenum_picked/genenum_whole)
}

###################################################
###################################################

	filename=paste0(sample_id,"_perc_gt6kb.csv")
	index_picked <- which(genelen_vector>=6000)
	gene_list_picked <- gene_list_rev_total[index_picked]
	genenum_picked <- length(index_picked)
	gene_matrix_picked <- matrix(nrow=genenum_picked,ncol=100)

	genenum_picked_vector <- c(genenum_picked_vector,genenum_picked)

	gene_list_picked_norm <- list()
	##############################

	for(i in 1:genenum_picked){
		gene_list_picked_norm[[i]] <- as.numeric(gene_list_picked[[i]][-1])/sum_forNormalize
	}
	#################################
	gene_list_picked_perc<- list()
	for(i in 1:genenum_picked){
		gene<- gene_list_picked_norm[[i]]
		genelen<- length(gene)
		gene_perc<- c()
		sum0<-0
		binlen=0.01*genelen
		for(j in 1:100){
			if(j < 100){
				sum_int<- sum(gene[1:floor(j*binlen)])
				sum_left <- (j*binlen-floor(j*binlen))*gene[(j*binlen)+1]
				sum1<- sum_int+sum_left

				sum=sum1-sum0
				gene_perc<- c(gene_perc,sum)
				sum0 <- sum1
			}
			if(j == 100){
				sum=sum(gene)-sum0
				gene_perc<- c(gene_perc,sum)
			}
		}
		gene_perc<- gene_perc/binlen    #normalize
		gene_list_picked_perc[[i]]<- gene_perc
	}

	##############################################

	for (i in 1:100){
		c<-0
		for(j in 1:genenum_picked){
			c <- c+gene_list_picked_perc[[j]][i]
		}
		gene_perc_base[i]<- c
	}

	gene_perc_base_smooth <- rep(0,100)
	for(i in 3:98){
		gene_perc_base_smooth[i]<- (gene_perc_base[i-2]+gene_perc_base[i-1]+gene_perc_base[i]+gene_perc_base[i+1]+gene_perc_base[i+2])/5
	}
	gene_perc_base_smooth[1] <- sum(gene_perc_base[1:3])/3
	gene_perc_base_smooth[2] <- sum(gene_perc_base[1:4])/4
	gene_perc_base_smooth[99] <- sum(gene_perc_base[97:100])/4
	gene_perc_base_smooth[100] <- sum(gene_perc_base[98:100])/3


	
	sum_base_smooth <- sum(gene_perc_base_smooth)
	sum_base <- sum(gene_perc_base)
	norm_gene_perc_base_smooth <- gene_perc_base_smooth/sum_base_smooth
	norm_gene_perc_base <- gene_perc_base/sum_base
	base <- c(1:100)
	frame_base_picked_norm <- data.frame(norm_gene_perc_base,norm_gene_perc_base_smooth,base)
	write.csv(frame_base_picked_norm,file=filename)
##############################
##############################
	filename=paste0(sample_id,"_perc_gt12kb.csv")
	index_picked <- which(genelen_vector>=12000)
	gene_list_picked <- gene_list_rev_total[index_picked]
	genenum_picked <- length(index_picked)
	gene_matrix_picked <- matrix(nrow=genenum_picked,ncol=100)

	genenum_picked_vector <- c(genenum_picked_vector,genenum_picked)

	gene_list_picked_norm <- list()
	##############################

	for(i in 1:genenum_picked){
		gene_list_picked_norm[[i]] <- as.numeric(gene_list_picked[[i]][-1])/sum_forNormalize
	}
	#################################
	gene_list_picked_perc<- list()
	for(i in 1:genenum_picked){
		gene<- gene_list_picked_norm[[i]]
		genelen<- length(gene)
		gene_perc<- c()
		sum0<-0
		binlen=0.01*genelen
		for(j in 1:100){
			if(j < 100){
				sum_int<- sum(gene[1:floor(j*binlen)])
				sum_left <- (j*binlen-floor(j*binlen))*gene[(j*binlen)+1]
				sum1<- sum_int+sum_left

				sum=sum1-sum0
				gene_perc<- c(gene_perc,sum)
				sum0 <- sum1
			}
			if(j == 100){
				sum=sum(gene)-sum0
				gene_perc<- c(gene_perc,sum)
			}
		}
		gene_perc<- gene_perc/binlen    #normalize
		gene_list_picked_perc[[i]]<- gene_perc
	}

	##############################################

	for (i in 1:100){
		c<-0
		for(j in 1:genenum_picked){
			c <- c+gene_list_picked_perc[[j]][i]
		}
		gene_perc_base[i]<- c
	}

	gene_perc_base_smooth <- rep(0,100)
	for(i in 3:98){
		gene_perc_base_smooth[i]<- (gene_perc_base[i-2]+gene_perc_base[i-1]+gene_perc_base[i]+gene_perc_base[i+1]+gene_perc_base[i+2])/5
	}
	gene_perc_base_smooth[1] <- sum(gene_perc_base[1:3])/3
	gene_perc_base_smooth[2] <- sum(gene_perc_base[1:4])/4
	gene_perc_base_smooth[99] <- sum(gene_perc_base[97:100])/4
	gene_perc_base_smooth[100] <- sum(gene_perc_base[98:100])/3


	
	sum_base_smooth <- sum(gene_perc_base_smooth)
	sum_base <- sum(gene_perc_base)
	norm_gene_perc_base_smooth <- gene_perc_base_smooth/sum_base_smooth
	norm_gene_perc_base <- gene_perc_base/sum_base
	base <- c(1:100)
	frame_base_picked_norm <- data.frame(norm_gene_perc_base,norm_gene_perc_base_smooth,base)
	write.csv(frame_base_picked_norm,file=filename)
#########################################################
#######################


