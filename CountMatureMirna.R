#!/usr/local/bin/R Rscript
args <- commandArgs(TRUE)



ReadsFile= args[1];
ReadsFile= as.character(ReadsFile);
# The input ReadsFile from the pipeline

miRNAFile= args[2];
miRNAFile= as.character(miRNAFile);
# The input miRNA location file

outputFile= args[3];
outputFile= as.character(outputFile); 
# The output file

##################################################################

Reads=read.table(ReadsFile,header=FALSE);
miRNA=read.table(miRNAFile,header=FALSE);
ReadsIndex=as.character(Reads[,1]);
miRNAIndex=as.character(miRNA[,1]);


miRNA=cbind(miRNA,0);

for (i in 1:length(ReadsIndex))
{
	print (i)
	a=which(miRNAIndex==ReadsIndex[i]);
	if (length(a)==0) next;
	for (j in 1:length(a))
	{
		if (miRNA[a[j],3]<(Reads[i,2]+8)&&miRNA[a[j],4]>(Reads[i,2]+8))
		{
			miRNA[a[j],5]=miRNA[a[j],5]+1;		
		}	
	}
}

write.table(miRNA,file=outputFile, row.names=FALSE, col.names=FALSE);


#####################################################################