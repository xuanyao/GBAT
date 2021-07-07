median.normalize.log <- function(matrix) {
    ## normalize the sequencing data
    log.matrix <- log2(matrix+1);

    ##browser();
    log.normalized.matrix <- apply(as.matrix(log.matrix), 2, function(x) { x-median(x)});

    return(log.normalized.matrix);
}
count_read=function(x){
 return(length(which(x>=0.5)))
 }

args=commandArgs(trailingOnly=TRUE)
output=args[1]
sample_dir=args[2]
count_dir=args[3]

samples=read.table(sample_dir,header=T)
N=nrow(samples)
all=matrix(0,26584,N)
for(i in 1:N){
 dat=read.table(paste(count_dir,samples[i,1],".count.txt",sep=""),header=T,as.is=T,sep="\t")
 all[,i]=dat[,7]
 }
chrom=sapply(strsplit(dat$Chr,";"),"[",1)
rm.flag=which(chrom=="chrX" | (chrom=="chrY" | chrom=="chrM"))
if(length(rm.flag)>0){
all=all[-rm.flag,]
dat=dat[-rm.flag,]
}
cpm=matrix(0,26584-length(rm.flag),N)
for(i in 1:ncol(all)){
 cpm[,i]=(all[,i]/(total_count[i]/1000000))
 }

ex_count=apply(cpm,1,count_read)
keep.flag=which(ex_count>=N/2)

all=all[keep.flag,]
dat=dat[keep.flag,]

gene_length=abs(dat$Length)
rpk=matrix(0,26584-length(rm.flag),N)

for(i in 1:ncol(all)){
rpk[,i]=(all[,i]/gene_length)
}
total_rpk=apply(rpk,2,sum)

for(i in 1:ncol(all)){
all[,i]=rpk[,i]/(total_rpk/1000000)
}

#all=median.normalize.log(all)
library(preprocessCore)
library(phenix)
#all <- normalize.quantiles(as.matrix(all))
all <- normalize.quantiles(as.matrix(all))
all.t=t(all)
all2=quantnorm(all.t)
all=t(all2)

#rm.flag=which(is.na(all.norm[,1]))
ex=t(all)
gnames=as.character(dat[,1])
#gnames=gnames[-rm.flag]
save(ex, gnames, file=output)
