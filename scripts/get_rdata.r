args=commandArgs(trailingOnly=TRUE)
output=args[1]
sample_dir=args[2]
count_dir=args[3]


median.normalize.log <- function(matrix) {
    ## normalize the sequencing data
    log.matrix <- log2(matrix+1);

    ##browser();
    log.normalized.matrix <- apply(as.matrix(log.matrix), 2, function(x) { x-median(x)});

    return(log.normalized.matrix);
}
count_read=function(x){
 return(length(which(x>=0.1)))
 }
args=commandArgs(trailingOnly=T)
sample=read.table(paste(sample_dir,sep=""),header=T)
nsample=nrow(sample)
file.names=paste(count_dir,sample[,1],".count.txt",sep="")
all=matrix(0,26584,nsample)
for(i in 1:length(file.names)){
 dat=read.table(file.names[i],header=T,as.is=T,sep="\t")
 all[,i]=dat[,7]
 }
chrom=sapply(strsplit(dat$Chr,";"),"[",1)
rm.flag=which(chrom=="chrX" | (chrom=="chrY" | chrom=="chrM"))
if(length(rm.flag)>0){
all=all[-rm.flag,]
dat=dat[-rm.flag,]
}
cpm=matrix(0,26584-length(rm.flag),nsample)
total_count=apply(all,2,sum)
for(i in 1:ncol(all)){
 cpm[,i]=(all[,i]/(total_count[i]/1000000))
 }

ex_count=apply(cpm,1,count_read)
keep.flag=which(ex_count>=nsample/2)

all=cpm[keep.flag,]
dat=dat[keep.flag,]

gene_length=abs(dat$Length)

for(i in 1:ncol(all)){
all[,i]=(all[,i]/gene_length)/(total_count[i]/1000000)
}

#all=median.normalize.log(all)
library(preprocessCore)
all <- normalize.quantiles(all)

all.norm=matrix(0,nrow(all),ncol(all))
ex_mean=apply(all,1,mean,na.rm=T)
for(i in 1:ncol(all)){
all.norm[,i]=all[,i]-ex_mean
}


ex_sd=apply(all,1,sd,na.rm=T)
for(i in 1:ncol(all)){
all.norm[,i]=all.norm[,i]/ex_sd
}
#rm.flag=which(is.na(all.norm[,1]))
all.t=t(all)
ex.t=all.norm
ex=t(ex.t)
gnames=as.character(dat[,1])
gnames=gnames[flag.row]
#gnames=gnames[-rm.flag]
save(all.t,ex, gnames, file=output)
