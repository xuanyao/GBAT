args=commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]
rdata_dir=args[3]
rsquared=function(p,a){
 r2=1-(sum((a-p)^2)/sum((a-mean(a))^2))
 return(r2)
}


load(rdata_dir)
exp_genes=gnames

dat=read.table(input,header=T,as.is=T,check.names=F)
#dat.names=read.table(paste("chr",chr,"_genes_175.txt",sep=""),sep="\t")

goos_cor=NULL
	for(i in 1:ncol(dat)){
	 namei=colnames(dat)[i]
	 ex.flag=which(exp_genes==namei)
	 goos_cor=c(goos_cor,cor(dat[,i],ex[,ex.flag]))
	 }
out=data.frame(genes=colnames(dat),rsquared_goos_sva=goos_cor)
write.table(out,output,row.names=F,quote=F,sep="\t")

