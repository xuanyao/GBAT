args=commandArgs(trailingOnly=T)
chr=args[1]
sub=args[2]
pr_dir=args[3]
output=args[4]
rdata_dir=args[5]
cvblup_dir=args[6]
tool_dir=args[7]
load(rdata_dir)
source(paste(tool_dir,'make_smartsva.R',sep=""))
exp_genes=gnames
library(phenix)
g=read.table(paste(cvblup_dir,sep=""),header=T,check.names=F)
pr=read.table(pr_dir,as.is=T,header=T,sep="\t")
cov1=read.table("/project2/gilad/xuanyao/GE_TWAS/DGN/covariates/Biological_and_hidden_factors.txt",header=T,sep="\t")
cov2=read.table("/project2/gilad/xuanyao/GE_TWAS/DGN/covariates/Technical_factors.txt",header=T,sep="\t")
samples=read.table("/project2/gilad/xuanyao/GE_TWAS/DGN/samples.txt",header=T)
keep.flag=which(is.element(cov1[,1],samples[,1]))
cov.keep=cbind(cov1[keep.flag,-1],cov2[keep.flag,-1])
cov=as.matrix(cov.keep[order(cov.keep[,1]),])

coln=NULL
allcor=NULL
allcoef=NULL
for(i in 1:ncol(g)){
		 if(!is.na(g[1,i]) & (pr[i,2]>0.0995)){
			 tryCatch({
			 g_sva=make_sva(ex,100,g[,i])
			 g_norm=(g[,i]-mean(g[,i]))/sd(g[,i])
			 if(!is.na(g_sva[1,1])){
			 out=as.matrix(cbind(g_sva[,1:20],cov))
			 rlm=lm(ex~out)
			 rQN=quantnorm(as.matrix(resid(rlm)))} else{rQN=ex}
			 test=lm(rQN~g_norm)
			 allcor=cbind(allcor,sapply(summary(test),function(x) x$coefficients[2,4]))
			 allcoef=cbind(allcoef,sapply(summary(test),function(x) x$coefficients[2,1]))
			 coln=c(coln,colnames(g)[i])

			  }, error=function(e){})
			 #test=apply(ex_norm,2,cor,g[,i],method='spearman')
		}
}
 write.table(allcor,output,row.names=F,col.names=coln)
#write.table(allcoef,paste(output,sep=""),row.names=F,col.names=coln)
 

