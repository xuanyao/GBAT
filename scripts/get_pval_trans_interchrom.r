args=commandArgs(trailingOnly=T)
chr=as.numeric(args[1])
allp=NULL
allrm=NULL
gene_pos=read.table("./cor/gene_pos.txt",header=T,as.is=T,fill=T,sep="\t")
pseudo=read.table("/ye/zaitlenlabstore/xuanyao/GE_cvp/scripts/mappability/hg19_GENCODE19/pseudogenes.txt",as.is=T)
rm.flag=which(is.element(gene_pos$gene,pseudo[,1]))
gene_pos=gene_pos[-rm.flag,]
  for(i in 1:10){
if(file.exists(paste("./cor/chr",chr,"_fix_p_rmext_20sva_sub",i,".txt",sep=""))){
 dat=read.table(paste("./cor/chr",chr,"_fix_p_rmext_20sva_sub",i,".txt",sep=""),header=T,as.is=T)
gene.names=colnames(dat)
dat=data.frame(dat[-rm.flag,])
for(nc in 1:ncol(dat)){
gene=gene.names[nc]
flag=which(gene_pos$gene==gene)
if(chr==1){
start=0
}else{start=max(which(gene_pos$chr==paste("chr",chr-1,sep="")))}

if(chr==22){
end=nrow(dat)}else{end=min(which(gene_pos$chr==paste("chr",chr+1,sep="")))}
 rm.flag=start:end
allrm=c(allrm,length(rm.flag))
 allp=c(allp,as.numeric(dat[-rm.flag,nc]))
}
}
#if(file.exists(paste("chr",chr,"_fix_p_permutation_sub",i,".txt",sep=""))){
# dat=read.table(paste("chr",chr,"_fix_p_permutation_sub",i,".txt",sep=""),header=T,as.is=T)
#for(nc in 1:ncol(dat)){
#gene=gene.names[nc]
#flag=which(gene_pos$gene==gene)
#start=max(0,flag-50)
#end=min(flag+50,nrow(dat))
# rm.flag=start:end
# allperm=c(allperm,as.numeric(dat[-rm.flag,nc]))
#}
#}
 }
write.table(allp,paste("./cor/chr",chr,"_allp_interchrom.txt",sep=""),quote=F,row.names=F,col.names=F)
#write.table(allperm,paste("chr",chr,"_perm_allp.txt",sep=""),quote=F,row.names=F,col.names=F)
