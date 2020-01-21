args=commandArgs(trailingOnly=T)
chr=args[1]
dat=read.table(paste("chr",chr,"_sig_genes_residQN_fdr10_p.txt",sep=""),header=T,sep="\t",as.is=T)
map=read.table("/ye/zaitlenlabstore/xuanyao/GE_cvp/scripts/mappability/hg19_GENCODE19/hg19_cross_mappability_strength_symmetric_mean_sorted_genename.txt",as.is=T)
gene1_all=NULL
gene2_all=NULL
 all=NULL
for(i in 1:nrow(dat)){
if(dat$num_reg[i]>1){
g2=strsplit(dat$reg_all,";")[[1]]
for(j in 1:length(g2)){
gene1=dat$genes[i]
 gene1_all=c(gene1_all,gene1)
gene2=g2[j]
 gene2_all=c(gene2_all,gene2)
 flag=which((map[,1]==gene1 &map[,2]==gene2) |(map[,1]==gene1 &map[,2]==gene2) )
if(length(flag)>0){all=c(all,map[flag,3])}
if(length(flag)==0){all=c(all,"NA")} 
 }
 }
 if(dat$num_reg[i]==1){
 gene1=dat$gene[i]
  gene1_all=c(gene1_all,gene1)
 gene2=dat$reg_all[i]
  gene2_all=c(gene2_all,gene2)
   flag=which(map[,1]==gene1 &map[,2]==gene2|(map[,1]==gene1 &map[,2]==gene2))
  if(length(flag)>0){all=c(all,map[flag,3])}
  if(length(flag)==0){all=c(all,"NA")}
 }
 }
 
 out=data.frame(gene1=gene1_all,gene2=gene2_all,mappability=all)
 write.table(out,paste("chr",chr,"_signal_mappability.txt",sep=""),quote=F,row.names=F,col.names=F)
