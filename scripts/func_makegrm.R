# Output a GCTA-format (zipped) GRM
# Input : fam = matix of sample IDs, out = prefix for file outputs, K = kinship matrix, M = # of markers

# Output id file
write.table( fam[,1:2] , file=paste(out,"grm.id",sep='.') , col.names=F , row.names=F , quote=F )

keep = which( lower.tri(K,diag=T) , arr.ind=T )
num = rep( M  , length(keep[,1]) )
grm = K[ keep ]
grm = grm[ order(keep[,1],keep[,2]) ]
writeBin( grm , con=paste(out,".grm.bin",sep="") , size=4 )
writeBin( num , con=paste(out,".grm.N.bin",sep="") , size=4 )
