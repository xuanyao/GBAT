# Output a GCTA-format (zipped) GRM
# Input : fam = matix of sample IDs, out = prefix for file outputs, K = kinship matrix, M = # of markers

# Output id file
write.table( fam[,1:2] , file=paste(out,"grm.id",sep='.') , col.names=F , row.names=F , quote=F )

keep = which( lower.tri(K,diag=T) , arr.ind=T )
non_mis = rep( M  , length(keep[,1]) )
# Bind row/column IDs with matrix
keep = cbind( keep , non_mis , K[ keep ] )
# Re-order in increasing order
keep = keep[ order(keep[,1],keep[,2]) , ]
write.table( keep , file=gzfile(paste(out,"grm.gz",sep='.')) , col.names=F , row.names=F , quote=F )
