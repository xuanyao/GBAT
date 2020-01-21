logfile=$1
cutoff=$(head -n1 ./log_residQN_h2g | cut -d" " -f2)
for i in {1..22}; do Rscript summarize_sig_genes_fromP.r $i $cutoff; done
for i in {1..22}; do Rscript summarize_gene_pair.r $i $cutoff; done

