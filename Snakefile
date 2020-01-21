import pandas as pd

configfile:"config.yaml"
samples=pd.read_table(config["samples"]).set_index("sample",drop=False)

CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22".split()
SUBS="1 2 3 4 5 6 7 8 9 10".split()

rule all:
	#input: expand(config['trans_DIR']+'cor_chr{chr}_sub{sub}.txt',chr=CHRS,sub=SUBS)
	input: expand('chr{chr}_sig_fdr10.txt',chr=CHRS)
	shell: 'rm *.nosex; rm *log; rm chr*sub*.map; rm chr*sub*.ped'

rule sig:
        input:'log_residQN_h2g'
        output:'chr{chr}_sig_fdr10.txt'
        params: chr="{chr}"
	shell: 'bash '+config['tool_DIR']+'run_sig.sh {input}'

rule qval:
        input:expand(config['trans_DIR']+'chr{chr}_allp_h2g_filtered_interchrom.txt',chr=CHRS)
        output: 'log_residQN_h2g'
        shell: 'Rscript '+config['tool_DIR']+'cal_qval.r '+config['trans_DIR']+'> {output}'

rule pval:
	input:expand(config['trans_DIR']+'cor_chr{{chr}}_sub{sub}.txt',sub=SUBS)
	output:config['trans_DIR']+'chr{chr}_allp_h2g_filtered_interchrom.txt'
	params: chr="{chr}"
	shell:config['Rscript_DIR']+'Rscript '+config['tool_DIR']+'get_pval_trans_perm_interchrom.r {params.chr} '+config['trans_DIR']

rule cor:
	input: config['trans_DIR']+'pearsonR_chr{chr}_sub{sub}.txt'
	output: config['trans_DIR']+'cor_chr{chr}_sub{sub}.txt'
	params: chr="{chr}",sub="{sub}"
	shell: config['Rscript_DIR']+'Rscript '+config['tool_DIR']+'cal_cor.r {params.chr} {params.sub} {input} {output} /project2/xuanyao/data/DGN/count/dgn_2qt.rdata '+config['trans_DIR']+'chr{params.chr}_sub{params.sub}.txt '+config['tool_DIR']
	
rule pearsonR:
	input: config['trans_DIR']+'chr{chr}_sub{sub}.txt'
	output: config['trans_DIR']+'pearsonR_chr{chr}_sub{sub}.txt'
	shell: config['Rscript_DIR']+'Rscript '+config['tool_DIR']+'cal_pearsonR.r {input} {output} /project2/xuanyao/data/DGN/count/dgn_2qt.rdata'
	

rule cvBLUP:
	input: config['trans_DIR']+'sva.txt'
	output: config['trans_DIR']+'chr{chr}_sub{sub}.txt'
	params: chr="{chr}",sub="{sub}"
	shell: config['Rscript_DIR']+'Rscript '+config['tool_DIR']+'cal_goos_fix_chr.r {params.chr} {params.sub} '+config['trans_DIR']+'gene_pos.txt '+config['genotype_DIR']+' /project2/xuanyao/data/DGN/count/dgn_2qt.rdata '+config['trans_DIR']+'sva.txt '+config['genpc_DIR']+' '+config['tool_DIR']+' {output} '+config['plink_samples']

rule sva:
	input: '/project2/xuanyao/data/DGN/count/dgn_2qt.rdata'
	output: config['trans_DIR']+'sva.txt'
	shell:config['Rscript_DIR']+'Rscript '+config['tool_DIR']+'get_gen_pos.r {input} '+config['count_DIR']+" "+config['trans_DIR']+';'+config['Rscript_DIR']+'Rscript '+config['tool_DIR']+'sva.r {input} {output} '+config['num_sv']+' '+config['tool_DIR']

rule rdata:
        input: expand(config['count_DIR']+'{sample}.count.txt',sample=samples.index)
        output: config['trans_DIR']+'my.rdata'
        shell: config['Rscript_DIR']+'Rscript '+config['tool_DIR']+'get_rdata.r {output} '+config['samples']+' '+config['count_DIR']


rule count:
        input: config['BAM_OUT']+'{sample}.bam'
        output: config['count_DIR']+'{sample}.count.txt'
        shell: config['tool_DIR']+'featureCounts -t exon -g gene_id -a '+config['gtf_DIR']+'hg19.mappability.filtered.gtf -o {output} {input}'

rule filter:
        input: config['BAM_DIR']+'{sample}.bam'
        output: config['BAM_OUT']+'{sample}.bam'
        shell: config['python_DIR']+'python '+config['tool_DIR']+'pysam_dgn_filter.py --bfile {input} --outdir '+config['BAM_OUT']
