import pathlib
from pathlib import Path


chromosomes = [f'chr{str(x)}' for x in range(1, 23)] + ['chrX']
#chromosomes = ['chr9', 'chr21', 'chr22']
hg38 = '/proj/ngi2016005/FU/FU164/analysis/REH_DeepVariant_Hybrid/reference/Homo_sapiens_assembly38.fasta'
loci = '/proj/ngi2016005/nobackup/private/reh_baf_figures/external_data/loci/G1000_loci_hg38'

with open('vcf_files.txt', 'r') as bf:
	samples = {x.split()[0]: {'vcf': x.split()[1]} for x in bf}


localrules: fix_loci

rule all:
	input:
		expand('allelecounter_vcf/{sample}.tsv', sample=samples),
		expand('figures/{sample}_baf_logr.pdf', sample=samples)

rule fix_loci:
	input:
		loci = lambda wildcards: f'{loci}_{wildcards.chr}.txt',

	output:
		'loci/{chr}.loci'

	run:
		with open(input.loci, 'r') as inp, open(output[0], 'w') as outp:
			for line in inp:
				outp.write(f'chr{line}')
		

rule split_vcf:
	input:
		vcf = lambda wildcards: samples[wildcards.sample]['vcf'],

	output:
		'vcf_by_chr/{sample}/{sample}_{chr}.vcf.gz'

	container:
		'container/allelecount_ggplot.0.1.sif'

	shell:
		"""
		tabix -h {input.vcf} {wildcards.chr} | bgzip > {output}
		tabix -p vcf {output}
		"""

rule allele_counter_vcf:
	input:
		vcf = 'vcf_by_chr/{sample}/{sample}_{chr}.vcf.gz',
		loci = 'loci/{chr}.loci'
	output:
		'allelecounter_vcf/{sample}_{chr}.tsv'

	container:
		'container/allelecount_ggplot.0.1.sif'

	script:
		'alleleCounter.R'

rule collect_counts:
	input:
		[f'allelecounter_vcf/{{sample}}_{chr}.tsv' for chr in chromosomes]

	output:
		'allelecounter_vcf/{sample}.tsv'

	shell:
		"""
		cat {input} > {output}
		"""

rule plot_baf_logr:
	input:
		'allelecounter_vcf/{sample}.tsv'

	output:
		'figures/{sample}_baf_logr.pdf'

	container:
		'container/allelecount_ggplot.0.1.sif'

	script:
		'baf_plot.R'
