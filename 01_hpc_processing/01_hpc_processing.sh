#####################################################################
# Commands used on HPC cluster.
# Singularity containers are built from Docker images available at:
# https://hub.docker.com/u/themariya
#####################################################################

### ILLUMINA WGS PROCESSING ###
# Nextflow v21.10.6 26 with nf-core/sarek pipeline v2.7
nextflow run sarek/workflow/ \
	-profile uppmax \
	--project REDACTED \
	-c conf/sarek_upps.config \
	--custom_config_base sw/sarek/configs/ \
	--igenomes_ignore \
	--genome GRCh38 \
	--input REH_WGS_Illumina_TruSeqDnaPcrFree/fastq_data/ \
	--outDir REH_PCRfree_analysis/sarek_out/ \
	--tools haplotypecaller,snpeff,tiddit \
	--no_gatk_spark \
	--resume

### PACBIO WGS MAPPING ###
# pbmm2 v1.9.0. files.fofn contains: pt_004_001.ccsreads.fastq.gz, pt_004_002.ccsreads.fastq.gz, pt_004_003.ccsreads.fastq.gz
singularity exec pbmm2.sif pbmm2 align \
	    REH_WGS_PacBio_HiFi_analysis/alignment/hg38.fa \
	    REH_WGS_PacBio_HiFi_analysis/alignment/files.fofn \
	    --sort --min-concordance-perc 70.0 --min-length 50 --sample REH --preset CCS -j 24 --log-level INFO

### ONT WGS PROCESSING ###
# Porechop 0.2.4 for trimming
for i in REH_WGS_ONT_PromethION_Ultralong/fastq/*.fastq.gz
do
    file=`basename $i`
    singularity exec porechop.sif porechop \
		-i "$i" \
		-o "REH_WGS_ONT_PromethION_Ultralong/porechopped/$file" \
		--discard_middle \
		--threads 16

# Merging porechopped files
ls -d REH_WGS_ONT_PromethION_Ultralong/porechopped/*.gz | \
    xargs cat >> REH_WGS_ONT_PromethION_Ultralong/merged/porechopped.merged.fastq.gz

# Minimap v2.24-r1122
singularity exec minimap2.sif minimap2 \
	    -Lax map-ont \
	    hg38.fa \
		REH_WGS_ONT_PromethION_Ultralong/merged/porechopped.merged.fastq.gz > \
	    REH_WGS_ONT_PromethION_Ultralong/alignment/ont.sam

# Sorting and indexing with samtools v1.15.1
singularity exec samtools.sif samtools sort \
	    REH_WGS_ONT_PromethION_Ultralong/alignment/ont.sam \
	    -o REH_WGS_ONT_PromethION_Ultralong/alignment/ont.sorted.bam

singularity exec samtools.sif samtools index \
	    REH_WGS_ONT_PromethION_Ultralong/alignment/ont.sorted.bam \
	    REH_WGS_ONT_PromethION_Ultralong/alignment/ont.sorted.bam.bai

### DEPTH OF COVERAGE STATISTICS ###
# samtools 1.15.1 coverage
# Same command was run on fastq_data.recal.bam (Illumina), REH_pb_mapped.bam (PB), ont.sorted.bam (ONT)
singularity exec samtools.sif samtools coverage \
		REH_WGS_Illumina_TruSeqDnaPcrFree/fastq_data/Recalibrated/fastq_data.recal.bam > \
	    analysis/coverage/samtools.coverage.illumina.txt

# Copycat v9e21f79 (https://github.com/marianattestad/copycat) to bin coverage
/bin/bash copycat \
	  REH_WGS_Illumina_TruSeqDnaPcrFree/fastq_data/Recalibrated/fastq_data.recal.bam \
	  hg38bundle/Homo_sapiens_assembly38.fasta \
	  copycat.pcrfree

/bin/bash copycat \
	  REH_WGS_PacBio_HiFi_analysis/alignment/REH_pb_mapped.bam \
	  hg38.fa \
	  copycat.pb

/bin/bash copycat \
	  REH_WGS_ONT_PromethION_Ultralong/alignment/ont.sorted.bam \
	  hg38.fa \
	  copycat.ont

### SV CALLING ON LONG READS ###
# Sniffles v 2.0.6. Same command run on ont.sorted.bam and REH_pb_mapped.bam.
singularity exec sniffles.sif /usr/local/bin/sniffles \
	    --threads 16 \
	    --non-germline \
	    --tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed \
	    --input REH_WGS_ONT_PromethION_Ultralong/alignment/ont.sorted.bam \
	    --vcf ont.sniffles.vcf

### FUSION GENE CALLING ON SHORT-READ RNA-SEQ ###
#  nf-core/rnafusion pipeline v2.0.0 on Illumina RNA-seq (both REH and GM12878 samples)
nextflow run rnafusion/nf-core-rnafusion-2.0.0/workflow \
	-profile uppmax \
	--outdir results \
	--genomes_base rnafusion/rnafusion-references \
	--input samples.csv \
	--project REDACTED \
	--ensembl_version 103 \
	--all

### FUSION GENE CALLING ON ISOSEQ ###
# JAFFAL pipeline v2.1, run on analysis-pt_047_001-6714-flnc.fasta and analysis-pt_047_002-6715-flnc.fasta
singularity exec jaffa-2.1.sif bpipe run \
	    -n32 \
	    -p readLayout=single \
	    -p refBase=REH_RnaSeq_HiFi_PacBio_IsoSeq_SMRT/jaffa/reference \
	    /opt/JAFFA/JAFFAL.groovy \
	    REH_RnaSeq_HiFi_PacBio_IsoSeq_SMRT/jaffa/*.fasta