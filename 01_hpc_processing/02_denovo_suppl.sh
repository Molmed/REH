#####################################################################
# Commands used on HPC cluster to generate supplementary data.
# Singularity containers are built from Docker images available at:
# https://hub.docker.com/u/themariya
#####################################################################

### PACBIO DE-NOVO ASSEMBLY ###
# hifiasm v0.16.1-r375
singularity exec hifiasm.sif hifiasm \
	    -o reh.pb.asm \
	    -t 48 \
	    pb.ccs.merged.fastq.gz

### ONT DE-NOVO ASSEMBLY ###
# NanoFilt v2.8.0
singularity exec nanofilt.sif NanoFilt \
	    -l 500 \
	    --headcrop 50 < \
	    porechopped.merged.fastq > \
	    trimmed.merged.porechopped.fastq

# flye v2.9.1-b1780
singularity exec flye.sif flye \
	    --nano-raw \
	    trimmed.merged.porechopped.fastq \
	    --genome-size 3.5g \
	    --threads 96 \
	    --scaffold \
	    --out-dir flye

### POLISHING FLYE WITH MEDAKA ###
# medaka v1.7.0
singularity exec medaka.sif medaka_consensus \
	    -i trimmed.merged.porechopped.fastq \
	    -d flye/00-assembly/draft_assembly.fasta \
	    -o medaka \
	    -t 32 \
	    -m r941_prom_hac_g507

### POLISHING FLYE WITH RACON ###
# Minimap v2.24-r1122
singularity exec minimap2.sif minimap2 \
	    -ax map-hifi \
	    flye/00-assembly/draft_assembly.fasta \
	    pb.ccs.merged.fastq.gz > \
	    racon/mapping/pb.overlaps.sam

### PHASING WITH HAPDU AND HAPDIFF ###
singularity minimap2.sif minimap2 \
	    -ax map-ont -t 30 \
	    medaka/consensus.fasta \
	    trimmed.merged.porechopped.fastq | \
        singularity exec samtools.sif samtools sort -@ 4 -m 4G > \
	    hapdup/ont_reads_to_flye_assembly.bam

singularity exec samtools.sif samtools index -@ 4 \
	    hapdup/ont_reads_to_flye_assembly.bam

singularity exec --bind ./hapdup \
	    hapdup.sif hapdup \
	    --assembly medaka/consensus.fasta \
	    --bam hapdup/ont_reads_to_flye_assembly.bam \
	    --out-dir hapdup/out \
	    -t 48 --rtype ont

singularity exec --bind ./hapdiff \
	    hapdiff.sif hapdiff.py \
	    --reference hapdiff/hg38.fa \
	    --pat hapdiff/hapdup_phased_1.fasta \
	    --mat hapdiff/hapdup_phased_2.fasta \
	    --out-dir ./hapdiff \
