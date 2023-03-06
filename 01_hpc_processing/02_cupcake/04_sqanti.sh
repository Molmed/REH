#! /bin/bash

###############
# DESCRIPTION #
###############
# Runs SQANTI3 fusion transcript classification.
#
#############
# VARIABLES #
#############
# SQANTI3_QC_PY=[path to `sqanti3_qc.py` executable]
# FUSION_GFF=[path to <prefix>.gff output from step 03]
# GENCODE_GTF=[path to gene annotation file gencode.<version>.annotation.gtf]
# REF_GENOME=[path to GRCh38.<version>.genome.fa]

$SQANTI3_QC_PY --gtf $FUSION_GFF $GENCODE_GTF $REF_GENOME --is_fusion
