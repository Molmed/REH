#! /bin/bash

###############
# DESCRIPTION #
###############
# Performs mapping and sorting on high-quality isoform file (<prefix>hq.fasta.gz)
#
#############
# VARIABLES #
#############
# MINIMAP2_BIN=[path to `minimap2` executable]
# REF_GENOME=[path to e.g. GRCh38.p13.genome.fa]
# HQ_ISOFORMS=[path to <prefix>hq.fasta.gz]
# NUM_THREADS
# OUTPUT_DIR
# OUTPUT_PREFIX

# Map
$MINIMAP2_BIN -ax splice:hq -uf --secondary=no -t$NUM_THREADS $REF_GENOME $HQ_ISOFORMS > "$OUTPUT_DIR/$OUTPUT_PREFIX.sam"

# Sort
sort -k 3,3 -k 4,4n "$OUTPUT_DIR/$OUTPUT_PREFIX.sam" > "$OUTPUT_DIR/$OUTPUT_PREFIX.sorted.sam"
