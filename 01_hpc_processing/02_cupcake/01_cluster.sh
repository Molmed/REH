#! /bin/bash

###############
# DESCRIPTION #
###############
# Performs isoseq3 clustering on FLNC BAM file
#
#############
# VARIABLES #
#############
# ISOSEQ3_BIN=[path to `isoseq3` executable]
# FLNC_BAM=[path to bam file]
# OUTPUT_BAM=[output name, e.g. sample.clustered.bam]
# NUM_THREADS

$ISOSEQ3_BIN cluster $FLNC_BAM $OUTPUT_BAM --verbose --use-qvs --num-threads $NUM_THREADS
