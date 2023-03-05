from pysam import VariantFile
import math
import csv
import re

CONSENSUS = "../01_consensus_sets/out/REH_consensus_3x.vcf"
CHROM_LENTHS = "data/chrom_lengths.tsv"
OUT_FILE="out/bins.tsv"

BIN_SIZE=1000000

bins = {}

def chromosome(str):
    # Possible formats: chr1, chrX, chr10
    # NOT chrUn_KI270742v1, N[chr14_GL000225v1_random:47528[
    # We exclude chrY since we know REH is from a female.
    m = re.search("chr([0-9X][0-9]?)$", str)
    if m:
        return m.group(1)
    

# Initialize bins
chrom_lengths = {}
with open(CHROM_LENTHS) as file:
    f = csv.reader(file, delimiter="\t")
    for row in f:
        chr = row[0]
        chrlen = int(row[2])
        chrom_lengths[chr] = chrlen
        bins[chromosome(chr)] = {}
        start = 1
        while start < chrlen:
            end = start + BIN_SIZE - 1
            end = min(end, chrlen)
            bins[chromosome(chr)][end] = 0
            start += BIN_SIZE

def process(vcf):
    for rec in vcf.fetch():
        chr = chromosome(rec.chrom)
        if chr:
            is_bnd = rec.info.get('SVTYPE') == "BND"
            pos = rec.pos
            chr_end = int(chrom_lengths["chr" + chr])

            # Optional: Ignore telomeres
            # is_telomere = pos < (BIN_SIZE * 2) or pos > chr_end - (BIN_SIZE * 2)
            
            if not is_bnd: # and not is_telomere:
                bin = min(math.ceil(pos / BIN_SIZE) * BIN_SIZE, chr_end)

                if not chr in bins:
                    bins[chr] = {}
                if not bin in bins[chr]:
                    bins[chr][bin] = 0

                bins[chr][bin] += 1

process(VariantFile(CONSENSUS))         

with open(OUT_FILE, 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    for chr in bins:
        for bin in bins[chr]:
            chr_name = "chr" + chr
            bin_name = chr_name + "_" + str(bin)
            start = bin - BIN_SIZE + 1  
            count = bins[chr][bin]
            row = [bin_name, chr_name, start, bin, count]
            tsv_output.writerow(row)