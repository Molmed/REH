import csv
import matplotlib.pyplot as plt

CONSENSUS = "out/trimmed.tsv"

afs = []

def append_af(r):
    dr = r.split(':')[3]
    reads = dr.split(',')
    ref_reads = float(reads[0])
    variant_reads = float(reads[1])
    total_reads = ref_reads + variant_reads
    af = variant_reads / total_reads
    if (af > 0):
        afs.append(af)

def plot(data):
    # Create a histogram
    plt.hist(data, bins=20, edgecolor='black')

    # Add labels and title
    plt.xlabel('AF (reads supporting SV / read depth at breakpoint)')
    plt.ylabel('Frequency')
    plt.title('Allele fractions (AF) in 3x consensus SV callset')

    # Save the histogram as a PNG image
    plt.savefig('out/histogram.png')

    # Show the histogram (optional)
    plt.show()


with open(CONSENSUS, 'r') as tsv_file:
    tsv_reader = csv.reader(tsv_file, delimiter='\t')
    for record in tsv_reader:
        ilmn = record[9]
        ont = record[10]
        pb = record[11]

        append_af(ilmn)
        append_af(ont)
        append_af(pb)

    plot(afs)
