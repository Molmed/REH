import pandas
import re
from gene_info import get_chromosome, set_chromosome, save_db

DATA_DIR="../data/fusion_callsets"
OUT_DIR="out"
GENE_WHITELIST="data/gene_whitelist.txt"
# NOTE: this list was generated AFTER data analysis
MANUALLY_CONFIRMED="data/manually_confirmed.txt"

# SHORT-READ RESULTS
SHORT_READ_DIR = DATA_DIR + "/short_read"
GM12878 = SHORT_READ_DIR + "/GM12878.fusionreport.txt"
ILLUMINA_ALL = SHORT_READ_DIR + "/illumina.all.txt"
ILLUMINA_FILTERED = SHORT_READ_DIR + "/illumina.filtered.csv"
ARRIBA = SHORT_READ_DIR + "/REH.arriba.fusions.tsv"
FUSIONCATCHER = SHORT_READ_DIR + "/REH.fusioncatcher.fusion-genes.txt"
PIZZLY = SHORT_READ_DIR + "/REH.pizzly.txt"
SQUID = SHORT_READ_DIR + "/REH.squid.fusions.annotated.txt"
STARFUSION = SHORT_READ_DIR + "/REH.starfusion.abridged.tsv"


# First get unique list of all fusions in format "GENE1--GENE2"
# illumina_all is already in this format, and illumina_filtered is a subset
with open(ILLUMINA_ALL) as f:
    illumina_all = f.read().splitlines()
illumina_filtered = pandas.read_csv(ILLUMINA_FILTERED)
arriba = pandas.read_csv(ARRIBA, sep='\t')
fusioncatcher = pandas.read_csv(FUSIONCATCHER, sep='\t')
pizzly = pandas.read_csv(PIZZLY, sep='\t')
squid = pandas.read_csv(SQUID, sep='\t')
starfusion = pandas.read_csv(STARFUSION, sep='\t')
# GM12878 normals to filter against
with open(GM12878) as f:
    gm12878 = f.read().splitlines()
# whitelist with known ALL-related genes
with open(GENE_WHITELIST) as f:
    whitelist = f.read().splitlines()
# the fusions that have been manually confirmed
with open(MANUALLY_CONFIRMED) as f:
   manually_confirmed_list = f.read().splitlines()

# LONG-READ RESULTS
LONG_READ_DIR = DATA_DIR + "/long_read"
JAFFA = LONG_READ_DIR + "/jaffa_results.csv"
CUPCAKE_STD = LONG_READ_DIR + "/cupcake.std.csv"
CUPCAKE_LONG = LONG_READ_DIR + "/cupcake.long.csv"

jaffa=pandas.read_csv(JAFFA)
cupcake_std=pandas.read_csv(CUPCAKE_STD)
cupcake_long=pandas.read_csv(CUPCAKE_LONG)

# OUTPUT COLUMNS
FUSION_GENE_SEPARATOR="--"
NOTES_SEPARATOR=","
YES = 'X'
MANUALLY_CONFIRMED_COL="Manually confirmed"
PASSED_AUTOMATED_FILTERING_COL="Passed automated filtering"
FUSION_GENE_COL="Fusion gene"
LEFT_GENE_COL="Left gene"
RIGHT_GENE_COL="Right gene"
LEFT_CHROM_COL="Left chromosome"
RIGHT_CHROM_COL="Right chromosome"
SHORT_READ_TOOLS_COL="Short-read tools"
LONG_READ_TOOLS_COL="Long-read tools"
TOTAL_TOOLS_COL="Total tools"
ILLUMINA_MAX_READS_COL="Illumina max supporting reads"
ISOSEQ_MAX_READS_COL="IsoSeq max supporting reads"
ISOSEQ_MAX_VARIANTS_COL="IsoSeq max variants"
KNOWN_ALL_COL="Known ALL gene"
OCCURS_IN_HEALTHY_COL="Occurs in healthy B cells"
NOTES_COL="Notes"

# PREPROCESS ARRIBA
# Genes can be of this format: "LINC02404(3758),AC090049.1(17041)"
# We remove the numbers in parantheses and convert to a list.
def preprocess_arriba_genes(col_name, new_col_name):
    lmb = lambda row: re.sub("\(\d+\)","",row[col_name]).split(",")
    arriba[new_col_name] = arriba.apply(lmb, axis=1)
preprocess_arriba_genes("#gene1", LEFT_GENE_COL)
preprocess_arriba_genes("gene2", RIGHT_GENE_COL)

# If a gene is a combination of two genes, split into two rows.
arriba = arriba.explode(LEFT_GENE_COL)
arriba = arriba.explode(RIGHT_GENE_COL)
arriba[FUSION_GENE_COL] = arriba[[LEFT_GENE_COL, RIGHT_GENE_COL]].apply(FUSION_GENE_SEPARATOR.join, axis=1)

# Use the max of split_reads1 and split_reads2 for each row
arriba["split_reads"] = arriba[["split_reads1", "split_reads2"]].apply(max, axis=1)

# Get the chromosomes
lmb = lambda row: "chr" + row.split(":")[0]
arriba[LEFT_CHROM_COL]=arriba['breakpoint1'].map(lmb)
arriba[RIGHT_CHROM_COL]=arriba['breakpoint2'].map(lmb)

# Misc notes
arriba[NOTES_COL] = arriba[["type", "reading_frame"]].apply(NOTES_SEPARATOR.join, axis=1)

# PREPROCESS FUSION CATCHER
fusioncatcher[FUSION_GENE_COL] = fusioncatcher[ \
    ["Gene_1_symbol(5end_fusion_partner)", "Gene_2_symbol(3end_fusion_partner)"] \
        ].apply(FUSION_GENE_SEPARATOR.join, axis=1)

# Get chromosomes
lmb = lambda row: "chr" + row.split(":")[0]
fusioncatcher[LEFT_CHROM_COL]=fusioncatcher["Fusion_point_for_gene_1(5end_fusion_partner)"].map(lmb)
fusioncatcher[RIGHT_CHROM_COL]=fusioncatcher["Fusion_point_for_gene_2(3end_fusion_partner)"].map(lmb)

# Misc notes
fusioncatcher["Fusion_description"] = fusioncatcher["Fusion_description"].astype(str)
fusioncatcher[NOTES_COL] = fusioncatcher[["Predicted_effect","Fusion_description"]].apply(NOTES_SEPARATOR.join, axis=1)

# PREPROCESS PIZZLY
pizzly[FUSION_GENE_COL] = pizzly[["geneA.name", "geneB.name"]].apply(FUSION_GENE_SEPARATOR.join, axis=1)

# PREPROCESS SQUID
# Replace separator then convert fused genes column into list
lmb = lambda row: row["FusedGenes"].replace(':',FUSION_GENE_SEPARATOR).split(",")
squid[FUSION_GENE_COL] = squid.apply(lmb, axis=1)

# If a gene is a combination of several genes, split into multiple rows.
squid = squid.explode(FUSION_GENE_COL)

# PREPROCESS STARFUSION
starfusion = starfusion.rename(columns={"#FusionName": FUSION_GENE_COL})

# Get chromosomes
lmb = lambda row: "chr" + row.split(":")[0]
starfusion[LEFT_CHROM_COL]=starfusion["LeftBreakpoint"].map(lmb)
starfusion[RIGHT_CHROM_COL]=starfusion["RightBreakpoint"].map(lmb)

# Misc notes
starfusion[NOTES_COL] = starfusion[["SpliceType", "annots"]].apply(NOTES_SEPARATOR.join, axis=1)

# TODO: ADD NOTES FROM LONG-READ
# PREPROCESS JAFFA
# JAFFA output is in the format "GENE1:GENE2"
jaffa[FUSION_GENE_COL]=jaffa['fusion genes'].map(
    lambda gene:FUSION_GENE_SEPARATOR.join(gene.split(":")))

# PREPROCESS CUPCAKE
def preprocess_cupcake_df(cupcake_df):
    # Cupcake output contains left and right gene in separate columns
    # Convert all to strings to avoid problems before fusing
    cupcake_df["LeftGeneID"] = cupcake_df["LeftGeneID"].astype(str)
    cupcake_df["RightGeneID"] = cupcake_df["RightGeneID"].astype(str)
    # If the gene id is a combination of two, just use the first one
    # TODO: explode them to two records instead
    lmb = lambda gene:gene.split("_")[0]
    cupcake_df["LeftGeneID"]=cupcake_df["LeftGeneID"].map(lmb)
    cupcake_df["RightGeneID"]=cupcake_df["RightGeneID"].map(lmb)
    cupcake_df[FUSION_GENE_COL] = cupcake_df[["LeftGeneID", "RightGeneID"]].apply(FUSION_GENE_SEPARATOR.join, axis=1)
    # Get the chromosomes
    lmb = lambda row:row.split(":")[0]
    cupcake_df[LEFT_CHROM_COL]=cupcake_df['LeftBreakpoint'].map(lmb)
    cupcake_df[RIGHT_CHROM_COL]=cupcake_df['RightBreakpoint'].map(lmb)

preprocess_cupcake_df(cupcake_std)
preprocess_cupcake_df(cupcake_long)

# Merge cupcake sets
cupcake_std["Read length"] = "standard"
cupcake_long["Read length"] = "long"
cupcake=pandas.concat([cupcake_std, cupcake_long])

# JOIN ALL FUSIONS IN COMPREHENSIVE LIST
fusion_genes = list(illumina_all) \
    + list(jaffa[FUSION_GENE_COL]) \
    + list(cupcake_std[FUSION_GENE_COL]) \
    + list(cupcake_long[FUSION_GENE_COL])

# GET UNIQUE SET
fusion_genes = list(set(fusion_genes))
print("{} unfiltered fusion candidates found.".format(len(fusion_genes)))

# DUMP TO NEW DATA FRAME
columns = [
    MANUALLY_CONFIRMED_COL, 
    PASSED_AUTOMATED_FILTERING_COL,
    FUSION_GENE_COL,
    LEFT_GENE_COL,
    RIGHT_GENE_COL,
    LEFT_CHROM_COL,
    RIGHT_CHROM_COL,
    SHORT_READ_TOOLS_COL,
    LONG_READ_TOOLS_COL,
    TOTAL_TOOLS_COL,
    ILLUMINA_MAX_READS_COL,
    ISOSEQ_MAX_READS_COL,
    ISOSEQ_MAX_VARIANTS_COL,
    'arriba',
    'fusioncatcher',
    'pizzly',
    'squid',
    'starfusion',
    'jaffa',
    'cupcake',
    KNOWN_ALL_COL,
    OCCURS_IN_HEALTHY_COL,
    NOTES_COL
]

data={FUSION_GENE_COL: fusion_genes}
fusion_genes = pandas.DataFrame(data=data, columns=columns)
fusion_genes[[LEFT_GENE_COL,RIGHT_GENE_COL]] = \
    fusion_genes[FUSION_GENE_COL].str.split(FUSION_GENE_SEPARATOR,expand=True)
fusion_genes.set_index(FUSION_GENE_COL, inplace=True)

for fusion_gene, row in fusion_genes.iterrows():
    manually_confirmed = ""
    passes_filtering = ""
    left_chromosome = ""
    right_chromosome = ""
    short_read_tools = 0
    long_read_tools = 0
    total_tools = 0
    illumina_reads = 0
    isoseq_reads = 0
    isoseq_variants = 0
    known_ALL = ""
    occurs_in_gm12878 = ""
    notes = ""

    def format_notes(notes):
        return NOTES_SEPARATOR + NOTES_SEPARATOR.join(list(notes)).replace('["', '').replace('"]', '')

    # Short-read
    if fusion_gene in illumina_all:
        short_read_tools = 1

    fusions_illumina_filtered = illumina_filtered.loc[illumina_filtered[FUSION_GENE_COL] == fusion_gene]
    if not fusions_illumina_filtered.empty:
        illumina_filtered_row = fusions_illumina_filtered.iloc[0]
        short_read_tools = illumina_filtered_row["Tools hits"]

    # Arriba
    fusions_arriba = arriba.loc[arriba[FUSION_GENE_COL] == fusion_gene]
    if not fusions_arriba.empty:
        first_arriba_fusion = fusions_arriba.iloc[0]
        left_chromosome = first_arriba_fusion[LEFT_CHROM_COL]
        right_chromosome = first_arriba_fusion[RIGHT_CHROM_COL]
        arriba_reads = sum(fusions_arriba['split_reads']) + sum(fusions_arriba['discordant_mates'])
        if arriba_reads > illumina_reads:
            illumina_reads = arriba_reads
        notes = notes + format_notes(fusions_arriba[NOTES_COL])
        row['arriba'] = YES

    # Fusioncatcher
    fusions_fusioncatcher = fusioncatcher.loc[fusioncatcher[FUSION_GENE_COL] == fusion_gene]
    if not fusions_fusioncatcher.empty:
        first_fusioncatcher_fusion = fusions_fusioncatcher.iloc[0]
        left_chromosome = first_fusioncatcher_fusion[LEFT_CHROM_COL]
        right_chromosome = first_fusioncatcher_fusion[RIGHT_CHROM_COL]
        split_reads_fusioncatcher = sum(fusions_fusioncatcher['Spanning_unique_reads'])
        # Multiply by 2 since fusioncatcher counts spanning pairs, not reads
        spanning_reads_fusioncatcher = sum(fusions_fusioncatcher['Spanning_pairs']) * 2
        fusioncatcher_reads = split_reads_fusioncatcher + spanning_reads_fusioncatcher
        if fusioncatcher_reads > illumina_reads:
            illumina_reads = fusioncatcher_reads
        notes = notes + format_notes(fusions_fusioncatcher[NOTES_COL])
        row['fusioncatcher'] = YES

    # Pizzly
    fusions_pizzly = pizzly.loc[pizzly[FUSION_GENE_COL] == fusion_gene]
    if not fusions_pizzly.empty:
        split_reads_pizzly = sum(fusions_pizzly['splitcount'])
        # Multiply by 2 since pizzly counts spanning pairs, not reads
        spanning_reads_pizzly = sum(fusions_pizzly['paircount']) * 2
        pizzly_reads = split_reads_pizzly + spanning_reads_pizzly
        if pizzly_reads > illumina_reads:
            illumina_reads = pizzly_reads
        row['pizzly'] = YES

    # Squid
    fusions_squid = squid.loc[squid[FUSION_GENE_COL] == fusion_gene]
    if not fusions_squid.empty:
        first_squid_fusion = fusions_squid.iloc[0]
        left_chromosome = "chr" + first_squid_fusion["# chrom1"]
        right_chromosome = "chr" + first_squid_fusion["chrom2"]
        # Squid uses just one metric
        squid_reads = sum(fusions_squid["score"])
        if squid_reads > illumina_reads:
            illumina_reads = squid_reads
        notes = notes + format_notes(fusions_squid["Type"])
        row['squid'] = YES

    # Starfusion
    fusions_starfusion = starfusion.loc[starfusion[FUSION_GENE_COL] == fusion_gene]
    if not fusions_starfusion.empty:
        first_starfusion_fusion = fusions_starfusion.iloc[0]
        left_chromosome = first_starfusion_fusion[LEFT_CHROM_COL]
        right_chromosome = first_starfusion_fusion[RIGHT_CHROM_COL]
        split_reads_starfusion = sum(fusions_starfusion['JunctionReadCount'])
        spanning_reads_starfusion = sum(fusions_starfusion['SpanningFragCount'])
        starfusion_reads = split_reads_starfusion + spanning_reads_starfusion
        if starfusion_reads > illumina_reads:
            illumina_reads = starfusion_reads
        notes = notes + format_notes(fusions_starfusion[NOTES_COL])
        row['starfusion'] = YES

    # Long-read
    fusions_jaffa = jaffa.loc[jaffa[FUSION_GENE_COL] == fusion_gene]
    if not fusions_jaffa.empty:
        first_jaffa_fusion = fusions_jaffa.iloc[0]
        left_chromosome = first_jaffa_fusion["chrom1"]
        right_chromosome = first_jaffa_fusion["chrom2"]
        long_read_tools += 1
        variants_jaffa = fusions_jaffa[['base1', 'base2']].value_counts().reset_index(name='count')
        isoseq_variants = len(variants_jaffa)
        reads_jaffa = sum(fusions_jaffa['spanning reads'])
        isoseq_reads = reads_jaffa
        known = list(fusions_jaffa['known'])
        if "Yes" in known:
            notes = notes + format_notes(["known"])
        row['jaffa'] = YES

    fusions_cupcake = cupcake.loc[cupcake[FUSION_GENE_COL] == fusion_gene]
    if not fusions_cupcake.empty:
        first_cupcake_fusion = fusions_cupcake.iloc[0]
        left_chromosome = first_cupcake_fusion[LEFT_CHROM_COL]
        right_chromosome = first_cupcake_fusion[RIGHT_CHROM_COL]
        long_read_tools += 1
        variants_cupcake = fusions_cupcake[['LeftBreakpoint', 'RightBreakpoint']].value_counts().reset_index(name='count')
        if len(variants_cupcake) > isoseq_variants:
            isoseq_variants = len(variants_cupcake)
        reads_cupcake = sum(fusions_cupcake['SpanningReads'])
        if reads_cupcake > isoseq_reads:
            isoseq_reads = reads_cupcake
        row['cupcake'] = YES

    # Get any missing chromosome info from Ensembl, or cache it for future use
    left_gene = row[LEFT_GENE_COL]
    right_gene = row[RIGHT_GENE_COL]
    if left_chromosome == "":
        left_chromosome = get_chromosome(left_gene)
    else:
        set_chromosome(left_gene, left_chromosome)

    if right_chromosome == "":
        right_chromosome = get_chromosome(right_gene)
    else:
        set_chromosome(right_gene, right_chromosome)

    # De-dupe and clean up notes
    def filter_notes(note):
        if len(note) == 0 or note == "." or note == "fusion-gene":
            return False
        return True

    notes = map(lambda str: str.strip('"'), list(set(notes.split(NOTES_SEPARATOR))))
    notes = NOTES_SEPARATOR.join(list(filter(filter_notes, notes)))

    # Is it cancer?
    if(left_gene in whitelist or
        right_gene in whitelist):
        known_ALL = YES

    # Is it a "healthy" fusion?
    if(fusion_gene in gm12878):
        occurs_in_gm12878 = YES

    # Has this fusion been manually confirmed?
    # NOTE: this is not part of the filtering process
    if(fusion_gene in manually_confirmed_list):
        manually_confirmed = YES

    # CONDITIONAL FILTERING
    # Filter genes that are not in GM12878, with different partner genes, and have one or more of the following:
    # - is known cancer gene with at least 5 short or long reads supporting
    # - 1 short-read and 1 long-read tool support
    # - 3 short-read tool support
    # - 10 supporting isoseq reads

    is_known_ALL_5_reads = known_ALL and (isoseq_reads >= 5 or illumina_reads >= 5)
    is_supported_by_short_and_long = short_read_tools >= 1 and long_read_tools >= 1
    is_supported_by_3_short_tools = short_read_tools >= 3
    is_supported_by_10_long_reads = isoseq_reads >= 10
    passes_min_support = is_known_ALL_5_reads or is_supported_by_short_and_long or is_supported_by_3_short_tools or is_supported_by_10_long_reads
    
    if not occurs_in_gm12878 and passes_min_support and left_gene != right_gene:
        passes_filtering = YES


    row[MANUALLY_CONFIRMED_COL] = manually_confirmed
    row[PASSED_AUTOMATED_FILTERING_COL] = passes_filtering
    row[LEFT_CHROM_COL] = left_chromosome
    row[RIGHT_CHROM_COL] = right_chromosome
    row[SHORT_READ_TOOLS_COL] = short_read_tools
    row[LONG_READ_TOOLS_COL] = long_read_tools
    row[TOTAL_TOOLS_COL] = short_read_tools + long_read_tools
    row[ILLUMINA_MAX_READS_COL] = illumina_reads
    row[ISOSEQ_MAX_READS_COL] = isoseq_reads
    row[ISOSEQ_MAX_VARIANTS_COL] = isoseq_variants
    row[KNOWN_ALL_COL] = known_ALL
    row[OCCURS_IN_HEALTHY_COL] = occurs_in_gm12878
    row[NOTES_COL] = notes

# Save
fusion_genes.to_csv(OUT_DIR + "/REH.svs.filtered.csv")
