import re
import pandas
from pysam import VariantFile

DATA_DIR="../data"
ILLUMINA = DATA_DIR + "/illumina.tiddit.vcf"
PB = DATA_DIR + "/pb.sniffles.vcf"
ONT = DATA_DIR + "/ont.sniffles.vcf"

OUT_DIR="out"
OUT_FILE=OUT_DIR + "/REH.svs.filtered.csv"

INVALID_POSITION = -1

# MIN SUPPORT
MIN_SUPPORT_READS = 5
MIN_SUPPORT_PCT = 0.20

# MAX COVERAGE, depends on dataset averages
# Very high coverage is indicative of chaotic regions
MAX_COV_PCT = 0.50
AVG_COV = {
	'ILLUMINA': {
		'ALL': 34.7,
		'16': 51.5,
		'X': 17.6
	},
	'PB': {
		'ALL': 14.6,
		'16': 20.8,
		'X': 7.5,
	},
	'ONT': {
		'ALL': 18.2,
		'16': 26.4,
		'X': 9.2
	}
}

illumina = VariantFile(ILLUMINA)
pb = VariantFile(PB)
ont = VariantFile(ONT)

# PREPARE NEW DATA FRAME
columns = [
	"dataset",
	"type",
	"chrom1",
	"pos1",
	"chrom2",
	"pos2",
	"size",
	"support",
	"coverage"
]
svs = pandas.DataFrame(columns=columns)

# COUNTS
# Raw
n_illumina = 0
n_pb = 0
n_ont = 0

# After size/type filtering
n_illumina_prefiltered = 0
n_pb_prefiltered = 0
n_ont_prefiltered = 0

# Filtered
n_illumina_filtered = 0
n_pb_filtered = 0
n_ont_filtered = 0

# HELPER FUNCTIOMS
def sv_length(rec):
	len = rec.info.get('SVLEN')
	try:
		return abs(len[0])
	except:
		try: 
			return abs(len)
		except:
			return 0

def chromosome(str):
	# Possible formats: chr1, chrX, chrUn_KI270742v1, N[chr14_GL000225v1_random:47528[
	# We exclude chrY since we know REH is from a female.
	m = re.search("chr([0-9X]*)", str)
	if m:
		return m.group(1)

def chromosomes(rec):
	chrom1 = chromosome(rec.chrom)
	chrom2 = chromosome(rec.alts[0]) or chrom1
	return (chrom1, chrom2)

def pos2(rec):
	p = rec.info.get('END')
	if p: 
		return p
	if is_bnd(rec):
		m = re.search("chr[0-9X]*:(\d*)", rec.alts[0])
		if not m:
			return INVALID_POSITION
		return m.group(1)

def passes_quality(rec):
	passes = True
	for filter in rec.filter.values():
		if filter.name != "PASS":
			passes = False
	return passes

def is_bnd(rec):
	is_bnd_type = rec.info.get('SVTYPE') == "BND"
	if not is_bnd_type:
		return False

	chrom1, chrom2 = chromosomes(rec)
	return chrom1 and chrom2 and chrom1 != chrom2

def is_big(rec):
	return sv_length(rec) >= 100000

def append_sv(rec, dict, svs):
	chrom1, chrom2 = chromosomes(rec)
	row = {
		"type": rec.info.get('SVTYPE'),
		"chrom1": chrom1,
		"pos1": rec.pos,
		"chrom2": chrom2,
		"pos2": pos2(rec),
		"size": sv_length(rec),
	}
	merged = {**row, **dict}
	df = pandas.DataFrame([merged])
	return pandas.concat([svs, df], ignore_index=True)

def max_cov_threshold(chr1, chr2, dataset):
	# If one of the chromosomes is 16, use chr16 DOC
	# If both chromosomes are X, use chrX DOC
	# Else use ALL DOC
	key = "ALL"
	chr1 = str(chr1)
	chr2 = str(chr2)
	if chr1 == "16" or chr2 == "16":
		key = "16"
	if chr1 == "X" and chr2 == "X":
		key = "X"

	# Get the coverage for the correct dataset and chr
	avg_cov = AVG_COV[dataset][key]
	return avg_cov + (avg_cov * MAX_COV_PCT)

# FILTERING CRITERIA:
# - all filters pass (FILTER="PASS")
# - SVTYPE="BND" (on different chromosomes) or SV len >= 100000
# - support is specified percentage of the coverage
# - coverage is below max threshold
def passes_prefilters(rec):
	return is_bnd(rec) or is_big(rec)

# TIDDIT FUNCTIONS
def passes_tiddit_filters(rec):
	# Coverage
	cov = get_tiddit_coverage(rec)
	chrom1, chrom2 = chromosomes(rec)
	cov_ok = cov <= max_cov_threshold(chrom1, chrom2, "ILLUMINA") 

	# Support
	support = rec.info.get('LTE')
	avg_cov = AVG_COV["ILLUMINA"]["ALL"]
	support_ok = support >= MIN_SUPPORT_READS and \
		support >= avg_cov * MIN_SUPPORT_PCT

	return passes_quality(rec) and cov_ok and support_ok

def get_tiddit_coverage(rec):
	return (rec.info.get('COVA') + rec.info.get('COVB')) / 2

def append_tiddit_row(rec, svs):
	row = {
		"dataset": "ILLUMINA",
		"support": rec.info.get('LTE'),
		"coverage": get_tiddit_coverage(rec) 
	}
	return append_sv(rec, row, svs)

# SNIFFLES FUNCTIONS
# - Filter out those labeled by Sniffles as "IMPRECISE"
def passes_sniffles_filters(rec, dataset):
	# Quality
	qual_ok = passes_quality(rec) and 'IMPRECISE' not in rec.info.keys()	

	# Coverage
	cov = get_sniffles_coverage(rec)
	chrom1, chrom2 = chromosomes(rec)
	cov_ok = cov <= max_cov_threshold(chrom1, chrom2, dataset)

	# Support
	support = rec.info.get('SUPPORT')
	avg_cov = AVG_COV[dataset]["ALL"]
	support_ok = support >= MIN_SUPPORT_READS and \
		support >= avg_cov * MIN_SUPPORT_PCT

	return qual_ok and cov_ok and support_ok

def get_sniffles_coverage(rec):
	cov_list = rec.info.get('COVERAGE')
	sum = 0
	for c in cov_list:
		if c is None:
			c = 0
		sum += c
	return sum / len(cov_list)

def append_sniffles_row(rec, svs, dataset):
	row = {
		"dataset": dataset,
		"support": rec.info.get('SUPPORT'),
		"coverage": get_sniffles_coverage(rec)
	}
	return append_sv(rec, row, svs)

# FILTER ILLUMINA
for rec in illumina.fetch():
	n_illumina += 1
	if(passes_prefilters(rec)):
		n_illumina_prefiltered += 1
		if (passes_tiddit_filters(rec)):
			n_illumina_filtered += 1
			svs = append_tiddit_row(rec, svs)

# FILTER PB
for rec in pb.fetch():
	n_pb += 1
	if(passes_prefilters(rec)):
		n_pb_prefiltered += 1
		if (passes_sniffles_filters(rec, "PB")):
			n_pb_filtered += 1
			svs = append_sniffles_row(rec, svs, "PB")


# FILTER ONT
for rec in ont.fetch():
	n_ont += 1
	if(passes_prefilters(rec)):
		n_ont_prefiltered += 1
		if (passes_sniffles_filters(rec, "ONT")):
			n_ont_filtered += 1
			svs = append_sniffles_row(rec, svs, "ONT")


print("{} unfiltered SV candidates found in Illumina data.".format(n_illumina))
print("{} unfiltered SV candidates found in PB data.".format(n_pb))
print("{} unfiltered SV candidates found in ONT data.".format(n_ont))

print("{} SVs filtered by type/size found in Illumina data.".format(n_illumina_prefiltered))
print("{} SVs filtered by type/size found in PB data.".format(n_pb_prefiltered))
print("{} SVs filtered by type/size found in ONT data.".format(n_ont_prefiltered))

print("{} SVs filtered by quality/support found in Illumina data.".format(n_illumina_filtered))
print("{} SVs filtered by quality/support found in PB data.".format(n_pb_filtered))
print("{} SVs filtered by quality/support found in ONT data.".format(n_ont_filtered))

# Save
svs.index = [x for x in range(1, len(svs.values)+1)]
svs.index.name = 'id'
svs.to_csv(OUT_FILE)