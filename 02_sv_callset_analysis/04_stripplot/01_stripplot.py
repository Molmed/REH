
import pandas
import matplotlib.pyplot as plt
import seaborn as sns
from pysam import VariantFile
import pickle

# Input files
DATA_DIR = "../../data/sv_callsets"
ILLUMINA = DATA_DIR + "/illumina.tiddit.vcf"
PB = DATA_DIR + "/pb.sniffles.vcf"
ONT = DATA_DIR + "/ont.sniffles.vcf"
CONSENSUS_ALL = "../01_consensus_sets/out/REH_consensus_3x.vcf"
CONSENSUS_LR = "../01_consensus_sets/out/REH_consensus_longread.vcf"

# Output files
OUT_DIR="out"
SIZES_DB=OUT_DIR + "/sizes.pickle"
OUT_FILE=OUT_DIR + "/stripplot.png"

MIN_BP = 100

do_process_sizes = False
try:
    with open(SIZES_DB, 'rb') as f:
        sizes_db = pickle.load(f)
except:
    do_process_sizes = True
    sizes_db = {}

def process_sizes():
	illumina = VariantFile(ILLUMINA)
	pb = VariantFile(PB)
	ont = VariantFile(ONT)
	consensus_all = VariantFile(CONSENSUS_ALL)
	consensus_lr = VariantFile(CONSENSUS_LR)

	# PREPARE NEW DATA FRAME
	columns = [
		"dataset",
		"size"
	]
	sizes = pandas.DataFrame(columns=columns)

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


	def append_size(size, dataset, sizes):
		if (size < MIN_BP):
			return sizes

		row = {
			"dataset": dataset,
			"size": size,
		}
		df = pandas.DataFrame([row])
		return pandas.concat([sizes, df], ignore_index=True)


	for rec in illumina.fetch():
		sizes = append_size(sv_length(rec), "ILLUMINA", sizes)


	for rec in pb.fetch():
		sizes = append_size(sv_length(rec), "PB", sizes)


	for rec in ont.fetch():
		sizes = append_size(sv_length(rec), "ONT", sizes)

	for rec in consensus_lr.fetch():
		sizes = append_size(sv_length(rec), "CONSENSUS (LONG-READ)", sizes)

	for rec in consensus_all.fetch():
		sizes = append_size(sv_length(rec), "CONSENSUS (ALL)", sizes)

	sizes = pandas.DataFrame(sizes.to_dict())
	with open(SIZES_DB, 'wb') as f:
		pickle.dump(sizes, f)

	return sizes

if do_process_sizes:
	sizes_db = process_sizes()

# Size plots
SMALL=1000
MEDIUM=10000
LARGE=1000000

plt.figure(figsize=(16, 11))
plt.subplot(4, 1, 1)
plt.xlim(MIN_BP,SMALL)
a=0.4
ax = sns.stripplot(x='size', y='dataset', data=sizes_db, color='#a73c52',
              alpha=a, size=4)
ax.set(xlabel='size (bp)', ylabel='')
ax.set_title(label="Small SVs (100 bp - 1 kb)", loc='right', fontweight="bold")

plt.subplot(4, 1, 2)
ax = sns.stripplot(x='size', y='dataset', data=sizes_db, color='#6b5f88',
              alpha=a, size=4)
plt.xlim(SMALL,MEDIUM)
ax.set(xlabel='size (bp)', ylabel='')
ax.set_title(label="Medium SVs (1 kbp - 10 kb)", loc='right', fontweight="bold")


plt.subplot(4, 1, 3)
ax = sns.stripplot(x='size', y='dataset', data=sizes_db, color='#3780b3',
              alpha=a, size=4)
plt.xlim(MEDIUM,LARGE)
ax.set(xlabel='size (Mb)', ylabel='')
ax.set_title(label="Large SVs (10 kbp - 1 Mb)", loc='right', fontweight="bold")

plt.subplot(4, 1, 4)
ax = sns.stripplot(x='size', y='dataset', data=sizes_db, color='#47a266',
              alpha=a, size=4)
plt.xlim(MIN_BP,8000000)
ax.set(xlabel='size (Mb)', ylabel='')
ax.set_title(label="All SVs (100 bp - 8 Mb)", loc='right', fontweight="bold")

plt.tight_layout()
plt.savefig(OUT_FILE)
