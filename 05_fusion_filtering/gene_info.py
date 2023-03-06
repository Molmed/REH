import requests
import pickle

GENE_API="https://rest.ensembl.org/lookup/symbol/homo_sapiens/{}?content-type=application/json"
GENE_DB="data/ensembl_gene_info.pickle"
try:
    with open(GENE_DB, 'rb') as f:
        genes_db = pickle.load(f)
except:
    genes_db = {}

def get_chromosome(symbol):
    if symbol in genes_db:
        return genes_db[symbol]

    api_url = GENE_API.format(symbol)
    gene = requests.get(api_url).json()
    try:
        chrom = "chr" + gene["seq_region_name"]
    except:
        chrom = ""
    genes_db[symbol] = chrom
    print("Fetched and returned {} from Ensembl gene {}".format(chrom, symbol))
    return chrom

def set_chromosome(symbol, chromosome):
    if symbol in genes_db and genes_db[symbol] == chromosome:
        return
    genes_db[symbol] = chromosome
    save_db()

def save_db():
    with open(GENE_DB, 'wb') as f:
        pickle.dump(genes_db, f)