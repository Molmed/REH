library(VariantAnnotation)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)

#save(snakemake, file=paste0(snakemake@rule, '_', snakemake@wildcards$sample, '_', snakemake@wildcards$chr, '.Rd'))
#load('allele_counter_vcf_Illumina_chr22.Rd')

dd <- fread(snakemake@input$loci)

rng <- GRanges(seqnames=dd[[1]], ranges=IRanges(start=dd[[2]], width=1))

tf <- open(TabixFile(snakemake@input$vcf, yieldSize=1000000))

vcf_tmp <- readVcf(tf, genome='BSgenome.Hsapiens.UCSC.hg38')

vcf_list <- list()
chunk <- 1
while(nrow(vcf_tmp)) {
	hits <- findOverlaps(vcf_tmp, rng)
	vcf_list[[chunk]] <- vcf_tmp[queryHits(hits)]
	chunk <- chunk + 1
	print(paste0('reading chunk ', chunk))
	vcf_tmp <- readVcf(tf, genome='BSgenome.Hsapiens.UCSC.hg38')
}
vcf <- do.call('rbind', vcf_list)
rm(vcf_list)
gc()

dd <- geno(vcf)$DP[,1]

# ONT calling only reports alt count in AD field, 
# the others ref and alt
if (geno(header(vcf))['AD','Number'] == 'R' ) {
	rr <- sapply(geno(vcf)$AD[,1], '[', 1)
	aa <- sapply(geno(vcf)$AD[,1], function(x) max(x[2:length(x)]))
} else {
	aa <- sapply(geno(vcf)$AD[,1], function(x) max(x))
	rr <- dd - aa
}	
	
dd <- data.table(chr=as.character(seqnames(vcf)), 
		pos=start(vcf),
		dp=dd,
		ref=rr,
		alt=aa,
		baf=aa/dd)

dd <- dd[!is.na(dp)]
		
fwrite(dd, file=snakemake@output[[1]], col.names=F, row.names=F, quote=F, sep='\t')
