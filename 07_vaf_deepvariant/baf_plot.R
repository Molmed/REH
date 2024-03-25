library(ggplot2)
library(RColorBrewer)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(ggrastr)

#save(snakemake, file=paste0(snakemake@rule, '_', snakemake@wildcards$sample, '.Rd'))
#load('plot_baf_logr_Illumina.Rd')

load_fun <- function(x) {
	dd <- fread(x)
	dd[['method']] <- sapply(strsplit(basename(x), '\\.'), '[', 1)
	return(dd)
}

tmp <- rbindlist(lapply(snakemake@input, load_fun))

chroms <- unique(tmp[[1]])
use_index <- T
if (use_index) {
	foo <- tmp[, .N, by='V1']
	sl <- foo[['N']]
	names(sl) <- foo[['V1']]
} else {
	sl <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[chroms]
}

corr <- c(0, cumsum(as.numeric(sl))[-length(sl)])
names(corr) <- names(sl)

if (use_index) {
	tmp[, genome_pos:=1:nrow(tmp)]
} else {
	tmp[, genome_pos:=V2 + corr[V1]]
}

tmp[[1]] <- factor(tmp[[1]], levels=chroms)
tmp[['odd']] <- factor(as.numeric(tmp[['V1']]) %% 2)


med_depth <- tmp[,.(median_depth=median(V3)), by='method']
tmp[['logR']] <- log2(tmp[['V3']]) - log2(med_depth[tmp[, method], on='method'][, median_depth])

plot_fun_baf <- function(x, xl, brk) {
	x <- x[V6 < 1]
	p <- plot_fun('V6', x, xl, c(0, 1), brk, hl=0.5)
	p <- p + ylab('Variant allele frequency')
	return(p)
}

plot_fun_lr <- function(x, xl, brk) {
	p <- plot_fun('logR', x, xl, c(-3,3), brk, hl=0)
	p <- p + ylab('LogR')
	return(p)
}

plot_fun <- function(v, x, xl, yl, brk, hl) {
	lab_pos <- sapply(seq_along(brk)[-1], function(i) (brk[i] - brk[i-1])/2 + brk[i-1])
	p <- ggplot(x, aes_string(x='genome_pos', y=v, color='odd')) + 
				geom_rug(aes(x=x, y=NULL, color=NULL), data=data.frame(x=brk), 
									sides='b', outside=T) +
				geom_hline(yintercept=hl, color='grey70') + 
				rasterize(geom_point(alpha=0.05, shape=19, size=0.1), dpi=300) + 
				scale_color_manual(values=c('red', 'blue')) + 
				guides(color='none') + 
				scale_x_continuous(breaks=lab_pos, labels=levels(x[['V1']]), 
								expand=expansion(mult=0.01)) + 
				scale_y_continuous(limit=yl) + 
				xlab('') + 
				coord_cartesian(clip='off', xlim=xl) + 
				theme(axis.text.x=element_text(angle=45, vjust=, hjust=1),
					axis.ticks.x=element_blank(),
					panel.grid=element_blank(),
					panel.background=element_rect(fill=NA, color='grey50'))
}

plot_panel <- function(x, xl, brk) {
	p_baf <- plot_fun_baf(x, xl, brk)
	p_lr <-  plot_fun_lr(x, xl, brk)
	return(wrap_plots(list(p_baf, p_lr), ncol=1))
}


spl_method <- split(tmp, tmp$method)
p_list <- list()
for (n in names(spl_method)) {
	dd_met <- spl_method[[n]]
	p_list[[n]] <- plot_panel(dd_met, range(dd_met[['genome_pos']]), c(0, cumsum(as.numeric(sl))))
}
pdf(snakemake@output[[1]],  width=16, height=8)
bar <- lapply(p_list, print)
dev.off()
