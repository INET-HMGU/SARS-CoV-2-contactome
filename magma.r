## ----setup, include=FALSE, message=FALSE, warning=FALSE----
.libPaths(c("~/packages/R/x86_64-redhat-linux-gnu-library/3.6/", .libPaths()))
knitr::opts_chunk$set(echo = TRUE, fig.width=8, fig.height=7)
knitr::opts_knit$set(root.dir=normalizePath(".."))
knitr::opts_knit$set(tidy=TRUE)

## knitr::opts_chunk$set(dev='pdf')
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(DT) # interactive html tables
library(knitr)
library(kableExtra)
library(Rgraphviz) # plot networks
library(Homo.sapiens) # annotation data
library(scales)
library(RColorBrewer)
library(ggpubr)

theme_set(theme_bw())


## -----------------------------------------------------
PPI <- read.csv("data/current/PPI/binary_edge_1126.csv", stringsAsFactors = FALSE)
virus_proteins <- unique(PPI$virus)
human_proteins <- unique(PPI$human)
all_proteins <- union(virus_proteins, human_proteins)


## -----------------------------------------------------
g <- new("graphNEL")
g <- graph::addNode(all_proteins, g)
g <- addEdge(PPI$virus, PPI$human, g)

g1 <- layoutGraph(g, layoutType="neato")
nodeRenderInfo(g1) <- list(shape="ellipse", height=15, widht=30)
renderGraph(g1)


## ---- echo=FALSE--------------------------------------
## define a wrapper function to call magma
magma.run <- function(...) { #function
    args <- list(...)
    cmd <- "./packages/magma/magma"
    for (arg in names(args)) {
      val <- args[[arg]]
      if (!is.null(val)) {
        cmd = paste0(cmd, " --", arg, " ", val)
      }
    }
    cat(cmd, "\n")
    system(cmd, ignore.stdout = TRUE)
}

## define a wrapper for preprocessing
preprocess.for.magma <- function(input.file, prefix, gene.loc.file, p.column, snp.column, chr.column, pos.column, sample.size) {
  ## files are very big, so we will do the parsing and reformating with awk and sed
  outdir <- dirname(prefix)
  dir.create(outdir, recursive=TRUE, showWarnings = FALSE)
  
  cat.cmd <- "cat"
  sed.cmd <- "sed" ## gsed needed on mac os
  if (length(grep(".gz$", input.file)) > 0) {
    cat.cmd <- "zcat" ## gzcat needed on mac os will fix later
  }
  
  system("pwd")
  cat(input.file, "\n")
  
  ## Rename the header column to P and SNP (from p and snptestid)
  cmd <- paste0(cat.cmd, " ", input.file, " | head -n 1 ")
  if (snp.column != "SNP") { 
    cmd <- paste0(cmd, " | ", 
                  sed.cmd, " 's/\\bSNP\\b/__old__SNP/g' | ", ## rename old SNP column (if exists)
                  sed.cmd, " 's/\\b", snp.column, "\\b/SNP/g' ") ## rename new SNP column 
  }
  if (p.column != "P") {
    cmd <- paste0(cmd, " |",
                  sed.cmd, " 's/\\bP\\s+/__old__P/g' | ", ## rename old Pvalue column (if exists)
                  sed.cmd, " 's/\\b", p.column, "\\b/P/g'") ## rename P val col
  }
  # cmd <- paste0(cmd, " > ", prefix, "_P.txt")
  # cat(cmd, "\n")
  # system(cmd)
  # cmd <- paste0(cat.cmd, " ", input.file, "| tail -n +2 >> ", prefix, "_P.txt")
  # cat(cmd, "\n")
  # system(cmd)
  
  ## Annotate the SNPs of our GWAS to genes. SNP locations are formated:
  ## rsid, chrom, bp
  ## find the indices of the rsid, chrom and bp
  cn <- colnames(read.csv(input.file, sep="\t", nrow=2))
  col.idx <- match(c(snp.column, chr.column, pos.column), cn)
  if (any(is.na(col.idx))) {
    cat("colname(s)", c(snp.column, chr.column, pos.column)[is.na(col.idx)], "not found in", cn, "\n")
    stop()
  }
  col.idx <- paste0("$", col.idx)
  cmd <- paste0(cat.cmd, " ", input.file, " | tail -n +2 ", " | awk 'BEGIN{OFS=\"\t\"}{gsub(/chr/, \"\", ", col.idx[2], "); print ", 
                paste(col.idx, collapse=", "), "}' | sed 's/$chr//g' > ", prefix, "_snppos.txt")
  cat(cmd, "\n")
  system(cmd)
  
  ## annotate SNPs to genes
  magma.run(annotate="", `snp-loc`=paste0(prefix, "_snppos.txt"), `gene-loc`=gene.loc.file, out=prefix)
  
  ## unzip the file if needed
  if (length(grep(".gz$", input.file)) > 0) {
    new.input <- paste0(prefix, "_P.txt")
    cmd <- paste(cat.cmd, input.file, ">", new.input)
    system(cmd)
    input.file <- new.input
  }
  
  ## run the gene level analysis
  if (is.na(as.numeric(sample.size))) {
    sample.size <- paste0("ncol=", sample.size)
  } else {
    sample.size <- paste0("N=", sample.size)
  }
  pval.arg <- paste0(input.file, " ", sample.size, " use=", paste(snp.column, p.column, sep=","))
  magma.run(bfile="data/current/references/magma/g1000_eur",  pval=pval.arg, `gene-annot`=paste0(prefix, ".genes.annot"), out=prefix)
}


## ---- message=FALSE-----------------------------------
symbol2entrez <- select(Homo.sapiens, columns=c("SYMBOL","ENTREZID"), keys=keys(Homo.sapiens, keytype="SYMBOL"), keytype="SYMBOL")
alias2entrez <- select(Homo.sapiens, columns=c("ALIAS","ENTREZID"), keys=keys(Homo.sapiens, keytype="ALIAS"), keytype="ALIAS")
colnames(alias2entrez)[1] <- "SYMBOL"
name2entrez <- unique(rbind(symbol2entrez, alias2entrez))

write.table(name2entrez, file="results/current/name2entrez.txt", sep="\t", quote=F, row.names=F)

table(human_proteins %in% name2entrez$SYMBOL)


## -----------------------------------------------------
PPI <- merge(PPI, name2entrez, by.x="human", by.y="SYMBOL")

sets <-c(all=list(unique(PPI$ENTREZID)), tapply(PPI$ENTREZID, PPI$virus, as.list))

setfile <- "results/current/gene_sets.txt"
for (set in names(sets)) {
  cat(set, "\t", paste(sets[[set]], collapse="\t"), "\n", sep="", file=setfile, append=(set != names(sets)[1]))
}


## ---- echo=FALSE--------------------------------------
run.magma.on.gwas.list <- function(magma.params, set.file, out.suffix, rerun=TRUE) {
  magma.res <- NULL
  for (i in 1:nrow(magma.params)) {
    ## print(magma.params[i,])
    prefix <- magma.params[i,"prefix"]
    gsa.file <- paste0(prefix, out.suffix, ".gsa.out")
    if (rerun || !file.exists(gsa.file)) {
      with(magma.params[i,], {
        print(genome.build)
        if (genome.build == "hg19") {
          gene.loc.file <- "data/current/references/magma/NCBI37.3.gene.loc"
        } else if (genome.build == "hg38") {
          gene.loc.file <- "data/current/references/magma/NCBI38.gene.loc"
        } else {
          cat("Only genome build hg19 and hg38 currently available!!\n")
          next  
        }
        gene.level.file <- paste0(prefix, ".genes.raw")
        if (!file.exists(gene.level.file)) {
          preprocess.for.magma(input.file, prefix, gene.loc.file, p.column, snp.column, chrom.column, pos.column, sample.size)
        }
        magma.run(`gene-results`=gene.level.file, `set-annot`=set.file, out=paste0(prefix, out.suffix)) # model="condition-hide=Average direction=greater",
        })
    }
    enrichment <- read.table(gsa.file, stringsAsFactors=F, header=TRUE, comment="#")
    enrichment <- data.frame(enrichment, 
                      prefix,
                      trait=basename(prefix),
                      FDR=p.adjust(enrichment$P, "BH"),
                      stringsAsFactors=FALSE)
    magma.res <- rbind(magma.res, enrichment)
  }
  write.table(magma.res, file=paste0("magma_enrichment", out.suffix, ".txt"), sep="\t", quote=F)
  return(magma.res)
}


## ---- message=FALSE, echo=FALSE-----------------------
gwas.info <- read_tsv("data/current/gtex_gwas_data/gwas_metadata.txt")

gwas.magma.params <- dplyr::rename(gwas.info[,c("new_abbreviation", "Sample_Size", "Tag")], prefix=new_abbreviation, sample.size=Sample_Size)
gwas.magma.params$input.file <- paste0("data/current/gtex_gwas_data/imputed_gwas_hg38_1.1/imputed_", gwas.magma.params$Tag, ".txt.gz")
gwas.magma.params$prefix <- paste0("results/current/magma/", gwas.magma.params$prefix)
## add genome build,  "snp.column"   "p.column"     "chrom.column" "pos.column"
gwas.magma.params <- data.frame(gwas.magma.params, genome.build="hg38", snp.column="variant_id", p.column="pvalue", chrom.column="chromosome", pos.column="position", stringsAsFactors = FALSE)

datatable(gwas.magma.params)                             


## -----------------------------------------------------
magma.res <- run.magma.on.gwas.list(gwas.magma.params, set.file=setfile, "_gtex_gwas_by_viral_protein", rerun=FALSE)
write.table(magma.res, "results/current/magma_gtex_gwas_by_viral_protein.txt", sep="\t", quote=FALSE, row.names = FALSE)


## ---- fig.width=5, fig.height=13----------------------
## ggplot(aes(x=celltype, y=region), data=magma.res) + geom_tile(aes(fill=-log10(P))) + facet_grid(trait~.)
global.cl.plot <- ggplot(aes(x=VARIABLE, y=trait), data=magma.res) + geom_tile(aes(fill=-log10(P))) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_gradient(low="white", high="purple") + geom_point(aes(x=VARIABLE, y=trait), data=magma.res[magma.res$FDR < 0.1,])
print(global.cl.plot)

pdf(file="results/current/magma_gtex_gwas_by_viral_protein.pdf", width=5, height=13)
print(global.cl.plot)
dev.off()


## ---- fig.width=5, fig.height=13----------------------
## ggplot(aes(x=celltype, y=region), data=magma.res) + geom_tile(aes(fill=-log10(P))) + facet_grid(trait~.)
global.cl.plot2 <- ggplot(aes(x=VARIABLE, y=trait), data=magma.res) + geom_tile(aes(fill=BETA/SE)) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_gradient2(low="blue", mid="white", high="red") + geom_point(aes(x=VARIABLE, y=trait), data=magma.res[magma.res$FDR < 0.1,])
print(global.cl.plot2)

pdf(file="results/current/magma_gtex_gwas_by_viral_protein_effect_size.pdf", width=5, height=13)
print(global.cl.plot2)
dev.off()


## ---- fig.width=5, fig.height=3-----------------------
sig.traits <- unique(magma.res$trait[magma.res$FDR < 0.1])
global.cl.plot3 <- ggplot(aes(x=VARIABLE, y=trait), data=filter(magma.res, trait %in% sig.traits)) + geom_tile(aes(fill=-log10(P))) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_gradient(low="white", high="purple") + geom_point(aes(x=VARIABLE, y=trait), data=filter(magma.res, FDR < 0.1 & trait %in% sig.traits))
print(global.cl.plot3)

pdf(file="results/current/magma_gtex_gwas_by_viral_protein_sig_traits_only.pdf", width=5, height=3)
print(global.cl.plot3)
dev.off()


## ---- fig.width=5, fig.height=3-----------------------
sig.traits <- unique(magma.res$trait[magma.res$FDR < 0.1])
global.cl.plot3 <- ggplot(aes(x=VARIABLE, y=trait), data=filter(magma.res, trait %in% sig.traits)) + geom_tile(aes(fill=BETA/SE)) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_gradient2(low="blue", mid="white", high="red") + geom_point(aes(x=VARIABLE, y=trait), data=filter(magma.res, FDR < 0.1 & trait %in% sig.traits))
print(global.cl.plot3)

pdf(file="results/current/magma_gtex_gwas_by_viral_protein_effect_size_sig_traits_only.pdf", width=5, height=3)
print(global.cl.plot3)
dev.off()


## ---- echo=FALSE--------------------------------------
cols <- c("Tag", "PUBMED_Paper_Link", "Phenotype", "Sample_Size", "Population", "Date", "Declared_Effect_Allele", "Genome_Reference", "new_abbreviation", "new_Phenotype", "Category")
magma.res <- merge(magma.res, gwas.info[,cols], by.x="trait", by.y="new_abbreviation")
magma.res <- magma.res[order(magma.res$P),]
write.table(magma.res, "results/current/magma_gtex_gwas_by_viral_protein_with_gwas_info.txt", sep="\t", quote=FALSE, row.names = FALSE)

datatable(magma.res %>% arrange(P))


## ---- echo=FALSE--------------------------------------
grep.cmd <- "grep"
## on my mac need to use gnu grep (faster)
if (system("hostname", intern=TRUE) == "MB080512") {
  grep.cmd <- "/opt/homebrew/bin/ggrep"
}


## ---- echo=FALSE--------------------------------------

get.gene.pvalues <- function(prefix, sets, name2entrez=NULL) {
  af.zscores <- read.table(paste0(prefix, ".genes.out"), header=T, stringsAsFactors=FALSE)
  af.zscores <- dplyr::rename(af.zscores, Nsamples=N)
  expr.cl <- sapply(sets, function(set) af.zscores$GENE %in% set)
  merged <- cbind(af.zscores, expr.cl)
  
  merged <- cbind(merged, leverage=NA, FDR=p.adjust(merged$P, "BH"))
  if (!is.null(name2entrez)) {
    merged <- merge(merged, name2entrez, by.x="GENE", by.y="ENTREZID")
  }
  return(merged)
}

## To get back to the SNP level we use the annotation files of magma. We select all genes in the interaction network. In a first step we reduce the size of the GWAS data.
get.snp.pvalues <- function(prefix, selected.entrez, grep.cmd="grep", redo=FALSE) {
  efile <- paste0(prefix, "_entrez_with_ppi.txt")
  sfile <- paste0(prefix, "_snps_with_ppi.txt")
  pfile <- paste0(prefix, "_P_with_ppi.txt")
  if (!file.exists(pfile) || redo) {
    cat(selected.entrez, file=efile, sep="\n")
    cmd <- paste0(grep.cmd, " -F -w -f ", efile, " ", prefix, ".genes.annot | cut -d '\t' -f 2- | tr '\\t' '\\n' | grep -v NA | sort -u > ", sfile)
    system(cmd)
    cmd <- paste0(grep.cmd, " -F -w -f ", sfile, " ", prefix, "_P.txt > ", pfile)
    system(cmd)
    snp_pval <- read.table(pfile, sep="\t", stringsAsFactors = FALSE)
    colnames(snp_pval) <- colnames(read.csv(paste0(prefix, "_P.txt"), sep="\t", nrows=3))
    
    #In the next step we read in the mapping of genes to SNPs
    ann <- readLines(paste0(prefix, ".genes.annot"))
    ann <- bind_rows(lapply(strsplit(ann, "\t"), function(x) {
      if (length(x) > 1) {
        snp <- x[-1]
      } else {
        snp <- NA
      }
      data.frame(gene=x[1], snp, stringsAsFactors = FALSE)
    }))
    ann <- ann[!is.na(ann$snp),]
    ann <- ann[ann$snp != "NA",]
    ann <- ann[ann$snp != ".",]
    
    # Annotate SNPs with genes and gene symbols
    snp_pval <- merge(snp_pval, ann, by.x="variant_id", by.y="snp")
    write.table(snp_pval, file=pfile, sep="\t", quote=F, row.names=F)
  } else {
    snp_pval <- read.csv(pfile, sep="\t", stringsAsFactors = FALSE)
  }
  return(snp_pval)
}

get.drivers <- function(prefix, sets, name2entrez=NULL, grep.cmd="grep", redo=FALSE) {
  selected.entrez <- unique(unlist(sets))
  gene_pvals <- get.gene.pvalues(prefix, sets, name2entrez)
  snp_pval <- get.snp.pvalues(prefix, selected.entrez, grep.cmd=grep.cmd)
  minp <- group_by(snp_pval, gene) %>% summarise(variant_id=variant_id[which.min(pvalue)], minp=min(pvalue))
  pvals <- merge(gene_pvals, minp, by.x="GENE", by.y="gene")
  return(pvals)
}


## -----------------------------------------------------
plist <- list()
sig.magma <- filter(magma.res, FDR < 0.1)
for (i in 1:nrow(sig.magma)) {
  set <- sig.magma[i, "VARIABLE"]
  prefix <- sig.magma[i, "prefix"]
  trait <- sig.magma[i, "trait"]
  gene_pvals <- get.gene.pvalues(prefix, sets)
  plist <- c(plist, list(ggplot(data=gene_pvals) + geom_boxplot(aes_string(x=set, y="ZSTAT")) + geom_point(aes_string(x=set, y="ZSTAT"), data=gene_pvals[gene_pvals[,set],]) + labs(title=trait)))
}
ggarrange(plotlist = plist, ncol=5, nrow=4)


## -----------------------------------------------------
gwas.hits <- data.frame()
for (prefix in unique(filter(magma.res, FDR < 0.1)[["prefix"]])) {
  snp_pval <- get.snp.pvalues(prefix, selected.entrez, grep.cmd=grep.cmd)
  snp_pval <- merge(snp_pval, name2entrez, by.x="gene", by.y="ENTREZID")
  minp <- group_by(snp_pval, SYMBOL) %>% summarise(gene=gene[which.min(pvalue)], variant_id=variant_id[which.min(pvalue)], gwas_pvalue=min(pvalue)) %>% filter(gwas_pvalue < 5.0e-8)
  if (nrow(minp) > 0) {
    gwas.hits <- rbind(gwas.hits, data.frame(prefix, minp))
  }
}
write.table(gwas.hits, file="results/current/magma_gtex_gwas_snps.txt")


## -----------------------------------------------------
hits <- merge(gwas.hits, PPI, by.x=c("SYMBOL", "gene"), by.y=c("human", "ENTREZID"))
hits <- merge(hits, magma.res, by.x=c("prefix", "virus"), by.y=c("prefix", "VARIABLE"))

filtered.hits <- filter(hits, FDR < 0.1) 
filtered.hits %>% dplyr::select(trait, virus, human=SYMBOL, variant_id, gwas_pvalue, P, FDR) %>% distinct() %>% arrange(human, FDR)  %>% datatable()


## -----------------------------------------------------
datatable(hits %>% dplyr::select(trait, virus, human=SYMBOL, variant_id, gwas_pvalue, P, FDR) %>% distinct() %>% arrange(trait, FDR))


## ---- echo=FALSE--------------------------------------
get_gwas_region <- function(prefix, chrom, start, end, awk.cmd="awk", cat.cmd="cat", redo=FALSE, chr.column="chromosome", pos.column="position") {
  out.file <- paste0(prefix, "_P_", chrom, "_", start, "_", end, ".txt")
  if (!file.exists(out.file) || redo) {
    input.file <- paste0(prefix, "_P.txt")
    ## find the indices of the chrom and bp
    cn <- colnames(read.csv(input.file, sep="\t", nrow=2))
    col.idx <- match(c(chr.column, pos.column), cn)
    if (any(is.na(col.idx))) {
      cat("colname(s)", c(chr.column, pos.column)[is.na(col.idx)], "not found in", cn, "\n")
      stop()
    }
    col.idx <- paste0("$", col.idx)
    ## build an awk script to filter the file
    cmd <- paste0(cat.cmd, " ", input.file, " | tail -n +2 ", " | ", awk.cmd, 
                  " 'BEGIN{OFS=\"\t\"}{gsub(/chr/, \"\", ", col.idx[1], "); ", ## replace chr prefix of chrom names
                  "if (", col.idx[1], ' == "', chrom, '" && ', ## match chrom name
                          col.idx[2], " > ", start, " && ", col.idx[2], " < ", end, ")", ## match position
                  "{print $0;}}' >> ", out.file)
    cat(cmd, "\n")
    ## write header to out file
    system(paste0("head -n 1 ", input.file, " > ", out.file))
    ## then extract the region
    system(cmd)
  }
  
  gwas_pvals <- read_tsv(out.file)
  return(gwas_pvals)
}

## convenience function that works with the gwas input parameter table
get_gwas_region_for_study <- function(prefix, magma.params, chrom, start, end, ...) {
  idx <- which(magma.params$prefix == prefix)
  get_gwas_region(prefix, chrom, start, end, chr.column=magma.params[idx,"chrom.column"], pos.column=magma.params[idx,"pos.column"])
}


## ---- message=FALSE-----------------------------------
for (i in 1:nrow(filtered.hits)) {
  prefix <- as.character(filtered.hits[i,"prefix"])
  gene_info <- get.drivers(prefix, sets)
  gene <- filtered.hits[i,"gene"]
  symbol <- filtered.hits[i,"SYMBOL"]
  chrom <- gene_info[gene_info$GENE == gene, "CHR"]
  start <- gene_info[gene_info$GENE == gene, "START"]
  end <- gene_info[gene_info$GENE == gene, "STOP"]
  gwas_pvals <- get_gwas_region_for_study(prefix, gwas.magma.params, chrom=chrom, start=start - 1e6, end=end + 1e6)
  p <- ggplot(aes(x=position, y=-log10(pvalue)), data=gwas_pvals) + geom_point() + 
    geom_rect(aes(xmin=start, xmax=end, ymin=-0.5, ymax=-0.4), col="red") + 
    # geom_label(aes(x=start + (end - start) / 2, y=-1), label=symbol) + 
    labs(title=paste(basename(prefix), "at", symbol, "Entrez id:", gene))
  print(p)
}


## ---- eval=FALSE--------------------------------------
## ## also add gene annotations -- does not work due to version conflicts
## library(ggbio)
## lim.gr <- GRanges(seqnames=chrom, ranges=IRanges(start, end))
## 
## p.txdb <- autoplot(Homo.sapiens, which=lim.gr)
## gr.gwas <- GRanges(seqnames=chrom, ranges=IRanges(test$position, width=1))
## mCol(gr.gwas) <- DataFrame(pvalue=test$pvalue)
## p.gwas <- autoplot(gr.gwas, aes(y=-log10(pvalue)))
## p <- p.txdb + p.gwas
## print(p)


## -----------------------------------------------------
hits.per.gene <- hits %>% group_by(SYMBOL) %>% summarise(ntraits=n_distinct(trait), nviralproteins=n_distinct(virus)) %>% arrange(-ntraits)
genes.w.multiple.hits <- filter(hits.per.gene, ntraits > 1)[["SYMBOL"]]
hits.per.gene 


## -----------------------------------------------------
filter(hits, SYMBOL %in% genes.w.multiple.hits)  %>% dplyr::select(trait, virus, human=SYMBOL, variant_id, gwas_pvalue, P, FDR) %>% distinct() %>% arrange(human, FDR)  %>% datatable()


## ---- echo=FALSE--------------------------------------
plot.graph.with.gwas <- function(g, pvals, virus_proteins, label.only.gwas=FALSE) {
  g1 <- g
  ## filter for genes in the graph and map p-values to a nice color scale
  pvals <- filter(pvals, SYMBOL %in% nodes(g1))
  breaks <- pretty(-log10(pvals$P), 9)
  pvals$bin <- findInterval(-log10(pvals$P), breaks)
  
  ## assign each bin a color
  pal <- brewer.pal(length(breaks) - 1, "Blues")
  pvals$fill <- pal[pvals$bin]
  
  fill <- pvals$fill
  names(fill) <- pvals$SYMBOL
  
  ## mark the viral proteins 
  red <- brewer.pal(3, "Set1")[1]
  viral.fill <- rep(red, length(virus_proteins))
  names(viral.fill) <- virus_proteins
  fill <- c(fill, viral.fill)
  
  ## shapes default: all ellipse
  ## if label only gwas: viral and gwas are elipse, rest is circle
  height <- 15
  width <- rep(height, numNodes(g1))
  names(width) <- nodes(g1)
  if (label.only.gwas) {
    
    ## labels
    label <- rep("", numNodes(g1))
    names(label) <- nodes(g1)
    label[virus_proteins] <- virus_proteins
    gwas_genes <- pvals[pvals$minp < 5e-8 | pvals$FDR < 0.05,"SYMBOL"]
    label[gwas_genes] <- gwas_genes
    
    ## shapes
    width[c(virus_proteins, gwas_genes)] <- 30
    shape <- rep("circle", numNodes(g1))
    names(shape) <- nodes(g1)
    shape[c(virus_proteins, gwas_genes)] <- "ellipse"
    
  } else {
    shape <- "ellipse"
    width <- 30
    label <- nodes(g1)
    names(label) <- label
  }

  graph.par(list(graph=list(overlap="scalexy", outputorder="edgesfirst"))) #overlap="false"
  g1 <- layoutGraph(g1, layoutType="neato")
  nodeRenderInfo(g1) <- list(label=label, fill=fill, shape=shape, height=height, widht=width)
  
  renderGraph(g1)
}

plot.color.legend <- function(pvals) {
  pvals <- filter(pvals, SYMBOL %in% nodes(g1))
  breaks <- pretty(-log10(pvals$P), 9)
  legend <- data.frame(bottom=breaks[-length(breaks)], top=breaks[-1], y=0)
  legend$mid <- with(legend, bottom + (top - bottom) / 2)
  legend$rev <- rev(legend$top)
  ggplot(data=legend) + geom_tile(aes(x=bottom, y=y, fill=top)) + geom_text(aes(x=bottom, y=0, label=rev)) + theme_void() + theme(legend.position="none")
}


## ---- include=FALSE-----------------------------------
out = NULL
for (i in 1:nrow(filtered.hits)) {
  ## set variables used inside the template child Rmd
  set <- filtered.hits[i, "virus"]
  trait <- filtered.hits[i, "trait"]
  # prefix <- sig.magma[i, "prefix"]
  prefix <- paste0("results/current/magma/", trait)
  
  
  out = c(out, knit_child('magma_template.Rmd'))
}


## ---- message=FALSE, eval=FALSE, echo=FALSE-----------
## ## these are the old communities
## community.huri <- read_tsv("data/current/PPI/HuRI_assays_community_HuRI.txt")
## community.huri.sig <- read_tsv("data/current/PPI/HuRI_signif_community-HuRI.txt")
## community.huri <- filter(community.huri, community %in% community.huri.sig$cluster) %>% inner_join(name2entrez, by=c("node"="SYMBOL")) %>% mutate(type="HuRI", cname=paste0("HuRI_", community))
## 
## community.hi.union <- read_tsv("data/current/PPI/HuRI_assays_community_HI-union.txt")
## community.hi.unionsig <- read_tsv("data/current/PPI/HuRI_signif_community-HuRI-HI-union.txt")
## community.hi.union <- filter(community.hi.union, community %in% community.hi.unionsig$cluster) %>% inner_join(name2entrez, by=c("node"="SYMBOL")) %>% mutate(type="HI-union", cname=paste0("HI-union_", community))
## 
## communities <- rbind(community.huri, community.hi.union)


## ---- message=FALSE-----------------------------------
## these are the new communities
communities <- read_tsv("data/current/PPI/statistics_community_membership.txt")
communities_sig <- read_tsv("data/current/PPI/statistics_community_significant.txt")
communities <- filter(communities, community %in% communities_sig$community) %>% inner_join(name2entrez, by=c("node"="SYMBOL")) %>% mutate(cname=community)


## -----------------------------------------------------
csets <- with(communities, tapply(ENTREZID, cname, as.list))

comm.setfile <- "results/current/community_gene_sets.txt"
for (set in names(csets)) {
  cat(set, "\t", paste(csets[[set]], collapse="\t"), "\n", sep="", file=comm.setfile, append=(set != names(csets)[1]))
}


## -----------------------------------------------------
magma.res.comm <- run.magma.on.gwas.list(gwas.magma.params, set.file=comm.setfile, "_gtex_gwas_by_community", rerun=FALSE)
write.table(magma.res.comm, "results/current/magma_gtex_gwas_by_community.txt", sep="\t", quote=FALSE, row.names = FALSE)


## -----------------------------------------------------
cols <- c("Tag", "PUBMED_Paper_Link", "Phenotype", "Sample_Size", "Population", "Date", "Declared_Effect_Allele", "Genome_Reference", "new_abbreviation", "new_Phenotype", "Category")
magma.res.comm <- merge(magma.res.comm, gwas.info[,cols], by.x="trait", by.y="new_abbreviation")
magma.res.comm <- magma.res.comm[order(magma.res.comm$P),]
write.table(magma.res.comm, "results/current/magma_gtex_gwas_by_community_with_gwas_info.txt", sep="\t", quote=FALSE, row.names = FALSE)


## -----------------------------------------------------
magma.res.comm %>% arrange(FDR) %>% datatable()


## -----------------------------------------------------

## we have many gene sets so only select the ones that are significant at least for one trait
sig.sets <- unique(magma.res.comm$VARIABLE[magma.res.comm$FDR < 0.1])
sig.traits <- unique(magma.res.comm$trait[magma.res.comm$FDR < 0.1])
plot.data <- filter(magma.res.comm, VARIABLE %in% sig.sets & trait %in% sig.traits)
plot.data$VARIABLE <- factor(plot.data$VARIABLE)

## ggplot(aes(x=celltype, y=region), data=magma.res) + geom_tile(aes(fill=-log10(P))) + facet_grid(trait~.)
global.cl.plot <- ggplot(aes(x=VARIABLE, y=trait), data=plot.data) + geom_tile(aes(fill=-log10(P))) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_gradient(low="white", high="purple") + geom_point(aes(x=VARIABLE, y=trait), data=filter(plot.data, FDR < 0.1)) + labs(x="Communitiy")
print(global.cl.plot)

pdf(file="results/current/magma_gtex_gwas_by_community.pdf", width=7, height=7)
print(global.cl.plot)
dev.off()


## -----------------------------------------------------

## we have many gene sets so only select the ones that are significant at least for one trait
sig.sets <- unique(magma.res.comm$VARIABLE[magma.res.comm$FDR < 0.1])
sig.traits <- unique(magma.res.comm$trait[magma.res.comm$FDR < 0.1])
plot.data <- filter(magma.res.comm, VARIABLE %in% sig.sets & trait %in% sig.traits)
plot.data$VARIABLE <- factor(plot.data$VARIABLE)

## ggplot(aes(x=celltype, y=region), data=magma.res) + geom_tile(aes(fill=-log10(P))) + facet_grid(trait~.)
global.cl.plot2 <- ggplot(aes(x=VARIABLE, y=trait), data=plot.data) + geom_tile(aes(fill=BETA/SE)) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_gradient2(low="blue", mid="white", high="red") + geom_point(aes(x=VARIABLE, y=trait), data=filter(plot.data, FDR < 0.1)) + labs(x="Communitiy")
print(global.cl.plot2)

pdf(file="results/current/magma_gtex_gwas_by_community_effect_size.pdf", width=7, height=7)
print(global.cl.plot2)
dev.off()


## ---- fig.width=9, fig.height=14----------------------
plist <- list()
sig.magma <- filter(magma.res.comm, FDR < 0.1)

for (i in 1:nrow(sig.magma)) {
  set <- sig.magma[i, "VARIABLE"]
  cset <- paste0("c", set)
  prefix <- sig.magma[i, "prefix"]
  trait <- sig.magma[i, "trait"]
  gene_pvals <- get.gene.pvalues(prefix, csets)
  ## replace "-" by "_" as ggplot does not like it in variable names
  colnames(gene_pvals)[10:ncol(gene_pvals)] <- paste0("c", colnames(gene_pvals)[10:ncol(gene_pvals)] )
  plist <- c(plist, list(ggplot(data=gene_pvals) + geom_boxplot(aes_string(x=cset, y="ZSTAT")) + geom_point(aes_string(x=cset, y="ZSTAT"), data=gene_pvals[gene_pvals[,cset],]) + labs(title=trait)))
}
ggarrange(plotlist = plist, ncol=6, nrow=8)


## -----------------------------------------------------
selected.entrez.comm <- unique(unlist(csets))
gwas.hits.comm <- data.frame()
for (prefix in unique(filter(magma.res.comm, FDR < 0.1)[["prefix"]])) {
  snp_pval <- get.snp.pvalues(prefix, selected.entrez.comm, grep.cmd=grep.cmd)
  ## snp_pval <- merge(snp_pval, name2entrez, by.x="gene", by.y="ENTREZID")
  minp <- group_by(snp_pval, gene) %>% summarise(variant_id=variant_id[which.min(pvalue)], gwas_pvalue=min(pvalue)) %>% filter(gwas_pvalue < 5.0e-8)
  if (nrow(minp) > 0) {
    gwas.hits.comm <- rbind(gwas.hits.comm, data.frame(prefix, minp))
  }
}
write.table(gwas.hits.comm, file="results/current/magma_gtex_gwas_by_community_snps.txt")


## -----------------------------------------------------
## here we do not merge by direct interaction but by community membership
hits.comm <- merge(gwas.hits.comm, communities, by.x=c("gene"), by.y=c("ENTREZID"))
hits.comm <- merge(hits.comm, magma.res.comm, by.x=c("prefix", "cname"), by.y=c("prefix", "VARIABLE"))

filtered.hits.comm <- filter(hits.comm, FDR < 0.1) 
filtered.hits.comm %>% dplyr::select(trait, cname, human=node, variant_id, gwas_pvalue, P, FDR) %>% distinct() %>% arrange(human, FDR)  %>% datatable()


## -----------------------------------------------------
filtered.hits.comm %>% group_by(node) %>% summarise(ntraits=length(unique(trait)), ncommunities=length(unique(cname))) %>% arrange(-ntraits, -ncommunities)


## -----------------------------------------------------
filtered.hits.comm %>% group_by(trait) %>% summarise(nhits=length(unique(node)), ncommunities=length(unique(cname))) %>% arrange(-nhits, -ncommunities)


## ---- include=FALSE-----------------------------------
out = NULL
for (i in 1:nrow(filtered.hits.comm)) {
  ## set variables used inside the template child Rmd
  set <- filtered.hits.comm[i, "cname"]
  trait <- filtered.hits.comm[i, "trait"]
  # prefix <- sig.magma[i, "prefix"]
  prefix <- paste0("results/current/magma/", trait)
  
  
  out = c(out, knit_child('magma_template_comm.Rmd'))
}

