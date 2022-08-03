## -----------------------------------------------------
pvals <- get.drivers(prefix, csets, name2entrez=name2entrez, grep.cmd=grep.cmd, redo=FALSE)
pvals <- filter(pvals, SYMBOL %in% nodes(g))

## community names are currently just numbers which are not good colnames!
colnames(pvals)[10:(ncol(pvals) - 3)] <- paste0("c", colnames(pvals)[10:(ncol(pvals) - 3)] )
cset <- paste0("c", set)

ggplot(aes_string(x=cset, y="ZSTAT"), data=pvals) + geom_boxplot() + geom_point(aes_string(x=cset, y="ZSTAT"), data=gene_pvals[pvals[,cset],]) + labs(title=trait)


## -----------------------------------------------------
datatable(dplyr::filter(pvals, !!(as.symbol(cset))) %>% dplyr::select(SYMBOL, ZSTAT, P, variant_id, gwas_pvalue=minp))


## -----------------------------------------------------
summary(dplyr::filter(pvals, !!(as.symbol(cset)))[["P"]])


## -----------------------------------------------------
table(dplyr::filter(pvals, !!(as.symbol(cset)))[["FDR"]] < 0.05)


## -----------------------------------------------------
summary(dplyr::filter(pvals, !!(as.symbol(cset)))[["minp"]])


## -----------------------------------------------------
table(dplyr::filter(pvals, !!(as.symbol(cset)))[["minp"]] < 5e-8)

