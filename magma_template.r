## -----------------------------------------------------
pvals <- get.drivers(prefix, sets, name2entrez=name2entrez, grep.cmd=grep.cmd, redo=FALSE)
pvals <- filter(pvals, SYMBOL %in% nodes(g))

ggplot(aes_string(x=set, y="ZSTAT"), data=pvals) + geom_boxplot() + geom_point(aes_string(x=set, y="ZSTAT"), data=gene_pvals[pvals[,set],]) + labs(title=trait)


## -----------------------------------------------------
datatable(dplyr::filter(pvals, !!(as.symbol(set))) %>% dplyr::select(SYMBOL, ZSTAT, P, variant_id, gwas_pvalue=minp))


## -----------------------------------------------------
summary(dplyr::filter(pvals, !!(as.symbol(set)))[["P"]])


## -----------------------------------------------------
table(dplyr::filter(pvals, !!(as.symbol(set)))[["FDR"]] < 0.05)


## -----------------------------------------------------
summary(dplyr::filter(pvals, !!(as.symbol(set)))[["minp"]])


## -----------------------------------------------------
table(dplyr::filter(pvals, !!(as.symbol(set)))[["minp"]] < 5e-8)


## ---- fig.width=20, fig.height=20---------------------
plot.graph.with.gwas(g, pvals, virus_proteins, label.only.gwas = TRUE)


## ---- fig.height=0.5, fig.width=3, echo=FALSE---------
plot.color.legend(pvals)

