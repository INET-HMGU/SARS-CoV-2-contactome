# HuRI community establishment with linkcomm package

# 20.04.2022 **11:11**-- AP-MS sets were included for clustering enrichment analysis

## environment ----
library(linkcomm)
library(openxlsx)
library(igraph)
library(dplyr)
library(stringr)
library(rstatix)

## data ----
huri <- read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = 2)
huri <- huri[, c(5:6)]
huri_graph <- simplify(graph_from_data_frame(huri, directed = FALSE), remove.loops = FALSE)
## community establishment ----
huri_ocg <- getOCG.clusters(as_data_frame(huri_graph))
# save.image(file = "../data/processed/HuRI_ocg.RData")

## save community data to XLSX ----
wb <- createWorkbook()
addWorksheet(wb, "HuRI")
writeData(wb, sheet = "HuRI", huri_ocg[4])
# saveWorkbook(wb, "../result/HuRI_assays_OCG.xlsx", overwrite = TRUE)

## generate community summary ----
clust <- huri_ocg$nodeclusters
clust2 <- clust %>% select(cluster) %>% group_by(cluster) %>% mutate(count = n()) %>% distinct() # have community size

ppi <- list(
    husci = read.csv("../data/HuSCI_node.csv", header = T),
    gordon = read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = 3),
    stukalov = read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = 4),
    li = read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = 5),
    nabeel = read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = 6)
)

ppi_node <- list(
    husci = unique(ppi[["husci"]][ppi[["husci"]]$category == "human", "node"]),
    gordon = unique(ppi[["gordon"]]$PreyGene),
    stukalov = unique(ppi[["stukalov"]]$human),
    li = unique(ppi[["li"]]$human),
    nabeel = unique(ppi[["nabeel"]]$PreyGene)
)

ppi_node[["husci"]] <- ppi_node[["husci"]][!ppi_node[["husci"]] == "TRAF1"] # TRAF1 was out!

ppi_node_df <- lapply(ppi_node, function(x) {
    df <- huri_ocg$nodeclusters[huri_ocg$nodeclusters$node %in% x, ] # viral target in community
    df2 <- df %>% select(cluster) %>% group_by(cluster) %>% mutate(count = n()) %>% distinct() # community viral target size
    return(df2)
})

ppi_node_all <- lapply(ppi_node_df, function(x) {
    merged <- merge(clust2, x, by = "cluster", all.x = TRUE) # merge community size and virtal target
    merged[is.na(merged)] <- 0
    merged$notTarget <- merged[, 2] - merged[, 3]
    row.names(merged) <- merged[, 1]
    return(merged)
})

ppi_node_toCheck <- lapply(ppi_node_all, function(x) {
    merged_checkP <- x[, c(3, 4)]
    names(merged_checkP) <- c("Target", "notTarget")
    return(merged_checkP)
})

## statistical analysis
total <- lapply(ppi_node, function(x) {
    len <- c(length(x), vcount(huri_graph))
    return(len)
})

ppi_node_stat <- list()
for (n in names(total)) {
    data <- data.frame()
    data <- apply(ppi_node_toCheck[[n]], 1, function(y) {
        tab <- rbind(y, total[[n]])
        fisher_test(tab)
        })

    pvalue <- c()
    for (i in 1:length(data)) {
        pvalue <- c(pvalue, data[[i]]$p)
    }
    padj <- p.adjust(pvalue, method = "fdr")

    enrich <- do.call(rbind.data.frame, data)
    enrich$cluster <- row.names(enrich)
    all <- cbind(ppi_node_all[[n]], enrich, padj)
    names(all) <- c("community", "community_size", "viral_target", "notTarget", "n", "p", "p.signif", "community_2", "p_adjust (fdr)")
    all$enriched <- ifelse(all$p < 0.05 & all$community_size >= 4, TRUE, FALSE)
    ppi_node_stat[[n]] <- all
}

ppi_node_stat <- lapply(ppi_node_stat, function(x) {
    x <- x[, c(1, 2, 3, 6, 7, 9, 10)]
})

write.xlsx(ppi_node_stat, file = "../result/HuRI_signif_community.xlsx", colNames = T, overwrite = T)
save.image("../result/community_stat.RData")