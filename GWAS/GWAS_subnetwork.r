# COVID-19 GWAS hit analysis
# Lin Chung-win

######
# load package
library(igraph)
library(rethinking)
library(gprofiler2)
library(gplots)
library(openxlsx)

source("~/Documents/INET-work/virus_network/src/combineNetwork.r")

huriRewire <- function(remove.loops = FALSE, ...) {
    huri_re <- rewire(huri_g, keeping_degseq(niter = gsize(huri_g) * 10))
    huri_sim <- simplify(huri_re, remove.loops = remove.loops)
    return(huri_sim)
}

huriRewireHusci <- function(node, remove.loops) {
    huri_re <- huriRewire(remove.loops)
    merged <- combineNetwork(huri_re, node)
    merged_inHuSCI <- as.numeric(table(V(merged)$name %in% husci_sym)["TRUE"])
    output <- c(merged_inHuSCI, gsize(merged))
    return(output)
}

######
# load dataset
huri <- read.xlsx("extended_table/Extended_Table_2_PPIs.xlsx", sheet = "HuRI")
husci <- read.csv("HuSCI_edge.csv", header = TRUE)
gwas <- read.csv("GWAS_fromNature_node.csv", header = T)

######
# 1. HuRI graph generation
huri_symbol <- huri[, c(5:6)]
huri_g_ori <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g_ori, remove.loops = FALSE) # V:8274, E:52558
# protein list filter
husci_sym <- husci[husci$group == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole
gwas_huri <- gwas$name[!is.na(gwas$ctl == 1)] # GWAS hit in HuRI

######
# 2. interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(huri_g, nodes = gwas_huri, order = 1, mode = "all")

######
# 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
# subnetwork of GWAS hit from HuRI
# inherit from above code
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(huri_g, names(V(gwas_all_g_merge))), remove.loops = F)
gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym]
gwas_all_husci_length <- length(gwas_all_husci)

# 3.1. do statistical analysis for all 5 GWAS hit candidate genes
gwas_rand_r2 <- c()
gwas_rand_r2 <- c(gwas_rand_r2, mcreplicate(10000, huriRewireHusci(gwas_huri, FALSE), mc.cores = detectCores()))
gwas_rand_df_r2 <- data.frame(matrix(gwas_rand_r2, ncol = 2, byrow = T))
names(gwas_rand_df_r2) <- c("viral_target", "interactions")
