# COVID-19 GWAS hit analysis with BioPlex 3.0
# meta-analysis: https://www.nature.com/articles/s41586-021-03767-x
# Mapping the human genetic architecture of COVID-19
# plus Gordon and Stukalov dataset
# Lin Chung-wen
# Date: 11.08.2021

######
# load package
library(igraph)
library(rethinking)
library(gprofiler2)
library(gplots)
library(openxlsx)
library(plotrix) # add table to plot

source("~/Documents/INET-work/virus_network/src/combineNetwork.r")

bioplexRewire <- function(remove.loops = TRUE, ...) {
    bioplex_re <- rewire(bioplex_g, keeping_degseq(niter = gsize(bioplex_g) * 10))
    bioplex_sim <- simplify(bioplex_re, remove.loops = remove.loops)
    return(bioplex_sim)
}

bioplexRewireMulti <- function(gwas, ctcl, hosp, infct, husci, gordon, stukalov, remove.loops = TRUE) {
    df <- c()
    re <- bioplexRewire(remove.loops) # rewire BioPlex
    merged <- combineNetwork(re, gwas) # get GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"]) # get viral targets from HuSCI in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"]) # get viral targets from Gordon in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"]) # get viral targets from Stukalov in GWAS+1 subnetwork
    df <- c(df, gsize(merged)) # get network size of GWAS+1 subnetwork

    merged <- combineNetwork(re, ctcl) # get GWAS+1 subnetwork only from critical illness candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))

    merged <- combineNetwork(re, hosp) # get GWAS+1 subnetwork only from hospitalized candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))

    merged <- combineNetwork(re, infct) # get GWAS+1 subnetwork only from reported infection candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))

    df <- c(df, mean(distances(re, v = gwas, to = gwas))[lower.tri(distances(re, v = gwas, to = gwas))])
    df <- c(df, mean(distances(re, v = ctcl, to = ctcl)[lower.tri(distances(re, v = ctcl, to = ctcl))]))
    df <- c(df, mean(distances(re, v = hosp, to = hosp)[lower.tri(distances(re, v = hosp, to = hosp))]))
    df <- c(df, mean(distances(re, v = infct, to = infct)[lower.tri(distances(re, v = infct, to = infct))]))
    return(df)
}

plotHist <- function(value, title, length, xmax, y1, y2) {
    dens_gwas <- hist(value, breaks = c(0:(max(value) + 1)), plot = FALSE, right = F)
    plot(dens_gwas, xlim = c(0, xmax), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", freq = FALSE, xlab = "Number of viral targets", ylab = "Frequency", main = "", cex.sub = 0.5)
    mytitle <- paste0("COVID19 GWAS subnetwork\nviral targets in ", title)
    mtext(side = 3, line = 1, cex = 1, mytitle)
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from BioPlex3.0")
    axis(side = 1, at = seq(0, xmax, by = 5) + 0.5, labels = seq(0, xmax, by = 5))
    arrows(length + 0.5, y1, length + 0.5, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value) + 4, max(dens_gwas$counts / 10000), paste0("median = ", median(value)), col = "grey", cex = 0.5)
    text(length - 2, y2, paste0("observed = ", length, "\np = ", table(value >= length)["TRUE"]/10000), cex = 0.4, pos = 4)
}

plotInteraction <- function(value, ymax, observe, phenotype) {
    dens_gwas <- hist(value, breaks = 20, plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), xlim = c(min(value), ifelse(observe > max(value), observe, max(value))), border = NA, las = 1, yaxt = "n", xlab = "Number of interactions", main = "", cex.sub = 0.5)
    mtext(side = 3, line = 1, cex = 1, paste0("COVID19 GWAS subnetwork: ", phenotype))
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from BioPlex3")
    axis(side = 2, at = seq(0, ymax, by = 500), labels = seq(0, ymax/10000, by = 0.05), las = 1)
    arrows(observe, 200, observe, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value), max(dens_gwas$counts), paste0("median = ", median(value)), col = "grey", cex = 0.5)
    text(observe - (observe / 9), 350, paste0("observed = ", observe, "\np ", ifelse(is.na(as.numeric(table(value >= observe)["TRUE"])), "< 0.0001", paste0(" = ", table(value >= observe)["TRUE"]/10000))), cex = 0.4, pos = 4)
}
plotDistance <- function(value, ymax, observe, phenotype) {
    dens_gwas <- hist(value, breaks = 8, plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xlim = c(2, 4), yaxt = "n", xlab = "Average shortest path", main = "", cex.sub = 0.5)
    mtext(side = 3, line = 1, cex = 1, paste0("COVID19 GWAS subnetwork: ", phenotype))
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from HuRI")
    axis(side = 2, at = seq(0, ymax, by = 500), labels = seq(0, ymax/10000, by = 0.05), las = 1)
    arrows(observe, 300, observe, 0, col = "#922687", lwd = 2, length = 0.1)
    text(round(median(value), 2), max(dens_gwas$counts), paste0("median = ", round(median(value), 2)), col = "grey", cex = 0.5)
    text(round(observe, 1) - 0.3, 400, paste0("observed = ", round(observe, 2), "\np = ", table(value >= round(observe, 2))["FALSE"]/10000), cex = 0.4, pos = 4)
}
######
# load dataset
bioplex <- read.delim("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/BioPlex.3.0_edge.tsv", header = T)
husci <- read.csv("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/HuSCI_node.csv", header = TRUE)
gwas <- read.xlsx("COVID_GWAS_hit_inHUSCI_v2.xlsx")
gordon <- read.xlsx("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Gordon")
stukalov <- read.xlsx("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Stukalov")

######
# 1. BioPlex graph generation
bioplex_symbol <- bioplex[, c(5:6)]
bioplex_g <- graph_from_data_frame(bioplex_symbol, directed = FALSE) # V:13957, E:118162

# GWAS list
gwas_bioplex <- gwas$All.LD[gwas$All.LD %in% V(bioplex_g)$name] # 24 of 42 candidates found in BioPlex3.0

ctcl <- gwas[, 1][gwas[, 6] == 1]
ctcl <- unique(ctcl[!is.na(ctcl)])
ctcl_bioplex <- ctcl[ctcl %in% V(bioplex_g)$name] # V:10
ctcl_1st <- combineNetwork(bioplex_g, ctcl_bioplex)

hosp <- gwas[, 1][gwas[, 7] == 1]
hosp <- unique(hosp[!is.na(hosp)])
hosp_bioplex <- hosp[hosp %in% V(bioplex_g)$name] # V:17
hosp_1st <- combineNetwork(bioplex_g, hosp_bioplex)

infct <- gwas[, 1][gwas[, 8] == 1]
infct <- unique(infct[!is.na(infct)])
infct_bioplex <- infct[infct %in% V(bioplex_g)$name] # V:10
infct_1st <- combineNetwork(bioplex_g, infct_bioplex)

# HuSCI, Gordon and Stukalov in BioPlex
husci_sym <- husci$node # V:171
husci_bioplex <- V(bioplex_g)$name[V(bioplex_g)$name %in% husci_sym] # HuSCI in BioPlex whole, V:132

gordon_sym <- unique(gordon$PreyGene) # V:384
gordon_bioplex <- V(bioplex_g)$name[V(bioplex_g)$name %in% gordon_sym] # V:346

stukalov_sym <- unique(stukalov$human) # V:876
stukalov_bioplex <- V(bioplex_g)$name[V(bioplex_g)$name %in% stukalov_sym] # V:723

######
# 2. interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(bioplex_g, nodes = gwas_bioplex, order = 1, mode = "all")

######
# 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
# subnetwork of GWAS hit from BioPlex
# inherit from above code
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(bioplex_g, names(V(gwas_all_g_merge))), remove.loops = TRUE)

# GWAS hit in HuSCI
gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym] # V:5
gwas_all_husci_length <- length(gwas_all_husci)
gwas_ctcl_husci <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% husci_sym] # V:2
gwas_ctcl_husci_length <- length(gwas_ctcl_husci)
gwas_hosp_husci <- V(hosp_1st)$name[V(hosp_1st)$name %in% husci_sym] # V:2
gwas_hosp_husci_length <- length(gwas_hosp_husci)
gwas_infct_husci <- V(infct_1st)$name[V(infct_1st)$name %in% husci_sym] # V:4
gwas_infct_husci_length <- length(gwas_infct_husci)

# GWAS hit in Gordon
gwas_all_gordon <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% gordon_sym] # V:10
gwas_all_gordon_length <- length(gwas_all_gordon)
gwas_ctcl_gordon <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% gordon_sym]
gwas_ctcl_gordon_length <- length(gwas_ctcl_gordon)
gwas_hosp_gordon <- V(hosp_1st)$name[V(hosp_1st)$name %in% gordon_sym]
gwas_hosp_gordon_length <- length(gwas_hosp_gordon)
gwas_infct_gordon <- V(infct_1st)$name[V(infct_1st)$name %in% gordon_sym]
gwas_infct_gordon_length <- length(gwas_infct_gordon)

# GWAS hit in Stukalov
gwas_all_stukalov <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% stukalov_sym] # V:19
gwas_all_stukalov_length <- length(gwas_all_stukalov)
gwas_ctcl_stukalov <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% stukalov_sym]
gwas_ctcl_stukalov_length <- length(gwas_ctcl_stukalov)
gwas_hosp_stukalov <- V(hosp_1st)$name[V(hosp_1st)$name %in% stukalov_sym]
gwas_hosp_stukalov_length <- length(gwas_hosp_stukalov)
gwas_infct_stukalov <- V(infct_1st)$name[V(infct_1st)$name %in% stukalov_sym]
gwas_infct_stukalov_length <- length(gwas_infct_stukalov)

# average shortest path of GWAS proteins
gwas_all_path <- mean(distances(bioplex_g, v = gwas_bioplex, to = gwas_bioplex)[lower.tri(distances(bioplex_g, v = gwas_bioplex, to = gwas_bioplex))])
gwas_ctcl_path <- mean(distances(bioplex_g, v = ctcl_bioplex, to = ctcl_bioplex)[lower.tri(distances(bioplex_g, v = ctcl_bioplex, to = ctcl_bioplex))])
gwas_hosp_path <- mean(distances(bioplex_g, v = hosp_bioplex, to = hosp_bioplex)[lower.tri(distances(bioplex_g, v = hosp_bioplex, to = hosp_bioplex))])
gwas_infct_path <- mean(distances(bioplex_g, v = infct_bioplex, to = infct_bioplex)[lower.tri(distances(bioplex_g, v = infct_bioplex, to = infct_bioplex))])

######
# permutation analysis
gwas_rand_r2 <- c()
gwas_rand_r2 <- c(gwas_rand_r2, mcreplicate(10000, bioplexRewireMulti(gwas_bioplex, ctcl_bioplex, hosp_bioplex, infct_bioplex, husci_sym, gordon_sym, stukalov_sym), mc.cores = detectCores()))
gwas_rand_r2[is.na(gwas_rand_r2)] <- 0

gwas_rand_df_r2 <- data.frame(matrix(gwas_rand_r2, ncol = 20, byrow = T))
gwas_rand_df_name <- c(
    "allGWAS_viral_target_inHuSCI",
    "allGWAS_viral_target_inGordon",
    "allGWAS_viral_target_inStukalov",
    "allGWAS_subnetworkSize",
    "ctclGWAS_viral_target_inHuSCI",
    "ctclGWAS_viral_target_inGordon",
    "ctclGWAS_viral_target_inStukalov",
    "ctclGWAS_subnetworkSize",
    "hospGWAS_viral_target_inHuSCI",
    "hospGWAS_viral_target_inGordon",
    "hospGWAS_viral_target_inStukalov",
    "hospGWAS_subnetworkSize",
    "infctGWAS_viral_target_inHuSCI",
    "infctGWAS_viral_target_inGordon",
    "infctGWAS_viral_target_inStukalov",
    "infctGWAS_subnetworkSize",
    "GWAS_average_shortest_path",
    "GWAS_ctcl_avg_path",
    "GWAS_hosp_avg_path",
    "GWAS_infct_avg_path"
)

names(gwas_rand_df_r2) <- gwas_rand_df_name

all_re_df_plot <- gwas_rand_df_r2[, c(1:3, 5:7, 9:11, 13:15)]
all_length <- c(
    gwas_all_husci_length,
    gwas_all_gordon_length,
    gwas_all_stukalov_length,
    gwas_ctcl_husci_length,
    gwas_ctcl_gordon_length,
    gwas_ctcl_stukalov_length,
    gwas_hosp_husci_length,
    gwas_hosp_gordon_length,
    gwas_hosp_stukalov_length,
    gwas_infct_husci_length,
    gwas_infct_gordon_length,
    gwas_infct_stukalov_length
)
title <- rep(c("HuSCI", "Gordon et al", "Stukalov et al"), 4)
phenotype <- rep(c(
    "all 24 genes",
    "critical illness, 10 genes",
    "hospitalization, 17 genes",
    "reported infection, 10 genes"
    ), each = 3)

phenotype2 <- rep(c(
    "all 24 genes",
    "critical illness, 10 genes",
    "hospitalization, 17 genes",
    "reported infection, 10 genes"
    ), each = 3)
xmax <- c(20, 25, 45, 15, 20, 30, 15, 20, 30, 15, 20, 30)

######
# plotting
pdf(file = "Nature2021b_3dataset_BioPlex3.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
for (i in 1:12) {
    plotHist(
        all_re_df_plot[, i],
        title[i],
        # phenotype2[i],
        all_length[i],
        xmax[i],
        0.03, 0.05
        )
}
# interaction
plotInteraction(gwas_rand_df_r2[, 4], 1500, gsize(gwas_all_final), "all GWAS")
plotInteraction(gwas_rand_df_r2[, 8], 1500, gsize(ctcl_1st), "critical illness")
plotInteraction(gwas_rand_df_r2[, 12], 1500, gsize(hosp_1st), "hospitalization")
plotInteraction(gwas_rand_df_r2[, 16], 1000, gsize(infct_1st), "reported infection")
# average shortest path
plotDistance(gwas_rand_df_r2[, 17], 3000, gwas_all_path, "all GWAS")
plotDistance(gwas_rand_df_r2[, 18], 3500, gwas_ctcl_path, "critical illness")
plotDistance(gwas_rand_df_r2[, 19], 4000, gwas_hosp_path, "hospitalization")
plotDistance(gwas_rand_df_r2[, 20], 4000, gwas_infct_path, "reported infection")
dev.off()

######################################################
# exclude paralogs
######################################################
gwas_bioplex2 <- gwas_bioplex[c(1:5, 7, 10:24)]
ctcl_bioplex2 <- ctcl_bioplex[c(1, 3, 6:10)]
ctcl_1st <- combineNetwork(bioplex_g, ctcl_bioplex2)
hosp_bioplex2 <- hosp_bioplex[c(1:4, 6, 9:17)]
hosp_1st <- combineNetwork(bioplex_g, hosp_bioplex2)
infct_bioplex2 <- infct_bioplex[c(1:10)]
infct_1st <- combineNetwork(bioplex_g, infct_bioplex2)

gwas_hit_1st <- make_ego_graph(bioplex_g, nodes = gwas_bioplex2, order = 1, mode = "all")
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
gwas_all_final <- simplify(induced_subgraph(bioplex_g, names(V(gwas_all_g_merge))), remove.loops = TRUE)

gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym] # V:5
gwas_all_husci_length <- length(gwas_all_husci)
gwas_ctcl_husci <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% husci_sym] # V:2
gwas_ctcl_husci_length <- length(gwas_ctcl_husci)
gwas_hosp_husci <- V(hosp_1st)$name[V(hosp_1st)$name %in% husci_sym] # V:2
gwas_hosp_husci_length <- length(gwas_hosp_husci)
gwas_infct_husci <- V(infct_1st)$name[V(infct_1st)$name %in% husci_sym] # V:4
gwas_infct_husci_length <- length(gwas_infct_husci)

gwas_all_gordon <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% gordon_sym] # V:10
gwas_all_gordon_length <- length(gwas_all_gordon)
gwas_ctcl_gordon <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% gordon_sym]
gwas_ctcl_gordon_length <- length(gwas_ctcl_gordon)
gwas_hosp_gordon <- V(hosp_1st)$name[V(hosp_1st)$name %in% gordon_sym]
gwas_hosp_gordon_length <- length(gwas_hosp_gordon)
gwas_infct_gordon <- V(infct_1st)$name[V(infct_1st)$name %in% gordon_sym]
gwas_infct_gordon_length <- length(gwas_infct_gordon)

gwas_all_stukalov <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% stukalov_sym] # V:19
gwas_all_stukalov_length <- length(gwas_all_stukalov)
gwas_ctcl_stukalov <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% stukalov_sym]
gwas_ctcl_stukalov_length <- length(gwas_ctcl_stukalov)
gwas_hosp_stukalov <- V(hosp_1st)$name[V(hosp_1st)$name %in% stukalov_sym]
gwas_hosp_stukalov_length <- length(gwas_hosp_stukalov)
gwas_infct_stukalov <- V(infct_1st)$name[V(infct_1st)$name %in% stukalov_sym]
gwas_infct_stukalov_length <- length(gwas_infct_stukalov)

gwas_all_path <- mean(distances(bioplex_g, v = gwas_bioplex, to = gwas_bioplex)[lower.tri(distances(bioplex_g, v = gwas_bioplex, to = gwas_bioplex))])
gwas_ctcl_path <- mean(distances(bioplex_g, v = ctcl_bioplex, to = ctcl_bioplex)[lower.tri(distances(bioplex_g, v = ctcl_bioplex, to = ctcl_bioplex))])
gwas_hosp_path <- mean(distances(bioplex_g, v = hosp_bioplex, to = hosp_bioplex)[lower.tri(distances(bioplex_g, v = hosp_bioplex, to = hosp_bioplex))])
gwas_infct_path <- mean(distances(bioplex_g, v = infct_bioplex, to = infct_bioplex)[lower.tri(distances(bioplex_g, v = infct_bioplex, to = infct_bioplex))])

gwas_rand_r3 <- c()
gwas_rand_r3 <- c(gwas_rand_r3, mcreplicate(10000, bioplexRewireMulti(gwas_bioplex2, ctcl_bioplex2, hosp_bioplex2, infct_bioplex2, husci_sym, gordon_sym, stukalov_sym), mc.cores = detectCores()))
gwas_rand_r3[is.na(gwas_rand_r3)] <- 0
gwas_rand_df_r3 <- data.frame(matrix(gwas_rand_r2, ncol = 20, byrow = T))
names(gwas_rand_df_r3) <- gwas_rand_df_name
all_re_df_plot <- gwas_rand_df_r3[, c(1:3, 5:7, 9:11, 13:15)]
all_length <- c(
    gwas_all_husci_length,
    gwas_all_gordon_length,
    gwas_all_stukalov_length,
    gwas_ctcl_husci_length,
    gwas_ctcl_gordon_length,
    gwas_ctcl_stukalov_length,
    gwas_hosp_husci_length,
    gwas_hosp_gordon_length,
    gwas_hosp_stukalov_length,
    gwas_infct_husci_length,
    gwas_infct_gordon_length,
    gwas_infct_stukalov_length
)

pdf(file = "Nature2021b_3dataset_BioPlex3_paralogs.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
for (i in 1:12) {
    plotHist(
        all_re_df_plot[, i],
        title[i],
        # phenotype2[i],
        all_length[i],
        xmax[i],
        0.03, 0.05
        )
}
# interaction
plotInteraction(gwas_rand_df_r3[, 4], 1500, gsize(gwas_all_final), "all GWAS")
plotInteraction(gwas_rand_df_r3[, 8], 1500, gsize(ctcl_1st), "critical illness")
plotInteraction(gwas_rand_df_r3[, 12], 1500, gsize(hosp_1st), "hospitalization")
plotInteraction(gwas_rand_df_r3[, 16], 1000, gsize(infct_1st), "reported infection")
# average shortest path
plotDistance(gwas_rand_df_r3[, 17], 3000, gwas_all_path, "all GWAS")
plotDistance(gwas_rand_df_r3[, 18], 3500, gwas_ctcl_path, "critical illness")
plotDistance(gwas_rand_df_r3[, 19], 4000, gwas_hosp_path, "hospitalization")
plotDistance(gwas_rand_df_r3[, 20], 4000, gwas_infct_path, "reported infection")
dev.off()

######
# save workarea data
save.image("Nature2021b_3dataset_BioPlex_noLoop.RData")
