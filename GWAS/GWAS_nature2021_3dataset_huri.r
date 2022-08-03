# COVID-19 GWAS hit analysis
# Lin Chung-win
# 11.08.2021
# result saved as "~/INET/GWAS_3data_HuRI.RData"

######
# load package
library(openxlsx)
library(igraph)
library(rethinking)
library(gprofiler2)
library(linkcomm)
library(gplots)

source("~/Documents/INET-work/virus_network/src/combineNetwork.r")

huriRewire <- function(remove.loops = TRUE, ...) {
    huri_re <- rewire(huri_g, keeping_degseq(niter = gsize(huri_g) * 10))
    huri_sim <- simplify(huri_re, remove.loops = remove.loops)
    return(huri_sim)
}

huriRewireMulti <- function(gwas, ctcl, hosp, infct, husci, gordon, stukalov, remove.loops = TRUE) {
    df <- c()
    re <- huriRewire(remove.loops) # rewire huri
    merged <- combineNetwork(re, gwas) # get GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"]) # get viral targets from HuSCI in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"]) # get viral targets from Gordon in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"]) # get viral targets from Stukalov in GWAS+1 subnetwork
    df <- c(df, gsize(merged)) # get network size of GWAS+1 subnetwork
    # average shortest path of GWAS subnetwork from the rewired HuRI
    df <- c(df, mean(distances(re, v = gwas, to = gwas)))

    merged <- combineNetwork(re, ctcl) # get GWAS+1 subnetwork only from critical illness candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    df <- c(df, mean(distances(re, v = ctcl, to = ctcl)))

    merged <- combineNetwork(re, hosp) # get GWAS+1 subnetwork only from hospitalized candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    df <- c(df, mean(distances(re, v = hosp, to = hosp)))

    merged <- combineNetwork(re, infct) # get GWAS+1 subnetwork only from reported infection candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    df <- c(df, mean(distances(re, v = infct, to = infct)))

    return(df)
}
plotHist <- function(value, title, length, y1, y2) {
    med <- median(value)
    dens_gwas <- hist(value, breaks = c(0:(max(value) + 1)), plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, 20), xaxt = "n", xlab = "Number of viral targets", ylab = "Frequency", main = "", cex.sub = 0.5)
    mytitle <- paste0("COVID19 GWAS subnetwork\nviral targets in ", title)
    mtext(side = 3, line = 1, cex = 1, mytitle)
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from HuRI")
    axis(side = 1, at = seq(0, 20, by = 5) + 0.5, labels = seq(0, 20, by = 5))
    arrows(length + 0.5, y1, length + 0.5, 0, col = "#922687", lwd = 2, length = 0.1)
    text(med + 2, max(dens_gwas$counts / 10000), paste0("median = ", med), col = "grey", cex = 0.5)
    text(length - 2, y2, paste0("observed = ", length, "\np = ", table(value >= length)["TRUE"] / 10000), cex = 0.4, pos = 4)
}

plotInteraction <- function(value, ymax, observe, phenotype) {
    dens_gwas <- hist(value, breaks = 20, plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of interactions", main = "", cex.sub = 0.5)
    mtext(side = 3, line = 1, cex = 1, paste0("COVID19 GWAS subnetwork: ", phenotype))
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from HuRI")
    axis(side = 2, at = seq(0, ymax, by = 500), labels = seq(0, ymax/10000, by = 0.05), las = 1)
    arrows(observe, 200, observe, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value), max(dens_gwas$counts), paste0("median = ", median(value)), col = "grey", cex = 0.5)
    text(observe - (observe / 10), 350, paste0("observed = ", observe, "\np = ", table(value >= observe)["TRUE"]/10000), cex = 0.4, pos = 4)
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
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
gwas <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits_v2.xlsx") # ref: Mapping the human genetic architecture of COVID-19 (url: https://doi.org/10.1038/s41586-021-03767-x)
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)

# add Gordon and Stukalov data
gordon <- read.xlsx("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Gordon")
stukalov <- read.xlsx("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Stukalov")

######
# 1. HuRI graph generation
huri_symbol <- huri[, c(5:6)]
huri_g_ori <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g_ori, remove.loops = TRUE) # V:8274, E:52558
# protein list filter
husci_sym <- husci[husci$group == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole

######
# all GWAS proteins
gwas_huri <- gwas$All.LD[gwas$All.LD %in% V(huri_g)$name] # GWAS hit in HuRI

# Gordon and Stukalov in HuRI
gordon_sym <- unique(gordon$PreyGene)
gordon_huri <- V(huri_g)$name[V(huri_g)$name %in% gordon_sym]

stukalov_sym <- unique(stukalov$human)
stukalov_huri <- V(huri_g)$name[V(huri_g)$name %in% stukalov_sym]

######
# 2. interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(huri_g, nodes = sort(gwas_huri2), order = 1, mode = "all") #14 of 42 in HuRI

# 2.1 interactors of GWAS hits with critical illness phenotypes
ctcl <- gwas[, 2][gwas[, 5] == 1]
ctcl <- unique(ctcl[!is.na(ctcl)])

ctcl_huri <- ctcl[ctcl %in% V(huri_g)$name]

ctcl_1st <- combineNetwork(huri_g, ctcl_huri)
gwas_ctcl_husci <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% husci_sym]
gwas_ctcl_husci_length <- length(gwas_ctcl_husci)

gwas_ctcl_gordon <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% gordon_sym]
gwas_ctcl_gordon_length <- length(gwas_ctcl_gordon)
gwas_ctcl_stukalov <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% stukalov_sym]
gwas_ctcl_stukalov_length <- length(gwas_ctcl_stukalov)

hosp <- gwas[, 2][gwas[, 6] == 1]
hosp <- unique(hosp[!is.na(hosp)])

hosp_huri <- hosp[hosp %in% V(huri_g)$name]

hosp_1st <- combineNetwork(huri_g, hosp_huri)
gwas_hosp_husci <- V(hosp_1st)$name[V(hosp_1st)$name %in% husci_sym]
gwas_hosp_husci_length <- length(gwas_hosp_husci)

gwas_hosp_gordon <- V(hosp_1st)$name[V(hosp_1st)$name %in% gordon_sym]
gwas_hosp_gordon_length <- length(gwas_hosp_gordon)
gwas_hosp_stukalov <- V(hosp_1st)$name[V(hosp_1st)$name %in% stukalov_sym]
gwas_hosp_stukalov_length <- length(gwas_hosp_stukalov)

infct <- gwas[, 2][gwas[, 7] == 1]
infct <- unique(infct[!is.na(infct)])

infct_huri <- infct[infct %in% V(huri_g)$name]

infct_1st <- combineNetwork(huri_g, infct_huri)
gwas_infct_husci <- V(infct_1st)$name[V(infct_1st)$name %in% husci_sym]
gwas_infct_husci_length <- length(gwas_infct_husci)

gwas_infct_gordon <- V(infct_1st)$name[V(infct_1st)$name %in% gordon_sym]
gwas_infct_gordon_length <- length(gwas_infct_gordon)
gwas_infct_stukalov <- V(infct_1st)$name[V(infct_1st)$name %in% stukalov_sym]
gwas_infct_stukalov_length <- length(gwas_infct_stukalov)

# average shortest path
all_path <- mean(distances(huri_g, v = gwas_huri, to = gwas_huri)[lower.tri(distances(huri_g, v = gwas_huri, to = gwas_huri))])
ctcl_path <- mean(distances(huri_g, v = ctcl_huri, to = ctcl_huri)[lower.tri(distances(huri_g, v = ctcl_huri, to = ctcl_huri))])
hosp_path <- mean(distances(huri_g, v = hosp_huri, to = hosp_huri)[lower.tri(distances(huri_g, v = hosp_huri, to = hosp_huri))])
infct_path <- mean(distances(huri_g, v = infct_huri, to = infct_huri)[lower.tri(distances(huri_g, v = infct_huri, to = infct_huri))])

# 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
# load gwas loci info, with 3 phenotype (critical illness, hospitalization and infection)
# subnetwork of GWAS hit from HuRI
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(huri_g, names(V(gwas_all_g_merge))), remove.loops = TRUE)

gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym] # 11 viral targets
gwas_all_husci_length <- length(gwas_all_husci)

gwas_all_gordon <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% gordon_sym] # 4 viral targets
gwas_all_gordon_length <- length(gwas_all_gordon)
gwas_all_stukalov <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% stukalov_sym] # 9 viral targets
gwas_all_stukalov_length <- length(gwas_all_stukalov)

# 3.1. do statistical analysis for all 17 GWAS hit candidate genes # time consuming
# HuSCI, Gordon and Stukalov
all_re <- c()
all_re <- c(all_re, mcreplicate(10000, huriRewireMulti(gwas_huri, ctcl_huri, hosp_huri, infct_huri, husci_sym, gordon_sym, stukalov_sym), mc.cores = detectCores()))
all_re[is.na(all_re)] <- 0

all_re_df <- data.frame(matrix(all_re, ncol = 20, byrow = T))
all_re_df_n <- c(
    "allGWAS_viral_target_inHuSCI",
    "allGWAS_viral_target_inGordon",
    "allGWAS_viral_target_inStukalov",
    "allGWAS_subnetworkSize",
    "allGWAS_shortestpath",
    "ctclGWAS_viral_target_inHuSCI",
    "ctclGWAS_viral_target_inGordon",
    "ctclGWAS_viral_target_inStukalov",
    "ctclGWAS_subnetworkSize",
    "ctclGWAS_shortestpath",
    "hospGWAS_viral_target_inHuSCI",
    "hospGWAS_viral_target_inGordon",
    "hospGWAS_viral_target_inStukalov",
    "hospGWAS_subnetworkSize",
    "hospGWAS_shortestpath",
    "infctGWAS_viral_target_inHuSCI",
    "infctGWAS_viral_target_inGordon",
    "infctGWAS_viral_target_inStukalov",
    "infctGWAS_subnetworkSize",
    "infctGWAS_shortestpath"
)
names(all_re_df) <- all_re_df_n

# write.xlsx(all_re_df, file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/Nature2021b_3dataset_paralog.xlsx", overwrite = T)
all_re_df_plot <- all_re_df[, c(1:3, 6:8, 11:13, 16:18)]
all_length <- c(gwas_all_husci_length,
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
# phenotype <- rep(c(
#     "all 17 genes and 1st interactors",
#     "critical illness, 10 genes and 1st interactors",
#     "hospitalization, 13 genes and 1st interactors",
#     "reported infection, 6 genes and 1st interactors"
#     ), each = 3)

# phenotype2 <- rep(c(
#     "all 14 genes",
#     "critical illness, 7 genes",
#     "hospitalization, 10 genes",
#     "reported infection, 5 genes"
#     ), each = 3)

y1 <- c(0.03, 0.03, 0.03, 0.03, 0.04, 0.03, 0.03, 0.03, 0.03, 0.04, 0.05, 0.05)
y2 <- c(0.05, 0.05, 0.05, 0.05, 0.06, 0.05, 0.05, 0.05, 0.05, 0.06, 0.07, 0.07)
# plotting
pdf("Nature2021b_3dataset_HuRI.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
for (i in 1:12) {
    plotHist(
        all_re_df_plot[, i],
        title[i],
        # phenotype2[i],
        all_length[i],
        y1[i], y2[i]
        )
}

# interaction, all GWAS
plotInteraction(all_re_df[, 4], 1500, gsize(gwas_all_final), "all GWAS")
# interaction, all Critical illness
plotInteraction(all_re_df[, 9], 1000, gsize(ctcl_1st), "critical")
# interaction, all hospitalization
plotInteraction(all_re_df[, 14], 1000, gsize(hosp_1st), "hospitalization")
# interaction, all reported infection
plotInteraction(all_re_df[, 19], 1000, gsize(infct_1st), "infection")

# average shortest path
plotDistance(all_re_df[, 5], 3000, all_path, "all GWAS")
plotDistance(all_re_df[, 10], 4000, ctcl_path, "critical illness")
plotDistance(all_re_df[, 15], 2500, hosp_path, "hospitalization")
plotDistance(all_re_df[, 20], 3500, infct_path, "infection")
dev.off()

############################################################
# exclude paralogs
# inherit from above code
############################################################
gwas_huri2 <- gwas_huri[c(1, 3, 5:9, 11:17)] # omit OAS2, ICAM1 and ICAM4
ctcl_huri2 <- ctcl_huri[c(1, 3, 5:7, 9, 10)] # only OAS1, ICAM3

ctcl_1st <- combineNetwork(huri_g, ctcl_huri2)
gwas_ctcl_husci <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% husci_sym]
gwas_ctcl_husci_length <- length(gwas_ctcl_husci)
gwas_ctcl_gordon <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% gordon_sym]
gwas_ctcl_gordon_length <- length(gwas_ctcl_gordon)
gwas_ctcl_stukalov <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% stukalov_sym]
gwas_ctcl_stukalov_length <- length(gwas_ctcl_stukalov)

hosp_huri2 <- hosp_huri[c(1, 3, 5:7, 9:13)] # only OAS1, ICAM3
hosp_1st <- combineNetwork(huri_g, hosp_huri2)
gwas_hosp_husci <- V(hosp_1st)$name[V(hosp_1st)$name %in% husci_sym]
gwas_hosp_husci_length <- length(gwas_hosp_husci)
gwas_hosp_gordon <- V(hosp_1st)$name[V(hosp_1st)$name %in% gordon_sym]
gwas_hosp_gordon_length <- length(gwas_hosp_gordon)
gwas_hosp_stukalov <- V(hosp_1st)$name[V(hosp_1st)$name %in% stukalov_sym]
gwas_hosp_stukalov_length <- length(gwas_hosp_stukalov)

infct_huri2 <- infct_huri[c(1:3, 5, 6)] # only OAS1
infct_1st <- combineNetwork(huri_g, infct_huri)
gwas_infct_husci <- V(infct_1st)$name[V(infct_1st)$name %in% husci_sym]
gwas_infct_husci_length <- length(gwas_infct_husci)

gwas_infct_gordon <- V(infct_1st)$name[V(infct_1st)$name %in% gordon_sym]
gwas_infct_gordon_length <- length(gwas_infct_gordon)
gwas_infct_stukalov <- V(infct_1st)$name[V(infct_1st)$name %in% stukalov_sym]
gwas_infct_stukalov_length <- length(gwas_infct_stukalov)

all_path <- mean(distances(huri_g, v = gwas_huri2, to = gwas_huri2)[lower.tri(distances(huri_g, v = gwas_huri2, to = gwas_huri2))])
ctcl_path <- mean(distances(huri_g, v = ctcl_huri2, to = ctcl_huri2)[lower.tri(distances(huri_g, v = ctcl_huri2, to = ctcl_huri2))])
hosp_path <- mean(distances(huri_g, v = hosp_huri2, to = hosp_huri2)[lower.tri(distances(huri_g, v = hosp_huri2, to = hosp_huri2))])
infct_path <- mean(distances(huri_g, v = infct_huri2, to = infct_huri2)[lower.tri(distances(huri_g, v = infct_huri2, to = infct_huri2))])

gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(huri_g, names(V(gwas_all_g_merge))), remove.loops = TRUE)

gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym] # 11 viral targets
gwas_all_husci_length <- length(gwas_all_husci)

gwas_all_gordon <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% gordon_sym] # 4 viral targets
gwas_all_gordon_length <- length(gwas_all_gordon)
gwas_all_stukalov <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% stukalov_sym] # 9 viral targets
gwas_all_stukalov_length <- length(gwas_all_stukalov)

all_re2 <- c()
all_re2 <- c(all_re2, mcreplicate(10000, huriRewireMulti(gwas_huri2, ctcl_huri2, hosp_huri2, infct_huri2, husci_sym, gordon_sym, stukalov_sym), mc.cores = detectCores()))
all_re2[is.na(all_re2)] <- 0

all_re_df <- data.frame(matrix(all_re2, ncol = 20, byrow = T))
names(all_re_df) <- all_re_df_n

all_re_df_plot <- all_re_df[, c(1:3, 6:8, 11:13, 16:18)]
all_length <- c(gwas_all_husci_length,
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

pdf("Nature2021b_3dataset_HuRI_paralogs.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
for (i in 1:12) {
    plotHist(
        all_re_df_plot[, i],
        title[i],
        # phenotype2[i],
        all_length[i],
        y1[i], y2[i]
        )
}

# interaction, all GWAS
plotInteraction(all_re_df[, 4], 1500, gsize(gwas_all_final), "all GWAS")
# interaction, all Critical illness
plotInteraction(all_re_df[, 9], 1000, gsize(ctcl_1st), "critical")
# interaction, all hospitalization
plotInteraction(all_re_df[, 14], 1000, gsize(hosp_1st), "hospitalization")
# interaction, all reported infection
plotInteraction(all_re_df[, 19], 1000, gsize(infct_1st), "infection")

# average shortest path
plotDistance(all_re_df[, 5], 3000, all_path, "all GWAS")
plotDistance(all_re_df[, 10], 4000, ctcl_path, "critical illness")
plotDistance(all_re_df[, 15], 2500, hosp_path, "hospitalization")
plotDistance(all_re_df[, 20], 3500, infct_path, "infection")
dev.off()

######
# save result
save.image("Nature2021b_3data_HuRI.RData")
