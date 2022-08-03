# COVID-19 GWAS hit analysis
# plus Gordon and Stukalov dataset
# Lin Chung-wen
# Date: 28.07.2021 **23:49**

######
# load package
library(igraph)
library(rethinking)
library(gprofiler2)
library(gplots)
library(openxlsx)
library(plotrix) # add table to plot

subnetwork <- function(network, node) {
    gwas_hit_1st <- make_ego_graph(network, nodes = node, order = 1, mode = "all")
    ######
    # 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
    # subnetwork of GWAS hit from HuRI
    # inherit from above code
    gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
    gwas_all_df <- do.call(rbind, gwas_all_list_df)
    gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
    # to have interaction between 1st interactors
    gwas_all_final <- simplify(induced_subgraph(network, names(V(gwas_all_g_merge))), remove.loops = TRUE)
    return(gwas_all_final)
}

rewire3Dataset <- function(network, node) {
    count <- c()
    re <- rewire(network, keeping_degseq(niter = gsize(network) * 10))
    merged <- subnetwork(re, node)
    # merged_inHuSCI
    count <- c(count, as.numeric(table(V(merged)$name %in% husci_sym)["TRUE"]))
    # merged_inGordon
    count <- c(count, as.numeric(table(V(merged)$name %in% gordon_sym)["TRUE"]))
    # merged_inStukalov
    count <- c(count, as.numeric(table(V(merged)$name %in% stukalov_sym)["TRUE"]))
    count <- c(count, gsize(merged))
    dist <- distances(re, v = node, to = node, mode = "all")
    count <- c(count, mean(dist[lower.tri(dist)]))
    return(count)
}

plotHist <- function(value, title, length, xmax, y1, y2) {
    dens_gwas <- hist(value, breaks = c(0:(max(value) + 1)), plot = FALSE, right = F)
    plot(dens_gwas, xlim = c(0, 25), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", freq = FALSE, xlab = "Number of viral targets", ylab = "Frequency", main = "", cex.sub = 0.5)
    mytitle <- paste0("viral targets in ", title)
    mtext(side = 3, line = 1, cex = 1, mytitle)
    mtext(side = 3, line = 0.2, cex = .8, "subnetwork extracted from HuRI")
    axis(side = 1, at = seq(0, xmax, by = 5) + 0.5, labels = seq(0, xmax, by = 5))
    arrows(length + 0.5, y1, length + 0.5, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value) + 4, max(dens_gwas$counts / 10000), paste0("median = ", median(value)), col = "grey", cex = 0.5)
    text(length - 2, y2, paste0("observed = ", length, "\np = ", table(value >= length)["TRUE"]/10000), cex = 0.4, pos = 4)
}
######
# load dataset
huri <- read.xlsx("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "HuRI")
husci <- read.csv("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/HuSCI_node.csv", header = TRUE)
gwas <- read.csv("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/GWAS_hits.csv", header = T)
gordon <- read.xlsx("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Gordon")
stukalov <- read.xlsx("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Stukalov")

######
# 1. HuRI graph generation
huri_symbol <- huri[, c(5:6)]
huri_g_ori <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g_ori, remove.loops = TRUE) # V:8274, E:52558
# protein list filter
husci_sym <- husci$node
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole

# GWAS hit in HuRI
gwas_huri <- gwas$name[!is.na(gwas$ctl == 1)]
gwas_huri_paralogs <- gwas_huri[c(1, 3:5)] # exclude paralogs: OAS2, 23.08.2021

# Gordon and Stukalov in HuRI
gordon_sym <- unique(gordon$PreyGene)
gordon_huri <- V(huri_g)$name[V(huri_g)$name %in% gordon_sym]

stukalov_sym <- unique(stukalov$human)
stukalov_huri <- V(huri_g)$name[V(huri_g)$name %in% stukalov_sym]

######
# 2. observation
# subnetwork establishment
observation_all <- subnetwork(huri_g, gwas_huri)
print("GWAS subnetwork")
observation_all

observation_paralogs <- subnetwork(huri_g, gwas_huri_paralogs)
print("GWAS subnetwork, without paralogs")
observation_paralogs

# a. viral targets within 3 dataset: HuSCI, Gordon et al, Stukalov et al
husci_viral_targets_all <- V(observation_all)$name[V(observation_all)$name %in% husci_sym]
print("Viral targets from HuSCI in GWAS subnetwork")
husci_viral_targets_all

gordon_viral_targets_all <- V(observation_all)$name[V(observation_all)$name %in% gordon_sym]
print("Viral targets from Gordon at al in GWAS subnetwork")
gordon_viral_targets_all

stukalov_viral_targets_all <- V(observation_all)$name[V(observation_all)$name %in% stukalov_sym]
print("Viral targets from Stukalov at al in GWAS subnetwork")
stukalov_viral_targets_all

husci_viral_targets_paralogs <- V(observation_paralogs)$name[V(observation_paralogs)$name %in% husci_sym]
print("Viral targets from HuSCI in GWAS subnetwork, without paralogs")
husci_viral_targets_paralogs

gordon_viral_targets_paralogs <- V(observation_paralogs)$name[V(observation_paralogs)$name %in% gordon_sym]
print("Viral targets from Gordon et al in GWAS subnetwork, without paralogs")
gordon_viral_targets_paralogs

stukalov_viral_targets_paralogs <- V(observation_paralogs)$name[V(observation_paralogs)$name %in% stukalov_sym]
print("Viral targets from Stukalov et al in GWAS subnetwork, without paralogs")
stukalov_viral_targets_paralogs

# b. interaction
interactions_all <- gsize(observation_all)
print("Interactions in GWAS subnetwork")
interactions_all

interactions_paralogs <- gsize(observation_paralogs)
print("Interactions in GWAS subnetwork, without paralogs")
interactions_paralogs

# c. average shortest path between GWAS proteins
dist_all <- distances(huri_g, gwas_huri, to = gwas_huri)
gwas_protein_shortest_path_all <- mean(dist_all[lower.tri(dist_all)])
print("Average shortest path between GWAS proteins")
gwas_protein_shortest_path_all

dist_paralogs <- distances(huri_g, gwas_huri_paralogs, to = gwas_huri_paralogs)
gwas_protein_shortest_path_paralogs <- mean(dist_paralogs[lower.tri(dist_paralogs)])
print("Average shortest path between GWAS proteins, without paralogs")
gwas_protein_shortest_path_paralogs

######
# 3. permutation analysis
# all GWAS proteins
permutation_all <- c()
permutation_all <- c(permutation_all, mcreplicate(10000, rewire3Dataset(huri_g, gwas_huri), mc.cores = detectCores()))
permutation_all[is.na(permutation_all)] <- 0
permutation_all_df <- data.frame(matrix(permutation_all, ncol = 5, byrow = T))
names(permutation_all_df) <- c("HuSCI_viral_target", "Gordon_viral_target", "Stukalov_viral_target", "interactions", "GWAS_average_shortest_path")

# without paralogs
permutation_paralogs <- c()
permutation_paralogs <- c(permutation_paralogs, mcreplicate(10000, rewire3Dataset(huri_g, gwas_huri_paralogs), mc.cores = detectCores()))
permutation_paralogs[is.na(permutation_paralogs)] <- 0
permutation_paralogs_df <- data.frame(matrix(permutation_paralogs, ncol = 5, byrow = T))
names(permutation_paralogs_df) <- c("HuSCI_viral_target", "Gordon_viral_target", "Stukalov_viral_target", "interactions", "GWAS_average_shortest_path")

######
# plot function
toPlot <- function(value, viral_husci, viral_gordon, viral_stukalov, interaction, distance) {
    # HuSCI viral target in GWAS subnetwork
    plotHist(value[, 1], "HuSCI", length(viral_husci), 25, 0.03, 0.05)
    # Gordon viral target in GWAS subnetwork
    plotHist(value[, 2], "Gordon et al", length(viral_gordon), 20, 0.03, 0.05)
    # Stukalov viral target in GWAS subnetwork
    plotHist(value[, 3], "Stukalov et al", length(viral_stukalov), 20, 0.03, 0.05)
    # interaction
    dens_gwas <- hist(value[, 4], breaks = 20, plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, xlim = c(200, 1200), las = 1, yaxt = "n", xlab = "Number of interactions", main = "", cex.sub = 0.5)
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from HuRI")
    axis(side = 2, at = seq(0, 1200, by = 200), labels = seq(0, 0.12, by = 0.02), las = 1)
    arrows(interaction, 200, interaction, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value[, 4]) + 100, max(dens_gwas$counts), paste0("median = ", median(value[, 4])), col = "grey", cex = 0.5)
    text(interaction - 200, 350, paste0("observed = ", interaction, "\np < 0.0001") , cex = 0.4, pos = 4)

    # mean distance
    dens_gwas <- hist(value[, 5], breaks = 8, plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xlim = c(1.5, 3.5), yaxt = "n", xlab = "Average shortest path", main = "", cex.sub = 0.5)
    axis(side = 2, at = seq(0, 5000, by = 500), labels = seq(0, 0.5, by = 0.05), las = 1)
    arrows(distance, 600, distance, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value[, 5]) + 0.1, max(dens_gwas$counts), paste0("median = ", round(median(value[, 5]), 2)), col = "grey", cex = 0.5)
    text(round(distance, 2), 1000, paste0("observed = ", round(distance, 2), "\np = ",
        round(
            as.numeric(table(value[value[, 5] != Inf, 5] >= distance)["FALSE"]) /
            as.numeric(table(value[, 5] != Inf)["TRUE"]),
        3)), cex = 0.4, pos = 4)
}

######
# plot
# all
pdf(file = "Nature2021a_3dataset_HuRI.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
toPlot(permutation_all_df, husci_viral_targets_all, gordon_viral_targets_all, stukalov_viral_targets_all, interactions_all, gwas_protein_shortest_path_all)
dev.off()
# without paralogs
pdf(file = "Nature2021a_3dataset_HuRI_paralogs.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
toPlot(permutation_paralogs_df, husci_viral_targets_paralogs, gordon_viral_targets_paralogs, stukalov_viral_targets_paralogs, interactions_paralogs, gwas_protein_shortest_path_paralogs)
dev.off()

######
# display degree of viral targets in HuRI
#------------------------------------------------
# NOT run
#------------------------------------------------
# husci_deg <- data.frame(degree(huri_g, v = husci_huri))
# names(husci_deg) <- "degree"
# gordon_deg <- data.frame(degree(huri_g, v = gordon_huri))
# names(gordon_deg) <- "degree"
# stukalov_deg <- data.frame(degree(huri_g, v = stukalov_huri))
# names(stukalov_deg) <- "degree"

# hist(husci_deg$degree, xlab = "degree", main = "Degree of HuSCI proteins in HuRI")
# addtable2plot(100, 50, summary(husci_deg), vlines = TRUE, bty = "l", cex = 2)

# hist(gordon_deg$degree, xlab = "degree", main = "Degree of Gordon proteins in HuRI")
# addtable2plot(150, 50, summary(gordon_deg), vlines = TRUE, bty = "l", cex = 2)

# hist(stukalov_deg$degree, xlab = "degree", main = "Degree of Stukalov proteins in HuRI")
# addtable2plot(70, 100, summary(stukalov_deg), vlines = TRUE, bty = "l", cex = 2)

# wb <- createWorkbook()
# addWorksheet(wb, "HuSCI")
# writeData(wb, "HuSCI", husci_deg, rowNames = TRUE)

# addWorksheet(wb, "Gordon et al")
# writeData(wb, "Gordon et al", gordon_deg, rowNames = TRUE)

# addWorksheet(wb, "Stukalov et al")
# writeData(wb, "Stukalov et al", stukalov_deg, rowNames = TRUE)

# saveWorkbook(wb, "~/Documents/INET-work/virus_network/statistic_results/GWAS/Nature2021a_3dataset_HuRI.xlsx", overwrite = TRUE)
#------------------------------------------------

######
# save workarea data
save.image("Nature2021a_3dataset_HuRI.RData")
