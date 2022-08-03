# COVID-19 GWAS hit analysis
# Lin Chung-win

######
# load packages
library(openxlsx)
library(igraph)
library(rethinking)

######
# load functions
source("combineNetwork.r")
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
# data
huri <- read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = 2)
gwas <- read.csv("../data/GWAS_hits.csv", header = T)
husci <- read.csv("../data/HuSCI_node.csv", header = T)

huri_symbol <- huri[, c(5:6)]
huri_g_ori <- graph_from_data_frame(huri_symbol, directed = FALSE)
huri_g <- simplify(huri_g_ori, remove.loops = FALSE)

husci_sym <- husci[husci$category == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole
gwas_huri <- gwas$name[!is.na(gwas$ctl == 1)] # GWAS hit in HuRI

gwas_hit_1st <- make_ego_graph(huri_g, nodes = gwas_huri, order = 1, mode = "all")

gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)

gwas_all_final <- simplify(induced_subgraph(huri_g, names(V(gwas_all_g_merge))), remove.loops = F) # V:169, E:1108
gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym] # 20 viral targets
gwas_all_husci_length <- length(gwas_all_husci)

gwas_rand_r2 <- c()
gwas_rand_r2 <- c(gwas_rand_r2, mcreplicate(10000, huriRewireHusci(gwas_huri, FALSE), mc.cores = detectCores()))
gwas_rand_df_r2 <- data.frame(matrix(gwas_rand_r2, ncol = 2, byrow = T))
names(gwas_rand_df_r2) <- c("viral_target", "interactions")
write.csv(gwas_rand_df_r2, file = "../result/10000_randomHuRI.csv", row.names = FALSE)

# plotting
pdf("../result/graph/random_GWAS_viral_target.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
dens_gwas <- hist(gwas_rand_df_r2[, "viral_target"], breaks = 15, plot = FALSE)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1 / 2), xlim = c(0, 25), freq = FALSE, border = NA, las = 1, xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes"
mtext(side = 3, line = 1, cex = 1, mytitle)
arrows(gwas_all_husci_length + 0.5, 0.03, gwas_all_husci_length + 0.5, 0.01, col = "#922687", lwd = 2, length = 0.1)
text(gwas_all_husci_length - 2, 0.05, paste0("observed = 20 \np = ", table(gwas_rand_df_r2[, "viral_target"] >= 20)["TRUE"] / 10000), cex = 0.4, pos = 4)

dens_gwas <- hist(gwas_rand_df_r2[, "interactions"], breaks = 25, plot = FALSE)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1 / 2), border = NA, las = 1, xlim = c(200, 1200), yaxt = "n", xlab = "Number of interactions", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes"
mtext(side = 3, line = 1, cex = 1, mytitle)
axis(side = 2, at = seq(0, 1400, by = 200), labels = seq(0, 0.14, by = 0.02), las = 2)
arrows(gsize(gwas_all_final) + 0.5, 200, gsize(gwas_all_final) + 0.5, 20, col = "#922687", lwd = 2, length = 0.1)
text(gsize(gwas_all_final) - 200, 400, paste0("observed = ", gsize(gwas_all_final), "\np < 0.0001"), cex = 0.4, pos = 4)

dev.off()
