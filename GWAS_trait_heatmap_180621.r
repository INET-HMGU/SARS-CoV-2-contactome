# Heatmap for GWAS hit traits and community correlation

library(openxlsx)
library(pheatmap)
library(circlize)
library(seriation)
library(dendextend)

file <- "../data/trait_matrix.xlsx"
data <- read.xlsx(file, sheet = "BETA_SE")
data_fdr <- read.xlsx(file, sheet = "FDR")
data_p <- read.xlsx(file, sheet = "pvalue")

plot_data <- data[, c(2:22)]
plot_fdr <- data_fdr[, c(2:22)]
plot_p <- data_p[, c(2:22)]

rownames(plot_fdr) <- data_fdr[, 1]
rownames(plot_data) <- data[, 1]
rownames(plot_p) <- data_p[, 1]

paletteLength <- 50
myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)
myBreaks <- c(seq(min(plot_data), 0, length.out = ceiling(paletteLength/2) + 1), seq(max(plot_data)/paletteLength, max(plot_data), length.out = floor(paletteLength/2)))
col_fun <- colorRamp2(c(-2, 0, 4), c("blue", "white", "red"))
col_fun_p <- colorRamp2(c(0, 1), c("white", "red"))
col_fun_p10 <- colorRamp2(c(0, 6), c("white", "red"))

plot_data_2 <- plot_data
plot_data_2[abs(plot_data_2) < 1]  <- 0
plot_p10 <- -log(plot_p, 10)
plot_fdr10 <- -log(plot_fdr, 10)
row_dim <- c(3, 5:13, 16:22, 26:29)
col_dim <- c(1:14, 16)
plot_p_binary0.1 <- ifelse(plot_fdr < 0.1, 1, 0)

row_dim2 <- c(16, 9, 11, 4, 10,  2, 5, 7, 12, 8, 13, 14, 6, 3, 1)
non_community_list <- c(891, 3205, 8, 3564, 364, 145, 1515, 3035, 1627)

gwas_comm <- read.xlsx("../data/GWAS_list.xlsx", sheet = "community")
gwas_trait <- read.xlsx("../data/GWAS_list.xlsx", sheet = "trait")
comm_position <- gwas_comm$index[!gwas_comm$community %in% non_community_list]

######
# plotting V1
phtmap <- pheatmap::pheatmap(t(plot_p10[comm_position, row_dim2]))
col_dend <- phtmap[[2]]
# reorder columns
col_dend <- rotate(col_dend, order = c("3682", "2769", "731", "316", "4227", "1353", "4160", "2831", "1833", "2545", "692", "926", "2398", "525", "1882", "1652", "3941", "415", "377", "571"))
pdf("../result/graph/gwas_10.pdf")
pheatmap::pheatmap(t(plot_p10[comm_position, row_dim2]),
    cluster_cols = as.hclust(col_dend), cluster_rows = F,
    color = colorRampPalette(c("white", "red"))(50), border_color = "white",
    display_numbers = ifelse(t(plot_fdr[comm_position, row_dim2]) < 0.05, "*", ""),
    fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(0, 1, 2, 3, 4, 5, 5.85),
    legend_labels = c("0", "1", "2", "3", "4", "5", "-log10(p)"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05) \n hierarchical cluster based on -log10(p)")
pheatmap::pheatmap(t(plot_data[comm_position, row_dim2]),
    cluster_cols = as.hclust(col_dend), cluster_rows = F,
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, border_color = "white",
    display_numbers = ifelse(t(plot_fdr[comm_position, row_dim2]) < 0.05, "*", ""),
    fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05) \n hierarchical cluster based on -log10(p)")
dev.off()
