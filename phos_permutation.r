# permutation analysis for phosphoregulated proteins
# Lin Chung-wen

######
# environment
google <- "/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/"
paper <- "/Volumes/GoogleDrive/My Drive/Paper_VirHostome_CoV2"

######
# load packages
library(rethinking)
library(openxlsx)
source(file.path(google, "src/random_sig.r"))

######
# load search space dataset
search_space <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_1_search_space.xlsx"), sheet = 2)
total_sym <- unique(search_space$ensembl_gene_name)
######
# load phospho protein dataset
phos_bouhaddou <- read.xlsx(file.path(google, "data/Phospho_Krogan.xlsx"), sheet = 1)
phos_stukalov <- read.xlsx(file.path(google, "data/Phospho_Pichlmair.xlsx"), sheet = 1)
phos <- unique(c(phos_stukalov$Gene.symbol, phos_bouhaddou$Gene_Name))
phos_inT <- phos[phos %in% total_sym]
######
# load HuSCI dataset
binary <- read.xlsx(file.path(paper, "04_Supplementary Information/Supplementary_Table_1.xlsx"), sheet = '1b - HuSCI', startRow = 4)
binary_human <- unique(binary[, "Host.protein_symbol"])
######
# permutation analysis
binary_phos <- table(binary_human %in% phos)[2]
random_phos <- simula(10000, total_sym, binary_human, phos_inT)
sign_phos <- round(1 - (abs(sig(as.numeric(random_phos), as.numeric(binary_phos))) / 10000), 4)
######
# plot generation
dens <- hist(random_phos, breaks = 40, right = TRUE, plot = FALSE)
pdf(file.path(google, "result/graph/random_phosphoprotein.pdf"), width = 3, height = 3)
par(mgp = c(0.1, 0.7, 0), ps = 8)
plot(dens, xlim = c(0, 50), col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlab = "", ylab = "", xaxt = "n", main = "Phosphoregulated protein")
axis(side = 1, at = seq(1, 50, 10) - 0.5, labels = seq(0, 49, 10), line = 0.3)
mtext(side = 1, line = 2.5, "Number of phosphoregulated host\ntargets")
mtext(side = 2, line = 2, "Fraction of\nrandom networks")
arrows(binary_phos + 0.5, 150/10000, binary_phos + 0.5, 10/10000, col = "#922687", lwd = 2, length = 0.1)
text(binary_phos - 5, 250/10000, paste0("observed = ", binary_phos, "\np = ", sign_phos), cex = 0.7, pos = 4)
dev.off()
