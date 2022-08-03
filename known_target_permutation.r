# previous identified viral targets
# Lin Chung-wen

######
# environment
google <- "/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/"
paper <- "/Volumes/GoogleDrive/My Drive/Paper_VirHostome_CoV2"

######
# load package
library(openxlsx)
library(rethinking)
source(file.path(google, "src/random_sig.r"))

######
# load search space dataset
search_space <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_1_search_space.xlsx"), sheet = 2)
total_sym <- unique(search_space$ensembl_gene_name)
######
# load known viral targets
previous <- read.xlsx(file.path(google, "data/known_targets.xlsx")) # IntAct
previous <- unique(previous$unique)
previous_inT <- previous[previous %in% total_sym]
# load HuSCI dataset
binary <- read.xlsx(file.path(paper, "04_Supplementary Information/Supplementary_Table_1.xlsx"), sheet = '1b - HuSCI', startRow = 4)
binary_node <- unique(binary[, "Host.protein_symbol"])
######
# permutation analysis
binary_pre <- table(binary_node %in% previous_inT)["TRUE"]

random_pre <- mcreplicate(10000, table(sample(total_sym, length(binary_node)) %in% previous_inT)[2], mc.cores = detectCores())
sign_pre <- round(1 - (abs(sig(as.numeric(random_pre), as.numeric(binary_pre))) / 10000), 4)

######
# message for checking
message(rep("#", 60), "\nHuSCI human targets: ", length(binary_node), "\n", rep("#", 60))
message(rep("#", 60), "\nPreviously known human targets: ", length(previous_inT), "\n", rep("#", 60))
message(rep("#", 60), "\nIntersection: ", binary_pre, "\n", rep("#", 60))
message(rep("#", 60), "\nPermutation significance: ", sign_pre, "\n", rep("#", 60))
######
# plotting
dens <- hist(random_pre, right = FALSE, plot = FALSE)
ymax <- round(max(dens$count)/1000, 1) * 1000
pdf(file.path(google, "result/graph/random_preIdentified.pdf"), width = 3, height = 3)
par(ps = 8)
plot(dens, xlim = c(0, 70), col = rgb(0.75, 0.75, 0.75, 1/2), freq = TRUE, border = NA, las = 1, xlab = "", ylab = "", yaxt = "n", main = "Known host targets")
axis(side = 2, at = seq(0, ymax, by = 500), labels = seq(0, ymax/10000, by = 0.05), las = 1)
mtext(side = 1, line = 2, "Number of previously\nknown host targets")
mtext(side = 2, line = 2, "Fraction of\nrandom networks")
arrows(binary_pre, 200, binary_pre, 0, col = "#922687", lwd = 2, length = 0.1)
text(binary_pre - 15, 400, label = paste("observed = ", binary_pre, "\np =", sign_pre), cex = 0.7, pos = 4)
dev.off()
