# Intraviral
# Lin Chung-wen

######
# environment
google <- "/Volumes/GoogleDrive/My Drive/VirHostome_CW"
paper <- "/Volumes/GoogleDrive/My Drive/Paper_VirHostome_CoV2"

######
# load packages
library(igraph)
library(openxlsx)
library(rethinking)
source(file.path(google, "GitHub/src/random_sig.r"))

######
# load data
sars2 <- read.xlsx(file.path(paper, "04_Supplementary Information/Supplementary_Table_1.xlsx"), sheet = '1a - Viral ORFs', startRow = 4)
sars2$"Viral.protein"[c(11, 12)] <- "NSP12"
sars2 <- toupper(unique(sars2$"Viral.protein"[c(1:30)]))

y2h_intra <- read.xlsx(file.path(paper, "04_Supplementary Information/Supplementary_Table_1.xlsx"), sheet = '1c - IntraSCI', startRow = 4)
y2h_intra <- data.frame(vAD = toupper(y2h_intra$"Prey.(AD)"), vDB = toupper(y2h_intra$"Bait.(DB)"))
y2h_graph <- graph_from_data_frame(y2h_intra, directed = FALSE)
y2h_graph <- simplify(y2h_graph, remove.loops = FALSE)

intra_viral <- read.csv(file.path(google, "cloning/data/Y2H_list.csv"), header = T)
intra_graph <- graph_from_data_frame(intra_viral, directed = FALSE)

binary_intra <- gsize(intersection(intra_graph, y2h_graph))

sars2_rep <- function(x) {
  df <- data.frame(from = sample(x, length(sars2), replace = TRUE), to = sample(x, length(sars2), replace = TRUE))
  sample_g <- simplify(graph_from_data_frame(df, directed = FALSE), remove.loops = FALSE)
  return(gsize(intersection(intra_graph, sample_g)))
}

random_intra <- mcreplicate(10000, sars2_rep(sars2), mc.cores = parallel::detectCores())

sign_intra <- round(1 - (abs(sig(as.numeric(random_intra), as.numeric(binary_intra))) / 10000), 4)

message(rep("#", 60), "\n", "Total intraSCI PPIs: ", gsize(y2h_graph), "\n", rep("#", 60))
message("\n", rep("#", 60), "\n", "Observed in intra_viral identified human target within VirHostome:", binary_intra, "\np = ", sign_intra, "\n", rep("#", 60))

######
# plotting
dens <- hist(random_intra, breaks = 6, right = FALSE, plot = FALSE)
ymax <- round(max(dens$count)/1000, 1) * 1000
pdf(file.path(google, "GitHub/result/graph/random_intraviral.pdf"), width = 3, height = 3)
par(ps = 8, mgp = c(0, 0.6, -0.1))
plot(dens, xlim = c(0, 8), col = rgb(0.75, 0.75, 0.75, 1/2), freq = TRUE, border = NA, las = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "IntraSCI PPI")
axis(side = 1, at = seq(0, 7, by = 1) + 0.5, labels = seq(0, 7, by = 1), las = 1)
axis(side = 2, at = seq(0, ymax, by = 500), labels = seq(0, ymax/10000, by = 0.05), las = 1)
mtext(side = 1, line = 1.5, "Number of PPIs in Li, et al.")
mtext(side = 2, line = 2, "Fraction of\nrandom networks")
arrows(binary_intra + 0.5, 500, binary_intra + 0.5, 0, col = "#922687", lwd = 2, length = 0.1)
text(binary_intra - 1, 1000, label = paste("observed = ", binary_intra, "\np =", sign_intra), cex = 0.7, pos = 4)
dev.off()
