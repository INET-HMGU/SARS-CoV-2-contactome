###### permutation analysis for RNA-binding proteins
# Lin Chung-wen

######
# environment
google <- "/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/"
paper <- "/Volumes/GoogleDrive/My Drive/Paper_VirHostome_CoV2"

######
# load packages
library(openxlsx)
library(rethinking)
source(file.path(google, "src/random_sig.r"))

######
# load reference
rbp_cRIC <- read.xlsx(file.path(google, "data/RBP_Kamel.xlsx"), sheet = 2) # https://www.biorxiv.org/content/biorxiv/early/2020/11/29/2020.11.25.398008/DC1/embed/media-1.xlsx
rbp_vRIC <- read.xlsx(file.path(google, "data/RBP_Kamel_vRIC.xlsx"), sheet = 1) # https://www.biorxiv.org/content/biorxiv/early/2020/11/29/2020.11.25.398008/DC5/embed/media-5.xlsx
calu3_HuH <- read.delim(file.path(google, "data/Calu3-HuH_RBPs.txt"), header = F, sep = "\t")
calu3_H1 <- calu3_HuH$V1
calu3_HuH_E6 <- read.delim(file.path(google, "data/Calu3-HuH-VeroE6_RBPs.txt"), header = F, sep = "\t")
calu3_H2 <- calu3_HuH_E6$V1
# load HuSCI ----
binary <- read.xlsx(file.path(paper, "04_Supplementary Information/Supplementary_Table_1.xlsx"), sheet = '1b - HuSCI', startRow = 4)
binary_human <- unique(binary[, "Host.protein_symbol"])
######
# load search space
search_space <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_1_search_space.xlsx"), sheet = 2)
total_sym <- unique(search_space$ensembl_gene_name)
######
# filter dataset
rbp_cRIC$adj.P.Val[is.na(rbp_cRIC$adj.P.Val)] <- 1 # necessary because NA will also be included
rbp_cRIC_sig <- rbp_cRIC[rbp_cRIC$adj.P.Val < 0.1, "Gene.name"][c(1:330)]
rbp_vRIC <- rbp_vRIC[rbp_vRIC$logFC > 0, ] # all
rbp_vRIC2 <- rbp_vRIC[rbp_vRIC$adj.P.Val < 0.1, ] # adjp < 0.1, fc > 0
rbp_c_inT <- unique(rbp_cRIC$Gene.name[c(1:773)])
rbp_c_sig_inT <- unique(rbp_cRIC_sig[rbp_cRIC_sig %in% total_sym])
rbp_v1_inT <- unique(rbp_vRIC$Gene.name[rbp_vRIC$Gene.name %in% total_sym])
rbp_v2_inT <- unique(rbp_vRIC2$Gene.name[rbp_vRIC2$Gene.name %in% total_sym])
calu3_H1_inT <- unique(calu3_H1[calu3_H1 %in% total_sym])
calu3_H2_inT <- unique(calu3_H2[calu3_H2 %in% total_sym])
######
# permutation analysis
## cRIC all
binary_rbp_cRIC <- table(binary_human %in% rbp_c_sig_inT)[2]
random_rbp_cRIC <- simula(10000, total_sym, binary_human, rbp_c_inT)
sign_rbpcRIC <- round(1 - (abs(sig(random_rbp_cRIC, binary_rbp_cRIC)) / 10000), 4)

## display protein numbers and items, cRIC
message("Number of unique RNA-binding proteins in cRIC: ", length(rbp_c_sig_inT)) # SARS-CoV-2 proteins excluded
message("Number of HuSCI identified proteins in cRIC: ", binary_rbp_cRIC)
message("HuSCI identified proteins in cRIC:", paste0(binary_human[binary_human %in% rbp_c_sig_inT], collapse = ", "))
message(rep("=", times = 60))

## cRIC significant
binary_rbp_cRIC2 <- table(binary_human %in% unique(rbp_cRIC_sig))[2]
random_rbp_cRIC2 <- simula(10000, total_sym, binary_human, rbp_c_sig_inT)
sign_rbp_cRIC2 <- round(1 - (abs(sig(random_rbp_cRIC2, binary_rbp_cRIC2)) / 10000), 4)

## vRIC all
binary_rbpvric <- table(binary_human %in% unique(rbp_vRIC$Gene.name))[2]
random_rbpvric <- simula(10000, total_sym, binary_human, rbp_v1_inT)
sign_rbpvric <- round(1 - (abs(sig(random_rbpvric, binary_rbpvric)) / 10000), 4)

## vRIC significant
binary_rbpvric2 <- table(binary_human %in% unique(rbp_vRIC2$Gene.name))[2]
random_rbpvric2 <- simula(10000, total_sym, binary_human, rbp_v2_inT)
sign_rbpvric2 <- round(1 - (abs(sig(random_rbpvric2, binary_rbpvric2)) / 10000), 4)

# Calu3_HuH
binary_calu3_H1 <- table(binary_human %in% unique(calu3_H1_inT))[2]
random_calu3_H1 <- simula(10000, total_sym, binary_human, calu3_H1_inT)
sign_calu3_H1 <- round(1 - (abs(sig(random_calu3_H1, binary_calu3_H1)) / 10000), 4)

# Calu3_HuH_E6
binary_calu3_H2 <- table(binary_human %in% unique(calu3_H2_inT))[2]
random_calu3_H2 <- simula(10000, total_sym, binary_human, calu3_H2_inT)
sign_calu3_H2 <- round(1 - (abs(sig(random_calu3_H2, binary_calu3_H2)) / 10000), 4)

######
# plotting
xlab <- "Number of identified RBPs"
ylab <- "Fraction of\n random networks"
## 1. cRIC
pdf(file.path(google, "result/graph/random_RBP_cRIC.pdf"), width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
barplot(prop.table(table(random_rbp_cRIC)), xlim = c(0, 20), ylim = c(0, 0.2), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", space = 0, main = "from cRIC\nall")
mtext(side = 1, line = 2, cex = 1, xlab)
mtext(side = 2, line = 2, cex = 1, ylab)
axis(side = 1, at = seq(1, 20, 4) - 0.5, labels = seq(0, 19, 4), line = 0.3)
arrows(binary_rbp_cRIC + 0.5, 270/10000, binary_rbp_cRIC + 0.5, 10/10000, col = "#922687", lwd = 2, length = 0.1)
text(binary_rbp_cRIC - 1, 450/10000, paste0("observed = ", binary_rbp_cRIC, "\np = ", sign_rbpcRIC), cex = 0.7, pos = 4)
dev.off()

## 2. cRIC significant
pdf(file.path(google, "result/graph/random_RBP_cRIC_p0.1.pdf"), width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
barplot(prop.table(table(random_rbp_cRIC2)), xlim = c(0, 20), ylim = c(0, 0.25), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", space = 0, main = "from cRIC\n(adj.p < 0.1)")
mtext(side = 1, line = 2, cex = 1, xlab)
mtext(side = 2, line = 2, cex = 1, ylab)
axis(side = 1, at = seq(1, 20, 4) - 0.5, labels = seq(0, 19, 4), line = 0.3)
arrows(binary_rbp_cRIC2 + 0.5, 250/10000, binary_rbp_cRIC2 + 0.5, 10/10000, col = "#922687", lwd = 2, length = 0.1)
text(binary_rbp_cRIC2 - 1, 450/10000, paste0("observed = ", binary_rbp_cRIC2, "\np = ", sign_rbp_cRIC2), cex = 0.7, pos = 4)
dev.off()

## 3. vRIC
pdf(file.path(google, "result/graph/random_RBP_vRICall.pdf"), width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
barplot(prop.table(table(random_rbpvric)), xlim = c(0, 15), ylim = c(0, 0.25), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", space = 0, main = "from vRIC RBP\n(all, logFC > 0)")
axis(side = 1, at = seq(1, 15, 3) - 0.5, labels = seq(0, 14, 3), line = 0.3)
mtext(side = 1, line = 2, cex = 1, xlab)
mtext(side = 2, line = 2, cex = 1, ylab)
arrows(binary_rbpvric + 0.5, 260/10000, binary_rbpvric + 0.5, 10/10000, col = "#922687", lwd = 2, length = 0.1)
text(binary_rbpvric - 1, 430/10000, paste0("observed = ", binary_rbpvric, "\np = ", sign_rbpvric), cex = 0.7, pos = 4)
dev.off()

## 4. vRIC significant
pdf(file.path(google, "result/graph/random_RBP_vRIC_p0.1.pdf"), width = 3, height = 3)
par(mgp = c(0.1, 0.7, 0), ps = 8)
barplot(prop.table(table(random_rbpvric2)), xlim = c(0, 10), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", space = 0, main = "from vRIC RBP\n(adj.p < 0.1, logFC > 0)")
axis(side = 1, at = seq(1, 10, 2) - 0.5, labels = seq(0, 9, 2), line = 0.3)
mtext(side = 1, line = 2, cex = 1, xlab)
mtext(side = 2, line = 2, cex = 1, ylab)
arrows(binary_rbpvric2 + 0.5, 500/10000, binary_rbpvric2 + 0.5, 10/10000, col = "#922687", lwd = 2, length = 0.1)
text(binary_rbpvric2 - 0.5, 750/10000, paste0("observed = ", binary_rbpvric2, "\np = ", sign_rbpvric2), cex = 0.7, pos = 4)
dev.off()

## 5. Calu3_HuH significant
pdf(file.path(google, "result/graph/random_RBP_calu3HuH.pdf"), width = 3, height = 3)
par(mgp = c(0.1, 0.7, 0), ps = 8)
barplot(prop.table(table(random_calu3_H1)), xlim = c(0, 15), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", space = 0, main = "from Calu3 HuH")
axis(side = 1, at = seq(1, 15, 3) - 0.5, labels = seq(0, 14, 3), line = 0.3)
mtext(side = 1, line = 2, cex = 1, xlab)
mtext(side = 2, line = 2, cex = 1, ylab)
arrows(binary_calu3_H1 + 0.5, 500/10000, binary_calu3_H1 + 0.5, 10/10000, col = "#922687", lwd = 2, length = 0.1)
text(binary_calu3_H1 - 0.5, 750/10000, paste0("observed = ", binary_calu3_H1, "\np = ", sign_calu3_H1), cex = 0.7, pos = 4)
dev.off()
## 6. Calu3 HuH VeroE6 significant
pdf(file.path(google, "result/graph/random_RBP_calu3HuH_VeroE6.pdf"), width = 3, height = 3)
par(mgp = c(0.1, 0.7, 0), ps = 8)
barplot(prop.table(table(random_calu3_H2)), xlim = c(0, 15), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", space = 0, main = "from Calu3 HuH VeroE6")
axis(side = 1, at = seq(1, 15, 3) - 0.5, labels = seq(0, 14, 3), line = 0.3)
mtext(side = 1, line = 2, cex = 1, xlab)
mtext(side = 2, line = 2, cex = 1, ylab)
arrows(binary_calu3_H2 + 0.5, 500/10000, binary_calu3_H2 + 0.5, 10/10000, col = "#922687", lwd = 2, length = 0.1)
text(binary_calu3_H2 - 0.5, 750/10000, paste0("observed = ", binary_calu3_H2, "\np = ", sign_calu3_H2), cex = 0.7, pos = 4)
dev.off()

# display protein numbers and items ----
## cRIC sign
cat(paste0("Number of unique RNA-binding proteins in cRIC (adjusted P < 0.1): ", length(unique(rbp_cRIC_sig))), "\n") # SARS-CoV-2 proteins excluded
cat(paste0("Number of HuSCI identified proteins in cRIC (adjP < 0.1): ", binary_rbp_cRIC2), "\n")
cat(paste(c("HuSCI identified proteins in cRIC:", binary_human[binary_human %in% unique(rbp_cRIC_sig)]), collapse = " "), "\n")
cat(paste(rep("=", times = 60), collapse = ""), "\n")
## vRIC
cat(paste0("Number of unique RNA-binding proteins in vRIC (logFC > 1): ", length(unique(rbp_vRIC$Gene.name))), "\n") # SARS-CoV-2 proteins excluded
cat(paste0("Number of HuSCI identified proteins in vRIC (logFC > 1): ", binary_rbpvric), "\n")
cat(paste(c("HuSCI identified proteins in vRIC:", binary_human[binary_human %in% unique(rbp_vRIC$Gene.name)]), collapse = " "), "\n")
cat(paste(rep("=", times = 60), collapse = ""), "\n")
## vRIC sign
cat(paste0("Number of unique RNA-binding proteins in vRIC (logFC > 1, adjP < 0.1): ", length(unique(rbp_vRIC2$Gene.name))), "\n") # SARS-CoV-2 proteins excluded
cat(paste0("Number of HuSCI identified proteins in vRIC (logFC > 1): ", binary_rbpvric2), "\n")
cat(paste(c("HuSCI identified proteins in vRIC:", binary_human[binary_human %in% unique(rbp_vRIC2$Gene.name)]), collapse = " "), "\n")
cat(paste(rep("=", times = 60), collapse = ""), "\n")

## Calu3 HuH7 sign
cat(paste0("Number of unique RNA-binding proteins in Calu3 HuH: ", length(unique(calu3_H1))), "\n") # SARS-CoV-2 proteins excluded
cat(paste0("Number of HuSCI identified proteins in Calu3 HuH: ", binary_calu3_H1), "\n")
cat(paste(c("HuSCI identified proteins in Calu3 HuH:", binary_human[binary_human %in% unique(calu3_H1)]), collapse = " "), "\n")
cat(paste(rep("=", times = 60), collapse = ""), "\n")
## Calu3 HuH7 VeroE6 sign
cat(paste0("Number of unique RNA-binding proteins in Calu3 HuH Vero E6: ", length(unique(calu3_H2))), "\n") # SARS-CoV-2 proteins excluded
cat(paste0("Number of HuSCI identified proteins in Calu3 HuH Vero E6: ", binary_calu3_H2), "\n")
cat(paste(c("HuSCI identified proteins in Calu3 HuH:", binary_human[binary_human %in% unique(calu3_H2)]), collapse = " "), "\n")
cat(paste(rep("=", times = 60), collapse = ""), "\n")
