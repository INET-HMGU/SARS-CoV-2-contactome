######
# Display protein subcellular location profiles
## Lin, Chung-wen

######
# environment
google <- "/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/"
paper <- "/Volumes/GoogleDrive/My Drive/Paper_VirHostome_CoV2"

######
# load packages
library(openxlsx)
library(reshape2)
library(ggplot2)
library(dplyr)
library(pastecs)
library(tidyverse)
source(file.path(google, "src/subLocation.r"))

######
# load dataset
## HPA
hpa <- read.xlsx(file.path(google, "data/HPA_20210409.xlsx")) # subtracted from https://www.proteinatlas.org/download/proteinatlas.tsv.zip
hpa_subcell <- hpa[, c(3, 1, 23)]

sub_major <- read.xlsx(file.path(google, "data/Extended_table/Extended_Table_3_subcellular.xlsx"), sheet = "lookup_table")
sub_major$major <- tolower(sub_major$major)
sub_major$subcellular <- tolower(sub_major$subcellular)
## HuSCI
binary <- read.xlsx(file.path(paper, "04_Supplementary Information/Supplementary_Table_1.xlsx"), sheet = '1b - HuSCI', startRow = 4)
## Gordon et al., 2020
gordon <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_2_PPIs.xlsx"), sheet = 3)
## Stukalov et al., 2021
stukalov <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_2_PPIs.xlsx"), sheet = 4)
## Li et al., 2020
li <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_2_PPIs.xlsx"), sheet = 5)
## Nabeel et al., 2020
nabeel <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_2_PPIs.xlsx"), sheet = 6)
## Laurent et al., 2020
bioid1 <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_2_PPIs.xlsx"), sheet = 7)
## St-Germain et al., 2020
bioid_st <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_2_PPIs.xlsx"), sheet = 8)
## loading Samavarchi et al., 2020
bioid_sama <- read.xlsx(file.path(google, "data/extended_table/Extended_Table_2_PPIs.xlsx"), sheet = 9)

######
# count ratio of gene subcellular location
total <- c(
    length(unique(hpa_subcell$Ensembl)),
    length(unique(binary$"Ensembl.gene.ID")),
    length(unique(gordon$Ensembl_uniprotIDmapping)),
    length(unique(stukalov$ensemblID)),
    length(unique(li$ensemblID)),
    length(unique(nabeel$Prey_ensembl)),
    length(unique(bioid1$ensemblID)),
    length(unique(bioid_st$ensemblID)),
    length(unique(bioid_sama$ensemblID))
    )

allPPIs <- list(
    HPA = hpa_subcell,
    HuSCI = unique(merge(binary, hpa_subcell, by.x = "Ensembl.gene.ID", by.y = "Ensembl", all.x = TRUE)[, c(1, 4, 17)]),
    Gordon = unique(merge(gordon, hpa_subcell, by.x = "Ensembl_uniprotIDmapping", by.y = "Ensembl", all.x = TRUE)[, c(1, 5, 8)]),
    Stukalov = unique(merge(stukalov, hpa_subcell, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1, 3, 6)]),
    Li = unique(merge(li, hpa_subcell, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1, 3, 8)]),
    Nabeel = unique(merge(nabeel, hpa_subcell, by.x = "Prey_ensembl", by.y = "Ensembl", all.x = TRUE)[, c(1, 4, 7)]),
    Laurent = unique(merge(bioid1[, c(4, 2)], hpa_subcell, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1, 2, 4)]),
    St_Germain = unique(merge(bioid_st[, c(2, 3)], hpa_subcell, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1, 2, 4)]),
    Samavarchi = unique(merge(bioid_sama[, c(2, 3)], hpa_subcell, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1, 2, 4)])
        )

######
# summary of raw data
## percentage calculation
subset_major_p <- lapply(allPPIs, subLocaPtoM)
subset_major_pt <- subset_major_p %>% reduce(left_join, by = "Var1") %>% setNames(., c("subcell", names(subset_major_p)))
subset_major_pt[is.na(subset_major_pt)] <- 0
subset_m_melt <- melt(subset_major_pt)
subset_m_melt$variable <- factor(subset_m_melt$variable, levels = c("HuSCI", "Gordon", "Stukalov", "Li", "Nabeel", "Laurent", "St_Germain", "Samavarchi",  "HPA"))

subset_major_c <- lapply(allPPIs, subLocaCtoM)
subset_major_ct <- subset_major_c %>% reduce(left_join, by = "Var1") %>% setNames(., c("subcell", names(subset_major_c)))
subset_major_ct[is.na(subset_major_ct)] <- 0

######
# fisher test
f_pvalue <- list()
count <- 1
for (i in 1:14) {
    for (j in 3:10) {
        dt_tmp <- matrix(as.numeric(c(subset_major_ct[i, c(j, 2)], total[c(j - 1, 1)] - subset_major_ct[i, c(j, 2)])), ncol = 2)
        f_pvalue[[count]] <- fisher.test(dt_tmp)$p.value
        count <- count + 1
    }
}
f_pvalue_df <- data.frame(matrix(unlist(f_pvalue), ncol = 8, byrow = TRUE))
names(f_pvalue_df) <- names(subset_major_ct)[c(3:10)]
f_pvalue_adj <- p.adjust(unlist(f_pvalue), method = "fdr")
f_pvalue_adj_df <- data.frame(matrix(unlist(f_pvalue_adj), ncol = 8, byrow = TRUE))
names(f_pvalue_adj_df) <- names(subset_major_ct)[c(3:10)]

anno <- data.frame(subcell = rep(subset_major_ct$subcell, each = 8),
    xstar = ifelse(unlist(f_pvalue) < 0.05, rep(1:8, 14), NA),
    ystar = ifelse(unlist(f_pvalue) < 0.05, as.vector(t(subset_major_pt[, c(3:10)])), NA),
    lab = ifelse(f_pvalue_adj < 0.05, "*", NA))

######
# plot generate
colour <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#FF7F00", "#FFFF00", "brown", "grey")

ggplot(subset_m_melt, aes(x = variable, y = value * 100)) +
    geom_bar(aes(fill = variable), stat = "identity") +
    scale_fill_manual(values = colour) +
    geom_text(data = anno, aes(x = xstar, y = ystar * 100, label = lab)) +
    labs(x = "", y = "% of genes", title = "proportion of protein major subcellular location") +
    facet_wrap(~ subcell, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color = "black"),
        strip.text = element_text(size = 8),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"))
ggsave(file.path(google, "result/graph/subcellular_major_location_all.pdf"), width = 12, height = 10)
