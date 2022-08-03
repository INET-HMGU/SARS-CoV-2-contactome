######
# Display protein tissue specificity profiles
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
library(RColorBrewer)
library(dplyr)
library(pastecs)

######
# load dataset
ext_table2 <- "data/extended_table/Extended_Table_2_PPIs.xlsx"
ext_table5 <- "data/extended_table/Extended_Table_5_infected_organ.xlsx"
## HPA
hpa <- read.xlsx(file.path(google, "data/HPA_20210409.xlsx")) # subtracted from https://www.proteinatlas.org/download/proteinatlas.tsv.zip
hpa_tissue <- hpa[, c(3, 1, 7, 8)]
names(hpa_tissue) <- c("ensemblID", "Gene", "RNA_tissue_specificity", "RNA_tissue_distribution")
hpa_dorward <- read.xlsx(file.path(google, ext_table5), sheet = "HPA")
## HuSCI
binary <- read.xlsx(file.path(paper, "04_Supplementary Information/Supplementary_Table_1.xlsx"), sheet = '1b - HuSCI', startRow = 4)
husci_tissue <- unique(hpa_tissue[hpa_tissue$ensemblID %in% binary$"Ensembl.gene.ID", ])
husci_dorward <- read.xlsx(file.path(google, ext_table5), sheet = "HuSCI")
## Gordon et al., 2020
gordon <- read.xlsx(file.path(google, ext_table2), sheet = 3)
gordon_tissue <- unique(hpa_tissue[hpa_tissue$ensemblID %in% gordon$Ensembl_uniprotIDmapping, ])
gordon_dorward <- read.xlsx(file.path(google, ext_table5), sheet = "Gordon")
## Stukalov et al., 2021
stukalov <- read.xlsx(file.path(google, ext_table2), sheet = 4)
stukalov_tissue <- unique(hpa_tissue[hpa_tissue$ensemblID %in% stukalov$ensemblID, ])
stukalov_dorward <- read.xlsx(file.path(google, ext_table5), sheet = "Stukalov")
## Li et al., 2020
li <- read.xlsx(file.path(google, ext_table2), sheet = 5)
li_tissue <- unique(hpa_tissue[hpa_tissue$ensemblID %in% li$ensemblID, ])
li_dorward <- read.xlsx(file.path(google, ext_table5), sheet = "Li")
## Nabeel et al., 2020
nabeel <- read.xlsx(file.path(google, ext_table2), sheet = 6)
nabeel_tissue <- unique(hpa_tissue[hpa_tissue$ensemblID %in% nabeel$Prey_ensembl,])
nabeel_dorward <- read.xlsx(file.path(google, ext_table5), sheet = "Nabeel")
## Laurent et al., 2020
bioid1 <- read.xlsx(file.path(google, ext_table2), sheet = 7)
laurent_tissue <- unique(hpa_tissue[hpa_tissue$ensemblID %in% bioid1$ensemblID, ])
laurent_dorward <- read.xlsx(file.path(google, ext_table5), sheet = "Laurent")
## St-Germain et al., 2020
bioid_st <- read.xlsx(file.path(google, ext_table2), sheet = 8)
st_germain_tissue <- unique(hpa_tissue[hpa_tissue$ensemblID %in% bioid_st$ensemblID, ])
st_germain_dorward <- read.xlsx(file.path(google, ext_table5), sheet = "St_Germain")
## loading Samavarchi et al., 2020
bioid_sama <- read.xlsx(file.path(google, ext_table2), sheet = 9)
samavarchi_tissue <- unique(hpa_tissue[hpa_tissue$ensemblID %in% bioid_sama$ensemblID, ])
samavarchi_dorward <- read.xlsx(file.path(google, ext_table5), sheet = "Samavarchi")

######
# proportion of protein tissue specificity
tissue_tables <- list(
    HPA = table(hpa_tissue$RNA_tissue_specificity),
    HuSCI = table(husci_tissue$RNA_tissue_specificity),
    Gordon = table(gordon_tissue$RNA_tissue_specificity),
    Stukalov = table(stukalov_tissue$RNA_tissue_specificity),
    Li = table(li_tissue$RNA_tissue_specificity),
    Nabeel = table(nabeel_tissue$RNA_tissue_specificity),
    Laurent = table(laurent_tissue$RNA_tissue_specificity),
    St_Germain = table(st_germain_tissue$RNA_tissue_specificity),
    Samavarchi = table(samavarchi_tissue$RNA_tissue_specificity)
)

tissue <- data.frame(
    interactome = names(tissue_tables)[c(2, 1, 3:9)],
    specific = c(
        sum(tissue_tables[["HuSCI"]][c(1, 4, 5)]),
        sum(tissue_tables[["HPA"]][c(1, 4, 5)]),
        sum(tissue_tables[["Gordon"]][c(1, 3:4)]),
        sum(tissue_tables[["Stukalov"]][c(1, 3:4)]),
        sum(tissue_tables[["Li"]][c(1, 3:4)]),
        sum(tissue_tables[["Nabeel"]][c(1, 3:4)]),
        sum(tissue_tables[["Laurent"]][c(1, 3:4)]),
        sum(tissue_tables[["St_Germain"]][c(1, 3:4)]),
        sum(tissue_tables[["Samavarchi"]][c(1, 3:4)])
        ),
    common = c(
        tissue_tables[["HuSCI"]][2],
        tissue_tables[["HPA"]][2],
        tissue_tables[["Gordon"]][2],
        tissue_tables[["Stukalov"]][2],
        tissue_tables[["Li"]][2],
        tissue_tables[["Nabeel"]][2],
        tissue_tables[["Laurent"]][2],
        tissue_tables[["St_Germain"]][2],
        tissue_tables[["Samavarchi"]][2]
    )
)
tissue$interactome <- factor(tissue$interactome, levels = c("HuSCI", "HPA", "Gordon", "Stukalov", "Li", "Nabeel", "Laurent", "St_Germain", "Samavarchi"))
tissue_melt <- melt(tissue)

table_dorward <- list(
    HPA = table(hpa_dorward$dorwie),
    HuSCI = table(husci_dorward$dorwie),
    Gordon = table(gordon_dorward$dorwie),
    Stukalov = table(stukalov_dorward$dorwie),
    Li = table(li_dorward$dorwie),
    Nabeel = table(nabeel_dorward$dorwie),
    Laurent = table(laurent_dorward$dorwie),
    St_Germain = table(st_germain_dorward$dorwie),
    Samavarchi = table(samavarchi_dorward$dorwie)
)
dorward <- as.data.frame(do.call(bind_rows, table_dorward))
dorward[is.na(dorward)] <- 0
row.names(dorward) <- names(table_dorward)
dorward_melt <- melt(t(dorward))
names(dorward_melt) <- c("Organ", "interactome", "value")
dorward_melt$interactome <- factor(dorward_melt$interactome, levels = c("HuSCI", "HPA", "Gordon", "Stukalov", "Li", "Nabeel", "Laurent", "St_Germain", "Samavarchi"))
######
# plotting
gp_style <- theme_bw() +
    theme(axis.title = element_text(color = "black", size = 14, face = "bold"),
        axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
colour <- c("#601B9E", "#7F7F7F", "#028CC1", "#028CC1", "#028CC1", "#028CC1", "#DB1517", "#DB1517", "#DB1517", "#9524E9", "#B3B3B3", "#3BA9ED", "#3BA9ED", "#3BA9ED", "#3BA9ED", "#F24141", "#F24141", "#F24141")
colour2 <- c("#670120", "#BB2934", "#E48267", "#FACDB5", "#F7F7F7", "#C1DDEC", "#69ACD1", "#2971B1", "#053061")
gp_percent <- ggplot(tissue_melt, aes(x = interactome, y = value)) +
    geom_bar(stat = "identity", position = "fill", fill = colour, color = "black") +
    scale_linetype_manual(values = c("common" = "dashed", "specific" = "blank")) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0%", "20%", "40%", "60%", "80%", "100%"), expand = expansion(mult = c(0, 0.1))) +
    labs(x = "", y = "Fraction") +
    gp_style
ggsave(gp_percent, file = file.path(google, "result/graph/tissue_specificity.pdf"))

gp_dorward <- ggplot(dorward_melt, aes(x = interactome, y = value)) +
    geom_bar(aes(fill = Organ), stat = "identity", position = "fill", color = "black") +
    scale_linetype_manual(values = c("common" = "dashed", "specific" = "blank")) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0%", "20%", "40%", "60%", "80%", "100%"), expand = expansion(mult = c(0, 0.1))) +
    scale_color_manual(colour) +
    labs(x = "", y = "Fraction") +
    gp_style
ggsave(gp_dorward, file = file.path(google, "result/graph/infected_organ.pdf"))
