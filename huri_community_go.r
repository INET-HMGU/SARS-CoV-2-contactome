# Gene Ontology analysis for HuRI communities

## load packages ----
library(gprofiler2)
library(dplyr)
library(rstatix)
library(stringr)

## custom functions ----
source("functionalEnrich.r")
source("exportGO.r")

## load community data ----
## and statistical results
load("../data/processed/HuRI_ocg.RData")
huri_all <- read.csv("../result/HuRI_signif_community.csv", header = TRUE)
check_list <- huri_all[huri_all[, 2] > 3 & huri_all[, 4] < 0.05, 1] # based on community size and fisher test p value

## load custom GMTs ----
bp_cust <- upload_GMT_file(gmtfile = "../data/processed/hsapien_HuRI_GOBP_EXP.gmt")
mf_cust <- upload_GMT_file(gmtfile = "../data/processed/hsapien_HuRI_GOMF_EXP.gmt")
cc_cust <- upload_GMT_file(gmtfile = "../data/processed/hsapien_HuRI_GOCC_EXP.gmt")

## GO enrichment analysis ----
huri_ocg_gobp <- funcEnrich(huri_ocg, bp_cust, check_list)
huri_ocg_gomf <- funcEnrich(huri_ocg, mf_cust, check_list)
huri_ocg_gocc <- funcEnrich(huri_ocg, cc_cust, check_list)

gobp <- export(huri_ocg_gobp, check_list)
gobp_df <- do.call(rbind.data.frame, gobp)
gomf <- export(huri_ocg_gomf, check_list)
gomf_df <- do.call(rbind.data.frame, gomf)
gocc <- export(huri_ocg_gocc, check_list)
gocc_df <- do.call(rbind.data.frame, gocc)

## save GO results ----
go_final <- rbind(gobp_df[, c(20, 2:13, 15:19)], gomf_df[, c(20, 2:13, 15:19)], gocc_df[, c(20, 2:13, 15:19)])

write.csv(go_final, file = "../result/HuRI_signif_community_GO.csv", row.names = FALSE)
