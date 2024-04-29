###############################################################################
# Usage:
# Rscript format_enrichment_result.R
# Function:
# formats the raw enrichment results (.RData) into readable file (.xlsx).
# Input:
# - "results/enrichmentAnalysis/RData/enrichment_result_{species}.RData"
# Output:
# - "results/enrichmentAnalysis/xlsx/{species}.xlsx"
###############################################################################

# install.packages("rlist")
# install.packages("tidyverse")
# install.packages("openxlsx")

options(tidyverse.quiet = TRUE)

library(rlist)
library(tidyverse)
library(openxlsx)

# Configuration file ----------------------------------------------------------

setwd("./enrichmentAnalysis")
config <- list.load("../config.yaml")

# Directory setting -----------------------------------------------------------

data_dir <- "../../data"

result_dir <- "../../results"

enrichment.dir <- file.path(result_dir, "enrichmentAnalysis")

raw_result.dir <- file.path(enrichment.dir, "RData")

result.dir <- file.path(enrichment.dir, "xlsx")
if (!dir.exists(result.dir)) dir.create(result.dir)

# Parameters setting ----------------------------------------------------------

# 在富集的原始计算结果中，筛选出通路和基因特征向量里公共基因个数在 [minCommonGeneNum, maxCommonGeneNum] 范围内的通路输出到 .xlsx 文件
minCommonGeneNum <- 1L
maxCommonGeneNum <- Inf

# 截取至多 trimGeneNum 个基因的名字存入字符串，防止字符串太长导致存入Excel单元格时溢出报错
trimGeneNum <- 50L

# signif_level <- 0.05
signif_level <- 0.10

from <- 1L
to <- 9L

species_list <- config$species_list
species_family_list <- config$species_info_list %>%
    list.filter(species %in% species_list) %>%
    list.map(list(species = species, family = family))
# family_list <- list.mapv(species_family_list, family) %>% unique()
family_list <- list.cases(species_family_list, family, sorted = FALSE)

gene_list <- str_c("gene_", from:to)
pathway_list <- c(
    "go.bp",
    "go.cc",
    "go.mf",
    "kegg",
    "reactome"
)

# load functions --------------------------------------------------------------

source("./enrichment_utils.R")

# read data -------------------------------------------------------------------

enrichment_result <- read_raw_result(raw_result.dir)

# format result and save ------------------------------------------------------

for (s in species_list) {
    for (g in gene_list) {
        for (p in pathway_list) {
            enrichment_result[[s]][[g]][[p]] <- enrichment_result[[s]][[g]][[p]] %>%
                filterGeneNum(minCommonGeneNum, maxCommonGeneNum) %>%
                sortGene() %>%
                trimLongGeneString(trimGeneNum)
        }
    }
}

outputCols <- c(
    "database",
    "id",
    str_c("pvalue.", rep(from:to, each = 2L), rep(c("+", "-"), times = to - from + 1L)),
    "description",
    "pathwayGeneNum",
    "num.of.targets",
    "num.of.non.targets",
    "pathwayGene_in_geneEigenvectors",
    "pathwayGene",
    NULL
)
bgFill <- list()
bgFill$pos <- RColorBrewer::brewer.pal(9, "Reds")[c(2, 4, 6)] %>% rev()
bgFill$neg <- RColorBrewer::brewer.pal(9, "Blues")[c(2, 4, 6)] %>% rev()
for (s in species_list) {
    result <- formatResult(enrichment_result, s)
    saveResult(result, s, result.dir)
}
