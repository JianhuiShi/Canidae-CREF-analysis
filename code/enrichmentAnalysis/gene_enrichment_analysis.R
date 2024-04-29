# Copyright © 2024, Jianhui Shi & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

###############################################################################
# Function:
# Perform the gene enrichment analysis by the Wilcoxon scoring method
# on the polarized gene-eigenvector.
# Input:
# - "data/pathway/{family}.{go.bp/go.cc/go.mf/reactome/kegg}.tsv"
# - "results/loadings/{species}_{gene_*}_loadings.txt"
# Output:
# - "results/enrichmentAnalysis/RData/enrichment_result_{species}.RData"
###############################################################################

# install.packages("rlist")
# install.packages("tidyverse")
# install.packages("snowfall")

options(tidyverse.quiet = TRUE)

library(rlist)
library(tidyverse)
library(snowfall)

# Configuration file ----------------------------------------------------------

setwd("./enrichmentAnalysis")
config <- list.load("../config.yaml")

# Directory setting -----------------------------------------------------------

data_dir <- "../../data"

result_dir <- "../../results"

database.dir <- file.path(data_dir, "pathway")

loading.dir <- file.path(result_dir, "loadings")

enrichment.dir <- file.path(result_dir, "enrichmentAnalysis")
if (!dir.exists(enrichment.dir)) dir.create(enrichment.dir, recursive = TRUE)

raw_result.dir <- file.path(enrichment.dir, "RData")
if (!dir.exists(raw_result.dir)) dir.create(raw_result.dir, recursive = TRUE)

# Parameters setting ----------------------------------------------------------

CPUs <- 120L

from <- 4L
to <- 4L

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

pathway_data.list <- list()
for (f in family_list) {
    for (p in pathway_list) {
        pathway_data.list[[f]][[p]] <- read_tsv(
            file.path(
                database.dir,
                str_glue("{f}.{p}.tsv")
            ),
            col_types = cols(
                id = col_character(),
                geneNum = col_integer(),
                .default = col_character()
            )
        )
    }
}

# pre-processing --------------------------------------------------------------

# 做富集前，筛选出基因个数在 [minGeneNum, maxGeneNum] 范围内的基因通路参与富集的计算
minGeneNum <- 1L
maxGeneNum <- Inf
for (f in family_list) {
    for (p in pathway_list) {
        pathway_data.list[[f]][[p]] <- pathway_data.list[[f]][[p]] %>%
            filter(geneNum %>% between(minGeneNum, maxGeneNum))
    }
}

# parallel computing ----------------------------------------------------------

outputCols <- c(
    "id",
    "pvalue.pos",
    "pvalue.neg",
    "description",
    "num.of.targets",
    "num.of.non.targets",
    "pathwayGeneNum",
    "pathwayGene_in_geneEigenvectors",
    "pathwayGene"
)

par.list <- expand_grid(
    species_family_list %>% bind_rows(),
    gene = gene_list
) %>%
    # transpose() %>%
    list.parse() %>%
    list.names(str_glue("{species} {gene}"))

time.start <- Sys.time()
sfInit(parallel = TRUE, cpus = CPUs)
sfLibrary(openxlsx)
sfLibrary(rlist)
sfLibrary(tidyverse)
sfExportAll()
result <- sfLapply(
    par.list,
    do_enrichment,
    loading.dir
)
sfStop()
time.end <- Sys.time()
print(time.end - time.start)

# transform result ------------------------------------------------------------

enrichment_result <- list()
for (s in species_list) {
    for (g in gene_list) {
        enrichment_result[[s]][[g]] <- result[[str_glue("{s} {g}")]]
    }
}

# write data ------------------------------------------------------------------

write_raw_result(enrichment_result, raw_result.dir)
