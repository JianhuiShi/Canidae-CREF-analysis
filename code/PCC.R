# Copyright © 2024, Jianhui Shi & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

###############################################################################
# Function:
# Compute the Pearson correlation coefficients between motif-eigenvectors
# across speices.
# Input:
# - "results/SVD/{species}_V0_1402.txt"
# Output (optional):
# - "results/[PCC/]cc.[RData/tsv/xlsx]"
###############################################################################

# install.packages("rlist")
# install.packages("tidyverse")

options(tidyverse.quiet = TRUE)

library(tidyverse)

# Configuration file ----------------------------------------------------------

setwd(".")
config <- rlist::list.load("./config.yaml")

## Directory setting ----------------------------------------------------------

data_dir <- "../data"

result_dir <- "../results"

SVD.dir <- file.path(result_dir, "SVD")

# cc.dir <- file.path(result_dir, "PCC")
cc.dir <- result_dir
if (!dir.exists(cc.dir)) dir.create(cc.dir, recursive = TRUE)

## Parameters setting ---------------------------------------------------------

species_list <- config$species_list
speciesNum <- length(species_list)

level_list <- 1L:10L

# read data -------------------------------------------------------------------

motif.v.list <- list()
for (s in species_list) {
    # cat(s, "\n")
    motif.v.list[[s]] <- read_tsv(
        file.path(SVD.dir, str_glue("{s}_V0_1402.txt")),
        col_select = num_range(prefix = "", range = level_list),
        col_names = as.character(level_list),
        col_types = str_dup("d", length(level_list))
    ) %>%
        as.list() %>%
        enframe(name = "level", value = "vec") %>%
        add_column(species = s, .before = 2)
}

motif.v.compact <- bind_rows(motif.v.list) %>%
    mutate(species = factor(species, levels = species_list)) %>%
    mutate(level = as.integer(level)) %>%
    filter(level %in% 2L:10L) %>%
    arrange(level, species)

# compute ---------------------------------------------------------------------

cc.tibble <- expand_grid(level = 2L:10L,
                         species_x = species_list,
                         species_y = species_list) %>%
    filter(species_x != species_y) %>%
    left_join(
        motif.v.compact,
        by = c("level" = "level", "species_x" = "species")
    ) %>%
    left_join(
        motif.v.compact,
        by = c("level" = "level", "species_y" = "species"),
        suffix = c("_x", "_y")
    ) %>%
    rowwise() %>%
    mutate(cc = cor(vec_x, vec_y)) %>%
    # mutate(cc = abs(cor(vec_x, vec_y))) %>%
    ungroup() %>%
    select(!c(vec_x, vec_y)) %>%
    mutate(level = level - 1L) %>%
    pivot_wider(
        names_from = level,
        values_from = cc
    )

# write data ------------------------------------------------------------------

# save(
#     cc.tibble,
#     file = file.path(cc.dir, "PCC.RData")
# )

# write_tsv(
#     x = cc.tibble,
#     file = file.path(cc.dir, "PCC.tsv")
# )

openxlsx::write.xlsx(
    x = cc.tibble %>%
        modify_at(.at = as.character(level_list), round, digits = 3),
    file = file.path(cc.dir, "PCC.xlsx"),
    overwrite = TRUE,
    colWidths = "auto"
)
