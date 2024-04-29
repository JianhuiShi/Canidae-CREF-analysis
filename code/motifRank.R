###############################################################################
# Function:
# Compute the rank of all motifs according to their loadings
# in polarized motif-eigenvectors.
# Input:
# - "data/motif_symbol.txt"
# - "results/loadings/{species}_{motif_*}_loadings.txt"
# Output (optional):
# - "results/[motifRank/]motifRank[.RData/_level_{l}.csv/.xlsx]"
###############################################################################

# install.packages("rlist")
# install.packages("tidyverse")
# install.packages("openxlsx")

options(tidyverse.quiet = TRUE)

library(rlist)
library(tidyverse)
library(openxlsx)

# Configuration file ----------------------------------------------------------

setwd(".")
config <- list.load("./config.yaml")

## Directory setting ----------------------------------------------------------

data_dir <- "../data"

result_dir <- "../results"

motif.path <- file.path(data_dir, "motif_symbol.txt")

loading.dir <- file.path(result_dir, "loadings")

# motifRank.dir <- file.path(result_dir, "motifRank")
motifRank.dir <- result_dir
if (!dir.exists(motifRank.dir)) dir.create(motifRank.dir, recursive = TRUE)

# Parameters setting ----------------------------------------------------------

species_list <- config$species_list
speciesNum <- length(species_list)

from <- 1L
to <- 9L
level_list <- str_c("level_", from:to)

# read data --------------------------------------------------------------------

loading.data <- list()
for (s in species_list) {
    for (l in level_list) {
        m <- str_replace(l, "level", "motif")
        loading.data[[s]][[l]] <-
            read_tsv(
                file.path(loading.dir, str_glue("{s}_{m}_loadings.txt")),
                col_types = "cd"
            )
    }
}
loading.data <- loading.data %>% transpose()

motif_list <- read_table(
    motif.path,
    col_names = FALSE,
    col_types = "c"
) %>%
    pull()

# compute motif rank ----------------------------------------------------------

# motif_rankTable.list <- list()
# for (l in level_list) {
#     motif_rankTable <- tibble(Motif = motif_list, .rows = 1403L)
#     for (s in species_list) {
#         motifLoadings <- loading.data[[l]][[s]] %>%
#             # arrange(desc(Loading)) %>%
#             # mutate(Symbol = str_sub(Symbol, start = 3L)) %>%
#             pull(Symbol)
#         motif_rankTable <- motif_rankTable %>%
#             add_column(!!s := match(motif_list, motifLoadings))
#     }
#     motif_rankTable <- motif_rankTable %>%
#         modify_at(
#             .at = -1,
#             ~ if_else(. > 702L, . - 1404L, .)
#         )
#     motif_rankTable.list[[l]] <- motif_rankTable
# }

motif_rankTable.list <- list()
for (l in level_list) {
    motif_rankTable.list[[l]] <- loading.data[[l]] %>%
        # arrange(desc(Loading)) %>%
        # map( ~ mutate(., Symbol = str_sub(Symbol, start = 3L))) %>%
        map(~ pull(., Symbol)) %>%
        # map( ~ match(motif_list, .)) %>%
        map(~ factor(., levels = motif_list)) %>%
        map(order) %>%
        map_dfc(~ if_else(. <= 702L, ., . - 1404L)) %>%
        add_column(Motif = motif_list, .before = 1)
}

# sort motif alphabetically
motif_rankTable.list <- motif_rankTable.list %>%
    map(~ arrange(., Motif))

# write data ------------------------------------------------------------------

# # save as RData file
# save(
#     motif_rankTable.list,
#     file = file.path(motifRank.dir, "motifRank.RData")
# )

# # write to separate csv files
# for (l in level_list) {
#     write_csv(
#         motif_rankTable.list[[l]],
#         file = file.path(
#             motifRank.dir,
#             str_glue("motifRank_{l}.csv")
#         )
#     )
# }

# write to one Excel file
wb <- createWorkbook()
for (l in level_list) {
    addWorksheet(
        wb = wb,
        sheetName = l
    )
    setColWidths(
        wb = wb,
        sheet = l,
        cols = seq_len(speciesNum + 1),
        widths = "auto"
    )
    writeData(
        wb = wb,
        sheet = l,
        x = motif_rankTable.list[[l]]
    )
}
saveWorkbook(
    wb = wb,
    file = file.path(motifRank.dir, "motifRank.xlsx"),
    overwrite = TRUE
)
