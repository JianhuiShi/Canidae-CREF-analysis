# Copyright © 2024, Jianhui Shi & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

###############################################################################
# Function:
# Compute the projections of the 4th and 5th motif-eigenvectors of one species
# onto the 2-D eigen space spanned by the 4th and 5th motif-eigenvectors of
# another species.
# Input:
# - "results/SVD/{species}_V0_1402.txt"
# Output (optional):
# - "results/[projection/]projection_4-5.[RData/tsv/xlsx]"
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

# projection.dir <- file.path(result_dir, "projection")
projection.dir <- result_dir
if (!dir.exists(projection.dir)) dir.create(projection.dir, recursive = TRUE)

## Parameters setting ---------------------------------------------------------

species_list <- config$species_list
speciesNum <- length(species_list)

# read data -------------------------------------------------------------------

motif.v.list <- list()
for (s in species_list) {
    # cat(s, "\n")
    motif.v.list[[s]] <- read_tsv(
        file.path(SVD.dir, str_glue("{s}_V0_1402.txt")),
        col_select = num_range(prefix = "", range = 1L:10L),
        col_names = as.character(1L:10L),
        col_types = str_dup("d", 10L)
    ) %>%
        as.list() %>%
        enframe(name = "level", value = "vec") %>%
        add_column(species = s, .before = 2)
}

motif.v.compact <- bind_rows(motif.v.list) %>%
    mutate(species = factor(species, levels = species_list)) %>%
    mutate(level = as.integer(level)) %>%
    # filter(level %in% level_list) %>%
    arrange(level, species)

# projection - ggplot (x-y coordinates) ----------------------------------------

for (level_list in list(c(5L, 6L))) {

    # 将物种 species 的第 level 层特征向量
    # 向物种 ref_species 的第4、5个特征向量张成的子空间作投影：
    # 若以物种 ref_species 的第4、5个特征向量为基，
    # 则投影向量的坐标为 x 和 y，模长 length，（4相对4的、5相对5的）旋转角度 angle。
    projection.tibble <-
        expand_grid(
            species = factor(species_list, level = species_list),
            level = level_list,
            ref_species = factor(species_list, level = species_list)
        ) %>%
        # filter(species != ref_species) %>%
        left_join(
            motif.v.compact,
            by = c("species", "level")
        ) %>%
        left_join(
            motif.v.compact %>%
                filter(level == level_list[1]) %>%
                select(!level) %>%
                rename(ref_vec1 = vec),
            by = c("ref_species" = "species")
        ) %>%
        left_join(
            motif.v.compact %>%
                filter(level == level_list[2]) %>%
                select(!level) %>%
                rename(ref_vec2 = vec),
            by = c("ref_species" = "species")
        ) %>%
        # compute x-y coordinates
        rowwise() %>%
        mutate(
            x = vec %*% ref_vec1 %>% as.double(),
            y = vec %*% ref_vec2 %>% as.double()
        ) %>%
        ungroup() %>%
        # compute polar coordinates
        rowwise() %>%
        mutate(
            angle = (atan2(y, x) * (180/pi)),
            length = norm(c(x, y), type = "2")
        ) %>%
        mutate(
            angle = if_else(level == level_list[2], angle - 90, angle)
        ) %>%
        ungroup() %>%
        mutate(level = level - 1L) %>%
        select(species, level, ref_species, x, y, angle, length) %>%
        arrange(species, ref_species, level)

    # write data ------------------------------------------------------------------

    # save(
    #     projection.tibble,
    #     file = file.path(
    #         projection.dir,
    #         str_glue("projection_{level_list[1] - 1}-{level_list[2] - 1}.RData")
    #     )
    # )

    # write_tsv(
    #     x = projection.tibble,
    #     file = file.path(
    #         projection.dir,
    #         str_glue("projection_{level_list[1] - 1}-{level_list[2] - 1}.tsv")
    #     )
    # )

    openxlsx::write.xlsx(
        x = projection.tibble %>%
            mutate(
                x = x %>% round(digits = 3),
                y = y %>% round(digits = 3),
                angle = angle %>% round(digits = 1),
                length = length %>% round(digits = 3)
            ),
        file = file.path(
            projection.dir,
            str_glue("projection_{level_list[1] - 1}-{level_list[2] - 1}.xlsx")
        ),
        overwrite = TRUE,
        colWidths = "auto"
    )

}
