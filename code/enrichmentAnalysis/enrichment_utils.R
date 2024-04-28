
enrichment_wilcoxon_test <- function(pathway.data, loading.data) {

    # 富集计算方法：
    # 以如下两组基因的载荷作为样本做 Wilcoxon 两样本秩和检验：
    # （1）基因特征向量全部基因和某个通路的交；
    # （2）基因特征向量全部基因减该通路；
    # 注意：每个通路正极的 p-value 和负极的 p-value 和为1。

    pathway.data <- pathway.data %>%
        select(
            id,
            geneNum,
            description,
            gene
        ) %>%
        mutate(
            gene = gene %>% str_split(",")
        )

    vector_gene_num <- nrow(loading.data)

    result <- pathway.data %>%
        rename(
            pathwayGeneNum = geneNum,
            pathwayGene = gene
        ) %>%
        rowwise() %>%
        mutate(
            index = list(loading.data$Symbol %in% pathwayGene),
            # index = list(toupper(loading.data$Symbol) %in% toupper(pathwayGene)),
            num.of.targets = sum(index),
            num.of.non.targets = vector_gene_num - num.of.targets,
            pathwayGene_in_geneEigenvectors = loading.data$Symbol[index] %>%
                # sort() %>%
                str_c(collapse = ",") %>%
                list(),
            pathwayGene = pathwayGene %>%
                str_c(collapse = ",") %>%
                list(),
            # 不可用if_else，ifelse短路但if_else不短路，
            # 在condition为FALSE的位置依然会计算true参数
            pvalue.pos = ifelse(
                (num.of.targets >= 1L) && (num.of.non.targets >= 1L),
                wilcox.test(
                    loading.data$Loading[index],
                    loading.data$Loading[!index],
                    "greater"
                )$p.val,
                NA_real_
            ),
            pvalue.neg = ifelse(
                (num.of.targets >= 1L) && (num.of.non.targets >= 1L),
                wilcox.test(
                    loading.data$Loading[index],
                    loading.data$Loading[!index],
                    "less"
                )$p.val,
                NA_real_
            )
        ) %>%
        ungroup() %>%
        select(all_of(outputCols))

    result

}

do_enrichment <- function(s_f_g, loading.dir) {

    s <- s_f_g$species
    f <- s_f_g$family
    g <- s_f_g$gene

    loading.data <- read_tsv(
        file.path(
            loading.dir,
            str_glue("{s}_{g}_loadings.txt")
        ),
        col_types = "cd"
    )

    result <- list()
    for (p in pathway_list) {
        result[[p]] <- enrichment_wilcoxon_test(
            pathway.data = pathway_data.list[[f]][[p]],
            loading.data = loading.data
        )
    }

    result

}

write_raw_result <- function(enrichment_result, raw_result.dir) {

    # # 单文件存储
    # list.save(
    #     enrichment_result,
    #     file = file.path(raw_result.dir, "enrichment_result.RData")
    # )

    # 分物种存储
    for (s in species_list) {
        list.save(
            enrichment_result[[s]],
            file = file.path(
                raw_result.dir,
                str_glue("enrichment_result_{s}.RData")
            )
        )
    }

}

read_raw_result <- function(raw_result.dir) {

    # 单文件读取
    # enrichment_result <- list.load(
    #     file = file.path(raw_result.dir, "enrichment_result.RData")
    # )

    # 分物种读取
    enrichment_result <- list()
    for (s in species_list) {
        enrichment_result[[s]] <- list.load(
            file = file.path(
                raw_result.dir,
                str_glue("enrichment_result_{s}.RData")
            )
        )
    }

    enrichment_result

}

mapColNames2Width <- function(colNames) {
    width <- case_when(
        colNames == "id" ~ 12L,
        colNames == "database" ~ 10L,
        colNames == "description" ~ 50L,
        colNames %>% str_detect("pvalue") ~ 10L,
        colNames == "pathwayGeneNum" ~ 16L,
        colNames == "num.of.targets" ~ 16L,
        colNames == "num.of.non.targets" ~ 16L,
        colNames == "pathwayGene" ~ 100L,
        colNames == "pathwayGene_in_geneEigenvectors" ~ 100L
    )
    width
}

filterGeneNum <- function(enrichment_result, minCommonGeneNum, maxCommonGeneNum) {
    # 筛选出通路和基因特征向量里公共基因个数在
    # [minCommonGeneNum, maxCommonGeneNum] 范围内的通路
    enrichment_result %>%
        filter(
            num.of.targets %>% between(minCommonGeneNum, maxCommonGeneNum)
        )
}

trimLongGeneString <- function(enrichment_result, trimGeneNum) {
    # 截取至多 trimGeneNum 个基因的名字存入字符串，
    # 防止字符串太长导致存入Excel单元格时溢出报错
    enrichment_result %>%
        mutate(
            pathwayGene = pathwayGene %>%
                str_split(",") %>%
                map(~ .[1:min(length(.), trimGeneNum)]) %>%
                map(~ str_c(., collapse = ",")),
            pathwayGene_in_geneEigenvectors = pathwayGene_in_geneEigenvectors %>%
                str_split(",") %>%
                map(~ .[1:min(length(.), trimGeneNum)]) %>%
                map(~ str_c(., collapse = ","))
        )
}

sortGene <- function(enrichment_result) {
    enrichment_result %>%
        mutate(
            # 对 pathwayGene 列按基因首字母排序
            pathwayGene = pathwayGene %>%
                str_split(",") %>%
                map(., sort) %>%
                map(., str_c, collapse = ","),
            # 对 pathwayGene_in_geneEigenvectors 列按基因首字母排序
            pathwayGene_in_geneEigenvectors = pathwayGene_in_geneEigenvectors %>%
                str_split(",") %>%
                map(., sort) %>%
                map(., str_c, collapse = ",")
        )
}

formatResult <- function(enrichment_result, s) {
    result.temp <- enrichment_result[[s]]
    for (g in gene_list) {
        for (p in pathway_list) {
            result.temp[[g]][[p]] <- result.temp[[g]][[p]] %>%
                arrange(id) %>%
                add_column(database = p, .before = 1)
        }
        result.temp[[g]] <- result.temp[[g]] %>%
            bind_rows()
    }
    result <- result.temp[[gene_list[1]]] %>%
        select(!c(pvalue.pos, pvalue.neg))
    for (g in gene_list) {
        result <- result %>%
            mutate(
                !!str_glue("pvalue.{str_sub(g, 6L, -1L)}+") := result.temp[[g]]$pvalue.pos,
                !!str_glue("pvalue.{str_sub(g, 6L, -1L)}-") := -result.temp[[g]]$pvalue.neg,
                .before = description
            )
    }
    result <- result %>%
        # filterGeneNum(minCommonGeneNum, maxCommonGeneNum) %>%
        # sortGene() %>%
        # trimLongGeneString(trimGeneNum) %>%
        select(all_of(outputCols))
    result
}

saveResult <- function(result, s, format_result.dir) {
    wb <- createWorkbook()
    addWorksheet(wb, "Sheet 1")
    setColWidths(
        wb = wb,
        sheet = "Sheet 1",
        cols = seq_along(outputCols),
        widths = mapColNames2Width(outputCols)
    )
    writeData(wb, "Sheet 1", result)

    colorCols <- str_which(outputCols, "pvalue")
    colorRows <- 1L + seq_len(nrow(result))
    # 0 < pos <= 0.01
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2>0,C2<=0.01)",
        style = createStyle(bgFill = bgFill$pos[1]),
        type = "expression"
    )
    # 0.01 < pos <= 0.05
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2>0.01,C2<=0.05)",
        style = createStyle(bgFill = bgFill$pos[2]),
        type = "expression"
    )
    # 0.05 < pos <= 0.10
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2>0.05,C2<=0.10)",
        style = createStyle(bgFill = bgFill$pos[3]),
        type = "expression"
    )
    # 0 < neg <= 0.01
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2<0,C2>=-0.01)",
        style = createStyle(bgFill = bgFill$neg[1]),
        type = "expression"
    )
    # 0.01 < neg <= 0.05
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2<-0.01,C2>=-0.05)",
        style = createStyle(bgFill = bgFill$neg[2]),
        type = "expression"
    )
    # 0.05 < neg <= 0.10
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2<-0.05,C2>=-0.10)",
        style = createStyle(bgFill = bgFill$neg[3]),
        type = "expression"
    )
    saveWorkbook(
        wb = wb,
        file = file.path(
            format_result.dir,
            str_glue("{s}.xlsx")
        ),
        overwrite = T
    )
}
