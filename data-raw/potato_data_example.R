library(tidyverse)
library(GWASpoly)


## ------------------------ ##
## Reading the GWAS results ##
## ------------------------ ##

load("W:/hrpoab/clever_culling/Potato/olivia_paper_1_2021/genomics_gwas/processed_data/data2_K_Qstructure_threshold.RData")

models <- colnames(data2_K_Qstructure_threshold@scores$bruising_score_mean)

gwas_data <- data2_K_Qstructure_threshold@scores$bruising_score_mean |>
  as_tibble(rownames = "id") |>
  mutate(score = base::pmax(!!!rlang::syms(models), na.rm = TRUE),
         pvalue = 10^(-score),
         chromosome = str_extract(id, "ST4.03ch\\d{2}"),
         position = as.numeric(str_extract(id, "(?<=_)\\d+"))) |>
  select(id, chromosome, position, score, pvalue)


## ----------------------------------------------------------- ##
## Reading the transcriptomics differential expression results ##
## ----------------------------------------------------------- ##

de_results <- read_csv("W:/hrpoab/clever_culling/Potato/olivia_paper_1_2021/transcriptomics_differential_expression/processed_data/res_DE_LFCshrinkage.csv")
genes_annotation <- read_csv("W:/hrpoab/clever_culling/Potato/olivia_paper_1_2021/transcriptomics_differential_expression/processed_data/descr_annot.csv")

de_data <- de_results |>
  left_join(genes_annotation, by = c("gene" = "gene_id")) |>
  mutate(score = -log10(padj)) |>
  select(gene,
         chromosome = chrom,
         position = pos,
         score,
         padj,
         log2FoldChange = log2FoldChange_shrinkage,
         start = tx_start,
         end = tx_end,
         label = description)


## --------------------------------------- ##
## Reading the candidate genes information ##
## --------------------------------------- ##

candidate_genes <- read_tsv("W:/hrpoab/clever_culling/Potato/olivia_paper_1_2021/genomics_gwas/known_qtls.tsv")

candidate_data <- candidate_genes |>
  filter(Trait == "Bruising") |>
  mutate(position = (Start + End) / 2) |>
  select(id = gene_id,
         chromosome = Chromosome,
         position,
         start = Start,
         end = End,
         name = gene_name_short,
         gene_name)


## --------------------------------------- ##
## Reading the candidate genes information ##
## --------------------------------------- ##

chromosome_info <- list(gwas_data, de_data, candidate_data) |>
  map(select, chromosome, position) |>
  reduce(bind_rows) |>
  group_by(chromosome) |>
  summarise(length = max(position)) |>
  mutate(length = ceiling(length))

usethis::use_data(gwas_data, de_data, candidate_data, chromosome_info, internal = TRUE, overwrite = TRUE)
