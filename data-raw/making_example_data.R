library(tidyverse)
library(GWASpoly)

set.seed(348)

## ------------------------ ##
## Reading the GWAS results ##
## ------------------------ ##

load("W:/hrpoab/clever_culling/Potato/olivia_paper_1_2021/genomics_gwas/processed_data/data2_K_Qstructure_threshold.RData")

models <- colnames(data2_K_Qstructure_threshold@scores$bruising_score_mean)

gwas_data <- data2_K_Qstructure_threshold@scores$bruising_score_mean |>
  as_tibble(rownames = "id") |>
  mutate(score = base::pmax(!!!rlang::syms(models), na.rm = TRUE),
         chromosome = str_extract(id, "ST4.03ch\\d{2}"),
         position = as.numeric(str_extract(id, "(?<=_)\\d+"))) |>
  select(id, chromosome, position, score)


## ----------------------------------------------------------- ##
## Reading the transcriptomics differential expression results ##
## ----------------------------------------------------------- ##

de_results <- read_csv("W:/hrpoab/clever_culling/Potato/olivia_paper_1_2021/transcriptomics_differential_expression/processed_data/res_DE_LFCshrinkage.csv")
genes_annotation <- read_csv("W:/hrpoab/clever_culling/Potato/olivia_paper_1_2021/transcriptomics_differential_expression/processed_data/descr_annot.csv")

de_data <- de_results |>
  left_join(genes_annotation, by = c("gene" = "gene_id")) |>
  select(gene,
         chromosome = chrom,
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
  select(id = gene_id,
         chromosome = Chromosome,
         start = Start,
         end = End,
         name = gene_name_short,
         gene_name) |>
  ## Get rid of all the "PPO" labels, otherwise too many for ggrepel
  mutate(name = case_when(name == "PPO" & start == 45676042 ~ "PPOs",
                          name == "PPO" & start != 45676042 ~ NA_character_,
                          TRUE ~ name)) |>
  ## removing some candidate genes from chromosome 3 to avoid ggrepel warning
  filter(!(chromosome == "ST4.03ch03") | name %in% c("PHO1A", "4CL", "KTPI", "4CL2"))

## ---------------------------------------- ##
## Constructing example QTL mapping results ##
## ---------------------------------------- ##

qtl_data <- tibble(
  chromosome = c("ST4.03ch01", "ST4.03ch01", "ST4.03ch04", "ST4.03ch05",
                 "ST4.03ch07", "ST4.03ch08", "ST4.03ch11"),
  start = c(5.6, 9.5, 12.7, 37.6, 49.9, 45, 21.4),
  end = c(11.3, 14, 17.6, 39.9, 51.7, 47, 26.7),
  score = c(4.1, 5.5, 4.3, 3.2, 5.1, 5.8, 2.1),
) |>
  mutate(
    id = paste0("qtl_", 1:n()),
    start = start * 1e6,
    end = end * 1e6,
    name = paste0("QTL ", 1:n())
  ) |>
  relocate(id, .before = everything())

## --------------------------------- ##
## Getting example data for GWASpoly ##
## --------------------------------- ##

## Code taken from the GWASpoly vignette (https://jendelman.github.io/GWASpoly/GWASpoly.html)
## with the first example data from GWASpoly

genofile <- system.file("extdata", "TableS1.csv", package = "GWASpoly")
phenofile <- system.file("extdata", "TableS2.csv", package = "GWASpoly")

data <- read.GWASpoly(ploidy = 4,
                      pheno.file = phenofile,
                      geno.file = genofile,
                      format = "ACGT",
                      n.traits = 13,
                      delim = ",")

data.original <- set.K(data,
                       LOCO = FALSE,
                       n.core = 2)

gwaspoly_res <- GWASpoly(data.original,
                         models = c("general", "additive", "1-dom"),
                         traits = c("tuber_eye_depth", "tuber_shape", "sucrose"),
                         n.core = 2)

gwaspoly_res_thr <- set.threshold(gwaspoly_res, method = "M.eff", level = 0.05)

## For users to access the example GWASpoly dataset without having it loaded
saveRDS(gwaspoly_res_thr, file = "inst/extdata/gwaspoly_res_thr.rda")

#manhattan.plot(gwaspoly_res_thr)


## ----------------- ##
## Saving everything ##
## ----------------- ##

## To reduce the size of the dataset, removing at random some of the non-significant
## markers and genes

## Preserving the end of each chromosome
gwas_data_ends <- gwas_data |>
  group_by(chromosome) |>
  slice_max(position, n = 1) |>
  ungroup()

de_data_ends <- de_data |>
  group_by(chromosome) |>
  slice_max(end, n = 1) |>
  ungroup()

## Keeping the significant markers/genes
gwas_data_sign <- gwas_data |>
  filter(score >= 3.5)

de_data_sign <- de_data |>
  filter(padj <= 0.05)

## Randomly subsetting non-significant markers/genes
gwas_data_subset <- gwas_data |>
  filter(score < 3.5) |>
  slice_sample(prop = 0.5)

de_data_subset <- de_data |>
  filter(padj > 0.05) |>
  slice_sample(prop = 0.5)

## Combining retained markers/genes
gwas_data <- bind_rows(
  gwas_data_ends,
  gwas_data_sign,
  gwas_data_subset
)

de_data <- bind_rows(
  de_data_ends,
  de_data_sign,
  de_data_subset
)

usethis::use_data(gwas_data, de_data, candidate_data, qtl_data, internal = TRUE, overwrite = TRUE)
