#' Checks whether some columns are present in a tibble
#'
#' @param x Tibble
#' @param col_names character vector of column names
#' @param param_name Character, name of the dataframe to use in the error message.
#' @returns invisible NULL
.check_cols <- function(x, col_names, param_name = "Input data-frame"){

  missing_cols <- setdiff(col_names, colnames(x))

  if(length(missing_cols)){
    stop("", param_name, " is missing the following columns: '", paste0(missing_cols, collapse = "', '"), "'.", call. = FALSE)
  }

  return(invisible(NULL))
}

#' Example dataset
#'
#' Returns a list of example datasets.
#'
#' The dataset used in this example is presented in Angelin-Bonnet et al., 2023
#' (submitted). In this study, tetraploid potato plants from a half-sibling
#' breeding population were used to assess the genetic components of tuber
#' bruising. Capture sequencing was used to obtain genomic information about the
#' individuals, and a genome-wide association study (GWAS) was performed on 72,847
#' genomic biallelic variants obtained from 158
#' plants for which a bruising score was measured. The GWAS analysis was carried
#' with the `GWASpoly` package. In addition, expression data was obtained for
#' 25,163 transcribed genes, and a differential expression (DE) analysis was
#' carried out between 41 low- and 33 high-bruising samples. Finally, a literature
#' search yielded a list of 42 candidate genes identified in previous studies as
#' involved in potato tuber bruising mechanisms. A subset of the GWAS and DE results,
#' as well as the list of candidate genes from the literature, are made available
#' in this function. From the
#' complete GWAS results table, half of the genomic variants with a GWAS score
#' < 3.5 were randomly selected and consequently discarded, yielding a dataset
#' with GWAS scores for 35,481 variants. Similarly, half of the transcribed
#' genes in the DE results table with an adjusted p-value > 0.05 were
#' randomly selected and discarded, yielding a dataset with DE results for
#' 10,671 transcribed genes. This filtering was performed to reduce the
#' size of the datasets (in accordance with CRAN policies), but ensures
#' that all significant markers and genes are retained in the datasets.
#' Finally, some of the candidate genes located on chromosome 3 were removed
#' from the example dataset for better clarity in the resulting HIDECAN
#' plot, leaving 32 candidate genes.
#'
#' @returns A list with the following elements:
#' * `GWAS`: a tibble of GWAS results, with columns `id`, `chromosome`,
#' `position` and `score`.
#' * `DE`: a tibble of differential expression results, with columns `gene`,
#' `chromosome`, `padj`, `log2FoldChange`, `start`, `end` and `label`.
#' * `CAN`: a tibble of candidate genes, with columns `id`, `chromosome`,
#' `start`, `end`, `name` and `gene_name`.
#' @export
get_example_data <- function(){
  list(
    GWAS = gwas_data,
    DE = de_data,
    CAN = candidate_data
  )
}

# #' GWASpoly example dataset
# #'
# #' Returns a list of GWASpoly example datasets. Not exported because otherwise
# #' would load the GWASpoly package which is not on CRAN.
# #'
# #' Install the GWASpoly package with `devtools::install_github('jendelman/GWASpoly', build_vignettes=FALSE)`.
# #'
# #' @param with_thresholds Logical, should the GWASpoly object returned
# #' contain significance threshold? Default value is `TRUE`.
# #' @returns
# #' * if `with_thresholds` is `TRUE`: a `GWASpoly.thresh` object
# #' (returned by the `GWASpoly::set.threshold()` function).
# #' if `with_thresholds` is `FALSE`: a `GWASpoly.fitted` object
# #' (returned by the `GWASpoly::GWASpoly()` function).
# #' @export
# get_gwaspoly_example_data <- function(with_thresholds = TRUE){
#
#   if(with_thresholds){
#     return(gwaspoly_res_thr)
#   } else {
#     return(gwaspoly_res)
#   }
# }
