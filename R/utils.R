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
