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

#' Get set of example data
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
