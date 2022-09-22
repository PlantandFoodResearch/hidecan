#' Checks whether some columns are present in a tibble
#'
#' @param x Tibble
#' @param col_names character vector of column names
#' @param param_name Character, name of the dataframe to use in the error message.
#' @returns invisible NULL
.check_cols <- function(x, col_names, param_name = "Input data-frame"){

  missing_cols <- setdiff(col_names, colnames(x))

  if(length(missing_cols)){
    stop("", param_name, " is missing the following columns: '", paste0(missing_cols, collapse = "', '"), "'.", .call = FALSE)
  }

  return(invisible(NULL))
}
