#' GWAS_data constructor
#'
#' @param dat Tibble, results from a GWAS analysis, with at least columns
#' `chromosome`. `position` and `score`.
#' @param title Character, the name of the analysis. Default value is `"GWAS"`.
#' @returns A `GWAS_data` object, i.e. a tibble.
#' @export
new_GWAS_data <- function(dat, title = "GWAS"){
  ## Making sure that the input is a tibble
  stopifnot(tibble::is_tibble(dat))
  stopifnot(is.character(title))

  class(dat) <- c("GWAS_data", class(dat))
  attr(dat, "title", title)

  dat
}

#' Check validity of input for GWAS_data constructor
#'
#' @param x, a GWAS_data object constructed via \link{new_GWAS_data}.
#' @returns A `GWAS_data` object, i.e. a tibble.
validate_GWAS_data <- function(x){
  ## A GWAS result table must at least contain these columns
  if(!all(c("chromosome", "position", "score") %in% colnames(x))) stop("Data frame should have the following columns: 'chromosome', 'position', 'score'.", call. = FALSE)

  if(!is.numeric(x[["position"]])) stop("'position' column should contain numeric values.", call. = FALSE)
  if(!is.numeric(x[["score"]])) stop("'score' column should contain numeric values.", call. = FALSE)

  x
}

#' Creates a GWAS_data object
#'
#' The input data should have one row per marker, and at least the
#' following columns:
#'
#' * `chromosome`: character column, chromosome on which the marker is located.
#'
#' * `position`: numeric, the physical position of the marker along the chromosome.
#'
#' * `score`: numeric, the GWAS score of the marker.
#'
#' @param dat Tibble, results from a GWAS analysis. See Details.
#' @param title Character, the name of the analysis. Default value is `"GWAS"`.
#' @returns A `GWAS_data` object, i.e. a tibble.
#' @export
GWAS_data <- function(dat, title = "GWAS"){
  ## If the data is not a tibble, transform it
  if(!tibble::is_tibble(dat)){

    ## Making sure that we are conserving
    if(is.null(rownames(dat))) rownames_arg <- NULL else rownames_arg <- "rowname"

    dat <- tibble::as_tibble(dat, rownames_arg)
  }

  dat <- new_GWAS_data(dat, title)
  validate_GWAS_data(dat)
}

