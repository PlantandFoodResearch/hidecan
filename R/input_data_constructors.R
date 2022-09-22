#' GWAS_data constructor
#'
#' @param dat Tibble, results from a GWAS analysis, with at least columns
#' `chromosome`, `position` and `score`.
#' @param title Character, the name of the analysis. Default value is `"GWAS"`.
#' @returns A `GWAS_data` object, i.e. a tibble.
#' @export
new_GWAS_data <- function(dat, title = "GWAS"){
  ## Making sure that the input is a tibble
  stopifnot(tibble::is_tibble(dat))
  stopifnot(is.character(title))

  class(dat) <- c("GWAS_data", class(dat))
  attr(dat, "title") <- title

  dat
}

#' Check validity of input for GWAS_data constructor
#'
#' @param x, a GWAS_data object constructed via \link{new_GWAS_data}.
#' @returns A `GWAS_data` object, i.e. a tibble.
validate_GWAS_data <- function(x){

  ## A GWAS result table must at least contain these columns
  .check_cols(x, c("chromosome", "position", "score"))

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
#' @param keep_rownames_as Character, the name of the column in which to save the
#' rownames of the input data-frame. Default value is `NULL`, i.e. rownames will
#' be discarded.
#' @returns A `GWAS_data` object, i.e. a tibble.
#' @export
GWAS_data <- function(dat, title = "GWAS", keep_rownames_as = NULL){
  ## If the data is not a tibble, transform it
  if(!tibble::is_tibble(dat)){

    dat <- tibble::as_tibble(dat, rownames = keep_rownames_as)
  }

  dat <- new_GWAS_data(dat, title)
  validate_GWAS_data(dat)
}






#' DE_data constructor
#'
#' @param dat Tibble, results from a differential expression analysis, with at least columns
#' `chromosome`, `score` and either `position` or both `start` and `end`.
#' @param title Character, the name of the analysis. Default value is `"DE"`.
#' @returns A `DE_data` object, i.e. a tibble.
#' @export
new_DE_data <- function(dat, title = "DE"){
  ## Making sure that the input is a tibble
  stopifnot(tibble::is_tibble(dat))
  stopifnot(is.character(title))

  class(dat) <- c("DE_data", class(dat))
  attr(dat, "title") <- title

  dat
}

#' Check validity of input for DE_data constructor
#'
#' @param x, a DE_data object constructed via \link{new_DE_data}.
#' @returns A `DE_data` object, i.e. a tibble.
validate_DE_data <- function(x){

  ## A GWAS result table must at least contain these columns
  .check_cols(x, c("chromosome", "score"))

  if(!( ("position" %in% colnames(x)) | all(c("start", "end") %in% colnames(x)) )){
    stop("Input data-frame should contain either a 'position' column or both a 'start' and 'end' columns.")
  }

  if("position" %in% colnames(x)){
    if(!is.numeric(x[["position"]])) stop("'position' column should contain numeric values.", call. = FALSE)

    if(all(c("start", "end") %in% colnames(x))) warning("'position' column takes precedence over 'start' and 'end' columns.")
  } else {
    if(!is.numeric(x[["start"]]) | !is.numeric(x[["end"]])) stop("'start' and 'end' columns should contain numeric values.", call. = FALSE)
  }


  if(!is.numeric(x[["score"]])) stop("'score' column should contain numeric values.", call. = FALSE)

  x
}

#' Creates a DE_data object
#'
#' The input data should have one row per gene or transcript, and at least the
#' following columns:
#'
#' * `chromosome`: character column, chromosome on which the gene/transcript is located.
#'
#' * `score`: numeric, the DE score of the gene/transcript.
#'
#' In addition, the input data should have either:
#'
#' * a `start` and an `end` column: numeric, start and end positions of the gene/transcript (in Mb).
#'
#' * a `position` column: numeric, the position of the gene/transcript (in Mb). Could represent
#' the middle of the gene or transcript.
#'
#' @param dat Tibble, results from a GWAS analysis. See Details.
#' @param title Character, the name of the analysis. Default value is `"GWAS"`.
#' @param keep_rownames_as Character, the name of the column in which to save the
#' rownames of the input data-frame. Default value is `NULL`, i.e. rownames will
#' be discarded.
#' @returns A `GWAS_data` object, i.e. a tibble.
#' @export

