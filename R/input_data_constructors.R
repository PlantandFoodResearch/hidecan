#' GWAS_data constructor
#'
#' @param dat Tibble, results from a GWAS analysis, with at least columns
#' `chromosome`, `position` and `score`.
#' @returns A `GWAS_data` object, i.e. a tibble.
#' @export
new_GWAS_data <- function(dat){
  ## Making sure that the input is a tibble
  stopifnot(tibble::is_tibble(dat))

  class(dat) <- c("GWAS_data", class(dat))

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
#' * `position`: numeric, the physical position of the marker along the chromosome (in bp).
#'
#' * `score`: numeric, the GWAS score of the marker.
#'
#' @param dat Tibble, results from a GWAS analysis. See Details.
#' @param keep_rownames_as Character, the name of the column in which to save the
#' rownames of the input data-frame. Default value is `NULL`, i.e. rownames will
#' be discarded.
#' @returns A `GWAS_data` object, i.e. a tibble.
#' @export
GWAS_data <- function(dat, keep_rownames_as = NULL){
  ## If the data is not a tibble, transform it
  if(!tibble::is_tibble(dat)){

    dat <- tibble::as_tibble(dat, rownames = keep_rownames_as)
  }

  dat <- new_GWAS_data(dat)
  validate_GWAS_data(dat)
}




#' DE_data constructor
#'
#' @param dat Tibble, results from a differential expression analysis, with at least columns
#' `chromosome`, `score`, `log2FoldChange`, `start` and `end`.
#' @returns A `DE_data` object, i.e. a tibble.
#' @export
new_DE_data <- function(dat){
  ## Making sure that the input is a tibble
  stopifnot(tibble::is_tibble(dat))

  class(dat) <- c("DE_data", class(dat))

  dat
}

#' Check validity of input for DE_data constructor
#'
#' @param x, a DE_data object constructed via \link{new_DE_data}.
#' @returns A `DE_data` object, i.e. a tibble.
validate_DE_data <- function(x){

  ## A GWAS result table must at least contain these columns
  .check_cols(x, c("chromosome", "start", "end", "score", "log2FoldChange"))


  if(!is.numeric(x[["start"]]) | !is.numeric(x[["end"]])) stop("'start' and 'end' columns should contain numeric values.", call. = FALSE)
  if(!is.numeric(x[["score"]])) stop("'score' column should contain numeric values.", call. = FALSE)
  if(!is.numeric(x[["log2FoldChange"]])) stop("'log2FoldChange' column should contain numeric values.", call. = FALSE)

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
#' * `log2FoldChange`: numeric, the log2(fold-change) of the gene/transcript.
#'
#' * `start`: numeric, starting position of the gene/transcript (in bp).
#'
#' * `end`: numeric, end position of the gene/transcript (in bp).
#'
#' @param dat Tibble, results from a GWAS analysis. See Details.
#' @param keep_rownames_as Character, the name of the column in which to save the
#' rownames of the input data-frame. Default value is `NULL`, i.e. rownames will
#' be discarded.
#' @returns A `DE_data` object, i.e. a tibble.
#' @export
DE_data <- function(dat, keep_rownames_as = NULL){
  ## If the data is not a tibble, transform it
  if(!tibble::is_tibble(dat)){

    dat <- tibble::as_tibble(dat, rownames = keep_rownames_as)
  }

  dat <- new_DE_data(dat)
  validate_DE_data(dat)
}



#' CAN_data constructor
#'
#' @param dat Tibble, containing information about genes of interest, with at least columns
#' `chromosome`, `start`, `end` and `name`.
#' @returns A `CAN_data` object, i.e. a tibble.
#' @export
new_CAN_data <- function(dat){
  ## Making sure that the input is a tibble
  stopifnot(tibble::is_tibble(dat))

  class(dat) <- c("CAN_data", class(dat))

  dat
}

#' Check validity of input for CAN_data constructor
#'
#' @param x, a CAN_data object constructed via \link{new_CAN_data}.
#' @returns A `CAN_data` object, i.e. a tibble.
validate_CAN_data <- function(x){

  ## A GWAS result table must at least contain these columns
  .check_cols(x, c("chromosome", "start", "end", "name"))


  if(!is.numeric(x[["start"]]) | !is.numeric(x[["end"]])) stop("'start' and 'end' columns should contain numeric values.", call. = FALSE)
  if(!is.character(x[["name"]])) stop("'name' column should contain character values.", call. = FALSE)

  x
}

#' Creates a CAN_data object
#'
#' The input data should have one row per gene, and at least the
#' following columns:
#'
#' * `chromosome`: character column, chromosome on which the gene is located.
#'
#' * `start`: numeric, starting position of the gene (in bp).
#'
#' * `end`: numeric, end position of the gene (in bp).
#'
#' * `name`: character, the name of the candidate genes to be displayed.
#'
#' @param dat Tibble, set of candidate genes of interest. See Details.
#' @param keep_rownames_as Character, the name of the column in which to save the
#' rownames of the input data-frame. Default value is `NULL`, i.e. rownames will
#' be discarded.
#' @returns A `CAN_data` object, i.e. a tibble.
#' @export
CAN_data <- function(dat, keep_rownames_as = NULL){
  ## If the data is not a tibble, transform it
  if(!tibble::is_tibble(dat)){

    dat <- tibble::as_tibble(dat, rownames = keep_rownames_as)
  }

  dat <- new_CAN_data(dat)
  validate_CAN_data(dat)
}
