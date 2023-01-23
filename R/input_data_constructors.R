#' `GWAS_data` constructor
#'
#' @param dat Tibble, results from a GWAS analysis, with at least columns
#' `chromosome`, `position` and `score`.
#' @returns A `GWAS_data` object, i.e. a tibble.
new_GWAS_data <- function(dat){
  ## Making sure that the input is a tibble
  stopifnot(tibble::is_tibble(dat))

  class(dat) <- c("GWAS_data", class(dat))

  dat
}

#' Checks validity of input for `GWAS_data` constructor
#'
#' @param x, a `GWAS_data` object constructed via \link{new_GWAS_data}.
#' @returns A `GWAS_data` object, i.e. a tibble.
validate_GWAS_data <- function(x){

  ## A GWAS result table must at least contain these columns
  .check_cols(x, c("chromosome", "position", "score"))

  if(!is.numeric(x[["position"]])) stop("'position' column should contain numeric values.", call. = FALSE)
  if(!is.numeric(x[["score"]])) stop("'score' column should contain numeric values.", call. = FALSE)

  x
}

#' Creates a `GWAS_data` object
#'
#' Creates a `GWAS_data` object from a tibble or data-frame of GWAS results.
#'
#' The input data should have one row per marker, and at least the
#' following columns:
#'
#' * `chromosome`: character column, chromosome on which the marker is located.
#'
#' * `position`: numeric, the physical position of the marker along the chromosome (in bp).
#'
#' * `score` or `padj`: numeric, the GWAS score or adjusted p-value of the marker.
#' If column `score` column is missing, will be constructed as `-log10(padj)`.
#'
#' @param dat Tibble, results from a GWAS analysis. See Details.
#' @param keep_rownames_as Character, the name of the column in which to save the
#' rownames of the input data-frame. Default value is `NULL`, i.e. rownames will
#' be discarded.
#' @returns A `GWAS_data` object, i.e. a tibble.
#' @examples
#' x <- get_example_data()
#'
#' GWAS_data(x[["GWAS"]])
#' @export
GWAS_data <- function(dat, keep_rownames_as = NULL){
  ## If the data is not a tibble, transform it
  if(!tibble::is_tibble(dat)){

    dat <- tibble::as_tibble(dat, rownames = keep_rownames_as)
  }

  ## If missing score column, construct it from adjusted p-value
  if(!("score" %in% colnames(dat))){

    if(!("padj" %in% colnames(dat))) stop("Input data-frame should have either a 'score' or a 'padj' column.")
    if(!is.numeric(dat[["padj"]])) stop("'padj' column in input data-frame should contain numeric values.")

    ## for devtools::check
    score <- padj <- NULL

    dat <- dat |>
      dplyr::mutate(score = -log10(padj))
  }

  dat <- new_GWAS_data(dat)
  validate_GWAS_data(dat)
}




#' `DE_data` constructor
#'
#' @param dat Tibble, results from a differential expression analysis, with at least columns
#' `chromosome`, `score`, `log2FoldChange`, `start`, `end` and `position`.
#' @returns A `DE_data` object, i.e. a tibble.
new_DE_data <- function(dat){
  ## Making sure that the input is a tibble
  stopifnot(tibble::is_tibble(dat))

  class(dat) <- c("DE_data", class(dat))

  dat
}

#' Checks validity of input for `DE_data` constructor
#'
#' @param x, a `DE_data` object constructed via \link{new_DE_data}.
#' @returns A `DE_data` object, i.e. a tibble.
validate_DE_data <- function(x){

  ## A GWAS result table must at least contain these columns
  .check_cols(x, c("chromosome", "start", "end", "position", "score", "log2FoldChange"))

  if(!is.numeric(x[["start"]]) | !is.numeric(x[["end"]])) stop("'start' and 'end' columns should contain numeric values.")
  if(!is.numeric(x[["position"]])) stop("'position' column should contain numeric values.", call. = FALSE)
  if(!is.numeric(x[["score"]])) stop("'score' column should contain numeric values.", call. = FALSE)
  if(!is.numeric(x[["log2FoldChange"]])) stop("'log2FoldChange' column should contain numeric values.", call. = FALSE)

  x
}

#' Creates a `DE_data` object
#'
#' Creates a `DE_data` object from a tibble or data-frame of differential expression
#' results.
#'
#' The input data should have one row per gene or transcript, and at least the
#' following columns:
#'
#' * `chromosome`: character column, chromosome on which the gene/transcript is located.
#'
#' * `start` and `end`: numeric, starting and end position of the gene/transcript
#' (in bp). A column `position` will be constructed as the middle value (mean) between
#' `start` and `end`.
#'
#' * `score` or `padj`: numeric, the DE score or adjusted p-value of the
#' gene/transcript. If column `score` column is missing, will be constructed
#' as `-log10(padj)`.
#'
#' * `foldChange` or `log2FoldChange`: numeric, the fold-change or log2(fold-change)
#' of the gene/transcript. If column `log2FoldChange` is missing, will be constructed
#' as `log2(foldChange)`.
#'
#' @param dat Tibble, results from a differential expression analysis. See Details.
#' @param keep_rownames_as Character, the name of the column in which to save the
#' rownames of the input data-frame. Default value is `NULL`, i.e. rownames will
#' be discarded.
#' @returns A `DE_data` object, i.e. a tibble.
#' @examples
#' x <- get_example_data()
#'
#' DE_data(x[["DE"]])
#' @export
DE_data <- function(dat, keep_rownames_as = NULL){
  ## If the data is not a tibble, transform it
  if(!tibble::is_tibble(dat)){

    dat <- tibble::as_tibble(dat, rownames = keep_rownames_as)
  }

  ## If missing score column, construct it from adjusted p-value
  if(!("score" %in% colnames(dat))){

    if(!("padj" %in% colnames(dat))) stop("Input data-frame should have either a 'score' or a 'padj' column.")
    if(!is.numeric(dat[["padj"]])) stop("'padj' column in input data-frame should contain numeric values.")

    ## for devtools::check
    score <- padj <- NULL

    dat <- dat |>
      dplyr::mutate(score = -log10(padj))
  }

  ## If missing log2FoldChange column, construct it from fold-change
  if(!("log2FoldChange" %in% colnames(dat))){

    if(!("foldChange" %in% colnames(dat))) stop("Input data-frame should have either a 'log2FoldChange' or a 'foldChange' column.")
    if(!is.numeric(dat[["foldChange"]])) stop("'foldChange' column in input data-frame should contain numeric values.")

    ## for devtools::check
    log2FoldChange <- foldChange <- NULL

    dat <- dat |>
      dplyr::mutate(log2FoldChange = log2(foldChange))
  }

  ## Construct the position column from start and end
  if(length(setdiff(c("start", "end"), colnames(dat)))) stop("Input data-frame should have a 'start' and an 'end' column.")
  if(!is.numeric(dat[["start"]]) | !is.numeric(dat[["end"]])) stop("'start' and 'end' columns should contain numeric values.")

  ## for devtools::check
  start <- end <- NULL

  dat <- dat |>
    dplyr::mutate(position = (start + end) / 2)

  dat <- new_DE_data(dat)
  validate_DE_data(dat)
}



#' `CAN_data` constructor
#'
#' @param dat Tibble, containing information about genes of interest, with at least columns
#' `chromosome`, `start`, `end`, `position` and `name`.
#' @returns A `CAN_data` object, i.e. a tibble.
new_CAN_data <- function(dat){
  ## Making sure that the input is a tibble
  stopifnot(tibble::is_tibble(dat))

  class(dat) <- c("CAN_data", class(dat))

  dat
}

#' Checks validity of input for `CAN_data` constructor
#'
#' @param x, a `CAN_data` object constructed via \link{new_CAN_data}.
#' @returns A `CAN_data` object, i.e. a tibble.
validate_CAN_data <- function(x){

  ## A GWAS result table must at least contain these columns
  .check_cols(x, c("chromosome", "start", "end", "position", "name"))


  if(!is.numeric(x[["start"]]) | !is.numeric(x[["end"]])) stop("'start' and 'end' columns should contain numeric values.", call. = FALSE)
  if(!is.numeric(x[["position"]])) stop("'position' column should contain numeric values.", call. = FALSE)
  if(!is.character(x[["name"]])) stop("'name' column should contain character values.", call. = FALSE)

  x
}

#' Creates a `CAN_data` object
#'
#' Creates a `CAN_data` object from a tibble or data-frame of candidate genes.
#'
#' The input data should have one row per gene, and at least the
#' following columns:
#'
#' * `chromosome`: character column, chromosome on which the gene is located.
#'
#' * `start` and `end`: numeric, starting and end position of the gene (in bp).
#' A column `position` will be constructed as the middle value (mean) between
#' `start` and `end`.
#'
#' * `name`: character, the name of the candidate genes to be displayed.
#'
#' @param dat Tibble, set of candidate genes of interest. See Details.
#' @param keep_rownames_as Character, the name of the column in which to save the
#' rownames of the input data-frame. Default value is `NULL`, i.e. rownames will
#' be discarded.
#' @returns A `CAN_data` object, i.e. a tibble.
#' @examples
#' x <- get_example_data()
#'
#' CAN_data(x[["CAN"]])
#' @export
CAN_data <- function(dat, keep_rownames_as = NULL){
  ## If the data is not a tibble, transform it
  if(!tibble::is_tibble(dat)){

    dat <- tibble::as_tibble(dat, rownames = keep_rownames_as)
  }

  ## Construct the position column from start and end
  if(length(setdiff(c("start", "end"), colnames(dat)))) stop("Input data-frame should have a 'start' and an 'end' column.")
  if(!is.numeric(dat[["start"]]) | !is.numeric(dat[["end"]])) stop("'start' and 'end' columns should contain numeric values.")

  ## for devtools::check
  start <- end <- NULL

  dat <- dat |>
    dplyr::mutate(position = (start + end) / 2)

  dat <- new_CAN_data(dat)
  validate_CAN_data(dat)
}
