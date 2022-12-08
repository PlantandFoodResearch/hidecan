#' Filters GWAS or DE results based on a threshold
#'
#' Filters markers or genes/transcripts based on a threshold applied to their
#' GWAS or DE score, and log2(fold-change) (if applicable). For a set of
#' candidate genes, simply returns the list. Note that markers or genes with
#' a missing score or log2(fold-change) will be removed from the dataset.
#'
#' @param x Either a `GWAS_data`, `DE_data` or `CAN_data` object.
#' @param score_thr Numeric, threshold to use on markers' or genes/transcripts' score.
#' Only markers or genes with a score equal to or higher than this threshold
#' will be retained. Default value is 0. Ignored for `CAN_data`.
#' @param log2fc_thr Numeric, threshold to use on the absolute value of genes/
#' transcripts' log2(fold-change). Only genes/transcripts with an absolute
#' log2(fold-change) equal to or higher than this threshold will be retained.
#' Ignored for `GWAS_data` and `CAN_data`.
#' @returns A filtered tibble (of class `GWAS_data_thr`, `DE_data_thr` or
#' `CAN_data_thr`).
#' @examples
#' x <- get_example_data()
#'
#' ## For GWAS results
#' apply_threshold(GWAS_data(x[["GWAS"]]), score_thr = 4)
#'
#' ## For DE results - in second line, no threshold is applied
#' ## on the log2(fold-change)
#' apply_threshold(DE_data(x[["DE"]]), score_thr = -log10(0.05), log2fc_thr = 1)
#' apply_threshold(DE_data(x[["DE"]]), score_thr = -log10(0.05), log2fc_thr = 0)
#'
#' ## No effect on the Candidate genes
#' apply_threshold(CAN_data(x[["CAN"]]))
#' @export
apply_threshold <- function(x, score_thr = 0, log2fc_thr = 0){
  UseMethod("apply_threshold")
}

#' @rdname apply_threshold
#' @export
apply_threshold.GWAS_data <- function(x, score_thr = 0, log2fc_thr = 0){

  score <- NULL

  res <- x |>
    dplyr::filter(score >= score_thr)

  class(res)[1] <- "GWAS_data_thr"

  return(res)
}

#' @rdname apply_threshold
#' @export
apply_threshold.DE_data <- function(x, score_thr = 0, log2fc_thr = 0){

  score <- log2FoldChange <- NULL

  res <- x |>
    dplyr::filter(score >= score_thr,
                  abs(log2FoldChange) >= log2fc_thr)

  class(res)[1] <- "DE_data_thr"

  return(res)
}

#' @rdname apply_threshold
#' @export
apply_threshold.CAN_data <- function(x, score_thr = 0, log2fc_thr = 0){

  res <- x

  class(res)[1] <- "CAN_data_thr"

  return(res)
}

#' @rdname apply_threshold
#' @export
apply_threshold.default <- function(x, score_thr = 0, log2fc_thr = 0){

  x

}
