#' Computes chromosomes' length
#'
#' Computes the length (in bp) of each chromosome as the maximum
#' position of markers or genes on the chromosome.
#'
#' @param x Either a `GWAS_data`, `DE_data` or `CAN_data` object.
#' @returns A tibble with two columns: `chromosome` (chromosome name) and
#' `length` (chromosome length in base pair).
#' @examples
#' x <- get_example_data()
#'
#' compute_chrom_length(GWAS_data(x[["GWAS"]]))
#' compute_chrom_length(DE_data(x[["DE"]]))
#' compute_chrom_length(CAN_data(x[["CAN"]]))
#' @export
compute_chrom_length <- function(x){
  UseMethod("compute_chrom_length")
}

#' @rdname compute_chrom_length
#' @export
compute_chrom_length.GWAS_data <- function(x){

  chromosome <- position <- NULL

  x |>
    dplyr::group_by(chromosome) |>
    dplyr::summarise(length = max(position),
                     .groups = "drop")

}

#' @rdname compute_chrom_length
#' @export
compute_chrom_length.DE_data <- function(x){

  .compute_chrom_length_genes(x)

}

#' @rdname compute_chrom_length
#' @export
compute_chrom_length.CAN_data <- function(x){

  .compute_chrom_length_genes(x)

}

#' Computes chromosomes' length for a tibble of genes
#'
#' Computes the length (in bp) of each chromosome as the maximum
#' position of genes on the chromosome.
#'
#' @param x Either a `DE_data` or `CAN_data` object.
#' @returns A tibble with two columns: `chromosome` (chromosome name) and
#' `length` (chromosome length in base pair).
#' @export
.compute_chrom_length_genes <- function(x){

  chromosome <- start <- end <- position <- NULL

  x |>
    ## To handle the case where start > end
    dplyr::mutate(position = purrr::map2_dbl(start, end, max)) |>
    dplyr::group_by(chromosome) |>
    dplyr::summarise(length = max(position),
                     .groups = "drop")
}


#' Computes chromosomes' length from list
#'
#' Computes the length (in bp) of each chromosome from a list of GWAS and
#' DE results as well as candidate gene lists.
#'
#' @param x A list of `GWAS_data`, `DE_data` or `CAN_data` objects.
#' @returns A tibble with two columns: `chromosome` (chromosome name) and
#' `length` (chromosome length in base pair).
#' @examples
#' x <- get_example_data()
#' y <- list("GWAS" = GWAS_data(x[["GWAS"]]),
#'           "DE" = DE_data(x[["DE"]]),
#'           "CAN" = CAN_data(x[["CAN"]]))
#'
#' combine_chrom_length(y)
#' @export
combine_chrom_length <- function(x){

  chromosome <- NULL

  x |>
    purrr::map(compute_chrom_length) |>
    purrr::reduce(dplyr::bind_rows) |>
    dplyr::group_by(chromosome) |>
    dplyr::summarise(length = max(length),
                     .groups = "drop")

}
