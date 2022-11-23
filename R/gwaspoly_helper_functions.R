#' Extract information from GWASpoly output
#'
#' Extracts GWAS results and chomosome length from GWASpoly output.
#'
#' @param gwaspoly_input A `GWASpoly.fitted` or `GWASpoly.thresh` object
#' (returned by `GWASpoly::GWASpoly()` or `GWASpoly::set.threshold()` functions).
#' @param traits Character vector, traits for which GWAS results should be
#' extracted. If `NULL` (default value), all traits present are considered.
#' @param models Character vector, genetic models for which GWAS results should be
#' extracted. If `NULL` (default value), all genetic models present are considered.
#' @returns A list with the following elements:
#' * `gwas_data_list`: A named list of `GWAS_data` object, giving the markers score for
#' each possible trait/genetic model combination. The names of the list are in
#' the form `trait (genetic model)`.
#' * `gwas_data_thr_list`: if the input data is a `GWASpoly.thresh` object
#' (from the `GWASpoly::set.threshold()` function), a named list of `GwAS_data_thr`,
#' with the significant markers score for each possible trait/genetic model
#' combination. The names of the list are in the form `trait (genetic model)`.
#' * `chrom_length`: A tibble with one row per chromosome, giving the length
#' (in bp) of each chromosome.
#' @export
GWAS_data_from_gwaspoly <- function(gwaspoly_input,
                                    traits = NULL,
                                    models = NULL){

  if(!inherits(gwaspoly_input, "GWASpoly.fitted") | !inherits(gwaspoly_input, "GWASpoly.thresh")){
    stop("'gwaspoly_input' should be a `GWASpoly.fitted` or `GWASpoly.thresh` object (returned by GWASpoly() or set.threshold() functions).")
  }

  ## For devtools::check
  trait <- model <- Marker <- Chrom <- chromosome <- position <- NULL

  ## Checking which traits to use

  traits_list <- names(gwaspoly_input@scores)

  if(!is.null(traits)){
    bad_traits <- setdiff(traits, traits_list)
    if(length(bad_traits)) stop("The following traits are not present in the input data: '",
                                paste0(bad_traits, collapse = "', '"),
                                "'.\nPossible values for traits argument are: '",
                                paste0(traits_list, collapse = "', '"),
                                "'.")
  } else {
    traits <- traits_list
  }

  ## Checking which models to use

  models_list <- colnames(gwaspoly_input@scores[[1]])

  if(!is.null(models)){
    bad_models <- setdiff(models, models_list)
    if(length(bad_models)) stop("The following models are not present in the input data: '",
                                paste0(bad_models, collapse = "', '"),
                                "'.\nPossible values for models argument are: '",
                                paste0(models_list, collapse = "', '"),
                                "'.")
  } else {
    models <- models_list
  }

  ## Get all possible combinations of traits and models
  grid_df <- tidyr::expand_grid(trait = traits,
                                model = models) |>
    dplyr::mutate(label = paste0(trait, " (", model, ")"))

  ## Get markers info
  markers_info <-   gwaspoly_input@map |>
    tibble::as_tibble() |>
    dplyr::select(marker = Marker,
                  chromosome = Chrom,
                  position = Position)

  ## Get chromosome length (need to do it only once)
  chrom_length <- markers_info |>
    dplyr::group_by(chromosome) |>
    dplyr::summarise(length = max(position),
                     .groups = "drop")

  ## Get the list of marker scores for each trait and model
  gwas_data_list <- purrr::map2(grid_df$trait,
                                grid_df$model,
                                ~ gwaspoly_input@scores[[.x]][, .y, drop = FALSE] |>
                                  tibble::as_tibble(rownames = "marker") |>
                                  dplyr::rename(score = !!dplyr::sym(.y))) |>
    purrr::map(~ markers_info |>
                 dplyr::full_join(.x, by = "marker") |>
                 dplyr::filter(!is.na(score))) |>
    purrr::map(GWAS_data)

    names(gwas_data_list) <- grid_df$label

  ## If possible, extract the threshold values
  if(inherits(gwaspoly_input, "GWASpoly.thresh")){

    thr_df <- gwaspoly_input@threshold |>
      tibble::as_tibble(rownames = "trait") |>
      tidyr::pivot_longer(cols = -trait,
                          names_to = "model",
                          values_to = "thr") |>
      dplyr::full_join(grid_df, by = c("trait", "model"))

    gwas_data_thr_list <- purrr::map2(
      thr_df$label,
      thr_df$thr,
      ~ apply_threshold(gwas_data_list[[.x]], score_thr = .y)
    )

    names(gwas_data_thr_list) <- grid_df$label

  } else {
    gwas_data_thr_list <- NULL
  }

  return(list(gwas_data_list = gwas_data_list,
              gwas_data_thr_list = gwas_data_thr_list,
              chrom_length = chrom_length))
}

#' Create a hidecan plot from a GWASpoly output
#'
#' Creates a hidecan plot from the GWAS results from GWASpoly.
#'
#' @param gwaspoly_input A `GWASpoly.thresh` object
#' (returned by the `GWASpoly::set.threshold()` function).
#' @param ... Further arguments passed to the \code{\link{create_hidecan_plot}()}
#' function.
#' @returns A ggplot.
#' @examples
#' \dontrun{
#' x <- get_gwaspoly_example_data(with_thresholds = TRUE)
#'
#' hidecan_plot_from_gwaspoly(x)
#' }
#' @export
hidecan_plot_from_gwaspoly <- function(gwaspoly_input, ...){

  if(!inherits(gwaspoly_input, "GWASpoly.fitted") | !inherits(gwaspoly_input, "GWASpoly.thresh")){
    stop("'gwaspoly_input' should be a `GWASpoly.thresh` object (returned by the set.threshold() function).")
  }

  input_list <- GWAS_data_from_gwaspoly(gwaspoly_input)

  create_hidecan_plot(
    input_list$gwas_data_thr_list,
    chrom_length = input_list$chrom_length,
    ...
  )
}
