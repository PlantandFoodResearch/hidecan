#' Creates a Manhattan plot
#'
#' Creates a Manhattan plot from a data-frame of GWAS results.
#'
#' @param gwas_list Data-frame or list of data-frames containing GWAS results,
#' each with at least a `chromosome`, `position` and either `padj` or `score`
#' columns. If a named list, the names will be used in the plot.
#' @param score_thr Numeric, the significance threshold on GWAS score. If a value
#' is provided, will be represented in the Manhattan plot as a red dashed line. If
#' `NULL` (default value), no significance threshold line will be drawn.
#' @param chroms Character vector, name of chromosomes to include in the plot.
#' If `NULL` (default value), all chromosomes will be included.
#' @param title Character, title of the plot. Default value is `NULL` (i.e.
#' no title will be added to the plot).
#' @param subtitle Character, subtitle of the plot. Default value is `NULL`
#' (i.e. no subtitle will be added to the plot).
#' @param size_range Numeric vector of length 2, the minimum and maximum point
#' size in the plot. Points size is proportional to their GWAS score. Default
#' value is `c(1, 3)`.
#' @param chrom_col Character vector of colour names or code, colours to use to draw
#' the points for each chromosome. Names will be ignored. If vector provided
#' contains less colours than the number of chromosomes to plot, the values will
#' be recycled. If `NULL`, default, colours will be automatically chosen from a
#' predefined palette.
#' @param ncol Integer, number of Manhattan plots per row when several GWAS results
#' are provided.
#' @returns A `ggplot`.
#' @examples
#' if (interactive()){
#' x <- get_example_data()[["GWAS"]]
#'
#' manhattan_plot(x)
#'
#' ## Adding a significance threshold line in the plot
#' manhattan_plot(x, score_thr = 4)
#'
#'
#' ## Use only two colours for the chromosomes
#' manhattan_plot(x,
#'                score_thr = 4,
#'                chrom_col = c("dodgerblue1", "dodgerblue4"))
#' }
#' @export
manhattan_plot <- function(gwas_list,
                           score_thr = NULL,
                           chroms = NULL,
                           title = NULL,
                           subtitle = NULL,
                           size_range = c(1, 3),
                           chrom_col = NULL,
                           ncol = NULL){

  ## For devtools::check
  facet <- NULL

  if(is.data.frame(gwas_list)){

    gwas_df <- GWAS_data(gwas_list)
    make_facets <- FALSE

  } else if(is.list(gwas_list)){

    if(is.null(names(gwas_list))) names(gwas_list) <- paste0("GWAS results ", seq_along(gwas_list))

    gwas_df <- gwas_list |>
      purrr::map_dfr(GWAS_data, .id = "facet") |>
      dplyr::mutate(facet = factor(facet, levels = names(gwas_list)))

    make_facets <- TRUE

  } else{
    stop("gwas_list argument should be either a list or a data-frame.")
  }



  chrom_length <- compute_chrom_length(gwas_df)

  ## for devtools::check
  chromosome <- offset <- position <- x_position <- score <- NULL

  ## Filtering chromosomes
  if(!is.null(chroms)){

    ## Checking that the chromosome names are valid
    bad_chroms <- setdiff(chroms, chrom_length$chromosome)
    if(length(bad_chroms)) stop("In chroms argument: '",
                                paste0(bad_chroms, collapse = "', '"),
                                "' are not valid chromosome names. Possible names are: '",
                                paste0(chrom_length$chromosome, collapse = "', '"),
                                "'.")

    ## Filters GWAS dataset and chrom length data
    gwas_df <- dplyr::filter(gwas_df, chromosome %in% chroms)
    chrom_length <- dplyr::filter(chrom_length, chromosome %in% chroms)

  } else {
    chroms <- chrom_length$chromosome
  }

  ## Selecting chromosome colours
  if(is.null(chrom_col)) {
    chrom_col <- RColorBrewer::brewer.pal(8, "Dark2")
  }

  chrom_col <- rep_len(chrom_col, length(chroms))
  names(chrom_col) <- chroms

  # Computing chromosomes x-axis position and middle (for x-axis labelling)
  chrom_length <- chrom_length |>
    dplyr::mutate(offset = c(0, cumsum(chrom_length$length)[-nrow(chrom_length)]),
                  middle = offset + ceiling(length / 2))

  ## Plot
  p <- gwas_df |>
    dplyr::left_join(chrom_length,
                     by = "chromosome") |>
    dplyr::mutate(x_position = position + offset) |>
    ggplot2::ggplot(ggplot2::aes(x = x_position,
                                 y = score,
                                 colour = chromosome,
                                 size = score)) +
    ggplot2::geom_point()

  ## Adding significance threshold horizontal line
  if(!is.null(score_thr)) {
    p <- p +
      ggplot2::geom_hline(yintercept = score_thr,
                          colour = "red",
                          linetype = 2)
  }

  if(make_facets){
    p <- p +
      ggplot2::facet_wrap(~facet, ncol = ncol)
  }

  p <- p +
    ggplot2::scale_colour_manual(values = chrom_col,
                                 guide = "none") +
    ggplot2::scale_size_continuous(range = size_range,
                                   guide = "none") +
    ggplot2::scale_x_continuous(breaks = chrom_length$middle,
                                labels = chrom_length$chromosome,
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL,
                  y = "GWAS score")


  return(p)

}

