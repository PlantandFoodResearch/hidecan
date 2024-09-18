#' Creates a HIDECAN plot
#'
#' Creates a HIDECAN plot from a list of filtered GWAS or DE results
#' and/or candidate genes.
#'
#' @param x A list of `GWAS_data_thr`, `DE_data_thr`, `CAN_data_thr` and/or
#'   `CUSTOM_data_thr` produced by the \code{\link{apply_threshold}()} function.
#'   If named, the names will be appended to the y-axis labels (use `' '` as
#'   empty name in the list).
#' @param chrom_length Tibble with columns `chromosome` and `length`, giving
#' for each chromosome its length in bp (see \code{\link{combine_chrom_length}()}
#' function).
#' @param remove_empty_chrom Logical, should chromosomes with no significant
#' markers/genes nor candidate genes be removed from the plot? Default value
#' if `FALSE`.
#' @param chroms Character vector, name of chromosomes to include in the plot.
#' If `NULL` (default value), all chromosomes will be included.
#' @param chrom_limits Integer vector of length 2, or named list where the
#' elements are integer vectors of length 2. If vector, gives the lower and upper
#' limit of the chromosomes (in bp) to use in the plot. If a named list, names
#' should correspond to chromosome names. Gives for each chromosome the lower
#' and upper limits (in bp) to use in the plot. Doesn't have to be specified
#' for all chromosomes. Default value is `NULL`, i.e. no limits are applied
#' to the chromosomes (they will be plotted in their entirety).
#' @param colour_genes_by_score Logical, whether to colour the genes by score
#' (`TRUE`) or by log2(fold-change) (`FALSE`). Default value is `TRUE`.
#' @param title Character, title of the plot. Default value is `NULL` (i.e.
#' no title will be added to the plot).
#' @param subtitle Character, subtitle of the plot. Default value is `NULL`
#' (i.e. no subtitle will be added to the plot).
#' @param n_rows Integer, number of rows of chromosomes to create in the plot.
#' Default value is `NULL`.
#' @param n_cols Integer, number of columns of chromosomes to create in the plot.
#' Default value is 2. Will be set to `NULL` if `n_rows` is not `NULL`.
#' @param legend_position Character, position of the legend in the plot. Can be
#' `bottom` (default value), `top`, `right`, `left` or `none`.
#' @param point_size Numeric, size of the points in the plot. Default value is 3.
#' @param label_size Numeric, size of the gene labels in the plot. Default value is
#' 3.5 (for \code{\link[ggrepel]{geom_label_repel}}).
#' @param label_padding Numeric, amount of padding around gene labels in the plot,
#' as unit or number. Default value is 0.15
#' (for \link[ggrepel]{geom_label_repel}).
#' @param custom_aes Named list of plot aesthetics for custom data types. See
#'   [hidecan_aes()] for information about the content of each element. Default
#'   is `NULL` (only needed if there are `CUSTOM_data_thr` objects in `x` or to
#'   customise the aesthetics for the default data tracks).
#' @returns A `ggplot`.
#' @examples
#' if (interactive()) {
#' x <- get_example_data()
#' y <- list("GWAS" = GWAS_data(x[["GWAS"]]),
#'           "DE" = DE_data(x[["DE"]]),
#'           "CAN" = CAN_data(x[["CAN"]]))
#'
#' chrom_length <- combine_chrom_length(y)
#'
#' z <- list(
#'   apply_threshold(y[["GWAS"]], score_thr = 4),
#'   apply_threshold(y[["DE"]], score_thr = 1.3, log2fc_thr = 0.5),
#'   apply_threshold(y[["CAN"]])
#' )
#'
#' create_hidecan_plot(z,
#'                     chrom_length,
#'                     label_size = 2)
#'
#' ## Colour genes according to their fold-change
#' create_hidecan_plot(z,
#'                     chrom_length,
#'                     colour_genes_by_score = FALSE,
#'                     label_size = 2)
#'
#' ## Add names to the datasets
#' create_hidecan_plot(setNames(z, c("Genomics", "RNAseq", "My list")),
#'                     chrom_length,
#'                     colour_genes_by_score = FALSE,
#'                     label_size = 2)
#'
#' ## Add names to some of the datasets only (e.g. not for GWAS results)
#' create_hidecan_plot(setNames(z, c(" ", "RNAseq", "My list")),
#'                     chrom_length,
#'                     colour_genes_by_score = FALSE,
#'                     label_size = 2)
#'
#' ## Set limits on all chromosomes (to "zoom in" to the 10-20Mb region)
#' create_hidecan_plot(z,
#'                     chrom_length,
#'                     label_size = 2,
#'                     chrom_limits = c(10e6, 20e6))
#'
#' ## Set limits on some chromosomes only
#' create_hidecan_plot(z,
#'                     chrom_length,
#'                     label_size = 2,
#'                     chrom_limits = list("ST4.03ch00" = c(10e6, 20e6),
#'                                         "ST4.03ch02" = c(15e6, 25e6)))
#' }
#' @export
create_hidecan_plot <- function(x,
                                chrom_length,
                                colour_genes_by_score = TRUE,
                                remove_empty_chrom = FALSE,
                                chroms = NULL,
                                chrom_limits = NULL,
                                title = NULL,
                                subtitle = NULL,
                                n_rows = NULL,
                                n_cols = 2,
                                legend_position = "bottom",
                                point_size = 3,
                                label_size = 3.5,
                                label_padding = 0.15,
                                custom_aes = NULL){

  ## for devtools::check()
  position <- dataset <- score <- chromosome <- position_mb <- NULL
  position_mb_end <- data_type <- name <- log2FoldChange <- NULL
  upper_limit_mb <- lower_limit_mb <- NULL

  ## Checking x input
  x_types <- purrr::map_chr(x, \(.x) class(.x)[[1]])

  accepted_types <- c("GWAS_data_thr", "DE_data_thr",
                      "CAN_data_thr", "CUSTOM_data_thr")
  if (!all(unique(x_types) %in% accepted_types)) {
    stop("Expecting a list of 'GWAS_data_thr', 'DE_data_thr', 'CAN_data_thr' ",
         "and/or 'CUSTOM_data_thr' objects (see apply_threshold() function).")
  }

  ## Checking chrom_length input
  if (length(setdiff(c("chromosome", "length"), names(chrom_length))) > 0) {
    stop("'chrom_length' argument should be a data-frame with columns 'chromosome' and 'length'.")
  }

  if (dplyr::n_distinct(chrom_length$chromosome) != nrow(chrom_length)) {
    stop("Duplicated chromosome names in 'chrom_length' data-frame.")
  }

  aes_types <- purrr::map_chr(x, .get_aes_type)

  ## Get aesthetics for the plots (line colours, point shape, etc) for each
  ## data type
  plot_aes <- .get_plot_aes(aes_types, colour_genes_by_score, custom_aes)

  data_type_labels <- purrr::map_chr(plot_aes, \(x) purrr::pluck(x, "y_label"))
  data_type_colours <- purrr::map_chr(plot_aes, \(x) purrr::pluck(x, "line_colour"))
  data_type_shapes <- purrr::map_int(plot_aes, \(x) purrr::pluck(x, "point_shape"))

  ## Making y labels
  if (is.null(names(x))) {
    names(x) <- make.unique(data_type_labels[aes_types], " - ")
  } else {
    names(x) <- paste0(names(x), " - ", data_type_labels[aes_types])
    names(x) <- sub("^ +- ", "", names(x))
    names(x) <- make.unique(names(x), ", ")
  }

  datasets_levels <- rev(names(x))

  ## Making sure that only one of n_rows and n_cols is not NULL
  if (!is.null(n_rows)) n_cols <- NULL

  ## Constructing main tibble for plot
  toplot <- x |>
    purrr::map(\(.x) dplyr::mutate(.x, data_type = .get_aes_type(.x))) |>
    purrr::list_rbind(names_to = "dataset") |>
    dplyr::mutate(
      position_mb = position / 1e6,
      dataset = factor(dataset, levels = datasets_levels)
    )

  if (!colour_genes_by_score) {
    toplot <- toplot |>
      dplyr::mutate(
        score = dplyr::case_when(
          data_type == "DE_data_thr" ~ log2FoldChange,
          TRUE ~ score
        )
      )
  }

  if ('score' %in% names(toplot)) toplot <- toplot |> dplyr::arrange(score)

  chroms <- .check_chroms(toplot, chrom_length, chroms, remove_empty_chrom)

  ## Now filter both toplot and chrom_lenght to only the desired chromosomes
  toplot <- toplot |>
    dplyr::filter(chromosome %in% chroms) |>
    dplyr::mutate(chromosome = factor(chromosome, levels = chroms, ordered = FALSE))

  chrom_length <- chrom_length |>
    dplyr::filter(chromosome %in% chroms) |>
    dplyr::mutate(chromosome = factor(chromosome, levels = chroms, ordered = FALSE))

  toplot_chroms <- chrom_length |>
    dplyr::mutate(position_mb = length / 1e6) |>
    dplyr::select(chromosome, position_mb) |>
    tidyr::expand_grid(dataset = factor(names(x), levels = datasets_levels)) |>
    dplyr::rename(position_mb_end = position_mb) |>
    dplyr::mutate(position_mb = 0)

  ## Handling chromosome limits
  if (!is.null(chrom_limits)) {
    chrom_limits_df <- .check_chrom_limits(chrom_limits, chrom_length, chroms)

    toplot <- toplot |>
      dplyr::left_join(chrom_limits_df, by = "chromosome") |>
      dplyr::filter(position_mb >= lower_limit_mb, position_mb <= upper_limit_mb) |>
      dplyr::select(-lower_limit_mb, -upper_limit_mb)

    toplot_chroms <- toplot_chroms |>
      dplyr::left_join(chrom_limits_df, by = "chromosome") |>
      dplyr::mutate(position_mb = lower_limit_mb,
                    position_mb_end = upper_limit_mb) |>
      dplyr::select(-lower_limit_mb, -upper_limit_mb)
  }

  ## Making sure that the order of the "Position of" legends matches the order
  ## in which the different data types appear in the y-axis
  unique_data_types <- unique(aes_types)

  p <- toplot |>
    ggplot2::ggplot(ggplot2::aes(x = position_mb, y = dataset)) +
    ## Chromosome segments
    ggplot2::geom_segment(data = toplot_chroms,
                          ggplot2::aes(xend = position_mb_end, yend = dataset)) +
    ggplot2::facet_wrap(~ chromosome, scales = "free_x", nrow = n_rows, ncol = n_cols) +
    ## vertical position indicators
    ggplot2::geom_vline(ggplot2::aes(xintercept = position_mb, colour = data_type),
                        alpha = 0.7) +
    ## General colours and shapes
    ggplot2::scale_colour_manual(values = data_type_colours[unique_data_types],
                                 labels = data_type_labels[unique_data_types],
                                 breaks = unique_data_types) +
    ggplot2::scale_shape_manual(values = data_type_shapes[unique_data_types],
                                labels = data_type_labels[unique_data_types],
                                breaks = unique_data_types) +
    ggplot2::guides(colour = ggplot2::guide_legend(title.position = "top",
                                                   title.hjust = 0.5,
                                                   override.aes = list(alpha = 1),
                                                   order = 1),
                    shape = ggplot2::guide_legend(order = 1))

  ## Check whether we need to add a new fill legend for each data type
  new_legend <- !(unique_data_types %in% c(unique_data_types[[1]], "CAN_data_thr"))
  names(new_legend) <- unique_data_types

  for (i in unique_data_types) {
    p <- .add_data_type(
      p,
      toplot,
      i,
      paes = plot_aes[[i]],
      add_new_legend = new_legend[[i]],
      point_size,
      label_size,
      label_padding
    )
  }

  p +
    ## Themes and labs
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = legend_position,
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(title = title,
                  subtitle = subtitle,
                  x = "Position (Mb)",
                  y = NULL,
                  colour = "Position of",
                  shape = "Position of")

}

#' Wrapper to create a HIDECAN plot
#'
#' Wrapper function to create a HIDECAN plot from GWAS results, DE results or
#' candidate genes.
#'
#' @param gwas_list Data-frame or list of data-frames containing GWAS results,
#'   each with at least a `chromosome`, `position` and either `padj` or `score`
#'   columns. If a named list, the names will be used in the plot.
#' @param de_list Data-frame or list of data-frames containing DE results, each
#'   with at least a `chromosome`, `start`, `end`, `log2FoldChange` and either
#'   `padj` or `score` columns. If a named list, the names will be used in the
#'   plot.
#' @param can_list Data-frame or list of data-frames containing candidate genes,
#'   each with at least a `chromosome`, `start`, `end` and `name` columns. If a
#'   named list, the names will be used in the plot.
#' @param score_thr_gwas Numeric, the score threshold for GWAS results that will
#'   be used to select which markers will be plotted. Default value is 4.
#' @param score_thr_de Numeric, the score threshold for DE results that will be
#'   used to select which markers will be plotted. Default value is 2.
#' @param log2fc_thr Numeric, the log2(fold-change) threshold that will be used
#'   to select which genes will be plotted. Default value is 1.
#' @param chrom_length Optional, tibble with columns `chromosome` and `length`,
#'   giving for each chromosome its length in bp. If `NULL` (the default), will
#'   be inferred from the GWAS, DE and candidate gene data.
#' @inheritParams create_hidecan_plot
#' @param custom_list Data-frame or list of data-frames containing custom
#'   genomic features, each with at least a `chromosome`, `position` and `score`
#'   columns. If a named list, the names will be used in the plot.
#' @param score_thr_custom Numeric, the score threshold for custom genomic
#'   features that will be used to select which markers will be plotted. Default
#'   value is 0.
#' @returns a `ggplot`.
#' @examples
#' if (interactive()) {
#' x <- get_example_data()
#'
#' ## Typical example with one GWAs result table, one DE result table and
#' ## one table of candidate genes
#' hidecan_plot(gwas_list = x[["GWAS"]],
#'              de_list = x[["DE"]],
#'              can_list = x[["CAN"]],
#'              score_thr_gwas = -log10(0.0001),
#'              score_thr_de = -log10(0.005),
#'              log2fc_thr = 0,
#'              label_size = 2)
#'
#' ## Example with two sets of GWAS results
#' hidecan_plot(gwas_list = list(x[["GWAS"]], x[["GWAS"]]),
#'              score_thr_gwas = 4)
#'
#' ## Example with two sets of DE results, with names
#' hidecan_plot(de_list = list("X vs Y" = x[["DE"]],
#'                             "X vs Z" = x[["DE"]]),
#'              score_thr_de = -log10(0.05),
#'              log2fc_thr = 0)
#'
#' ## Set limits on all chromosomes (to "zoom in" to the 10-20Mb region)
#' hidecan_plot(gwas_list = x[["GWAS"]],
#'              de_list = x[["DE"]],
#'              can_list = x[["CAN"]],
#'              score_thr_gwas = -log10(0.0001),
#'              score_thr_de = -log10(0.005),
#'              log2fc_thr = 0,
#'              label_size = 2,
#'              chrom_limits = c(10e6, 20e6))
#'
#' ## Set limits on some chromosomes only
#' hidecan_plot(gwas_list = x[["GWAS"]],
#'              de_list = x[["DE"]],
#'              can_list = x[["CAN"]],
#'              score_thr_gwas = -log10(0.0001),
#'              score_thr_de = -log10(0.005),
#'              log2fc_thr = 0,
#'              label_size = 2,
#'              chrom_limits = list("ST4.03ch00" = c(10e6, 20e6),
#'                                   "ST4.03ch02" = c(15e6, 25e6)))
#' }
#' @export
hidecan_plot <- function(gwas_list = NULL,
                         de_list = NULL,
                         can_list = NULL,
                         score_thr_gwas = 4,
                         score_thr_de = 2,
                         log2fc_thr = 1,
                         chrom_length = NULL,
                         colour_genes_by_score = TRUE,
                         remove_empty_chrom = FALSE,
                         chroms = NULL,
                         chrom_limits = NULL,
                         title = NULL,
                         subtitle = NULL,
                         n_rows = NULL,
                         n_cols = 2,
                         legend_position = "bottom",
                         point_size = 3,
                         label_size = 3.5,
                         label_padding = 0.15,
                         custom_aes = NULL,
                         custom_list = NULL,
                         score_thr_custom = 0
){

  ## Input should either be NULL, a data-frame or a list
  list_or_null <- list(gwas_list, de_list, can_list) |>
    purrr::map_lgl(\(.x) (is.list(.x) | is.null(.x)))

  if (!all(list_or_null)) {
    stop("Arguments 'gwas_list', 'de_list' or 'can_list' should either be a data-frame or a list, or NULL.")
  }

  ## Separate error message for custom input because don't want to confuse
  ## people using it with only intended data types
  if (!(is.list(custom_list) | is.null(custom_list))) {
    stop("Argument 'custom_list' should either be a data-frame or a list, or NULL.")
  }

  ## If providing a single data-frame, turn into a list
  if (is.data.frame(gwas_list)) gwas_list <- list(gwas_list)
  if (is.data.frame(de_list)) de_list <- list(de_list)
  if (is.data.frame(can_list)) can_list <- list(can_list)
  if (is.data.frame(custom_list)) custom_list <- list(custom_list)

  ## create the input tibbles from ordinary matrices or data-frames
  error_func <- function(arg){
    res <- substitute(function(err){
      msg <- conditionMessage(err)
      stop("In '", arg, "' argument: ", msg, call. = FALSE)
    })

    return(res)
  }

  gwas_list <- gwas_list |>
    purrr::map(~ tryCatch(GWAS_data(.x), error = eval(error_func("gwas_list"))))

  de_list <- de_list |>
    purrr::map(~ tryCatch(DE_data(.x), error = eval(error_func("de_list"))))

  can_list <- can_list |>
    purrr::map(~ tryCatch(CAN_data(.x), error = eval(error_func("can_list"))))

  custom_list <- custom_list |>
    purrr::map(~ tryCatch(CUSTOM_data(.x), error = eval(error_func("custom_list"))))

  x <- c(gwas_list, de_list, can_list, custom_list)

  if(length(x) == 0) {
    stop("Should provide at least one non-empty list for 'gwas_list', 'de_list',",
         " 'can_list' or 'custom_list' argument.")
  }

  ## if not provided, compute the chromosome length
  if(is.null(chrom_length)){
    chrom_length <- combine_chrom_length(x)
  }

  ## Check threshold values
  if(length(score_thr_gwas) != 1) stop("'score_thr_gwas' argument should be a single numeric value.")
  if(!is.numeric(score_thr_gwas)) stop("'score_thr_gwas' argument should be a numeric value.")

  if(length(score_thr_de) != 1) stop("'score_thr_de' argument should be a single numeric value.")
  if(!is.numeric(score_thr_de)) stop("'score_thr_de' argument should be a numeric value.")

  if(length(log2fc_thr) != 1) stop("'log2fc_thr' argument should be a single numeric value.")
  if(!is.numeric(log2fc_thr)) stop("'log2fc_thr' argument should be a numeric value.")

  if(length(score_thr_custom) != 1) stop("'score_thr_custom' argument should be a single numeric value.")
  if(!is.numeric(score_thr_custom)) stop("'score_thr_custom' argument should be a numeric value.")

  ## Apply threshold to datasets
  score_thr <- c(
    "GWAS_data" = score_thr_gwas,
    "DE_data" = score_thr_de,
    "CAN_data" = 0,
    "CUSTOM_data" = score_thr_custom
  )

  x <- purrr::map2(
    x,
    purrr::map_chr(x, ~ class(.x)[[1]]),
    ~ apply_threshold(.x, score_thr = score_thr[[.y]], log2fc_thr = log2fc_thr)
  )

  ## create the plot
  create_hidecan_plot(x,
                      chrom_length,
                      colour_genes_by_score,
                      remove_empty_chrom,
                      chroms,
                      chrom_limits,
                      title,
                      subtitle,
                      n_rows,
                      n_cols,
                      legend_position,
                      point_size,
                      label_size,
                      label_padding,
                      custom_aes)
}


#' Default aesthetics for HIDECAN plot
#'
#' Generates a list of the default aesthetics used for HIDECAN plots.
#'
#' @param colour_genes_by_score Logical, whether to colour the genes by score
#' (`TRUE`) or by log2(fold-change) (`FALSE`). Default value is `TRUE`.
#' @returns A named list, with one element per type of data (e.g. GWAS, DE, etc).
#' Each element is itself a list with the following elements:
#' * `y_label`: prefix added to the name of a track on the y-axis of the plot.
#' * `line_colour`: colour of the vertical line used to show the position of
#'  elements of this type.
#' * `point_shape`: shape used for the points of this type.
#' * `show_name`: whether a label with `name` value should be added to points of
#'  this type.
#' * `fill_scale`: fill scale to use for the points of this type.
#' @export
hidecan_aes <- function(colour_genes_by_score = TRUE) {
  if (colour_genes_by_score) {
    de_fill_scale <- viridis::scale_fill_viridis(
      "DE gene score",
      option = "viridis",
      guide = ggplot2::guide_colourbar(title.position = "top",
                                       title.hjust = 0.5,
                                       order = 3)
    )
  } else {
    de_fill_scale <- ggplot2::scale_fill_gradient2(
      "DE gene log2(fold-change)",
      low = "darkblue",
      mid = "white",
      high = "firebrick",
      guide = ggplot2::guide_colourbar(title.position = "top",
                                       title.hjust = 0.5,
                                       order = 3)
    )
  }

  list(
    "GWAS_data_thr" = list(
      "y_label" = "GWAS peaks",
      "line_colour" = "red",
      "point_shape" = 21,
      "show_name" = FALSE,
      "fill_scale" = viridis::scale_fill_viridis(
        "GWAS marker score",
        option = "plasma",
        guide = ggplot2::guide_colourbar(title.position = "top",
                                         title.hjust = 0.5,
                                         order = 2))
    ),

    "DE_data_thr" = list(
      "y_label" = "DE genes",
      "line_colour" = "darkcyan",
      "point_shape" = 23,
      "show_name" = FALSE,
      "fill_scale" = de_fill_scale
    ),

    "CAN_data_thr" = list(
      "y_label" = "Candidate genes",
      "line_colour" = "grey30",
      "point_shape" = 4,
      "show_name" = TRUE,
      "fill_scale" = NULL
    ),


    "CUSTOM_data_thr" = list(
      "y_label" = "Custom track",
      "line_colour" = "darkgoldenrod2",
      "point_shape" = 21,
      "show_name" = FALSE,
      "fill_scale" = viridis::scale_fill_viridis(
        "Custom score",
        option = "rocket",
        guide = ggplot2::guide_colourbar(title.position = "top",
                                         title.hjust = 0.5,
                                         order = 4))
    )
  )
}

#' Returns either "aes_type" attribute or object class
#'
#' @param x A `GWAS_data_thr`, `DE_data_thr`, `CAN_data_thr` or
#'   `CUSTOM_data_thr` object.
#' @returns Either "aes_type" attribute or object class.
.get_aes_type <- function(x) {
  res <- attr(x, "aes_type")
  ifelse(
    is.null(res),
    class(x)[[1]],
    res
  )
}

#' Computes plot aesthetics
#'
#' Get the list of aesthetics for each data types in the plot.
#'
#' @param aes_types Character vector of data types (one per dataset to plot),
#' should be computed with [.get_aes_type()].
#' @param colour_genes_by_score Logical, whether to colour the genes by score
#' (`TRUE`) or by log2(fold-change) (`FALSE`).
#' @param custom_aes Named list of plot aesthetics for custom data types. See
#' [hidecan_aes()] for information about the content of each element.
#' @returns Named list of plot aesthetics.
.get_plot_aes <- function(aes_types, colour_genes_by_score, custom_aes) {

  default_aes <- hidecan_aes(colour_genes_by_score)

  if (!is.null(custom_aes)) {
    default_aes <- c(
      custom_aes,
      default_aes[setdiff(names(default_aes), names(custom_aes))]
    )
  }

  missing_aes <- setdiff(aes_types, names(default_aes))

  if (length(missing_aes) > 0) {
    stop("'custom_aes' function should contain the following elements: '",
         paste0(missing_aes, collapse = "', "), "'.",
         call. = FALSE)
  }

  return(default_aes)
}

#' Check chromosomes to plot
#'
#' @param toplot Tibble, main data-frame for plot.
#' @param chrom_length Tibble of chromosome length.
#' @param chroms Character vector, name of chromosomes to include in the plot.
#' If `NULL` (default value), all chromosomes will be included.
#' @param remove_empty_chrom Logical, should chromosomes with no significant
#' markers/genes nor candidate genes be removed from the plot? Default value
#' if `FALSE`.
#' @returns Character vector of chromosomes to plot.
.check_chroms <- function(toplot, chrom_length, chroms, remove_empty_chrom) {
  ## chromosomes present in dataset and chrom_length
  cl_chroms <- unique(chrom_length$chromosome)
  tp_chroms <- unique(toplot$chromosome)
  all_chrom <- union(cl_chroms, tp_chroms)

  if (is.null(chroms)) {
    chroms <- all_chrom
  } else {
    wrong_chroms <- setdiff(chroms, all_chrom)
    if (length(wrong_chroms)) {
      stop("In 'chroms' argument, the following values are not valid chromosome names: '",
           paste0(wrong_chroms, collapse = "', '"),
           ". Possible values are: '",
           paste0(all_chrom, collapse = "', '"),
           ".")
    }
  }

  ## If want to remove empty chromosomes, look only at the ones that are in toplot
  if(remove_empty_chrom) chroms <- intersect(chroms, tp_chroms)

  ## checking whether some chromosomes in toplot are missing from chrom_length
  ## check happens after selecting the specified chromosomes (so then if we don't have info about a chromosome
  ## but don't want to plot it it's fine), and also after having removed the empty chromosomes if needed
  ## for the same reason
  missing_chroms <- setdiff(chroms, cl_chroms)
  if (length(missing_chroms)) {
    stop("The following chromosomes are present in the data but missing from 'chrom_length' data-frame: '",
         paste0(missing_chroms, collapse = "', '"), ".")
  }

  chroms
}

#' Check chromosome limits
#'
#' @param chrom_limits Integer vector of length 2, or named list where the
#'   elements are integer vectors of length 2. If vector, gives the lower and
#'   upper limit of the chromosomes (in bp) to use in the plot. If a named list,
#'   names should correspond to chromosome names. Gives for each chromosome the
#'   lower and upper limits (in bp) to use in the plot. Doesn't have to be
#'   specified for all chromosomes. Default value is `NULL`, i.e. no limits are
#'   applied to the chromosomes (they will be plotted in their entirety).
#' @param chrom_length Tibble of chromosome length.
#' @param chroms Character vector, name of chromosomes to include in the plot.
#' @returns Tibble of chromosome limits to use in the plot.
.check_chrom_limits <- function(chrom_limits, chrom_length, chroms) {
  ## for devtools::check()
  chromosome <- upper_limit_mb <- lower_limit_mb <- NULL

  if (is.numeric(chrom_limits) & length(chrom_limits) == 2) {

    ## Apply the limits to all chromosomes
    chrom_limits <- purrr::map(chroms, ~ chrom_limits) |>
      stats::setNames(chroms)

  } else if (is.list(chrom_limits)) {

    if (is.null(names(chrom_limits))) {
      stop("The chrom_limits argument should be a named list, with the names corresponding to chromosomes' name.")
    }
    if (!all(purrr::map_lgl(chrom_limits, is.numeric)) | !all(purrr::map_lgl(chrom_limits, ~length(.x) == 2))) {
      stop("The chrom_limits argument should be a named list where each element is an integer vector of length 2.")
    }

    bad_chroms <- setdiff(names(chrom_limits), chroms)
    if (length(bad_chroms)) {
      stop("In chrom_limits argument: '",
           paste0(bad_chroms, collapse = "', '"),
           "' are not valid chromosome names. Possible names are: '",
           paste0(chroms, collapse = "', '"),
           "'.")
    }
    chrom_limits <- chrom_limits[intersect(names(chrom_limits), chroms)]

  } else {
    stop("The chrom_limits argument should either be an integer vector of length 2 or a named list where each element is an integer vector of length 2.")
  }

  chrom_length |>
    dplyr::left_join(
      purrr::map_dfr(
        chrom_limits,
        ~ tibble::tibble(lower_limit_mb = .x[[1]],
                         upper_limit_mb = .x[[2]]),
        .id = "chromosome"
      ) |>
        dplyr::mutate(chromosome = factor(chromosome, levels = chroms)),
      by = "chromosome"
    ) |>
    tidyr::replace_na(list(lower_limit_mb = 0)) |>
    dplyr::mutate(upper_limit_mb = dplyr::coalesce(upper_limit_mb, length),
                  upper_limit_mb = purrr::map2_dbl(length, upper_limit_mb, ~ min(c(.x, .y))),
                  lower_limit_mb = lower_limit_mb / 1e6,
                  upper_limit_mb = upper_limit_mb / 1e6) |>
    dplyr::select(-length)
}

#' Add data type to HIDECAN plot
#'
#' @param p Current ggplot.
#' @param toplot Tibble of data to plot.
#' @param i Character, data type to add in the plot.
#' @param paes Named list, aesthetics to use for this data type (see
#'   [hidecan_aes()] for information about the required values).
#' @param add_new_legend Logical, whether the current data type requires a new
#'   legend for fill.
#' @param point_size Numeric, size of the points in the plot.
#' @param label_size Numeric, size of the gene labels in the plot.
#' @param label_padding Numeric, amount of padding around gene labels in the
#'   plot, as unit or number.
#' @returns Ggplot `p` to which the new data type has been added.
.add_data_type <- function(p,
                           toplot,
                           i,
                           paes,
                           add_new_legend,
                           point_size,
                           label_size,
                           label_padding) {
  ## for devtools::check()
  data_type <- name <- score <- NULL

  sub_data <- dplyr::filter(toplot, data_type == i)

  if (add_new_legend) p <- p + ggnewscale::new_scale_fill()

  if (paes$show_name) {
    p <- p +
      ggrepel::geom_label_repel(ggplot2::aes(label = name),
                                data = sub_data,
                                nudge_y = 0.5,
                                na.rm = TRUE,
                                size = label_size,
                                label.padding = label_padding,
                                alpha = 0.5)
  }

  if (is.null(paes$fill_scale)) {
    scatter_mapping <- ggplot2::aes(shape = data_type)
  } else {
    scatter_mapping <- ggplot2::aes(shape = data_type, fill = score)
  }

  p <- p +
    ggplot2::geom_point(data = sub_data,
                        mapping = scatter_mapping,
                        size = point_size) +
    paes$fill_scale
}

# create_hidecan_plot_interactive <- function(x,
#                                             chrom_length,
#                                             colour_genes_by_score = TRUE,
#                                             title = NULL,
#                                             subtitle = NULL,
#                                             n_rows = NULL,
#                                             n_cols = 2,
#                                             legend_position = "bottom",
#                                             point_size = 3,
#                                             label_size = 3.5,
#                                             label_padding = 0.15){
#
#   position <- dataset <- score <- chromosome <- position_mb <- position_mb_end <- data_type <- name <- log2FoldChange <- NULL
#
#   ## Labels, colours and shapes
#   data_type_labels <- c("GWAS_data_thr" = "GWAS peaks",
#                         "DE_data_thr" = "DE genes",
#                         "CAN_data_thr" = "Candidate genes")
#
#   data_type_colours <- c("GWAS_data_thr" = "red",
#                          "DE_data_thr" = "darkcyan",
#                          "CAN_data_thr" = "grey30")
#
#   data_type_shapes <- c("GWAS_data_thr" = 21,
#                         "DE_data_thr" = 23,
#                         "CAN_data_thr" = 4)
#
#
#   x_types <- purrr::map_chr(x, ~ class(.x)[[1]])
#
#   if(!all(unique(x_types) %in% c("GWAS_data_thr", "DE_data_thr", "CAN_data_thr"))) stop("Expecting a list of 'GWAS_data_thr', 'DE_data_thr' and/or 'CAN_data_thr' objects (see apply_threshold() function).")
#
#   if(is.null(names(x))){
#     names(x) <- make.unique(data_type_labels[x_types], " - ")
#   } else {
#     names(x) <- paste0(names(x), " - ", data_type_labels[x_types])
#     names(x) <- sub("^  - ", "", names(x))
#   }
#
#   datasets_levels <- rev(names(x))
#
#   toplot <- x |>
#     purrr::map_dfr(~ dplyr::mutate(.x, data_type = class(.x)[[1]]),
#                    .id = "dataset") |>
#     dplyr::mutate(position_mb = position / 1e6,
#                   dataset = factor(dataset, levels = datasets_levels)) |>
#     dplyr::arrange(score)
#
#   toplot_chroms <- chrom_length |>
#     dplyr::mutate(position_mb = length / 1e6) |>
#     dplyr::select(chromosome, position_mb) |>
#     tidyr::expand_grid(dataset = factor(names(x), levels = datasets_levels)) |>
#     dplyr::rename(position_mb_end = position_mb) |>
#     dplyr::mutate(position_mb = 0)
#
#   p <- toplot |>
#     ggplot2::ggplot(ggplot2::aes(x = position_mb, y = dataset)) +
#     ggplot2::facet_wrap(~ chromosome, scales = "free_x", nrow = n_rows, ncol = n_cols) +
#     ## General colours and shapes
#     ggplot2::scale_shape_manual(values = data_type_shapes, labels = data_type_labels, guide = "none") +
#     ggplot2::geom_point(ggplot2::aes(shape = data_type, fill = score),
#                         size = point_size) +
#     ggplot2::geom_label(ggplot2::aes(label = name), na.rm = TRUE) +
#     ## Themes and labs
#     viridis::scale_fill_viridis(option = "plasma") +
#     ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
#     ggplot2::theme_bw() +
#     ggplot2::theme(legend.position = legend_position) +
#     ggplot2::labs(title = title,
#                   subtitle = subtitle,
#                   x = "Position (Mb)",
#                   y = NULL,
#                   colour = "Position of",
#                   shape = "Position of")
#
#   ggplotly(p)
#
# }
#
# chrom_midpoints <- chrom_length |>
#   dplyr::mutate(midpoint = length / 2e6) |>
#   dplyr::select(chromosome, midpoint) |>
#   tibble::deframe()
#
# toplot |>
#   dplyr::group_by(chromosome) |>
#   tidyr::nest() |>
#   dplyr::left_join(toplot_chroms |>
#                      group_by(chromosome) |>
#                      tidyr::nest() |>
#                      dplyr::rename(data_chrom = data),
#                    by = "chromosome") |>
#   dplyr::mutate(subplot = purrr::pmap(list(chromosome, data, data_chrom),
#                                       function(chrom, data, data_chrom){
#                                         plot_ly(
#                                           type = "scatter",
#                                           mode = "markers"
#                                         ) |>
#                                           add_markers(x = data$position_mb,
#                                                       y = data$dataset,
#                                                       symbol = data$dataset) |>
#                                         add_annotations(x = chrom_midpoints[[chrom]],
#                                                         y = 3,
#                                                         text = chrom,
#                                                         showarrow = FALSE,
#                                                         font = list(size = 16),
#                                                         #xref = "paper",
#                                                         #yref = "paper",
#                                                         xanchor = "center",
#                                                         yanchor = "bottom") |>
#                                           add_segments(x = data_chrom[["position_mb"]],
#                                                        y = data_chrom[["dataset"]],
#                                                        xend = data_chrom[["position_mb_end"]],
#                                                        yend = data_chrom[["dataset"]],
#                                                        color = "black")
#                                       })) |>
#   dplyr::pull(subplot) |>
#   subplot(shareX = FALSE, shareY = TRUE)
