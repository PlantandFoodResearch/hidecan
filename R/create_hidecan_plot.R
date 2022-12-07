#' Creates a HIDECAN plot
#'
#' Creates a HIDECAN plot from a list of filtered GWAS or DE results
#' and/or candidate genes.
#'
#' @param x A list of `GWAS_data_thr`, `DE_data_thr` and/or `CAN_data_thr`
#' produced by the \code{\link{apply_threshold}()} function. If named, the names
#' will be appended to the y-axis labels (use `' '` as empty name in the list).
#' @param chrom_length Tibble with columns `chromosome` and `length`, giving
#' for each chromosome its length in bp (see \code{\link{combine_chrom_length}()}
#' function).
#' @param remove_empty_chrom Logical, should chromosomes with no significant
#' markers/genes nor candidate genes be removed from the plot? Default value
#' if `FALSE`.
#' @param chroms Character vector, name of chromosomes to include in the plot.
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
#' @returns A ggplot.
#' @examples
#' \dontrun{
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
#' }
#' @export
create_hidecan_plot <- function(x,
                                chrom_length,
                                colour_genes_by_score = TRUE,
                                remove_empty_chrom = FALSE,
                                chroms = NULL,
                                title = NULL,
                                subtitle = NULL,
                                n_rows = NULL,
                                n_cols = 2,
                                legend_position = "bottom",
                                point_size = 3,
                                label_size = 3.5,
                                label_padding = 0.15){

  ## for devtools::check()
  position <- dataset <- score <- chromosome <- position_mb <- position_mb_end <- data_type <- name <- log2FoldChange <- NULL

  ## Labels, colours and shapes
  data_type_labels <- c("GWAS_data_thr" = "GWAS peaks",
                        "DE_data_thr" = "DE genes",
                        "CAN_data_thr" = "Candidate genes")

  data_type_colours <- c("GWAS_data_thr" = "red",
                         "DE_data_thr" = "darkcyan",
                         "CAN_data_thr" = "grey30")

  data_type_shapes <- c("GWAS_data_thr" = 21,
                        "DE_data_thr" = 23,
                        "CAN_data_thr" = 4)


  ## Checking x input
  x_types <- purrr::map_chr(x, ~ class(.x)[[1]])

  if(!all(unique(x_types) %in% c("GWAS_data_thr", "DE_data_thr", "CAN_data_thr"))) stop("Expecting a list of 'GWAS_data_thr', 'DE_data_thr' and/or 'CAN_data_thr' objects (see apply_threshold() function).")

  ## Chromosomes present in the input data
  x_chroms <- x |>
    purrr::map(~ unique(.x[["chromosome"]])) |>
    unlist() |>
    unique()

  if(!is.null(chroms)){
    ## Checking whether the chroms input is valid
    wrong_chroms <- setdiff(chroms, x_chroms)
    if(length(wrong_chroms)) stop("In 'chroms' argument, the following values are not valid chromosome names: '",
                                  paste0(wrong_chroms, collapse = "', '"),
                                  ". Possible values are: '",
                                  paste0(x_chroms, collapse = "', '"),
                                  ".")

    ## Selecting only requested chromosomes in data
    x_chroms <- chroms

    x <- x |>
      purrr::map(
        ~ dplyr::filter(.x, chromosome %in% x_chroms)
      )
  }

  ## Checking chrom_length input
  if(length(setdiff(c("chromosome", "length"), names(chrom_length)))) stop("'chrom_length' argument should be a data-frame with columns 'chromosome' ad length.")
  if(dplyr::n_distinct(chrom_length$chromosome) != nrow(chrom_length)) stop("Duplicated chromosome names in 'chrom_length' data-frame.")

  chrom_length <- chrom_length |>
    dplyr::filter(chromosome %in% x_chroms)

  missing_chroms <- setdiff(x_chroms, chrom_length$chromosome)
  if(length(missing_chroms)) stop("The following chromosomes are present in the data but missing from 'chrom_length' data-frame: '", paste0(missing_chroms, collapse = "', '"), ".")

  ## Removing chromosomes with no significant markers or genes
  if(remove_empty_chrom){
    chrom_length <- chrom_length |>
      dplyr::filter(chromosome %in% x_chroms)

    if(is.factor(chrom_length$chromosome)) chrom_length$chromosome <- droplevels(chrom_length$chromosome)
  }

  ## Making y labels
  if(is.null(names(x))){
    names(x) <- make.unique(data_type_labels[x_types], " - ")
  } else {
    names(x) <- paste0(names(x), " - ", data_type_labels[x_types])
    names(x) <- sub("^  - ", "", names(x))
    names(x) <- make.unique(names(x), ", ")
  }

  ## Making sure that only one of n_rows and n_cols is not NULL
  if(!is.null(n_rows)) n_cols <- NULL

  datasets_levels <- rev(names(x))

  toplot <- x |>
    purrr::map_dfr(~ dplyr::mutate(.x, data_type = class(.x)[[1]]),
                   .id = "dataset") |>
    dplyr::mutate(position_mb = position / 1e6,
                  dataset = factor(dataset, levels = datasets_levels))

  if('score' %in% names(toplot)) toplot <- toplot |> dplyr::arrange(score)

  toplot_chroms <- chrom_length |>
    dplyr::mutate(position_mb = length / 1e6) |>
    dplyr::select(chromosome, position_mb) |>
    tidyr::expand_grid(dataset = factor(names(x), levels = datasets_levels)) |>
    dplyr::rename(position_mb_end = position_mb) |>
    dplyr::mutate(position_mb = 0)

  ## Making sure that the order of the "Position of" legends matches the order
  ## in which the different data types appear in the y-axis
  unique_data_types <- toplot |>
    dplyr::select(dataset, data_type) |>
    dplyr::distinct() |>
    dplyr::arrange(dplyr::desc(dataset)) |>
    dplyr::pull(data_type) |>
    unique()

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
                                 labels = data_type_labels[unique_data_types]) +
    ggplot2::scale_shape_manual(values = data_type_shapes[unique_data_types],
                                labels = data_type_labels[unique_data_types]) +
    ggplot2::guides(colour = ggplot2::guide_legend(title.position = "top",
                                                   title.hjust = 0.5,
                                                   override.aes = list(alpha = 1),
                                                   order = 1),
                    shape = ggplot2::guide_legend(order = 1))

  if("CAN_data_thr" %in% x_types){
    p <- p +
      ## Candidate genes labels
      ggrepel::geom_label_repel(ggplot2::aes(label = name),
                                nudge_y = 0.5,
                                na.rm = TRUE,
                                size = label_size,
                                label.padding = label_padding,
                                alpha = 0.5) +
      ## Candidate points
      ggplot2::geom_point(data = dplyr::filter(toplot, data_type == "CAN_data_thr"),
                          mapping = ggplot2::aes(shape = data_type),
                          size = point_size)
  }

  if("DE_data_thr" %in% x_types){

    if(colour_genes_by_score){

      p <- p +
        ## DE points
        ggplot2::geom_point(data = dplyr::filter(toplot, data_type == "DE_data_thr"),
                            mapping = ggplot2::aes(shape = data_type, fill = score),
                            size = point_size) +
        viridis::scale_fill_viridis("DE gene score",
                                    option = "viridis",
                                    guide = ggplot2::guide_colourbar(title.position="top",
                                                                     title.hjust = 0.5,
                                                                     order = 3))
    } else {

      p <- p  +
        ## DE points
        ggplot2::geom_point(data = dplyr::filter(toplot, data_type == "DE_data_thr"),
                            mapping = ggplot2::aes(shape = data_type, fill = log2FoldChange),
                            size = point_size) +
        ggplot2::scale_fill_gradient2("DE gene log2(fold-change)",
                                      low = "darkblue",
                                      mid = "white",
                                      high = "firebrick",
                                      guide = ggplot2::guide_colourbar(title.position="top",
                                                                       title.hjust = 0.5,
                                                                       order = 3))
    }

  }

  if("GWAS_data_thr" %in% x_types){

    if("DE_data_thr" %in% x_types) p <- p + ggnewscale::new_scale_fill()

    p <- p +
      ## GWAS points
      ggplot2::geom_point(data = dplyr::filter(toplot, data_type == "GWAS_data_thr"),
                          mapping = ggplot2::aes(shape = data_type, fill = score),
                          size = point_size) +
      viridis::scale_fill_viridis("GWAS marker score",
                                  option = "plasma",
                                  guide = ggplot2::guide_colourbar(title.position="top",
                                                                   title.hjust = 0.5,
                                                                   order = 2))
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
#' Wrapper function to create a HIDECAN plot from GWAS results,
#' DE results or candidate genes.
#'
#' @param gwas_list Data-frame or list of data-frames containing GWAS results,
#' each with at least a `chromosome`, `position` and either `padj` or `score`
#' columns. If a named list, the names will be used in the plot.
#' @param de_list Data-frame or list of data-frames containing DE results,
#' each with at least a `chromosome`, `start`, `end`, `log2FoldChange` and
#' either `padj` or `score` columns. If a named list, the names will be used in
#' the plot.
#' @param can_list Data-frame or list of data-frames containing candidate genes,
#' each with at least a `chromosome`, `start`, `end` and `name` columns. If a
#' named list, the names will be used in the plot.
#' @param score_thr_gwas Numeric, the score threshold for GWAS results that will be used to
#' select which markers will be plotted. Default value is 4.
#' @param score_thr_de Numeric, the score threshold for DE results that will be used to
#' select which markers will be plotted. Default value is 2.
#' @param log2fc_thr Numeric, the log2(fold-change) threshold that will be used
#' to select which genes will be plotted. Default value is 1.
#' @param chrom_length Optional, tibble with columns `chromosome` and `length`,
#' giving for each chromosome its length in bp. If `NULL` (the default), will
#' be inferred from the GWAS, DE and candidate gene data.
#' @inheritParams create_hidecan_plot
#' @returns a ggplot.
#' @examples
#' \dontrun{
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
                         title = NULL,
                         subtitle = NULL,
                         n_rows = NULL,
                         n_cols = 2,
                         legend_position = "bottom",
                         point_size = 3,
                         label_size = 3.5,
                         label_padding = 0.15
){

  ## Input should either be NULL, a data-frame or a list
  if(!all(purrr::map_lgl(list(gwas_list, de_list, can_list), ~ (is.list(.x) | is.null(.x))))) stop("Arguments 'gwas_list', 'de_list' or 'can_list' should either be a data-frame or a list, or NULL.")

  ## If providing a single data-frame, turn into a list
  if(is.data.frame(gwas_list)) gwas_list <- list(gwas_list)
  if(is.data.frame(de_list)) de_list <- list(de_list)
  if(is.data.frame(can_list)) can_list <- list(can_list)

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

  x <- c(gwas_list, de_list, can_list)

  if(!length(x)) stop("Should provide at least one non-empty list for 'gwas_list', 'de_list' or 'can_list' argument.")

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


  ## Apply threshold to datasets
  score_thr <- c("GWAS_data" = score_thr_gwas, "DE_data" = score_thr_de, "CAN_data" = 0)

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
                      title,
                      subtitle,
                      n_rows,
                      n_cols,
                      legend_position,
                      point_size,
                      label_size,
                      label_padding)
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
