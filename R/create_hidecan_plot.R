#' creates a HIDECAN plot
#'
#' Creates a HIDECAN plot from a list of filtered GWAS or DE results
#' and/or candidate genes.
#'
#' @param x A list of `GWAS_data_thr`, `DE_data_thr` and/or `CAN_data_thr`
#' produced by the \link{apply_threshold}() function. If named, the names
#' will be integrated to the y-axis labels (use `' '` to exclude the name
#' for any element of the list).
#' @param chrom_length Tibble with columns `chromosome` and `length`, giving
#' for each chromosome its length in bp (see \link{combine_chrom_length})()
#' function.
#' @param colour_genes_by_score Logical, whether to colour the genes by score
#' (`TRUE`) or by log2(fold-change) (`FALSE`). Default value is `TRUE`.
#' @param title Character, title of the plot. Default value is `NULL` (i.e.
#' no title will be added to the plot).
#' @param subtitle Character, subtitle of the plot. Default value is `NULL`
#' (i.e. no subtitle will be added to the plot).
#' @param n_rows Integer, number of rows of facets to create in the plot.
#' Default value is `NULL`.
#' @param n_cols Integer, number of columns of facets to create in the plot.
#' Default value is 2.
#' @param legend_position Character, position of the legend in the plot. Can be
#' `bottom` (default value), `top`, `right` or `left`.
#' @param point_size Numeric, size of the points in the plot. Default value is 3.
#' @param label_size Numeric, size of the gene labels in the plot. Default value is
#' 3.5 (for \link[ggrepel]{geom_label_repel}).
#' @param label_padding Numeric, amount of padding around gene labels in the plot,
#' as unit or number. Default value is 0.15
#' (for \link[ggrepel]{geom_label_repel}).
#' @returns A ggplot.
#' @export
create_hidecan_plot <- function(x,
                                chrom_length,
                                colour_genes_by_score = TRUE,
                                title = NULL,
                                subtitle = NULL,
                                n_rows = NULL,
                                n_cols = 2,
                                legend_position = "bottom",
                                point_size = 3,
                                label_size = 3.5,
                                label_padding = 0.15){

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


  x_types <- purrr::map_chr(x, ~ class(.x)[[1]])

  if(!all(unique(x_types) %in% c("GWAS_data_thr", "DE_data_thr", "CAN_data_thr"))) stop("Expecting a list of 'GWAS_data_thr', 'DE_data_thr' and/or 'CAN_data_thr' objects (see apply_threshold() function).")

  if(is.null(names(x))){
    names(x) <- make.unique(data_type_labels[x_types], " - ")
  } else {
    names(x) <- paste0(names(x), " - ", data_type_labels[x_types])
    names(x) <- sub("^  - ", "", names(x))
  }

  datasets_levels <- rev(names(x))

  toplot <- x |>
    purrr::map_dfr(~ dplyr::mutate(.x, data_type = class(.x)[[1]]),
                   .id = "dataset") |>
    dplyr::mutate(position_mb = position / 1e6,
                  dataset = factor(dataset, levels = datasets_levels)) |>
    dplyr::arrange(score)

  toplot_chroms <- chrom_length |>
    dplyr::mutate(position_mb = length / 1e6) |>
    dplyr::select(chromosome, position_mb) |>
    tidyr::expand_grid(dataset = factor(names(x), levels = datasets_levels)) |>
    dplyr::rename(position_mb_end = position_mb) |>
    dplyr::mutate(position_mb = 0)

  p <- toplot |>
    ggplot2::ggplot(ggplot2::aes(x = position_mb, y = dataset)) +
    ## Chromosome segments
    ggplot2::geom_segment(data = toplot_chroms,
                          ggplot2::aes(xend = position_mb_end, yend = dataset)) +
    ggplot2::facet_wrap(~ chromosome, scales = "free_x", nrow = n_rows, ncol = n_cols) +
    ## General colours and shapes
    ggplot2::scale_colour_manual(values = data_type_colours, labels = data_type_labels) +
    ggplot2::scale_shape_manual(values = data_type_shapes, labels = data_type_labels) +
    ggplot2::guides(colour = ggplot2::guide_legend(title.position = "top",
                                                   title.hjust = 0.5,
                                                   override.aes = list(alpha = 1),
                                                   order = 1),
                    shape = ggplot2::guide_legend(order = 1)) +
    ## vertical position indicators
    ggplot2::geom_vline(ggplot2::aes(xintercept = position_mb, colour = data_type),
                        alpha = 0.7)

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
                          mapping = ggplot2::aes(shape = data_type, fill = score),
                          size = point_size)
  }

  if("GWAS_data_thr" %in% x_types){
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

  if("DE_data_thr" %in% x_types){

    if(colour_genes_by_score){
      p <- p  +
        ## DE points
        ggnewscale::new_scale_fill() +
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
        ggnewscale::new_scale_fill() +
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

  p +
    ## Themes and labs
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = legend_position) +
    ggplot2::labs(title = title,
                  subtitle = subtitle,
                  x = "Position (Mb)",
                  y = NULL,
                  colour = "Position of",
                  shape = "Position of")

}
