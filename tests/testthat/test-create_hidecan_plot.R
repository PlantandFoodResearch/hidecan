test_that("create_hidecan_plot works", {
  gwas_res <- test_get_gwas()
  de_res <- test_get_de()
  can_res <- test_get_can()

  x <- list(
    apply_threshold(gwas_res, score_thr = 4),
    apply_threshold(de_res, score_thr = -log10(0.05), log2fc_thr = 0.5),
    apply_threshold(can_res)
  )

  chrom_length <- combine_chrom_length(list(gwas_res, de_res, can_res))

  expect_error(
    create_hidecan_plot(list(LETTERS[1:5]), chrom_length),
    "Expecting a list of 'GWAS_data_thr', 'DE_data_thr', 'CAN_data_thr' and/or 'CUSTOM_data_thr' objects (see apply_threshold() function).",
    fixed = TRUE
  )

  expect_no_error(create_hidecan_plot(x, chrom_length))
  expect_no_error(create_hidecan_plot(x, chrom_length, remove_empty_chrom = TRUE))

  p <- create_hidecan_plot(x, chrom_length)
  expect_s3_class(p, "ggplot")

  ## Checking problems with chrom_length
  expect_error(
    create_hidecan_plot(x, dplyr::rename(chrom_length, chrom = chromosome)),
    "'chrom_length' argument should be a data-frame with columns 'chromosome' and 'length'."
  )
  expect_error(
    create_hidecan_plot(x, dplyr::bind_rows(chrom_length, chrom_length[1, ])),
    "Duplicated chromosome names in 'chrom_length' data-frame."
  )
  expect_error(
    create_hidecan_plot(x, chrom_length[-13, ]),
    "The following chromosomes are present in the data but missing from 'chrom_length' data-frame: '.+"
  )

  ## Checking names on y-axis
  expect_equal(p$data[["dataset"]] |>
                 unique() |>
                 sort(),
               c("Candidate genes", "DE genes", "GWAS peaks") |>
                 factor())

  p <- create_hidecan_plot(list(x[[1]], x[[2]], x[[1]]), chrom_length)
  expect_equal(p$data[["dataset"]] |>
                 unique() |>
                 sort(),
               c("GWAS peaks - 1", "DE genes", "GWAS peaks") |>
                 factor(levels = c("GWAS peaks - 1", "DE genes", "GWAS peaks")))

  p <- create_hidecan_plot(list("GWAS 1" = x[[1]], " " =  x[[2]], "GWAS 2" =  x[[1]]), chrom_length)
  expect_equal(p$data[["dataset"]] |>
                 unique() |>
                 sort(),
               c("GWAS 2 - GWAS peaks", "DE genes", "GWAS 1 - GWAS peaks") |>
                 factor(levels = c("GWAS 2 - GWAS peaks", "DE genes", "GWAS 1 - GWAS peaks")))

  p <- create_hidecan_plot(list("GWAS 1" = x[[1]], " " =  x[[2]], "GWAS 1" =  x[[1]]), chrom_length)
  expect_equal(p$data[["dataset"]] |>
                 unique() |>
                 sort(),
               c("GWAS 1 - GWAS peaks, 1", "DE genes", "GWAS 1 - GWAS peaks") |>
                 factor(levels = c("GWAS 1 - GWAS peaks, 1", "DE genes", "GWAS 1 - GWAS peaks")))

  ## Checking chrom_limits arguments
  expect_error(
    create_hidecan_plot(x, chrom_length, chrom_limits = "TEST"),
    "The chrom_limits argument should either be an integer vector of length 2 or a named list where each element is an integer vector of length 2."
  )
  expect_no_error(create_hidecan_plot(x, chrom_length, chrom_limits = c(3e6, 5e6)))
  expect_error(
    create_hidecan_plot(x, chrom_length, chrom_limits = list(c(3e6, 5e6))),
    "The chrom_limits argument should be a named list, with the names corresponding to chromosomes' name."
  )
  expect_error(
    create_hidecan_plot(x, chrom_length, chrom_limits = list("TEST" = c(3e6, 5e6))),
    "In chrom_limits argument: 'TEST' are not valid chromosome names. Possible names are: '.+"
  )
  expect_no_error(create_hidecan_plot(x, chrom_length, chrom_limits = list("ST4.03ch06" = c(3e6, 5e6))))
  expect_error(
    create_hidecan_plot(x, chrom_length, chrom_limits = list("ST4.03ch06" = c(3e6, 5e6, 5e6))),
    "The chrom_limits argument should be a named list where each element is an integer vector of length 2."
  )
  expect_error(
    create_hidecan_plot(x, chrom_length, chrom_limits = list("ST4.03ch06" = c("a", "B"))),
    "The chrom_limits argument should be a named list where each element is an integer vector of length 2."
  )

  p <- create_hidecan_plot(x, chrom_length, chrom_limits = c(3e6, 5e6))
  cl_df <- get_chrom_limits_plot(p)
  expect_true(all(cl_df$position_mb == 3))
  expect_true(all(cl_df$position_mb_end == 5))

  p <- create_hidecan_plot(x, chrom_length, chrom_limits = list("ST4.03ch06" = c(3e6, 5e6)))
  cl_df <- get_chrom_limits_plot(p)
  expect_equal(dplyr::filter(cl_df, chromosome == "ST4.03ch06")$position_mb, 3)
  expect_equal(dplyr::filter(cl_df, chromosome == "ST4.03ch06")$position_mb_end, 5)
  expect_true(all(dplyr::filter(cl_df, chromosome != "ST4.03ch06")$position_mb == 0))
  expect_true(all(dplyr::filter(cl_df, chromosome != "ST4.03ch06")$position_mb_end > 5))

  p <- create_hidecan_plot(x, chrom_length, chrom_limits = list("ST4.03ch06" = c(3e6, 500e6)))
  cl_df <- get_chrom_limits_plot(p)
  expect_equal(dplyr::filter(cl_df, chromosome == "ST4.03ch06")$position_mb, 3)
  expect_equal(dplyr::filter(cl_df, chromosome == "ST4.03ch06")$position_mb_end |> round(2), 59.36)
  expect_true(all(dplyr::filter(cl_df, chromosome != "ST4.03ch06")$position_mb == 0))
  expect_true(all(dplyr::filter(cl_df, chromosome != "ST4.03ch06")$position_mb_end > 5))
})


test_that("create_hidecan_plot works with custom tracks", {
  gwas_res <- test_get_gwas()
  de_res <- test_get_de()
  can_res <- test_get_can()
  custom_res <- test_get_custom()

  x <- list(
    apply_threshold(gwas_res, score_thr = 4),
    apply_threshold(de_res, score_thr = -log10(0.05), log2fc_thr = 0.5),
    apply_threshold(can_res),
    apply_threshold(custom_res, score_thr = 4.95)
  )

  chrom_length <- combine_chrom_length(list(gwas_res, de_res, can_res, custom_res))

  expect_no_error(create_hidecan_plot(x, chrom_length))

  custom_aes <- list(
    "CUSTOM_data_thr" = list(
      "y_label" = "QTL peaks",
      "line_colour" = "purple",
      "point_shape" = 21,
      "show_name" = FALSE,
      fill_scale = viridis::scale_fill_viridis(
        "QTL marker score",
        option = "inferno",
        guide = ggplot2::guide_colourbar(title.position = "top",
                                         title.hjust = 0.5,
                                         order = 4)
      )
    )
  )

  expect_no_error(create_hidecan_plot(x, chrom_length, custom_aes = custom_aes))
})

test_that("hidecan_plot works", {

  gwas_res <- test_get_gwas()
  de_res <- test_get_de()
  can_res <- test_get_can()

  ## getting input as tibbles
  gwas_data <- gwas_res; class(gwas_data) <- class(gwas_data)[-1]
  de_data <- de_res; class(de_data) <- class(de_data)[-1]
  can_data <- can_res; class(can_data) <- class(can_data)[-1]

  ## Minimal graph should work
  expect_error(hidecan_plot(gwas_list = gwas_res), NA)
  expect_error(hidecan_plot(de_list = de_res), NA)
  expect_message(hidecan_plot(de_list = de_res), NA) ## got rid of the duplicated x scale
  expect_error(hidecan_plot(can_list = can_res), NA)

  expect_error(hidecan_plot(gwas_list = gwas_data,
                            de_list = de_data,
                            can_list = can_data), NA)

  ## Checking GWAS, DE and CAN input
  expect_error(hidecan_plot(), "Should provide at least one non-empty list for 'gwas_list', 'de_list', 'can_list' or 'custom_list' argument.")
  expect_error(hidecan_plot(gwas_list = "TEST"), "Arguments 'gwas_list', 'de_list' or 'can_list' should either be a data-frame or a list, or NULL.")

  expect_error(hidecan_plot(gwas_list = dplyr::rename(gwas_data, chrom = chromosome)), "In 'gwas_list' argument: Input data-frame is missing the following columns: 'chromosome'.")
  expect_error(hidecan_plot(de_list = dplyr::mutate(de_data, score = paste0(score))), "In 'de_list' argument: 'score' column should contain numeric values.")
  expect_error(hidecan_plot(can_list = dplyr::rename(can_data, chrom = chromosome)), "In 'can_list' argument: Input data-frame is missing the following columns: 'chromosome'.")

  ## Check score and log2fc thresholds
  expect_error(hidecan_plot(gwas_list = gwas_data, score_thr_gwas = "TEST"), "'score_thr_gwas' argument should be a numeric value.")
  expect_error(hidecan_plot(gwas_list = gwas_data, score_thr_gwas = 1:2), "'score_thr_gwas' argument should be a single numeric value.")
  expect_error(hidecan_plot(gwas_list = gwas_data, score_thr_de = "TEST"), "'score_thr_de' argument should be a numeric value.")
  expect_error(hidecan_plot(gwas_list = gwas_data, score_thr_de = 1:2), "'score_thr_de' argument should be a single numeric value.")
  expect_error(hidecan_plot(gwas_list = gwas_data, log2fc_thr = "TEST"),"'log2fc_thr' argument should be a numeric value.")
  expect_error(hidecan_plot(gwas_list = gwas_data, log2fc_thr = c(1, 2)),"'log2fc_thr' argument should be a single numeric value.")

})
