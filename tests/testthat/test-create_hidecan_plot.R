load("data-test/test_input.rda")

test_that("create_hidecan_plot works", {

  x <- list(
    apply_threshold(gwas_res, score_thr = 4),
    apply_threshold(de_res, score_thr = -log10(0.05), log2fc_thr = 0.5),
    apply_threshold(can_res)
  )

  chrom_length <- combine_chrom_length(list(gwas_res, de_res, can_res))

  expect_error(create_hidecan_plot(list(LETTERS[1:5]),
                                   "Expecting a list of 'GWAS_data_thr', 'DE_data_thr' and/or 'CAN_data_thr' objects (see apply_threshold() function)."))

  expect_error(create_hidecan_plot(x, chrom_length), NA)

  p <- create_hidecan_plot(x, chrom_length)
  expect_s3_class(p, "ggplot")

  ## Checking problems with chrom_length
  expect_error(create_hidecan_plot(x, dplyr::rename(chrom_length, chrom = chromosome)), "'chrom_length' argument should be a data-frame with columns 'chromosome' ad length.")
  expect_error(create_hidecan_plot(x, dplyr::bind_rows(chrom_length, chrom_length[1, ])), "Duplicated chromosome names in 'chrom_length' data-frame.")
  expect_error(create_hidecan_plot(x, chrom_length[-1, ]), "The following chromosomes are present in the data but missing from 'chrom_length' data-frame: '.+")

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
})


test_that("hidecan_plot works", {

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
  expect_error(hidecan_plot(), "Should provide at least one non-empty list for 'gwas_list', 'de_list' or 'ca_list' argument.")
  expect_error(hidecan_plot(gwas_list = "TEST"), "Arguments 'gwas_list', 'de_list' or 'can_list' should either be a data-frame or a list, or NULL.")

  expect_error(hidecan_plot(gwas_list = dplyr::rename(gwas_data, chrom = chromosome)), "In 'gwas_list' argument: Input data-frame is missing the following columns: 'chromosome'.")
  expect_error(hidecan_plot(de_list = dplyr::mutate(de_data, score = paste0(score))), "In 'de_list' argument: 'score' column should contain numeric values.")
  expect_error(hidecan_plot(can_list = dplyr::rename(can_data, chrom = chromosome)), "In 'can_list' argument: Input data-frame is missing the following columns: 'chromosome'.")

  ## Check score and log2fc thresholds
  expect_error(hidecan_plot(gwas_list = gwas_data, score_thr = "TEST"), "'score_thr' argument should be a named vector of length 2 with names 'GWAS' and 'DE'.")
  expect_error(hidecan_plot(gwas_list = gwas_data, score_thr = c("Test" = 1, "DE" = 2)), "'score_thr' argument should be a named vector of length 2 with names 'GWAS' and 'DE'.")
  expect_error(hidecan_plot(gwas_list = gwas_data, score_thr = c(2, 2)), "'score_thr' argument should be a named vector of length 2 with names 'GWAS' and 'DE'.")
  expect_error(hidecan_plot(gwas_list = gwas_data, score_thr = c("GWAS" = "1", "DE" = "2")), "'score_thr' argument should be a numeric vector.")
  expect_error(hidecan_plot(gwas_list = gwas_data, log2fc_thr = "TEST"),"'log2fc_thr' argument should be a numeric value.")
  expect_error(hidecan_plot(gwas_list = gwas_data, log2fc_thr = c(1, 2)),"'log2fc_thr' argument should be a single numeric value.")

})
