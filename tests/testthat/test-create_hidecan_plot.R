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
})
