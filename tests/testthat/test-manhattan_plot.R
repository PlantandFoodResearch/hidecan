test_that("manhattan_plot works", {

  ## Setup
  x <- get_example_data()[["GWAS"]]

  expect_error(manhattan_plot(x), NA)
  expect_error(manhattan_plot(x, 4), NA)

  expect_error(manhattan_plot(list(x, x)), NA)
  expect_error(manhattan_plot(list("Trait A" = x, "Trait B" = x)), NA)

  expect_error(manhattan_plot(x, chroms = "TEST"), "In chroms argument: 'TEST' are not valid chromosome names. Possible names are:.+")
  expect_error(manhattan_plot(x, chroms = "ST4.03ch00"), NA)

  expect_error(manhattan_plot(x, chrom_col = c("dodgerblue1", "dodgerblue4")), NA)
  expect_error(manhattan_plot(x, chrom_col = rep_len(c("dodgerblue1", "dodgerblue4"), 20)), NA)

  p <- manhattan_plot(x)
  expect_s3_class(p, "ggplot")
})
