test_that("GWAS_data_from_gwaspoly works", {

  ## because the example data is a S4 object from the GWASpoly package,
  ## Cannot load it without using GWASpoly
  skip_if_not_installed("GWASpoly")

  x_thr <- readRDS("data-test/test_gwaspoly_res_thr.rda")
  x_nothr <- readRDS("data-test/test_gwaspoly_res.rda")

  expect_error(GWAS_data_from_gwaspoly("TEST"), "'gwaspoly_output' should be a `GWASpoly.fitted` or `GWASpoly.thresh` object (returned by GWASpoly() or set.threshold() functions).", fixed = TRUE)

  expect_error(GWAS_data_from_gwaspoly(x_nothr), NA)
  expect_error(GWAS_data_from_gwaspoly(x_thr), NA)

  expect_error(GWAS_data_from_gwaspoly(x_thr, traits = "TEST"), "The following traits are not present in the input data: 'TEST'.\nPossible values for traits argument are: 'tuber_eye_depth', 'tuber_shape'.")
  expect_error(GWAS_data_from_gwaspoly(x_thr, models = "TEST"), "The following models are not present in the input data: 'TEST'.\nPossible values for models argument are: 'general', 'additive'.")

  expect_error(GWAS_data_from_gwaspoly(x_thr, traits = "tuber_eye_depth"), NA)
  expect_error(GWAS_data_from_gwaspoly(x_thr, models = "additive"), NA)
  expect_error(GWAS_data_from_gwaspoly(x_thr, traits = "tuber_eye_depth", models = "additive"), NA)

  res_nothr <- GWAS_data_from_gwaspoly(x_nothr)
  res_thr <- GWAS_data_from_gwaspoly(x_thr)
  expect_equal(names(res_nothr), c("gwas_data_list", "gwas_data_thr_list", "chrom_length"))
  expect_equal(names(res_thr), c("gwas_data_list", "gwas_data_thr_list", "chrom_length"))

  tm_labels <- tidyr::expand_grid(
    trait = c('tuber_eye_depth', 'tuber_shape'),
    model = c('general', 'additive')
  ) |>
    dplyr::mutate(label = paste0(trait, " (", model, ")")) |>
    dplyr::pull(label)

  expect_equal(names(res_thr[["gwas_data_list"]]), tm_labels)
  expect_equal(purrr::map_chr(res_thr[["gwas_data_list"]], ~ class(.x)[[1]]), rep("GWAS_data", length(tm_labels)) |> setNames(tm_labels))
  expect_equal(res_nothr[["gwas_data_list"]], res_thr[["gwas_data_list"]])

  expect_equal(res_nothr[["gwas_data_thr_list"]], NULL)
  expect_equal(names(res_thr[["gwas_data_thr_list"]]), tm_labels)
  expect_equal(purrr::map_chr(res_thr[["gwas_data_thr_list"]], ~ class(.x)[[1]]), rep("GWAS_data_thr", length(tm_labels)) |> setNames(tm_labels))

  expect_equal(names(res_thr[["chrom_length"]]), c("chromosome", "length"))
  expect_equal(nrow(res_thr[["chrom_length"]]), 13)
  expect_equal(res_nothr[["chrom_length"]], res_thr[["chrom_length"]])
})
