test_that("apply_threshold.GWAS_data works", {

  gwas_res <- test_get_gwas()

  ## If no threshold is specified, returns the dataset without rows
  ## that have a missing score
  expect_identical(unclass(apply_threshold(gwas_res)),
                   gwas_res |>
                     dplyr::filter(!is.na(score)) |>
                     unclass())

  ## applying the threshold should get rid of the rows with a score
  ## below the threshold
  res <- apply_threshold(gwas_res, score_thr = 2)
  expect_true(min(res[["score"]]) >= 2)

  ## the log2fc_thr parameter should have no effect on the output
  expect_identical(apply_threshold(gwas_res, score_thr = 2, log2fc_thr = 0),
                   apply_threshold(gwas_res, score_thr = 2, log2fc_thr = 1))

})

test_that("apply_threshold.DE_data works", {

  de_res <- test_get_de()

  ## If no threshold is specified, returns the dataset without rows
  ## that have a missing score
  expect_identical(unclass(apply_threshold(de_res)),
                   de_res |>
                     dplyr::filter(!is.na(score)) |>
                     unclass())

  ## applying the threshold should get rid of the rows with a score
  ## below the threshold
  res <- apply_threshold(de_res, score_thr = 1)
  expect_true(min(res[["score"]]) >= 1)

  res <- apply_threshold(de_res, log2fc_thr = 0.5)
  expect_true(min(abs(res[["log2FoldChange"]])) >= 0.5)
  expect_true(min(res[["log2FoldChange"]]) <= -0.5)

  ## Testing threshold on both score and log2FC at the same time
  res <- apply_threshold(de_res, score_thr = 1, log2fc_thr = 0.5)
  expect_true(min(res[["score"]]) >= 1)
  expect_true(min(abs(res[["log2FoldChange"]])) >= 0.5)
})

test_that("apply_threshold.CAN_data works", {

  can_res <- test_get_can()

  ## Should return the dataset
  expect_identical(unclass(apply_threshold(can_res)),
                   unclass(can_res))

  expect_identical(unclass(apply_threshold(can_res, score_thr = 2)),
                   unclass(can_res))

  expect_identical(unclass(apply_threshold(can_res, log2fc_thr = 2)),
                   unclass(can_res))
})

test_that("apply_threshold.QTL_data works", {

  qtl_res <- test_get_qtl()

  ## If no threshold is specified, returns the dataset without rows
  ## that have a missing score
  expect_identical(
    unclass(apply_threshold(qtl_res)),
    qtl_res |>
      dplyr::filter(!is.na(score)) |>
      unclass()
  )

  ## applying the threshold should get rid of the rows with a score
  ## below the threshold
  res <- apply_threshold(qtl_res, score_thr = 2)
  expect_true(min(res[["score"]]) >= 2)

  ## the log2fc_thr parameter should have no effect on the output
  expect_identical(
    apply_threshold(qtl_res, score_thr = 2, log2fc_thr = 0),
    apply_threshold(qtl_res, score_thr = 2, log2fc_thr = 1)
  )
})

test_that("apply_threshold.CUSTOM_data works", {

  custom_res <- test_get_custom()

  ## If no threshold is specified, returns the dataset without rows
  ## that have a missing score
  expect_identical(unclass(apply_threshold(custom_res)),
                   custom_res |>
                     dplyr::filter(!is.na(score)) |>
                     unclass())

  ## applying the threshold should get rid of the rows with a score
  ## below the threshold
  res <- apply_threshold(custom_res, score_thr = 2)
  expect_true(min(res[["score"]]) >= 2)

  ## the log2fc_thr parameter should have no effect on the output
  expect_identical(apply_threshold(custom_res, score_thr = 2, log2fc_thr = 0),
                   apply_threshold(custom_res, score_thr = 2, log2fc_thr = 1))

})
