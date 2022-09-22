test_that("new_GWAS_data works", {

  data <- tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = 1:5)

  expect_error(new_GWAS_data(as.data.frame(data)))
  expect_error(new_GWAS_data(data, 2))

  res <- new_GWAS_data(data)

  expect_equal(unclass(res), unclass(data), ignore_attr = TRUE)

  expect_equal(class(res)[1], "GWAS_data")
  expect_equal(attr(res, "title"), "GWAS")

  res <-  new_GWAS_data(data, "test")
  expect_equal(attr(res, "title"), "test")
})


test_that("GWAS_data works", {

  data <- tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = 1:5)

  expect_error(GWAS_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, score = 1:5)),
               "Input data-frame is missing the following columns: 'chromosome', 'position'.")

  expect_error(GWAS_data(data), NA)

  expect_equal(unclass(GWAS_data(data)), unclass(data), ignore_attr = TRUE)

  data_2 <- as.data.frame(data)
  rownames(data_2) <- letters[1:5]

  res <- GWAS_data(data_2)
  expect_equal(colnames(res),  colnames(data_2))

  res <- GWAS_data(data_2, keep_rownames_as = "rowname")
  expect_equal(res[["rowname"]], rownames(data_2))
})



test_that("new_DE_data works", {

  data <- tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = 1:5)

  expect_error(new_DE_data(as.data.frame(data)))
  expect_error(new_DE_data(data, 2))

  res <- new_DE_data(data)

  expect_equal(unclass(res), unclass(data), ignore_attr = TRUE)

  expect_equal(class(res)[1], "DE_data")
  expect_equal(attr(res, "title"), "DE")

  res <-  new_DE_data(data, "test")
  expect_equal(attr(res, "title"), "test")
})
