test_that("new_GWAS_data works", {

  data <- tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = 1:5)

  ## Expecting a tibble as input
  expect_error(new_GWAS_data(as.data.frame(data)))

  res <- new_GWAS_data(data)

  ## output should be equal to input
  expect_equal(unclass(res), unclass(data))

  ## class should be GWAS_data
  expect_equal(class(res)[1], "GWAS_data")

})


test_that("GWAS_data works", {

  data <- tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = 1:5)

  ## Need the correct columns
  expect_error(GWAS_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, score = 1:5)),
               "Input data-frame is missing the following columns: 'chromosome', 'position'.")

  ## Checking correct input type
  expect_error(GWAS_data(tibble::tibble(chromosome = LETTERS[1:5], position = LETTERS[1:5], score = 1:5)),
               "'position' column should contain numeric values.")
  expect_error(GWAS_data(tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = LETTERS[1:5])),
               "'score' column should contain numeric values.")

  ## Correct input should trigger no error
  expect_error(GWAS_data(data), NA)

  ## output should be equal to input
  expect_equal(unclass(GWAS_data(data)), unclass(data))


  data_2 <- as.data.frame(data)
  rownames(data_2) <- letters[1:5]

  ## By default, rownames of a data-frame should be ignored
  res <- GWAS_data(data_2)
  expect_equal(colnames(res),  colnames(data_2))

  ## If specified, rownames should be saved in new column
  res <- GWAS_data(data_2, keep_rownames_as = "rowname")
  expect_equal(res[["rowname"]], rownames(data_2))
})



test_that("new_DE_data works", {

  data <- tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5, log2FoldChange = 1:5)

  ## Expecting a tibble as input
  expect_error(new_DE_data(as.data.frame(data)))

  res <- new_DE_data(data)

  ## output should be equal to input
  expect_equal(unclass(res), unclass(data))

  ## class should be DE_data
  expect_equal(class(res)[1], "DE_data")

})

test_that("DE_data works", {

  data <- tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5, log2FoldChange = 1:5)

  ## Need the correct columns
  expect_error(DE_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, score = 1:5)),
               "Input data-frame is missing the following columns: 'chromosome', 'start', 'end', 'log2FoldChange'.")

  ## Checking correct input type
  expect_error(DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = LETTERS[1:5], end = 2:6, score = 1:5, log2FoldChange = 1:5)),
               "'start' and 'end' columns should contain numeric values.")
  expect_error(DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = LETTERS[1:5], score = 1:5, log2FoldChange = 1:5)),
               "'start' and 'end' columns should contain numeric values.")
  expect_error(DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = LETTERS[1:5], log2FoldChange = 1:5)),
               "'score' column should contain numeric values.")
  expect_error(DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5, log2FoldChange = letters[1:5])),
               "'log2FoldChange' column should contain numeric values.")

  ## Correct input should trigger no error
  expect_error(DE_data(data), NA)

  ## output should be equal to input (except for the position column)
  expect_equal(unclass(dplyr::select(DE_data(data), -position)), unclass(data))

  ## position column should be the average of start and end position
  expect_equal(DE_data(data)[["position"]], purrr::map2_dbl(data[["start"]], data[["end"]], mean))

  data_2 <- as.data.frame(data)

  rownames(data_2) <- letters[1:5]

  ## By default, rownames of a data-frame should be ignored
  res <- DE_data(data_2)
  expect_equal(colnames(res),  c(colnames(data_2), "position"))

  ## If specified, rownames should be saved in new column
  res <- DE_data(data_2, keep_rownames_as = "rowname")
  expect_equal(res[["rowname"]], rownames(data_2))
})



test_that("new_CAN_data works", {

  data <- tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, name = letters[1:5])

  ## Expecting a tibble as input
  expect_error(new_CAN_data(as.data.frame(data)))

  res <- new_CAN_data(data)

  ## output should be equal to input
  expect_equal(unclass(res), unclass(data))

  ## class should be CAN_data
  expect_equal(class(res)[1], "CAN_data")
})

test_that("CAN_data works", {

  data <- tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, name = letters[1:5])

  ## Need the correct columns
  expect_error(CAN_data(tibble::tibble(chromosome = LETTERS[1:5], b = 1:5, score = 1:5)),
               "Input data-frame is missing the following columns: 'start', 'end', 'name'.")

  ## Checking correct input type
  expect_error(CAN_data(tibble::tibble(chromosome = LETTERS[1:5], start = LETTERS[1:5], end = 2:6, name = letters[1:5])),
               "'start' and 'end' columns should contain numeric values.")
  expect_error(CAN_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = LETTERS[1:5], name = letters[1:5])),
               "'start' and 'end' columns should contain numeric values.")
  expect_error(CAN_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, name = 1:5)),
               "'name' column should contain character values.")

  ## Correct input should trigger no error
  expect_error(CAN_data(data), NA)

  ## output should be equal to input, except for the position column
  expect_equal(unclass(dplyr::select(CAN_data(data), -position)), unclass(data))

  ## position column should be the average of start and end position
  expect_equal(CAN_data(data)[["position"]], purrr::map2_dbl(data[["start"]], data[["end"]], mean))

  data_2 <- as.data.frame(data)
  rownames(data_2) <- letters[1:5]

  ## By default, rownames of a data-frame should be ignored
  res <- CAN_data(data_2)
  expect_equal(colnames(res),  c(colnames(data_2), "position"))

  ## If specified, rownames should be saved in new column
  res <- CAN_data(data_2, keep_rownames_as = "rowname")
  expect_equal(res[["rowname"]], rownames(data_2))
})
