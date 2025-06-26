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
  expect_error(
    GWAS_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, score = 1:5)),
    "Input data-frame is missing the following columns: 'chromosome', 'position'."
  )
  expect_error(
    GWAS_data(tibble::tibble(a = LETTERS[1:5], b = 1:5)),
    "Input data-frame should have either a 'score' or a 'padj' column."
  )
  expect_error(
    GWAS_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, padj = LETTERS[1:5])),
    "'padj' column in input data-frame should contain numeric values."
  )
  expect_error(
    GWAS_data(tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, padj = 1:5)),
    NA
  ) ## no error if no score but padj
  expect_equal(
    GWAS_data(tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, padj = 1:5))[["score"]],
    -log10(1:5)
  )

  ## Checking correct input type
  expect_error(
    GWAS_data(tibble::tibble(chromosome = LETTERS[1:5], position = LETTERS[1:5], score = 1:5)),
    "'position' column should contain numeric values."
  )
  expect_error(
    GWAS_data(tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = LETTERS[1:5])),
    "'score' column should contain numeric values."
  )

  ## Correct input should trigger no error
  expect_error(GWAS_data(data), NA)

  ## output should be equal to input
  expect_equal(unclass(GWAS_data(data)), unclass(data))


  data_2 <- as.data.frame(data)
  rownames(data_2) <- letters[1:5]

  ## By default, rownames of a data-frame should be ignored
  res <- GWAS_data(data_2)
  expect_equal(colnames(res), colnames(data_2))

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
  expect_error(
    DE_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, start = 1:5, end = 2:6, log2FoldChange = 1:5, score = 1:5)),
    "Input data-frame is missing the following columns: 'chromosome'."
  )
  expect_error(
    DE_data(tibble::tibble(a = LETTERS[1:5], b = 1:5)),
    "Input data-frame should have either a 'score' or a 'padj' column."
  )
  expect_error(
    DE_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, padj = LETTERS[1:5])),
    "'padj' column in input data-frame should contain numeric values."
  )
  expect_error(
    DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, padj = 1:5, log2FoldChange = 1:5)),
    NA
  ) ## no error if no score but padj
  expect_equal(
    DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, padj = 1:5, log2FoldChange = 1:5))[["score"]],
    -log10(1:5)
  )
  expect_error(
    DE_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, score = 1:5)),
    "Input data-frame should have either a 'log2FoldChange' or a 'foldChange' column."
  )
  expect_error(
    DE_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, score = 1:5, foldChange = LETTERS[1:5])),
    "'foldChange' column in input data-frame should contain numeric values."
  )
  expect_error(
    DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5, foldChange = 1:5)),
    NA
  ) ## no error if no score but padj
  expect_equal(
    DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5, foldChange = 1:5))[["log2FoldChange"]],
    log2(1:5)
  )
  expect_error(
    DE_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, score = 1:5, log2FoldChange = 1:5)),
    "Input data-frame should have a 'start' and an 'end' column."
  )
  expect_error(
    DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = LETTERS[1:5], end = 2:6, score = 1:5, log2FoldChange = 1:5)),
    "'start' and 'end' columns should contain numeric values."
  )
  expect_error(
    DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = LETTERS[1:5], score = 1:5, log2FoldChange = 1:5)),
    "'start' and 'end' columns should contain numeric values."
  )
  expect_equal(
    DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5, foldChange = 1:5))[["position"]],
    seq(1.5, 5.5, 1)
  )

  ## Checking correct input type
  expect_error(
    DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = LETTERS[1:5], log2FoldChange = 1:5)),
    "'score' column should contain numeric values."
  )
  expect_error(
    DE_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5, log2FoldChange = letters[1:5])),
    "'log2FoldChange' column should contain numeric values."
  )

  ## Correct input should trigger no error
  expect_error(DE_data(data), NA)

  ## output should be equal to input (except for the position column)
  expect_equal(unclass(dplyr::select(DE_data(data), -position)), unclass(data))

  data_2 <- as.data.frame(data)

  rownames(data_2) <- letters[1:5]

  ## By default, rownames of a data-frame should be ignored
  res <- DE_data(data_2)
  expect_equal(colnames(res), c(colnames(data_2), "position"))

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
  expect_error(
    CAN_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, start = 1:5, end = 2:6)),
    "Input data-frame is missing the following columns: 'chromosome', 'name'."
  )
  expect_error(
    CAN_data(tibble::tibble(a = LETTERS[1:5], b = 1:5)),
    "Input data-frame should have a 'start' and an 'end' column."
  )
  expect_error(
    CAN_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, start = LETTERS[1:5], end = LETTERS[1:5])),
    "'start' and 'end' columns should contain numeric values."
  )
  expect_equal(
    CAN_data(tibble::tibble(chromosome = LETTERS[1:5], name = letters[1:5], start = 1:5, end = 2:6))[["position"]],
    seq(1.5, 5.5, 1)
  )

  ## Checking correct input type
  expect_error(
    CAN_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, name = 1:5)),
    "'name' column should contain character values."
  )

  ## Correct input should trigger no error
  expect_error(CAN_data(data), NA)

  ## output should be equal to input, except for the position column
  expect_equal(unclass(dplyr::select(CAN_data(data), -position)), unclass(data))

  data_2 <- as.data.frame(data)
  rownames(data_2) <- letters[1:5]

  ## By default, rownames of a data-frame should be ignored
  res <- CAN_data(data_2)
  expect_equal(colnames(res), c(colnames(data_2), "position"))

  ## If specified, rownames should be saved in new column
  res <- CAN_data(data_2, keep_rownames_as = "rowname")
  expect_equal(res[["rowname"]], rownames(data_2))
})

test_that("new_QTL_data works", {
  data <- tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5)

  ## Expecting a tibble as input
  expect_error(new_QTL_data(as.data.frame(data)))

  res <- new_QTL_data(data)

  ## output should be equal to input
  expect_equal(unclass(res), unclass(data))

  ## class should be DE_data
  expect_equal(class(res)[1], "QTL_data")
})

test_that("QTL_data works", {
  data <- tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5)

  ## Need the correct columns
  expect_error(
    QTL_data(
      tibble::tibble(
        a = LETTERS[1:5],
        b = 1:5,
        start = 1:5,
        end = 2:6,
        score = 1:5
      )
    ),
    "Input data-frame is missing the following columns: 'chromosome'."
  )
  expect_error(
    QTL_data(tibble::tibble(a = LETTERS[1:5], b = 1:5)),
    "Input data-frame should have either a 'score' or a 'padj' column."
  )
  expect_error(
    QTL_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, padj = LETTERS[1:5])),
    "'padj' column in input data-frame should contain numeric values."
  )
  expect_no_error(
    QTL_data(
      tibble::tibble(
        chromosome = LETTERS[1:5],
        start = 1:5,
        end = 2:6,
        padj = 1:5
      )
    )
  ) ## no error if no score but padj
  expect_equal(
    QTL_data(
      tibble::tibble(
        chromosome = LETTERS[1:5],
        start = 1:5,
        end = 2:6,
        padj = 1:5
      )
    )[["score"]],
    -log10(1:5)
  )

  expect_error(
    QTL_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, score = 1:5)),
    "Input data-frame should have a 'start' and an 'end' column."
  )
  expect_error(
    QTL_data(tibble::tibble(chromosome = LETTERS[1:5], start = LETTERS[1:5], end = 2:6, score = 1:5)),
    "'start' and 'end' columns should contain numeric values."
  )
  expect_error(
    QTL_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = LETTERS[1:5], score = 1:5)),
    "'start' and 'end' columns should contain numeric values."
  )
  expect_equal(
    QTL_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = 1:5))[["position"]],
    seq(1.5, 5.5, 1)
  )

  ## Checking correct input type
  expect_error(
    QTL_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 2:6, score = LETTERS[1:5])),
    "'score' column should contain numeric values."
  )

  ## Correct input should trigger no error
  expect_no_error(QTL_data(data))

  ## output should be equal to input (except for the position column)
  expect_equal(unclass(dplyr::select(QTL_data(data), -position)), unclass(data))

  data_2 <- as.data.frame(data)

  rownames(data_2) <- letters[1:5]

  ## By default, rownames of a data-frame should be ignored
  res <- QTL_data(data_2)
  expect_equal(colnames(res), c(colnames(data_2), "position"))

  ## If specified, rownames should be saved in new column
  res <- QTL_data(data_2, keep_rownames_as = "rowname")
  expect_equal(res[["rowname"]], rownames(data_2))
})


test_that("new_CUSTOM_data works", {
  data <- tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = 1:5)

  ## Expecting a tibble as input
  expect_error(new_CUSTOM_data(as.data.frame(data)))

  res <- new_CUSTOM_data(data)

  ## output should be equal to input
  expect_equal(unclass(res), unclass(data))

  ## class should be CUSTOM_data
  expect_equal(class(res)[1], "CUSTOM_data")
})


test_that("CUSTOM_data works", {
  ## Need the correct columns
  expect_error(
    CUSTOM_data(tibble::tibble(a = LETTERS[1:5], b = 1:5)),
    "Input data-frame should have a 'start' and an 'end' column, as there is no 'position' column."
  )

  expect_error(
    CUSTOM_data(tibble::tibble(a = LETTERS[1:5], b = 1:5, position = 1:5)),
    "Input data-frame is missing the following columns: 'chromosome', 'score'."
  )

  expect_no_error(
    CUSTOM_data(tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = 1:5)),
  )
  expect_no_error(
    CUSTOM_data(tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 1:5, score = 1:5)),
  )

  ## Checking correct input type
  expect_error(
    CUSTOM_data(tibble::tibble(chromosome = LETTERS[1:5], position = LETTERS[1:5], score = 1:5)),
    "'position' column should contain numeric values."
  )
  expect_error(
    CUSTOM_data(tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = LETTERS[1:5])),
    "'score' column should contain numeric values."
  )

  ## Correct input should trigger no error
  data_pos <- tibble::tibble(chromosome = LETTERS[1:5], position = 1:5, score = 1:5)
  data_se <- tibble::tibble(chromosome = LETTERS[1:5], start = 1:5, end = 3:7, score = 1:5)
  data_full <- tibble::tibble(
    chromosome = LETTERS[1:5], position = 2:6, start = 1:5, end = 3:7, score = 1:5
  )
  expect_no_error(CUSTOM_data(data_pos))
  expect_no_error(CUSTOM_data(data_se))
  expect_no_error(CUSTOM_data(data_full))

  ## output should be equal to input
  expect_equal(unclass(CUSTOM_data(data_full)), unclass(data_full))
  expect_equal(
    unclass(CUSTOM_data(data_pos)),
    data_pos |> dplyr::mutate(start = position, end = position) |> unclass()
  )
  expect_equal(
    unclass(CUSTOM_data(data_se)),
    data_se |> dplyr::mutate(position = (start + end) / 2) |> unclass()
  )

  data_2 <- as.data.frame(data_full)
  rownames(data_2) <- letters[1:5]

  ## By default, rownames of a data-frame should be ignored
  res <- CUSTOM_data(data_2)
  expect_equal(colnames(res), colnames(data_2))

  ## If specified, rownames should be saved in new column
  res <- CUSTOM_data(data_2, keep_rownames_as = "rowname")
  expect_equal(res[["rowname"]], rownames(data_2))
})
