test_that(".check_cols works", {

  data <- tibble::tibble(a = 1:5, b = 1:5, c = LETTERS[1:5])

  expect_error(.check_cols(data, c("a", "b", "d")), "Input data-frame is missing the following columns: 'd'.")
  expect_error(.check_cols(data, c("a", "e", "d")), "Input data-frame is missing the following columns: 'e', 'd'.")
  expect_error(.check_cols(data, c("d", "e", "f")), "Input data-frame is missing the following columns: 'd', 'e', 'f'.")
  expect_error(.check_cols(data, c("a", "b", "d"), "test"), "test is missing the following columns: 'd'.")

  expect_error(.check_cols(data, c("a", "b", "c")), NA)
})
