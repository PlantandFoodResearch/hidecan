load("data-test/test_input.rda")

test_that("compute_chrom_length.GWAS_data works", {

  chr_length <- compute_chrom_length(gwas_res)

  ## One value per chromosome
  expect_equal(chr_length[["chromosome"]], sort(unique(gwas_res[["chromosome"]])))

  ## Each value is the max position within the chromosome
  manual_vals <- gwas_res[["chromosome"]] |>
    unique() |>
    sort() |>
    purrr::map_dbl(~ max(gwas_res$position[gwas_res$chromosome == .x]))

  expect_equal(chr_length[["length"]],
               manual_vals)

})


test_that("compute_chrom_length.DE_data works", {

  chr_length <- compute_chrom_length(de_res)

  ## One value per chromosome
  expect_equal(chr_length[["chromosome"]], sort(unique(de_res[["chromosome"]])))

  ## Each value is the max position within the chromosome
  manual_vals <- de_res[["chromosome"]] |>
    unique() |>
    sort() |>
    purrr::map_dbl(~ max(c(max(de_res$start[de_res$chromosome == .x]),
                           max(de_res$end[de_res$chromosome == .x]))))

  expect_equal(chr_length[["length"]],
               manual_vals)

})

test_that("compute_chrom_length.CAN_data works", {

  chr_length <- compute_chrom_length(can_res)

  ## One value per chromosome
  expect_equal(chr_length[["chromosome"]], sort(unique(can_res[["chromosome"]])))

  ## Each value is the max position within the chromosome
  manual_vals <- can_res[["chromosome"]] |>
    unique() |>
    sort() |>
    purrr::map_dbl(~ max(c(max(can_res$start[can_res$chromosome == .x]),
                           max(can_res$end[can_res$chromosome == .x]))))

  expect_equal(chr_length[["length"]],
               manual_vals)

})
