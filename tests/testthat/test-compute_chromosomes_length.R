test_that("compute_chrom_length.GWAS_data works", {
  gwas_res <- test_get_gwas()

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
  de_res <- test_get_de()

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
  can_res <- test_get_can()

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


test_that("compute_chrom_length.QTL_data works", {
  qtl_res <- test_get_qtl()

  chr_length <- compute_chrom_length(qtl_res)

  ## One value per chromosome
  expect_equal(chr_length[["chromosome"]], sort(unique(qtl_res[["chromosome"]])))

  ## Each value is the max position within the chromosome
  manual_vals <- qtl_res[["chromosome"]] |>
    unique() |>
    sort() |>
    purrr::map_dbl(~ max(c(max(qtl_res$start[qtl_res$chromosome == .x]),
                           max(qtl_res$end[qtl_res$chromosome == .x]))))

  expect_equal(chr_length[["length"]],
               manual_vals)

})

test_that("compute_chrom_length.CUSTOM_data works", {
  custom_res <- test_get_custom()

  chr_length <- compute_chrom_length(custom_res)

  ## One value per chromosome
  expect_equal(chr_length[["chromosome"]], sort(unique(custom_res[["chromosome"]])))

  ## Each value is the max position within the chromosome
  manual_vals <- custom_res[["chromosome"]] |>
    unique() |>
    sort() |>
    purrr::map_dbl(~ max(custom_res$position[custom_res$chromosome == .x]))

  expect_equal(chr_length[["length"]],
               manual_vals)

})

test_that("combine_chrom_length works", {
  gwas_res <- test_get_gwas()
  de_res <- test_get_de()
  can_res <- test_get_can()

  res <- combine_chrom_length(list(gwas_res, de_res, can_res))

  expect_equal(res[["chromosome"]],
               c(gwas_res[["chromosome"]], de_res[["chromosome"]], can_res[["chromosome"]]) |>
                 unique() |>
                 sort())

  expect_equal(res[["length"]],
               c(gwas_res[["chromosome"]], de_res[["chromosome"]], can_res[["chromosome"]]) |>
                 unique() |>
                 sort() |>
                 purrr::map_dbl(function(x){
                  max(
                    gwas_res[["position"]][gwas_res[["chromosome"]] == x],
                    de_res[["start"]][de_res[["chromosome"]] == x],
                    de_res[["end"]][de_res[["chromosome"]] == x],
                    can_res[["start"]][can_res[["chromosome"]] == x],
                    can_res[["end"]][can_res[["chromosome"]] == x]
                  )
                 }))
})
