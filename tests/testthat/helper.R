test_get_gwas <- function() {
  readRDS(test_path("fixtures", "test_gwas.rds"))
}

test_get_de <- function() {
  readRDS(test_path("fixtures", "test_de.rds"))
}

test_get_can <- function() {
  readRDS(test_path("fixtures", "test_can.rds"))
}

test_get_qtl <- function() {
  readRDS(test_path("fixtures", "test_qtl.rds"))
}

test_get_custom <- function() {
  readRDS(test_path("fixtures", "test_custom.rds"))
}

test_get_gwaspoly <- function() {
  readRDS(test_path("fixtures", "test_gwaspoly_res.rds"))
}

test_get_gwaspoly_thr <- function() {
  readRDS(test_path("fixtures", "test_gwaspoly_res_thr.rds"))
}

get_chrom_limits_plot <- function(p){
  p$layers[[1]]$data |>
    dplyr::select(chromosome, position_mb, position_mb_end) |>
    dplyr::distinct()
}
