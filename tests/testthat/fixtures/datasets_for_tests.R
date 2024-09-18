library(dplyr)
load("R/sysdata.rda")

set.seed(834)
gwas_res <- hidecan::GWAS_data(gwas_data |> dplyr::slice_sample(n = 1000))
de_res <- hidecan::DE_data(de_data|> dplyr::slice_sample(n = 1000))
can_res <- hidecan::CAN_data(candidate_data)
custom_res <- gwas_data |>
  dplyr::slice_sample(n = 1000) |>
  dplyr::mutate(score = runif(1000, 0, 5)) |>
  hidecan::CUSTOM_data()

saveRDS(gwas_res, file = "tests/testthat/fixtures/test_gwas.rds")
saveRDS(de_res, file = "tests/testthat/fixtures/test_de.rds")
saveRDS(can_res, file = "tests/testthat/fixtures/test_can.rds")
saveRDS(custom_res, file = "tests/testthat/fixtures/test_custom.rds")

library(GWASpoly)
library(readr)
library(dplyr)

genofile <- system.file("extdata", "TableS1.csv", package = "GWASpoly")
phenofile <- system.file("extdata", "TableS2.csv", package = "GWASpoly")

set.seed(10)
geno_data <- read_csv(genofile) |>
  slice_sample(n = 500)

temp_file <- tempfile()
write_csv(geno_data, temp_file)

data <- read.GWASpoly(ploidy = 4,
                      pheno.file = phenofile,
                      geno.file = temp_file,
                      format = "ACGT",
                      n.traits = 13,
                      delim = ",")

data.original <- set.K(data,
                       LOCO = FALSE,
                       n.core = 2)


gwaspoly_res <- GWASpoly(data.original,
                         models = c("general", "additive"),
                         traits = c("tuber_eye_depth", "tuber_shape"),
                         n.core = 2)

gwaspoly_res_thr <- set.threshold(gwaspoly_res, method = "M.eff", level = 0.05)

## For users to access the example GWASpoly dataset without having it loaded
saveRDS(gwaspoly_res, file = "tests/testthat/fixtures/test_gwaspoly_res.rds")
saveRDS(gwaspoly_res_thr, file = "tests/testthat/fixtures/test_gwaspoly_res_thr.rds")
