library(dplyr)
load("R/sysdata.rda")

set.seed(834)
gwas_res <- hidecan::GWAS_data(gwas_data |>
                                 dplyr::slice_sample(n = 1000))
de_res <- hidecan::DE_data(de_data|>
                             dplyr::slice_sample(n = 1000))
can_res <- hidecan::CAN_data(candidate_data)

save(gwas_res, de_res, can_res, file = "tests/testthat/data-test/test_input.rda")
