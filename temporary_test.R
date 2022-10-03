load_all()

load("R/sysdata.rda")


gwas_res <- GWAS_data(gwas_data)
de_res <- DE_data(de_data)
can_res <- CAN_data(candidate_data)

list("GWAS" = gwas_res,
     "DE" = de_res,
     "CAN" = can_res) |>
  purrr::map(compute_chrom_length)

list("GWAS" = gwas_res,
     "DE" = de_res,
     "CAN" = can_res) |>
  purrr::map(apply_threshold)
