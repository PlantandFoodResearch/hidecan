load_all()

gwas_res <- GWAS_data(gwas_data)
de_res <- DE_data(de_data)
can_res <- CAN_data(candidate_data)

# test <- list("GWAS" = gwas_res,
#      "DE" = de_res,
#      "CAN" = can_res) |>
#   purrr::map(compute_chrom_length)
#
# purrr::map(test, class)

chrom_length <- list("GWAS" = gwas_res,
     "DE" = de_res,
     "CAN" = can_res) |>
  combine_chrom_length()


gwas_res_thr <- apply_threshold(gwas_res, score_thr = 4)
de_res_thr <- apply_threshold(de_res, score_thr = -log10(0.05), log2fc_thr = 0)
can_res_thr <- apply_threshold(can_res)

x <- list(gwas_res_thr, de_res_thr, can_res_thr)

create_hidecan_plot(x, chrom_length, legend_position = "left", label_size = 1)
