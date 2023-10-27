library(hexSticker)
library(hidecan)
library(tibble)
library(ggplot2)

z <- list()

z[["GWAS"]] <- tibble(
  chromosome = "chromosome 1",
  position = c(0, 15, 80, 100),
  score = c(0.1, 5, 6.5, 0.1)
)

z[["DE"]] <- tibble(
  chromosome = "chromosome 1",
  start = c(0, 15.5, 40, 50, 52, 79, 100),
  end = start,
  score = c(0.1, 5, 4.5, 8, 6.5, 3.5, 0.1),
  log2FoldChange = 3
)

z[["CAN"]] <- tibble(
  chromosome = "chromosome 1",
  start = c(50.5, 79.5),
  end = start,
  name = c("", "")
)

hidecan_plot <- hidecan_plot(
  gwas_list = z$GWAS,
  de_list = z$DE,
  can_list = z$CAN,
  title = "HIDECAN plot",
  point_size = 1.5,
  label_size = 3.5
) +
  labs(
    title = NULL,
    x = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    strip.background.x = element_rect(fill = "#6cc3c6"),
    strip.text.x = element_text(
      margin = margin(0, 0, 0, 0),
      color = NA
    ),
    panel.grid = element_blank()
  ) +
  theme_transparent()

col_border <- "#399093"
col_background <- "#ceebec"
sticker(
  hidecan_plot,
  package = "hidecan",
  p_size = 20,
  p_color = col_border,
  s_x = 1,
  s_y = 0.8,
  s_width = 1.3,
  s_height = 0.7,
  h_color = col_border,
  h_fill = col_background,
  filename = "man/figures/logo.png"
)
