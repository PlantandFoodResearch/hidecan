#' Launch HIDECAN shiny app
#'
#' Starts the HIDECAN shiny app. The app reads in csv data to produce
#' a HIDECAN plot.
#'
#' @export
run_hidecan_shiny <- function(){
  app_dir <- system.file("shiny_apps/hidecan_shiny.R", package = "hidecan")

  shiny::runApp(app_dir, display.mode = "normal")
}
