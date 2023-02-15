#' Launches the HIDECAN shiny app
#'
#' Starts the HIDECAN shiny app. The app reads in csv data to produce
#' a HIDECAN plot.
#'
#' @returns No return value, called for side effects (launching the shiny app).
#'
#' @export
run_hidecan_shiny <- function(){
  app_dir <- system.file("shiny_apps/hidecan_shiny.R", package = "hidecan")

  shiny::runApp(app_dir, display.mode = "normal")
}

## To get rid of the check NOTE
ignore_unused_imports <- function(){
  vroom::vroom
}
