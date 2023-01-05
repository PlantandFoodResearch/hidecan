library(shiny)

################################################################################
##                           HELPER FUNCTIONS                                 ##
################################################################################

process_upload_data <- function(input_list, func){

  ## Extract files extensions
  extensions <- tools::file_ext(input_list$name)

  ## Safe version of the function
  safe_func <- purrr::safely(func)

  ## Read each file
  res <- purrr::map2(
    extensions,
    input_list$datapath,
    ~ switch(.x,
             csv = vroom::vroom(.y, show_col_types = FALSE, progress = FALSE) |> safe_func(),
             validate("Invalid file; Please upload a .csv file."))

  ) |>
    purrr::transpose()

  ## Checking whether any file returned an error
  no_error <- res$error |>
    purrr::map_lgl(is.null) |>
    all()

  if(no_error){
    res$error <- NULL
  } else {

    ## Extract error message
    error_msg <- res$error |>
      setNames(input_list$name) |>
      purrr::map(purrr::pluck, "message") |>
      purrr::imap(~ paste0("Input file ", .y, ": ", .x)) |>
      purrr::reduce(paste0, collapse = "\n")

    res$error <- error_msg

    ## Remove NULL elements from the list or results
    ## If all are NULL, will return an empty list()
    res$result <- purrr::discard(res$result,
                                 is.null)
  }

  return(res)
}


get_input_data_reactive <- function(input_id, func){

  substitute(reactive({

    if(!is.null(input$input_id)){
      ## Read in the data
      res <- process_upload_data(input$input_id, func)

      ## Print error message for all files that are not in
      ## the right format
      shinyFeedback::feedbackDanger(input_id,
                                    !is.null(res$error),
                                    res$error)

      res$result
    } else{
      list()
    }
  }))
}


render_custom_names_textInput <- function(input_id, suffix){
  substitute(
    renderUI({
      req(input$input_id)

      input_files_df <- input$input_id

      purrr::map(seq_len(nrow(input_files_df)),
                 ~ textInput(inputId = paste0("custom_name_", suffix, "_", .x),
                             label = input_files_df$name[.x],
                             value = " "))

    })

  )
}


gwas_threshold_tabs <- tabsetPanel(
  id = "thr_selection_gwas",
  type = "hidden",
  tabPanel(
    "gwas_thr",
    sliderInput("score_thr_gwas", "Score threshold for GWAS results", value = 4, min = 0, max = 10, step = 0.1)
  )
)


de_threshold_tabs <- tabsetPanel(
  id = "thr_selection_de",
  type = "hidden",
  tabPanel(
    "de_thr",
    sliderInput("score_thr_de", "Score threshold for DE results", value = 1.3, min = 0, max = 10, step = 0.1),
    sliderInput("log2fc_thr_de", "log2(fold-change) threshold for DE results", value = 1, min = 0, max = 10, step = 0.1)
  )
)


# download_params <- tabsetPanel(
#   id = "download_param",
#   type = "hidden",
#   tabPanel(
#     "download_param",
#     tags$div(selectInput("download_format", "Format", c("PNG", "PDF")),  style="display:inline-block; width: 20%; margin-bottom: 0px"),
#     tags$div(numericInput("download_plot_width", "Width (in)", value = 16, min = 0, max = Inf),  style="display:inline-block; width: 20%"),
#     tags$div(numericInput("download_plot_height", "Height (in)", value = 10, min = 0, max = Inf),  style="display:inline-block; width: 20%"),
#     tags$div(downloadButton("download_plot"),  style="display:inline-block")
#   )
# )

################################################################################
##                                  UI                                        ##
################################################################################

ui <- fluidPage(

  ## User feedback
  shinyFeedback::useShinyFeedback(),

  titlePanel(
    "HIDECAN plot"
  ),
  sidebarLayout(

    ## ---------------------------------------------------------------------- ##
    ##                              Inputs                                    ##
    ## ---------------------------------------------------------------------- ##

    sidebarPanel(

      tabsetPanel(
        id = "sidebar_panel",

        tabPanel(
          title = "Import data",

          ## Select GWAS results
          fileInput("upload_gwas",
                    "Upload one or more GWAS results file",
                    accept = ".csv",
                    multiple = TRUE),


          ## Select DE results
          fileInput("upload_de",
                    "Upload one or more DE results file",
                    accept = ".csv",
                    multiple = TRUE),


          ## Select CAN results
          fileInput("upload_can",
                    "Upload one or more candidate genes list file",
                    accept = ".csv",
                    multiple = TRUE),

          ## Choose custom names for the different input files
          uiOutput("input_custom_names_gwas"),
          uiOutput("input_custom_names_de"),
          uiOutput("input_custom_names_can")
        ),

        tabPanel(
          title = "Thresholds",

          ## Input sliders for GWAS score threshold
          gwas_threshold_tabs,

          ## Input sliders for DE score and log2FC threshold
          de_threshold_tabs
        ),

        tabPanel(
          title = "Plot options",

          textInput("title", "Title"),
          textInput("subtitle", "Subtitle"),
          fluidRow(
            column(6,
                   numericInput("nrows", "Number of rows", value = NULL, min = 1, max = Inf)),
            column(6,
                   numericInput("ncols", "Number of columns", value = 2, min = 1, max = Inf))
          ),
          selectInput("legend_position", "Legend position", c("bottom", "top", "left", "right", "none")),
          fluidRow(
            column(6,
                   numericInput("point_size", "Point size", value = 3, min = 0, max = Inf))
          ),
          fluidRow(
            column(6,
                   numericInput("label_size", "Label size", value = 3.5, min = 0, max = Inf)),
            column(6,
                   numericInput("label_padding", "Label padding", value = 0.15, min = 0, max = Inf))
          ),
          checkboxInput("colour_genes_by_score", "Colour genes by score?", value = TRUE),
          checkboxInput("remove_empty_chrom", "Remove empty chromosomes?", value = TRUE)
        ),

        tabPanel(
          title = "Download",

          uiOutput("download_params")
        )


      ),

      #uiOutput("download_params")
    ),

    ## ---------------------------------------------------------------------- ##
    ##                              Outputs                                   ##
    ## ---------------------------------------------------------------------- ##

    mainPanel(

      # ## Preview results
      # verbatimTextOutput("head_all"),
      # tableOutput("chrom_length")

      plotOutput("hidecan_plot", height = "800px")

    )
  )
)


################################################################################
##                                 SERVER                                     ##
################################################################################


server <- function(input, output, session){

  ## ------------------------------------------------------------------------ ##
  ##                               Reactives                                  ##
  ## ------------------------------------------------------------------------ ##

  ## Reading in gwas_data
  gwas_data_list <- eval(get_input_data_reactive("upload_gwas", GWAS_data))

  ## Reading in de_data
  de_data_list <- eval(get_input_data_reactive("upload_de", DE_data))

  ## Reading in can_data
  can_data_list <- eval(get_input_data_reactive("upload_can", CAN_data))


  ## Adding the buttons for the custom names
  output$input_custom_names_gwas <- eval(render_custom_names_textInput("upload_gwas", "gwas"))
  output$input_custom_names_de <- eval(render_custom_names_textInput("upload_de", "de"))
  output$input_custom_names_can <- eval(render_custom_names_textInput("upload_can", "can"))

  get_custom_names <- reactive({
    c(
      purrr::map(seq_along(gwas_data_list()),
                 ~ eval(str2expression(paste0("input$custom_name_gwas_", .x)))),
      purrr::map(seq_along(de_data_list()),
                 ~ eval(str2expression(paste0("input$custom_name_de_", .x)))),
      purrr::map(seq_along(can_data_list()),
                 ~ eval(str2expression(paste0("input$custom_name_can_", .x))))
    )
  })

  ## Getting all the data
  all_data_list <- reactive({
    c(gwas_data_list(),
      de_data_list(),
      can_data_list())
  })

  ## Checking whether the threshold inputs should be displayed
  observeEvent(gwas_data_list(), {
    if(length(gwas_data_list())){
      showTab(inputId = "thr_selection_gwas", target = "gwas_thr", select = TRUE)
    } else {
      hideTab(inputId = "thr_selection_gwas", target = "gwas_thr")
    }
  })

  observeEvent(de_data_list(), {
    if(length(de_data_list())){
      showTab(inputId = "thr_selection_de", target = "de_thr", select = TRUE)
    } else {
      hideTab(inputId = "thr_selection_de", target = "de_thr")
    }
  })

  # observeEvent(any_data(), {
  #   if(any_data()){
  #     showTab(inputId = "download_param", target = "download_param", select = TRUE)
  #   } else {
  #     hideTab(inputId = "download_param", target = "download_param")
  #   }
  # })


  ## Check whether there is any input data
  any_data <- reactive({
    length(all_data_list()) > 0
  })

  chrom_length <- reactive({

    ## wait for at list one input dataset
    req(any_data())

    all_data_list() |>
      combine_chrom_length()
  })

  plot_nrows <- reactive({
    res <- input$nrows
    print(res)
    if(missing(res) | is.na(res)) res <- NULL

    res
  })

  hidecan_plot <- reactive({
    ## wait for at list one input dataset
    req(any_data())

    x <- c(
      gwas_data_list() |>
        purrr::map(apply_threshold, input$score_thr_gwas),
      de_data_list() |>
        purrr::map(apply_threshold, input$score_thr_de, input$log2fc_thr_de),
      can_data_list() |>
        purrr::map(apply_threshold)
    ) |>
      setNames(get_custom_names())

    create_hidecan_plot(x,
                        chrom_length(),
                        colour_genes_by_score = input$colour_genes_by_score,
                        remove_empty_chrom = input$remove_empty_chrom,
                        title = input$title,
                        subtitle = input$subtitle,
                        n_rows = plot_nrows(),
                        n_cols = input$ncols,
                        legend_position = input$legend_position,
                        point_size = input$point_size,
                        label_size = input$label_size,
                        label_padding = input$label_padding)
  })

  ## Display the download button
  output$download_params <- renderUI({
    req(any_data())

    #h3("Download"),

    ## Parameters for downloading the data
    #download_params
    fluidRow(
      column(
        3,
        selectInput("download_format", "Format", c("png", "pdf"))
      ),
      column(
        3,
        numericInput("download_plot_width", "Width (in)", value = 16, min = 0, max = Inf)
      ),
      column(
        3,
        numericInput("download_plot_height", "Height (in)", value = 10, min = 0, max = Inf)
      ),
      column(
        3,
        style = "margin-top: 25px;",
        downloadButton("download_plot")
      )

    )

  })

  ## ------------------------------------------------------------------------ ##
  ##                                Outputs                                   ##
  ## ------------------------------------------------------------------------ ##

  output$head_all <- renderPrint({
    c(gwas_data_list(),
      de_data_list(),
      can_data_list()) |>
      purrr::map(head)
  })

  output$chrom_length <- renderTable(
    chrom_length()
  )

  output$hidecan_plot <- renderPlot({
    hidecan_plot()

  })

  output$download_plot <- downloadHandler(
    filename = function(){
      paste0("hidecan_plot_", Sys.Date(), ".", input$download_format)
    },
    content = function(file){
      ggplot2::ggsave(file,
                      hidecan_plot(),
                      device = input$download_format,
                      width = input$download_plot_width,
                      height = input$download_plot_height)
    }
  )

}

shinyApp(ui, server)
