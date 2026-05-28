#------------------------
# Section 0: Packages
#------------------------
library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(GenErlangFit)

#------------------------
# UI
#------------------------

ui <- navbarPage(

  title = "GenErlangFit",

  theme = bs_theme(
    version = 5,
    bootswatch = "minty",
    navbar_bg = "#90C0AE"
  ),

  header = tags$head(
    tags$style(HTML("

    .format-btn {
      font-size: 11px !important;
      padding: 4px 8px !important;
      margin-bottom: 8px;
      background-color: #F3969A !important;
      border-color: #F3969A !important;
      color: white !important;
    }

    .format-btn:hover {
      background-color: #e67c86 !important;
      border-color: #e67c86 !important;
      filter: brightness(90%);
      transition: 0.2s;
    }

    .format-box {
      margin-top: 10px;
      padding: 12px;
      border: 1px solid #ddd;
      border-radius: 8px;
      background: #f8f9fa;
    }

    .card {
      padding: 15px;
      margin-bottom: 15px;
      border-radius: 10px;
      border: 1px solid #ddd;
      background: #ffffff;
    }

    .card-title {
      font-weight: bold;
      font-size: 18px;
      margin-bottom: 10px;
    }

    .error-text {
      color: #d9534f;
      font-weight: bold;
      margin-top: 5px;
    }

  "))
  ),

  # =========================================================
  # 1. DATA ENTRY
  # =========================================================
  tabPanel(
    "Data Entry",

    sidebarLayout(
      sidebarPanel(

        h4("Enter Data"),

        actionButton(
          "toggle_format",
          "Show Required Format",
          class = "btn btn-primary format-btn"
        ),

        conditionalPanel(
          condition = "input.toggle_format % 2 == 1",
          div(
            class = "format-box",
            strong("Required Format"),
            HTML("
Upload CSV or manual input:<br><br>
A single column of non-negative integers.<br><br>
Example:<br>
1<br>2<br>3<br>4<br>5")
          )
        ),

        fileInput("file", "Upload CSV"),
        div(textOutput("csv_error"), class = "error-text"),

        br(),

        h4("Manual Data Entry"),
        textAreaInput("manual_data", "Manual Data Entry", rows = 5),
        actionButton("submit_manual", "Submit Manual Data"),
        div(textOutput("manual_error"), class = "error-text"),

        br(),

        actionButton("clear", "Clear Data")
      ),

      mainPanel(
        h3("Dataset"),
        plotOutput("data_histogram", height = "400px")
      )
    )
  ),

  # =========================================================
  # 2. COMPUTE FIT
  # =========================================================
  tabPanel(
    "Compute Fit",

    sidebarLayout(

      # =========================
      # SIDEBAR (3 CARDS)
      # =========================
      sidebarPanel(

        # -------------------------
        # CARD 1: FIT SETTINGS
        # -------------------------
        div(class = "card",

            div(class = "card-title", "Fit Settings"),

            radioButtons(
              "fit_type",
              "Select Fit Type",
              choices = c("Default", "Erlang", "Erlang-Exp"),
              selected = "Default"
            ),

            conditionalPanel(
              condition = "input.fit_type == 'Erlang'",
              numericInput("initial_k", "Initial Guess for K", value = NA),
              checkboxInput("find_smallest_erlang", "Find Erlang Smallest K", value = FALSE)
            ),

            conditionalPanel(
              condition = "input.fit_type == 'Erlang-Exp'",
              numericInput("initial_k_exp", "Initial Guess for K", value = NA),

              radioButtons(
                "search_type",
                "Search Type",
                choices = c("Fixed K", "Search over a Window")
              ),

              conditionalPanel(
                condition = "input.search_type == 'Search over a Window'",
                numericInput("window_size", "Window Size", value = NA)
              ),

              checkboxInput(
                "find_smallest_erlang_exp",
                "Find Erlang-Exp Smallest K",
                value = FALSE
              )
            ),

            uiOutput("run_fit_button")
        ),

        # -------------------------
        # CARD 2: COMPUTE GOF
        # -------------------------
        div(class = "card",

            div(class = "card-title", "Compute GOF"),

            radioButtons(
              "gof_mode",
              "GOF Mode",
              choices = c("Default", "User Selection")
            ),

            conditionalPanel(
              condition = "input.gof_mode == 'User Selection'",
              numericInput("alpha_value", "Alpha Value", value = NA),
              numericInput("num_bootstraps", "Bootstraps", value = NA),
              checkboxGroupInput(
                "test_statistics",
                "Test Statistics",
                choices = c(
                  "KS",
                  "AD",
                  "CvM",
                  "Likelihood Ratio"
                )
              )
            ),

            uiOutput("run_gof_button"),
            verbatimTextOutput("gof_output")
        ),

        # -------------------------
        # CARD 3: SMALLEST K
        # -------------------------
        div(class = "card",

            div(class = "card-title", "Smallest K"),

            checkboxInput(
              "find_smallest_erlang",
              "Erlang Smallest K",
              value = FALSE
            ),

            checkboxInput(
              "find_smallest_erlang_exp",
              "Erlang-Exp Smallest K",
              value = FALSE
            )
        )
      ),

      # =========================
      # MAIN PANEL
      # =========================
      mainPanel(

        h3("Model Fit Output"),
        plotOutput("fit_plot", height = "500px"),
        DTOutput("fit_table")
      )
    )
  ),

  # =========================================================
  # 3. DATA VIEW (UNCHANGED LOGIC SIMPLIFIED)
  # =========================================================
  tabPanel(
    "Data View",
    h3("Data Preview"),
    plotOutput("data_histogram", height = "500px")
  )
)

#------------------------
# SERVER (UNCHANGED LOGIC)
#------------------------

server <- function(input, output, session) {

  data <- reactiveVal(NULL)

  csv_error <- reactiveVal("")
  manual_error <- reactiveVal("")

  output$csv_error <- renderText(csv_error())
  output$manual_error <- renderText(manual_error())

  observeEvent(input$file, {

    req(input$file)
    csv_error("")

    df <- read.csv(input$file$datapath, header = FALSE)

    colnames(df) <- "Value"
    data(df)
  })

  observeEvent(input$submit_manual, {

    manual_error("")

    tokens <- unlist(strsplit(input$manual_data, "[,\n\r]+"))
    tokens <- trimws(tokens)
    tokens <- tokens[tokens != ""]

    if (!all(grepl("^\\d+$", tokens))) {
      manual_error("Invalid input")
      return()
    }

    data(data.frame(Value = as.numeric(tokens)))
  })

  observeEvent(input$clear, {
    data(NULL)
  })

  output$data_histogram <- renderPlot({

    req(data())

    ggplot(data(), aes(Value)) +
      geom_histogram(fill = "#90C0AE", color = "black") +
      theme_minimal()
  })

  fit_results <- eventReactive(input$run_fit, {

    req(data())
    empiricaldata <- data()[[1]]

    if (input$fit_type == "Default") {
      GenErlang_Fit("QuickFitAllModels", empiricaldata)
    } else {
      GenErlang_Fit("Erlang", empiricaldata)
    }
  })

  output$run_fit_button <- renderUI({
    actionButton("run_fit", "Run Fit", class = "btn-primary")
  })

  output$fit_table <- renderDT({
    req(fit_results())
    fit_results()$ResultsTable
  })

  output$fit_plot <- renderPlot({
    req(data())
    ggplot(data(), aes(Value)) +
      geom_histogram(fill = "#90C0AE", color = "black") +
      theme_minimal()
  })

  output$gof_output <- renderPrint({
    req(input$run_gof)
    cat("GOF Mode:", input$gof_mode)
  })

  output$run_gof_button <- renderUI({
    actionButton("run_gof", "Compute GOF", class = "btn-primary")
  })
}

#------------------------
# RUN APP
#------------------------
shinyApp(ui, server)
