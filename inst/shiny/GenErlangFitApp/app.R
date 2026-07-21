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


  # =========================================================
  # ADD CSS (hover + small button styling) - Just a custom UI Feature.
  # =========================================================
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

    .smallest-k-box {
      margin-top: 15px;
      padding: 12px;
      border: 2px solid #78C2AD;
      border-radius: 8px;
      background: #f0f9f7;
    }

  "))
  ),

  # =========================================================
  # 0. ABOUT
  # =========================================================
  tabPanel(
    "About",

    fluidPage(
      h2("About GenErlangFit"),

      br(),

      p("Welcome to the GenErlangFit Shiny application."),

      p("This app provides an interactive interface for fitting Erlang and
        Erlang-Exponential distributions to empirical data using the GenErlangFit R package."),

      br(),

      h4("Features"),
      tags$ul(
        tags$li("Upload or manually enter data"),
        tags$li("Fit Erlang and Erlang-Exp distributions"),
        tags$li("Compute goodness-of-fit tests with bootstrap p-values"),
        tags$li("Find the smallest K that passes the goodness-of-fit test")
      ),

      br(),

      h4("Resources"),
      p(
        "Package documentation: ",
        tags$a(
          href = "https://vivekmurals.github.io/GenErlangFit/",
          target = "_blank",
          "https://vivekmurals.github.io/GenErlangFit/"
        )
      ),
      p(
        "GitHub repository: ",
        tags$a(
          href = "https://github.com/vivekmurals/GenErlangFit",
          target = "_blank",
          "https://github.com/vivekmurals/GenErlangFit"
        )
      )
    )
  ),

  # =========================================================
  # 1. DATA ENTRY
  # =========================================================
  tabPanel(
    "Data Entry",

    sidebarLayout(

      # -----------------------------------------------------
      # SIDEBAR
      # -----------------------------------------------------
      sidebarPanel(

        h4("Enter Data"),
        # =====================================================
        # FORMAT TOGGLE (MOVED UP + SMALL BUTTON)
        # =====================================================
        actionButton(
          "toggle_format",
          "Show Required Format",
          class = "btn btn-primary format-btn"
        ),

        # FORMAT BOX
        conditionalPanel(
          condition = "input.toggle_format % 2 == 1",

          div(
            class = "format-box",

            strong("Required Format"),

            div(
              HTML("
Upload a CSV file or manually enter data in the following format:<br>

A single column of non-negative integers (no header) and each row contains one value.<br><br>

Example CSV / Manual Entry:<br>
1<br>
2<br>
3<br>
4<br>
5<br>
6<br>
7<br>
8<br>
9
")
            )
          )
        ),



        fileInput(
          "file",
          "Upload CSV"
        ),
        textOutput("csv_error"),
        br(),

        h4("Manual Data Entry"),

        textAreaInput(
          "manual_data",
          "Paste or type data here:",
          rows = 5
        ),

        actionButton(
          "submit_manual",
          "Submit Manual Data"
        ),

        textOutput("manual_error"),

        br(), br(),

        actionButton(
          "clear",
          "Clear Data"
        ),

      ),

      # -----------------------------------------------------
      # MAIN PANEL
      # -----------------------------------------------------
      mainPanel(

        tabsetPanel(

          # ONLY ONE TAB NOW
          tabPanel(
            "View Data",

            h3("Uploaded / Entered Dataset"),

            br(),

            plotOutput(
              "data_histogram",
              height = "500px"
            )
          )

        )
      )
    )
  ),


  # =========================================================
  # 2. COMPUTE FIT
  # =========================================================
  tabPanel(
    "Compute Fit",

    sidebarLayout(

      # -----------------------------------------------------
      # SIDEBAR
      # -----------------------------------------------------
      sidebarPanel(

        h4("Fit Settings"),


        # ----------------------------------------------
        # FIT TYPE
        # ----------------------------------------------
        radioButtons(
          "fit_type",
          "Select Fit Type",
          choices = c(
            "Default",
            "Erlang",
            "Erlang-Exp"
          ),
          selected = "Default"
        ),


        # =================================================
        # ERLANG
        # =================================================
        conditionalPanel(
          condition = "
            input.fit_type == 'Erlang'
          ",


          # ------------------------------------------
          # OPTIONAL K
          # ------------------------------------------
          numericInput(
            "initial_k",
            "Initial Guess for K (optional)",
            value = NA,
            min = 1,
            step = 1
          ),


          # ------------------------------------------
          # OPTIONAL SMALLEST K
          # ------------------------------------------
          div(
            style = "display: flex; align-items: flex-start; gap: 5px;",
            checkboxInput(
              "find_smallest_erlang",
              "Find Erlang Smallest K",
              value = FALSE
            ),
            tooltip(
              icon("info-circle", class = "text-muted"),
              "Quick check to identify the smallest K consistent with the data based on Kolmogorov-Smirnov test statistic (α = 0.050). Provides immediate exploratory insight into whether simpler models may be adequate. For alternative test statistics, see Goodness of Fit tab.",
              placement = "right"
            )
          )
        ),


        # =================================================
        # ERLANG-EXP
        # =================================================
        conditionalPanel(
          condition = "
            input.fit_type == 'Erlang-Exp'
          ",


          # ------------------------------------------
          # REQUIRED K
          # ------------------------------------------
          numericInput(
            "initial_k_exp",
            "Initial Guess for K",
            value = NA,
            min = 1,
            step = 1
          ),


          br(),


          # ------------------------------------------
          # SEARCH TYPE
          # ------------------------------------------
          radioButtons(
            "search_type",
            "Search Type",
            choices = c(
              "Fixed K",
              "Search over a Window"
            ),
            selected = character(0)
          ),


          # ------------------------------------------
          # WINDOW SIZE
          # ------------------------------------------
          conditionalPanel(
            condition = "
              input.search_type == 'Search over a Window'
            ",

            numericInput(
              "window_size",
              "Window Size",
              value = NA,
              min = 1,
              step = 1
            )
          ),


          br(),


          # ------------------------------------------
          # OPTIONAL SMALLEST K
          # ------------------------------------------
          div(
            style = "display: flex; align-items: flex-start; gap: 5px;",
            checkboxInput(
              "find_smallest_erlang_exp",
              "Find Erlang-Exp Smallest K",
              value = FALSE
            ),
            tooltip(
              icon("info-circle", class = "text-muted"),
              "Quick check to identify the smallest K consistent with the based on Kolmogorov-Smirnov test statistic (α = 0.050). Provides immediate exploratory insight into whether simpler models may be adequate. For alternative test statistics, see Goodness of Fit tab.",
              placement = "right"
            )
          )
        ),



        br(),


        # ----------------------------------------------
        # RUN FIT BUTTON
        # ----------------------------------------------
        uiOutput("run_fit_button")
      ),


      # -----------------------------------------------------
      # MAIN PANEL
      # -----------------------------------------------------
      mainPanel(

        h3("Model Fit Output"),

        br(),

        plotOutput(
          "fit_plot",
          height = "550px"
        ),

        br(),

        DT::DTOutput("fit_table")
      )
    )
  ),


  # =========================================================
  # 3. GOODNESS OF FIT
  # =========================================================
  tabPanel(
    "Compute Goodness of Fit",

    sidebarLayout(

      # -----------------------------------------------------
      # SIDEBAR
      # -----------------------------------------------------
      sidebarPanel(

        h4("GOF Options"),

        # ---------------------------------------------------
        # DISPLAY CURRENT FIT INFO
        # ---------------------------------------------------
        uiOutput("current_fit_info"),

        hr(),

        # ---------------------------------------------------
        # GOF MODE
        # ---------------------------------------------------
        radioButtons(
          "gof_mode",
          "Select GOF Mode",
          choices = c(
            "Default",
            "User Selection"
          ),
          selected = "Default"
        ),


        # ===================================================
        # USER SELECTION OPTIONS
        # ===================================================
        conditionalPanel(
          condition = "
            input.gof_mode == 'User Selection'
          ",


          # -----------------------------------------------
          # ALPHA VALUE
          # -----------------------------------------------
          numericInput(
            "alpha_value",
            "Alpha Value",
            value = 0.05,
            min = 0.001,
            max = 1,
            step = 0.01
          ),


          # -----------------------------------------------
          # NUMBER OF BOOTSTRAPS
          # -----------------------------------------------
          numericInput(
            "num_bootstraps",
            "Number of Bootstraps",
            value = 200,
            min = 10,
            step = 10
          ),


          # -----------------------------------------------
          # TEST STATISTICS
          # -----------------------------------------------
          radioButtons(
            "test_statistic",
            "Choice of Test Statistic",
            choices = c(
              "KS" = "KS",
              "Anderson-Darling" = "AD",
              "Cramer-von Mises" = "CvM"
            ),
            selected = "KS"
          )
        ),


        br(),


        # ---------------------------------------------------
        # COMPUTE GOF BUTTON
        # ---------------------------------------------------
        uiOutput("run_gof_button"),


        # ---------------------------------------------------
        # SMALLEST K SECTION (appears after GOF is run)
        # ---------------------------------------------------
        uiOutput("smallest_k_section")
      ),


      # -----------------------------------------------------
      # MAIN PANEL
      # -----------------------------------------------------
      mainPanel(

        h3("Goodness of Fit Results"),

        br(),

        # ---------------------------------------------------
        # DYNAMIC GOF OUTPUT
        # ---------------------------------------------------
        uiOutput("gof_main_output")
      )
    )
  )
)



#------------------------
# SERVER
#------------------------

server <- function(input, output, session) {

  # =========================================================
  # DATA STORAGE
  # =========================================================

  data <- reactiveVal(NULL)


  # =========================================================
  # CSV UPLOAD
  # =========================================================

  observeEvent(input$file, {

    req(input$file)

    uploaded_data <- read.csv(
      input$file$datapath,
      header = FALSE
    )

    colnames(uploaded_data) <- "Value"

    data(uploaded_data)
  })


  # =========================================================
  # MANUAL UPLOAD AND ERROR
  # =========================================================
  manual_error <- reactiveVal("")

  output$manual_error <- renderText({
    manual_error()
  })

  observeEvent(input$submit_manual, {

    req(input$manual_data)

    raw <- input$manual_data

    manual_error("")

    # -------------------------------------------------------
    # Split by commas OR newlines
    # -------------------------------------------------------
    tokens <- unlist(strsplit(raw, "[,\n\r]+"))

    tokens <- trimws(tokens)
    tokens <- tokens[tokens != ""]

    # -------------------------------------------------------
    # Validate format
    # -------------------------------------------------------
    is_valid <- all(grepl("^\\d+$", tokens))

    if (!is_valid || length(tokens) == 0) {

      manual_error("Input data does not match the required format.")
      return(NULL)
    }

    values <- as.numeric(tokens)

    df <- data.frame(Value = values)

    data(df)

    manual_error("")  # clear error on success
  })



  # =========================================================
  # FORMAT TOGGLE BUTTON
  # =========================================================

  format_state <- reactiveVal(FALSE)

  observeEvent(input$toggle_format, {

    format_state(!format_state())

    updateActionButton(
      session,
      "toggle_format",
      label = if (format_state()) {
        "Hide Required Format"
      } else {
        "Show Required Format"
      }
    )
  })

  # =========================================================
  # FORMAT CLEAR BUTTON
  # =========================================================

  observeEvent(input$clear, {

    data(NULL)

  })


  # =========================================================
  # DYNAMIC RUN FIT BUTTON
  # =========================================================

  output$run_fit_button <- renderUI({

    valid <- FALSE


    # -------------------------------------------------------
    # DEFAULT
    # -------------------------------------------------------
    if (input$fit_type == "Default") {

      valid <- TRUE
    }


    # -------------------------------------------------------
    # ERLANG
    # -------------------------------------------------------
    if (input$fit_type == "Erlang") {

      valid <- TRUE
    }


    # -------------------------------------------------------
    # ERLANG-EXP
    # -------------------------------------------------------
    if (input$fit_type == "Erlang-Exp") {

      if (!is.na(input$initial_k_exp)) {

        if (input$search_type == "Fixed K") {

          valid <- TRUE
        }


        if (
          input$search_type == "Search over a Window" &&
          !is.na(input$window_size)
        ) {

          valid <- TRUE
        }
      }
    }


    # -------------------------------------------------------
    # BUTTON STATE
    # -------------------------------------------------------
    if (valid) {

      actionButton(
        "run_fit",
        "Run Fit",
        class = "btn-primary"
      )

    } else {

      actionButton(
        "run_fit",
        "Run Fit",
        class = "btn-primary disabled"
      )
    }
  })


  # =========================================================
  # FIT COMPUTATION
  # =========================================================

  fit_results <- eventReactive(input$run_fit, {

    req(data())

    empiricaldata <- data()[[1]]


    # =====================================================
    # DEFAULT
    # =====================================================
    if (input$fit_type == "Default") {

      results <- GenErlang_Fit("QuickFitAllModels", empiricaldata)
    }


    # =====================================================
    # ERLANG
    # =====================================================
    if (input$fit_type == "Erlang") {

      fit_args <- list(
        mode = "Erlang",
        empiricaldata = empiricaldata
      )

      fit_args$pvaloption = "nil" # By default, does not compute p-value.

      if (!is.na(input$initial_k)) {
        fit_args$K <- input$initial_k
      }


      if (isTRUE(input$find_smallest_erlang)) {
        fit_args$SmallestK <- TRUE
        fit_args$pvaloption = "KS"
      }



      results <- do.call(
        GenErlang_Fit,
        fit_args
      )
    }


    # =====================================================
    # ERLANG-EXP
    # =====================================================
    if (input$fit_type == "Erlang-Exp") {

      fit_args <- list(
        mode = "ErlangExp",
        empiricaldata = empiricaldata,
        K = input$initial_k_exp
      )

      fit_args$pvaloption = "nil" # By default, does not compute p-value.

      if (isTRUE(input$find_smallest_erlang_exp)) {
        fit_args$SmallestK <- TRUE
        fit_args$pvaloption = "KS"
      }


      if (input$search_type == "Fixed K") {

        fit_args$FixedK <- TRUE

      } else if (input$search_type == "Search over a Window") {

        fit_args$FixedK <- FALSE
        fit_args$KWindowSize <- input$window_size

      }


      results <- do.call(
        GenErlang_Fit,
        fit_args
      )
    }


    return(results)
  })


  # =========================================================
  # FIT OUTPUT
  # =========================================================

  output$fit_table <- DT::renderDT({

    req(fit_results())

    # Get the results table
    fit_table_df <- fit_results()$ResultsTable

    # Round numeric columns to 4 decimal places
    numeric_cols <- sapply(fit_table_df, is.numeric)
    fit_table_df[numeric_cols] <- lapply(fit_table_df[numeric_cols], function(x) round(x, 4))

    # Rename columns based on what's present
    current_names <- colnames(fit_table_df)
    new_names <- current_names

    new_names[current_names == "ErlangLambda"] <- "Erlang λ"
    new_names[current_names == "ExpLambda"] <- "Exp λ"
    new_names[current_names == "Lambda"] <- "Erlang λ"
    new_names[current_names == "LogLikelihood"] <- "Neg LogLikelihood"

    colnames(fit_table_df) <- new_names

    # Replace 0 with NA in Exp λ column (shows as blank)
    if ("Exp λ" %in% colnames(fit_table_df)) {
      fit_table_df[["Exp λ"]][fit_table_df[["Exp λ"]] == 0] <- NA
    }

    fit_table_df
  },
  options = list(
    pageLength = 5,
    searching = FALSE,
    lengthChange = FALSE,
    dom = "t"
  ))


  # =========================================================
  # DYNAMIC GOF BUTTON
  # =========================================================

  output$run_gof_button <- renderUI({

    # Check if fit results exist
    fit_exists <- !is.null(tryCatch(fit_results(), error = function(e) NULL))

    valid <- FALSE

    # Check for Erlang, Erlang-Exp, or Default fit types
    if (fit_exists && input$fit_type %in% c("Erlang", "Erlang-Exp", "Default")) {

      if (input$gof_mode == "Default") {
        valid <- TRUE
      }

      if (input$gof_mode == "User Selection") {
        if (
          !is.na(input$alpha_value) &&
          !is.na(input$num_bootstraps) &&
          input$alpha_value > 0 &&
          input$alpha_value <= 1 &&
          input$num_bootstraps >= 10
        ) {
          valid <- TRUE
        }
      }
    }


    if (valid) {

      actionButton(
        "run_gof",
        "Compute GOF",
        class = "btn-primary"
      )

    } else {

      actionButton(
        "run_gof",
        "Compute GOF",
        class = "btn-primary disabled"
      )
    }
  })


  # =========================================================
  # CURRENT FIT INFO DISPLAY
  # =========================================================

  output$current_fit_info <- renderUI({

    fit_exists <- !is.null(tryCatch(fit_results(), error = function(e) NULL))

    if (!fit_exists) {
      return(
        div(
          style = "color: #856404; background-color: #fff3cd; padding: 10px; border-radius: 5px;",
          icon("exclamation-triangle"),
          " No fit computed yet. Please run a fit first."
        )
      )
    }

    # Get fit type
    fit_type <- input$fit_type

    if (fit_type == "Erlang") {

      k_star <- fit_results()$Best$K_star
      lambda_star <- fit_results()$Best$Lambda_star

      return(
        div(
          style = "background-color: #d4edda; padding: 10px; border-radius: 5px;",
          strong("Current Fit: Erlang"),
          br(),
          sprintf("K* = %d, λ* = %.4f", k_star, lambda_star)
        )
      )

    } else if (fit_type == "Erlang-Exp") {

      k_star <- fit_results()$Best$K_star
      erlang_lambda_star <- fit_results()$Best$ErlangLambda_star
      exp_lambda_star <- fit_results()$Best$ExpLambda_star

      return(
        div(
          style = "background-color: #d4edda; padding: 10px; border-radius: 5px;",
          strong("Current Fit: Erlang-Exp"),
          br(),
          sprintf("K* = %d", k_star),
          br(),
          sprintf("Erlang λ* = %.4f", erlang_lambda_star),
          br(),
          sprintf("Exp λ* = %.4f", exp_lambda_star)
        )
      )

    } else if (fit_type == "Default") {

      # Get Erlang results
      k_erlang <- fit_results()$Erlang_Results$Best$K_star
      lambda_erlang <- fit_results()$Erlang_Results$Best$Lambda_star

      # Get Erlang-Exp results
      k_erlang_exp <- fit_results()$ErlangExp_Results$Best$K_star
      erlang_lambda_exp <- fit_results()$ErlangExp_Results$Best$ErlangLambda_star
      exp_lambda_exp <- fit_results()$ErlangExp_Results$Best$ExpLambda_star

      return(
        div(
          style = "background-color: #d4edda; padding: 10px; border-radius: 5px;",
          strong("Current Fit: Default (Both Models)"),
          br(), br(),
          tags$span(style = "color: red;", strong("Erlang:")),
          sprintf(" K* = %d, λ Erlang* = %.4f", k_erlang, lambda_erlang),
          br(),
          tags$span(style = "color: blue;", strong("Erlang-Exp:")),
          sprintf(" K* = %d, λ Erlang* = %.4f, λ Exp* = %.4f", k_erlang_exp, erlang_lambda_exp, exp_lambda_exp)
        )
      )
    }
  })


  # =========================================================
  # GOF COMPUTATION
  # =========================================================

  gof_results <- eventReactive(input$run_gof, {

    req(fit_results())
    req(data())

    empiricaldata <- data()[[1]]

    # =====================================================
    # DEFAULT GOF (Both Erlang and Erlang-Exp)
    # =====================================================
    if (input$fit_type == "Default") {

      # Get Erlang fitted parameters
      k_erlang <- fit_results()$Erlang_Results$Best$K_star
      lambda_erlang <- fit_results()$Erlang_Results$Best$Lambda_star

      # Get Erlang-Exp fitted parameters
      k_erlang_exp <- fit_results()$ErlangExp_Results$Best$K_star
      erlang_lambda_exp <- fit_results()$ErlangExp_Results$Best$ErlangLambda_star
      exp_lambda_exp <- fit_results()$ErlangExp_Results$Best$ExpLambda_star

      # Set GOF parameters based on mode
      if (input$gof_mode == "Default") {
        alpha <- 0.05
        n_bootstraps <- 200
        pvaloption <- "KS"
      } else {
        alpha <- input$alpha_value
        n_bootstraps <- input$num_bootstraps
        pvaloption <- input$test_statistic
      }

      # Compute Erlang GOF
      gof_erlang <- GenErlangFit:::Erlang_Fit_v2_Pvalue(
        empiricaldata = empiricaldata,
        k_star = k_erlang,
        lambda_star = 1/lambda_erlang,
        s = length(empiricaldata),
        n = n_bootstraps,
        alpha = alpha,
        pvaloption = pvaloption,
        ShowFigures = FALSE
      )

      # Compute Erlang-Exp GOF
      gof_erlang_exp <- GenErlangFit:::ErlangExp_Fit_v2_Pvalue(
        empiricaldata = empiricaldata,
        k_star = k_erlang_exp,
        erlambda_star = erlang_lambda_exp,
        explambda_star = exp_lambda_exp,
        s = length(empiricaldata),
        n = n_bootstraps,
        alpha = alpha,
        pvaloption = pvaloption,
        ShowFigures = FALSE
      )

      # Return combined results
      return(list(
        fit_type = "Default",
        alpha = alpha,
        n_bootstraps = n_bootstraps,
        pvaloption = pvaloption,
        empiricaldata = empiricaldata,
        # Erlang results
        erlang = list(
          k_star = k_erlang,
          lambda_star = lambda_erlang,
          p_value = gof_erlang$p_value,
          q_value = gof_erlang$q_value,
          metric_star = gof_erlang$metric_star,
          sample_stats = gof_erlang$sample_stats
        ),
        # Erlang-Exp results
        erlang_exp = list(
          k_star = k_erlang_exp,
          erlang_lambda_star = erlang_lambda_exp,
          exp_lambda_star = exp_lambda_exp,
          p_value = gof_erlang_exp$p_value,
          q_value = gof_erlang_exp$q_value,
          metric_star = gof_erlang_exp$metric_star,
          sample_stats = gof_erlang_exp$sample_stats
        )
      ))
    }

    # =====================================================
    # ERLANG GOF
    # =====================================================
    if (input$fit_type == "Erlang") {

      # Get fitted parameters
      k_star <- fit_results()$Best$K_star
      lambda_star <- fit_results()$Best$Lambda_star

      # Set GOF parameters based on mode
      if (input$gof_mode == "Default") {
        alpha <- 0.05
        n_bootstraps <- 200
        pvaloption <- "KS"
      } else {
        alpha <- input$alpha_value
        n_bootstraps <- input$num_bootstraps
        pvaloption <- input$test_statistic
      }

      # Call the GOF function
      gof_res <- GenErlangFit:::Erlang_Fit_v2_Pvalue(
        empiricaldata = empiricaldata,
        k_star = k_star,
        lambda_star = 1/lambda_star,
        s = length(empiricaldata),
        n = n_bootstraps,
        alpha = alpha,
        pvaloption = pvaloption,
        ShowFigures = FALSE
      )

      # Return results with metadata
      return(list(
        fit_type = "Erlang",
        k_star = k_star,
        lambda_star = lambda_star,
        erlang_lambda_star = NULL,
        exp_lambda_star = NULL,
        alpha = alpha,
        n_bootstraps = n_bootstraps,
        pvaloption = pvaloption,
        p_value = gof_res$p_value,
        q_value = gof_res$q_value,
        metric_star = gof_res$metric_star,
        sample_stats = gof_res$sample_stats,
        empiricaldata = empiricaldata
      ))
    }

    # =====================================================
    # ERLANG-EXP GOF
    # =====================================================
    if (input$fit_type == "Erlang-Exp") {

      # Get fitted parameters
      k_star <- fit_results()$Best$K_star
      erlang_lambda_star <- fit_results()$Best$ErlangLambda_star
      exp_lambda_star <- fit_results()$Best$ExpLambda_star

      # Set GOF parameters based on mode
      if (input$gof_mode == "Default") {
        alpha <- 0.05
        n_bootstraps <- 200
        pvaloption <- "KS"
      } else {
        alpha <- input$alpha_value
        n_bootstraps <- input$num_bootstraps
        pvaloption <- input$test_statistic
      }

      # Call the GOF function
      gof_res <- GenErlangFit:::ErlangExp_Fit_v2_Pvalue(
        empiricaldata = empiricaldata,
        k_star = k_star,
        erlambda_star = erlang_lambda_star,
        explambda_star = exp_lambda_star,
        s = length(empiricaldata),
        n = n_bootstraps,
        alpha = alpha,
        pvaloption = pvaloption,
        ShowFigures = FALSE
      )

      # Return results with metadata
      return(list(
        fit_type = "Erlang-Exp",
        k_star = k_star,
        lambda_star = NULL,
        erlang_lambda_star = erlang_lambda_star,
        exp_lambda_star = exp_lambda_star,
        alpha = alpha,
        n_bootstraps = n_bootstraps,
        pvaloption = pvaloption,
        p_value = gof_res$p_value,
        q_value = gof_res$q_value,
        metric_star = gof_res$metric_star,
        sample_stats = gof_res$sample_stats,
        empiricaldata = empiricaldata
      ))
    }

    return(NULL)
  })


  # =========================================================
  # GOF SUMMARY TABLE DATA (WITH LOG-LIKELIHOOD) - V1
  # =========================================================

  gof_summary_data <- reactive({

    # Check if GOF results exist
    gof_exists <- !is.null(tryCatch(gof_results(), error = function(e) NULL))

    if (!gof_exists) {
      return(NULL)
    }

    res <- gof_results()
    results <- list()
    empiricaldata <- res$empiricaldata

    # ---------------------------------------------------
    # Helper function to calculate Erlang log-likelihood
    # ---------------------------------------------------
    calc_erlang_loglik <- function(k, lambda, data) {
      if (is.null(k) || is.null(lambda) || is.na(k) || is.na(lambda)) return(NA_real_)
      if (k <= 0 || lambda <= 0) return(NA_real_)
      tryCatch({
        sum(dgamma(data, shape = k, scale = 1/lambda, log = TRUE))
      }, error = function(e) NA_real_)
    }

    # ---------------------------------------------------
    # Helper function to calculate Erlang-Exp log-likelihood
    # ---------------------------------------------------
    calc_erlangexp_loglik <- function(k, erlang_lambda, exp_lambda, data) {
      if (is.null(k) || is.null(erlang_lambda) || is.null(exp_lambda)) return(NA_real_)
      if (is.na(k) || is.na(erlang_lambda) || is.na(exp_lambda)) return(NA_real_)
      if (k <= 0 || erlang_lambda <= 0 || exp_lambda <= 0) return(NA_real_)
      tryCatch({
        res <- GenErlangFit:::ErlangExp_Func(data, ErK = k, Erlam = erlang_lambda, Explam = exp_lambda)
        res$Likelihood  # This is the log-likelihood (sum of log PDFs)
      }, error = function(e) NA_real_)
    }

    # =====================================================
    # DEFAULT MODE - Both models
    # =====================================================
    if (res$fit_type == "Default") {

      # Erlang results
      erl_loglik <- calc_erlang_loglik(res$erlang$k_star, res$erlang$lambda_star, empiricaldata)
      results$erlang <- data.frame(
        Model = "Erlang",
        K = res$erlang$k_star,
        `Erlang λ` = round(res$erlang$lambda_star, 4),
        `Exp λ` = NA,
        `Neg LogLikelihood` = if (!is.na(erl_loglik)) round(erl_loglik, 2) else NA_real_,
        `p-value` = round(res$erlang$p_value, 4),
        Result = ifelse(res$erlang$p_value >= res$alpha, "PASS", "FAIL"),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      # Erlang-Exp results
      erlexp_loglik <- calc_erlangexp_loglik(
        res$erlang_exp$k_star,
        res$erlang_exp$erlang_lambda_star,
        res$erlang_exp$exp_lambda_star,
        empiricaldata
      )
      results$erlang_exp <- data.frame(
        Model = "Erlang-Exp",
        K = res$erlang_exp$k_star,
        `Erlang λ` = round(res$erlang_exp$erlang_lambda_star, 4),
        `Exp λ` = round(res$erlang_exp$exp_lambda_star, 4),
        `Neg LogLikelihood` = if (!is.na(erlexp_loglik)) round(erlexp_loglik, 2) else NA_real_,
        `p-value` = round(res$erlang_exp$p_value, 4),
        Result = ifelse(res$erlang_exp$p_value >= res$alpha, "PASS", "FAIL"),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      # Check for Smallest K results (Default mode)
      smallest_k_res <- tryCatch(smallest_k_results_default(), error = function(e) NULL)


      if (isTRUE(smallest_k_res$erlang_selected) && !is.null(smallest_k_res$erlang)) {

        # Smallest K Erlang
        if (!is.null(smallest_k_res$erlang)) {
          sk_erl_loglik <- calc_erlang_loglik(
            smallest_k_res$erlang$smallest_k,
            smallest_k_res$erlang$smallest_lambda,
            empiricaldata
          )
          results$smallest_k_erlang <- data.frame(
            Model = "Erlang (Smallest K)",
            K = smallest_k_res$erlang$smallest_k,
            `Erlang λ` = round(smallest_k_res$erlang$smallest_lambda, 4),
            `Exp λ` = NA,
            `Neg LogLikelihood` = if (!is.na(sk_erl_loglik)) round(sk_erl_loglik, 2) else NA_real_,
            `p-value` = round(smallest_k_res$erlang$smallest_p_value, 4),
            Result = ifelse(smallest_k_res$erlang$smallest_q_value == 1, "PASS", "FAIL"),
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        }

        # Smallest K Erlang-Exp
        if (isTRUE(smallest_k_res$erlang_exp_selected) && !is.null(smallest_k_res$erlang_exp)){
          sk_erlexp_loglik <- calc_erlangexp_loglik(
            smallest_k_res$erlang_exp$smallest_k,
            smallest_k_res$erlang_exp$smallest_erlang_lambda,
            smallest_k_res$erlang_exp$smallest_exp_lambda,
            empiricaldata
          )
          results$smallest_k_erlang_exp <- data.frame(
            Model = "Erlang-Exp (Smallest K)",
            K = smallest_k_res$erlang_exp$smallest_k,
            `Erlang λ` = round(smallest_k_res$erlang_exp$smallest_erlang_lambda, 4),
            `Exp λ` = round(smallest_k_res$erlang_exp$smallest_exp_lambda, 4),
            `Neg LogLikelihood` = if (!is.na(sk_erlexp_loglik)) round(sk_erlexp_loglik, 2) else NA_real_,
            `p-value` = round(smallest_k_res$erlang_exp$smallest_p_value, 4),
            Result = ifelse(smallest_k_res$erlang_exp$smallest_q_value == 1, "PASS", "FAIL"),
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        }
      }

    } else {
      # =====================================================
      # SINGLE MODEL MODE
      # =====================================================

      if (res$fit_type == "Erlang") {
        erl_loglik <- calc_erlang_loglik(res$k_star, res$lambda_star, empiricaldata)
        results$single <- data.frame(
          Model = "Erlang",
          K = res$k_star,
          `Erlang λ` = round(res$lambda_star, 4),
          `Exp λ` = NA,
          `Neg LogLikelihood` = if (!is.na(erl_loglik)) round(erl_loglik, 2) else NA_real_,
          `p-value` = round(res$p_value, 4),
          Result = ifelse(res$p_value >= res$alpha, "PASS", "FAIL"),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      } else {
        erlexp_loglik <- calc_erlangexp_loglik(
          res$k_star,
          res$erlang_lambda_star,
          res$exp_lambda_star,
          empiricaldata
        )
        results$single <- data.frame(
          Model = "Erlang-Exp",
          K = res$k_star,
          `Erlang λ` = round(res$erlang_lambda_star, 4),
          `Exp λ` = round(res$exp_lambda_star, 4),
          `Neg LogLikelihood` = if (!is.na(erlexp_loglik)) round(erlexp_loglik, 2) else NA_real_,
          `p-value` = round(res$p_value, 4),
          Result = ifelse(res$p_value >= res$alpha, "PASS", "FAIL"),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }

      # Smallest K for single model
      smallest_k_single <- tryCatch(smallest_k_results(), error = function(e) NULL)

      if (!is.null(smallest_k_single)) {
        if (res$fit_type == "Erlang") {
          sk_erl_loglik <- calc_erlang_loglik(
            smallest_k_single$smallest_k,
            smallest_k_single$smallest_lambda,
            empiricaldata
          )
          results$smallest_k_single <- data.frame(
            Model = "Erlang (Smallest K)",
            K = smallest_k_single$smallest_k,
            `Erlang λ` = round(smallest_k_single$smallest_lambda, 4),
            `Exp λ` = NA,
            `Neg LogLikelihood` = if (!is.na(sk_erl_loglik)) round(sk_erl_loglik, 2) else NA_real_,
            `p-value` = round(smallest_k_single$smallest_p_value, 4),
            Result = ifelse(smallest_k_single$smallest_q_value == 1, "PASS", "FAIL"),
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        } else {
          sk_erlexp_loglik <- calc_erlangexp_loglik(
            smallest_k_single$smallest_k,
            smallest_k_single$smallest_erlang_lambda,
            smallest_k_single$smallest_exp_lambda,
            empiricaldata
          )
          results$smallest_k_single <- data.frame(
            Model = "Erlang-Exp (Smallest K)",
            K = smallest_k_single$smallest_k,
            `Erlang λ` = round(smallest_k_single$smallest_erlang_lambda, 4),
            `Exp λ` = round(smallest_k_single$smallest_exp_lambda, 4),
            `Neg LogLikelihood` = if (!is.na(sk_erlexp_loglik)) round(sk_erlexp_loglik, 2) else NA_real_,
            `p-value` = round(smallest_k_single$smallest_p_value, 4),
            Result = ifelse(smallest_k_single$smallest_q_value == 1, "PASS", "FAIL"),
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        }
      }
    }

    # Combine all results
    if (length(results) > 0) {
      # DEBUG - add these lines
      print("=== DEBUG results list before rbind ===")
      print(paste("Length of results:", length(results)))
      print(paste("Names of results:", paste(names(results), collapse=", ")))
      for (nm in names(results)) {
        print(paste("  ", nm, "- class:", class(results[[nm]])))
        print(paste("  ", nm, "- is.data.frame:", is.data.frame(results[[nm]])))
        if (is.data.frame(results[[nm]])) {
          print(paste("  ", nm, "- ncol:", ncol(results[[nm]])))
          print(paste("  ", nm, "- colnames:", paste(colnames(results[[nm]]), collapse=", ")))
        }
      }
      print("=== END DEBUG ===")
      do.call(rbind, results)
    } else {
      NULL
    }
  })


  # =========================================================
  # GOF SUMMARY TABLE RENDER
  # =========================================================

  output$gof_summary_table <- DT::renderDT({

    df <- gof_summary_data()
    req(df)
    # Get alpha value for caption
    gof_res <- tryCatch(gof_results(), error = function(e) NULL)
    alpha_val <- if (!is.null(gof_res)) gof_res$alpha else 0.05

    DT::datatable(
      df,
      rownames = FALSE,
      options = list(
        dom = 't',
        ordering = FALSE,
        pageLength = 10,
        columnDefs = list(
          list(className = 'dt-center', targets = '_all')
        )
      ),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left; font-size: 0.95em; color: #555;',
        sprintf('PASS: p-value ≥ %.3f (consistent with data) | FAIL: p-value < %.3f (not consistent with data)', alpha_val, alpha_val)
      )
    ) %>%
      DT::formatStyle(
        'Result',
        backgroundColor = DT::styleEqual(
          c('PASS', 'FAIL'),
          c('#d4edda', '#f8d7da')
        ),
        color = DT::styleEqual(
          c('PASS', 'FAIL'),
          c('#155724', '#721c24')
        ),
        fontWeight = 'bold'
      ) %>%
      DT::formatStyle(
        'Model',
        fontWeight = 'bold'
      )
  })


  # =========================================================
  # DYNAMIC GOF MAIN OUTPUT
  # =========================================================

  output$gof_main_output <- renderUI({

    # Check if GOF results exist
    gof_exists <- !is.null(tryCatch(gof_results(), error = function(e) NULL))

    if (!gof_exists) {
      return(
        div(
          class = "alert alert-info",
          "Click 'Compute GOF' to see results."
        )
      )
    }

    res <- gof_results()

    # =====================================================
    # DEFAULT (BOTH MODELS) OUTPUT
    # =====================================================
    if (res$fit_type == "Default") {

      return(
        tagList(
          # ---------------------------------------------------
          # SUMMARY TABLE AT TOP
          # ---------------------------------------------------
          card(
            card_header(
              class = "bg-light",
              tags$h4(icon("table"), "Results Summary", style = "margin: 0;")
            ),
            card_body(
              DT::DTOutput("gof_summary_table")
            )
          ),

          br(),

          # ---------------------------------------------------
          # PLOTS ROW - Combined CDF
          # ---------------------------------------------------
          fluidRow(
            column(
              width = 12,
              plotOutput(
                "gof_cdf_plot_combined",
                height = "400px"
              )
            )
          ),

          br(),

          # ---------------------------------------------------
          # BOOTSTRAP HISTOGRAMS ROW
          # ---------------------------------------------------
          fluidRow(
            column(
              width = 6,
              h4("Erlang Bootstrap Distribution", style = "color: red;"),
              plotOutput(
                "gof_bootstrap_plot_erlang",
                height = "350px"
              )
            ),
            column(
              width = 6,
              h4("Erlang-Exp Bootstrap Distribution", style = "color: blue;"),
              plotOutput(
                "gof_bootstrap_plot_erlang_exp",
                height = "350px"
              )
            )
          ),

          br(),

          hr(),

          # ---------------------------------------------------
          # TEXT OUTPUT
          # ---------------------------------------------------
          h4("Detailed Results"),

          fluidRow(
            column(
              width = 6,
              h5("Erlang GOF Results", style = "color: red;"),
              verbatimTextOutput("gof_output_erlang")
            ),
            column(
              width = 6,
              h5("Erlang-Exp GOF Results", style = "color: blue;"),
              verbatimTextOutput("gof_output_erlang_exp")
            )
          ),

          # ---------------------------------------------------
          # SMALLEST K OUTPUT FOR DEFAULT (appears after run)
          # ---------------------------------------------------
          uiOutput("smallest_k_output_section_default")
        )
      )
    }

    # =====================================================
    # SINGLE MODEL OUTPUT (Erlang or Erlang-Exp)
    # =====================================================
    return(
      tagList(
        # ---------------------------------------------------
        # SUMMARY TABLE AT TOP
        # ---------------------------------------------------
        card(
          card_header(
            class = "bg-light",
            tags$h4(icon("table"), "Results Summary", style = "margin: 0;")
          ),
          card_body(
            DT::DTOutput("gof_summary_table")
          )
        ),

        br(),

        # ---------------------------------------------------
        # PLOTS ROW
        # ---------------------------------------------------
        fluidRow(
          column(
            width = 6,
            plotOutput(
              "gof_cdf_plot",
              height = "400px"
            )
          ),
          column(
            width = 6,
            plotOutput(
              "gof_bootstrap_plot",
              height = "400px"
            )
          )
        ),

        br(),

        hr(),

        # ---------------------------------------------------
        # TEXT OUTPUT
        # ---------------------------------------------------
        h4("Detailed Results"),

        verbatimTextOutput("gof_output"),

        # ---------------------------------------------------
        # SMALLEST K OUTPUT (appears after Smallest K is run)
        # ---------------------------------------------------
        uiOutput("smallest_k_output_section")
      )
    )
  })


  # =========================================================
  # SMALLEST K SECTION (appears after GOF is run)
  # =========================================================

  output$smallest_k_section <- renderUI({

    # Check if GOF results exist
    gof_exists <- !is.null(tryCatch(gof_results(), error = function(e) NULL))

    if (!gof_exists) {
      return(NULL)
    }

    res <- gof_results()

    # =====================================================
    # DEFAULT MODE - Show checkboxes for model selection
    # =====================================================
    if (res$fit_type == "Default") {

      stat_name <- switch(
        toupper(res$pvaloption),
        "KS" = "Kolmogorov-Smirnov",
        "AD" = "Anderson-Darling",
        "CVM" = "Cramér-von Mises",
        res$pvaloption
      )

      return(
        div(
          class = "smallest-k-box",

          strong("Find Smallest K"),
          br(),
          br(),

          p(
            style = "font-size: 12px; color: #555;",
            sprintf(
              "Select model(s) to find Smallest K based on %s test statistic (α = %.3f)",
              stat_name,
              res$alpha
            )
          ),

          # Checkboxes for model selection
          div(
            style = "margin-bottom: 10px;",
            checkboxInput(
              "smallest_k_erlang_check",
              tags$span(style = "color: red;", "Erlang"),
              value = FALSE
            ),
            checkboxInput(
              "smallest_k_erlang_exp_check",
              tags$span(style = "color: blue;", "Erlang-Exp"),
              value = FALSE
            )
          ),

          actionButton(
            "run_smallest_k_default",
            "Compute Smallest K",
            class = "btn-success btn-sm"
          )
        )
      )
    }

    # =====================================================
    # SINGLE MODEL MODE
    # =====================================================
    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    div(
      class = "smallest-k-box",

      strong("Find Smallest K"),
      br(),
      br(),

      p(
        style = "font-size: 12px; color: #555;",
        sprintf(
          "Compute Smallest K based on %s test statistic (α = %.3f)",
          stat_name,
          res$alpha
        )
      ),

      actionButton(
        "run_smallest_k",
        "Compute Smallest K",
        class = "btn-success btn-sm"
      )
    )
  })


  # =========================================================
  # SMALLEST K COMPUTATION (Single Model)
  # =========================================================

  smallest_k_results <- eventReactive(input$run_smallest_k, {

    req(gof_results())
    req(data())

    empiricaldata <- data()[[1]]
    gof_res <- gof_results()

    # Get the pvaloption and alpha from GOF results
    pvaloption <- gof_res$pvaloption
    alpha <- gof_res$alpha

    # =====================================================
    # ERLANG SMALLEST K
    # =====================================================
    if (input$fit_type == "Erlang") {

      fit_args <- list(
        mode = "Erlang",
        empiricaldata = empiricaldata,
        SmallestK = TRUE,
        pvaloption = pvaloption,
        Alpha = alpha
      )

      # Include initial K if it was specified
      if (!is.na(input$initial_k)) {
        fit_args$K <- input$initial_k
      }

      results <- do.call(
        GenErlang_Fit,
        fit_args
      )

      return(list(
        fit_type = "Erlang",
        pvaloption = pvaloption,
        alpha = alpha,
        smallest_k = results$Smallest$K_star,
        smallest_lambda = results$Smallest$Lambda_star,
        smallest_p_value = results$Smallest$P_star,
        smallest_q_value = results$Smallest$Q_Value,
        smallest_metric = results$Smallest$metric_star,
        smallest_sample_stats = results$Smallest$samplestats_star,
        best_k = results$Best$K_star,
        best_lambda = results$Best$Lambda_star,
        empiricaldata = empiricaldata
      ))
    }

    # =====================================================
    # ERLANG-EXP SMALLEST K
    # =====================================================
    if (input$fit_type == "Erlang-Exp") {

      fit_args <- list(
        mode = "ErlangExp",
        empiricaldata = empiricaldata,
        K = input$initial_k_exp,
        SmallestK = TRUE,
        pvaloption = pvaloption,
        Alpha = alpha
      )

      # Set FixedK or Window based on original selection
      if (input$search_type == "Fixed K") {
        fit_args$FixedK <- TRUE
      } else if (input$search_type == "Search over a Window") {
        fit_args$FixedK <- FALSE
        fit_args$KWindowSize <- input$window_size
      }

      results <- do.call(
        GenErlang_Fit,
        fit_args
      )

      return(list(
        fit_type = "Erlang-Exp",
        pvaloption = pvaloption,
        alpha = alpha,
        smallest_k = results$Smallest$K_star,
        smallest_erlang_lambda = results$Smallest$ErlangLambda_star,
        smallest_exp_lambda = results$Smallest$ExpLambda_star,
        smallest_p_value = results$Smallest$P_star,
        smallest_q_value = results$Smallest$Q_Value,
        smallest_metric = results$Smallest$metric_star,
        smallest_sample_stats = results$Smallest$samplestats_star,
        best_k = results$Best$K_star,
        best_erlang_lambda = results$Best$ErlangLambda_star,
        best_exp_lambda = results$Best$ExpLambda_star,
        empiricaldata = empiricaldata
      ))
    }

    return(NULL)
  })


  # =========================================================
  # SMALLEST K COMPUTATION (Default Mode - Both Models)
  # =========================================================

  smallest_k_results_default <- eventReactive(input$run_smallest_k_default, {

    req(gof_results())
    req(data())

    # Check if at least one checkbox is selected
    erlang_selected <- isTRUE(input$smallest_k_erlang_check)
    erlang_exp_selected <- isTRUE(input$smallest_k_erlang_exp_check)

    if (!erlang_selected && !erlang_exp_selected) {
      return(list(error = "Please select at least one model."))
    }

    empiricaldata <- data()[[1]]
    gof_res <- gof_results()

    # Get the pvaloption and alpha from GOF results
    pvaloption <- gof_res$pvaloption
    alpha <- gof_res$alpha

    results <- list(
      pvaloption = pvaloption,
      alpha = alpha,
      erlang_selected = erlang_selected,
      erlang_exp_selected = erlang_exp_selected,
      empiricaldata = empiricaldata
    )

    # =====================================================
    # ERLANG SMALLEST K
    # =====================================================
    if (erlang_selected) {

      fit_args <- list(
        mode = "Erlang",
        empiricaldata = empiricaldata,
        SmallestK = TRUE,
        pvaloption = pvaloption,
        Alpha = alpha
      )

      erlang_results <- do.call(
        GenErlang_Fit,
        fit_args
      )

      results$erlang <- list(
        smallest_k = erlang_results$Smallest$K_star,
        smallest_lambda = erlang_results$Smallest$Lambda_star,
        smallest_p_value = erlang_results$Smallest$P_star,
        smallest_q_value = erlang_results$Smallest$Q_Value,
        smallest_metric = erlang_results$Smallest$metric_star,
        smallest_sample_stats = erlang_results$Smallest$samplestats_star,
        best_k = erlang_results$Best$K_star,
        best_lambda = erlang_results$Best$Lambda_star
      )
    }

    # =====================================================
    # ERLANG-EXP SMALLEST K
    # =====================================================
    if (erlang_exp_selected) {

      # Get K from the original fit results
      k_erlang_exp <- gof_res$erlang_exp$k_star

      fit_args <- list(
        mode = "ErlangExp",
        empiricaldata = empiricaldata,
        K = k_erlang_exp,
        SmallestK = TRUE,
        pvaloption = pvaloption,
        Alpha = alpha,
        FixedK = TRUE  # Use Fixed K for Default mode
      )

      erlang_exp_results <- do.call(
        GenErlang_Fit,
        fit_args
      )

      results$erlang_exp <- list(
        smallest_k = erlang_exp_results$Smallest$K_star,
        smallest_erlang_lambda = erlang_exp_results$Smallest$ErlangLambda_star,
        smallest_exp_lambda = erlang_exp_results$Smallest$ExpLambda_star,
        smallest_p_value = erlang_exp_results$Smallest$P_star,
        smallest_q_value = erlang_exp_results$Smallest$Q_Value,
        smallest_metric = erlang_exp_results$Smallest$metric_star,
        smallest_sample_stats = erlang_exp_results$Smallest$samplestats_star,
        best_k = erlang_exp_results$Best$K_star,
        best_erlang_lambda = erlang_exp_results$Best$ErlangLambda_star,
        best_exp_lambda = erlang_exp_results$Best$ExpLambda_star
      )
    }

    return(results)
  })


  # =========================================================
  # SMALLEST K OUTPUT SECTION (Single Model)
  # =========================================================

  output$smallest_k_output_section <- renderUI({

    # Check if smallest K results exist
    smallest_k_exists <- !is.null(tryCatch(smallest_k_results(), error = function(e) NULL))

    if (!smallest_k_exists) {
      return(NULL)
    }

    tagList(
      hr(),

      h4("Smallest K Results"),

      br(),

      # ---------------------------------------------------
      # ROW 1: PDF Plot (full width)
      # ---------------------------------------------------
      fluidRow(
        column(
          width = 12,
          plotOutput(
            "smallest_k_pdf_plot",
            height = "450px"
          )
        )
      ),

      br(),

      # ---------------------------------------------------
      # ROW 2: CDF Plot and Bootstrap Histogram
      # ---------------------------------------------------
      fluidRow(
        column(
          width = 6,
          plotOutput(
            "smallest_k_cdf_plot",
            height = "400px"
          )
        ),
        column(
          width = 6,
          plotOutput(
            "smallest_k_bootstrap_plot",
            height = "400px"
          )
        )
      ),

      br(),

      hr(),

      h4("Detailed Smallest K Results"),

      verbatimTextOutput("smallest_k_output")
    )
  })


  # =========================================================
  # SMALLEST K OUTPUT SECTION (Default Mode)
  # =========================================================

  output$smallest_k_output_section_default <- renderUI({

    # Check if smallest K results exist for default mode
    smallest_k_exists <- !is.null(tryCatch(smallest_k_results_default(), error = function(e) NULL))

    if (!smallest_k_exists) {
      return(NULL)
    }

    res <- smallest_k_results_default()

    # Check for error
    if (!is.null(res$error)) {
      return(
        div(
          class = "alert alert-warning",
          icon("exclamation-triangle"),
          res$error
        )
      )
    }

    erlang_selected <- res$erlang_selected
    erlang_exp_selected <- res$erlang_exp_selected

    tagList(
      hr(),

      h4("Smallest K Results"),

      br(),

      # ---------------------------------------------------
      # PDF PLOTS
      # ---------------------------------------------------
      if (erlang_selected && erlang_exp_selected) {
        # Both models - side by side
        fluidRow(
          column(
            width = 6,
            h5("Erlang: Best K vs Smallest K", style = "color: red;"),
            plotOutput(
              "smallest_k_pdf_plot_erlang_default",
              height = "400px"
            )
          ),
          column(
            width = 6,
            h5("Erlang-Exp: Best K vs Smallest K", style = "color: blue;"),
            plotOutput(
              "smallest_k_pdf_plot_erlang_exp_default",
              height = "400px"
            )
          )
        )
      } else if (erlang_selected) {
        # Only Erlang
        fluidRow(
          column(
            width = 12,
            plotOutput(
              "smallest_k_pdf_plot_erlang_default",
              height = "450px"
            )
          )
        )
      } else {
        # Only Erlang-Exp
        fluidRow(
          column(
            width = 12,
            plotOutput(
              "smallest_k_pdf_plot_erlang_exp_default",
              height = "450px"
            )
          )
        )
      },

      br(),

      # ---------------------------------------------------
      # CDF AND BOOTSTRAP PLOTS
      # ---------------------------------------------------
      if (erlang_selected && erlang_exp_selected) {
        # Both models
        tagList(
          # Erlang CDF and Bootstrap
          fluidRow(
            column(
              width = 6,
              h5("Erlang CDF Comparison", style = "color: red;"),
              plotOutput(
                "smallest_k_cdf_plot_erlang_default",
                height = "350px"
              )
            ),
            column(
              width = 6,
              h5("Erlang Bootstrap Distribution", style = "color: red;"),
              plotOutput(
                "smallest_k_bootstrap_plot_erlang_default",
                height = "350px"
              )
            )
          ),
          br(),
          # Erlang-Exp CDF and Bootstrap
          fluidRow(
            column(
              width = 6,
              h5("Erlang-Exp CDF Comparison", style = "color: blue;"),
              plotOutput(
                "smallest_k_cdf_plot_erlang_exp_default",
                height = "350px"
              )
            ),
            column(
              width = 6,
              h5("Erlang-Exp Bootstrap Distribution", style = "color: blue;"),
              plotOutput(
                "smallest_k_bootstrap_plot_erlang_exp_default",
                height = "350px"
              )
            )
          )
        )
      } else if (erlang_selected) {
        # Only Erlang
        fluidRow(
          column(
            width = 6,
            plotOutput(
              "smallest_k_cdf_plot_erlang_default",
              height = "400px"
            )
          ),
          column(
            width = 6,
            plotOutput(
              "smallest_k_bootstrap_plot_erlang_default",
              height = "400px"
            )
          )
        )
      } else {
        # Only Erlang-Exp
        fluidRow(
          column(
            width = 6,
            plotOutput(
              "smallest_k_cdf_plot_erlang_exp_default",
              height = "400px"
            )
          ),
          column(
            width = 6,
            plotOutput(
              "smallest_k_bootstrap_plot_erlang_exp_default",
              height = "400px"
            )
          )
        )
      },

      br(),

      hr(),

      # ---------------------------------------------------
      # TEXT OUTPUT
      # ---------------------------------------------------
      h4("Detailed Smallest K Results"),

      if (erlang_selected && erlang_exp_selected) {
        fluidRow(
          column(
            width = 6,
            h5("Erlang", style = "color: red;"),
            verbatimTextOutput("smallest_k_output_erlang_default")
          ),
          column(
            width = 6,
            h5("Erlang-Exp", style = "color: blue;"),
            verbatimTextOutput("smallest_k_output_erlang_exp_default")
          )
        )
      } else if (erlang_selected) {
        verbatimTextOutput("smallest_k_output_erlang_default")
      } else {
        verbatimTextOutput("smallest_k_output_erlang_exp_default")
      }
    )
  })


  # =========================================================
  # SMALLEST K PDF PLOT - ERLANG (Default Mode)
  # =========================================================

  output$smallest_k_pdf_plot_erlang_default <- renderPlot({

    req(smallest_k_results_default())

    res <- smallest_k_results_default()

    if (is.null(res$erlang)) {
      return(NULL)
    }

    empiricaldata <- res$empiricaldata
    erlang <- res$erlang

    # Calculate bin width
    bin_width <- 2 * IQR(empiricaldata) / (length(empiricaldata)^(1/3))

    # For cases where all data points are identical or nearly identical
    if (bin_width < .Machine$double.eps || !is.finite(bin_width)) {
      bin_width <- 1  # Default fallback bin width
    }

    # Create x grid for density curves
    x_grid <- seq(0, 1.2 * max(empiricaldata), length.out = 1000)

    # Best K density
    best_density <- dgamma(
      x_grid,
      shape = erlang$best_k,
      scale = 1/erlang$best_lambda
    )

    # Smallest K density
    smallest_density <- dgamma(
      x_grid,
      shape = erlang$smallest_k,
      scale = 1/erlang$smallest_lambda
    )

    # Create labels
    best_label <- sprintf("Best K (MLE): K=%d, λ=%.3f", erlang$best_k, erlang$best_lambda)
    smallest_label <- sprintf("Smallest K: K=%d, λ=%.3f", erlang$smallest_k, erlang$smallest_lambda)

    # Create data frames
    df_best <- data.frame(x = x_grid, density = best_density)
    df_smallest <- data.frame(x = x_grid, density = smallest_density)

    # Base plot
    p <- ggplot(
      data.frame(Value = empiricaldata),
      aes(x = Value)
    ) +
      geom_histogram(
        aes(y = after_stat(density)),
        binwidth = bin_width,
        color = "black",
        fill = "#90C0AE",
        alpha = 0.7
      ) +
      geom_line(
        data = df_best,
        aes(x = x, y = density, color = "Best K (MLE)"),
        linewidth = 1.2
      ) +
      geom_line(
        data = df_smallest,
        aes(x = x, y = density, color = "Smallest K"),
        linewidth = 1.2,
        linetype = "dashed"
      ) +
      scale_color_manual(
        values = c(
          "Best K (MLE)" = "red",
          "Smallest K" = "darkgreen"
        ),
        labels = c(
          "Best K (MLE)" = best_label,
          "Smallest K" = smallest_label
        ),
        name = NULL
      ) +
      labs(
        title = "Erlang: Data Histogram with Fitted Distributions",
        subtitle = "Comparing Best K (MLE) vs Smallest K",
        x = "Observed Values",
        y = "Density"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "bottom"
      )

    p
  })


  # =========================================================
  # SMALLEST K PDF PLOT - ERLANG-EXP (Default Mode)
  # =========================================================

  output$smallest_k_pdf_plot_erlang_exp_default <- renderPlot({

    req(smallest_k_results_default())

    res <- smallest_k_results_default()

    if (is.null(res$erlang_exp)) {
      return(NULL)
    }

    empiricaldata <- res$empiricaldata
    erlang_exp <- res$erlang_exp

    # Calculate bin width
    bin_width <- 2 * IQR(empiricaldata) / (length(empiricaldata)^(1/3))

    # For cases where all data points are identical or nearly identical
    if (bin_width < .Machine$double.eps || !is.finite(bin_width)) {
      bin_width <- 1  # Default fallback bin width
    }

    # Create x grid for density curves
    x_grid <- seq(0, 1.2 * max(empiricaldata), length.out = 1000)

    # Best K density
    best_density <- GenErlangFit:::ErlangExp_Func(
      x_grid,
      ErK = erlang_exp$best_k,
      Erlam = erlang_exp$best_erlang_lambda,
      Explam = erlang_exp$best_exp_lambda
    )$Probability

    # Smallest K density
    smallest_density <- GenErlangFit:::ErlangExp_Func(
      x_grid,
      ErK = erlang_exp$smallest_k,
      Erlam = erlang_exp$smallest_erlang_lambda,
      Explam = erlang_exp$smallest_exp_lambda
    )$Probability

    # Create labels
    best_label <- sprintf(
      "Best K (MLE): K=%d, λEr=%.3f, λExp=%.3f",
      erlang_exp$best_k, erlang_exp$best_erlang_lambda, erlang_exp$best_exp_lambda
    )
    smallest_label <- sprintf(
      "Smallest K: K=%d, λEr=%.3f, λExp=%.3f",
      erlang_exp$smallest_k, erlang_exp$smallest_erlang_lambda, erlang_exp$smallest_exp_lambda
    )

    # Create data frames
    df_best <- data.frame(x = x_grid, density = best_density)
    df_smallest <- data.frame(x = x_grid, density = smallest_density)

    # Base plot
    p <- ggplot(
      data.frame(Value = empiricaldata),
      aes(x = Value)
    ) +
      geom_histogram(
        aes(y = after_stat(density)),
        binwidth = bin_width,
        color = "black",
        fill = "#90C0AE",
        alpha = 0.7
      ) +
      geom_line(
        data = df_best,
        aes(x = x, y = density, color = "Best K (MLE)"),
        linewidth = 1.2
      ) +
      geom_line(
        data = df_smallest,
        aes(x = x, y = density, color = "Smallest K"),
        linewidth = 1.2,
        linetype = "dashed"
      ) +
      scale_color_manual(
        values = c(
          "Best K (MLE)" = "blue",
          "Smallest K" = "purple"
        ),
        labels = c(
          "Best K (MLE)" = best_label,
          "Smallest K" = smallest_label
        ),
        name = NULL
      ) +
      labs(
        title = "Erlang-Exp: Data Histogram with Fitted Distributions",
        subtitle = "Comparing Best K (MLE) vs Smallest K",
        x = "Observed Values",
        y = "Density"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "bottom"
      )

    p
  })


  # =========================================================
  # SMALLEST K CDF PLOT - ERLANG (Default Mode)
  # =========================================================

  output$smallest_k_cdf_plot_erlang_default <- renderPlot({

    req(smallest_k_results_default())

    res <- smallest_k_results_default()

    if (is.null(res$erlang)) {
      return(NULL)
    }

    empiricaldata <- res$empiricaldata
    erlang <- res$erlang

    # Compute empirical CDF
    ecdf_data <- ecdf(empiricaldata)
    x_vals <- sort(empiricaldata)
    ecdf_vals <- ecdf_data(x_vals)

    # Best K CDF
    best_cdf <- pgamma(x_vals, shape = erlang$best_k, scale = 1/erlang$best_lambda)

    # Smallest K CDF
    smallest_cdf <- pgamma(x_vals, shape = erlang$smallest_k, scale = 1/erlang$smallest_lambda)

    # Create labels
    best_label <- sprintf("Best K: K=%d, λ=%.3f", erlang$best_k, erlang$best_lambda)
    smallest_label <- sprintf("Smallest K: K=%d, λ=%.3f", erlang$smallest_k, erlang$smallest_lambda)

    # Build data frame
    df_cdf <- data.frame(
      x = rep(x_vals, 3),
      cdf = c(ecdf_vals, best_cdf, smallest_cdf),
      Type = factor(
        c(
          rep("Empirical CDF", length(x_vals)),
          rep(best_label, length(x_vals)),
          rep(smallest_label, length(x_vals))
        ),
        levels = c("Empirical CDF", best_label, smallest_label)
      )
    )

    # Create plot
    p <- ggplot(df_cdf, aes(x = x, y = cdf, color = Type, linetype = Type)) +
      geom_step(
        data = subset(df_cdf, Type == "Empirical CDF"),
        linewidth = 1.2
      ) +
      geom_line(
        data = subset(df_cdf, Type == best_label),
        linewidth = 1.2
      ) +
      geom_line(
        data = subset(df_cdf, Type == smallest_label),
        linewidth = 1.2
      ) +
      scale_color_manual(
        values = c(
          "Empirical CDF" = "black",
          setNames("red", best_label),
          setNames("darkgreen", smallest_label)
        )
      ) +
      scale_linetype_manual(
        values = c(
          "Empirical CDF" = "solid",
          setNames("solid", best_label),
          setNames("dashed", smallest_label)
        )
      ) +
      labs(
        title = "Erlang CDF Comparison",
        x = "x",
        y = "CDF",
        color = NULL,
        linetype = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = -10),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      guides(
        color = guide_legend(ncol = 1),
        linetype = guide_legend(ncol = 1)
      )

    p
  })


  # =========================================================
  # SMALLEST K CDF PLOT - ERLANG-EXP (Default Mode)
  # =========================================================

  output$smallest_k_cdf_plot_erlang_exp_default <- renderPlot({

    req(smallest_k_results_default())

    res <- smallest_k_results_default()

    if (is.null(res$erlang_exp)) {
      return(NULL)
    }

    empiricaldata <- res$empiricaldata
    erlang_exp <- res$erlang_exp

    # Compute empirical CDF
    ecdf_data <- ecdf(empiricaldata)
    x_vals <- sort(empiricaldata)
    ecdf_vals <- ecdf_data(x_vals)

    # Best K CDF
    best_params <- c(erlang_exp$best_k, erlang_exp$best_erlang_lambda, erlang_exp$best_exp_lambda)
    best_cdf <- GenErlangFit:::ErlangExpCDF_Func(best_params, x_vals, interval = 0.01)

    # Smallest K CDF
    smallest_params <- c(erlang_exp$smallest_k, erlang_exp$smallest_erlang_lambda, erlang_exp$smallest_exp_lambda)
    smallest_cdf <- GenErlangFit:::ErlangExpCDF_Func(smallest_params, x_vals, interval = 0.01)

    # Create labels
    best_label <- sprintf(
      "Best K: K=%d, λEr=%.2f, λExp=%.2f",
      erlang_exp$best_k, erlang_exp$best_erlang_lambda, erlang_exp$best_exp_lambda
    )
    smallest_label <- sprintf(
      "Smallest K: K=%d, λEr=%.2f, λExp=%.2f",
      erlang_exp$smallest_k, erlang_exp$smallest_erlang_lambda, erlang_exp$smallest_exp_lambda
    )

    # Build data frame
    df_cdf <- data.frame(
      x = rep(x_vals, 3),
      cdf = c(ecdf_vals, best_cdf, smallest_cdf),
      Type = factor(
        c(
          rep("Empirical CDF", length(x_vals)),
          rep(best_label, length(x_vals)),
          rep(smallest_label, length(x_vals))
        ),
        levels = c("Empirical CDF", best_label, smallest_label)
      )
    )

    # Create plot
    p <- ggplot(df_cdf, aes(x = x, y = cdf, color = Type, linetype = Type)) +
      geom_step(
        data = subset(df_cdf, Type == "Empirical CDF"),
        linewidth = 1.2
      ) +
      geom_line(
        data = subset(df_cdf, Type == best_label),
        linewidth = 1.2
      ) +
      geom_line(
        data = subset(df_cdf, Type == smallest_label),
        linewidth = 1.2
      ) +
      scale_color_manual(
        values = c(
          "Empirical CDF" = "black",
          setNames("blue", best_label),
          setNames("purple", smallest_label)
        )
      ) +
      scale_linetype_manual(
        values = c(
          "Empirical CDF" = "solid",
          setNames("solid", best_label),
          setNames("dashed", smallest_label)
        )
      ) +
      labs(
        title = "Erlang-Exp CDF Comparison",
        x = "x",
        y = "CDF",
        color = NULL,
        linetype = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = -10),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      guides(
        color = guide_legend(ncol = 1),
        linetype = guide_legend(ncol = 1)
      )

    p
  })


  # =========================================================
  # SMALLEST K BOOTSTRAP PLOT - ERLANG (Default Mode)
  # =========================================================

  output$smallest_k_bootstrap_plot_erlang_default <- renderPlot({

    req(smallest_k_results_default())

    res <- smallest_k_results_default()

    if (is.null(res$erlang)) {
      return(NULL)
    }

    erlang <- res$erlang

    # Get test statistic name for labels
    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    # Create data frame for histogram
    df_bootstrap <- data.frame(
      Statistic = erlang$smallest_sample_stats
    )

    # Calculate appropriate bin width
    bin_width <- diff(range(erlang$smallest_sample_stats)) / 30

    # Create plot
    p <- ggplot(df_bootstrap, aes(x = Statistic)) +
      geom_histogram(
        binwidth = bin_width,
        fill = "#98D8AA",
        color = "black",
        alpha = 0.7
      ) +
      geom_vline(
        xintercept = erlang$smallest_metric,
        linetype = "dashed",
        color = "black",
        linewidth = 1.2
      ) +
      annotate(
        "text",
        x = erlang$smallest_metric,
        y = Inf,
        label = sprintf("Observed = %.4f", erlang$smallest_metric),
        hjust = -0.1,
        vjust = 2,
        size = 4,
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = sprintf("p-value = %.4f", erlang$smallest_p_value),
        hjust = 1.1,
        vjust = 2,
        size = 4,
        fontface = "bold",
        color = if (erlang$smallest_q_value == 1) "darkgreen" else "red"
      ) +
      labs(
        title = sprintf("Erlang Bootstrap (%s)", stat_name),
        subtitle = sprintf("Smallest K = %d", erlang$smallest_k),
        x = sprintf("%s Statistic", stat_name),
        y = "Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40")
      )

    p
  })


  # =========================================================
  # SMALLEST K BOOTSTRAP PLOT - ERLANG-EXP (Default Mode)
  # =========================================================

  output$smallest_k_bootstrap_plot_erlang_exp_default <- renderPlot({

    req(smallest_k_results_default())

    res <- smallest_k_results_default()

    if (is.null(res$erlang_exp)) {
      return(NULL)
    }

    erlang_exp <- res$erlang_exp

    # Get test statistic name for labels
    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    # Create data frame for histogram
    df_bootstrap <- data.frame(
      Statistic = erlang_exp$smallest_sample_stats
    )

    # Calculate appropriate bin width
    bin_width <- diff(range(erlang_exp$smallest_sample_stats)) / 30

    # Create plot
    p <- ggplot(df_bootstrap, aes(x = Statistic)) +
      geom_histogram(
        binwidth = bin_width,
        fill = "#B4A7D6",
        color = "black",
        alpha = 0.7
      ) +
      geom_vline(
        xintercept = erlang_exp$smallest_metric,
        linetype = "dashed",
        color = "black",
        linewidth = 1.2
      ) +
      annotate(
        "text",
        x = erlang_exp$smallest_metric,
        y = Inf,
        label = sprintf("Observed = %.4f", erlang_exp$smallest_metric),
        hjust = -0.1,
        vjust = 2,
        size = 4,
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = sprintf("p-value = %.4f", erlang_exp$smallest_p_value),
        hjust = 1.1,
        vjust = 2,
        size = 4,
        fontface = "bold",
        color = if (erlang_exp$smallest_q_value == 1) "darkgreen" else "red"
      ) +
      labs(
        title = sprintf("Erlang-Exp Bootstrap (%s)", stat_name),
        subtitle = sprintf("Smallest K = %d", erlang_exp$smallest_k),
        x = sprintf("%s Statistic", stat_name),
        y = "Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40")
      )

    p
  })


  # =========================================================
  # SMALLEST K TEXT OUTPUT - ERLANG (Default Mode)
  # =========================================================

  output$smallest_k_output_erlang_default <- renderPrint({

    req(smallest_k_results_default())

    res <- smallest_k_results_default()

    if (is.null(res$erlang)) {
      return(NULL)
    }

    erlang <- res$erlang

    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    cat("==============================================\n")
    cat("       SMALLEST K RESULTS (ERLANG)            \n")
    cat("==============================================\n\n")

    cat("--- Search Settings ---\n")
    cat(sprintf("  Test Statistic : %s (%s)\n", res$pvaloption, stat_name))
    cat(sprintf("  Alpha          : %.4f\n", res$alpha))
    cat("\n")

    cat("--- Best K (MLE) ---\n")
    cat(sprintf("  K*      : %d\n", erlang$best_k))
    cat(sprintf("  Lambda* : %.6f\n", erlang$best_lambda))
    cat("\n")

    cat("--- Smallest K (passing GOF) ---\n")
    cat(sprintf("  K*      : %d\n", erlang$smallest_k))
    cat(sprintf("  Lambda* : %.6f\n", erlang$smallest_lambda))
    cat("\n")

    cat("--- GOF Results for Smallest K ---\n")
    cat(sprintf("  Observed Statistic : %.6f\n", erlang$smallest_metric))
    cat(sprintf("  P-value            : %.6f\n", erlang$smallest_p_value))
    cat(sprintf("  Q-value            : %d\n", erlang$smallest_q_value))
    cat("\n")

    cat("--- Interpretation ---\n")
    if (erlang$smallest_q_value == 1) {
      cat(sprintf("  Smallest K = %d passes the %s test at alpha = %.4f\n",
                  erlang$smallest_k, stat_name, res$alpha))
    } else {
      cat(sprintf("  No K found that passes the %s test at alpha = %.4f\n",
                  stat_name, res$alpha))
    }
    cat("\n")

    cat("==============================================\n")
  })


  # =========================================================
  # SMALLEST K TEXT OUTPUT - ERLANG-EXP (Default Mode)
  # =========================================================

  output$smallest_k_output_erlang_exp_default <- renderPrint({

    req(smallest_k_results_default())

    res <- smallest_k_results_default()

    if (is.null(res$erlang_exp)) {
      return(NULL)
    }

    erlang_exp <- res$erlang_exp

    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    cat("==============================================\n")
    cat("     SMALLEST K RESULTS (ERLANG-EXP)          \n")
    cat("==============================================\n\n")

    cat("--- Search Settings ---\n")
    cat(sprintf("  Test Statistic : %s (%s)\n", res$pvaloption, stat_name))
    cat(sprintf("  Alpha          : %.4f\n", res$alpha))
    cat("\n")

    cat("--- Best K (MLE) ---\n")
    cat(sprintf("  K*             : %d\n", erlang_exp$best_k))
    cat(sprintf("  Erlang Lambda* : %.6f\n", erlang_exp$best_erlang_lambda))
    cat(sprintf("  Exp Lambda*    : %.6f\n", erlang_exp$best_exp_lambda))
    cat("\n")

    cat("--- Smallest K (passing GOF) ---\n")
    cat(sprintf("  K*             : %d\n", erlang_exp$smallest_k))
    cat(sprintf("  Erlang Lambda* : %.6f\n", erlang_exp$smallest_erlang_lambda))
    cat(sprintf("  Exp Lambda*    : %.6f\n", erlang_exp$smallest_exp_lambda))
    cat("\n")

    cat("--- GOF Results for Smallest K ---\n")
    cat(sprintf("  Observed Statistic : %.6f\n", erlang_exp$smallest_metric))
    cat(sprintf("  P-value            : %.6f\n", erlang_exp$smallest_p_value))
    cat(sprintf("  Q-value            : %d\n", erlang_exp$smallest_q_value))
    cat("\n")

    cat("--- Interpretation ---\n")
    if (erlang_exp$smallest_q_value == 1) {
      cat(sprintf("  Smallest K = %d passes the %s test at alpha = %.4f\n",
                  erlang_exp$smallest_k, stat_name, res$alpha))
    } else {
      cat(sprintf("  No K found that passes the %s test at alpha = %.4f\n",
                  stat_name, res$alpha))
    }
    cat("\n")

    cat("==============================================\n")
  })


  # =========================================================
  # SMALLEST K PDF PLOT (Single Model)
  # =========================================================

  output$smallest_k_pdf_plot <- renderPlot({

    req(smallest_k_results())
    req(gof_results())

    res <- smallest_k_results()
    gof_res <- gof_results()
    empiricaldata <- res$empiricaldata

    # Calculate bin width
    bin_width <- 2 * IQR(empiricaldata) / (length(empiricaldata)^(1/3))
    if (bin_width < .Machine$double.eps || !is.finite(bin_width)) {
      bin_width <- 1  # Default fallback bin width
    }

    # Create x grid for density curves
    x_grid <- seq(0, 1.2 * max(empiricaldata), length.out = 1000)

    # Base plot with histogram
    p <- ggplot(
      data.frame(Value = empiricaldata),
      aes(x = Value)
    ) +
      geom_histogram(
        aes(y = after_stat(density)),
        binwidth = bin_width,
        color = "black",
        fill = "#90C0AE",
        alpha = 0.7
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "bottom"
      )

    # =====================================================
    # ERLANG PDF
    # =====================================================
    if (res$fit_type == "Erlang") {

      # Best K density
      best_density <- dgamma(
        x_grid,
        shape = res$best_k,
        scale = 1/res$best_lambda
      )

      # Smallest K density
      smallest_density <- dgamma(
        x_grid,
        shape = res$smallest_k,
        scale = 1/res$smallest_lambda
      )

      # Create labels
      best_label <- sprintf("Best K (MLE): K=%d, λ=%.3f", res$best_k, res$best_lambda)
      smallest_label <- sprintf("Smallest K: K=%d, λ=%.3f", res$smallest_k, res$smallest_lambda)

      # Create data frames
      df_best <- data.frame(x = x_grid, density = best_density)
      df_smallest <- data.frame(x = x_grid, density = smallest_density)

      p <- p +
        geom_line(
          data = df_best,
          aes(x = x, y = density, color = "Best K (MLE)"),
          linewidth = 1.2
        ) +
        geom_line(
          data = df_smallest,
          aes(x = x, y = density, color = "Smallest K"),
          linewidth = 1.2,
          linetype = "dashed"
        ) +
        scale_color_manual(
          values = c(
            "Best K (MLE)" = "red",
            "Smallest K" = "darkgreen"
          ),
          labels = c(
            "Best K (MLE)" = best_label,
            "Smallest K" = smallest_label
          ),
          name = NULL
        ) +
        labs(
          title = "Data Histogram with Fitted Erlang Distributions",
          subtitle = "Comparing Best K (MLE) vs Smallest K",
          x = "Observed Values",
          y = "Density"
        )
    }

    # =====================================================
    # ERLANG-EXP PDF
    # =====================================================
    if (res$fit_type == "Erlang-Exp") {

      # Best K density
      best_density <- GenErlangFit:::ErlangExp_Func(
        x_grid,
        ErK = res$best_k,
        Erlam = res$best_erlang_lambda,
        Explam = res$best_exp_lambda
      )$Probability

      # Smallest K density
      smallest_density <- GenErlangFit:::ErlangExp_Func(
        x_grid,
        ErK = res$smallest_k,
        Erlam = res$smallest_erlang_lambda,
        Explam = res$smallest_exp_lambda
      )$Probability

      # Create labels
      best_label <- sprintf(
        "Best K (MLE): K=%d, λEr=%.3f, λExp=%.3f",
        res$best_k, res$best_erlang_lambda, res$best_exp_lambda
      )
      smallest_label <- sprintf(
        "Smallest K: K=%d, λEr=%.3f, λExp=%.3f",
        res$smallest_k, res$smallest_erlang_lambda, res$smallest_exp_lambda
      )

      # Create data frames
      df_best <- data.frame(x = x_grid, density = best_density)
      df_smallest <- data.frame(x = x_grid, density = smallest_density)

      p <- p +
        geom_line(
          data = df_best,
          aes(x = x, y = density, color = "Best K (MLE)"),
          linewidth = 1.2
        ) +
        geom_line(
          data = df_smallest,
          aes(x = x, y = density, color = "Smallest K"),
          linewidth = 1.2,
          linetype = "dashed"
        ) +
        scale_color_manual(
          values = c(
            "Best K (MLE)" = "blue",
            "Smallest K" = "purple"
          ),
          labels = c(
            "Best K (MLE)" = best_label,
            "Smallest K" = smallest_label
          ),
          name = NULL
        ) +
        labs(
          title = "Data Histogram with Fitted Erlang-Exp Distributions",
          subtitle = "Comparing Best K (MLE) vs Smallest K",
          x = "Observed Values",
          y = "Density"
        )
    }

    p
  })


  # =========================================================
  # SMALLEST K CDF PLOT (Single Model)
  # =========================================================

  output$smallest_k_cdf_plot <- renderPlot({

    req(smallest_k_results())
    req(gof_results())

    res <- smallest_k_results()
    gof_res <- gof_results()
    empiricaldata <- res$empiricaldata

    # Compute empirical CDF
    ecdf_data <- ecdf(empiricaldata)
    x_vals <- sort(empiricaldata)
    ecdf_vals <- ecdf_data(x_vals)

    # =====================================================
    # ERLANG CDF
    # =====================================================
    if (res$fit_type == "Erlang") {

      # Best K CDF
      best_cdf <- pgamma(x_vals, shape = res$best_k, scale = 1/res$best_lambda)

      # Smallest K CDF
      smallest_cdf <- pgamma(x_vals, shape = res$smallest_k, scale = 1/res$smallest_lambda)

      # Create labels
      best_label <- sprintf("Best K: K=%d, λ=%.3f", res$best_k, res$best_lambda)
      smallest_label <- sprintf("Smallest K: K=%d, λ=%.3f", res$smallest_k, res$smallest_lambda)

      # Build data frame
      df_cdf <- data.frame(
        x = rep(x_vals, 3),
        cdf = c(ecdf_vals, best_cdf, smallest_cdf),
        Type = factor(
          c(
            rep("Empirical CDF", length(x_vals)),
            rep(best_label, length(x_vals)),
            rep(smallest_label, length(x_vals))
          ),
          levels = c("Empirical CDF", best_label, smallest_label)
        )
      )

      # Create plot
      p <- ggplot(df_cdf, aes(x = x, y = cdf, color = Type, linetype = Type)) +
        geom_step(
          data = subset(df_cdf, Type == "Empirical CDF"),
          linewidth = 1.2
        ) +
        geom_line(
          data = subset(df_cdf, Type == best_label),
          linewidth = 1.2
        ) +
        geom_line(
          data = subset(df_cdf, Type == smallest_label),
          linewidth = 1.2
        ) +
        scale_color_manual(
          values = c(
            "Empirical CDF" = "black",
            setNames("red", best_label),
            setNames("darkgreen", smallest_label)
          )
        ) +
        scale_linetype_manual(
          values = c(
            "Empirical CDF" = "solid",
            setNames("solid", best_label),
            setNames("dashed", smallest_label)
          )
        ) +
        labs(
          title = "CDF Comparison (Erlang)",
          x = "x",
          y = "CDF",
          color = NULL,
          linetype = NULL
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin(t = -10),
          plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        guides(
          color = guide_legend(ncol = 1),
          linetype = guide_legend(ncol = 1)
        )
    }

    # =====================================================
    # ERLANG-EXP CDF
    # =====================================================
    if (res$fit_type == "Erlang-Exp") {

      # Best K CDF
      best_params <- c(res$best_k, res$best_erlang_lambda, res$best_exp_lambda)
      best_cdf <- GenErlangFit:::ErlangExpCDF_Func(best_params, x_vals, interval = 0.01)

      # Smallest K CDF
      smallest_params <- c(res$smallest_k, res$smallest_erlang_lambda, res$smallest_exp_lambda)
      smallest_cdf <- GenErlangFit:::ErlangExpCDF_Func(smallest_params, x_vals, interval = 0.01)

      # Create labels
      best_label <- sprintf(
        "Best K: K=%d, λEr=%.2f, λExp=%.2f",
        res$best_k, res$best_erlang_lambda, res$best_exp_lambda
      )
      smallest_label <- sprintf(
        "Smallest K: K=%d, λEr=%.2f, λExp=%.2f",
        res$smallest_k, res$smallest_erlang_lambda, res$smallest_exp_lambda
      )

      # Build data frame
      df_cdf <- data.frame(
        x = rep(x_vals, 3),
        cdf = c(ecdf_vals, best_cdf, smallest_cdf),
        Type = factor(
          c(
            rep("Empirical CDF", length(x_vals)),
            rep(best_label, length(x_vals)),
            rep(smallest_label, length(x_vals))
          ),
          levels = c("Empirical CDF", best_label, smallest_label)
        )
      )

      # Create plot
      p <- ggplot(df_cdf, aes(x = x, y = cdf, color = Type, linetype = Type)) +
        geom_step(
          data = subset(df_cdf, Type == "Empirical CDF"),
          linewidth = 1.2
        ) +
        geom_line(
          data = subset(df_cdf, Type == best_label),
          linewidth = 1.2
        ) +
        geom_line(
          data = subset(df_cdf, Type == smallest_label),
          linewidth = 1.2
        ) +
        scale_color_manual(
          values = c(
            "Empirical CDF" = "black",
            setNames("blue", best_label),
            setNames("purple", smallest_label)
          )
        ) +
        scale_linetype_manual(
          values = c(
            "Empirical CDF" = "solid",
            setNames("solid", best_label),
            setNames("dashed", smallest_label)
          )
        ) +
        labs(
          title = "CDF Comparison (Erlang-Exp)",
          x = "x",
          y = "CDF",
          color = NULL,
          linetype = NULL
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin(t = -10),
          plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        guides(
          color = guide_legend(ncol = 1),
          linetype = guide_legend(ncol = 1)
        )
    }

    p
  })


  # =========================================================
  # SMALLEST K BOOTSTRAP PLOT (Single Model)
  # =========================================================

  output$smallest_k_bootstrap_plot <- renderPlot({

    req(smallest_k_results())

    res <- smallest_k_results()

    # Get test statistic name for labels
    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    # Create data frame for histogram
    df_bootstrap <- data.frame(
      Statistic = res$smallest_sample_stats
    )

    # Calculate appropriate bin width
    bin_width <- diff(range(res$smallest_sample_stats)) / 30

    # Determine color based on fit type
    hist_fill <- if (res$fit_type == "Erlang") "#98D8AA" else "#B4A7D6"

    # Create plot
    p <- ggplot(df_bootstrap, aes(x = Statistic)) +
      geom_histogram(
        binwidth = bin_width,
        fill = hist_fill,
        color = "black",
        alpha = 0.7
      ) +
      geom_vline(
        xintercept = res$smallest_metric,
        linetype = "dashed",
        color = "black",
        linewidth = 1.2
      ) +
      annotate(
        "text",
        x = res$smallest_metric,
        y = Inf,
        label = sprintf("Observed = %.4f", res$smallest_metric),
        hjust = -0.1,
        vjust = 2,
        size = 4,
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = sprintf("p-value = %.4f", res$smallest_p_value),
        hjust = 1.1,
        vjust = 2,
        size = 4,
        fontface = "bold",
        color = if (res$smallest_q_value == 1) "darkgreen" else "red"
      ) +
      labs(
        title = sprintf("Bootstrap Distribution (%s)", stat_name),
        subtitle = sprintf(
          "Smallest K = %d | %s",
          res$smallest_k,
          res$fit_type
        ),
        x = sprintf("%s Statistic", stat_name),
        y = "Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40")
      )

    p
  })


  # =========================================================
  # SMALLEST K TEXT OUTPUT (Single Model)
  # =========================================================

  output$smallest_k_output <- renderPrint({

    req(smallest_k_results())

    res <- smallest_k_results()

    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    # =====================================================
    # ERLANG SMALLEST K OUTPUT
    # =====================================================
    if (res$fit_type == "Erlang") {

      cat("==============================================\n")
      cat("       SMALLEST K RESULTS (ERLANG)            \n")
      cat("==============================================\n\n")

      cat("--- Search Settings ---\n")
      cat(sprintf("  Test Statistic : %s (%s)\n", res$pvaloption, stat_name))
      cat(sprintf("  Alpha          : %.4f\n", res$alpha))
      cat("\n")

      cat("--- Best K (MLE) ---\n")
      cat(sprintf("  K*      : %d\n", res$best_k))
      cat(sprintf("  Lambda* : %.6f\n", res$best_lambda))
      cat("\n")

      cat("--- Smallest K (passing GOF) ---\n")
      cat(sprintf("  K*      : %d\n", res$smallest_k))
      cat(sprintf("  Lambda* : %.6f\n", res$smallest_lambda))
      cat("\n")

      cat("--- GOF Results for Smallest K ---\n")
      cat(sprintf("  Observed Statistic : %.6f\n", res$smallest_metric))
      cat(sprintf("  P-value            : %.6f\n", res$smallest_p_value))
      cat(sprintf("  Q-value            : %d\n", res$smallest_q_value))
      cat("\n")

      cat("--- Interpretation ---\n")
      if (res$smallest_q_value == 1) {
        cat(sprintf("  Smallest K = %d passes the %s test at alpha = %.4f\n",
                    res$smallest_k, stat_name, res$alpha))
      } else {
        cat(sprintf("  No K found that passes the %s test at alpha = %.4f\n",
                    stat_name, res$alpha))
      }
      cat("\n")

      cat("==============================================\n")
    }

    # =====================================================
    # ERLANG-EXP SMALLEST K OUTPUT
    # =====================================================
    if (res$fit_type == "Erlang-Exp") {

      cat("==============================================\n")
      cat("     SMALLEST K RESULTS (ERLANG-EXP)          \n")
      cat("==============================================\n\n")

      cat("--- Search Settings ---\n")
      cat(sprintf("  Test Statistic : %s (%s)\n", res$pvaloption, stat_name))
      cat(sprintf("  Alpha          : %.4f\n", res$alpha))
      cat("\n")

      cat("--- Best K (MLE) ---\n")
      cat(sprintf("  K*             : %d\n", res$best_k))
      cat(sprintf("  Erlang Lambda* : %.6f\n", res$best_erlang_lambda))
      cat(sprintf("  Exp Lambda*    : %.6f\n", res$best_exp_lambda))
      cat("\n")

      cat("--- Smallest K (passing GOF) ---\n")
      cat(sprintf("  K*             : %d\n", res$smallest_k))
      cat(sprintf("  Erlang Lambda* : %.6f\n", res$smallest_erlang_lambda))
      cat(sprintf("  Exp Lambda*    : %.6f\n", res$smallest_exp_lambda))
      cat("\n")

      cat("--- GOF Results for Smallest K ---\n")
      cat(sprintf("  Observed Statistic : %.6f\n", res$smallest_metric))
      cat(sprintf("  P-value            : %.6f\n", res$smallest_p_value))
      cat(sprintf("  Q-value            : %d\n", res$smallest_q_value))
      cat("\n")

      cat("--- Interpretation ---\n")
      if (res$smallest_q_value == 1) {
        cat(sprintf("  Smallest K = %d passes the %s test at alpha = %.4f\n",
                    res$smallest_k, stat_name, res$alpha))
      } else {
        cat(sprintf("  No K found that passes the %s test at alpha = %.4f\n",
                    stat_name, res$alpha))
      }
      cat("\n")

      cat("==============================================\n")
    }
  })


  # =========================================================
  # GOF COMBINED CDF PLOT (for Default fit type)
  # =========================================================

  output$gof_cdf_plot_combined <- renderPlot({

    req(gof_results())

    res <- gof_results()

    if (res$fit_type != "Default") {
      return(NULL)
    }

    empiricaldata <- res$empiricaldata

    # Compute empirical CDF
    ecdf_data <- ecdf(empiricaldata)
    x_vals <- sort(empiricaldata)
    ecdf_vals <- ecdf_data(x_vals)

    # Compute Erlang CDF
    erlang_cdf <- pgamma(
      x_vals,
      shape = res$erlang$k_star,
      scale = 1/res$erlang$lambda_star
    )

    # Compute Erlang-Exp CDF
    params_exp <- c(res$erlang_exp$k_star, res$erlang_exp$erlang_lambda_star, res$erlang_exp$exp_lambda_star)
    erlang_exp_cdf <- GenErlangFit:::ErlangExpCDF_Func(
      params_exp,
      x_vals,
      interval = 0.01
    )

    # Create labels
    erlang_label <- sprintf(
      "Erlang: K=%d, λ=%.3f",
      res$erlang$k_star,
      res$erlang$lambda_star
    )

    erlang_exp_label <- sprintf(
      "Erlang-Exp: K=%d, λEr=%.3f, λExp=%.3f",
      res$erlang_exp$k_star,
      res$erlang_exp$erlang_lambda_star,
      res$erlang_exp$exp_lambda_star
    )

    # Build data frame
    df_cdf <- data.frame(
      x = rep(x_vals, 3),
      cdf = c(ecdf_vals, erlang_cdf, erlang_exp_cdf),
      Type = factor(
        c(
          rep("Empirical CDF", length(x_vals)),
          rep(erlang_label, length(x_vals)),
          rep(erlang_exp_label, length(x_vals))
        ),
        levels = c("Empirical CDF", erlang_label, erlang_exp_label)
      )
    )

    # Create plot
    p <- ggplot(df_cdf, aes(x = x, y = cdf, color = Type, linetype = Type)) +
      geom_step(
        data = subset(df_cdf, Type == "Empirical CDF"),
        linewidth = 1.2
      ) +
      geom_line(
        data = subset(df_cdf, Type == erlang_label),
        linewidth = 1.2
      ) +
      geom_line(
        data = subset(df_cdf, Type == erlang_exp_label),
        linewidth = 1.2
      ) +
      scale_color_manual(
        values = c(
          "Empirical CDF" = "black",
          setNames("red", erlang_label),
          setNames("blue", erlang_exp_label)
        )
      ) +
      scale_linetype_manual(
        values = c(
          "Empirical CDF" = "solid",
          setNames("solid", erlang_label),
          setNames("dashed", erlang_exp_label)
        )
      ) +
      labs(
        title = "Empirical vs Fitted CDFs (Both Models)",
        x = "x",
        y = "CDF",
        color = NULL,
        linetype = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = -10),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      guides(
        color = guide_legend(ncol = 1),
        linetype = guide_legend(ncol = 1)
      )

    p
  })


  # =========================================================
  # GOF BOOTSTRAP PLOT - ERLANG (for Default fit type)
  # =========================================================

  output$gof_bootstrap_plot_erlang <- renderPlot({

    req(gof_results())

    res <- gof_results()

    if (res$fit_type != "Default") {
      return(NULL)
    }

    # Get test statistic name for labels
    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    # Create data frame for histogram
    df_bootstrap <- data.frame(
      Statistic = res$erlang$sample_stats
    )

    # Calculate appropriate bin width
    bin_width <- diff(range(res$erlang$sample_stats)) / 30

    # Create plot
    p <- ggplot(df_bootstrap, aes(x = Statistic)) +
      geom_histogram(
        binwidth = bin_width,
        fill = "#F8766D",
        color = "black",
        alpha = 0.7
      ) +
      geom_vline(
        xintercept = res$erlang$metric_star,
        linetype = "dashed",
        color = "black",
        linewidth = 1.2
      ) +
      annotate(
        "text",
        x = res$erlang$metric_star,
        y = Inf,
        label = sprintf("Observed = %.4f", res$erlang$metric_star),
        hjust = -0.1,
        vjust = 2,
        size = 4,
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = sprintf("p-value = %.4f", res$erlang$p_value),
        hjust = 1.1,
        vjust = 2,
        size = 4,
        fontface = "bold",
        color = if (res$erlang$q_value == 1) "darkgreen" else "red"
      ) +
      labs(
        title = sprintf("Erlang Bootstrap (%s)", stat_name),
        subtitle = sprintf("K = %d | n = %d bootstraps", res$erlang$k_star, res$n_bootstraps),
        x = sprintf("%s Statistic", stat_name),
        y = "Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40")
      )

    p
  })


  # =========================================================
  # GOF BOOTSTRAP PLOT - ERLANG-EXP (for Default fit type)
  # =========================================================

  output$gof_bootstrap_plot_erlang_exp <- renderPlot({

    req(gof_results())

    res <- gof_results()

    if (res$fit_type != "Default") {
      return(NULL)
    }

    # Get test statistic name for labels
    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    # Create data frame for histogram
    df_bootstrap <- data.frame(
      Statistic = res$erlang_exp$sample_stats
    )

    # Calculate appropriate bin width
    bin_width <- diff(range(res$erlang_exp$sample_stats)) / 30

    # Create plot
    p <- ggplot(df_bootstrap, aes(x = Statistic)) +
      geom_histogram(
        binwidth = bin_width,
        fill = "#619CFF",
        color = "black",
        alpha = 0.7
      ) +
      geom_vline(
        xintercept = res$erlang_exp$metric_star,
        linetype = "dashed",
        color = "black",
        linewidth = 1.2
      ) +
      annotate(
        "text",
        x = res$erlang_exp$metric_star,
        y = Inf,
        label = sprintf("Observed = %.4f", res$erlang_exp$metric_star),
        hjust = -0.1,
        vjust = 2,
        size = 4,
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = sprintf("p-value = %.4f", res$erlang_exp$p_value),
        hjust = 1.1,
        vjust = 2,
        size = 4,
        fontface = "bold",
        color = if (res$erlang_exp$q_value == 1) "darkgreen" else "red"
      ) +
      labs(
        title = sprintf("Erlang-Exp Bootstrap (%s)", stat_name),
        subtitle = sprintf("K = %d | n = %d bootstraps", res$erlang_exp$k_star, res$n_bootstraps),
        x = sprintf("%s Statistic", stat_name),
        y = "Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40")
      )

    p
  })


  # =========================================================
  # GOF TEXT OUTPUT - ERLANG (for Default fit type)
  # =========================================================

  output$gof_output_erlang <- renderPrint({

    req(gof_results())

    res <- gof_results()

    if (res$fit_type != "Default") {
      return(NULL)
    }

    cat("==============================================\n")
    cat("       GOODNESS OF FIT RESULTS (ERLANG)       \n")
    cat("==============================================\n\n")

    cat("--- Fitted Model Parameters ---\n")
    cat(sprintf("  K*      : %d\n", res$erlang$k_star))
    cat(sprintf("  Lambda* : %.6f\n", res$erlang$lambda_star))
    cat("\n")

    cat("--- GOF Test Settings ---\n")
    cat(sprintf("  Test Statistic : %s\n", res$pvaloption))
    cat(sprintf("  Alpha          : %.4f\n", res$alpha))
    cat(sprintf("  Bootstraps     : %d\n", res$n_bootstraps))
    cat("\n")

    cat("--- GOF Test Results ---\n")
    cat(sprintf("  Observed Statistic : %.6f\n", res$erlang$metric_star))
    cat(sprintf("  P-value            : %.6f\n", res$erlang$p_value))
    cat(sprintf("  Q-value            : %d\n", res$erlang$q_value))
    cat("\n")

    cat("--- Decision ---\n")
    if (res$erlang$q_value == 1) {
      cat(sprintf("  FAIL TO REJECT null hypothesis at alpha = %.4f\n", res$alpha))
      cat("  Interpretation: Data is CONSISTENT with the Erlang model.\n")
    } else {
      cat(sprintf("  REJECT null hypothesis at alpha = %.4f\n", res$alpha))
      cat("  Interpretation: Data is NOT consistent with the Erlang model.\n")
    }
    cat("\n")

    cat("==============================================\n")
  })


  # =========================================================
  # GOF TEXT OUTPUT - ERLANG-EXP (for Default fit type)
  # =========================================================

  output$gof_output_erlang_exp <- renderPrint({

    req(gof_results())

    res <- gof_results()

    if (res$fit_type != "Default") {
      return(NULL)
    }

    cat("==============================================\n")
    cat("     GOODNESS OF FIT RESULTS (ERLANG-EXP)     \n")
    cat("==============================================\n\n")

    cat("--- Fitted Model Parameters ---\n")
    cat(sprintf("  K*             : %d\n", res$erlang_exp$k_star))
    cat(sprintf("  Erlang Lambda* : %.6f\n", res$erlang_exp$erlang_lambda_star))
    cat(sprintf("  Exp Lambda*    : %.6f\n", res$erlang_exp$exp_lambda_star))
    cat("\n")

    cat("--- GOF Test Settings ---\n")
    cat(sprintf("  Test Statistic : %s\n", res$pvaloption))
    cat(sprintf("  Alpha          : %.4f\n", res$alpha))
    cat(sprintf("  Bootstraps     : %d\n", res$n_bootstraps))
    cat("\n")

    cat("--- GOF Test Results ---\n")
    cat(sprintf("  Observed Statistic : %.6f\n", res$erlang_exp$metric_star))
    cat(sprintf("  P-value            : %.6f\n", res$erlang_exp$p_value))
    cat(sprintf("  Q-value            : %d\n", res$erlang_exp$q_value))
    cat("\n")

    cat("--- Decision ---\n")
    if (res$erlang_exp$q_value == 1) {
      cat(sprintf("  FAIL TO REJECT null hypothesis at alpha = %.4f\n", res$alpha))
      cat("  Interpretation: Data is CONSISTENT with the Erlang-Exp model.\n")
    } else {
      cat(sprintf("  REJECT null hypothesis at alpha = %.4f\n", res$alpha))
      cat("  Interpretation: Data is NOT consistent with the Erlang-Exp model.\n")
    }
    cat("\n")

    cat("==============================================\n")
  })


  # =========================================================
  # GOF CDF COMPARISON PLOT (for single model)
  # =========================================================

  output$gof_cdf_plot <- renderPlot({

    req(gof_results())

    res <- gof_results()

    # Skip if Default fit type (handled by combined plot)
    if (res$fit_type == "Default") {
      return(NULL)
    }

    empiricaldata <- res$empiricaldata

    # Compute empirical CDF
    ecdf_data <- ecdf(empiricaldata)
    x_vals <- sort(empiricaldata)
    ecdf_vals <- ecdf_data(x_vals)

    # =====================================================
    # ERLANG CDF
    # =====================================================
    if (res$fit_type == "Erlang") {

      # Compute theoretical Erlang CDF
      theoretical_cdf <- pgamma(
        x_vals,
        shape = res$k_star,
        scale = 1/res$lambda_star
      )

      # Create label for legend
      fit_label <- sprintf(
        "Erlang CDF: K=%d, λ=%.3f",
        res$k_star,
        res$lambda_star
      )

      # Build data frame for plotting
      df_cdf <- data.frame(
        x = rep(x_vals, 2),
        cdf = c(ecdf_vals, theoretical_cdf),
        Type = factor(
          c(
            rep("Empirical CDF", length(x_vals)),
            rep(fit_label, length(x_vals))
          ),
          levels = c("Empirical CDF", fit_label)
        )
      )

      # Create plot
      p <- ggplot(df_cdf, aes(x = x, y = cdf, color = Type)) +
        geom_step(
          data = subset(df_cdf, Type == "Empirical CDF"),
          linewidth = 1.2
        ) +
        geom_line(
          data = subset(df_cdf, Type == fit_label),
          linewidth = 1.2
        ) +
        scale_color_manual(
          values = c(
            "Empirical CDF" = "black",
            setNames("red", fit_label)
          )
        ) +
        labs(
          title = "Empirical vs Fitted CDF (Erlang)",
          x = "x",
          y = "CDF",
          color = NULL
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin(t = -10),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
    }

    # =====================================================
    # ERLANG-EXP CDF
    # =====================================================
    if (res$fit_type == "Erlang-Exp") {

      # Compute theoretical Erlang-Exp CDF using custom function
      params <- c(res$k_star, res$erlang_lambda_star, res$exp_lambda_star)
      theoretical_cdf <- GenErlangFit:::ErlangExpCDF_Func(
        params,
        x_vals,
        interval = 0.01
      )

      # Create label for legend
      fit_label <- sprintf(
        "Erlang-Exp CDF: K=%d, λEr=%.3f, λExp=%.3f",
        res$k_star,
        res$erlang_lambda_star,
        res$exp_lambda_star
      )

      # Build data frame for plotting
      df_cdf <- data.frame(
        x = rep(x_vals, 2),
        cdf = c(ecdf_vals, theoretical_cdf),
        Type = factor(
          c(
            rep("Empirical CDF", length(x_vals)),
            rep(fit_label, length(x_vals))
          ),
          levels = c("Empirical CDF", fit_label)
        )
      )

      # Create plot
      p <- ggplot(df_cdf, aes(x = x, y = cdf, color = Type)) +
        geom_step(
          data = subset(df_cdf, Type == "Empirical CDF"),
          linewidth = 1.2
        ) +
        geom_line(
          data = subset(df_cdf, Type == fit_label),
          linewidth = 1.2
        ) +
        scale_color_manual(
          values = c(
            "Empirical CDF" = "black",
            setNames("blue", fit_label)
          )
        ) +
        labs(
          title = "Empirical vs Fitted CDF (Erlang-Exp)",
          x = "x",
          y = "CDF",
          color = NULL
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin(t = -10),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
    }

    p
  })


  # =========================================================
  # GOF BOOTSTRAP HISTOGRAM PLOT (for single model)
  # =========================================================

  output$gof_bootstrap_plot <- renderPlot({

    req(gof_results())

    res <- gof_results()

    # Skip if Default fit type (handled by separate plots)
    if (res$fit_type == "Default") {
      return(NULL)
    }

    # Get test statistic name for labels
    stat_name <- switch(
      toupper(res$pvaloption),
      "KS" = "Kolmogorov-Smirnov",
      "AD" = "Anderson-Darling",
      "CVM" = "Cramér-von Mises",
      res$pvaloption
    )

    # Create data frame for histogram
    df_bootstrap <- data.frame(
      Statistic = res$sample_stats
    )

    # Calculate appropriate bin width
    bin_width <- diff(range(res$sample_stats)) / 30

    # Determine color based on fit type
    hist_fill <- if (res$fit_type == "Erlang") "#F8766D" else "#619CFF"

    # Create plot
    p <- ggplot(df_bootstrap, aes(x = Statistic)) +
      geom_histogram(
        binwidth = bin_width,
        fill = hist_fill,
        color = "black",
        alpha = 0.7
      ) +
      geom_vline(
        xintercept = res$metric_star,
        linetype = "dashed",
        color = "black",
        linewidth = 1.2
      ) +
      annotate(
        "text",
        x = res$metric_star,
        y = Inf,
        label = sprintf("Observed = %.4f", res$metric_star),
        hjust = -0.1,
        vjust = 2,
        size = 4,
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = sprintf("p-value = %.4f", res$p_value),
        hjust = 1.1,
        vjust = 2,
        size = 4,
        fontface = "bold",
        color = if (res$q_value == 1) "darkgreen" else "red"
      ) +
      labs(
        title = sprintf("Bootstrap Distribution (%s)", stat_name),
        subtitle = sprintf(
          "%s | n = %d bootstraps",
          res$fit_type,
          res$n_bootstraps
        ),
        x = sprintf("%s Statistic", stat_name),
        y = "Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40")
      )

    p
  })


  # =========================================================
  # GOF TEXT OUTPUT (for single model)
  # =========================================================

  output$gof_output <- renderPrint({

    req(gof_results())

    res <- gof_results()

    # Skip if Default fit type (handled by separate outputs)
    if (res$fit_type == "Default") {
      return(NULL)
    }

    # =====================================================
    # ERLANG OUTPUT
    # =====================================================
    if (res$fit_type == "Erlang") {

      cat("==============================================\n")
      cat("       GOODNESS OF FIT RESULTS (ERLANG)       \n")
      cat("==============================================\n\n")

      cat("--- Fitted Model Parameters ---\n")
      cat(sprintf("  K*      : %d\n", res$k_star))
      cat(sprintf("  Lambda* : %.6f\n", res$lambda_star))
      cat("\n")

      cat("--- GOF Test Settings ---\n")
      cat(sprintf("  Test Statistic : %s\n", res$pvaloption))
      cat(sprintf("  Alpha          : %.4f\n", res$alpha))
      cat(sprintf("  Bootstraps     : %d\n", res$n_bootstraps))
      cat("\n")

      cat("--- GOF Test Results ---\n")
      cat(sprintf("  Observed Statistic : %.6f\n", res$metric_star))
      cat(sprintf("  P-value            : %.6f\n", res$p_value))
      cat(sprintf("  Q-value            : %d\n", res$q_value))
      cat("\n")

      cat("--- Decision ---\n")
      if (res$q_value == 1) {
        cat(sprintf("  FAIL TO REJECT null hypothesis at alpha = %.4f\n", res$alpha))
        cat("  Interpretation: Data is CONSISTENT with the Erlang model.\n")
      } else {
        cat(sprintf("  REJECT null hypothesis at alpha = %.4f\n", res$alpha))
        cat("  Interpretation: Data is NOT consistent with the Erlang model.\n")
      }
      cat("\n")

      cat("--- Bootstrap Sample Statistics Summary ---\n")
      cat(sprintf("  Min    : %.6f\n", min(res$sample_stats)))
      cat(sprintf("  Q1     : %.6f\n", quantile(res$sample_stats, 0.25)))
      cat(sprintf("  Median : %.6f\n", median(res$sample_stats)))
      cat(sprintf("  Q3     : %.6f\n", quantile(res$sample_stats, 0.75)))
      cat(sprintf("  Max    : %.6f\n", max(res$sample_stats)))
      cat("\n")

      cat("==============================================\n")
    }

    # =====================================================
    # ERLANG-EXP OUTPUT
    # =====================================================
    if (res$fit_type == "Erlang-Exp") {

      cat("==============================================\n")
      cat("     GOODNESS OF FIT RESULTS (ERLANG-EXP)     \n")
      cat("==============================================\n\n")

      cat("--- Fitted Model Parameters ---\n")
      cat(sprintf("  K*             : %d\n", res$k_star))
      cat(sprintf("  Erlang Lambda* : %.6f\n", res$erlang_lambda_star))
      cat(sprintf("  Exp Lambda*    : %.6f\n", res$exp_lambda_star))
      cat("\n")

      cat("--- GOF Test Settings ---\n")
      cat(sprintf("  Test Statistic : %s\n", res$pvaloption))
      cat(sprintf("  Alpha          : %.4f\n", res$alpha))
      cat(sprintf("  Bootstraps     : %d\n", res$n_bootstraps))
      cat("\n")

      cat("--- GOF Test Results ---\n")
      cat(sprintf("  Observed Statistic : %.6f\n", res$metric_star))
      cat(sprintf("  P-value            : %.6f\n", res$p_value))
      cat(sprintf("  Q-value            : %d\n", res$q_value))
      cat("\n")

      cat("--- Decision ---\n")
      if (res$q_value == 1) {
        cat(sprintf("  FAIL TO REJECT null hypothesis at alpha = %.4f\n", res$alpha))
        cat("  Interpretation: Data is CONSISTENT with the Erlang-Exp model.\n")
      } else {
        cat(sprintf("  REJECT null hypothesis at alpha = %.4f\n", res$alpha))
        cat("  Interpretation: Data is NOT consistent with the Erlang-Exp model.\n")
      }
      cat("\n")

      cat("--- Bootstrap Sample Statistics Summary ---\n")
      cat(sprintf("  Min    : %.6f\n", min(res$sample_stats)))
      cat(sprintf("  Q1     : %.6f\n", quantile(res$sample_stats, 0.25)))
      cat(sprintf("  Median : %.6f\n", median(res$sample_stats)))
      cat(sprintf("  Q3     : %.6f\n", quantile(res$sample_stats, 0.75)))
      cat(sprintf("  Max    : %.6f\n", max(res$sample_stats)))
      cat("\n")

      cat("==============================================\n")
    }
  })


  # =========================================================
  # DATA HISTOGRAM
  # =========================================================

  output$data_histogram <- renderPlot({

    req(data())

    empiricaldata <- data()[[1]]


    bin_width <- 2 * IQR(empiricaldata) /
      (length(empiricaldata)^(1/3))

    if (bin_width < .Machine$double.eps || !is.finite(bin_width)) {
      bin_width <- 1  # Default fallback bin width
    }

    ggplot(
      data.frame(Value = empiricaldata),
      aes(x = Value)
    ) +
      geom_histogram(
        binwidth = bin_width,
        color = "black",
        fill = "#90C0AE"
      ) +
      labs(
        title = "Empirical Data Histogram",
        x = "Observed Values",
        y = "Frequency"
      ) +
      theme_minimal(base_size = 15)
  })


  # =========================================================
  # FIT PLOT
  # =========================================================

  output$fit_plot <- renderPlot({

    req(data())

    empiricaldata <- data()[[1]]


    bin_width <- 2 * IQR(empiricaldata) /
      (length(empiricaldata)^(1/3))

    if (bin_width < .Machine$double.eps || !is.finite(bin_width)) {
      bin_width <- 1  # Default fallback bin width
    }

    p <- ggplot(
      data.frame(Value = empiricaldata),
      aes(x = Value)
    ) +
      geom_histogram(
        aes(y = after_stat(density)),
        binwidth = bin_width,
        color = "black",
        fill = "#90C0AE"
      ) +
      labs(
        title = "Empirical Data with Fitted Distribution",
        x = "Observed Values",
        y = "Density"
      ) +
      theme_minimal(base_size = 15)


    if (input$run_fit > 0) {

      req(fit_results())


      x_grid <- seq(
        0,
        1.2 * max(empiricaldata),
        length.out = 1000
      )


      # ===================================================
      # DEFAULT FIT
      # ===================================================
      if (input$fit_type == "Default") {

        # Erlang Results
        K_erlang <- fit_results()$Erlang_Results$Best$K_star

        lambda_erlang <- fit_results()$Erlang_Results$Best$Lambda_star


        erlang_density <- dgamma(
          x_grid,
          shape = K_erlang,
          scale = 1/lambda_erlang
        )


        erlang_df <- data.frame(
          x = x_grid,
          density = erlang_density
        )

        # Erlang-Exp Results
        K_ErExp <- fit_results()$ErlangExp_Results$Best$K_star

        erlang_lambda_ErExp <- fit_results()$ErlangExp_Results$Best$ErlangLambda_star

        exp_lambda_ErExp <- fit_results()$ErlangExp_Results$Best$ExpLambda_star


        density_exp <-
          GenErlangFit:::ErlangExp_Func(
            x_grid,
            ErK = K_ErExp,
            Erlam = erlang_lambda_ErExp,
            Explam = exp_lambda_ErExp
          )$Probability

        df_exp <- data.frame(
          x = x_grid,
          density = density_exp,
          Model = paste0(
            "Erlang-Exp: K = ",
            K_ErExp,
            ", λEr = ",
            round(erlang_lambda_ErExp, 3),
            ", λExp = ",
            round(exp_lambda_ErExp, 3)
          )
        )



        p <- p +

          geom_line(
            data = erlang_df,
            aes(x = x, y = density, color = "Erlang"),
            linewidth = 1.5
          ) +

          geom_line(
            data = df_exp,
            aes(x = x, y = density, color = "Erlang-Exp"),
            linewidth = 1.5
          ) +

          scale_color_manual(
            values = c(
              "Erlang" = "red",
              "Erlang-Exp" = "blue"
            ),
            name = "Fit Type"
          )
      }



      # ===================================================
      # ERLANG FIT
      # ===================================================
      if (input$fit_type == "Erlang") {

        K_fit <- fit_results()$Best$K_star

        lambda_fit <- fit_results()$Best$Lambda_star


        fitted_density <- dgamma(
          x_grid,
          shape = K_fit,
          scale = 1/lambda_fit
        )


        fit_df <- data.frame(
          x = x_grid,
          density = fitted_density
        )


        p <- p +
          geom_line(
            data = fit_df,
            aes(x = x, y = density, color = "Erlang"),
            linewidth = 1.5
          )


        if (isTRUE(input$find_smallest_erlang)) {

          K_small <- fit_results()$Smallest$K_star

          lambda_small <- fit_results()$Smallest$Lambda_star


          smallest_density <- dgamma(
            x_grid,
            shape = K_small,
            scale = 1/lambda_small
          )


          smallest_df <- data.frame(
            x = x_grid,
            density = smallest_density
          )


          p <- p +
            geom_line(
              data = smallest_df,
              aes(x = x, y = density, color = "Erlang Smallest K"),
              linewidth = 1.5
            )
        }


        p <- p +
          scale_color_manual(
            values = c(
              "Erlang" = "red",
              "Erlang Smallest K" = "darkgreen"
            ),
            name = "Fit Type"
          )
      }



      # ===================================================
      # ERLANG-EXP FIT
      # ===================================================
      if (input$fit_type == "Erlang-Exp") {

        K_ErExp <- fit_results()$Best$K_star

        erlang_lambda_ErExp <- fit_results()$Best$ErlangLambda_star

        exp_lambda_ErExp <- fit_results()$Best$ExpLambda_star


        density_exp <-
          GenErlangFit:::ErlangExp_Func(
            x_grid,
            ErK = K_ErExp,
            Erlam = erlang_lambda_ErExp,
            Explam = exp_lambda_ErExp
          )$Probability

        df_exp <- data.frame(
          x = x_grid,
          density = density_exp,
          Model = paste0(
            "Erlang-Exp: K = ",
            K_ErExp,
            ", λEr = ",
            round(erlang_lambda_ErExp, 3),
            ", λExp = ",
            round(exp_lambda_ErExp, 3)
          )
        )


        p <- p +
          geom_line(
            data = df_exp,
            aes(x = x, y = density, color = "Erlang-Exp"),
            linewidth = 1.5
          )


        if (isTRUE(input$find_smallest_erlang_exp)) {

          K_small <- fit_results()$Smallest$K_star

          erlang_lambda_small <- fit_results()$Smallest$ErlangLambda_star

          exp_lambda_small <- fit_results()$Smallest$ExpLambda_star


          density_small <-
            GenErlangFit:::ErlangExp_Func(
              x_grid,
              ErK = K_small,
              Erlam = erlang_lambda_small,
              Explam = exp_lambda_small
            )$Probability

          smallest_df <- data.frame(
            x = x_grid,
            density = density_small
          )


          p <- p +
            geom_line(
              data = smallest_df,
              aes(x = x, y = density, color = "Erlang-Exp Smallest K"),
              linewidth = 1.5
            )
        }


        p <- p +
          scale_color_manual(
            values = c(
              "Erlang-Exp" = "blue",
              "Erlang-Exp Smallest K" = "purple"
            ),
            name = "Fit Type"
          )
      }
    }


    p
  })

}



#------------------------
# RUN APP
#------------------------

shinyApp(ui = ui, server = server)
