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
          checkboxInput(
            "find_smallest_erlang",
            "Find Erlang Smallest K",
            value = FALSE
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
          checkboxInput(
            "find_smallest_erlang_exp",
            "Find Erlang-Exp Smallest K",
            value = FALSE
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

    fit_results()$ResultsTable
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

    # Check for Erlang or Erlang-Exp fit types
    if (fit_exists && input$fit_type %in% c("Erlang", "Erlang-Exp")) {

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

      return(
        div(
          style = "color: #856404; background-color: #fff3cd; padding: 10px; border-radius: 5px;",
          icon("info-circle"),
          " Please select Erlang or Erlang-Exp fit type for GOF computation."
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
        lambda_star = lambda_star,
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
  # SMALLEST K SECTION (appears after GOF is run)
  # =========================================================

  output$smallest_k_section <- renderUI({

    # Check if GOF results exist
    gof_exists <- !is.null(tryCatch(gof_results(), error = function(e) NULL))

    if (!gof_exists) {
      return(NULL)
    }

    # Get the test statistic used
    res <- gof_results()

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
  # SMALLEST K COMPUTATION
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
  # SMALLEST K OUTPUT SECTION
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
  # SMALLEST K PDF PLOT (Plot 1)
  # =========================================================

  output$smallest_k_pdf_plot <- renderPlot({

    req(smallest_k_results())
    req(gof_results())

    res <- smallest_k_results()
    gof_res <- gof_results()
    empiricaldata <- res$empiricaldata

    # Calculate bin width
    bin_width <- 2 * IQR(empiricaldata) / (length(empiricaldata)^(1/3))

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
        scale = res$best_lambda
      )

      # Smallest K density
      smallest_density <- dgamma(
        x_grid,
        shape = res$smallest_k,
        scale = res$smallest_lambda
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
        "Best K (MLE): K=%d, λE=%.3f, λX=%.3f",
        res$best_k, res$best_erlang_lambda, res$best_exp_lambda
      )
      smallest_label <- sprintf(
        "Smallest K: K=%d, λE=%.3f, λX=%.3f",
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
  # SMALLEST K CDF PLOT (Plot 2)
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
      best_cdf <- pgamma(x_vals, shape = res$best_k, scale = res$best_lambda)

      # Smallest K CDF
      smallest_cdf <- pgamma(x_vals, shape = res$smallest_k, scale = res$smallest_lambda)

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
        "Best K: K=%d, λE=%.2f, λX=%.2f",
        res$best_k, res$best_erlang_lambda, res$best_exp_lambda
      )
      smallest_label <- sprintf(
        "Smallest K: K=%d, λE=%.2f, λX=%.2f",
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
  # SMALLEST K BOOTSTRAP PLOT (Plot 3)
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
  # SMALLEST K TEXT OUTPUT
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
  # GOF CDF COMPARISON PLOT
  # =========================================================

  output$gof_cdf_plot <- renderPlot({

    req(gof_results())

    res <- gof_results()
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
        scale = res$lambda_star
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
        "Erlang-Exp CDF: K=%d, λE=%.3f, λX=%.3f",
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
  # GOF BOOTSTRAP HISTOGRAM PLOT
  # =========================================================

  output$gof_bootstrap_plot <- renderPlot({

    req(gof_results())

    res <- gof_results()

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
  # GOF TEXT OUTPUT
  # =========================================================

  output$gof_output <- renderPrint({

    req(gof_results())

    res <- gof_results()

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
          scale = lambda_erlang
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
            ", λE = ",
            round(erlang_lambda_ErExp, 3),
            ", λX = ",
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
          scale = lambda_fit
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
            scale = lambda_small
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
            ", λE = ",
            round(erlang_lambda_ErExp, 3),
            ", λX = ",
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
