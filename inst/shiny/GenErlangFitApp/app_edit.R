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
        uiOutput("run_gof_button")
      ),


      # -----------------------------------------------------
      # MAIN PANEL
      # -----------------------------------------------------
      mainPanel(

        h3("Goodness of Fit Results"),

        br(),

        verbatimTextOutput("gof_output")
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
        sample_stats = gof_res$sample_stats
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
        sample_stats = gof_res$sample_stats
      ))
    }

    return(NULL)
  })


  # =========================================================
  # GOF OUTPUT
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
