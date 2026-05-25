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
  # 1. DATA ENTRY
  # =========================================================
  tabPanel(
    "Data Entry",

    sidebarLayout(

      # -----------------------------------------------------
      # SIDEBAR
      # -----------------------------------------------------
      sidebarPanel(

        h4("Data Controls"),

        fileInput(
          "file",
          "Upload CSV"
        ),

        actionButton(
          "clear",
          "Clear Data"
        )
      ),


      # -----------------------------------------------------
      # MAIN PANEL
      # -----------------------------------------------------
      mainPanel(

        tabsetPanel(


          # =================================================
          # ENTER DATA
          # =================================================
          tabPanel(
            "Enter Data",

            h3("Manual Data Entry"),

            textAreaInput(
              "manual_data",
              "Paste or type data here:",
              rows = 10
            ),

            br(),

            actionButton(
              "toggle_format",
              "Show Required Format"
            ),


            # ----------------------------------------------
            # REQUIRED FORMAT BOX
            # ----------------------------------------------
            conditionalPanel(
              condition = "input.toggle_format % 2 == 1",

              div(
                style = "
                  margin-top:10px;
                  padding:12px;
                  border:1px solid #ddd;
                  border-radius:8px;
                  background:#f8f9fa;
                ",

                strong("Required Format"),

                tags$pre(
                  "Manually enter rows or upload a CSV file in the following format:

Dataset name,Time units,Case counts"
                )
              )
            )
          ),


          # =================================================
          # VIEW DATA
          # =================================================
          tabPanel(
            "View Data",

            h3("Uploaded Dataset"),

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
            value = NA,
            min = 0,
            max = 1,
            step = 0.01
          ),


          # -----------------------------------------------
          # NUMBER OF BOOTSTRAPS
          # -----------------------------------------------
          numericInput(
            "num_bootstraps",
            "Number of Bootstraps",
            value = NA,
            min = 1,
            step = 1
          ),


          # -----------------------------------------------
          # TEST STATISTICS
          # -----------------------------------------------
          checkboxGroupInput(
            "test_statistics",
            "Choice of Test Statistic",
            choices = c(
              "Likelihood-Ratio Test",
              "Distance Measure - KS",
              "Distance Measure - AD",
              "Distance Measure - CvM"
            ),
            selected = NULL
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

    valid <- FALSE


    if (input$gof_mode == "Default") {

      valid <- TRUE
    }


    if (input$gof_mode == "User Selection") {

      if (
        !is.na(input$alpha_value) &&
        !is.na(input$num_bootstraps) &&
        length(input$test_statistics) > 0
      ) {

        valid <- TRUE
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


  # =========================================================
  # GOF OUTPUT
  # =========================================================

  output$gof_output <- renderPrint({

    input$run_gof

    isolate({

      cat("Selected GOF Mode:\n\n")

      cat("-", input$gof_mode, "\n")


      if (input$gof_mode == "User Selection") {

        cat(
          "\nAlpha Value:",
          input$alpha_value,
          "\n"
        )

        cat(
          "Number of Bootstraps:",
          input$num_bootstraps,
          "\n"
        )

        cat(
          "\nSelected Test Statistics:\n"
        )

        for (i in input$test_statistics) {

          cat("-", i, "\n")
        }
      }
    })
  })

}



#------------------------
# RUN APP
#------------------------

shinyApp(ui = ui, server = server)
