#------------------------
# Section 0: Packages
#------------------------
library(shiny)
library(bslib)
library(GenErlangFit)

#------------------------
# UI
#------------------------

ui <- navbarPage(
  
  title = "GenErlangFit",
  
  theme = bs_theme(
    version = 5,
    bootswatch = "minty",
    navbar_bg = "#79c2ad"
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
            
            h3("Current Dataset"),
            
            tableOutput("data_table")
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
        
        verbatimTextOutput("fit_output")
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
        
        # Fixed K
        if (input$search_type == "Fixed K") {
          
          valid <- TRUE
        }
        
        
        # Search Window
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
  # FIT OUTPUT
  # =========================================================
  
  output$fit_output <- renderPrint({
    
    input$run_fit
    
    isolate({
      
      cat("Selected Fit Type:\n\n")
      
      cat("-", input$fit_type, "\n")
      
      
      # -----------------------------------------------------
      # ERLANG
      # -----------------------------------------------------
      if (input$fit_type == "Erlang") {
        
        if (!is.na(input$initial_k)) {
          
          cat(
            "\nInitial K Guess:",
            input$initial_k,
            "\n"
          )
        }
        
        cat(
          "Find Erlang Smallest K:",
          input$find_smallest_erlang,
          "\n"
        )
      }
      
      
      # -----------------------------------------------------
      # ERLANG-EXP
      # -----------------------------------------------------
      if (input$fit_type == "Erlang-Exp") {
        
        cat(
          "\nInitial K Guess:",
          input$initial_k_exp,
          "\n"
        )
        
        cat(
          "Search Type:",
          input$search_type,
          "\n"
        )
        
        
        if (
          input$search_type ==
          "Search over a Window"
        ) {
          
          cat(
            "Window Size:",
            input$window_size,
            "\n"
          )
        }
        
        
        cat(
          "Find Erlang-Exp Smallest K:",
          input$find_smallest_erlang_exp,
          "\n"
        )
      }
    })
  })
  
  
  # =========================================================
  # DYNAMIC GOF BUTTON
  # =========================================================
  
  output$run_gof_button <- renderUI({
    
    valid <- FALSE
    
    
    # -------------------------------------------------------
    # DEFAULT
    # -------------------------------------------------------
    if (input$gof_mode == "Default") {
      
      valid <- TRUE
    }
    
    
    # -------------------------------------------------------
    # USER SELECTION
    # -------------------------------------------------------
    if (input$gof_mode == "User Selection") {
      
      if (
        !is.na(input$alpha_value) &&
        !is.na(input$num_bootstraps) &&
        length(input$test_statistics) > 0
      ) {
        
        valid <- TRUE
      }
    }
    
    
    # -------------------------------------------------------
    # BUTTON STATE
    # -------------------------------------------------------
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
  # GOF OUTPUT
  # =========================================================
  
  output$gof_output <- renderPrint({
    
    input$run_gof
    
    isolate({
      
      cat("Selected GOF Mode:\n\n")
      
      cat("-", input$gof_mode, "\n")
      
      
      # -----------------------------------------------------
      # USER SELECTION SETTINGS
      # -----------------------------------------------------
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
  
  
  # =========================================================
  # DATA TABLE
  # =========================================================
  
  output$data_table <- renderTable({
    data()
  })
  
}



#------------------------
# RUN APP
#------------------------

shinyApp(ui = ui, server = server)