#' Fit Empirical Data Using Erlang + Exponential Mixture Model (Version 2)
#'
#' Fits empirical data to a mixture of Erlang and Exponential distributions
#' using Maximum Likelihood Estimation (MLE), with options for goodness-of-fit
#' testing via bootstrapping and adaptive search over shape parameter `K`.
#'
#' @param empiricaldata Numeric vector of observed data.
#' @param K Integer initial guess for the Erlang shape parameter.
#' @param ... Optional named arguments to override defaults:
#' \describe{
#'   \item{Alpha}{Significance level for goodness-of-fit test. Default = 0.05.}
#'   \item{pvaloption}{Goodness-of-fit metric: `"KS"`, `"CvM"`, `"AD"`, or `"NIL"`. Default = `"KS"`.}
#'   \item{InitialguessErLam}{Initial guess for Erlang rate parameter. Default = 3.}
#'   \item{SmallestK}{Logical. If TRUE, searches for smallest acceptable K under alpha. Default = FALSE.}
#'   \item{SmallestKValue}{Initial smallest K value for search. Default = -1 (auto).}
#'   \item{ShowFigures}{Logical. If TRUE, generates diagnostic plots. Default = TRUE.}
#'   \item{FixedK}{Logical. If TRUE, fits only the specified K without adaptive search. Default = FALSE.}
#'   \item{KWindowSize}{Window size around K for adaptive search. Default = 1.}
#'   \item{NumBootstraps}{Number of bootstrap samples for p-value estimation. Default = round(10 * 10 / Alpha).}
#' }
#'
#' @details
#' This function fits an Erlang + Exponential mixture model to the input data,
#' optionally performs bootstrap-based goodness-of-fit tests, and can generate
#' diagnostic plots of the fit. It supports an adaptive search over the shape
#' parameter `K` or fitting with fixed `K`.
#'
#' @return A list containing:
#' \describe{
#'   \item{Best}{A list with the best fit parameters, p-value, and goodness-of-fit metrics.}
#'   \item{Smallest}{If \code{SmallestK = TRUE}, results for the smallest acceptable K.}
#'   \item{AllFits}{All fit results across the searched K window (unless \code{FixedK = TRUE}).}
#' }
#'
#' @seealso [ErlangExp_Fit_v2_FixedK()], [ErlangExp_Fit_v2_Pvalue()], [ErlangExp_NegLogLikelihood()]
#'
#' @examples
#' \dontrun{
#' data <- rexp(200, rate = 0.5)
#' fit <- ErlangExp_Fit_v2(data, K = 2, Alpha = 0.05, ShowFigures = TRUE)
#' print(fit$Best)
#' }
#'
#' @keywords internal
ErlangExp_Fit_v2 <- function(empiricaldata, K, ...) {


  defaultAlpha = 0.05;
  defaultSmallestKValue = -1;

  # -----------------------------
  # 1. Defaults
  # -----------------------------
  defaults <- list(
    Alpha = defaultAlpha,
    pvaloption = "KS",
    InitialguessErLam = 3,
    SmallestK = FALSE,
    SmallestKValue = defaultSmallestKValue,
    ShowFigures = TRUE,
    FixedK = FALSE,
    KWindowSize = 1,
    NumBootstraps = round(10 * 10 / defaultAlpha)
  )

  # -----------------------------
  # 2. Parse name–value arguments (case-insensitive) - i.e., override default parameter values when provided as input options
  # -----------------------------
  user_args <- list(...)
  if (length(user_args) > 0) {
    names_lower <- tolower(names(user_args))
    defaults_lower <- tolower(names(defaults))
    invalid <- setdiff(names_lower, defaults_lower)
    if (length(invalid) > 0)
      stop(paste("Unknown parameter(s):", paste(invalid, collapse = ", ")))
    matched <- match(names_lower, defaults_lower)
    names(user_args) <- names(defaults)[matched]
  }

  # Override default parameter values when provided as input options
  options <- modifyList(defaults, user_args)

  # Update Numbootstraps to whatever Alpha was provided
  if (is.null(user_args$NumBootstraps)) {
    options$NumBootstraps <- round(10 * 10 / options$Alpha)
  }

  OGPVal <- options$pvaloption # Store this because we will turn this off when searching through multiple K values.

  # Parse out options for use in THIS code.
  FixedK <- options$FixedK
  KWindowSize <- options$KWindowSize
  SmallestK = options$SmallestK
  SmallestKValue = options$SmallestKValue
  ShowFigures = options$ShowFigures
  pvaloption <- options$pvaloption
  cat("Options:\n"); print(options)

  # Completely Remove options that ErlangExp_Fit doesn’t use
  suboptions <- options
  suboptions$FixedK <- NULL
  suboptions$KWindowSize <- NULL
  suboptions$SmallestKValue <- NULL

  cat("Sub-Options:\n"); print(options)

  # -----------------------------
  # 3. If FixedK is specified, fit is done for that K
  # -----------------------------
  if (FixedK) {
    res_fixed <- do.call(ErlangExp_Fit_v2_FixedK, c(list(empiricaldata, K), suboptions))
    best_row = res_fixed
    combined = res_fixed
    best_row_list = res_fixed

  }

  else {
    # -----------------------------
    # 4. Otherwise adaptive multi-K search around given K
    # -----------------------------
    ErExpK <- K
    WindowK <- seq(ErExpK - KWindowSize, ErExpK + KWindowSize)
    WindowK <- WindowK[WindowK >= 0] # Ensures that ErExp Ks are >= 0.

    suboptions$pvaloption <- "NIL"  # speed up inner fits by not computing the P-value.

    ErlangExpBest <- matrix(NA, nrow = length(WindowK), ncol = 4)
    colnames(ErlangExpBest) <- c("K_star", "ErlangLambda_star", "ExpLambda_star", "LogLikelihood")

    for (i in seq_along(WindowK)) {
      cat(sprintf("Running K Window %d / %d ...\n", i, length(WindowK)))
      k_try <- WindowK[i]
      # suboptions$InitialguessErLam <- k_try / mean(empiricaldata) # MODIFIED
      res <- do.call(ErlangExp_Fit_v2_FixedK, c(list(empiricaldata, k_try), suboptions))
      ErlangExpBest[i, ] <- c(k_try, res$ErlangLambda, res$ExpLambda, res$LogLikelihood)
    }

    LLs <- ErlangExpBest[, "LogLikelihood"]
    maxIdx <- which.max(LLs)
    maxLL <- LLs[maxIdx]

    # -----------------------------
    # Adaptive expansion logic
    # -----------------------------
    if (maxIdx > 1 && maxIdx < nrow(ErlangExpBest)) {
      cat("Peak found in middle. Search complete.\n")
      combined <- ErlangExpBest
      best_row <- combined[maxIdx, ]

    } else if (maxIdx == 1 && WindowK[maxIdx] == 0) {
      cat("Peak is an exponential function (K = 0). Search complete.\n")
      combined <- ErlangExpBest
      best_row <- combined[maxIdx, ]

    } else if (maxIdx == 1 && WindowK[maxIdx] > 0) {
      cat("Extending window to the left ...\n")
      criteria <- 1
      count <- 0
      currentK <- WindowK[maxIdx]
      currentML <- maxLL
      shift <- 0
      ErlangExpBest_LeftWindow <- NULL

      while (criteria > 0) {
        count <- count + 1
        newK <- currentK - 1
        if (newK < 0) {
          shift <- 1
          break
        }
        # suboptions$InitialguessErLam <- newK / mean(empiricaldata) # MODIFIED
        res <- do.call(ErlangExp_v2_Fit_FixedK, c(list(empiricaldata, newK), suboptions))
        new_row <- c(newK, res$ErlangLambda, res$ExpLambda, res$LogLikelihood)
        ErlangExpBest_LeftWindow <- rbind(ErlangExpBest_LeftWindow, new_row)
        criteria <- new_row[4] > currentML
        currentML <- new_row[4]
        currentK <- newK
      }

      combined <- rbind(ErlangExpBest_LeftWindow[nrow(ErlangExpBest_LeftWindow):1, ], ErlangExpBest)
      best_row <- combined[1 + 1 - shift, ]

    } else if (maxIdx == nrow(ErlangExpBest)) {
      cat("Extending window to the right ...\n")
      criteria <- 1
      count <- 0
      currentK <- WindowK[maxIdx]
      currentML <- maxLL
      ErlangExpBest_RightWindow <- NULL

      while (criteria > 0) {
        count <- count + 1
        newK <- currentK + 1
        # suboptions$InitialguessErLam <- newK / mean(empiricaldata) # MODIFIED
        res <- do.call(ErlangExp_Fit_v2_FixedK, c(list(empiricaldata, newK), suboptions))
        new_row <- c(newK, res$ErlangLambda, res$ExpLambda, res$LogLikelihood)
        ErlangExpBest_RightWindow <- rbind(ErlangExpBest_RightWindow, new_row)
        criteria <- new_row[4] > currentML
        currentML <- new_row[4]
        currentK <- newK
      }

      combined <- rbind(ErlangExpBest, ErlangExpBest_RightWindow)
      best_row <- combined[nrow(combined) - 1, ]
    }


    # Compute P-Value for best row identified in the window
    if (tolower(OGPVal) == "nil"){ # User does not want p-Value

      bestK <- best_row[1]
      bestERLambda <- best_row[2]
      bestEXPLambda <- best_row[3]
      maxLL <- best_row[4]


      cat(sprintf(
        "Best K = %d | Erlang λ = %.4f | Exp λ = %.4f | Log-Likelihood = %.4f\n",
        bestK, bestERLambda, bestEXPLambda, maxLL
      ))

      best_row_list <- list(
        K_star = best_row[[1]],
        ErlangLambda_star = best_row[[2]],
        ExpLambda_star = best_row[[3]],
        LogLikelihood = best_row[[4]]
      )

    }

    else{ # User wants p-Value for best case

      # Compute p-value and q-value
      PQres <- ErlangExp_Fit_v2_Pvalue(empiricaldata,
                                       best_row[1],
                                       best_row[2],
                                       best_row[3],
                                       s = length(empiricaldata),
                                       n = options$NumBootstraps,
                                       alpha = options$Alpha,
                                       pvaloption = OGPVal,
                                       ShowFigures = options$ShowFigures)
      p_value <- PQres$p_value
      q_value <- PQres$q_value
      gof_metric <- PQres$metric_star
      gof_Sample_Stats <- PQres$sample_stats

      # Append p and q to best_row
      best_row <- c(best_row, p_value, q_value)
      names(best_row) <- c("K_star", "ErlangLambda_star", "ExpLambda_star", "LogLikelihood", "P_Value", "Q_Value")

      cat(sprintf(
        "Best K = %d | Erlang λ = %.4f | Exp λ = %.4f | P = %.4f | Q = %.4f | Log-Likelihood = %.4f\n",
        best_row[1], best_row[2], best_row[3], best_row[5], best_row[6], best_row[4]
      ))

      best_row_list <- list(
        K_star = best_row[[1]],
        ErlangLambda_star = best_row[[2]],
        ExpLambda_star = best_row[[3]],
        P_star = PQres$p_value,
        Q_Value = PQres$q_value,
        LogLikelihood = best_row[[4]],
        metric_star = PQres$metric_star,
        samplestats_star = PQres$sample_stats
      )

    }

    print(best_row_list)


  }

  if (SmallestK) {

    # If SmallestKValue is not provided, compute it via Erlang_Fit
    if (SmallestKValue < 0) {
      # Run Erlang_Fit to get the initial smallest K
      ErlangSmallest <- Erlang_Fit_v2(empiricaldata, SmallestK = TRUE, ShowFigures = FALSE)
      SmallestKValue <- ErlangSmallest$Smallest$K_star
    }


    # Quick check: does the current best fit pass? If it doesnt then return current value.
    if (length(best_row) >= 5 && best_row[5] == 0) {
      SmallestKOutput <- best_row
      cat("No K smaller than the current Best Fit K Found.\n")
    } else {

      # Compute P-value for the initial guess for Smallest K Value
      suboptions$ShowFigures <- FALSE  # avoid intermediate plots
      suboptions$pvaloption <- OGPVal  # avoid intermediate plots
      args_to_call <- c(list(empiricaldata), SmallestKValue, suboptions)

      ErlangExpSmallest <- list()
      ErlangExpSmallest[[1]] <- do.call(ErlangExp_Fit_v2_FixedK, args_to_call)

      print( ErlangExpSmallest[[1]])
      # Quick check: does the current smallest fit pass? If it doesnt then return current value. If not run iteratively down.
      if (ErlangExpSmallest[[1]]$Q_Value == 0) {
        SmallestKOutput <- ErlangExpSmallest[[1]]
        cat("No K smaller than initial guess for Smallest K Value Found. Try a higher value of K\n")

      }  else {

        k_temp <- ErlangExpSmallest[[1]]$K_star
        q_temp <- ErlangExpSmallest[[1]]$Q_Value
        metric_star_store <- c(ErlangExpSmallest[[1]]$metric_star)
        sampleStatsStore <- as.matrix(ErlangExpSmallest[[1]]$samplestats_star)
        p_temp <- c(ErlangExpSmallest[[1]]$p_value)

        count <- 2

        while (q_temp > 0) {
          k_temp <- k_temp - 1
          if (k_temp < 1) break
          args_to_call[[2]] <- k_temp
          ErlangExpSmallest[[count]] <- do.call(ErlangExp_Fit_v2_FixedK, args_to_call)
          q_temp <- ErlangExpSmallest[[count]]$Q_Value

          # Wanted to store these to be able to plot but since i am alr storing the whole thing will just call smallestkoutput
          # p_temp <- c(p_temp, res_temp$P_star)
          # q_temp <- res_temp$Q_Value
          # metric_star_store <- c(metric_star_store, res_temp$metric_star)
          # sampleStatsStore <- cbind(sampleStatsStore, res_temp$samplestats_star)

          count <- count + 1
        }

        # Select the last valid K
        if (k_temp < 1) {
          SmallestKOutput <- ErlangExpSmallest[[length(ErlangExpSmallest)]]
          # gof_metric_smallest <- metric_star_store[1]
          # gof_Sample_Stats_smallest <- sampleStatsStore[, 1, drop = FALSE]
        } else {
          SmallestKOutput <- ErlangExpSmallest[[length(ErlangExpSmallest) - 1]] # ERROR flagged
          # gof_metric_smallest <- metric_star_store[length(metric_star_store) - 1]
          # gof_Sample_Stats_smallest <- sampleStatsStore[, ncol(sampleStatsStore) - 1, drop = FALSE]
        }

      }





    }


    cat("Smallest K search complete.\n")


    SmallestKOutput_list <- list(
      K_star = SmallestKOutput$K_star,
      ErlangLambda_star = SmallestKOutput$ErlangLambda_star,
      ExpLambda_star = SmallestKOutput$ExpLambda_star,
      P_star = SmallestKOutput$P_star,
      Q_Value = SmallestKOutput$Q_Value,
      LogLikelihood = SmallestKOutput$LogLikelihood,
      metric_star = SmallestKOutput$metric_star,
      samplestats_star = SmallestKOutput$samplestats_star
    )


  }




  if (ShowFigures) {


    K = best_row_list$K_star;
    ErlangLambda =  best_row_list$ErlangLambda_star
    ExpLambda = best_row_list$ExpLambda_star
    metric_star = best_row_list$metric_star
    sampleStats = best_row_list$samplestats_star

    if (SmallestK) {


      k_smallest = SmallestKOutput$K_star
      lambda_smallest = SmallestKOutput$ErlangLambda_star
      ExpLambda_small = SmallestKOutput$ExpLambda_star



    }



    pvaloption_upper <- toupper(pvaloption)

    # Plot main PDF

    bin_width <- 2 * IQR(empiricaldata) / (length(empiricaldata)^(1/3))

    # Plotting Function
    xax <- seq(0, 1.2*max(empiricaldata), by = 0.01)

    resplot <- ErlangExp_Func(xax, ErK = K, Erlam = ErlangLambda, Explam = ExpLambda)

    # Create dataframe for plotting
    df_main <- data.frame(
      x = xax,
      y = resplot$Probability,
      Type = paste0("ErlangExp PDF: K* = ", K, ", λ* = ", round(ErlangLambda, 2),
                    ", λ_exp = ", round(ExpLambda, 2))
    )

    main_label <- df_main$Type[1]

    # Base plot: histogram + main fit
    P1 <- ggplot() +
      geom_histogram(data = data.frame(x = empiricaldata),
                     aes(x = x, y = after_stat(density), fill = "Observed Data"),
                     binwidth = bin_width,
                     color = "black", alpha = 0.6) +
      geom_line(data = df_main, aes(x = x, y = y, color = main_label),
                linewidth = 1.2) +
      scale_fill_manual(name = "", values = c("Observed Data" = "#78A5A3")) +
      scale_color_manual(name = "", values = setNames(c("red"), main_label)) +
      labs(title = "Histogram with Erlang-Exponential Fits",
           x = "x", y = "Density") +
      theme_gray() +
      theme(
        legend.position = "bottom",        # Move legend below the plot
        legend.box = "horizontal",         # Horizontal layout
        legend.title = element_blank(),    # Remove legend title
        legend.justification = "center",   # Center the legend
        plot.title = element_text(hjust = 0.5) # Center title
      )


    # Overlay smallest K if it exists
    if (SmallestK) {
      res_small <- ErlangExp_Func(xax, ErK = k_smallest, Erlam = lambda_smallest, Explam = ExpLambda_small)

      df_small <- data.frame(
        x = xax,
        y = res_small$Probability,
        Type = paste0("Smallest K PDF: K = ", k_smallest,
                      ", λ = ", round(lambda_smallest, 2),
                      ", λ_exp = ", round(ExpLambda_small, 2))
      )

      P1 <- P1 +
        geom_line(data = df_small, aes(x = x, y = y, color = Type),
                  linewidth = 1.2) +
        scale_color_manual(
          name = "",
          values = setNames(c("red", "grey50"), c(main_label, df_small$Type[1]))
        )
    }

    print(P1)


    # Plot main CDF

    # Base CDF data
    ecdf_data <- ecdf(empiricaldata)
    x_vals <- sort(empiricaldata)
    datY <- ecdf_data(x_vals)

    # Best-fit Erlang+Exponential CDF
    params <- c(K, ErlangLambda, ExpLambda)
    gammaY <- ErlangExpCDF_Func(params, x_vals, interval = 0.01)

    df_main_cdf <- data.frame(
      x = x_vals,
      y = gammaY,
      Type = paste0("ErlangExp CDF: K* = ", K,
                    ", λ* = ", round(ErlangLambda, 2),
                    ", λ_exp = ", round(ExpLambda, 2))
    )

    main_label_CDF <- df_main_cdf$Type[1]

    # Base plot: Empirical CDF + main fit
    P2 <- ggplot() +
      geom_step(aes(x = x_vals, y = datY, color = "Empirical CDF"), linewidth = 1.1) +
      geom_line(data = df_main_cdf, aes(x = x, y = y, color = main_label_CDF),
                linewidth = 1.2) +
      scale_color_manual(name = "",
                         values = setNames(c("red", "black"), c(main_label_CDF, "Empirical CDF"))) +
      labs(title = "Empirical vs Erlang-Exponential CDF",
           x = "x", y = "CDF") +
      theme_gray() +
      theme(
        legend.position = "bottom",        # Move legend below the plot
        legend.box = "horizontal",         # Horizontal layout
        legend.title = element_blank(),    # Remove legend title
        legend.justification = "center",   # Center the legend
        plot.title = element_text(hjust = 0.5) # Center title
      )
    # Overlay smallest K if it exists
    if (SmallestK) {
      params_small <- c(k_smallest, lambda_smallest, ExpLambda_small)
      gammaY_small <- ErlangExpCDF_Func(params_small, x_vals, interval = 0.01)

      df_small_cdf <- data.frame(
        x = x_vals,
        y = gammaY_small,
        Type = paste0("Smallest K CDF: K = ", k_smallest,
                      ", λ = ", round(lambda_smallest, 2),
                      ", λ_exp = ", round(ExpLambda_small, 2))
      )

      P2 <- P2 +
        geom_line(data = df_small_cdf, aes(x = x, y = y, color = Type), linewidth = 1.2) +
        scale_color_manual(name = "",
                           values = setNames(c("red", "black", "grey50"),
                                             c(main_label_CDF, "Empirical CDF", df_small_cdf$Type[1])))
    }
    print(P2)


    # Plot main bootstaps if pvaloption is on.

    if (tolower(OGPVal) != "nil") {


      df_stats <- data.frame(Statistic = sampleStats)

      P3 <- ggplot(df_stats, aes(x = Statistic)) +
        geom_histogram(binwidth = diff(range(sampleStats)) / 30,
                       fill = "#78A5A3", alpha = 0.6, color = "black") +
        geom_vline(xintercept = metric_star, linetype = "dashed", color = "black", linewidth = 1) +
        labs(title = paste("Bootstrap Distribution of", pvaloption_upper, "Statistic"),
             x = paste(pvaloption_upper, "Statistic"), y = "Count") +
        theme_gray()
      # Optionally print
      print(P3)

      # Smallest K fit (if it exists)
      if (SmallestK) {
        sampleStats_small <- SmallestKOutput$samplestats_star
        metric_star_small <- SmallestKOutput$metric_star
        df_stats_small <- data.frame(Statistic = sampleStats_small)

        P4 <- ggplot(df_stats_small, aes(x = Statistic)) +
          geom_histogram(binwidth = diff(range(sampleStats_small)) / 30,
                         fill = "#78A5A3", alpha = 0.6, color = "black") +
          geom_vline(xintercept = metric_star_small, linetype = "dashed", color = "black", linewidth = 1) +
          labs(title = paste("Bootstrap Distribution of", pvaloption_upper, "Statistic (Smallest K)"),
               x = paste(pvaloption_upper, "Statistic"), y = "Count") +
          theme_gray()


        print(P4)

      }

    }


  }

  return(list(
    Best = best_row_list,
    Smallest = if (SmallestK) SmallestKOutput else NULL,
    AllFits = if (FixedK) NULL else combined
  ))

}


#' Fit Erlang + Exponential Mixture Model with Fixed Shape Parameter K
#'
#' Fits empirical data to an Erlang + Exponential mixture model using
#' Maximum Likelihood Estimation (MLE), with a fixed Erlang shape parameter `K`.
#' Optionally computes bootstrap-based goodness-of-fit p-values.
#'
#' @param empiricaldata Numeric vector of observed data.
#' @param FixedKValue Integer specifying the fixed Erlang shape parameter `K`.
#' @param ... Optional named arguments to override defaults:
#' \describe{
#'   \item{Alpha}{Significance level for goodness-of-fit test. Default = 0.05.}
#'   \item{pvaloption}{Goodness-of-fit metric: `"KS"`, `"CvM"`, `"AD"`, or `"NIL"` (skip p-value). Default = `"KS"`.}
#'   \item{InitialguessErlam}{Initial guess for Erlang rate parameter. Default = 3.}
#'   \item{SmallestK}{Logical. Not used in this function but included for consistency. Default = FALSE.}
#'   \item{ShowFigures}{Logical. If TRUE, generates diagnostic plots. Default = TRUE.}
#'   \item{NumBootstraps}{Number of bootstrap samples for p-value estimation. Default = round(10 * 10 / Alpha).}
#' }
#'
#' @details
#' This function performs MLE to estimate the Erlang and Exponential rate
#' parameters, fixing the Erlang shape parameter `K`. Goodness-of-fit testing
#' via bootstrap is optional and controlled by the `pvaloption` argument.
#'
#' @return A list containing:
#' \describe{
#'   \item{K_star}{Estimated Erlang shape parameter (fixed input).}
#'   \item{ErlangLambda_star}{Estimated Erlang rate parameter.}
#'   \item{ExpLambda_star}{Estimated Exponential rate parameter.}
#'   \item{P_star}{Bootstrap p-value for goodness-of-fit (if computed).}
#'   \item{Q_Value}{Bootstrap q-value (if computed).}
#'   \item{LogLikelihood}{Log-likelihood of the fitted model.}
#'   \item{metric_star}{Goodness-of-fit test statistic (if computed).}
#'   \item{samplestats_star}{Bootstrap sample statistics (if computed).}
#' }
#'
#' @seealso [ErlangExp_Fit_v2()], [ErlangExp_Fit_v2_Pvalue()], [ErlangExp_NegLogLikelihood()]
#'
#' @examples
#' \dontrun{
#' data <- rexp(200, rate = 0.5)
#' fit_fixed <- ErlangExp_Fit_v2_FixedK(data, FixedKValue = 2, Alpha = 0.05)
#' print(fit_fixed$ErlangLambda_star)
#' }
#'
#' @keywords internal
ErlangExp_Fit_v2_FixedK <- function(empiricaldata, FixedKValue, ...) {

  # -----------------------------
  # 1. Default Input Arguments
  # -----------------------------

  defaultAlpha = 0.05;

  defaults <- list(
    Alpha = defaultAlpha,
    pvaloption = "KS",
    InitialguessErlam = 3,
    SmallestK = FALSE,
    ShowFigures = TRUE,
    NumBootstraps = round(10 * 10 / defaultAlpha)
  )

  # -----------------------------
  # 2. Parse User Inputs
  # -----------------------------
  user_args <- list(...)

  if (length(user_args) > 0) {
    # Convert all user-supplied names to lowercase - so case insensitive
    user_names_lower <- tolower(names(user_args))
    default_names_lower <- tolower(names(defaults))

    # Check for invalid names
    invalid <- setdiff(user_names_lower, default_names_lower)
    if (length(invalid) > 0) {
      stop(paste0("Unknown parameter(s): ", paste(invalid, collapse = ", ")))
    }

    # Match lowercased user args to defaults
    matched_names <- match(user_names_lower, default_names_lower)
    names(user_args) <- names(defaults)[matched_names]
  }

  # -----------------------------
  # 3. Override defaults with user-supplied values
  # -----------------------------
  opts <- modifyList(defaults, user_args)

  # -----------------------------
  # 4. Post-processing or dependent defaults
  # -----------------------------
  if (is.null(user_args$NumBootstraps)) {
    opts$NumBootstraps <- round(10 * 10 / opts$Alpha)
  }

  # -----------------------------
  # 5. Assign to variables (for convenience)
  # -----------------------------
  Alpha <- opts$Alpha
  pvaloption <- opts$pvaloption
  InitialguessErlam <- opts$InitialguessErlam
  SmallestK <- opts$SmallestK
  ShowFigures <- opts$ShowFigures
  NumBootstraps <- opts$NumBootstraps

  # -----------------------------
  # 6. Display configuration summary
  # -----------------------------
  message("ErlangExp Fit initialized with user-parsed parameters:")
  print(opts)



  # -----------------------------
  # Optimize to find MLE
  # -----------------------------

  data_vec = empiricaldata

  # Initial guesses
  Erlam_init <- InitialguessErlam

  # Use 1D optimizer since we're only searching over Erlam for a given K
  opt <- optimize(f = ErlangExp_NegLogLikelihood,
                  interval = c(FixedKValue/mean(data_vec)+1e-2, 10*mean(data_vec)),  # adjust bounds if needed
                  tol = 1e-8,
                  data = data_vec,
                  K = FixedKValue)

  k_mle <- FixedKValue
  erlam_mle <- opt$minimum
  explam_mle <- 1 / (mean(data_vec) - (FixedKValue / erlam_mle))
  cat("MLE shape (k):", FixedKValue, "\n")
  cat("MLE scale (Erlang lambda):", erlam_mle, "\n")
  cat("MLE scale (Exp lambda):", explam_mle, "\n")
  resu <- ErlangExp_Func(data_vec, ErK = k_mle, Erlam = erlam_mle, Explam = explam_mle)
  cat("MLE Log-Likelihood:", resu$Likelihood, "\n")

  # -----------------------------
  # P. Value computation
  # -----------------------------
  if (tolower(pvaloption) == "nil") {
    # Skip p-value computation (used for internal multi-K searches)
    cat("Skipping p-value computation (pvaloption = NIL)\n")

  } else {

    # Using the function
    res <- ErlangExp_Fit_v2_Pvalue(empiricaldata,
                                   k_mle,
                                   erlam_mle,
                                   explam_mle,
                                   s = length(data_vec),
                                   n = NumBootstraps,
                                   alpha = Alpha,
                                   pvaloption = pvaloption,
                                   ShowFigures)
    p_value <- res$p_value
    q_value <- res$q_value
    gof_metric <- res$metric_star
    gof_Sample_Stats <- res$sample_stats
  }



  # # -----------------------------
  # # 8. Plot necessary figures  - Diagnostic
  # # -----------------------------
  #
  # bin_width <- 2 * IQR(empiricaldata) / (length(empiricaldata)^(1/3))
  #
  # # Plotting Function
  # xax <- seq(0, 1.2*max(empiricaldata), by = 0.01)
  # resplot <- ErlangExp_Func(xax, ErK = k_mle, Erlam = erlam_mle, Explam = explam_mle)
  #
  # # Create dataframe for plotting
  # df_plot <- data.frame(
  #   x = xax,
  #   pdf = resplot$Probability
  # )
  #
  # # Plot histogram + PDF
  # fit_fig1 <- ggplot() +
  #   geom_histogram(data = data.frame(x = empiricaldata),
  #                  aes(x = x, y = after_stat(density)),
  #                  binwidth = bin_width,
  #                  fill = "#78A5A3", color = "black", alpha = 0.6) +
  #   geom_line(data = df_plot, aes(x = x, y = pdf),
  #             color = "black", linewidth = 1.2) +
  #   labs(
  #     title = "Erlang+Exponential PDF",
  #     x = "x",
  #     y = "Density"
  #   ) +
  #   theme_minimal(base_size = 14)
  # print(fit_fig1)
  #
  #
  # if (ShowFigures) {
  #   print(fit_fig1)  # always available
  #
  #   if (tolower(pvaloption) != "nil") {
  #     # Only print GOF figures if p-value was computed
  #     print(gof_fig1)
  #     print(gof_fig2)
  #   }
  # }



  # -----------------------------
  # 10. Return the proper outputs
  # -----------------------------

  # return(list(
  #   K_star = k_star,
  #   Lambda_star = lambda_star,
  #   P_star = p_value,
  #   Q_Value = q_value,
  #   Loglikelihood = log_likelihood_star,
  #   Figures_star = list(Fig1 = fit_fig1, Fig2 = gof_fig1, Fig3 = gof_fig2)
  # ))


  if (tolower(pvaloption) == "nil")  {
    output_list <- list(
      K_star = k_mle,
      ErlangLambda_star = erlam_mle,
      ExpLambda_star = explam_mle,
      LogLikelihood = resu$Likelihood
    )
  } else {
    output_list <- list(
      K_star = k_mle,
      ErlangLambda_star = erlam_mle,
      ExpLambda_star = explam_mle,
      P_star = p_value,
      Q_Value = q_value,
      LogLikelihood = resu$Likelihood,
      metric_star = gof_metric,
      samplestats_star = gof_Sample_Stats

    )
  }

  return(output_list)

}



#' Bootstrap-based Goodness-of-Fit Test for Erlang + Exponential Fit
#'
#' Performs a bootstrap hypothesis test to assess goodness-of-fit of empirical data
#' to an Erlang + Exponential mixture model, using the specified test statistic
#' (Kolmogorov-Smirnov, Cramér-von Mises, or Anderson-Darling).
#'
#' @param empiricaldata Numeric vector of observed data.
#' @param k_star Integer Erlang shape parameter (fixed).
#' @param erlambda_star Numeric Erlang rate parameter.
#' @param explambda_star Numeric Exponential rate parameter.
#' @param s Integer sample size for bootstrap samples. Default is length of empiricaldata.
#' @param n Integer number of bootstrap samples to generate. Default is 1000.
#' @param alpha Significance level for hypothesis testing. Default is 0.05.
#' @param pvaloption Character specifying goodness-of-fit metric to use: `"KS"` (Kolmogorov-Smirnov), `"CvM"` (Cramér-von Mises), or `"AD"` (Anderson-Darling). Default is `"KS"`.
#' @param ShowFigures Logical indicating whether to generate diagnostic plots. Default is TRUE (currently commented out in code).
#'
#' @details
#' This function computes the empirical CDF and theoretical Erlang+Exponential CDF,
#' calculates the specified goodness-of-fit test statistic, and then uses bootstrap
#' resampling to estimate the p-value for the test. It also provides a simple
#' hypothesis test decision based on the significance level `alpha`.
#'
#' @return A list containing:
#' \describe{
#'   \item{p_value}{Bootstrap p-value for the goodness-of-fit test.}
#'   \item{q_value}{Binary decision indicator: 1 = fail to reject null, 0 = reject null.}
#'   \item{metric_star}{Observed test statistic value from empirical data.}
#'   \item{sample_stats}{Vector of test statistics from bootstrap samples.}
#' }
#'
#' @seealso [ErlangExp_Fit_v2()], [ErlangExp_Fit_v2_FixedK()], [ErlangExpCDF_Func()]
#'
#' @examples
#' \dontrun{
#' data <- rexp(200, rate = 0.5)
#' fit <- ErlangExp_Fit_v2(data)
#' pval_res <- ErlangExp_Fit_v2_Pvalue(data, fit$Best$K_star, fit$Best$ErlangLambda_star, fit$Best$ExpLambda_star)
#' print(pval_res$p_value)
#' }
#'
#' @keywords internal
ErlangExp_Fit_v2_Pvalue <- function(empiricaldata, k_star, erlambda_star, explambda_star,
                                    s = length(empiricaldata),
                                    n = 1000,            # number of bootstrap samples
                                    alpha = 0.05,
                                    pvaloption = "KS",   # only KS implemented for now
                                    ShowFigures = TRUE) {

  data_vec <- empiricaldata
  pvaloption_upper <- toupper(pvaloption)

  # 1. Empirical CDF
  ecdf_data <- ecdf(data_vec)
  x_vals <- sort(data_vec)
  datY <- ecdf_data(x_vals)
  n_points <- length(datY) # for AD

  # 2. Erlang CDF (theoretical)
  params <- c(k_star, erlambda_star, explambda_star)
  gammaY <- ErlangExpCDF_Func(params, x_vals, interval = 0.01)

  # 3. Observed test statistic
  if (pvaloption_upper == "KS") {
    metric_star <- max((datY - gammaY)^2)
  } else if (pvaloption_upper == "CVM") {
    metric_star <- sum((datY - gammaY)^2)
  } else if (pvaloption_upper == "AD") {
    weights <- ((2 * (1:n_points) - 1) / n_points)
    gammaY <- pmin(pmax(gammaY, 1e-10), 1 - 1e-10)
    metric_star <- -n_points - sum(weights * (log(gammaY) + log(1 - rev(gammaY))))
  } else {
    stop("Invalid p-value option. Choose 'KS', 'CvM', or 'AD'.")
  }

  # 4. Bootstrap loop
  sampleStats <- numeric(n)
  for (i in 1:n) {
    # sample_data <- rgamma(s, shape = k_star, scale = lambda_star) # CHANGE THIS TO GENERATE RANDOM ERLANG-EXP DATA ############
    sample_data <- ErlangExpRnd(params, s)
    ecdf_sample <- ecdf(sample_data)
    x_sample <- sort(sample_data)
    sampleY <- ecdf_sample(x_sample)
    # gammaY_sample <- pgamma(x_sample, shape = k_star, scale = lambda_star)
    gammaY_sample <- ErlangExpCDF_Func(params, x_sample, interval = 0.01) # CHANGE THIS TO GET ERLANG-EXP CDF FOR SAMPLE DATA ############

    if (pvaloption_upper == "KS") {
      sampleStats[i] <- max((sampleY - gammaY_sample)^2)
    } else if (pvaloption_upper == "CVM") {
      sampleStats[i] <- sum((sampleY - gammaY_sample)^2)
    } else if (pvaloption_upper == "AD") {
      weights_sample <- ((2 * (1:length(sampleY)) - 1) / length(sampleY))
      gammaY_sample <- pmin(pmax(gammaY_sample, 1e-10), 1 - 1e-10)
      sampleStats[i] <- -length(sampleY) - sum(weights_sample * (log(gammaY_sample) + log(1 - rev(gammaY_sample))))
    }
  }

  # 5. Compute p-value
  p_star <- mean(sampleStats >= metric_star)

  # 6. Compute Q-value
  if (p_star > alpha) {
    q_star <- 1  # fail to reject null
    message("Fail to reject the null hypothesis at alpha = ", alpha,
            ". Data consistent with the model.")
  } else {
    q_star <- 0  # reject null
    message("Reject the null hypothesis at alpha = ", alpha,
            ". Data not consistent with the model.")
  }

  # -----------------------------
  # 6. Figures
  # -----------------------------
  # P2 <- P3 <- NULL
  # if (ShowFigures) {
  #   # Empirical vs Erlang CDF
  #   df_cdf <- data.frame(
  #     x = x_vals,
  #     Empirical = datY,
  #     Erlang = gammaY
  #   )
  #   P2 <- ggplot(df_cdf, aes(x = x)) +
  #     geom_step(aes(y = Empirical, color = "Empirical CDF")) +
  #     geom_line(aes(y = Erlang, color = "Erlang CDF"), linewidth = 1.2) +
  #     scale_color_manual(name = "", values = c("Empirical CDF" = "black", "Erlang CDF" = "red")) +
  #     labs(title = "Empirical vs Erlang CDF",
  #          x = "x", y = "CDF") +
  #     theme_minimal()
  #
  #   # Histogram of bootstrap statistics
  #   df_stats <- data.frame(Statistic = sampleStats)
  #   P3 <- ggplot(df_stats, aes(x = Statistic)) +
  #     geom_histogram(binwidth = diff(range(sampleStats)) / 30,
  #                    fill = "#78A5A3", alpha = 0.6, color = "black") +
  #     geom_vline(xintercept = metric_star, linetype = "dashed", color = "black", linewidth = 1) +
  #     labs(title = paste("Bootstrap Distribution of", pvaloption_upper, "Statistic"),
  #          x = paste(pvaloption_upper, "Statistic"), y = "Count") +
  #     theme_minimal()
  # }

  # -----------------------------
  # 7. Return
  # -----------------------------
  return(list(
    p_value = p_star,
    q_value = q_star,
    metric_star = metric_star,
    sample_stats = sampleStats
  ))

}


#' Negative Log-Likelihood for Erlang + Exponential Mixture
#'
#' Computes the negative log-likelihood of the Erlang + Exponential mixture
#' distribution given the Erlang rate parameter, data, and Erlang shape parameter.
#'
#' @param ErlangLambda Numeric rate parameter for the Erlang component.
#' @param data Numeric vector of observed data.
#' @param K Integer shape parameter for the Erlang distribution.
#'
#' @details
#' Calculates the exponential rate parameter from the mean constraint and
#' evaluates the negative log-likelihood using the Erlang + Exponential PDF
#' implemented in \code{ErlangExp_Func}. Returns \code{Inf} if the exponential
#' parameter is invalid (e.g., negative or infinite) to penalize invalid parameter sets.
#'
#' @return Numeric value of the negative log-likelihood.
#'
#' @seealso [ErlangExp_Func()], [ErlangExp_Fit_v2_FixedK()]
#'
#' @examples
#' \dontrun{
#' data <- rexp(100, rate = 0.5)
#' negLL <- ErlangExp_NegLogLikelihood(3, data, 2)
#' }
#'
#' @keywords internal
ErlangExp_NegLogLikelihood <- function(ErlangLambda, data, K) {
  # Compute ExpLambda from mean constraint
  ExpLambda <- 1 / (mean(data) - (K / ErlangLambda))

  # If ExpLambda is invalid (e.g. negative or infinite), return Inf
  if (ExpLambda <= 0 || is.infinite(ExpLambda)) {
    return(Inf)
  }

  # Call the direct integral PDF
  res <- ErlangExp_Func(data, ErK = K, Erlam = ErlangLambda, Explam = ExpLambda)

  # Return the negative log-likelihood
  return(-res$Likelihood)
}


#' Erlang + Exponential Mixture PDF and Log-Likelihood
#'
#' Computes the probability density function (PDF) values of the Erlang plus
#' Exponential mixture distribution for the given data, and calculates the
#' total log-likelihood.
#'
#' @param data Numeric vector of observed data points.
#' @param ErK Integer shape parameter for the Erlang distribution.
#' @param Erlam Numeric rate parameter for the Erlang distribution.
#' @param Explam Numeric rate parameter for the Exponential distribution.
#'
#' @details
#' For each data point, numerically integrates the Erlang + Exponential PDF
#' using the given parameters. Returns a list containing the vector of PDF
#' values for each data point and the sum of their log-likelihoods.
#'
#' @return A list with components:
#' \describe{
#'   \item{Probability}{Numeric vector of PDF values for each data point.}
#'   \item{Likelihood}{Sum of log of PDF values (total log-likelihood).}
#' }
#'
#' @seealso [ErlangExp_NegLogLikelihood()], [ErlangExp_Fit_v2_FixedK()]
#'
#' @examples
#' \dontrun{
#' data <- rexp(100, rate = 0.5)
#' res <- ErlangExp_Func(data, ErK = 3, Erlam = 2, Explam = 1)
#' res$Likelihood
#' }
#'
#' @keywords internal
ErlangExp_Func <- function(data, ErK, Erlam, Explam) {
  q <- numeric(length(data))  # store PDF values

  for (i in seq_along(data)) {

    # constant multiplier
    const <- exp(-Explam * data[i]) * (Erlam^ErK) * Explam * (1 / (factorial(ErK-1))) # Can use gamma(ErK) over factorial too
    # integrand function
    integrand <- function(x) {
      exp((Explam - Erlam) * x) * (x^(ErK - 1)) * const
    }

    # numerical integration [0, data[i]]
    # integral_val <- integrate(integrand, lower = 0, upper = data[i])$value

    integral_val <- integrate(
      integrand, 0, data[i],
      rel.tol = 1e-10, abs.tol = 1e-10
    )$value




    # PDF at data[i]
    q[i] <- integral_val
  }

  Probability <- q
  Likelihood <- sum(log(q))

  return(list(Probability = Probability, Likelihood = Likelihood))
}

#' Compute CDF of Erlang + Exponential Mixture Distribution
#'
#' Calculates the cumulative distribution function (CDF) values for the Erlang
#' plus Exponential mixture distribution at specified data points using numerical integration.
#'
#' @param params Numeric vector of length 3 containing the parameters:
#' \itemize{
#'   \item Erlang shape parameter `K` (integer)
#'   \item Erlang rate parameter `lambda`
#'   \item Exponential rate parameter `lambda`
#' }
#' @param DatX Numeric vector of data points where CDF is evaluated.
#' @param interval Numeric step size for numerical integration grid (default 0.01).
#'
#' @details
#' The function computes the PDF values on a grid using `ErlangExp_Func` and
#' numerically integrates them using the trapezoidal rule to obtain the CDF.
#' The results are then interpolated at the requested points.
#'
#' @return Numeric vector of CDF values at `DatX`.
#'
#' @seealso [ErlangExp_Func()], [ErlangExp_Fit_v2_Pvalue()]
#'
#' @examples
#' \dontrun{
#' params <- c(3, 2, 1)
#' x <- seq(0, 10, by = 0.1)
#' cdf_vals <- ErlangExpCDF_Func(params, x)
#' plot(x, cdf_vals, type = "l")
#' }
#'
#' @keywords internal
ErlangExpCDF_Func <- function(params, DatX, interval = 0.01) {
  # Unpack parameters
  ErlangK <- params[1]
  ErlangLam <- params[2]
  ExpLam <- params[3]

  # Define grid from 0 to max(DatX)
  t_grid <- seq(0, max(DatX), by = interval)

  # Compute PDF values over grid using your ErlangExp_DirectInt
  pdf_vals <- ErlangExp_Func(t_grid, ErK = ErlangK, Erlam = ErlangLam, Explam = ExpLam)$Probability

  # Compute cumulative integral using trapezoidal rule
  F_cum <- cumsum(c(0, diff(t_grid) * (head(pdf_vals, -1) + tail(pdf_vals, -1)) / 2))

  # Interpolate to match DatX points
  ErlangExpCDF <- approx(x = t_grid, y = F_cum, xout = DatX, method = "linear", rule = 2)$y

  # Return CDF values
  return(ErlangExpCDF)
}


#' Generate Random Samples from Erlang + Exponential Mixture Distribution
#'
#' Simulates random variates from a distribution defined as the sum of
#' an Erlang-distributed random variable and an independent Exponential random variable.
#'
#' @param params Numeric vector of length 3 containing the parameters:
#' \itemize{
#'   \item Erlang shape parameter `K` (integer)
#'   \item Erlang rate parameter `lambda`
#'   \item Exponential rate parameter `lambda`
#' }
#' @param rows Integer, number of rows for the output sample matrix.
#' @param cols Integer, number of columns for the output sample matrix (default 1).
#'
#' @details
#' Random samples are generated by independently sampling from the Erlang and
#' Exponential components and summing them. The output is reshaped as a matrix.
#'
#' @return Numeric matrix of random samples with dimensions `rows` × `cols`.
#'
#' @seealso [ErlangExp_Func()], [ErlangExp_Fit_v2_Pvalue()]
#'
#' @examples
#' \dontrun{
#' params <- c(3, 2, 1)
#' samples <- ErlangExpRnd(params, rows = 100, cols = 1)
#' hist(samples, breaks = 30)
#' }
#'
#' @keywords internal
ErlangExpRnd <- function(params, rows, cols = 1) {
  # Unpack parameters
  ErlangK <- params[1]
  ErlangLam <- params[2]
  ExpLam <- params[3]

  n <- rows * cols

  # Generate random samples
  ErlangRnd <- rgamma(n, shape = ErlangK, scale = 1/ErlangLam)
  ExpRnd <- rexp(n, rate = ExpLam)

  # Sum component-wise
  rndarray <- ErlangRnd + ExpRnd

  # Reshape into matrix form (rows × cols)
  matrix(rndarray, nrow = rows, ncol = cols)
}
