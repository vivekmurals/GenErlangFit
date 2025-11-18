#' Fit Empirical Data Using Erlang Distribution (Version 2)
#'
#' Fits empirical data to an Erlang (Gamma with integer shape) distribution
#' using Maximum Likelihood Estimation (MLE), and evaluates goodness-of-fit
#' using bootstrap-based hypothesis testing.
#'
#' @param empiricaldata Numeric vector of observed data.
#' @param ... Optional named arguments to override defaults:
#' \describe{
#'   \item{Alpha}{Significance level for goodness-of-fit test. Default = 0.05.}
#'   \item{pvaloption}{Goodness-of-fit metric: `"KS"`, `"CvM"`, or `"AD"`. Default = `"KS"`.}
#'   \item{InitialguessK}{Initial guess for Erlang shape parameter. Default = 3.}
#'   \item{SmallestK}{Logical. If TRUE, searches for smallest acceptable K under alpha. Default = FALSE.}
#'   \item{ShowFigures}{Logical. If TRUE, generates plots. Default = TRUE.}
#'   \item{NumBootstraps}{Number of bootstrap samples for p-value estimation. Default = round(10 * 10 / Alpha).}
#' }
#'
#' @details
#' This function estimates the Erlang parameters (shape `K` and scale `lambda`)
#' by maximizing the log-likelihood function, and tests whether the data are
#' consistent with the fitted model via bootstrap-based Kolmogorov–Smirnov
#' (or other) tests.
#' Optionally, it can also search for the smallest integer `K` that still passes
#' the chosen significance threshold (`SmallestK = TRUE`).
#'
#' @return A list containing:
#' \describe{
#'   \item{Best}{A list with MLE estimates, p-value, log-likelihood, and test metrics.}
#'   \item{Smallest}{If \code{SmallestK = TRUE}, results for the smallest acceptable K.}
#' }
#'
#' @seealso [Erlang_Fit_v2_Pvalue()], [Erlang_NegLogLikelihood()]
#'
#' @examples
#' \dontrun{
#' data <- rexp(200, rate = 0.5)
#' fit <- Erlang_Fit_v2(data, Alpha = 0.05, pvaloption = "KS", ShowFigures = FALSE)
#' fit$Best$K_star
#' }
#'
#' @keywords internal
Erlang_Fit_v2 <- function(empiricaldata, ...) {

  # -----------------------------
  # 1. Default Input Arguments
  # -----------------------------

  defaultAlpha = 0.05;

  defaults <- list(
    Alpha = defaultAlpha,
    pvaloption = "KS",
    InitialguessK = 3,
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
  # 3. Override defaults with user-supplied values & Post-processing dependent defaults
  # -----------------------------
  opts <- modifyList(defaults, user_args)

  if (is.null(user_args$NumBootstraps)) {
    opts$NumBootstraps <- round(10 * 10 / opts$Alpha)
  }

  # -----------------------------
  # 4. Assign to variables (for convenience)
  # -----------------------------
  Alpha <- opts$Alpha
  pvaloption <- opts$pvaloption
  InitialguessK <- opts$InitialguessK
  SmallestK <- opts$SmallestK
  ShowFigures <- opts$ShowFigures
  NumBootstraps <- opts$NumBootstraps

  # -----------------------------
  # 5. Display configuration summary
  # -----------------------------
  message("Erlang_Fit_v2 initialized with user-parsed parameters:")
  print(opts)

  # -----------------------------
  # 6. Fitting data to Erlang via MLE
  # -----------------------------

  # -----------------------------
  # Find MLE for Gamma (continuous K)
  # -----------------------------
  data_vec = empiricaldata
  # Initial guesses
  k_init <- InitialguessK

  # Use 1D optimizer since we're only searching over k
  opt <- optimize(f = Erlang_NegLogLikelihood,
                  interval = c(1e-5, 1e5),  # adjust bounds if needed
                  tol = 1e-8,
                  data = data_vec)

  k_mle <- opt$minimum
  lambda_mle <- mean(data_vec) / k_mle
  cat("MLE shape (k):", k_mle, "\n")
  cat("MLE scale (lambda):", lambda_mle, "\n")

  # -----------------------------
  # Find MLE for Erlang (integer K)
  # -----------------------------
  mu <- mean(data_vec)
  floor_k <- floor(k_mle)
  ceil_k <- ceiling(k_mle)

  if (floor_k == 0 || ceil_k == 0) {
    # Special case: K = 0 not admissible
    k_star <- 1
    lambda_star <- mu / k_star
    log_likelihood_star <- sum(dgamma(data_vec, shape = k_star, scale = lambda_star, log = TRUE))
    cat("Special case: K=0 not allowed.\n")
    cat("Setting K* =", k_star, ", lambda* =", lambda_star, "\n")
    cat("Log-likelihood at MLE:", log_likelihood_star, "\n")
  } else {
    # Candidate integer K values
    Kint <- c(floor_k, ceil_k)

    # Compute log-likelihoods for each candidate K
    log_likelihoods <- sapply(Kint, function(K) {
      lambda <- mu / K
      sum(dgamma(data_vec, shape = K, scale = lambda, log = TRUE))
    })

    # Select K* with max log-likelihood
    idx <- which.max(log_likelihoods)
    k_star <- Kint[idx]
    lambda_star <- mu / k_star
    log_likelihood_star <- log_likelihoods[idx]

    cat("Erlang MLE K*:", k_star, "\n")
    cat("Erlang MLE Lambda*:", lambda_star, "\n")
    cat("Log-likelihood at MLE:", log_likelihood_star, "\n")
  }

  # -----------------------------
  # 7. Goodness of Fit Evaluation
  # -----------------------------


  if (tolower(pvaloption) == "nil") {
    # Skip p-value computation (used for internal multi-K searches)
    cat("Skipping p-value computation (pvaloption = NIL)\n")

  } else {

    # Using the function
    res <- Erlang_Fit_v2_Pvalue(empiricaldata,
                                k_star,
                                lambda_star,
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



  # -----------------------------
  # 8. Smallest K
  # -----------------------------

  if (SmallestK) {

    if (q_value == 0) {
      message("No smaller K found: data already fails goodness-of-fit.")
      # Store results for smallest K same as best fit results
      k_smallest <- k_star
      lambda_smallest <- lambda_star
      p_smallest <- p_value
      q_smallest <- q_value
      log_likelihood_smallest <- log_likelihood_star
      gof_metric_smallest <- gof_metric
      gof_Sample_Stats_smallest <- gof_Sample_Stats


    } else if (q_value > 0) {

      count <- 1
      k_temp <- k_star
      p_temp <- p_value
      q_temp <- q_value
      metric_star_store <- c(gof_metric)
      sampleStatsStore <- as.matrix(gof_Sample_Stats)

      # ShowFigures_backup <- ShowFigures  # temporarily disable plots
      # ShowFigures <- FALSE

      while (q_temp > 0) {
        k_temp <- k_temp - 1
        if (k_temp == 0) break
        lambda_temp <- mean(empiricaldata) / k_temp
        res_temp <- Erlang_Fit_v2_Pvalue(empiricaldata,
                                         k_temp,
                                         lambda_temp,
                                         s = length(data_vec),
                                         n = NumBootstraps,
                                         alpha = Alpha,
                                         pvaloption = pvaloption,
                                         ShowFigures)
        p_temp <- c(p_temp, res_temp$p_value)
        q_temp <- res_temp$q_value
        metric_star_store <- c(metric_star_store, res_temp$metric_star)
        sampleStatsStore <- cbind(sampleStatsStore, res_temp$sample_stats)
        count <- count + 1

      }

      # Smallest acceptable K
      k_smallest <- k_temp + 1
      lambda_smallest <- mean(empiricaldata) / k_smallest


      if (count > 1) {
        p_smallest <- p_temp[length(p_temp)-1]  # last acceptable K same as main
        gof_metric_smallest <- metric_star_store[length(metric_star_store) - 1]
        gof_Sample_Stats_smallest <- sampleStatsStore[, ncol(sampleStatsStore) - 1, drop = FALSE]


      } else {
        p_smallest <- p_temp

        gof_metric_smallest <- metric_star_store[1]
        gof_Sample_Stats_smallest <- sampleStatsStore[, 1, drop = FALSE]

      }

      q_smallest <- as.numeric(p_smallest > Alpha)
      log_likelihood_smallest <- sum(dgamma(empiricaldata, shape = k_smallest, scale = lambda_smallest, log = TRUE))

      message("Smallest K found: ", k_smallest,
              ", Lambda: ", round(lambda_smallest, 3),
              ", p-value: ", round(p_smallest, 4))
    }


  }


  # # -----------------------------
  # # 8. Plot necessary figures
  # # -----------------------------


  if (ShowFigures) {

    # Plot Empirical Data Histogram and Best Fit PDF

    bin_width <- 2 * IQR(empiricaldata) / (length(empiricaldata)^(1/3))


    x <- seq(0, max(data_vec) + 5, length.out = 500)

    df_erlang <- data.frame(
      x = x,
      y = dgamma(x, shape = k_star, scale = lambda_star),
      Type = paste0("Erlang PDF: K* = ", k_star,
                    ", λ* = ", round(lambda_star, 2))
    )


    erlang_label <- df_erlang$Type[1]

    P1 <- ggplot() +
      geom_histogram(aes(x = data_vec, y = ..density.., fill = "Observed Data"),
                     binwidth = bin_width, color = "black", alpha = 0.8) +
      geom_line(data = df_erlang, aes(x = x, y = y, color = erlang_label),
                linewidth = 1.2, linetype = "solid") +
      scale_fill_manual(name = "", values = c("Observed Data" = "#78A5A3")) +
      scale_color_manual(name = "", values = setNames(c("red"), c(erlang_label))) +
      labs(title = paste("Histogram with Erlang Fits"),
           x = "x", y = "Density") +
      theme_gray() +
      theme(
        legend.position = "bottom",        # Move legend below the plot
        legend.box = "horizontal",         # Horizontal layout
        legend.title = element_blank(),    # Remove legend title
        legend.justification = "center",   # Center the legend
        plot.title = element_text(hjust = 0.5) # Center title
      )

    if (SmallestK){

      y_small <- dgamma(x, shape = k_smallest, scale = lambda_smallest)
      df_small <- data.frame(
        x = x,
        y = y_small,
        Type = paste0("Smallest K PDF: K = ", k_smallest,
                      ", λ = ", round(lambda_smallest, 2))
      )

      # Add as a new geom_line with fixed color
      P1 <- P1 +
        geom_line(data = df_small, aes(x = x, y = y, color = Type),
                  linewidth = 1.2, linetype = "solid")
      P1 <- P1 +
        scale_color_manual(
          name = "",
          values = setNames(c("red", "grey50"),
                            c(erlang_label, df_small$Type[1]))
        )
    }


    # Base CDF data
    ecdf_data <- ecdf(data_vec)
    x_vals <- sort(data_vec)
    datY <- ecdf_data(x_vals)

    # Best-fit Erlang CDF
    gammaY <- pgamma(x_vals, shape = k_star, scale = lambda_star)

    df_erlang_cdf <- data.frame(
      x = x_vals,
      y = gammaY,
      Type = paste0("Erlang CDF: K* = ", k_star,
                    ", λ* = ", round(lambda_star, 2))
    )

    erlang_label_CDF <- df_erlang_cdf$Type[1]

    # Base plot
    P2 <- ggplot() +
      geom_step(aes(x = x_vals, y = datY, color = "Empirical CDF")) +
      geom_line(data = df_erlang_cdf, aes(x = x, y = y, color = erlang_label_CDF),
                linewidth = 1.2, linetype = "solid") +
      scale_color_manual(name = "", values = setNames(c("red", "black"),
                                                      c(erlang_label_CDF, "Empirical CDF"))) +
      labs(title = "Empirical vs Erlang CDF",
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
      y_small <- pgamma(x_vals, shape = k_smallest, scale = lambda_smallest)
      df_small <- data.frame(
        x = x_vals,
        y = y_small,
        Type = paste0("Smallest K CDF: K = ", k_smallest,
                      ", λ = ", round(lambda_smallest, 2))
      )

      P2 <- P2 +
        geom_line(data = df_small, aes(x = x, y = y, color = Type),
                  linewidth = 1.2, linetype = "solid") +
        scale_color_manual(
          name = "",
          values = setNames(c("red", "black", "grey50"),
                            c(erlang_label_CDF, "Empirical CDF", df_small$Type[1]))
        )
    }

    print(P1)
    print(P2)


    if (tolower(pvaloption) != "nil") {

      # Plot Bootstrap Statistics
      pvaloption_upper <- toupper(pvaloption)
      # Histogram of bootstrap statistics
      sampleStats = gof_Sample_Stats
      metric_star = gof_metric
      df_stats <- data.frame(Statistic = sampleStats)
      P3 <- ggplot(df_stats, aes(x = Statistic)) +
        geom_histogram(binwidth = diff(range(sampleStats)) / 30,
                       fill = "#78A5A3", alpha = 0.6, color = "black") +
        geom_vline(xintercept = metric_star, linetype = "dashed", color = "black", linewidth = 1) +
        labs(title = paste("Bootstrap Distribution of", pvaloption_upper, "Statistic"),
             x = paste(pvaloption_upper, "Statistic"), y = "Count") +
        theme_gray()

      print(P3)

      if (SmallestK){

        sampleStats = gof_Sample_Stats_smallest
        metric_star = gof_metric_smallest
        df_stats <- data.frame(Statistic = sampleStats)
        P4 <- ggplot(df_stats, aes(x = Statistic)) +
          geom_histogram(binwidth = diff(range(sampleStats)) / 30,
                         fill = "#78A5A3", alpha = 0.6, color = "black") +
          geom_vline(xintercept = metric_star, linetype = "dashed", color = "black", linewidth = 1) +
          labs(title = paste("Bootstrap Distribution of", pvaloption_upper, "Statistic (Smallest K)"),
               x = paste(pvaloption_upper, "Statistic"), y = "Count") +
          theme_minimal()        # Bootstrap for Smallest K

        print(P4)

      }


    }


  }



  return(list(
    Best = if (tolower(pvaloption) != "nil") {
      list(
        K_star = k_star,
        Lambda_star = lambda_star,
        P_star = p_value,
        Q_Value = q_value,
        Loglikelihood = log_likelihood_star,
        metric_star = gof_metric,
        samplestats_star = gof_Sample_Stats
      )
    } else {
      list(
        K_star = k_star,
        Lambda_star = lambda_star,
        Loglikelihood = log_likelihood_star
      )
    },

    Smallest = if (SmallestK) list(
      K_star = k_smallest,
      Lambda_star = lambda_smallest,
      P_star = p_smallest,
      Q_Value = q_smallest,
      Loglikelihood = log_likelihood_smallest,
      metric_star = gof_metric_smallest,
      samplestats_star = gof_Sample_Stats_smallest
    ) else NULL
  ))

}



#' Bootstrap-Based Goodness-of-Fit for Erlang Distribution
#'
#' Computes a bootstrap-based p-value for the Erlang model fit
#' using Kolmogorov–Smirnov, Cramér–von Mises, or Anderson–Darling statistics.
#'
#' @param empiricaldata Numeric vector of observed data.
#' @param k_star Numeric, fitted Erlang shape parameter.
#' @param lambda_star Numeric, fitted Erlang scale parameter.
#' @param s Integer, number of observations. Defaults to \code{length(empiricaldata)}.
#' @param n Integer, number of bootstrap samples. Default = 1000.
#' @param alpha Numeric, significance level. Default = 0.05.
#' @param pvaloption Character. Choice of test statistic: `"KS"`, `"CvM"`, or `"AD"`.
#' @param ShowFigures Logical. If TRUE, produces optional plots (currently suppressed).
#'
#' @details
#' The function compares the empirical CDF of the observed data to the
#' theoretical Erlang CDF and computes a test statistic. Bootstrap samples
#' from the fitted Erlang model are used to estimate the null distribution
#' and compute the corresponding p-value.
#'
#' @return A list containing:
#' \describe{
#'   \item{p_value}{Bootstrap-estimated p-value.}
#'   \item{q_value}{Binary indicator (1 = fail to reject, 0 = reject).}
#'   \item{metric_star}{Observed test statistic.}
#'   \item{sample_stats}{Bootstrap sample statistics.}
#' }
#'
#' @seealso [Erlang_Fit_v2()], [Erlang_NegLogLikelihood()]
#'
#' @examples
#' \dontrun{
#' data <- rexp(200, rate = 0.5)
#' fit <- Erlang_Fit_v2(data, ShowFigures = FALSE)
#' pvals <- Erlang_Fit_v2_Pvalue(data, fit$Best$K_star, fit$Best$Lambda_star)
#' }
#'
#' @keywords internal
Erlang_Fit_v2_Pvalue <- function(empiricaldata, k_star, lambda_star,
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
  gammaY <- pgamma(x_vals, shape = k_star, scale = lambda_star)

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
    sample_data <- rgamma(s, shape = k_star, scale = lambda_star)
    ecdf_sample <- ecdf(sample_data)
    x_sample <- sort(sample_data)
    sampleY <- ecdf_sample(x_sample)
    gammaY_sample <- pgamma(x_sample, shape = k_star, scale = lambda_star)


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

  # # -----------------------------
  # # 6. Figures
  # # -----------------------------
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




#' Negative Log-Likelihood for Erlang Distribution
#'
#' Computes the negative log-likelihood (to be minimized) for a given
#' Erlang (Gamma) shape parameter `k` and data vector. The scale parameter
#' is internally estimated as \eqn{lambda = mean(data) / k}.
#'
#' @param k Numeric, Erlang shape parameter.
#' @param data Numeric vector of observed data.
#'
#' @return Numeric scalar representing the negative log-likelihood.
#'
#' @seealso [Erlang_Fit_v2()], [Erlang_Fit_v2_Pvalue()]
#'
#' @examples
#' \dontrun{
#' data <- rexp(200, rate = 0.5)
#' Erlang_NegLogLikelihood(3, data)
#' }
#'
#' @keywords internal
Erlang_NegLogLikelihood <- function(k, data) {
  if (k <= 0) return(Inf)
  lambda <- mean(data) / k
  -sum(dgamma(data, shape = k, scale = lambda, log = TRUE))
}
