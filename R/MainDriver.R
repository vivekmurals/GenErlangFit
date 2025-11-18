#' Generalized Erlang and Erlang-Exponential Distribution Fitting
#'
#' Fits empirical data to Erlang or Erlang-Exponential mixture models using maximum likelihood estimation.
#' This wrapper function calls the appropriate specialized fitting routines based on the \code{mode} parameter.
#'
#' @param mode Character string specifying the model to fit. Options are \code{"Erlang"} or \code{"ErlangExp"}.
#'   Alternatively, if \code{mode} is numeric with no additional arguments, it is interpreted as empirical data, and both models are fit.
#' @param empiricaldata Numeric vector of empirical data to be fitted. Required when \code{mode} is \code{"Erlang"} or \code{"ErlangExp"}.
#' @param K Integer specifying the initial guess for the shape parameter \code{K} (number of Erlang phases) for \code{"ErlangExp"} mode. Required if \code{mode = "ErlangExp"}.
#' @param ... Additional optional parameters passed down to the underlying fitting functions. For example:
#'   \itemize{
#'     \item \code{SmallestK} Logical. Whether to search for the smallest acceptable \code{K} satisfying goodness-of-fit.
#'     \item Other model-specific options.
#'   }
#'
#' @return A list containing fitting results:
#' \itemize{
#'   \item For \code{"Erlang"} mode, returns the result list from \code{\link{Erlang_Fit_v2}}.
#'   \item For \code{"ErlangExp"} mode, returns the result list from \code{\link{ErlangExp_Fit_v2}}.
#'   \item If only empirical data is provided (numeric \code{mode}), returns a list with both Erlang and Erlang-Exp fit results.
#' }
#'
#' @examples
#' \dontrun{
#' data <- rexp(100, rate = 0.1)
#' GenErlang_Fit("Erlang", data, SmallestK = TRUE)
#' GenErlang_Fit("ErlangExp", data, K = 3)
#' GenErlang_Fit(data) # fits both models and returns results
#' }
#'
#' @seealso \code{\link{Erlang_Fit_v2}}, \code{\link{ErlangExp_Fit_v2}}
#'
#' @export
GenErlang_Fit <- function(mode, empiricaldata = NULL, K = NULL, ...) {

  # Parse additional arguments
  user_args <- list(...)

  # ---------------------------
  # 1. Erlang fit
  # ---------------------------
  if (is.character(mode) && length(mode) == 1) {

    if (tolower(mode) == "erlang") {

      if (isTRUE(user_args$SmallestK)) {
        Erlang_Results <- Erlang_Fit_v2(empiricaldata, SmallestK = TRUE, ...)
        ModelNames <- c("Erlang", "Erlang Smallest K")
        K_vals <- c(Erlang_Results$Best$K_star, Erlang_Results$Smallest$K_star)
        Lambda <- c(Erlang_Results$Best$Lambda_star, Erlang_Results$Smallest$Lambda_star)
        LogL <- c(Erlang_Results$Best$Loglikelihood, Erlang_Results$Smallest$Loglikelihood)
      } else {
        Erlang_Results <- Erlang_Fit_v2(empiricaldata, ...)
        ModelNames <- "Erlang"
        K_vals <- Erlang_Results$Best$K_star
        Lambda <- Erlang_Results$Best$Lambda_star
        LogL <- Erlang_Results$Best$Loglikelihood
      }

      ResultsTable <- data.frame(
        Model = ModelNames,
        K = K_vals,
        ErlangLambda = Lambda,
        LogLikelihood = LogL,
        row.names = NULL
      )
      print(ResultsTable, row.names = FALSE)
      return(Erlang_Results)

    } else if (tolower(mode) == "erlangexp") {

      if (is.null(K)) stop("Provide initial K for Erlang-Exp fit")

      if (isTRUE(user_args$SmallestK)) {
        ErlangExp_Results <- ErlangExp_Fit_v2(empiricaldata, K, SmallestK = TRUE, ...)

        ModelNames <- c("Erlang–Exp", "Erlang–Exp (Smallest K)")
        K_vals <- c(ErlangExp_Results$Best$K_star, ErlangExp_Results$Smallest$K_star)
        ErlangLambda <- c(ErlangExp_Results$Best$ErlangLambda_star, ErlangExp_Results$Smallest$ErlangLambda_star)
        ExpLambda <- c(ErlangExp_Results$Best$ExpLambda_star, ErlangExp_Results$Smallest$ExpLambda_star)
        LogL <- c(ErlangExp_Results$Best$LogLikelihood, ErlangExp_Results$Smallest$LogLikelihood)

      } else {
        ErlangExp_Results <- ErlangExp_Fit_v2(empiricaldata, K, ...)

        ModelNames <- "Erlang–Exp"
        K_vals <- ErlangExp_Results$Best$K_star
        ErlangLambda <- ErlangExp_Results$Best$ErlangLambda_star
        ExpLambda <- ErlangExp_Results$Best$ExpLambda_star
        LogL <- ErlangExp_Results$Best$LogLikelihood
      }

      ResultsTable <- data.frame(
        Model = ModelNames,
        K = K_vals,
        ErlangLambda = ErlangLambda,
        ExpLambda = ExpLambda,
        LogLikelihood = LogL,
        row.names = NULL
      )
      print(ResultsTable, row.names = FALSE)
      return(ErlangExp_Results)

    } else {
      stop('Unknown mode. Use "Erlang" or "ErlangExp".')
    }

    # ---------------------------
    # 2. Data-only call
    # ---------------------------
  } else if (is.numeric(mode) && length(list(...)) == 0) {

    data <- mode  # `mode` actually holds empiricaldata

    # Run both Erlang and Erlang–Exp fits
    Erlang_Results <- Erlang_Fit_v2(data, SmallestK = TRUE)
    ErlangExp_Results <- ErlangExp_Fit_v2(
      data,
      Erlang_Results$Best$K_star - 1,
      SmallestK = TRUE,
      SmallestKValue = Erlang_Results$Smallest$K_star - 1
    )

    # Store results
    GenErlang_Results <- list(
      Erlang_Results = Erlang_Results,
      ErlangExp_Results = ErlangExp_Results
    )

    # Build summary table
    Model <- c(
      "Erlang",
      "Erlang Smallest K",
      "Erlang–Exp (Local Best)",
      "Erlang–Exp Smallest K"
    )

    K <- c(
      Erlang_Results$Best$K_star,
      Erlang_Results$Smallest$K_star,
      ErlangExp_Results$Best$K_star,
      ErlangExp_Results$Smallest$K_star
    )

    ErlangLambda <- c(
      Erlang_Results$Best$Lambda_star,
      Erlang_Results$Smallest$Lambda_star,
      ErlangExp_Results$Best$ErlangLambda_star,
      ErlangExp_Results$Smallest$ErlangLambda_star
    )

    ExpLambda <- c(
      0,
      0,
      ErlangExp_Results$Best$ExpLambda_star,
      ErlangExp_Results$Smallest$ExpLambda_star
    )

    LogLikelihood <- c(
      Erlang_Results$Best$LogLikelihood,
      Erlang_Results$Smallest$LogLikelihood,
      ErlangExp_Results$Best$LogLikelihood,
      ErlangExp_Results$Smallest$LogLikelihood
    )

    ResultsTable <- data.frame(
      Model,
      K,
      ErlangLambda,
      ExpLambda,
      LogLikelihood,
      row.names = NULL
    )


    print(ResultsTable, row.names = FALSE)

    return(GenErlang_Results)

    # ---------------------------
    # 3. Invalid input
    # ---------------------------
  } else {
    stop('Unknown mode. Use "Erlang", "ErlangExp", or provide only empirical data.')
  }

}

