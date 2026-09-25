# tests/testthat/test-model-selection-logic.R
#
# Tests for model-selection / orchestration logic (SmallestK selection,
# window/peak logic, quick-fit-all-modes wiring), as distinct from the
# core statistical engine validated in test-core-model-fitting.R.


# ---- Fixed reference dataset ----
set.seed(123)
ref_data <- rgamma(200, shape = 5, rate = 0.8)

# ---- Custom CDF wrapper for Erlang-Exp (needed for base R/goftest cross-checks) ----
pErlangExp <- function(q, ErlangK, ErlangLam, ExpLam) {
  ErlangExpCDF_Func(c(ErlangK, ErlangLam, ExpLam), q, interval = 0.01)
}



# ===========================================================
# ERLANG-Smallest K — Smallest K Selection Logic
# ===========================================================

test_that("[Erlang] SmallestK selection logic is internally consistent", {
  set.seed(123)
  fit <- GenErlang_Fit('Erlang', ref_data, SmallestK = TRUE, pvaloption = 'KS',
                       NumBootstraps = 200, ShowFigures = FALSE)

  best <- fit$Best
  smallest <- fit$Smallest

  expect_false(is.null(smallest))
  expect_true(smallest$K_star <= best$K_star)
  expect_equal(smallest$Q_Value, 1)
  expect_true(smallest$P_star > 0.05)

  manual_ll_smallest <- sum(dgamma(ref_data, shape = smallest$K_star,
                                   rate = smallest$Lambda_star, log = TRUE))
  expect_equal(smallest$Loglikelihood, manual_ll_smallest, tolerance = 1e-4)

  if (smallest$K_star < best$K_star) {
    k_below <- smallest$K_star - 1
    lambda_below_scale <- mean(ref_data) / k_below

    set.seed(123)
    res_below <- Erlang_Fit_v2_Pvalue(
      ref_data,
      k_below,
      lambda_below_scale,
      s = length(ref_data),
      n = 200,
      alpha = 0.05,
      pvaloption = "KS",
      ShowFigures = FALSE
    )

    expect_equal(res_below$q_value, 0)
  } else {
    succeed("Smallest K equals Best K; no decrementing occurred, boundary check not applicable.")
  }
})



# ===========================================================
# ERLANG-Exp Smallest K — Smallest K Selection Logic
# ===========================================================

test_that("[ErlangExp] SmallestK selection logic is internally consistent", {
  set.seed(123)
  fit <- GenErlang_Fit('erlangexp', ref_data, K = 4, SmallestK = TRUE,
                       pvaloption = 'KS', NumBootstraps = 200, ShowFigures = FALSE)

  best <- fit$Best
  smallest <- fit$Smallest

  expect_false(is.null(smallest))
  expect_true(smallest$K_star <= best$K_star)
  expect_equal(smallest$Q_Value, 1)
  expect_true(smallest$P_star > 0.05)

  # Manual log-likelihood cross-check via ErlangExp_Func
  manual_ll_smallest <- ErlangExp_Func(ref_data,
                                       ErK = smallest$K_star,
                                       Erlam = smallest$ErlangLambda_star,
                                       Explam = smallest$ExpLambda_star)$Likelihood
  expect_equal(smallest$LogLikelihood, manual_ll_smallest, tolerance = 1e-4)

  # Boundary check: confirm K_star - 1 independently fails GoF,
  # only if that value is a valid K (>= 1) and decrementing actually occurred.
  if (smallest$K_star < best$K_star && smallest$K_star > 1) {
    k_below <- smallest$K_star - 1

    set.seed(123)
    res_below <- ErlangExp_Fit_v2_FixedK(
      ref_data,
      k_below,
      pvaloption = "KS",
      NumBootstraps = 200,
      ShowFigures = FALSE
    )

    expect_equal(res_below$Q_Value, 0)
  } else {
    succeed("No valid K below Smallest K to test (either no decrement occurred, or Smallest K == 1).")
  }
})


# ===========================================================
# ERLANG-Exp Windowed K — Peak Selection Logic
# ===========================================================

test_that("[ErlangExp] Window K peak-finding correctly identifies local maximum in log-likelihood", {
  fit <- GenErlang_Fit('erlangexp', ref_data, 3, FixedK = FALSE, KWindowSize = 10,
                       pvaloption = 'NIL', ShowFigures = FALSE)

  peak <- fit$Best
  peak_K <- peak$K_star

  # Independently recompute LL at peak, peak-1, and peak+1 via direct FixedK calls
  ll_at <- function(k) {
    res <- ErlangExp_Fit_v2_FixedK(ref_data, k, pvaloption = "NIL", ShowFigures = FALSE)
    res$LogLikelihood
  }

  ll_peak  <- ll_at(peak_K)
  ll_below <- if (peak_K > 1) ll_at(peak_K - 1) else -Inf
  ll_above <- ll_at(peak_K + 1)

  # Confirm reported LL matches independent recomputation
  expect_equal(ll_peak, peak$LogLikelihood, tolerance = 1e-4)

  # Confirm it's actually a local maximum: neighbors should not exceed it
  expect_true(ll_peak >= ll_below)
  expect_true(ll_peak >= ll_above)
})


# ===========================================================
# GenErlangFit Quick Fit all Models Logic
# ===========================================================
