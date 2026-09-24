library(ggplot2)
library(goftest)

# ---- Fixed reference dataset ----
set.seed(123)
ref_data <- rgamma(200, shape = 5, rate = 0.8)

# ---- Custom CDF wrapper for Erlang-Exp (needed for base R/goftest cross-checks) ----
pErlangExp <- function(q, ErlangK, ErlangLam, ExpLam) {
  ErlangExpCDF_Func(c(ErlangK, ErlangLam, ExpLam), q, interval = 0.01)
}


# ===========================================================
# ERLANG — 4 tests: LogL, KS, CVM, AD
# ===========================================================

test_that("[Erlang] Log-likelihood matches manual dgamma computation", {
  fit <- GenErlang_Fit('Erlang', ref_data, pvaloption = 'KS')
  k_hat <- fit$Best$K_star; lam_hat <- fit$Best$Lambda_star

  manual_ll <- sum(dgamma(ref_data, shape = k_hat, rate = lam_hat, log = TRUE))

  expect_equal(fit$Best$Loglikelihood, manual_ll, tolerance = 1e-4)
})

test_that("[Erlang] KS statistic matches base R", {
  fit <- GenErlang_Fit('Erlang', ref_data, pvaloption = 'KS')
  k_hat <- fit$Best$K_star; lam_hat <- fit$Best$Lambda_star

  baseR_ks <- suppressWarnings(ks.test(ref_data, "pgamma", shape = k_hat, rate = lam_hat))

  expect_equal(as.numeric(fit$Best$metric_star), as.numeric(baseR_ks$statistic), tolerance = 1e-6)
})

test_that("[Erlang] CvM statistic matches goftest implementation", {
  fit <- GenErlang_Fit('Erlang', ref_data, pvaloption = 'CVM')
  k_hat <- fit$Best$K_star; lam_hat <- fit$Best$Lambda_star

  baseR_cvm <- goftest::cvm.test(ref_data, "pgamma", shape = k_hat, rate = lam_hat)

  expect_equal(as.numeric(fit$Best$metric_star), as.numeric(baseR_cvm$statistic), tolerance = 1e-6)
})

test_that("[Erlang] AD statistic matches goftest implementation", {
  fit <- GenErlang_Fit('Erlang', ref_data, pvaloption = 'AD')
  k_hat <- fit$Best$K_star; lam_hat <- fit$Best$Lambda_star

  baseR_ad <- goftest::ad.test(ref_data, "pgamma", shape = k_hat, rate = lam_hat)

  expect_equal(as.numeric(fit$Best$metric_star), as.numeric(baseR_ad$statistic), tolerance = 1e-6)
})


# ===========================================================
# FOUNDATIONAL CHECK — validates ErlangExpCDF_Func against an
# independent oracle (base R simulation), so the 4 tests below
# that rely on it are not circular.
# ===========================================================

test_that("[ErlangExp] Helper CDF function matches independent Monte Carlo simulation", {
  set.seed(456)
  K_true <- 3; erlam_true <- 0.5; explam_true <- 1.2

  # Independent simulation: sum of Erlang(K, rate) and Exp(rate),
  # built only from base R primitives (rgamma + rexp) — no dependence
  # on ErlangExp_Func, ErlangExpCDF_Func, or GenErlang_Fit at all.
  sim_data <- rgamma(5000, shape = K_true, rate = erlam_true) + rexp(5000, rate = explam_true)

  # Compare the simulated sample against the package helper's CDF
  # via a one-sample KS test.
  ks_result <- suppressWarnings(
    ks.test(sim_data, pErlangExp, ErlangK = K_true, ErlangLam = erlam_true, ExpLam = explam_true)
  )

  # If ErlangExpCDF_Func correctly represents the Erlang + Exponential
  # convolution, the simulated data should be statistically consistent
  # with it (i.e., no significant mismatch).
  expect_gt(ks_result$p.value, 0.01)
})


# ===========================================================
# ERLANG-EXP — 4 tests: LogL, KS, CVM, AD
# ===========================================================

test_that("[ErlangExp] Log-likelihood matches manual density computation", {
  fit <- GenErlang_Fit('ErlangExp', ref_data, 3, FixedK = TRUE, pvaloption = 'KS')
  k_hat <- fit$Best$K_star
  erlam <- fit$Best$ErlangLambda_star
  explam <- fit$Best$ExpLambda_star

  pdf_vals <- ErlangExp_Func(ref_data, ErK = k_hat, Erlam = erlam, Explam = explam)$Probability
  manual_ll <- sum(log(pmax(pdf_vals, 1e-300)))

  # NOTE: field name is "LogLikelihood" (capital L) for ErlangExp,
  # vs "Loglikelihood" (lowercase l) for Erlang. This is a naming
  # inconsistency across model branches in the package -- flagged
  # separately in CODE_REVIEW.md. Confirmed via testing that this
  # field, unlike the Erlang one, already returns a value directly
  # comparable to manual_ll with no sign flip needed.
  expect_equal(fit$Best$LogLikelihood, manual_ll, tolerance = 1e-3)
})

test_that("[ErlangExp] KS statistic matches base R via custom CDF wrapper", {
  fit <- GenErlang_Fit('ErlangExp', ref_data, 3, FixedK = TRUE, pvaloption = 'KS')
  k_hat <- fit$Best$K_star
  erlam <- fit$Best$ErlangLambda_star
  explam <- fit$Best$ExpLambda_star

  baseR_ks <- suppressWarnings(
    ks.test(ref_data, pErlangExp, ErlangK = k_hat, ErlangLam = erlam, ExpLam = explam)
  )

  expect_equal(as.numeric(fit$Best$metric_star), as.numeric(baseR_ks$statistic), tolerance = 1e-4)
})

test_that("[ErlangExp] CvM statistic matches goftest implementation", {
  fit <- GenErlang_Fit('ErlangExp', ref_data, 3, FixedK = TRUE, pvaloption = 'CVM')
  k_hat <- fit$Best$K_star
  erlam <- fit$Best$ErlangLambda_star
  explam <- fit$Best$ExpLambda_star

  baseR_cvm <- goftest::cvm.test(ref_data, pErlangExp,
                                 ErlangK = k_hat, ErlangLam = erlam, ExpLam = explam)

  expect_equal(as.numeric(fit$Best$metric_star), as.numeric(baseR_cvm$statistic), tolerance = 1e-4)
})

test_that("[ErlangExp] AD statistic matches goftest implementation", {
  fit <- GenErlang_Fit('ErlangExp', ref_data, 3, FixedK = TRUE, pvaloption = 'AD')
  k_hat <- fit$Best$K_star
  erlam <- fit$Best$ErlangLambda_star
  explam <- fit$Best$ExpLambda_star

  baseR_ad <- goftest::ad.test(ref_data, pErlangExp,
                               ErlangK = k_hat, ErlangLam = erlam, ExpLam = explam)

  expect_equal(as.numeric(fit$Best$metric_star), as.numeric(baseR_ad$statistic), tolerance = 1e-4)
})
