#######################################################################
# R Code for Real Data Application: Bladder Cancer Data
# Progressive Type-II censoring with three removal schemes:
#   - TP (Truncated Poisson)
#   - TDW (Truncated Discrete Weibull)
#   - BIN (Binomial with p=0.5)
# Effective sample size: m=32 (out of n=128)
# 
# This code produces:
#   - MLE for complete data (baseline)
#   - MLE and Bayes estimates for each censoring scheme
#   - Point estimates and interval estimates (CI for MLE, HPD for Bayes)
#   - Comparison table similar to Table 5 in Pathak et al. (2020)
#   - FULL display of observed X and removal R vectors for each scheme
#######################################################################

# --------------------- Load required packages ---------------------
if (!requireNamespace("maxLik", quietly = TRUE)) install.packages("maxLik")
if (!requireNamespace("numDeriv", quietly = TRUE)) install.packages("numDeriv")
if (!requireNamespace("coda", quietly = TRUE)) install.packages("coda")

library(maxLik)
library(numDeriv)
library(coda)

set.seed(2026)

# --------------------- Load and prepare data ---------------------
# Bladder cancer remission times (n=128) from Lee & Wang (2003)
bladder_data <- c(
  4.50, 32.15, 3.88, 13.80, 19.13, 4.87, 5.85, 14.24, 5.71, 7.09,
  7.87, 7.59, 20.28, 5.32, 5.49, 3.02, 46.12, 2.02, 4.51, 5.17,
  2.83, 9.22, 1.05, 0.20, 8.37, 3.82, 9.47, 36.66, 14.77, 26.31,
  79.05, 10.06, 8.53, 2.02, 4.98, 11.98, 2.62, 4.26, 5.06, 1.76,
  0.90, 11.25, 16.62, 4.40, 21.73, 10.34, 12.07, 34.26, 10.66, 6.97,
  2.07, 0.51, 12.03, 0.08, 17.12, 3.36, 2.64, 1.40, 12.63, 43.01,
  14.76, 2.75, 7.66, 0.81, 1.19, 7.32, 4.18, 3.36, 8.66, 1.26,
  13.29, 1.46, 14.83, 6.76, 23.63, 5.62, 3.25, 18.10, 7.62, 7.63,
  17.14, 25.74, 3.52, 2.87, 15.96, 17.36, 9.74, 3.31, 7.28, 1.35,
  0.40, 2.26, 4.33, 9.02, 5.41, 2.69, 22.69, 6.94, 2.54, 11.79,
  2.46, 7.26, 2.69, 5.34, 3.48, 8.26, 6.93, 4.23, 3.70, 0.50,
  10.75, 6.54, 3.64, 5.32, 13.11, 8.65, 3.57, 5.09, 7.39, 5.41,
  11.64, 2.09, 2.23, 6.25, 7.93, 4.34, 25.82, 12.02
)

n <- length(bladder_data)  # 128
m <- 32                    # effective sample size (25% of n)
p_bin <- 0.5               # Binomial removal probability

# Sort data (for progressive censoring)
X_sorted <- sort(bladder_data)

# --------------------- Weibull-Poisson distribution functions ---------------------

# Stable log-density
log_expm1_pos <- function(x) {
  ifelse(x < log(2), log(expm1(x)), x + log1p(-exp(-x)))
}

log_dwp <- function(x, a, t, l) {
  if (any(x <= 0) || any(c(a, t, l) <= 0)) return(-Inf)
  z <- t * x^a
  log(a) + log(t) + log(l) + (a - 1) * log(x) - z + l * exp(-z) - log_expm1_pos(l)
}

log_swp <- function(x, a, t, l) {
  if (any(x <= 0) || any(c(a, t, l) <= 0)) return(0)
  z <- t * x^a
  log_expm1_pos(l * exp(-z)) - log_expm1_pos(l)
}

# --------------------- Removal schemes (conditional on observed data) ---------------------

remaining_units <- function(R, i, n, m) {
  n - m - if (i == 1) 0 else sum(R[1:(i - 1)])
}

# Generate TP removal vector
sample_R_TP <- function(n, m, l, data_sorted) {
  R <- integer(m)
  rem <- n - m
  for (i in seq_len(m - 1)) {
    avail <- remaining_units(R, i, n, m)
    if (avail <= 0) break
    w <- dpois(0:avail, l)
    w <- w / sum(w)
    R[i] <- sample(0:avail, 1, prob = w)
    rem <- rem - R[i]
  }
  R[m] <- rem
  R
}

# Generate TDW removal vector
sample_R_TDW <- function(n, m, a, t, data_sorted) {
  R <- integer(m)
  rem <- n - m
  for (i in seq_len(m - 1)) {
    avail <- remaining_units(R, i, n, m)
    if (avail <= 0) break
    k <- 0:avail
    w <- exp(-t * k^a) - exp(-t * (k + 1)^a)
    w <- pmax(w, 0)
    w <- w / sum(w)
    R[i] <- sample(k, 1, prob = w)
    rem <- rem - R[i]
  }
  R[m] <- rem
  R
}

# Generate BIN removal vector
sample_R_BIN <- function(n, m, p, data_sorted) {
  R <- integer(m)
  rem <- n - m
  for (i in seq_len(m - 1)) {
    avail <- remaining_units(R, i, n, m)
    if (avail <= 0) break
    R[i] <- rbinom(1, avail, p)
    R[i] <- min(R[i], rem)
    rem <- rem - R[i]
  }
  R[m] <- rem
  R
}

# --------------------- Progressive Type-II censoring from real data ---------------------

# Given a sorted full sample and a removal vector, produce the observed sample
apply_removal_to_data <- function(data_sorted, R, n, m) {
  observed <- numeric(m)
  remaining_indices <- 1:n
  for (i in 1:m) {
    # Observe the smallest remaining
    obs_idx <- remaining_indices[1]
    observed[i] <- data_sorted[obs_idx]
    # Remove the observed unit
    remaining_indices <- remaining_indices[-1]
    # Remove R[i] units from remaining (if any)
    if (R[i] > 0 && length(remaining_indices) > 0) {
      remove_idx <- sample(1:length(remaining_indices), 
                           size = min(R[i], length(remaining_indices)), 
                           replace = FALSE)
      remaining_indices <- remaining_indices[-remove_idx]
    }
    if (length(remaining_indices) == 0) break
  }
  return(observed)
}

# Generate a censored sample with a specific scheme
generate_censored_sample <- function(scheme, data_sorted, n, m, true_par, p_bin = 0.5) {
  if (scheme == "TP") {
    R <- sample_R_TP(n, m, true_par[3], data_sorted)
  } else if (scheme == "TDW") {
    R <- sample_R_TDW(n, m, true_par[1], true_par[2], data_sorted)
  } else if (scheme == "BIN") {
    R <- sample_R_BIN(n, m, p_bin, data_sorted)
  } else {
    stop("Unknown scheme")
  }
  x_obs <- apply_removal_to_data(data_sorted, R, n, m)
  return(list(X = x_obs, R = R))
}

# --------------------- Log-likelihood functions ---------------------

loglik_L1 <- function(par, x, R) {
  a <- par[1]; t <- par[2]; l <- par[3]
  if (any(par <= 0)) return(-Inf)
  sum(log_dwp(x, a, t, l)) + sum(R * log_swp(x, a, t, l))
}

loglik_L2_TP <- function(par, R, n, m) {
  l <- par[3]
  if (l <= 0) return(-Inf)
  ans <- 0
  for (i in 1:(m - 1)) {
    avail <- remaining_units(R, i, n, m)
    vals <- 0:avail
    den <- sum(dpois(vals, l))
    ans <- ans + log(dpois(R[i], l) / den)
  }
  ans
}

loglik_L2_TDW <- function(par, R, n, m) {
  a <- par[1]; t <- par[2]
  if (a <= 0 || t <= 0) return(-Inf)
  ans <- 0
  for (i in 1:(m - 1)) {
    avail <- remaining_units(R, i, n, m)
    ri <- R[i]
    num <- exp(-t * ri^a) - exp(-t * (ri + 1)^a)
    den <- 1 - exp(-t * (avail + 1)^a)
    if (num <= 0 || den <= 0) return(-Inf)
    ans <- ans + log(num) - log(den)
  }
  ans
}

# Full log-likelihood
loglik_wp <- function(par, x, R, n, m, scheme) {
  ans <- loglik_L1(par, x, R)
  if (!is.finite(ans)) return(-Inf)
  if (scheme == "TP") ans <- ans + loglik_L2_TP(par, R, n, m)
  else if (scheme == "TDW") ans <- ans + loglik_L2_TDW(par, R, n, m)
  # For BIN, L2 is constant (0)
  ans
}

# --------------------- MLE (on log scale) ---------------------

mle_estim <- function(x, R, n, m, scheme, start_par) {
  obj <- function(lp) {
    par <- exp(lp)
    if (any(par > 50)) return(1e30)
    val <- loglik_wp(par, x, R, n, m, scheme)
    ifelse(is.finite(val), -val, 1e30)
  }
  
  starts <- list(log(start_par),
                 log(c(1, 0.01, 4)),  # plausible for bladder data
                 log(c(0.5, 0.005, 2)),
                 log(c(2, 0.02, 6)))
  
  best <- NULL; best_val <- Inf
  for (st in starts) {
    opt <- try(optim(st, obj, method = "Nelder-Mead",
                     control = list(maxit = 3000, reltol = 1e-8)), silent = TRUE)
    if (!inherits(opt, "try-error") && opt$convergence %in% c(0,1) && opt$value < best_val) {
      best <- opt; best_val <- opt$value
    }
  }
  if (is.null(best)) return(NULL)
  est <- exp(best$par)
  names(est) <- c("alpha", "theta", "lambda")
  est
}

# Confidence intervals for MLE (using Hessian)
mle_ci <- function(x, R, n, m, scheme, est) {
  # Compute Hessian on log scale
  H <- tryCatch(numDeriv::hessian(function(lp) {
    par <- exp(lp)
    -loglik_wp(par, x, R, n, m, scheme)
  }, log(est)), error = function(e) NULL)
  
  if (is.null(H) || !all(is.finite(H))) return(matrix(NA, 3, 2))
  V <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V)) return(matrix(NA, 3, 2))
  se <- est * sqrt(diag(V))
  ci <- cbind(lower = est - qnorm(0.975) * se,
              upper = est + qnorm(0.975) * se)
  rownames(ci) <- c("alpha", "theta", "lambda")
  return(ci)
}

# --------------------- Bayesian estimation (MCMC) ---------------------

logprior <- function(par, type) {
  if (any(par <= 0)) return(-Inf)
  if (type == "Gamma22") {
    sum(dgamma(par, shape = 2, rate = 2, log = TRUE))
  } else {  # Jeffreys
    -sum(log(par))
  }
}

bayes_estim <- function(x, R, n, m, scheme, init_par, prior_type = "Gamma22",
                        n_iter = 15000, burn = 5000, thin = 2,
                        proposal_sd = c(0.1, 0.1, 0.1)) {
  eta <- log(init_par)
  target <- function(e) {
    par <- exp(e)
    loglik_wp(par, x, R, n, m, scheme) + logprior(par, prior_type) + sum(e)
  }
  cur <- target(eta)
  out <- matrix(NA, n_iter, 3)
  ac <- 0
  for (i in seq_len(n_iter)) {
    pr <- eta + rnorm(3, 0, proposal_sd)
    v <- target(pr)
    if (is.finite(v) && log(runif(1)) < v - cur) {
      eta <- pr; cur <- v; ac <- ac + 1
    }
    out[i, ] <- exp(eta)
  }
  keep <- out[seq(burn + 1, n_iter, by = thin), , drop = FALSE]
  colnames(keep) <- c("alpha", "theta", "lambda")
  est <- colMeans(keep)
  hpd <- t(apply(keep, 2, function(v) as.numeric(HPDinterval(mcmc(v), 0.95))))
  colnames(hpd) <- c("lower", "upper")
  list(est = est, hpd = hpd, acceptance = ac / n_iter)
}

# --------------------- Complete data MLE (baseline) ---------------------

# MLE for complete data (using maxLik package)
logL_complete <- function(par) {
  a <- par[1]; t <- par[2]; l <- par[3]
  if (any(par <= 0)) return(-Inf)
  sum(log_dwp(bladder_data, a, t, l))
}

mle_complete <- maxLik(logLik = logL_complete, start = c(1, 0.01, 4), method = "BFGS")
true_par <- coef(mle_complete)
names(true_par) <- c("alpha", "theta", "lambda")

cat("\n========== COMPLETE DATA MLE ==========\n")
cat(sprintf("alpha = %.5f, theta = %.5f, lambda = %.5f\n",
            true_par[1], true_par[2], true_par[3]))

# --------------------- Generate censored samples for each scheme ---------------------

set.seed(2026)

schemes <- c("TP", "TDW", "BIN")
results <- list()

for (sch in schemes) {
  cat("\n\n========== Scheme:", sch, "==========\n")
  
  # Generate censored sample
  censored <- generate_censored_sample(sch, X_sorted, n, m, true_par, p_bin)
  x <- censored$X
  R <- censored$R
  
  # --- FULL DISPLAY of observed sample and removal vector ---
  cat("\nObserved X (m = 32):\n")
  # Print as a row of values, but for readability we can print in columns
  df_sample <- data.frame(i = 1:m, X = round(x, 4), R = R)
  print(df_sample, row.names = FALSE)
  
  # Alternatively, if you want a compact vector print:
  # cat("X =", paste(round(x,4), collapse=", "), "\n")
  # cat("R =", paste(R, collapse=", "), "\n")
  
  # MLE
  mle_est <- mle_estim(x, R, n, m, sch, true_par)
  if (is.null(mle_est)) {
    cat("MLE did not converge.\n")
    next
  }
  ci_mle <- mle_ci(x, R, n, m, sch, mle_est)
  
  cat("\nMLE estimates:\n")
  cat(sprintf("  alpha = %.5f\n", mle_est[1]))
  cat(sprintf("  theta = %.5f\n", mle_est[2]))
  cat(sprintf("  lambda = %.5f\n", mle_est[3]))
  
  cat("\n95% CI (MLE):\n")
  cat(sprintf("  alpha: [%.5f, %.5f]\n", ci_mle[1,1], ci_mle[1,2]))
  cat(sprintf("  theta: [%.5f, %.5f]\n", ci_mle[2,1], ci_mle[2,2]))
  cat(sprintf("  lambda: [%.5f, %.5f]\n", ci_mle[3,1], ci_mle[3,2]))
  
  # Bayes (Gamma22)
  init_par <- if (all(is.finite(mle_est)) && all(mle_est > 0)) mle_est else true_par
  bayes_gamma <- bayes_estim(x, R, n, m, sch, init_par, "Gamma22",
                             n_iter = 15000, burn = 5000, thin = 2,
                             proposal_sd = c(0.1, 0.01, 0.1))
  
  cat("\nBayes estimates (Gamma22 prior):\n")
  cat(sprintf("  alpha = %.5f\n", bayes_gamma$est[1]))
  cat(sprintf("  theta = %.5f\n", bayes_gamma$est[2]))
  cat(sprintf("  lambda = %.5f\n", bayes_gamma$est[3]))
  cat(sprintf("  Acceptance rate = %.3f\n", bayes_gamma$acceptance))
  
  cat("\n95% HPD (Gamma22):\n")
  cat(sprintf("  alpha: [%.5f, %.5f]\n", bayes_gamma$hpd[1,1], bayes_gamma$hpd[1,2]))
  cat(sprintf("  theta: [%.5f, %.5f]\n", bayes_gamma$hpd[2,1], bayes_gamma$hpd[2,2]))
  cat(sprintf("  lambda: [%.5f, %.5f]\n", bayes_gamma$hpd[3,1], bayes_gamma$hpd[3,2]))
  
  # Bayes (Jeffreys)
  bayes_jeff <- bayes_estim(x, R, n, m, sch, init_par, "Jeffreys",
                            n_iter = 15000, burn = 5000, thin = 2,
                            proposal_sd = c(0.1, 0.01, 0.1))
  
  cat("\nBayes estimates (Jeffreys prior):\n")
  cat(sprintf("  alpha = %.5f\n", bayes_jeff$est[1]))
  cat(sprintf("  theta = %.5f\n", bayes_jeff$est[2]))
  cat(sprintf("  lambda = %.5f\n", bayes_jeff$est[3]))
  cat(sprintf("  Acceptance rate = %.3f\n", bayes_jeff$acceptance))
  
  cat("\n95% HPD (Jeffreys):\n")
  cat(sprintf("  alpha: [%.5f, %.5f]\n", bayes_jeff$hpd[1,1], bayes_jeff$hpd[1,2]))
  cat(sprintf("  theta: [%.5f, %.5f]\n", bayes_jeff$hpd[2,1], bayes_jeff$hpd[2,2]))
  cat(sprintf("  lambda: [%.5f, %.5f]\n", bayes_jeff$hpd[3,1], bayes_jeff$hpd[3,2]))
  
  # Store results
  results[[sch]] <- list(
    x = x, R = R,
    mle = mle_est, ci_mle = ci_mle,
    bayes_gamma = bayes_gamma,
    bayes_jeff = bayes_jeff
  )
}

# --------------------- Summary Table (like Table 5 in Pathak et al.) ---------------------

cat("\n\n========== SUMMARY TABLE (similar to Table 5) ==========\n")
cat(sprintf("Complete data MLE: alpha=%.5f, theta=%.5f, lambda=%.5f\n\n",
            true_par[1], true_par[2], true_par[3]))

cat("Scheme  Parameter    MLE      Gamma22_Bayes   Jeffreys_Bayes\n")
cat("------------------------------------------------------------\n")

for (sch in schemes) {
  if (is.null(results[[sch]])) next
  res <- results[[sch]]
  
  # alpha
  cat(sprintf("%-7s %-10s %-8.4f   %-14.4f   %-14.4f\n",
              sch, "alpha", res$mle[1], res$bayes_gamma$est[1], res$bayes_jeff$est[1]))
  cat(sprintf("%-7s %-10s %-8.4f   %-14.4f   %-14.4f\n",
              "", "theta", res$mle[2], res$bayes_gamma$est[2], res$bayes_jeff$est[2]))
  cat(sprintf("%-7s %-10s %-8.4f   %-14.4f   %-14.4f\n",
              "", "lambda", res$mle[3], res$bayes_gamma$est[3], res$bayes_jeff$est[3]))
  cat("\n")
}

# --------------------- Comparison with complete data MLE (absolute errors) ---------------------

cat("\n========== ABSOLUTE ERRORS (vs complete data MLE) ==========\n")
for (sch in schemes) {
  if (is.null(results[[sch]])) next
  res <- results[[sch]]
  
  err_mle <- abs(res$mle - true_par)
  err_gamma <- abs(res$bayes_gamma$est - true_par)
  err_jeff <- abs(res$bayes_jeff$est - true_par)
  
  cat(sprintf("\n%s:\n", sch))
  cat(sprintf("  MLE:         alpha=%.5f, theta=%.5f, lambda=%.5f\n",
              err_mle[1], err_mle[2], err_mle[3]))
  cat(sprintf("  Gamma22:     alpha=%.5f, theta=%.5f, lambda=%.5f\n",
              err_gamma[1], err_gamma[2], err_gamma[3]))
  cat(sprintf("  Jeffreys:    alpha=%.5f, theta=%.5f, lambda=%.5f\n",
              err_jeff[1], err_jeff[2], err_jeff[3]))
}

cat("\n*** Analysis completed successfully ***\n")
