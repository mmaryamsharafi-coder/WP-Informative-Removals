# Revised Weibull-Poisson simulation: ML vs informative Bayes vs noninformative Bayes
# n=60, m=15, B=2000; Table 4 logic preserved
# True parameters: alpha=1.5, theta=1.0, lambda=3.0
# Priors: informative Gamma(2,2) for alpha, theta, lambda
#          noninformative independent Jeffreys-type prior 1/(alpha*theta*lambda)

if (!requireNamespace("coda", quietly = TRUE)) install.packages("coda", repos = "https://cloud.r-project.org")
library(coda)
if (!requireNamespace("numDeriv", quietly = TRUE)) install.packages("numDeriv")
library(numDeriv)
set.seed(2026)

# --------------------- User settings ---------------------
B       <- 2000        # number of replications (reduce for testing)
n       <- 60          # changed from 40
m       <- 15           # changed from 10
p_bin   <- 0.5
true_par <- c(alpha = 1.5, theta = 1.0, lambda = 3.0)   # NEW true values

mcmc_iter <- 6000
burn      <- 2000
thin      <- 5

# --------------------- Core functions (unchanged) ---------------------
log_expm1_pos <- function(x) {
  ifelse(x < log(2), log(expm1(x)), x + log1p(-exp(-x)))
}

log_dwp <- function(x, a, t, l) {
  z <- t * x^a
  log(a) + log(t) + log(l) + (a - 1) * log(x) - z + l * exp(-z) - log_expm1_pos(l)
}

log_swp <- function(x, a, t, l) {
  z <- t * x^a
  log_expm1_pos(l * exp(-z)) - log_expm1_pos(l)
}

qwp <- function(u, a, t, l) {
  q <- log1p((1 - u) * expm1(l)) / l
  (-log(q) / t)^(1 / a)
}

remaining_units <- function(R, i, n, m) {
  n - m - sum(R[seq_len(max(0, i - 1))])
}

# ---------- Censoring scheme generators ----------
sample_R_TP <- function(n, m, l) {
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

sample_R_TDW <- function(n, m, a, t) {
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

sample_R_BIN <- function(n, m, p = 0.5) {
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

# ---------- Progressive Type-II sample generator ----------
rprog_wp <- function(n, m, R, a, t, l) {
  u <- runif(m)
  y <- numeric(m)
  surv <- 1
  for (i in seq_len(m)) {
    surv <- surv * runif(1)^(1 / (n - sum(R[seq_len(i - 1)]) - i + 1))
    y[i] <- qwp(1 - surv, a, t, l)
  }
  sort(y)
}

# ---------- Log-likelihood of removal scheme ----------
log_rem <- function(R, scheme, a, t, l, p = 0.5) {
  if (scheme == "BIN") return(0)
  if (scheme == "TP") {
    rem <- n - m
    v <- 0
    for (i in seq_len(m - 1)) {
      av <- remaining_units(R, i, n, m)
      k <- R[i]
      w <- dpois(0:av, l)
      v <- v + log(w[k + 1] / sum(w))
      rem <- rem - k
    }
    return(v)
  }
  # TDW
  rem <- n - m
  v <- 0
  for (i in seq_len(m - 1)) {
    av <- remaining_units(R, i, n, m)
    k <- R[i]
    den <- 1 - exp(-t * (av + 1)^a)
    num <- exp(-t * k^a) - exp(-t * (k + 1)^a)
    v <- v + log(pmax(num, 1e-300)) - log(pmax(den, 1e-300))
    rem <- rem - k
  }
  v
}

# ---------- Full log-likelihood (log scale parameters) ----------
loglik <- function(eta, x, R, scheme) {
  q <- exp(eta)
  if (any(!is.finite(q)) || any(q > 50)) return(-Inf)
  sum(log_dwp(x, q[1], q[2], q[3])) +
    sum(R * log_swp(x, q[1], q[2], q[3])) +
    log_rem(R, scheme, q[1], q[2], q[3], p_bin)
}

# ---------- MLE fit ----------
mle_fit <- function(x, R, scheme) {
  # UPDATED starting points for new true parameters
  st <- rbind(log(c(1.5, 1.0, 3.0)),   # true values
              log(c(1, 1, 2)),
              log(c(2, 2, 5)))
  z <- lapply(seq_len(nrow(st)), function(i) {
    tryCatch(
      optim(st[i, ], function(e) -loglik(e, x, R, scheme),
            method = "BFGS", hessian = TRUE,
            control = list(maxit = 3000)),
      error = function(e) NULL
    )
  })
  z <- Filter(Negate(is.null), z)
  if (!length(z)) return(NULL)
  fit <- z[[which.min(sapply(z, `[[`, "value"))]]
  if (fit$convergence > 1) return(NULL)
  est <- exp(fit$par)
  if (any(!is.finite(est)) || any(est > 50)) return(NULL)

  H <- tryCatch(numDeriv::hessian(function(e) -loglik(e, x, R, scheme), fit$par),
                error = function(e) NULL)
  se <- rep(NA, 3)
  if (!is.null(H) && all(is.finite(H))) {
    V <- tryCatch(solve(H), error = function(e) NULL)
    if (!is.null(V) && all(diag(V) > 0))
      se <- est * sqrt(diag(V))
  }
  ci <- cbind(lower = est - qnorm(0.975) * se,
              upper = est + qnorm(0.975) * se)
  rownames(ci) <- names(true_par)
  list(est = est, ci = ci)
}

# ---------- Log-prior ----------
logprior <- function(q, type) {
  if (any(q <= 0) || any(!is.finite(q))) return(-Inf)
  if (type == "Gamma22") {
    sum(dgamma(q, 2, rate = 2, log = TRUE))
  } else {  # Jeffreys
    -sum(log(q))
  }
}

# ---------- Bayesian fit (MCMC) ----------
bayes_fit <- function(x, R, scheme, type, init) {
  eta <- log(init)
  target <- function(e) loglik(e, x, R, scheme) + logprior(exp(e), type) + sum(e)
  cur <- target(eta)
  out <- matrix(NA, mcmc_iter, 3)
  ac <- 0
  for (i in seq_len(mcmc_iter)) {
    pr <- eta + rnorm(3, 0, c(0.1, 0.1, 0.1))
    v <- target(pr)
    if (is.finite(v) && log(runif(1)) < v - cur) {
      eta <- pr
      cur <- v
      ac <- ac + 1
    }
    out[i, ] <- exp(eta)
  }
  keep <- out[seq(burn + 1, mcmc_iter, by = thin), , drop = FALSE]
  colnames(keep) <- names(true_par)
  h <- t(apply(keep, 2, function(v) as.numeric(HPDinterval(mcmc(v), 0.95))))
  colnames(h) <- c("lower", "upper")
  list(est = colMeans(keep), hpd = h)
}

# ---------- One replication ----------
one_rep <- function(s) {
  R <- switch(s,
              TP  = sample_R_TP(n, m, true_par[3]),
              TDW = sample_R_TDW(n, m, true_par[1], true_par[2]),
              BIN = sample_R_BIN(n, m, p_bin))
  x <- rprog_wp(n, m, R, true_par[1], true_par[2], true_par[3])
  ml <- mle_fit(x, R, s)
  init <- if (is.null(ml)) c(1.5, 1.0, 3.0) else ml$est
  g <- bayes_fit(x, R, s, "Gamma22", init)
  u <- bayes_fit(x, R, s, "Jeffreys", init)
  list(
    ml   = if (is.null(ml)) rep(NA, 3) else ml$est,
    ci   = if (is.null(ml)) matrix(NA, 3, 2) else ml$ci,
    g    = g,
    u    = u,
    time = max(x)
  )
}

# ---------- Run all replications for each scheme ----------
res <- setNames(
  lapply(c("TP", "TDW", "BIN"), function(s) {
    lapply(seq_len(B), function(i) one_rep(s))
  }),
  c("TP", "TDW", "BIN")
)

# ---------- Table 1: Bias and RMSE ----------
T1 <- do.call(rbind, lapply(names(res), function(s) {
  z <- res[[s]]
  do.call(rbind, lapply(1:3, function(j) {
    a <- sapply(z, function(q) q$ml[j])
    g <- sapply(z, function(q) q$g$est[j])
    u <- sapply(z, function(q) q$u$est[j])
    data.frame(
      Scheme = s,
      Parameter = names(true_par)[j, drop = TRUE],
      ML_Bias = mean(a - true_par[j], na.rm = TRUE),
      ML_RMSE = sqrt(mean((a - true_par[j])^2, na.rm = TRUE)),
      Gamma22_Bias = mean(g - true_par[j], na.rm = TRUE),
      Gamma22_RMSE = sqrt(mean((g - true_par[j])^2, na.rm = TRUE)),
      Jeffreys_Bias = mean(u - true_par[j], na.rm = TRUE),
      Jeffreys_RMSE = sqrt(mean((u - true_par[j])^2, na.rm = TRUE))
    )
  }))
}))

# ---------- Table 3: CP and AW for intervals ----------
metric <- function(A, j, t) {
  lo <- A[, j, 1]
  hi <- A[, j, 2]
  ok <- is.finite(lo) & is.finite(hi)
  # ordering band for avoiding  lower > upper
  lo_ok <- pmin(lo[ok], hi[ok])
  hi_ok <- pmax(lo[ok], hi[ok])
  c(CP = mean(lo_ok <= t & hi_ok >= t),
    AW = mean(hi_ok - lo_ok))
}

T3 <- do.call(rbind, lapply(names(res), function(s) {
  z <- res[[s]]
  do.call(rbind, lapply(1:3, function(j) {
    # MLE CI
    A <- array(unlist(lapply(z, function(q) q$ci)), c(2, 3, length(z)))
    # Gamma22 HPD
    G <- array(unlist(lapply(z, function(q) q$g$hpd)), c(2, 3, length(z)))
    # Jeffreys HPD
    U <- array(unlist(lapply(z, function(q) q$u$hpd)), c(2, 3, length(z)))
    a <- metric(aperm(A, c(3, 2, 1)), j, true_par[j])
    g <- metric(aperm(G, c(3, 2, 1)), j, true_par[j])
    u <- metric(aperm(U, c(3, 2, 1)), j, true_par[j])
    data.frame(
      Scheme = s,
      Parameter = names(true_par)[j],
      ML_CP = a[1],
      ML_AW = a[2],
      Gamma22_HPD_CP = g[1],
      Gamma22_HPD_AW = g[2],
      Jeffreys_HPD_CP = u[1],
      Jeffreys_HPD_AW = u[2]
    )
  }))
}))

# ---------- Table 4: Expected test time ----------
T4 <- data.frame(
  Scheme = names(res),
  Expected_Test_Time = sapply(res, function(z) mean(sapply(z, `[[`, "time"), na.rm = TRUE))
)

# ---------- Write CSV files (updated names) ----------
write.csv(T1, "Table1_n60_m15_alpha15_theta10_lambda30.csv", row.names = FALSE)
write.csv(T3, "Table3_n60_m15_alpha15_theta10_lambda30.csv", row.names = FALSE)
write.csv(T4, "Table4_n60_m15_alpha15_theta10_lambda30.csv", row.names = FALSE)

# ---------- Print results ----------
print(T1)
print(T3)
print(T4)