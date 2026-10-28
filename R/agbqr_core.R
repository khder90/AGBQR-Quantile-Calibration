ALPHA <- 0.05

# Estimation, covariance, and posterior simulation
check_loss <- function(u, tau) u * (tau - as.numeric(u < 0))
q_loss <- function(beta, y, X, tau) sum(check_loss(y - as.vector(X %*% beta), tau))
positive_matrix <- function(M, floor_rel = 1e-9) {
  M <- (M + t(M)) / 2
  eg <- eigen(M, symmetric = TRUE)
  floor_val <- max(max(abs(eg$values)), 1e-9) * floor_rel
  vals <- pmax(eg$values, floor_val)
  out <- eg$vectors %*% (vals * t(eg$vectors))
  (out + t(out)) / 2
}
inverse_spd <- function(A) chol2inv(chol(positive_matrix(A)))
find_qr <- function(y, X, tau) {
  if (nrow(X) <= ncol(X) + 2) stop("Not enough observations for QR.")
  fit <- quantreg::rq.fit.br(x = X, y = as.numeric(y), tau = tau)
  beta <- as.numeric(fit$coefficients)
  if (any(!is.finite(beta))) stop("Non-finite QR coefficients")
  names(beta) <- colnames(X)
  beta
}
default_lag <- function(n) max(0L, min(n - 2L, as.integer(floor(4 * (n/100)^(2/9)))))
score_hac <- function(scores, lag = NULL) {
  scores <- as.matrix(scores)
  n <- nrow(scores)
  if (is.null(lag)) lag <- default_lag(n)
  lag <- as.integer(min(max(0, lag), n - 2L))
  # Centered score autocovariances. Bartlett/Newey-West construction.
  S0 <- sweep(scores, 2, colMeans(scores), "-")
  V <- crossprod(S0) / n
  if (lag >= 1) {
    for (ll in seq_len(lag)) {
      Cll <- crossprod(S0[(ll + 1):n, , drop = FALSE],
                       S0[seq_len(n - ll), , drop = FALSE]) / n
      w <- 1 - ll / (lag + 1)
      V <- V + w * (Cll + t(Cll))
    }
  }
  positive_matrix(V)
}
hac_components <- function(beta_hat, y, X, tau, lag = NULL,
                           density = "weighted", bandwidth_multiplier = 1) {
  if (!density %in% c("pooled", "weighted")) stop("density must be pooled or weighted")
  n <- length(y)
  p <- ncol(X)
  u <- as.vector(y - X %*% beta_hat)
  psi <- tau - as.numeric(u < 0)
  I_hac <- score_hac(X * psi, lag = lag)
  h <- 1.06 * stats::sd(u) * n^(-1/5) * bandwidth_multiplier
  if (!is.finite(h) || h <= 1e-8) h <- max(1e-3, stats::mad(y) * 0.2)
  kw <- stats::dnorm(u / h) / h
  if (density == "pooled") {
    J <- mean(kw) * crossprod(X) / n
  } else {
    # Local residual-kernel estimator of E[f_{U|X}(0|X)XX'].
    # It does NOT estimate each conditional density separately.
    J <- crossprod(X, X * kw) / n
  }
  J <- positive_matrix(J)
  Jinv <- inverse_spd(J)
  V_asym <- Jinv %*% I_hac %*% Jinv
  V_coef <- positive_matrix(V_asym / n)
  dimnames(V_coef) <- list(colnames(X), colnames(X))
  list(J = J, I_hac = I_hac, V_asym = V_asym, V_coef = V_coef,
       bandwidth = h, lag = if (is.null(lag)) default_lag(n) else as.integer(lag),
       density = density, n = n)
}
# Prior-aware local covariance calibration.
calibrate_eta <- function(hac, prior_var = 100, bounds = c(1e-4, 1e4)) {
  target_diag <- pmax(diag(hac$V_coef), 1e-15)
  p <- length(target_diag)
  prior_precision <- if (is.infinite(prior_var)) 0 else 1/prior_var
  D <- function(eta) inverse_spd(eta * hac$n * hac$J + diag(prior_precision, p))
  objective <- function(log_eta) {
    d <- diag(D(exp(log_eta)))
    sum((d / target_diag - 1)^2)
  }
  opt <- optimize(objective, interval = log(bounds), tol = 1e-8)
  eta <- exp(opt$minimum)
  laplace <- D(eta)
  list(eta = eta, D = laplace, objective = opt$objective,
       near_boundary = eta < 1.001*bounds[1] || eta > .999*bounds[2])
}
# Cholesky proposal; proposal is symmetric, so standard Metropolis ratio suffices.
log_gibbs <- function(beta, y, X, tau, eta, prior_var) {
  lp <- if (is.infinite(prior_var)) 0 else -sum(beta^2)/(2*prior_var)
  -eta * q_loss(beta, y, X, tau) + lp
}
posterior_mode <- function(beta_start, y, X, tau, eta, prior_var) {
  # The check loss is nonsmooth; a derivative-free optimizer is used.
  objective <- function(b) -log_gibbs(b, y, X, tau, eta, prior_var)
  fit <- tryCatch(stats::optim(as.numeric(beta_start), objective,
                               method = "Nelder-Mead",
                               control = list(maxit = 1600, reltol = 1e-9)),
                  error = function(e) NULL)
  if (is.null(fit) || any(!is.finite(fit$par))) return(rep(NA_real_, length(beta_start)))
  out <- as.numeric(fit$par)
  names(out) <- names(beta_start)
  out
}
sample_gibbs <- function(y, X, tau, eta, beta_start, J, prior_var = 100,
                         M = 3000L, burn = 1000L, seed = NULL,
                         compute_mode = TRUE) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  n <- length(y)
  p <- ncol(X)
  if (M <= burn + 30) stop("M must exceed burn by at least 30")
  prior_prec <- if (is.infinite(prior_var)) 0 else 1/prior_var
  proposal_cov <- inverse_spd(eta * n * J + diag(prior_prec, p))
  proposal_chol <- chol(proposal_cov)
  beta <- as.numeric(beta_start)
  log_old <- log_gibbs(beta, y, X, tau, eta, prior_var)
  draws <- matrix(NA_real_, nrow = M - burn, ncol = p)
  colnames(draws) <- colnames(X)
  delta <- 1
  accepted <- 0L
  retained_accepted <- 0L
  for (m in seq_len(M)) {
    proposal <- beta + delta * as.vector(crossprod(proposal_chol, rnorm(p)))
    log_new <- log_gibbs(proposal, y, X, tau, eta, prior_var)
    log_alpha <- min(0, log_new - log_old)
    accept <- as.integer(log(runif(1)) < log_alpha)
    if (accept == 1L) {
      beta <- proposal
      log_old <- log_new
      accepted <- accepted + 1L
    }
    if (m <= burn) {
      gamma_m <- min(.01, m^(-.5))
      delta <- exp(max(-10, min(10, log(delta) + gamma_m*(accept - .234))))
    } else {
      draws[m - burn, ] <- beta
      retained_accepted <- retained_accepted + accept
    }
  }
  post_mean <- colMeans(draws)
  post_cov <- positive_matrix(stats::cov(draws))
  post_mode <- if (compute_mode) posterior_mode(beta_start, y, X, tau, eta, prior_var) else rep(NA_real_, p)
  names(post_mode) <- colnames(X)
  ess <- tryCatch(as.numeric(coda::effectiveSize(coda::as.mcmc(draws))),
                  error = function(e) rep(NA_real_, p))
  list(draws = draws, mean = post_mean, mode = post_mode, cov = post_cov,
       lower = apply(draws, 2, stats::quantile, probs = ALPHA / 2),
       upper = apply(draws, 2, stats::quantile, probs = 1-ALPHA / 2),
       ess = ess, acceptance = retained_accepted/(M-burn),
       scale = delta, eta = eta)
}
# Matrix comparisons explicitly expose the limits of scalar-rate calibration.
compare_covariances <- function(post_cov, hac_target, laplace_cov = NULL) {
  safe_norm <- function(A) sqrt(sum(A^2))
  diag_rmse <- sqrt(mean((diag(post_cov)/pmax(diag(hac_target),1e-12)-1)^2))
  full_gap <- safe_norm(post_cov - hac_target) / max(safe_norm(hac_target),1e-12)
  od <- function(A) { diag(A) <- 0; A }
  off_gap <- safe_norm(od(post_cov)-od(hac_target)) / max(safe_norm(od(hac_target)),1e-12)
  lap_gap <- if (is.null(laplace_cov)) NA_real_ else
    safe_norm(post_cov-laplace_cov)/max(safe_norm(laplace_cov),1e-12)
  c(Diagonal_RelRMSE = diag_rmse, Full_Frobenius_Gap = full_gap,
    OffDiagonal_Gap = off_gap, MCMC_to_Laplace_Gap = lap_gap)
}
intervals_wald <- function(beta, V) {
  se <- sqrt(pmax(diag(V), 0))
  z <- stats::qnorm(1-ALPHA/2)
  list(lower = as.numeric(beta)-z*se, upper = as.numeric(beta)+z*se)
}
# Conventional non-iid quantile regression intervals.
qr_nid_intervals <- function(y, X, tau) {
  df <- data.frame(y = y, X[, -1, drop = FALSE], check.names = TRUE)
  fit <- quantreg::rq(y ~ ., data = df, tau = tau)
  ss <- summary(fit, se = "nid")$coefficients
  s <- as.numeric(ss[, 2])
  est <- as.numeric(stats::coef(fit))
  z <- stats::qnorm(1-ALPHA/2)
  list(lower = est-z*s, upper = est+z*s)
}

# Numerical consistency checks.
stopifnot(all(check_loss(c(-2,-1,0,1,2),.5)>=0))
stopifnot(abs(sum(diag(inverse_spd(diag(1,3))))-3)<1e-8)
stopifnot(default_lag(60) >= 0, default_lag(60) <= 58)


# Simulation data-generating processes
ar_predictor <- function(n, rho) {
  x <- numeric(n)
  x[1] <- stats::rnorm(1)
  for (t in 2:n) x[t] <- rho*x[t-1] + stats::rnorm(1)
  as.numeric(scale(x))
}
simulate_data <- function(n, dgp) {
  X <- cbind(Intercept = 1, x1 = ar_predictor(n,.4),
             x2 = ar_predictor(n,.5), x3 = ar_predictor(n,.3))
  beta <- c(Intercept = .5, x1 = 1, x2 = -.7, x3 = .4)
  if (dgp == "Gaussian") {
    u <- rnorm(n)
    target <- function(tau) beta + c(Intercept = qnorm(tau), x1 = 0, x2 = 0, x3 = 0)
    y <- as.vector(X %*% beta + u)
  } else if (dgp == "HeavyTail") {
    u <- stats::rt(n, df = 3) / sqrt(3)
    target <- function(tau) beta + c(Intercept = stats::qt(tau, 3)/sqrt(3),x1=0,x2=0,x3=0)
    y <- as.vector(X %*% beta + u)
  } else if (dgp == "Heteroskedastic") {
    X <- cbind(Intercept = 1, x1 = runif(n,-1,1), x2 = runif(n,-1,1), x3 = runif(n,-1,1))
    delta <- c(Intercept = 1, x1 = .3, x2 = .2, x3 = -.1)
    u <- rnorm(n)
    y <- as.vector(X %*% beta + (X %*% delta) * u)
    target <- function(tau) beta + delta * stats::qnorm(tau)
  } else if (dgp == "SerialAR") {
    rho <- .65
    u <- numeric(n)
    u[1] <- rnorm(1)
    for (t in 2:n) u[t] <- rho*u[t-1] + sqrt(1-rho^2)*rnorm(1)
    y <- as.vector(X %*% beta + u)
    # Unconditional conditional-on-X target. AR error dependent through time.
    target <- function(tau) beta + c(Intercept = qnorm(tau),x1=0,x2=0,x3=0)
  } else stop("Unknown DGP: ", dgp)
  list(y = y, X = X, truth = target)
}

