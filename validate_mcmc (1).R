# Four-chain numerical validation for selected quantile-regression settings.
# Run from the repository root after installing the packages listed in README.md.

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_ROOT, "R", "agbqr_core.R"))

required <- c("readxl", "quantreg", "coda", "posterior", "ggplot2", "tidyr", "openxlsx")
for (pkg in required) if (!requireNamespace(pkg, quietly = TRUE)) {
  stop("Install required package: ", pkg)
}

N_ITER <- 20000L
BURN_IN <- 5000L
N_CHAINS <- 4L
PRIOR_VARIANCE <- 100
OUT <- file.path(PROJECT_ROOT, "results", "validation")
for (folder in c(OUT, file.path(OUT, "chains"), file.path(OUT, "figures"))) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

input <- as.data.frame(readxl::read_excel(file.path(PROJECT_ROOT, "data", "Syrian_Inflation_1961_2025.xlsx")))
stopifnot(all(c("YEAR", "EXR", "TRDGDP", "INF", "M2GDP") %in% names(input)))
input <- input[order(input$YEAR), ]
scaled <- function(x) as.numeric((x - mean(x)) / stats::sd(x))
empirical <- list(
  y = scaled(input$INF),
  X = cbind(Intercept = 1, LEXR = scaled(log(input$EXR)),
            M2GDP = scaled(input$M2GDP), TRDGDP = scaled(input$TRDGDP))
)

cases <- rbind(
  do.call(rbind, lapply(c("Gaussian", "HeavyTail", "Heteroskedastic", "SerialAR"),
                        function(dgp) data.frame(DGP = dgp, N = 60L, Tau = c(0.10, 0.90)))),
  data.frame(DGP = c("Gaussian", "HeavyTail", "Heteroskedastic", "SerialAR"), N = 120L, Tau = 0.50),
  data.frame(DGP = "Empirical", N = nrow(empirical$X), Tau = c(0.10, 0.25, 0.50, 0.75, 0.90))
)
stopifnot(nrow(cases) == 17L)

summaries <- vector("list", nrow(cases))
parameters <- vector("list", nrow(cases))
for (i in seq_len(nrow(cases))) {
  setting <- cases[i, ]
  label <- sprintf("%s_n%d_tau%02d", setting$DGP, setting$N, round(setting$Tau * 100))
  message("Validating ", label, " (", i, "/", nrow(cases), ")")
  stored <- file.path(OUT, "chains", paste0(label, ".rds"))
  if (file.exists(stored)) {
    chains <- readRDS(stored)
  } else {
    set.seed(510000L + i)
    sample <- if (setting$DGP == "Empirical") empirical else simulate_data(setting$N, setting$DGP)
    beta_qr <- find_qr(sample$y, sample$X, setting$Tau)
    components <- hac_components(beta_qr, sample$y, sample$X, setting$Tau,
                                 density = "weighted")
    rate <- calibrate_eta(components, prior_var = PRIOR_VARIANCE)$eta
    chains <- lapply(seq_len(N_CHAINS), function(j) {
      fit <- sample_gibbs(sample$y, sample$X, setting$Tau,
                          eta = rate, beta_start = beta_qr, J = components$J,
                          prior_var = PRIOR_VARIANCE,
                          M = N_ITER, burn = BURN_IN,
                          seed = 520000L + i * 100L + j,
                          compute_mode = FALSE)
      list(draws = fit$draws, acceptance = fit$acceptance)
    })
    saveRDS(chains, stored)
  }
  stopifnot(length(chains) == N_CHAINS,
            all(vapply(chains, function(x) nrow(x$draws) == N_ITER - BURN_IN, logical(1))))
  npar <- ncol(chains[[1]]$draws)
  arrays <- array(NA_real_, dim = c(N_ITER - BURN_IN, N_CHAINS, npar),
                  dimnames = list(NULL, paste0("chain_", 1:N_CHAINS), colnames(chains[[1]]$draws)))
  for (j in seq_len(N_CHAINS)) arrays[, j, ] <- chains[[j]]$draws
  diagnostic <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws_array(arrays),
    "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"
  ))
  diagnostic$MCSE_over_SD <- diagnostic$mcse_mean / diagnostic$sd
  diagnostic$Case <- label
  diagnostic$Pass <- with(diagnostic, is.finite(rhat) & rhat < 1.01 &
    ess_bulk >= 400 & ess_tail >= 400 & MCSE_over_SD <= 0.05)
  parameters[[i]] <- diagnostic
  summaries[[i]] <- data.frame(
    Case = label, DGP = setting$DGP, N = setting$N, Tau = setting$Tau,
    Max_Rhat = max(diagnostic$rhat), Min_Bulk_ESS = min(diagnostic$ess_bulk),
    Min_Tail_ESS = min(diagnostic$ess_tail), Max_MCSE_SD = max(diagnostic$MCSE_over_SD),
    Mean_Acceptance = mean(vapply(chains, `[[`, numeric(1), "acceptance")),
    All_Pass = all(diagnostic$Pass)
  )
  if (setting$DGP == "Empirical" && setting$Tau %in% c(0.10, 0.50, 0.90)) {
    plot_data <- do.call(rbind, lapply(seq_len(N_CHAINS), function(j) {
      keep <- seq.int(1L, N_ITER - BURN_IN, by = 10L)
      wide <- as.data.frame(chains[[j]]$draws[keep, , drop = FALSE])
      wide$Iteration <- keep
      wide$Chain <- paste0("Chain ", j)
      tidyr::pivot_longer(wide, cols = -c(Iteration, Chain),
                          names_to = "Coefficient", values_to = "Value")
    }))
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(Iteration, Value, colour = Chain)) +
      ggplot2::geom_line(linewidth = 0.3, alpha = 0.7) +
      ggplot2::facet_wrap(~Coefficient, ncol = 2, scales = "free_y") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(x = "Retained iteration", y = "Posterior coefficient")
    ggplot2::ggsave(file.path(OUT, "figures", paste0("Trace_", label, ".png")),
                    p, width = 12, height = 8, dpi = 600, units = "in", bg = "white")
  }
}
openxlsx::write.xlsx(
  list(Case_Summary = do.call(rbind, summaries),
       Parameter_Diagnostics = do.call(rbind, parameters)),
  file.path(OUT, "MCMC_Diagnostics.xlsx"), overwrite = TRUE
)
writeLines(capture.output(sessionInfo()), file.path(OUT, "Session_Info.txt"))
message("Completed independent-chain validation; see results/validation/MCMC_Diagnostics.xlsx")
