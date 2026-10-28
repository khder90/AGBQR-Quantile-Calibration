# Adaptive Gibbs posteriors for quantile regression
# Numerical analysis and Syrian inflation application
#
# From the repository root:
#   source("R/run_full_analysis.R")
#
# Data: data/Syrian_Inflation_1961_2025.xlsx
# Output: results/

RUN_MODE <- "full"
BASE_SEED <- 20261021L
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
INPUT_FILE <- file.path(PROJECT_ROOT, "data", "Syrian_Inflation_1961_2025.xlsx")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "results")
RUN_MONTE_CARLO <- TRUE
RUN_SENSITIVITY <- TRUE
RUN_EMPIRICAL <- TRUE
RUN_ROLLING_VALIDATION <- TRUE

if (!RUN_MODE %in% c("pilot", "full")) stop("RUN_MODE must be pilot or full")
IS_FULL <- identical(RUN_MODE, "full")
DGP_SET <- if (IS_FULL) c("Gaussian", "HeavyTail", "Heteroskedastic", "SerialAR") else c("Gaussian", "Heteroskedastic", "SerialAR")
SAMPLE_SIZES <- if (IS_FULL) c(60L, 120L, 250L, 500L) else c(60L, 120L)
TAUS <- c(0.10, 0.25, 0.50, 0.75, 0.90)
MC_REPS <- if (IS_FULL) 300L else 8L
MC_M <- if (IS_FULL) 3000L else 700L
MC_BURN <- if (IS_FULL) 1000L else 250L
MC_CHUNK_REPS <- if (IS_FULL) 10L else 4L
SENS_REPS <- if (IS_FULL) 75L else 5L
SENS_M <- if (IS_FULL) 1800L else 650L
SENS_BURN <- if (IS_FULL) 600L else 200L
EMP_M <- if (IS_FULL) 30000L else 4000L
EMP_BURN <- if (IS_FULL) 10000L else 1200L
OOS_M <- if (IS_FULL) 3000L else 900L
OOS_BURN <- if (IS_FULL) 1000L else 300L
OOS_INITIAL <- 40L
PRIOR_VAR <- 100
DIFFUSE_VAR <- 1e6
DENSITY_DEFAULT <- "weighted"  
BW_MULT_DEFAULT <- 1
ALPHA <- 0.05


# Package dependencies and output directories
options(repos = c(CRAN = "https://cloud.r-project.org"),
        stringsAsFactors = FALSE, scipen = 999, warn = 1)
required <- c("readxl", "quantreg", "dplyr", "tidyr", "ggplot2", "openxlsx", "coda")
for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
# Attach dplyr/tidyr for the pipe operator and explicit dataframe verbs.
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
if (!file.exists(INPUT_FILE)) stop("Missing input file: ", INPUT_FILE)
for (sub in c("data", "paper_tables", "paper_figures", "appendix_tables",
              "appendix_figures", "mc_data", "mc_data/blocks_full", "mc_data/blocks_pilot",
              "sensitivity", "oos", "logs", "notes")) {
  dir.create(file.path(OUTPUT_DIR, sub), recursive = TRUE, showWarnings = FALSE)
}
MODE_TAG <- if (IS_FULL) "FULL" else "PILOT"
LOG_FILE <- file.path(OUTPUT_DIR, "logs", paste0("Execution_", MODE_TAG, ".log"))
log_line <- function(...) {
  message <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(message, "\n")
  cat(message, "\n", file = LOG_FILE, append = TRUE)
}
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
write_xlsx <- function(x, path) openxlsx::write.xlsx(x, path, overwrite = TRUE)
fmt4 <- function(x) if (is.numeric(x)) round(x, 4) else x
ROUND4 <- function(df) dplyr::mutate(df, dplyr::across(where(is.numeric), fmt4))
plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 16) + ggplot2::theme(
    plot.title = ggplot2::element_text(size = 17, face = "bold"),
    axis.title = ggplot2::element_text(size = 16),
    axis.text = ggplot2::element_text(size = 16),
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    legend.text = ggplot2::element_text(size = 16),
    legend.title = ggplot2::element_text(size = 16),
    panel.grid.minor = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(size = 13))
}
save_figure <- function(plot, file, width = 12, height = 8) {
  # 600 DPI = 600 dots per inch; 16-point text for axes, labels and ticks.
  ggplot2::ggsave(file, plot, width = width, height = height,
                  units = "in", dpi = 600, bg = "white", limitsize = FALSE)
}
COLORS <- c("QR-Wald-HAC" = "#1B9E77", "QR-NID" = "#7570B3",
            "Fixed-Normal" = "#D95F02", "Adaptive-Normal" = "#E7298A",
            "Fixed-Diffuse" = "#66A61E", "Adaptive-Diffuse" = "#1F78B4")
log_line("Started mode=", MODE_TAG, ". Results root: ", OUTPUT_DIR)
writeLines(c(paste("Mode:", RUN_MODE), paste("MC repetitions:", MC_REPS),
             paste("DGPs:", paste(DGP_SET, collapse=", ")),
             paste("Sample sizes:", paste(SAMPLE_SIZES, collapse=", ")),
             paste("MC MCMC:", MC_M, "burn:", MC_BURN),
             paste("Empirical MCMC:", EMP_M, "burn:", EMP_BURN),
             "Figures: 600 dpi; axis labels/ticks 16 pt",
             "Simulation estimates are conditional on the stated design and sampling settings."),
           file.path(OUTPUT_DIR, paste0("RUN_SETTINGS_", MODE_TAG, ".txt")))

source(file.path(PROJECT_ROOT, "R", "agbqr_core.R"), local = FALSE)

# ====================== INPUT AND PREPROCESSING ==============================
raw_input <- as.data.frame(readxl::read_excel(INPUT_FILE))
# Normalize input column names.
canonical_header <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
keys <- canonical_header(names(raw_input))
map <- c(YEAR="YEAR", EX="EXR", EXR="EXR", OPEN="TRDGDP",
         TRDGDP="TRDGDP", INF="INF", M2GDP="M2GDP")
for (i in seq_along(keys)) if (keys[i] %in% names(map)) names(raw_input)[i] <- map[[keys[i]]]
if (anyDuplicated(names(raw_input))) stop("Ambiguous duplicate Excel column names after normalization")
needed <- c("YEAR", "EXR", "TRDGDP", "INF", "M2GDP")
if (!all(needed %in% names(raw_input))) {
  stop("The Excel file is missing columns: ", paste(setdiff(needed, names(raw_input)), collapse = ", "))
}
raw_input <- raw_input[, needed]
for (nm in needed) raw_input[[nm]] <- as.numeric(raw_input[[nm]])
raw_input <- raw_input[order(raw_input$YEAR), , drop = FALSE]
if (any(!is.finite(as.matrix(raw_input)))) stop("Input contains non-finite values")
if (anyDuplicated(raw_input$YEAR)) stop("YEAR is not unique")
if (any(raw_input$EXR <= 0)) stop("EXR must be strictly positive for log transform")
raw_input$LEXR <- log(raw_input$EXR)
standardize <- function(v) {
  s <- stats::sd(v)
  if (!is.finite(s) || s < 1e-12) stop("Predictor or outcome has near-zero SD")
  (v - mean(v)) / s
}
raw_input$INF_s <- standardize(raw_input$INF)
raw_input$LEXR_s <- standardize(raw_input$LEXR)
raw_input$M2GDP_s <- standardize(raw_input$M2GDP)
raw_input$TRDGDP_s <- standardize(raw_input$TRDGDP)
y_emp <- raw_input$INF_s
X_emp <- cbind(Intercept = 1, LEXR = raw_input$LEXR_s,
               M2GDP = raw_input$M2GDP_s, TRDGDP = raw_input$TRDGDP_s)
write_csv(raw_input, file.path(OUTPUT_DIR, "data", "Clean_Inflation_Data.csv"))
log_line("Data imported: n=", nrow(raw_input), ", years ", min(raw_input$YEAR), " to ", max(raw_input$YEAR))
# Optional descriptive time-series diagnostics (not proof of model validity).
if (!requireNamespace("tseries", quietly = TRUE)) install.packages("tseries")
stationarity <- lapply(c("INF_s","LEXR_s","M2GDP_s","TRDGDP_s"),function(nm) {
  vec <- raw_input[[nm]]
  a <- tryCatch(tseries::adf.test(vec),error=function(e)NULL)
  k <- tryCatch(tseries::kpss.test(vec),error=function(e)NULL)
  data.frame(Variable=nm,
    ADF_Statistic=if(is.null(a)) NA_real_ else unname(a$statistic),
    ADF_p=if(is.null(a)) NA_real_ else a$p.value,
    KPSS_Statistic=if(is.null(k)) NA_real_ else unname(k$statistic),
    KPSS_p=if(is.null(k)) NA_real_ else k$p.value)
})
write_xlsx(ROUND4(dplyr::bind_rows(stationarity)),
 file.path(OUTPUT_DIR,"appendix_tables","Appendix_Stationarity_Diagnostics.xlsx"))

# Descriptives and empirical plots, available in both modes.
descriptive <- data.frame(Variable = c("INF", "LEXR", "M2GDP", "TRDGDP"))
for (metric in c("N", "Mean", "Median", "SD", "Min", "Max")) descriptive[[metric]] <- NA_real_
for (i in seq_len(nrow(descriptive))) {
  v <- raw_input[[descriptive$Variable[i]]]
  descriptive[i, 2:7] <- list(length(v), mean(v), stats::median(v),
                              stats::sd(v), min(v), max(v))
}
write_xlsx(ROUND4(descriptive), file.path(OUTPUT_DIR, "paper_tables", "Table_E1_Descriptive_Statistics.xlsx"))
plot_ts <- tidyr::pivot_longer(raw_input[, c("YEAR", "INF_s", "LEXR_s", "M2GDP_s", "TRDGDP_s")],
                               cols = -YEAR, names_to = "Variable", values_to = "Value")
p_ts <- ggplot2::ggplot(plot_ts, ggplot2::aes(YEAR, Value, colour = Variable)) +
  ggplot2::geom_line(linewidth = .9) + ggplot2::facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  ggplot2::scale_colour_brewer(palette = "Dark2") +
  ggplot2::labs(x = "Year", y = "Standardized value", title = "Syrian inflation and covariates") +
  plot_theme() + ggplot2::theme(legend.position = "none")
save_figure(p_ts, file.path(OUTPUT_DIR, "paper_figures", "Figure_E1_Standardized_Time_Series_600dpi.png"))

# Each replication uses the same HAC covariance for Wald inference and calibration.
METHODS_POSTERIOR <- c("Fixed-Normal", "Adaptive-Normal", "Fixed-Diffuse", "Adaptive-Diffuse")
run_rep_tau <- function(y, X, tau, true_beta, mc_seed, M, burn,
                        density = DENSITY_DEFAULT, bandwidth_multiplier = 1,
                        lag = NULL, include_nid = TRUE, compute_modes = TRUE) {
  n <- length(y)
  qr_t0 <- proc.time()[["elapsed"]]
  qr_beta <- find_qr(y, X, tau)
  hh <- hac_components(qr_beta, y, X, tau, lag = lag,
                       density = density, bandwidth_multiplier = bandwidth_multiplier)
  wald <- intervals_wald(qr_beta, hh$V_coef)
  qr_hac_seconds <- proc.time()[["elapsed"]] - qr_t0
  p <- ncol(X)
  truth <- as.numeric(true_beta[colnames(X)])
  wald_joint <- as.numeric(t(qr_beta - truth) %*% inverse_spd(hh$V_coef) %*% (qr_beta - truth) <= qchisq(.95,p))
  rows <- list()
  add_method <- function(method, estimate, lower, upper, mode = NULL,
                         eta = NA_real_, accept = NA_real_, ess = NULL,
                         cov_gaps = NULL, joint_covered = NA_real_, time_seconds = NA_real_) {
    if (is.null(mode)) mode <- rep(NA_real_,p)
    if (is.null(ess)) ess <- rep(NA_real_,p)
    if (is.null(cov_gaps)) cov_gaps <- rep(NA_real_,4)
    dat <- data.frame(Method = method, Variable = colnames(X),
                      Estimate = as.numeric(estimate), Mode = as.numeric(mode),
                      True_Beta = truth, Lower = as.numeric(lower), Upper = as.numeric(upper),
                      Covered = as.integer(truth >= lower & truth <= upper),
                      Interval_Width = as.numeric(upper-lower),
                      Eta = eta, Acceptance = accept, ESS = ess,
                      Diagonal_RelRMSE = cov_gaps[1], Full_Frobenius_Gap = cov_gaps[2],
                      OffDiagonal_Gap = cov_gaps[3], MCMC_to_Laplace_Gap = cov_gaps[4],
                      Joint_Covered = joint_covered,
                      HAC_Lag = hh$lag, Density_Estimator = density,
                      Bandwidth = hh$bandwidth, Time_Seconds = time_seconds)
    rows[[length(rows)+1L]] <<- dat
  }
  add_method("QR-Wald-HAC",qr_beta,wald$lower,wald$upper,
             joint_covered = wald_joint,time_seconds = qr_hac_seconds)
  if (include_nid) {
    nid_t0 <- proc.time()[["elapsed"]]
    nid <- tryCatch(qr_nid_intervals(y,X,tau),error = function(e) NULL)
    nid_elapsed <- proc.time()[["elapsed"]]-nid_t0
    if (!is.null(nid) && all(is.finite(c(nid$lower,nid$upper))))
      add_method("QR-NID",qr_beta,nid$lower,nid$upper,time_seconds=nid_elapsed)
  }
  normal_cal <- calibrate_eta(hh, PRIOR_VAR)
  diffuse_cal <- calibrate_eta(hh, DIFFUSE_VAR)
  meta <- data.frame(Tau = tau, Eta_Normal = normal_cal$eta,
                     Eta_Diffuse = diffuse_cal$eta,
                     Objective_Normal = normal_cal$objective,
                     Objective_Diffuse = diffuse_cal$objective,
                     Near_Boundary_Normal = normal_cal$near_boundary,
                     Near_Boundary_Diffuse = diffuse_cal$near_boundary,
                     HAC_Lag = hh$lag, Bandwidth = hh$bandwidth)
  for (j in seq_along(METHODS_POSTERIOR)) {
    method <- METHODS_POSTERIOR[j]
    pv <- if (grepl("Diffuse", method)) DIFFUSE_VAR else PRIOR_VAR
    cal <- if (grepl("Diffuse", method)) diffuse_cal else normal_cal
    eta <- if (grepl("^Fixed",method)) 1 else cal$eta
    mcmc_t0 <- proc.time()[["elapsed"]]
    post <- sample_gibbs(y,X,tau,eta,qr_beta,hh$J,prior_var = pv,
                         M = M,burn = burn, seed = mc_seed + j,
                         compute_mode = compute_modes)
    mcmc_elapsed <- proc.time()[["elapsed"]]-mcmc_t0
    gaps <- compare_covariances(post$cov, hh$V_coef,
                                if (grepl("^Adaptive",method)) cal$D else NULL)
    joint <- as.numeric(t(post$mean-truth) %*% inverse_spd(post$cov) %*%
                          (post$mean-truth) <= qchisq(.95,p))
    add_method(method, post$mean, post$lower,post$upper, mode = post$mode,
               eta = eta,accept=post$acceptance,ess=post$ess,
               cov_gaps=gaps,joint_covered=joint,time_seconds=mcmc_elapsed)
  }
  list(rows = dplyr::bind_rows(rows), metadata = meta)
}

# ====================== MAIN MONTE CARLO ====================================
mc_path <- file.path(OUTPUT_DIR,"mc_data",paste0("MC_",MODE_TAG,"_AllRows.rds"))
errors <- list()
if (RUN_MONTE_CARLO) {
  log_line("Starting Monte Carlo (restarts resume from completed chunks)")
  chunk_folder <- file.path(OUTPUT_DIR,"mc_data",if (IS_FULL) "blocks_full" else "blocks_pilot")
  dgp_count <- length(DGP_SET)
  for (di in seq_along(DGP_SET)) for (ni in seq_along(SAMPLE_SIZES)) {
    dgp <- DGP_SET[di]; n <- SAMPLE_SIZES[ni]
    blocks <- split(seq_len(MC_REPS), ceiling(seq_len(MC_REPS) / MC_CHUNK_REPS))
    for (block in blocks) {
      block_file <- file.path(chunk_folder, sprintf("%s_n%03d_rep%04d_%04d.rds",
                                                   dgp,n,min(block),max(block)))
      if (file.exists(block_file)) next
      chunk_results <- list(); chunk_errors <- list(); idx <- 1L
      for (r in block) {
        # Design replication seeds remain fixed even if jobs are restarted.
        seed_base <- BASE_SEED + ((di-1)*length(SAMPLE_SIZES)+ni)*100000L + r*100L
        set.seed(seed_base)
        sim <- simulate_data(n,dgp)
        for (ti in seq_along(TAUS)) {
          tau <- TAUS[ti]
          true_b <- sim$truth(tau)
          result <- tryCatch(run_rep_tau(sim$y,sim$X,tau,true_b,
                                         mc_seed=seed_base+ti*10L,
                                         M=MC_M,burn=MC_BURN,
                                         include_nid=TRUE, compute_modes=TRUE),
                             error=function(e)e)
          if (inherits(result,"error")) {
            chunk_errors[[length(chunk_errors)+1L]] <- data.frame(
              DGP=dgp,N=n,Replication=r,Tau=tau,Error=conditionMessage(result))
          } else {
            tab <- result$rows
            tab$DGP <- dgp;tab$N <- n;tab$Replication <- r;tab$Tau <- tau
            tab$Bias <- tab$Estimate-tab$True_Beta
            tab$Squared_Error <- tab$Bias^2
            tab$Mode_Bias <- tab$Mode-tab$True_Beta
            tab$Mode_Squared_Error <- tab$Mode_Bias^2
            chunk_results[[idx]] <- tab; idx <- idx+1L
          }
        }
      }
      saveRDS(list(data=dplyr::bind_rows(chunk_results),
                   errors=dplyr::bind_rows(chunk_errors)),block_file)
      log_line("MC complete ",dgp," n=",n," rep ",min(block),"–",max(block),
               " failures=",length(chunk_errors))
    }
  }
  block_files <- list.files(chunk_folder, pattern="\\.rds$",full.names=TRUE)
  all_blocks <- lapply(block_files,readRDS)
  mc <- dplyr::bind_rows(lapply(all_blocks, `[[`, "data"))
  fail <- dplyr::bind_rows(lapply(all_blocks, `[[`, "errors"))
  if (nrow(fail)>0) write_csv(fail,file.path(OUTPUT_DIR,"logs",paste0("Failures_",MODE_TAG,".csv")))
  if (nrow(mc)==0) stop("No successful Monte Carlo results. See execution log")
  saveRDS(mc,mc_path)
  write_csv(mc,file.path(OUTPUT_DIR,"mc_data",paste0("MC_",MODE_TAG,"_AllRows.csv")))
} else {
  if (!file.exists(mc_path)) stop("No MC cache available at ",mc_path)
  mc <- readRDS(mc_path)
  log_line("Loaded cached MC rows: ", nrow(mc))
}

# ====================== MC TABLES AND FIGURES ===============================
mc_coef <- mc[mc$Variable != "Intercept", ,drop=FALSE]
mc_core <- mc_coef %>%
  dplyr::group_by(DGP,N,Tau,Method,Variable) %>%
  dplyr::summarise(Replications = dplyr::n_distinct(Replication),
                   Bias=mean(Bias,na.rm=TRUE),
                   RMSE=sqrt(mean(Squared_Error,na.rm=TRUE)),
                   Mode_Bias=if(all(is.na(Mode_Bias))) NA_real_ else mean(Mode_Bias,na.rm=TRUE),
                   Mode_RMSE=if(all(is.na(Mode_Squared_Error))) NA_real_ else sqrt(mean(Mode_Squared_Error,na.rm=TRUE)),
                   Coverage=mean(Covered,na.rm=TRUE),
                   Coverage_Error=abs(Coverage-.95),
                   Mean_Width=mean(Interval_Width,na.rm=TRUE),
                   Median_Width=stats::median(Interval_Width,na.rm=TRUE),
                   Mean_Eta=if(all(is.na(Eta))) NA_real_ else mean(Eta,na.rm=TRUE),
                   Mean_ESS=if(all(is.na(ESS))) NA_real_ else mean(ESS,na.rm=TRUE),
                   Mean_Acceptance=if(all(is.na(Acceptance))) NA_real_ else mean(Acceptance,na.rm=TRUE),
                   Mean_Time_Seconds=if(all(is.na(Time_Seconds))) NA_real_ else mean(Time_Seconds,na.rm=TRUE),
                   .groups="drop")
# Aggregate unrounded estimates before preparing display tables.
mc_aggregate <- mc_coef %>%
  dplyr::group_by(DGP,N,Method) %>%
  dplyr::summarise(Replications=dplyr::n_distinct(Replication),
                   Bias=mean(Bias,na.rm=TRUE),
                   RMSE=sqrt(mean(Squared_Error,na.rm=TRUE)),
                   Coverage=mean(Covered,na.rm=TRUE),
                   Coverage_Error=abs(Coverage-.95),
                   Mean_Width=mean(Interval_Width,na.rm=TRUE),
                   Mean_Mode_RMSE=if(all(is.na(Mode_Squared_Error))) NA_real_ else sqrt(mean(Mode_Squared_Error,na.rm=TRUE)),
                   Mean_Time_Seconds=mean(Time_Seconds,na.rm=TRUE),
                   .groups="drop")
mc_diagnostic <- mc %>% dplyr::filter(Method %in% METHODS_POSTERIOR) %>%
  dplyr::group_by(DGP,N,Tau,Method) %>%
  dplyr::summarise(Mean_Diag_Discrepancy=mean(Diagonal_RelRMSE,na.rm=TRUE),
                   Mean_Full_Discrepancy=mean(Full_Frobenius_Gap,na.rm=TRUE),
                   Mean_Offdiag_Discrepancy=mean(OffDiagonal_Gap,na.rm=TRUE),
                   Mean_MCMC_Laplace_Discrepancy=if(all(is.na(MCMC_to_Laplace_Gap))) NA_real_ else mean(MCMC_to_Laplace_Gap,na.rm=TRUE),
                   Mean_ESS=mean(ESS,na.rm=TRUE),
                   Acceptance=mean(Acceptance,na.rm=TRUE), .groups="drop")
mc_joint <- mc %>% dplyr::distinct(DGP,N,Replication,Tau,Method,Joint_Covered) %>%
  dplyr::group_by(DGP,N,Tau,Method) %>%
  dplyr::summarise(Replications=sum(is.finite(Joint_Covered)),
                   Joint_Coverage=if(all(is.na(Joint_Covered))) NA_real_ else mean(Joint_Covered,na.rm=TRUE),
                   .groups="drop")
mc_design <- expand.grid(DGP=DGP_SET,N=SAMPLE_SIZES,Tau=TAUS) %>%
  dplyr::mutate(Replication_Target=MC_REPS, MCMC_Iterations=MC_M,
                Burnin=MC_BURN, Estimators="QR-HAC; QR-NID; fixed/adaptive × normal/diffuse")
mc_tables <- list(Design=mc_design, Aggregate_Results=ROUND4(mc_aggregate),
                  By_Quantile_Coefficient=ROUND4(mc_core),
                  Covariance_Diagnostics=ROUND4(mc_diagnostic),
                  Joint_Regions=ROUND4(mc_joint))
write_xlsx(mc_tables,file.path(OUTPUT_DIR,"paper_tables",paste0("MC_Main_Tables_",MODE_TAG,".xlsx")))
write_csv(ROUND4(mc_core),file.path(OUTPUT_DIR,"appendix_tables",paste0("MC_Quantile_Coefficient_",MODE_TAG,".csv")))
write_csv(ROUND4(mc_joint),file.path(OUTPUT_DIR,"appendix_tables",paste0("MC_Joint_Coverage_",MODE_TAG,".csv")))
# Preserve the completeness record for every design cell.
completeness <- mc %>%
  dplyr::group_by(DGP,N,Tau,Method) %>%
  dplyr::summarise(Replications_Available=dplyr::n_distinct(Replication),.groups="drop") %>%
  dplyr::mutate(Expected=MC_REPS,Complete=(Replications_Available==Expected))
write_csv(completeness,file.path(OUTPUT_DIR,"logs",paste0("MC_Completeness_",MODE_TAG,".csv")))
if(any(!completeness$Complete))
  log_line("ATTENTION: Incomplete Monte Carlo method cells. See MC_Completeness CSV.")
p_cov <- ggplot2::ggplot(mc_aggregate,ggplot2::aes(x=factor(N),y=Coverage,colour=Method,group=Method))+
  ggplot2::geom_hline(yintercept=.95, linetype="dashed",linewidth=.9)+
  ggplot2::geom_line(linewidth=1.1)+ggplot2::geom_point(size=3)+
  ggplot2::facet_wrap(~DGP,ncol=2)+
  ggplot2::scale_colour_manual(values=COLORS)+
  ggplot2::labs(x="Sample size",y="Empirical 95% coverage",colour="Method",
                title="Coverage relative to nominal 95%")+plot_theme()
save_figure(p_cov,file.path(OUTPUT_DIR,"paper_figures",paste0("Figure_MC_Coverage_",MODE_TAG,"_600dpi.png")),14,10)
p_width <- ggplot2::ggplot(mc_aggregate,ggplot2::aes(x=factor(N),y=Mean_Width,colour=Method,group=Method))+
  ggplot2::geom_line(linewidth=1.1)+ggplot2::geom_point(size=3)+
  ggplot2::facet_wrap(~DGP,scales="free_y",ncol=2)+
  ggplot2::scale_colour_manual(values=COLORS)+
  ggplot2::labs(x="Sample size",y="Mean 95% interval width",colour="Method",
                title="Interval width across methods")+plot_theme()
save_figure(p_width,file.path(OUTPUT_DIR,"paper_figures",paste0("Figure_MC_Widths_",MODE_TAG,"_600dpi.png")),14,10)
p_rmse <- ggplot2::ggplot(mc_aggregate,ggplot2::aes(x=factor(N),y=RMSE,colour=Method,group=Method))+
  ggplot2::geom_line(linewidth=1.1)+ggplot2::geom_point(size=3)+
  ggplot2::facet_wrap(~DGP,scales="free_y",ncol=2)+
  ggplot2::scale_colour_manual(values=COLORS)+
  ggplot2::labs(x="Sample size",y="Root mean squared error",colour="Method",
                title="Estimation accuracy: posterior mean and QR")+plot_theme()
save_figure(p_rmse,file.path(OUTPUT_DIR,"paper_figures",paste0("Figure_MC_RMSE_",MODE_TAG,"_600dpi.png")),14,10)
mc_modes_plot <- mc_aggregate %>%
  dplyr::filter(Method %in% METHODS_POSTERIOR) %>%
  tidyr::pivot_longer(c(RMSE,Mean_Mode_RMSE),names_to="Point_Estimator",values_to="Error")
p_modes <- ggplot2::ggplot(mc_modes_plot,ggplot2::aes(x=factor(N),y=Error,colour=Point_Estimator,group=Point_Estimator))+
  ggplot2::geom_line(linewidth=1)+ggplot2::geom_point(size=2.8)+
  ggplot2::facet_grid(DGP~Method,scales="free_y")+
  ggplot2::scale_colour_manual(values=c(RMSE="#006D77",Mean_Mode_RMSE="#E76F51"),
                                labels=c(RMSE="Posterior mean", Mean_Mode_RMSE="Posterior mode"))+
  ggplot2::labs(x="Sample size",y="RMSE",colour="Point estimator",
                title="Prior shrinkage: posterior mean versus mode")+plot_theme()
save_figure(p_modes,file.path(OUTPUT_DIR,"appendix_figures",paste0("Figure_MC_Posterior_Mean_Mode_",MODE_TAG,"_600dpi.png")),19,11)
p_gaps <- ggplot2::ggplot(mc_diagnostic %>% dplyr::filter(Method %in% c("Adaptive-Normal","Adaptive-Diffuse")),
                          ggplot2::aes(x=factor(N),y=Mean_Full_Discrepancy,colour=Method,group=Method))+
  ggplot2::geom_line(linewidth=1)+ggplot2::geom_point(size=3)+
  ggplot2::facet_wrap(~DGP,ncol=2)+ggplot2::scale_colour_manual(values=COLORS)+
  ggplot2::labs(x="Sample size",y="Relative full-matrix discrepancy",colour="Method",
                title="Posterior covariance vs identical HAC target")+plot_theme()
save_figure(p_gaps,file.path(OUTPUT_DIR,"appendix_figures",paste0("Figure_MC_Covariance_Gap_",MODE_TAG,"_600dpi.png")),14,10)
log_line("Main Monte Carlo tables and graphics saved")

# ====================== SENSITIVITY OF HAC / KERNEL DENSITY ==================
# Density and HAC lag sensitivity
sens_file <- file.path(OUTPUT_DIR,"sensitivity",paste0("Sensitivity_",MODE_TAG,".rds"))
if (RUN_SENSITIVITY) {
  if (!file.exists(sens_file)) {
    sens_results <- list(); sens_errors <- list(); ii <- 1L
    sens_designs <- c("Heteroskedastic", "SerialAR")
    for (di in seq_along(sens_designs)) for (n in c(60L,120L)) for (r in seq_len(SENS_REPS)) {
      seed <- BASE_SEED+900000L+di*100000L+n*100L+r
      set.seed(seed)
      dat <- simulate_data(n,sens_designs[di])
      for (tau in c(.10,.50,.90)) {
        base_lag <- default_lag(n)
        variants <- list(
          Pooled_Default = list(density="pooled",bw=1,lag=base_lag),
          Weighted_Default = list(density="weighted",bw=1,lag=base_lag),
          Weighted_NarrowBW = list(density="weighted",bw=.75,lag=base_lag),
          Weighted_WideBW = list(density="weighted",bw=1.5,lag=base_lag),
          Weighted_Lag0 = list(density="weighted",bw=1,lag=0),
          Weighted_Lag1 = list(density="weighted",bw=1,lag=1),
          Weighted_Lag5 = list(density="weighted",bw=1,lag=5))
        for (vi in seq_along(variants)) {
          v <- variants[[vi]];vname <- names(variants)[vi]
          ans <- tryCatch({
            bb <- find_qr(dat$y,dat$X,tau)
            hc <- hac_components(bb,dat$y,dat$X,tau,
                                 lag=v$lag,density=v$density,bandwidth_multiplier=v$bw)
            cc <- calibrate_eta(hc,PRIOR_VAR)
            ps <- sample_gibbs(dat$y,dat$X,tau,cc$eta,bb,hc$J,PRIOR_VAR,
                               M=SENS_M,burn=SENS_BURN,seed=seed+vi*10L+round(tau*100),
                               compute_mode=FALSE)
            target <- as.numeric(dat$truth(tau)[colnames(dat$X)])
            data.frame(DGP=sens_designs[di],N=n,Replication=r,Tau=tau,
                       Variant=vname,Density=v$density,Lag=v$lag,BW_Multiplier=v$bw,
                       Eta=cc$eta, Variable=colnames(dat$X),
                       Covered=as.numeric(target>=ps$lower & target<=ps$upper),
                       Width=as.numeric(ps$upper-ps$lower),
                       Cov_Gap=compare_covariances(ps$cov,hc$V_coef)[2])
          }, error=function(e)e)
          if (inherits(ans,"error")) {
            sens_errors[[length(sens_errors)+1L]] <- data.frame(DGP=sens_designs[di],
              N=n,Replication=r,Tau=tau,Variant=vname,Error=conditionMessage(ans))
          } else {
            sens_results[[ii]] <- ans; ii <- ii+1L
          }
        }
      }
      if (r %% 10L == 0) log_line("Sensitivity ",sens_designs[di]," n=",n," r=",r)
    }
    sens_raw <- dplyr::bind_rows(sens_results)
    sens_fail <- dplyr::bind_rows(sens_errors)
    saveRDS(sens_raw,sens_file)
    if(nrow(sens_fail)>0) write_csv(sens_fail,file.path(OUTPUT_DIR,"logs",paste0("Sensitivity_Failures_",MODE_TAG,".csv")))
  } else sens_raw <- readRDS(sens_file)
  if(nrow(sens_raw)>0) {
    sens_tab <- sens_raw %>% dplyr::filter(Variable != "Intercept") %>%
      dplyr::group_by(DGP,N,Tau,Variant,Density,Lag,BW_Multiplier) %>%
      dplyr::summarise(Replications=dplyr::n_distinct(Replication),
                       Coverage=mean(Covered),Mean_Width=mean(Width),
                       Mean_Eta=mean(Eta),Full_Cov_Gap=mean(Cov_Gap),.groups="drop")
    write_xlsx(ROUND4(sens_tab),file.path(OUTPUT_DIR,"appendix_tables",paste0("Table_Sensitivity_",MODE_TAG,".xlsx")))
    p_sens <- ggplot2::ggplot(sens_tab,
      ggplot2::aes(x=factor(Tau),y=Coverage,fill=Variant))+
      ggplot2::geom_col(position=ggplot2::position_dodge(width=.85),width=.78)+
      ggplot2::geom_hline(yintercept=.95,linetype="dashed")+
      ggplot2::facet_grid(DGP~N)+
      ggplot2::scale_fill_brewer(palette="Set2")+
      ggplot2::labs(x="Quantile",y="Coverage",fill="Estimator setting",
        title="Sensitivity to density estimation and HAC lag")+plot_theme()
    save_figure(p_sens,file.path(OUTPUT_DIR,"appendix_figures",paste0("Figure_Sensitivity_",MODE_TAG,"_600dpi.png")),18,11)
  }
  log_line("Sensitivity analysis saved")
}

# ====================== EMPIRICAL IN-SAMPLE COMPARISONS ======================
empirical_cache <- file.path(OUTPUT_DIR,"data",paste0("Empirical_",MODE_TAG,".rds"))
if (RUN_EMPIRICAL) {
  if (!file.exists(empirical_cache)) {
    empirical_rows <- list();empirical_summary <- list();empirical_cov <- list()
    for (ti in seq_along(TAUS)) {
      tau <- TAUS[ti]
      bb <- find_qr(y_emp,X_emp,tau)
      hc <- hac_components(bb,y_emp,X_emp,tau)
      ca <- calibrate_eta(hc,PRIOR_VAR)
      cd <- calibrate_eta(hc,DIFFUSE_VAR)
      empirical_cov[[ti]] <- list(Tau=tau, HAC=hc$V_coef)
      summary_rows <- list()
      add_emp <- function(method, est, lower,upper,mode=rep(NA_real_,ncol(X_emp)),
                          eta=NA_real_, accept=NA_real_,ess=rep(NA_real_,ncol(X_emp)),
                          cv=rep(NA_real_,4), postcov=NULL) {
        summary_rows[[length(summary_rows)+1L]] <<- data.frame(
          Tau=tau,Method=method,Variable=colnames(X_emp),
          Estimate=as.numeric(est),Posterior_Mode=as.numeric(mode),
          Lower=as.numeric(lower),Upper=as.numeric(upper),
          Width=as.numeric(upper-lower),Eta=eta,
          Acceptance=accept,ESS=ess,Diag_Gap=cv[1],Full_Gap=cv[2],
          Offdiag_Gap=cv[3],MCMC_Laplace_Gap=cv[4])
      }
      ww <- intervals_wald(bb,hc$V_coef)
      add_emp("QR-Wald-HAC",bb,ww$lower,ww$upper)
      nid <- tryCatch(qr_nid_intervals(y_emp,X_emp,tau),error=function(e)NULL)
      if(!is.null(nid)) add_emp("QR-NID",bb,nid$lower,nid$upper)
      for (j in seq_along(METHODS_POSTERIOR)) {
        method <- METHODS_POSTERIOR[j]
        pv <- if (grepl("Diffuse",method)) DIFFUSE_VAR else PRIOR_VAR
        calibrated <- if (grepl("Diffuse",method)) cd else ca
        eta <- if (grepl("^Fixed",method)) 1 else calibrated$eta
        fit <- sample_gibbs(y_emp,X_emp,tau,eta,bb,hc$J,pv,
                            M=EMP_M,burn=EMP_BURN,seed=BASE_SEED+ti*1000L+j,
                            compute_mode=TRUE)
        saveRDS(fit,file.path(OUTPUT_DIR,"data",sprintf("Empirical_Draws_%s_tau%03d_%s.rds",
                                                        method,round(tau*100),MODE_TAG)))
        if (method == "Adaptive-Normal") {
          trace_df <- as.data.frame(fit$draws)
          trace_df$Iteration <- seq_len(nrow(trace_df))
          trace_df <- tidyr::pivot_longer(trace_df, cols = -Iteration,
                            names_to="Variable", values_to="Draw")
          p_trace <- ggplot2::ggplot(trace_df,ggplot2::aes(Iteration,Draw,colour=Variable))+
            ggplot2::geom_line(linewidth=.3)+
            ggplot2::facet_wrap(~Variable,scales="free_y",ncol=2)+
            ggplot2::scale_colour_brewer(palette="Dark2")+
            ggplot2::labs(x="Retained MCMC iteration",y="Coefficient draw",
                  title=paste("Adaptive Metropolis trace, tau =",tau))+plot_theme()+
            ggplot2::theme(legend.position="none")
          save_figure(p_trace,file.path(OUTPUT_DIR,"appendix_figures",
                  sprintf("Figure_MCMC_Trace_tau%03d_%s_600dpi.png",round(tau*100),MODE_TAG)),14,10)
          p_den <- ggplot2::ggplot(trace_df,ggplot2::aes(Draw,fill=Variable))+
            ggplot2::geom_density(alpha=.6,linewidth=.7)+
            ggplot2::facet_wrap(~Variable,scales="free",ncol=2)+
            ggplot2::scale_fill_brewer(palette="Set2")+
            ggplot2::labs(x="Posterior coefficient",y="Density",
                  title=paste("Gibbs posterior density, tau =",tau))+plot_theme()+
            ggplot2::theme(legend.position="none")
          save_figure(p_den,file.path(OUTPUT_DIR,"appendix_figures",
                  sprintf("Figure_MCMC_Density_tau%03d_%s_600dpi.png",round(tau*100),MODE_TAG)),14,10)
        }
        cv <- compare_covariances(fit$cov,hc$V_coef,
                                  if(grepl("^Adaptive",method)) calibrated$D else NULL)
        add_emp(method,fit$mean,fit$lower,fit$upper,fit$mode,
                eta,fit$acceptance,fit$ess,cv)
        # Save full covariance matrices, not just diagonal discrepancies.
        covariance_table <- as.data.frame(fit$cov)
        names(covariance_table) <- colnames(X_emp)
        covariance_table$Variable <- colnames(X_emp)
        covariance_table <- covariance_table[,c("Variable",colnames(X_emp))]
        write_csv(covariance_table,
                  file.path(OUTPUT_DIR,"appendix_tables",
                    sprintf("Empirical_MCMC_Cov_%s_tau%03d_%s.csv",method,round(tau*100),MODE_TAG)))
      }
      covariance_target_table <- as.data.frame(hc$V_coef)
      names(covariance_target_table) <- colnames(X_emp)
      covariance_target_table$Variable <- colnames(X_emp)
      covariance_target_table <- covariance_target_table[,c("Variable",colnames(X_emp))]
      write_csv(covariance_target_table,
                file.path(OUTPUT_DIR,"appendix_tables",
                  sprintf("Empirical_HAC_Target_tau%03d.csv",round(tau*100))))
      empirical_rows[[ti]] <- dplyr::bind_rows(summary_rows)
      empirical_summary[[ti]] <- data.frame(Tau=tau,
          Normal_Eta=ca$eta,Diffuse_Eta=cd$eta,
          HAC_Lag=hc$lag,Density_Estimator=DENSITY_DEFAULT,
          Kernel_Bandwidth=hc$bandwidth,
          Normal_Objective=ca$objective, Diffuse_Objective=cd$objective,
          Normal_Boundary=ca$near_boundary,Diffuse_Boundary=cd$near_boundary)
      log_line("Empirical τ=",tau," eta_normal=",round(ca$eta,5))
    }
    emp_tab <- dplyr::bind_rows(empirical_rows)
    emp_meta <- dplyr::bind_rows(empirical_summary)
    saveRDS(list(tab=emp_tab,meta=emp_meta),empirical_cache)
  } else {
    tmp <- readRDS(empirical_cache)
    emp_tab <- tmp$tab; emp_meta <- tmp$meta
  }
  write_xlsx(list(Estimates=ROUND4(emp_tab),Learning_Rates=ROUND4(emp_meta)),
             file.path(OUTPUT_DIR,"paper_tables",paste0("Empirical_Coefficients_",MODE_TAG,".xlsx")))
  fitted_emp <- emp_tab %>% dplyr::group_by(Tau,Method) %>%
    dplyr::summarise(Coefficients=list(setNames(Estimate,Variable)),.groups="drop")
  empirical_fit_metrics <- list(); k <- 1L
  for (i in seq_len(nrow(fitted_emp))) {
    b <- fitted_emp$Coefficients[[i]]
    qhat <- as.vector(X_emp %*% b[colnames(X_emp)])
    tau <- fitted_emp$Tau[i]
    empirical_fit_metrics[[k]] <- data.frame(Tau=tau,Method=fitted_emp$Method[i],
      In_Sample_Pinball=mean(check_loss(y_emp-qhat,tau)),
      In_Sample_MAE=mean(abs(y_emp-qhat)),
      In_Sample_RMSE=sqrt(mean((y_emp-qhat)^2)))
    k <- k+1L
  }
  empirical_fit_metrics <- dplyr::bind_rows(empirical_fit_metrics)
  write_xlsx(ROUND4(empirical_fit_metrics),
             file.path(OUTPUT_DIR,"paper_tables",paste0("Empirical_In_Sample_All_Methods_",MODE_TAG,".xlsx")))
  p_coef <- ggplot2::ggplot(emp_tab %>% dplyr::filter(Variable!="Intercept",
                                                      Method %in% c("QR-Wald-HAC","Fixed-Normal","Adaptive-Normal")),
           ggplot2::aes(x=Tau,y=Estimate,colour=Method,group=Method))+
    ggplot2::geom_hline(yintercept=0,linetype="dashed",colour="#555555")+
    ggplot2::geom_line(linewidth=1)+ggplot2::geom_point(size=2.8)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower,ymax=Upper),width=.015,linewidth=.6)+
    ggplot2::facet_wrap(~Variable,scales="free_y",ncol=1)+
    ggplot2::scale_colour_manual(values=COLORS)+
    ggplot2::labs(x="Quantile level",y="Standardized coefficient",colour="Method",
                  title="Syrian inflation: HAC confidence vs Gibbs credible intervals")+plot_theme()
  save_figure(p_coef,file.path(OUTPUT_DIR,"paper_figures",paste0("Figure_E2_Intervals_",MODE_TAG,"_600dpi.png")),12,14)
  p_eta <- ggplot2::ggplot(emp_meta,ggplot2::aes(Tau,Normal_Eta))+
    ggplot2::geom_line(colour="#1B9E77",linewidth=1.3)+
    ggplot2::geom_point(colour="#D95F02",size=3.5)+
    ggplot2::labs(x="Quantile level",y="Calibrated learning rate",
                  title="Quantile-specific learning-rate calibration")+plot_theme()
  save_figure(p_eta,file.path(OUTPUT_DIR,"paper_figures",paste0("Figure_E3_Eta_",MODE_TAG,"_600dpi.png")),10,7)
  p_emp_gap <- ggplot2::ggplot(emp_tab %>% dplyr::filter(Variable=="Intercept",
           Method %in% METHODS_POSTERIOR),
           ggplot2::aes(x=Tau,y=Full_Gap,colour=Method))+
    ggplot2::geom_line(linewidth=1)+ggplot2::geom_point(size=3)+
    ggplot2::scale_colour_manual(values=COLORS)+
    ggplot2::labs(x="Quantile level",y="Relative Frobenius discrepancy",
                  colour="Method",title="Real-data posterior covariance versus HAC")+plot_theme()
  save_figure(p_emp_gap,file.path(OUTPUT_DIR,"appendix_figures",paste0("Figure_E4_Covariance_Discrepancy_",MODE_TAG,"_600dpi.png")),12,8)
  # Empirical comparisons and estimated quantile functions
  emp_width_plot <- emp_tab %>% dplyr::filter(Variable != "Intercept",
      Method %in% c("QR-Wald-HAC","Fixed-Normal","Adaptive-Normal"))
  p_emp_width <- ggplot2::ggplot(emp_width_plot,
      ggplot2::aes(x=factor(Tau),y=Width,fill=Method))+
    ggplot2::geom_col(position="dodge",width=.75)+
    ggplot2::facet_wrap(~Variable,scales="free_y")+
    ggplot2::scale_fill_manual(values=COLORS)+
    ggplot2::labs(x="Quantile level",y="95% interval width",
                   title="Empirical interval widths under matched HAC benchmark")+plot_theme()
  save_figure(p_emp_width,file.path(OUTPUT_DIR,"paper_figures",
      paste0("Figure_E5_Interval_Widths_",MODE_TAG,"_600dpi.png")),13,9)
  fitted_lines <- lapply(TAUS,function(tau) {
    dat <- emp_tab %>% dplyr::filter(Tau==tau,Method=="Adaptive-Normal")
    betas <- setNames(dat$Estimate,dat$Variable)
    data.frame(YEAR=raw_input$YEAR,Quantile=paste0("tau = ",tau),
               Fitted=as.vector(X_emp %*% betas[colnames(X_emp)]))
  })
  fitted_lines <- dplyr::bind_rows(fitted_lines)
  p_emp_fit <- ggplot2::ggplot()+
    ggplot2::geom_point(data=raw_input,ggplot2::aes(YEAR,INF_s),size=2.2,
                        colour="#252525")+
    ggplot2::geom_line(data=fitted_lines,
       ggplot2::aes(YEAR,Fitted,colour=Quantile),linewidth=1.1)+
    ggplot2::scale_colour_brewer(palette="Dark2")+
    ggplot2::labs(x="Year",y="Standardized inflation",colour="Quantile",
                  title="Observed inflation and estimated conditional quantiles")+plot_theme()
  save_figure(p_emp_fit,file.path(OUTPUT_DIR,"paper_figures",
      paste0("Figure_E6_Fitted_Quantiles_",MODE_TAG,"_600dpi.png")),14,9)
  log_line("Empirical in-sample estimates saved")
}

# ====================== ROLLING ORIGIN OUT-OF-SAMPLE =========================
# Expanding-window conditional quantile prediction.
# Models use prior-year outcomes, while test-year regressors are realized values.
oos_cache <- file.path(OUTPUT_DIR,"oos",paste0("Rolling_OOS_",MODE_TAG,".rds"))
if (RUN_ROLLING_VALIDATION) {
  if (!file.exists(oos_cache)) {
    forecasts <- list(); errors_oos <- list(); ii <- 1L
    for (train_end in seq.int(OOS_INITIAL,nrow(raw_input)-1L)) {
      tr <- raw_input[seq_len(train_end),,drop=FALSE]
      te <- raw_input[train_end+1L,,drop=FALSE]
      cols <- c("LEXR","M2GDP","TRDGDP")
      mu <- vapply(tr[,cols,drop=FALSE],mean,0.0)
      sd_x <- vapply(tr[,cols,drop=FALSE],stats::sd,0.0)
      if(any(sd_x <= 1e-12)) stop("Near-constant training predictor")
      mean_y <- mean(tr$INF); sd_y <- stats::sd(tr$INF)
      yt <- (tr$INF-mean_y)/sd_y
      XX <- cbind(Intercept=1,sweep(sweep(as.matrix(tr[,cols]),2,mu,"-"),2,sd_x,"/"))
      XX_test <- c(Intercept=1,(as.numeric(unlist(te[1,cols],use.names=FALSE))-mu)/sd_x)
      names(XX_test) <- colnames(XX)
      for (ti in seq_along(TAUS)) {
        tau <- TAUS[ti]
        result <- tryCatch({
          bb <- find_qr(yt,XX,tau)
          hc <- hac_components(bb,yt,XX,tau)
          ca <- calibrate_eta(hc,PRIOR_VAR)
          cd <- calibrate_eta(hc,DIFFUSE_VAR)
          methods <- list("QR-Wald-HAC" = bb)
          eta_settings <- list(Fixed_Normal=1,Adaptive_Normal=ca$eta,
                                Fixed_Diffuse=1,Adaptive_Diffuse=cd$eta)
          labels <- c(Fixed_Normal="Fixed-Normal",Adaptive_Normal="Adaptive-Normal",
                      Fixed_Diffuse="Fixed-Diffuse",Adaptive_Diffuse="Adaptive-Diffuse")
          for(j in seq_along(eta_settings)) {
            setting <- names(eta_settings)[j]
            pv <- if(grepl("Diffuse",setting)) DIFFUSE_VAR else PRIOR_VAR
            pp <- sample_gibbs(yt,XX,tau,eta_settings[[j]],bb,hc$J,pv,
                    M=OOS_M,burn=OOS_BURN,
                    seed=BASE_SEED+train_end*1000L+ti*100L+j,
                    compute_mode=FALSE)
            methods[[labels[[setting]]]] <- pp$mean
          }
          dplyr::bind_rows(lapply(names(methods),function(method) {
            q_std <- as.numeric(sum(XX_test*methods[[method]]))
            pred <- mean_y+sd_y*q_std
            data.frame(Year=te$YEAR,Tau=tau,Method=method,
               Observed=te$INF,Prediction=pred,
               Pinball=check_loss(te$INF-pred,tau),
               AE=abs(te$INF-pred),SE=(te$INF-pred)^2,
               Training_N=train_end)
          }))
        },error=function(e)e)
        if(inherits(result,"error")) {
          errors_oos[[length(errors_oos)+1L]] <- data.frame(Train_End=train_end,Tau=tau,
                                                            Error=conditionMessage(result))
        } else { forecasts[[ii]] <- result;ii <- ii+1L }
      }
      if(train_end %% 5L == 0) log_line("Rolling OOS: completed training n=",train_end)
    }
    forecast_rows <- dplyr::bind_rows(forecasts)
    saveRDS(forecast_rows,oos_cache)
    if(length(errors_oos)>0) write_csv(dplyr::bind_rows(errors_oos),
                    file.path(OUTPUT_DIR,"logs",paste0("OOS_Failures_",MODE_TAG,".csv")))
  } else forecast_rows <- readRDS(oos_cache)
  if(nrow(forecast_rows)>0) {
    oos_metrics <- forecast_rows %>% dplyr::group_by(Tau,Method) %>%
      dplyr::summarise(N_Forecasts=dplyr::n(),
                       Pinball_Loss=mean(Pinball),MAE=mean(AE),
                       RMSE=sqrt(mean(SE)),.groups="drop")
    oos_total <- forecast_rows %>% dplyr::group_by(Method) %>%
      dplyr::summarise(N_Forecasts=dplyr::n(),
                       Mean_Pinball_Loss=mean(Pinball),MAE=mean(AE),
                       RMSE=sqrt(mean(SE)),.groups="drop")
    write_xlsx(list(By_Quantile=ROUND4(oos_metrics),Overall=ROUND4(oos_total)),
      file.path(OUTPUT_DIR,"paper_tables",paste0("OOS_Rolling_Forecast_",MODE_TAG,".xlsx")))
    write_csv(forecast_rows,file.path(OUTPUT_DIR,"oos",paste0("Rolling_Forecasts_",MODE_TAG,".csv")))
    p_oos <- ggplot2::ggplot(oos_metrics,ggplot2::aes(x=Tau,y=Pinball_Loss,colour=Method,group=Method))+
      ggplot2::geom_line(linewidth=1.1)+ggplot2::geom_point(size=3)+
      ggplot2::scale_colour_manual(values=COLORS)+
      ggplot2::labs(x="Quantile level",y="Out-of-sample pinball loss",
          colour="Method",title="Expanding-window conditional quantile comparison")+plot_theme()
    save_figure(p_oos,file.path(OUTPUT_DIR,"paper_figures",paste0("Figure_OOS_Pinball_",MODE_TAG,"_600dpi.png")),13,8)
    # Quantile crossing, separate from predictive error.
    crossings <- forecast_rows %>% dplyr::group_by(Year,Method) %>%
      dplyr::arrange(Tau,.by_group=TRUE) %>%
      dplyr::summarise(Crossing=as.integer(any(diff(Prediction) < 0)),.groups="drop") %>%
      dplyr::group_by(Method) %>%
      dplyr::summarise(Years=dplyr::n(),Fraction_With_Crossing=mean(Crossing),.groups="drop")
    write_csv(crossings,file.path(OUTPUT_DIR,"appendix_tables",paste0("OOS_Quantile_Crossing_",MODE_TAG,".csv")))
    log_line("Rolling out-of-sample comparison saved")
  }
}

# ====================== EMPIRICAL CALIBRATION SENSITIVITY ====================
# This is inexpensive (Laplace-based calibration only; no repeated MCMC).
emp_sens <- list(); ee <- 1L
for (tau in TAUS) {
  bb <- find_qr(y_emp,X_emp,tau)
  configs <- expand.grid(Density=c("pooled","weighted"),
                         BandwidthMultiplier=c(.75,1,1.5),
                         Lag=unique(c(0L,1L,default_lag(length(y_emp)),5L)))
  for (j in seq_len(nrow(configs))) {
    hc <- hac_components(bb,y_emp,X_emp,tau,density=configs$Density[j],
                         bandwidth_multiplier=configs$BandwidthMultiplier[j],
                         lag=configs$Lag[j])
    ca <- calibrate_eta(hc,PRIOR_VAR)
    cd <- calibrate_eta(hc,DIFFUSE_VAR)
    emp_sens[[ee]] <- data.frame(Tau=tau,Density=configs$Density[j],
      BandwidthMultiplier=configs$BandwidthMultiplier[j],Lag=configs$Lag[j],
      Eta_Normal=ca$eta,Eta_Diffuse=cd$eta,
      Objective_Normal=ca$objective,Objective_Diffuse=cd$objective,
      Normal_Near_Boundary=ca$near_boundary,Diffuse_Near_Boundary=cd$near_boundary)
    ee <- ee+1L
  }
}
write_xlsx(ROUND4(dplyr::bind_rows(emp_sens)),
           file.path(OUTPUT_DIR,"appendix_tables","Empirical_Learning_Rate_Sensitivity.xlsx"))

# Execution environment and settings
writeLines(c(
  paste0("R version: ",R.version.string),
  paste0("Execution mode: ",MODE_TAG),
  paste0("Created: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),
  paste0("Source dataset: ",INPUT_FILE),
  paste0("Source observations: ",nrow(raw_input)),
  paste0("Years: ",min(raw_input$YEAR),"-",max(raw_input$YEAR)),
  paste0("Simulation replications per DGP × n: ",MC_REPS),
  paste0("Core quantiles: ",paste(TAUS,collapse=", ")),
  "Posterior kernel: exponentiated quantile check loss and a Gaussian prior.",
  "Same HAC covariance target for Wald and adaptive Gibbs comparisons.",
  paste0("Default prior variances: ",PRIOR_VAR," and ",DIFFUSE_VAR),
  paste0("Default HAC lag rule: floor[4(n/100)^(2/9)]"),
  "HAC: centered scores, Bartlett kernel.",
  "Density: weighted residual Gaussian kernel; bandwidth 1.06*SD(u)*n^(-1/5).",
  "MCMC: adaptive random-walk Metropolis, target acceptance 0.234 during burn-in.",
  "MCMC post-burn proposal frozen; 95% marginal credible intervals and joint ellipsoids.",
  "Learning rate optimized on prior-aware Laplace covariance diagonal; scalar target.",
  "Prediction: expanding window; realized test-year predictors; train-only normalization.",
  "Image output: PNG 600 DPI; axis tick labels and titles 16 pt."
),file.path(OUTPUT_DIR,"notes",paste0("Reproducibility_Record_",MODE_TAG,".txt")))
saveRDS(sessionInfo(),file.path(OUTPUT_DIR,"logs",paste0("R_Session_Info_",MODE_TAG,".rds")))
log_line("FINISHED. Tables, figures, Monte Carlo, and conditional predictions saved to ", OUTPUT_DIR)
cat("\nAnalysis completed. Output directory:\n", normalizePath(OUTPUT_DIR), "\n")
