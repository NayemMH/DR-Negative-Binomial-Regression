suppressWarnings(suppressMessages({library(MASS); library(parallel)}))

geti <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.numeric(v) else d }
gets <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }

CFG <- list(
  N       = as.integer(geti("SIM_N",     5000000)),
  TOTAL   = as.integer(geti("SIM_TOTAL", 1000)),
  BATCH   = as.integer(geti("SIM_BATCH", 100)),
  THETA   = geti("SIM_THETA", 2),
  SEED    = as.integer(geti("SIM_SEED",  123)),
  CMLE_N  = as.integer(geti("SIM_CMLE_N", 50)),
  KSET    = as.integer(strsplit(gets("SIM_K", "25,50,100"), ",")[[1]]),
  BETA    = c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80),
  OUTDIR  = "results_output",
  BATCHDIR= file.path("results_output", "sim_batches")
)
CFG$BETA_NAMES <- c("(Intercept)", paste0("x", seq_len(length(CFG$BETA) - 1L)))

rep_file <- function(r) file.path(CFG$BATCHDIR, sprintf("rep_%04d.rds", r))
done_reps <- function() {
  fs <- list.files(CFG$BATCHDIR, pattern = "^rep_\\d+\\.rds$")
  if (!length(fs)) return(integer(0))
  as.integer(sub("^rep_(\\d+)\\.rds$", "\\1", fs))
}

run_rep <- function(r, CFG) {
  suppressWarnings(suppressMessages(library(MASS)))
  source("00_estimator/drglm_nb.R")
  source("01_simulation/sgd_nb.R")

  seed <- CFG$SEED + r
  set.seed(seed)
  n <- CFG$N; beta <- CFG$BETA; ncov <- length(beta) - 1L
  X  <- matrix(runif(n * ncov), n, ncov); colnames(X) <- paste0("x", seq_len(ncov))
  mu <- exp(as.numeric(cbind(1, X) %*% beta))
  y  <- rnbinom(n, mu = mu, size = CFG$THETA)
  dat <- data.frame(y = y, X)
  f <- as.formula(paste("y ~", paste(colnames(X), collapse = " + ")))

  res <- list(rep = r, seed = seed, n = n)

  g <- NULL
  if (r <= CFG$CMLE_N) {
    t <- system.time(g <- suppressWarnings(glm.nb(f, data = dat)))["elapsed"]
    res$cmle <- list(beta = coef(g), se = sqrt(diag(vcov(g))), theta = g$theta, time = unname(t))
  }
  ref <- if (!is.null(g)) coef(g) else NULL
  rm(g); gc()

  res$dr <- list()
  for (k in CFG$KSET) {
    t <- system.time(fit <- {
      capture.output(o <- drglm(f, family = "NB", data = dat, k = k)); o
    })["elapsed"]
    b  <- setNames(fit$coef, rownames(fit$Estimates))
    se <- setNames(fit$Estimates[, "Std. Error"], rownames(fit$Estimates))
    res$dr[[as.character(k)]] <- list(
      beta = b, se = se, theta = fit$theta, time = unname(t),
      maxdiff = if (!is.null(ref)) max(abs(b[names(ref)] - ref)) else NA_real_)
    gc()
  }

  t <- system.time({
    idx <- sample.int(n, floor(0.1 * n))
    gs  <- suppressWarnings(glm.nb(f, data = dat[idx, ]))
  })["elapsed"]
  res$sub <- list(beta = coef(gs), se = sqrt(diag(vcov(gs))), theta = gs$theta, time = unname(t))
  rm(gs); gc()

  t <- system.time({
    Xm <- cbind(`(Intercept)` = 1, as.matrix(X))
    sg <- sgd_nb(Xm, y, theta_init = 1, epochs = 8, batch_size = 10000L, seed = seed)
  })["elapsed"]
  res$sgd <- list(beta = sg$coef, theta = sg$theta, time = unname(t))
  rm(Xm, dat, X, y); gc()

  saveRDS(res, rep_file(r))
  drmax <- suppressWarnings(max(vapply(res$dr, function(z) z$maxdiff, 0), na.rm = TRUE))
  sprintf("rep %d done (theta_dr=%.4f, cmle=%s, dr_maxdiff=%s)", r,
          res$dr[[1]]$theta,
          if (is.null(res$cmle)) "skip" else sprintf("%.4f", res$cmle$theta),
          if (is.finite(drmax)) sprintf("%.1e", drmax) else "NA")
}

cmd_status <- function() {
  d <- done_reps()
  cat(sprintf("Simulation progress: %d / %d replications complete.\n", length(d), CFG$TOTAL))
  cat(sprintf("  N per rep = %d | theta = %g | K = {%s}\n",
              CFG$N, CFG$THETA, paste(CFG$KSET, collapse = ", ")))
  cat(sprintf("  Batch dir: %s\n", CFG$BATCHDIR))
  if (length(d) && length(d) < CFG$TOTAL) {
    nxt <- setdiff(seq_len(CFG$TOTAL), d)
    cat(sprintf("  Next reps to run: %s%s\n", paste(head(nxt, 10), collapse = ", "),
                if (length(nxt) > 10) ", ..." else ""))
  }
}

cmd_run <- function(batch_n, workers) {
  dir.create(CFG$BATCHDIR, recursive = TRUE, showWarnings = FALSE)
  todo <- setdiff(seq_len(CFG$TOTAL), done_reps())
  if (!length(todo)) { cat("All replications already complete.\n"); return(invisible()) }
  batch <- head(todo, batch_n)
  cat(sprintf("Running %d replications (%s) on %d worker(s)...\n",
              length(batch), paste(range(batch), collapse = "-"), workers))
  t0 <- Sys.time()
  if (workers <= 1) {
    msgs <- lapply(batch, run_rep, CFG = CFG)
  } else {
    cl <- makeCluster(workers, type = "PSOCK")
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl, c("CFG", "run_rep", "rep_file"), envir = environment())
    msgs <- parLapplyLB(cl, batch, function(r) run_rep(r, CFG))
  }
  for (m in msgs) cat("  ", m, "\n", sep = "")
  cat(sprintf("Batch done in %.1f min. Progress: %d / %d.\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              length(done_reps()), CFG$TOTAL))
}

cmd_reset <- function() {
  if (dir.exists(CFG$BATCHDIR))
    file.remove(list.files(CFG$BATCHDIR, pattern = "^rep_\\d+\\.rds$", full.names = TRUE))
  cat("Cleared all replication files.\n")
}

cmd_combine <- function() {
  d <- sort(done_reps())
  if (!length(d)) { cat("No replications to combine.\n"); return(invisible()) }
  reps <- lapply(d, function(r) readRDS(rep_file(r)))
  M <- length(reps); bn <- CFG$BETA_NAMES; true <- setNames(CFG$BETA, bn)

  pick <- function(v) {
    if (is.null(v)) return(rep(NA_real_, length(bn)))
    if (is.null(names(v)) && length(v) == length(bn)) names(v) <- bn
    v[bn]
  }
  gather <- function(getter) {
    grab <- function(z, field) {
      m <- getter(z)
      if (is.null(m)) pick(NULL) else pick(m[[field]])
    }
    B <- t(vapply(reps, grab, numeric(length(bn)), field = "beta"))
    S <- t(vapply(reps, grab, numeric(length(bn)), field = "se"))
    colnames(B) <- bn
    colnames(S) <- bn
    list(B = B, S = S)
  }

  dr_getter <- function(k) {
    key <- as.character(k)
    function(z) z$dr[[key]]
  }
  methods <- list(CMLE = function(z) z$cmle,
                  Subsample = function(z) z$sub,
                  SGD = function(z) z$sgd)
  for (k in CFG$KSET) methods[[paste0("DR_K", k)]] <- dr_getter(k)

  mean_tab <- data.frame(parameter = bn, true = unname(true))
  for (mn in names(methods)) {
    g <- gather(methods[[mn]])
    mean_tab[[mn]] <- colMeans(g$B, na.rm = TRUE)
  }
  write.csv(mean_tab, file.path(CFG$OUTDIR, "sim_mean_estimates.csv"), row.names = FALSE)

  detail <- list()
  for (mn in names(methods)) {
    g <- gather(methods[[mn]]); B <- g$B; S <- g$S
    bias <- colMeans(B, na.rm = TRUE) - true
    emp  <- apply(B, 2, sd, na.rm = TRUE)
    est  <- colMeans(S, na.rm = TRUE)
    rmse <- sqrt(colMeans(sweep(B, 2, true)^2, na.rm = TRUE))
    cov  <- if (all(is.na(S))) rep(NA_real_, length(bn)) else
      colMeans(abs(sweep(B, 2, true)) <= 1.96 * S, na.rm = TRUE)
    detail[[mn]] <- data.frame(method = mn, parameter = bn, bias = bias,
                               se_emp = emp, se_est = est, rmse = rmse, coverage = cov,
                               row.names = NULL)
  }
  detail <- do.call(rbind, detail)
  write.csv(detail, file.path(CFG$OUTDIR, "sim_detailed_properties.csv"), row.names = FALSE)

  cmle_reps <- which(vapply(reps, function(z) !is.null(z$cmle), logical(1)))
  if (!length(cmle_reps)) cmle_reps <- seq_along(reps)
  mse_of <- function(mn) {
    B <- gather(methods[[mn]])$B[cmle_reps, , drop = FALSE]
    mean(sweep(B, 2, true)^2, na.rm = TRUE)
  }
  mse_cmle <- mse_of("CMLE")
  eff <- data.frame(method = names(methods),
                    mse = vapply(names(methods), mse_of, 0),
                    rel_efficiency = mse_cmle / vapply(names(methods), mse_of, 0),
                    n_reps_used = length(cmle_reps),
                    row.names = NULL)
  write.csv(eff, file.path(CFG$OUTDIR, "sim_relative_efficiency.csv"), row.names = FALSE)

  rep_time <- function(z, getter) {
    v <- getter(z)$time
    if (is.null(v)) NA_real_ else v
  }
  mean_time <- function(getter) {
    mean(vapply(reps, rep_time, 0, getter = getter), na.rm = TRUE)
  }
  tim <- data.frame(
    method = names(methods),
    mean_time_s = vapply(methods, mean_time, 0),
    row.names = NULL)

  max_diff_for_k <- function(k) {
    key <- as.character(k)
    max(vapply(reps, function(z) z$dr[[key]]$maxdiff, 0), na.rm = TRUE)
  }
  drdiff <- vapply(CFG$KSET, max_diff_for_k, 0)
  write.csv(tim, file.path(CFG$OUTDIR, "sim_timing.csv"), row.names = FALSE)
  write.csv(data.frame(K = CFG$KSET, max_abs_coef_diff_vs_CMLE = drdiff),
            file.path(CFG$OUTDIR, "sim_dr_equivalence.csv"), row.names = FALSE)

  writeLines(as.character(M), file.path(CFG$OUTDIR, "sim_n_replications.txt"))
  cat(sprintf("Combined %d replications -> results_output/sim_*.csv\n", M))
  cat(sprintf("  Max D&R-vs-CMLE coef diff across all reps/K: %.2e\n", max(drdiff)))
  cat("  Mean coverage (CMLE):",
      sprintf("%.3f", mean(detail$coverage[detail$method == "CMLE"], na.rm = TRUE)), "\n")
}

main <- function() {
  a <- commandArgs(TRUE)
  cmd <- if (length(a)) a[1] else "status"
  if (cmd == "run") {
    batch_n <- if (length(a) >= 2) as.integer(a[2]) else CFG$BATCH
    workers <- if (length(a) >= 3) as.integer(a[3]) else max(1L, detectCores() - 1L)
    cmd_run(batch_n, workers)
  } else if (cmd == "status")  cmd_status()
    else if (cmd == "combine") cmd_combine()
    else if (cmd == "reset")   cmd_reset()
    else cat("Unknown command. Use: run | status | combine | reset\n")
}
main()
