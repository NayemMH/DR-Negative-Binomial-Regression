geti <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.numeric(v) else d }

peak_mb <- function() {
  g <- gc(); cn <- colnames(g); j <- which(cn == "max used")
  if (length(j)) sum(g[, j[1] + 1L]) else sum(g[, ncol(g)])
}

BETA  <- c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80)

run_worker <- function(config, n, theta, seed) {
  gc(reset = TRUE)
  suppressWarnings(suppressMessages(library(MASS)))
  if (config == "cmle") {
    set.seed(seed); ncov <- length(BETA) - 1L
    X  <- matrix(runif(n * ncov), n, ncov); colnames(X) <- paste0("x", seq_len(ncov))
    mu <- exp(as.numeric(cbind(1, X) %*% BETA))
    y  <- rnbinom(n, mu = mu, size = theta)
    dat <- data.frame(y = y, X)
    f <- as.formula(paste("y ~", paste(colnames(X), collapse = " + ")))
    invisible(suppressWarnings(glm.nb(f, data = dat)))
  } else {
    k <- as.integer(sub("^dr:", "", config))
    source("00_estimator/drglm_nb_stream.R")
    dir <- file.path(tempdir(), sprintf("membench_%d", k))
    nb_partition_sim(n = n, k = k, beta = BETA, theta = theta, dir = dir, seed = seed)
    invisible(drglm_nb_stream(dir))
    unlink(dir, recursive = TRUE)
  }
  cat(sprintf("PEAK_MB=%.1f\n", peak_mb()))
}

driver <- function() {
  n     <- as.integer(geti("MEM_N", 3000000))
  theta <- geti("MEM_THETA", 2)
  seed  <- as.integer(geti("MEM_SEED", 1))
  kset  <- as.integer(strsplit(Sys.getenv("MEM_K", "25,50,100"), ",")[[1]])
  configs <- c("cmle", paste0("dr:", kset))
  rows <- list()
  for (cfg in configs) {
    cat(sprintf("Measuring peak memory: %s (n=%d) ...\n", cfg, n))
    out <- system2("Rscript",
                   c("04_computational/memory_benchmark.R", "worker", cfg, n, theta, seed),
                   stdout = TRUE, stderr = FALSE)
    pk  <- suppressWarnings(as.numeric(sub(".*PEAK_MB=", "", grep("PEAK_MB", out, value = TRUE))))
    rows[[cfg]] <- data.frame(
      method = ifelse(cfg == "cmle", "CMLE (glm.nb)", paste0("D&R stream K=", sub("dr:", "", cfg))),
      K = ifelse(cfg == "cmle", 1L, as.integer(sub("dr:", "", cfg))),
      peak_heap_mb = if (length(pk)) pk[1] else NA_real_)
  }
  res <- do.call(rbind, rows); rownames(res) <- NULL
  base <- res$peak_heap_mb[res$method == "CMLE (glm.nb)"]
  res$reduction_vs_cmle_pct <- round(100 * (base - res$peak_heap_mb) / base, 1)
  dir.create("results_output", showWarnings = FALSE)
  write.csv(res, "results_output/mem_benchmark.csv", row.names = FALSE)
  cat("\n"); print(res); cat("\nWrote results_output/mem_benchmark.csv\n")
}

a <- commandArgs(TRUE)
if (length(a) && a[1] == "worker") {
  run_worker(a[2], as.integer(a[3]), as.numeric(a[4]), as.integer(a[5]))
} else {
  driver()
}
