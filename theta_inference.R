suppressWarnings(suppressMessages({library(MASS)}))
source("00_estimator/drglm_nb.R")

geti <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.numeric(v) else d }

NN    <- as.integer(geti("TH_N", 200000))
MM    <- as.integer(geti("TH_M", 300))
KK    <- as.integer(geti("TH_K", 20))
BIGN  <- as.integer(geti("TH_BIGN", 5000000))
BIGM  <- as.integer(geti("TH_BIGM", 3))
THETA <- 2
BETA  <- c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80)
OUT   <- "results_output"
SEED0 <- 770000

theta_info <- function(y, mu, th, chunk_size = 500000L) {
  n <- length(y); acc <- 0
  for (lo in seq(1L, n, by = chunk_size)) {
    hi <- min(lo + chunk_size - 1L, n)
    yi <- y[lo:hi]; mi <- mu[lo:hi]
    acc <- acc + sum(trigamma(th) - trigamma(yi + th)
                     - mi / (th * (th + mi))
                     + (mi - yi) / (th + mi)^2)
  }
  acc
}

gen <- function(n, seed) {
  set.seed(seed)
  ncov <- length(BETA) - 1L
  X <- matrix(runif(n * ncov), n, ncov); colnames(X) <- paste0("x", seq_len(ncov))
  mu <- exp(as.numeric(cbind(1, X) %*% BETA))
  y <- rnbinom(n, mu = mu, size = THETA)
  list(dat = data.frame(y = y, X),
       f = as.formula(paste("y ~", paste(colnames(X), collapse = " + "))))
}

one_rep <- function(n, k, seed) {
  d <- gen(n, seed)
  g <- suppressWarnings(glm.nb(d$f, data = d$dat))
  capture.output(o <- drglm(d$f, family = "NB", data = d$dat, k = k))
  Xm  <- model.matrix(d$f, data = d$dat)
  mu  <- exp(as.numeric(Xm %*% o$coef))
  J   <- theta_info(d$dat$y, mu, o$theta)
  se  <- if (is.finite(J) && J > 0) 1 / sqrt(J) else NA_real_
  c(theta_dr = o$theta, se_dr = se,
    theta_glm = g$theta, se_glm = g$SE.theta,
    cover = as.numeric(abs(o$theta - THETA) <= qnorm(0.975) * se),
    iters = o$iterations)
}

cat(sprintf("[theta-inf] stage 1: n=%d, K=%d, M=%d\n", NN, KK, MM))
S1 <- matrix(NA_real_, MM, 6,
             dimnames = list(NULL, c("theta_dr", "se_dr", "theta_glm",
                                     "se_glm", "cover", "iters")))
for (m in seq_len(MM)) {
  S1[m, ] <- one_rep(NN, KK, SEED0 + m)
  if (m %% 25 == 0) cat(sprintf("  rep %d/%d  theta=%.5f se=%.5f\n",
                                m, MM, S1[m, "theta_dr"], S1[m, "se_dr"]))
  gc()
}

s1 <- data.frame(
  stage = "coverage study", n = NN, K = KK, M = MM,
  mean_theta = mean(S1[, "theta_dr"]),
  bias = mean(S1[, "theta_dr"]) - THETA,
  se_emp = sd(S1[, "theta_dr"]),
  se_est_mean = mean(S1[, "se_dr"], na.rm = TRUE),
  se_glmnb_mean = mean(S1[, "se_glm"], na.rm = TRUE),
  coverage_95 = mean(S1[, "cover"], na.rm = TRUE),
  max_abs_se_diff_vs_glmnb = max(abs(S1[, "se_dr"] - S1[, "se_glm"]), na.rm = TRUE),
  max_abs_theta_diff_vs_glmnb = max(abs(S1[, "theta_dr"] - S1[, "theta_glm"]), na.rm = TRUE),
  mean_outer_iters = mean(S1[, "iters"]))

cat(sprintf("\n[theta-inf] stage 2: n=%d, M=%d\n", BIGN, BIGM))
S2 <- matrix(NA_real_, BIGM, 6, dimnames = list(NULL, colnames(S1)))
for (m in seq_len(BIGM)) {
  S2[m, ] <- one_rep(BIGN, 50L, SEED0 + 5000 + m)
  cat(sprintf("  rep %d/%d  theta=%.6f  se_dr=%.6f  se_glm=%.6f\n",
              m, BIGM, S2[m, "theta_dr"], S2[m, "se_dr"], S2[m, "se_glm"]))
  gc()
}
s2 <- data.frame(
  stage = "large-n confirmation", n = BIGN, K = 50, M = BIGM,
  mean_theta = mean(S2[, "theta_dr"]),
  bias = mean(S2[, "theta_dr"]) - THETA,
  se_emp = if (BIGM > 1) sd(S2[, "theta_dr"]) else NA_real_,
  se_est_mean = mean(S2[, "se_dr"], na.rm = TRUE),
  se_glmnb_mean = mean(S2[, "se_glm"], na.rm = TRUE),
  coverage_95 = mean(S2[, "cover"], na.rm = TRUE),
  max_abs_se_diff_vs_glmnb = max(abs(S2[, "se_dr"] - S2[, "se_glm"]), na.rm = TRUE),
  max_abs_theta_diff_vs_glmnb = max(abs(S2[, "theta_dr"] - S2[, "theta_glm"]), na.rm = TRUE),
  mean_outer_iters = mean(S2[, "iters"]))

res <- rbind(s1, s2)
dir.create(OUT, showWarnings = FALSE)
write.csv(res, file.path(OUT, "rev_theta_inference.csv"), row.names = FALSE)
saveRDS(list(stage1 = S1, stage2 = S2), file.path(OUT, "rev_theta_inference_raw.rds"))
cat("\n"); print(t(res), quote = FALSE)
cat("\nWrote rev_theta_inference.csv\n")
