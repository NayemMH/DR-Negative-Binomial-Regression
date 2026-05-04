
  library(data.table)
  library(MASS)


# Configuration

INPUT_CSV <- "Hospital_Inpatient_Discharges_(SPARCS_De-Identified)__2024_20260419.csv"
OUTPUT_DIR <- "results_output"
DRGLM_CORE_PATH <- "NB drglm version 2.R"

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Switches for running sections. The full SPARCS model and the original
RUN_MAIN_MONTE_CARLO <- FALSE
RUN_SMALL_EFFICIENCY_SIMULATION <- FALSE
RUN_SPARCS_FULL_MODEL <- FALSE
RUN_SPARCS_BENCHMARK <- FALSE
RUN_SPARCS_DIAGNOSTICS <- FALSE
RUN_EXTRA_STRESS_SIMULATIONS <- FALSE
RUN_ROOTOGRAM <- FALSE

set.seed(20260423)

---------------------------------------------------------------------
load_drglm <- function(path = DRGLM_CORE_PATH) {
  txt <- readLines(path, warn = FALSE)
  cut1 <- grep("set.seed(", txt, fixed = TRUE)
  cut2 <- grep("Demonstration on simulated data", txt, fixed = TRUE)
  cut <- c(cut1, cut2)
  if (length(cut)) txt <- txt[seq_len(min(cut) - 1L)]
  env <- new.env(parent = globalenv())
  eval(parse(text = txt), envir = env)
  get("drglm", envir = env)
}

drglm_nb <- load_drglm(DRGLM_CORE_PATH)

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- "Unknown"
  x[is.na(x)] <- "Unknown"
  x
}

make_reflevel <- function(x, ref) {
  x <- factor(x)
  if (ref %in% levels(x)) x <- relevel(x, ref = ref)
  x
}

fmt_mb <- function(x) as.numeric(x) / 1024^2

summarise_counts <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75))
  data.table(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    median = unname(q[2]),
    q1 = unname(q[1]),
    q3 = unname(q[3]),
    min = min(x),
    max = max(x),
    variance = var(x),
    variance_to_mean = var(x) / mean(x)
  )
}

extract_glmnb_table <- function(fit_object) {
  sm <- summary(fit_object)$coefficients
  ci <- suppressMessages(confint.default(fit_object))
  ci <- ci[rownames(sm), , drop = FALSE]
  data.table(
    term = rownames(sm),
    estimate = sm[, 1],
    irr = exp(sm[, 1]),
    std_error = sm[, 2],
    z_value = sm[, 3],
    p_value = sm[, 4],
    ci_95 = sprintf("[%.6f, %.6f]", ci[, 1], ci[, 2])
  )
}

build_top_table <- function(coef_dt, n_top = 25L) {
  dt <- copy(coef_dt)
  dt <- dt[term != "(Intercept)"]
  dt[, abs_log_irr := abs(log(irr))]
  dt[order(p_value, -abs_log_irr)][1:min(.N, n_top),
    .(term, irr = round(irr, 3), std_error = round(std_error, 4),
      z_value = round(z_value, 2), p_value = signif(p_value, 3), ci_95)]
}

sparcs_keep_cols <- c(
  "Length of Stay",
  "Health Service Area",
  "Age Group",
  "Gender",
  "Type of Admission",
  "APR Severity of Illness Description",
  "APR Risk of Mortality",
  "APR Medical Surgical Description",
  "Emergency Department Indicator"
)

load_sparcs_analysis_data <- function(input_csv = INPUT_CSV, permute = FALSE) {
  dt <- fread(input_csv, select = sparcs_keep_cols)
  dt[, los := suppressWarnings(as.integer(`Length of Stay`))]
  dt[`Length of Stay` == "120+", los := 120L]
  dt <- dt[!is.na(los) & los > 0]

  for (nm in setdiff(sparcs_keep_cols, "Length of Stay")) {
    set(dt, j = nm, value = clean_chr(dt[[nm]]))
  }

  analysis_dt <- data.table(
    los = dt$los,
    hsa = make_reflevel(dt[["Health Service Area"]], "New York City"),
    age_group = make_reflevel(dt[["Age Group"]], "18-29"),
    gender = make_reflevel(dt[["Gender"]], "F"),
    admission_type = make_reflevel(dt[["Type of Admission"]], "Elective"),
    severity = make_reflevel(dt[["APR Severity of Illness Description"]], "Minor"),
    mortality_risk = make_reflevel(dt[["APR Risk of Mortality"]], "Minor"),
    med_surg = make_reflevel(dt[["APR Medical Surgical Description"]], "Medical"),
    ed_indicator = make_reflevel(dt[["Emergency Department Indicator"]], "N")
  )

  if (permute) analysis_dt <- analysis_dt[sample(.N)]
  analysis_dt
}

nb_formula <- los ~ hsa + age_group + gender +
  admission_type + severity + mortality_risk +
  med_surg + ed_indicator

# Monte Carlo simulation

make_main_mc_data <- function(nobs = 5000000L, theta = 2) {
  X <- replicate(10, runif(nobs))
  colnames(X) <- paste0("x", 1:10)
  beta <- c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80)
  eta <- drop(cbind(1, X) %*% beta)
  mu <- exp(eta)
  y <- rnbinom(nobs, mu = mu, size = theta)
  data.frame(y = y, X)
}

fit_main_mc_once <- function(k = 1L, nobs = 5000000L) {
  dat <- make_main_mc_data(nobs)
  f <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10

  if (k == 1L) {
    fit <- MASS::glm.nb(f, data = dat)
    list(beta = coef(fit), se = sqrt(diag(vcov(fit))), theta = fit$theta)
  } else {
    fit <- drglm_nb(f, family = "NB", data = dat, k = k)
    list(beta = fit$Estimate, se = fit$`Std. Error`, theta = NA_real_)
  }
}

run_main_monte_carlo <- function(m = 100L, nobs = 5000000L, k_values = c(1L, 25L, 50L, 100L)) {
  for (k in k_values) {
    out <- replicate(m, fit_main_mc_once(k = k, nobs = nobs), simplify = FALSE)
    save(out, file = file.path(OUTPUT_DIR, sprintf("main_mc_k%s.RData", k)))
  }
}

# Smaller efficiency simulation used for reproducible timing

make_efficiency_data <- function(n, theta = 2.5) {
  x1 <- round(rnorm(n, 50, 10))
  x2 <- round(rnorm(n, 7.5, 2.1))
  x3 <- factor(sample(c("0", "1"), n, TRUE))
  x4 <- factor(sample(c("0", "1", "2"), n, TRUE))
  x6 <- round(rnorm(n, 60, 5))
  eta <- -0.2 + 0.01 * x1 - 0.03 * x2 + 0.15 * (x3 == "1") +
    0.12 * (x4 == "1") - 0.08 * (x4 == "2") + 0.02 * x6
  mu <- exp(eta)
  y <- rnegbin(n, mu = mu, theta = theta)
  data.frame(pred_1 = x1, pred_2 = x2, pred_3 = x3, pred_4 = x4, pred_5 = y, pred_6 = x6)
}

fit_efficiency_once <- function(dat, k) {
  f <- pred_5 ~ pred_1 + pred_2 + pred_3 + pred_4 + pred_6

  gc()
  t_dr <- system.time(fit_dr <- drglm_nb(f, family = "NB", data = dat, k = k, max_iter = 100))

  gc()
  t_glm <- system.time(fit_glm <- glm.nb(f, data = dat))

  coef_dr <- setNames(fit_dr$Estimate, rownames(fit_dr))
  coef_glm <- coef(fit_glm)
  all_terms <- union(names(coef_dr), names(coef_glm))

  data.table(
    n = nrow(dat),
    k = k,
    drglm_seconds = unname(t_dr["elapsed"]),
    glmnb_seconds = unname(t_glm["elapsed"]),
    time_reduction_pct = 100 * (unname(t_glm["elapsed"]) - unname(t_dr["elapsed"])) / unname(t_glm["elapsed"]),
    speedup_factor = unname(t_glm["elapsed"]) / unname(t_dr["elapsed"]),
    drglm_object_mb = fmt_mb(object.size(fit_dr)),
    glmnb_object_mb = fmt_mb(object.size(fit_glm)),
    size_reduction_pct = 100 * (fmt_mb(object.size(fit_glm)) - fmt_mb(object.size(fit_dr))) / fmt_mb(object.size(fit_glm)),
    max_abs_coef_diff = max(abs(coef_dr[all_terms] - coef_glm[all_terms]), na.rm = TRUE),
    theta_glmnb = fit_glm$theta
  )
}

run_small_efficiency_simulation <- function(n = 100000L, k_values = c(10L, 20L, 40L), reps = 3L) {
  results <- rbindlist(lapply(k_values, function(k) {
    rbindlist(lapply(seq_len(reps), function(rep_id) {
      out <- fit_efficiency_once(make_efficiency_data(n), k)
      out[, replication := rep_id]
      out
    }))
  }))

  summary_dt <- results[, .(
    mean_drglm_seconds = mean(drglm_seconds),
    mean_glmnb_seconds = mean(glmnb_seconds),
    mean_time_reduction_pct = mean(time_reduction_pct),
    mean_speedup_factor = mean(speedup_factor),
    mean_drglm_object_mb = mean(drglm_object_mb),
    mean_glmnb_object_mb = mean(glmnb_object_mb),
    mean_size_reduction_pct = mean(size_reduction_pct),
    mean_max_abs_coef_diff = mean(max_abs_coef_diff)
  ), by = .(n, k)]

  fwrite(results, file.path(OUTPUT_DIR, "simulation_efficiency_raw.csv"))
  fwrite(summary_dt, file.path(OUTPUT_DIR, "simulation_efficiency_summary.csv"))
  summary_dt
}

# SPARCS descriptive analysis and full NB model

run_sparcs_full_model <- function() {
  analysis_dt <- load_sparcs_analysis_data(permute = TRUE)

  desc_all <- summarise_counts(analysis_dt$los)
  desc_by_admission <- analysis_dt[, .(
    n = .N,
    mean_los = mean(los),
    median_los = median(los),
    sd_los = sd(los)
  ), by = admission_type][order(-n)]
  desc_by_severity <- analysis_dt[, .(
    n = .N,
    mean_los = mean(los),
    median_los = median(los),
    sd_los = sd(los)
  ), by = severity][order(match(severity, c("Minor", "Moderate", "Major", "Extreme", "Undetermined")))]

  fwrite(desc_all, file.path(OUTPUT_DIR, "los_summary.csv"))
  fwrite(desc_by_admission, file.path(OUTPUT_DIR, "los_by_admission.csv"))
  fwrite(desc_by_severity, file.path(OUTPUT_DIR, "los_by_severity.csv"))

  full_time <- system.time(full_fit <- MASS::glm.nb(nb_formula, data = as.data.frame(analysis_dt)))
  coef_dt <- extract_glmnb_table(full_fit)
  top_dt <- build_top_table(coef_dt, n_top = 25L)

  fwrite(coef_dt, file.path(OUTPUT_DIR, "nb_full_coefficients.csv"))
  fwrite(top_dt, file.path(OUTPUT_DIR, "nb_top_effects.csv"))
  writeLines(capture.output(print(full_time)), file.path(OUTPUT_DIR, "nb_runtime.txt"))
  fwrite(data.table(metric = c("theta", "aic"), value = c(full_fit$theta, AIC(full_fit))),
         file.path(OUTPUT_DIR, "nb_model_fit_stats.csv"))

  val_n <- min(200000L, nrow(analysis_dt))
  val_data <- as.data.frame(analysis_dt[sample.int(nrow(analysis_dt), val_n)])
  pilot_fit <- MASS::glm.nb(nb_formula, data = val_data)
  drglm_ok <- TRUE
  v_time <- system.time({
    val_drglm <- tryCatch(
      drglm_nb(nb_formula, family = "NB", data = val_data, k = 20,
               tol = 1e-7, max_iter = 100,
               start_beta = coef(pilot_fit), start_theta = pilot_fit$theta),
      error = function(e) {
        drglm_ok <<- FALSE
        NULL
      }
    )
  })

  if (drglm_ok) {
    val_drglm_coef <- setNames(val_drglm$Estimate, rownames(val_drglm))
    val_glm_coef <- coef(pilot_fit)
    all_terms <- union(names(val_drglm_coef), names(val_glm_coef))
    max_abs_diff <- max(abs(val_drglm_coef[all_terms] - val_glm_coef[all_terms]), na.rm = TRUE)
  } else {
    max_abs_diff <- NA_real_
  }

  validation_dt <- data.table(
    metric = c("validation_n", "drglm_seconds", "glm_nb_seconds", "max_abs_coef_diff", "drglm_ok"),
    value = c(val_n, unname(v_time["elapsed"]), 0, max_abs_diff, drglm_ok)
  )
  fwrite(validation_dt, file.path(OUTPUT_DIR, "validation_comparison.csv"))

  invisible(full_fit)
}

# SPARCS small-sample efficiency benchmark

benchmark_sparcs_one <- function(dat, k = 10L) {
  gc()
  t_glm <- system.time(glm_fit <- MASS::glm.nb(nb_formula, data = dat))
  size_glm <- object.size(glm_fit)

  gc()
  t_dr <- system.time(
    dr_fit <- tryCatch(
      drglm_nb(nb_formula, family = "NB", data = dat, k = k,
               tol = 1e-7, max_iter = 100,
               start_beta = coef(glm_fit), start_theta = glm_fit$theta),
      error = function(e) e
    )
  )

  dr_ok <- !inherits(dr_fit, "error")
  if (dr_ok) {
    size_dr <- object.size(dr_fit)
    dr_coef <- setNames(dr_fit$Estimate, rownames(dr_fit))
    glm_coef <- coef(glm_fit)
    all_terms <- union(names(dr_coef), names(glm_coef))
    max_diff <- max(abs(dr_coef[all_terms] - glm_coef[all_terms]), na.rm = TRUE)
  } else {
    size_dr <- NA_real_
    max_diff <- NA_real_
  }

  data.table(
    n = nrow(dat),
    k = k,
    glm_seconds = unname(t_glm["elapsed"]),
    drglm_seconds = unname(t_dr["elapsed"]),
    glm_object_mb = fmt_mb(size_glm),
    drglm_object_mb = ifelse(is.na(size_dr), NA_real_, fmt_mb(size_dr)),
    time_reduction_pct = ifelse(dr_ok, 100 * (unname(t_glm["elapsed"]) - unname(t_dr["elapsed"])) / unname(t_glm["elapsed"]), NA_real_),
    size_reduction_pct = ifelse(dr_ok, 100 * (fmt_mb(size_glm) - fmt_mb(size_dr)) / fmt_mb(size_glm), NA_real_),
    max_abs_coef_diff = max_diff,
    drglm_ok = dr_ok
  )
}

run_sparcs_benchmark <- function(sample_sizes = c(1000L, 2000L, 5000L)) {
  analysis_dt <- load_sparcs_analysis_data(permute = TRUE)
  bench_dt <- rbindlist(lapply(sample_sizes, function(n) {
    benchmark_sparcs_one(as.data.frame(analysis_dt[1:n]), k = 10L)
  }), fill = TRUE)
  fwrite(bench_dt, file.path(OUTPUT_DIR, "efficiency_benchmark.csv"))
  bench_dt
}

# Post-review SPARCS diagnostics and forest plot

run_sparcs_diagnostics <- function() {
  dt <- fread(INPUT_CSV, select = sparcs_keep_cols)
  initial_n <- nrow(dt)
  dt[, los := suppressWarnings(as.integer(`Length of Stay`))]
  dt[`Length of Stay` == "120+", los := 120L]
  invalid_los_n <- dt[is.na(los) | los <= 0, .N]
  dt <- dt[!is.na(los) & los > 0]
  blank_admission_n <- dt[trimws(as.character(`Type of Admission`)) == "" | is.na(`Type of Admission`), .N]
  dt <- dt[!(trimws(as.character(`Type of Admission`)) == "" | is.na(`Type of Admission`))]

  attrition <- data.table(
    step = c("Raw SPARCS inpatient records", "Removed invalid or missing LOS",
             "Removed blank admission type", "Final analytical sample"),
    removed = c(0L, invalid_los_n, blank_admission_n, 0L),
    remaining = c(initial_n, initial_n - invalid_los_n, nrow(dt), nrow(dt))
  )
  fwrite(attrition, file.path(OUTPUT_DIR, "sample_attrition.csv"))

  for (nm in setdiff(sparcs_keep_cols, "Length of Stay")) {
    set(dt, j = nm, value = clean_chr(dt[[nm]]))
  }

  analysis_dt <- data.table(
    los = dt$los,
    hsa = make_reflevel(dt[["Health Service Area"]], "New York City"),
    age_group = make_reflevel(dt[["Age Group"]], "18-29"),
    gender = make_reflevel(dt[["Gender"]], "F"),
    admission_type = make_reflevel(dt[["Type of Admission"]], "Elective"),
    severity = make_reflevel(dt[["APR Severity of Illness Description"]], "Minor"),
    mortality_risk = make_reflevel(dt[["APR Risk of Mortality"]], "Minor"),
    med_surg = make_reflevel(dt[["APR Medical Surgical Description"]], "Medical"),
    ed_indicator = make_reflevel(dt[["Emergency Department Indicator"]], "N")
  )

  adm_sev <- dcast(analysis_dt[, .N, by = .(admission_type, severity)],
                   admission_type ~ severity, value.var = "N", fill = 0)
  adm_tot <- analysis_dt[, .N, by = admission_type]
  adm_extreme <- analysis_dt[severity == "Extreme", .N, by = admission_type]
  adm_major_extreme <- analysis_dt[severity %chin% c("Major", "Extreme"), .N, by = admission_type]
  adm_summary <- merge(adm_tot, adm_extreme, by = "admission_type", all.x = TRUE,
                       suffixes = c("_total", "_extreme"))
  adm_summary <- merge(adm_summary, adm_major_extreme, by = "admission_type", all.x = TRUE)
  setnames(adm_summary, "N", "N_major_extreme")
  adm_summary[is.na(N_extreme), N_extreme := 0L]
  adm_summary[is.na(N_major_extreme), N_major_extreme := 0L]
  adm_summary[, `:=`(
    pct_extreme = 100 * N_extreme / N_total,
    pct_major_extreme = 100 * N_major_extreme / N_total
  )]
  fwrite(adm_sev, file.path(OUTPUT_DIR, "admission_by_severity_counts.csv"))
  fwrite(adm_summary[order(-N_total)], file.path(OUTPUT_DIR, "admission_severity_summary.csv"))

  make_sparcs_forest_plot()
}

make_sparcs_forest_plot <- function(coef_path = file.path(OUTPUT_DIR, "nb_full_coefficients.csv")) {
  if (!file.exists(coef_path)) return(invisible(NULL))
  coef_dt <- fread(coef_path)
  coef_dt <- coef_dt[term != "(Intercept)"]
  coef_dt[, abs_z := abs(z_value)]
  top <- coef_dt[order(-abs_z)][1:min(15L, .N)]
  top[, term_label := gsub("^hsa", "HSA: ", term)]
  top[, term_label := gsub("^age_group", "Age: ", term_label)]
  top[, term_label := gsub("^gender", "Gender: ", term_label)]
  top[, term_label := gsub("^admission_type", "Admission: ", term_label)]
  top[, term_label := gsub("^severity", "Severity: ", term_label)]
  top[, term_label := gsub("^mortality_risk", "Mortality risk: ", term_label)]
  top[, term_label := gsub("^med_surg", "Med/Surg: ", term_label)]
  top[, term_label := gsub("^ed_indicator", "ED indicator: ", term_label)]
  top[, lower := exp(estimate - 1.96 * std_error)]
  top[, upper := exp(estimate + 1.96 * std_error)]
  top <- top[order(irr)]

  png(file.path(OUTPUT_DIR, "sparcs_irr_forest_top15.png"), width = 2200, height = 1500, res = 220)
  op <- par(mar = c(5, 12, 2, 2))
  y <- seq_len(nrow(top))
  plot(top$irr, y, xlim = range(c(top$lower, top$upper, 1)), yaxt = "n",
       xlab = "Incidence rate ratio (95% Wald CI)", ylab = "",
       pch = 19, log = "x")
  segments(top$lower, y, top$upper, y, lwd = 2)
  abline(v = 1, lty = 2, col = "gray45")
  axis(2, at = y, labels = top$term_label, las = 1, cex.axis = 0.75)
  par(op)
  dev.off()
}

# Additional heavy-overdispersion and sparse simulations

stress_beta <- c(0.25, 0.6, -0.5, 0.35, -0.25, 0.45)
names(stress_beta) <- c("(Intercept)", paste0("x", 1:5))

make_stress_data <- function(n, theta, sparse = FALSE) {
  x <- replicate(5, rnorm(n))
  b <- stress_beta
  if (sparse) b[1] <- -1.3
  eta <- drop(cbind(1, x) %*% b)
  mu <- exp(eta)
  y <- rnbinom(n, mu = mu, size = theta)
  data.frame(y = y, x1 = x[, 1], x2 = x[, 2], x3 = x[, 3], x4 = x[, 4], x5 = x[, 5])
}

run_stress_scenario <- function(label, n, theta, k, sparse = FALSE, m = 30L) {
  f <- y ~ x1 + x2 + x3 + x4 + x5
  coef_rows <- vector("list", m)
  timing_rows <- vector("list", m)
  zero_rate <- numeric(m)

  for (r in seq_len(m)) {
    dat <- make_stress_data(n, theta, sparse)
    zero_rate[r] <- mean(dat$y == 0)

    t_glm <- system.time(fit_glm <- glm.nb(f, data = dat))
    t_dr <- system.time(
      fit_dr <- drglm_nb(f, family = "NB", data = dat, k = k, tol = 1e-7,
                         max_iter = 100, start_beta = coef(fit_glm),
                         start_theta = fit_glm$theta)
    )

    dr_coef <- setNames(fit_dr$Estimate, rownames(fit_dr))
    glm_coef <- coef(fit_glm)
    terms <- names(glm_coef)

    coef_rows[[r]] <- data.table(
      scenario = label, replication = r, term = terms,
      cmle = as.numeric(glm_coef[terms]),
      dr = as.numeric(dr_coef[terms]),
      abs_diff = abs(as.numeric(glm_coef[terms]) - as.numeric(dr_coef[terms]))
    )
    timing_rows[[r]] <- data.table(
      scenario = label, replication = r,
      cmle_seconds = unname(t_glm["elapsed"]),
      dr_seconds = unname(t_dr["elapsed"]),
      theta_cmle = fit_glm$theta,
      theta_dr = NA_real_
    )
  }

  list(coef = rbindlist(coef_rows), timing = rbindlist(timing_rows)[, zero_rate := zero_rate])
}

run_extra_stress_simulations <- function() {
  set.seed(20260503)
  scenarios <- list(
    run_stress_scenario("Heavy overdispersion", n = 100000L, theta = 0.5, k = 20L),
    run_stress_scenario("Sparse small-n", n = 50000L, theta = 0.7, k = 10L, sparse = TRUE)
  )

  coef_all <- rbindlist(lapply(scenarios, `[[`, "coef"))
  timing_all <- rbindlist(lapply(scenarios, `[[`, "timing"))

  coef_summary <- coef_all[, .(
    cmle_mean = mean(cmle),
    dr_mean = mean(dr),
    mean_abs_diff = mean(abs_diff),
    max_abs_diff = max(abs_diff)
  ), by = .(scenario, term)]

  timing_summary <- timing_all[, .(
    n_replications = .N,
    mean_zero_rate = mean(zero_rate),
    cmle_seconds = mean(cmle_seconds),
    dr_seconds = mean(dr_seconds),
    speedup = mean(cmle_seconds) / mean(dr_seconds),
    theta_cmle = mean(theta_cmle),
    theta_dr = mean(theta_dr, na.rm = TRUE),
    theta_abs_diff = mean(abs(theta_cmle - theta_dr), na.rm = TRUE)
  ), by = scenario]
  timing_summary[is.nan(theta_dr), theta_dr := NA_real_]
  timing_summary[is.nan(theta_abs_diff), theta_abs_diff := NA_real_]

  fwrite(coef_all, file.path(OUTPUT_DIR, "extra_simulation_coefficients_raw.csv"))
  fwrite(timing_all, file.path(OUTPUT_DIR, "extra_simulation_timing_raw.csv"))
  fwrite(coef_summary, file.path(OUTPUT_DIR, "extra_simulation_coefficients_summary.csv"))
  fwrite(timing_summary, file.path(OUTPUT_DIR, "extra_simulation_timing_summary.csv"))
  list(coef_summary = coef_summary, timing_summary = timing_summary)
}

# SPARCS hanging rootogram

run_sparcs_rootogram <- function(max_count = 40L) {
  analysis_dt <- load_sparcs_analysis_data(permute = FALSE)
  coef_path <- file.path(OUTPUT_DIR, "nb_full_coefficients.csv")
  fit_path <- file.path(OUTPUT_DIR, "nb_model_fit_stats.csv")
  if (!file.exists(coef_path) || !file.exists(fit_path)) {
    stop("Run the SPARCS full model first, or provide nb_full_coefficients.csv and nb_model_fit_stats.csv.")
  }

  coef_dt <- fread(coef_path)
  beta <- setNames(coef_dt$estimate, coef_dt$term)
  theta <- as.numeric(fread(fit_path)[metric == "theta", value])[1]
  X <- model.matrix(delete.response(terms(nb_formula, data = analysis_dt)), data = analysis_dt)

  missing_terms <- setdiff(colnames(X), names(beta))
  if (length(missing_terms)) {
    warning("Missing coefficient terms set to zero for diagnostic plot: ",
            paste(missing_terms, collapse = ", "))
    beta <- c(beta, setNames(rep(0, length(missing_terms)), missing_terms))
  }
  beta <- beta[colnames(X)]
  mu <- as.numeric(exp(X %*% beta))

  counts <- 1:max_count
  observed <- tabulate(analysis_dt$los, nbins = max_count)[counts]
  expected <- vapply(counts, function(k) sum(dnbinom(k, mu = mu, size = theta)), numeric(1))
  root_dt <- data.table(
    los = counts,
    observed = observed,
    expected = expected,
    hanging = sqrt(observed) - sqrt(expected)
  )
  fwrite(root_dt, file.path(OUTPUT_DIR, "sparcs_rootogram_values.csv"))

  png(file.path(OUTPUT_DIR, "sparcs_nb_hanging_rootogram.png"), width = 2200, height = 1400, res = 220)
  op <- par(mar = c(5, 5, 2, 1))
  plot(root_dt$los, root_dt$hanging, type = "h", lwd = 5, lend = "butt",
       xlab = "Length of stay (days)",
       ylab = expression(sqrt(Observed) - sqrt(Expected)),
       main = "", col = ifelse(root_dt$hanging >= 0, "#2f6f9f", "#b24a3b"))
  abline(h = 0, lty = 2, col = "gray35")
  par(op)
  dev.off()

  summary_dt <- root_dt[, .(
    max_abs_hanging = max(abs(hanging)),
    mean_abs_hanging = mean(abs(hanging)),
    total_observed_1_40 = sum(observed),
    total_expected_1_40 = sum(expected)
  )]
  fwrite(summary_dt, file.path(OUTPUT_DIR, "sparcs_rootogram_summary.csv"))
  summary_dt
}

if (RUN_MAIN_MONTE_CARLO) {
  run_main_monte_carlo()
}

if (RUN_SMALL_EFFICIENCY_SIMULATION) {
  print(run_small_efficiency_simulation())
}

if (RUN_SPARCS_FULL_MODEL) {
  run_sparcs_full_model()
}

if (RUN_SPARCS_BENCHMARK) {
  print(run_sparcs_benchmark())
}

if (RUN_SPARCS_DIAGNOSTICS) {
  run_sparcs_diagnostics()
}

if (RUN_EXTRA_STRESS_SIMULATIONS) {
  print(run_extra_stress_simulations())
}

if (RUN_ROOTOGRAM) {
  print(run_sparcs_rootogram())
}

