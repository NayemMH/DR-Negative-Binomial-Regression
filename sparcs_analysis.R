library(data.table)
library(MASS)

load_drglm <- function(path) {
  txt <- readLines(path, warn = FALSE)
  cut1 <- grep("set.seed(", txt, fixed = TRUE)
  cut2 <- grep("Demonstration on simulated data", txt, fixed = TRUE)
  cut <- c(cut1, cut2)
  if (length(cut)) {
    txt <- txt[seq_len(min(cut) - 1L)]
  }
  env <- new.env(parent = globalenv())
  eval(parse(text = txt), envir = env)
  get("drglm", envir = env)
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- "Unknown"
  x[is.na(x)] <- "Unknown"
  x
}

summarise_counts <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75))
  data.frame(
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

make_reflevel <- function(x, ref) {
  x <- factor(x)
  if (ref %in% levels(x)) {
    x <- relevel(x, ref = ref)
  }
  x
}

extract_coef_table <- function(fit_object) {
  out <- as.data.table(fit_object, keep.rownames = "term")
  names(out) <- c("term", "estimate", "irr", "std_error", "z_value", "p_value", "ci_95")
  out
}

extract_glmnb_table <- function(fit_object) {
  sm <- summary(fit_object)$coefficients
  ci <- suppressMessages(confint.default(fit_object))
  ci <- ci[rownames(sm), , drop = FALSE]
  out <- data.table(
    term = rownames(sm),
    estimate = sm[, 1],
    irr = exp(sm[, 1]),
    std_error = sm[, 2],
    z_value = sm[, 3],
    p_value = sm[, 4],
    ci_95 = sprintf("[%.6f, %.6f]", ci[, 1], ci[, 2])
  )
  out
}

build_top_table <- function(coef_dt, n_top = 20L) {
  dt <- copy(coef_dt)
  dt <- dt[term != "(Intercept)"]
  dt[, abs_log_irr := abs(log(irr))]
  dt[order(p_value, -abs_log_irr)][1:min(.N, n_top),
    .(term, irr = round(irr, 3), std_error = round(std_error, 4),
      z_value = round(z_value, 2), p_value = signif(p_value, 3), ci_95)]
}

set.seed(20260423)

input_csv <- "Hospital_Inpatient_Discharges_(SPARCS_De-Identified)__2024_20260419.csv"
drglm_path <- "NB drglm version 2.R"
output_dir <- "results_output"
dir.create(output_dir, showWarnings = FALSE)

drglm_nb <- load_drglm(drglm_path)

keep_cols <- c(
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

dt <- fread(input_csv, select = keep_cols)

dt[, los := suppressWarnings(as.integer(`Length of Stay`))]
dt[`Length of Stay` == "120+", los := 120L]
dt <- dt[!is.na(los) & los > 0]

cat_vars <- setdiff(keep_cols, "Length of Stay")
for (nm in cat_vars) {
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

analysis_dt <- analysis_dt[sample(.N)]

rm(dt)
gc()

desc_all <- as.data.table(summarise_counts(analysis_dt$los))
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

write.csv(desc_all, file.path(output_dir, "los_summary.csv"), row.names = FALSE)
fwrite(desc_by_admission, file.path(output_dir, "los_by_admission.csv"))
fwrite(desc_by_severity, file.path(output_dir, "los_by_severity.csv"))

nb_formula <- los ~ hsa + age_group + gender +
  admission_type + severity + mortality_risk +
  med_surg + ed_indicator

full_time <- system.time({
  full_fit <- MASS::glm.nb(nb_formula, data = as.data.frame(analysis_dt))
})

coef_dt <- extract_glmnb_table(full_fit)
top_dt <- build_top_table(coef_dt, n_top = 25L)

fwrite(coef_dt, file.path(output_dir, "nb_full_coefficients.csv"))
fwrite(top_dt, file.path(output_dir, "nb_top_effects.csv"))
writeLines(capture.output(print(full_time)), file.path(output_dir, "nb_runtime.txt"))
fwrite(data.table(metric = c("theta", "aic"),
                  value = c(full_fit$theta, AIC(full_fit))),
       file.path(output_dir, "nb_model_fit_stats.csv"))

val_n <- min(10000L, nrow(analysis_dt))
val_idx <- sample.int(nrow(analysis_dt), val_n)
val_data <- as.data.frame(analysis_dt[val_idx])

pilot_fit <- MASS::glm.nb(nb_formula, data = val_data)

drglm_ok <- TRUE
v_time <- system.time({
  val_drglm <- tryCatch(
    drglm_nb(
      formula = nb_formula,
      family = "NB",
      data = val_data,
      k = 20,
      tol = 1e-7,
      max_iter = 100,
      start_beta = coef(pilot_fit),
      start_theta = pilot_fit$theta
    ),
    error = function(e) {
      drglm_ok <<- FALSE
      NULL
    }
  )
})

g_time <- system.time({
  val_glmnb <- pilot_fit
})

val_glm_coef <- coef(val_glmnb)
if (drglm_ok) {
  val_drglm_coef <- setNames(val_drglm$Estimate, rownames(val_drglm))
  all_terms <- union(names(val_drglm_coef), names(val_glm_coef))
  max_abs_diff <- max(abs(val_drglm_coef[all_terms] - val_glm_coef[all_terms]), na.rm = TRUE)
} else {
  max_abs_diff <- NA_real_
}
validation_dt <- data.table(
  metric = c("validation_n", "drglm_seconds", "glm_nb_seconds", "max_abs_coef_diff", "drglm_ok"),
  value = c(
    val_n,
    unname(v_time["elapsed"]),
    unname(g_time["elapsed"]),
    max_abs_diff,
    drglm_ok
  )
)
fwrite(validation_dt, file.path(output_dir, "validation_comparison.csv"))

cat("Analysis complete.\n")
cat(sprintf("Observations used: %d\n", nrow(analysis_dt)))
cat(sprintf("LOS mean: %.3f\n", mean(analysis_dt$los)))
cat(sprintf("LOS variance-to-mean ratio: %.3f\n", var(analysis_dt$los) / mean(analysis_dt$los)))
cat(sprintf("Full-data theta: %.6f\n", full_fit$theta))
cat("Top effects table written to results_output/nb_top_effects.csv\n")
