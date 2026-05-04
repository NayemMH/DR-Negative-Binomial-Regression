library(MASS)


drglm_residuals <- function(model, type = "response") {
  y   <- model$y
  mu  <- model$fitted.values
  th  <- model$theta
  
  if (type == "response") {
    y - mu
  } else if (type == "pearson") {

    (y - mu) / sqrt(mu + mu^2 / th)
  } else if (type == "deviance") {
   
    sign(y - mu) * sqrt(2 * abs(
      ifelse(y > 0, y * log(y / mu), 0) -
        (y + th) * log((y + th) / (mu + th))
    ))
  } else {
    stop("Unsupported type. Choose from: 'response', 'pearson', 'deviance'.")
  }
}


drglm <- function(formula, family, data, k,
                  tol        = 1e-7,
                  max_iter   = 200,
                  start_beta = NULL,
                  start_theta = NULL) {
  
  if (family != "NB") stop("Unsupported family. Currently supported: 'NB'")
  
  n              <- nrow(data)
  rows_per_chunk <- ceiling(n / k)
  chunk_bounds   <- lapply(seq_len(k), function(i) {
    lo <- (i - 1L) * rows_per_chunk + 1L
    hi <- min(i * rows_per_chunk, n)
    if (lo > hi) return(NULL)
    lo:hi
  })
  chunk_bounds <- Filter(Negate(is.null), chunk_bounds)
  k_eff        <- length(chunk_bounds)
  
  response_name <- all.vars(formula)[1]
  rhs_terms     <- delete.response(terms(formula, data = data))
  y_source      <- data[[response_name]]
  
  get_chunk <- function(i) data[chunk_bounds[[i]], , drop = FALSE]
  x_cache   <- lapply(seq_len(k_eff), function(i) {
    model.matrix(rhs_terms, data = get_chunk(i))
  })
  
  nms <- colnames(x_cache[[1L]])
  p   <- length(nms)
  
  Y_counts <- numeric(0)
  for (i in seq_len(k_eff)) {
    chunk_y      <- y_source[chunk_bounds[[i]]]
    chunk_counts <- table(chunk_y)
    new_levels   <- setdiff(names(chunk_counts), names(Y_counts))
    if (length(new_levels) > 0) Y_counts[new_levels] <- 0
    Y_counts[names(chunk_counts)] <-
      Y_counts[names(chunk_counts)] + as.numeric(chunk_counts)
  }
  y_vals <- as.integer(names(Y_counts))
  y_cnts <- as.integer(Y_counts)
  
  if (!is.null(start_beta) && !is.null(start_theta)) {
    B_cur  <- as.numeric(start_beta)
    th_cur <- as.numeric(start_theta)[1]
  } else {
    pilot_pois <- glm(formula, family = poisson(), data = get_chunk(1L))
    B_cur      <- as.numeric(coef(pilot_pois))
    
    mu_pilot <- fitted(pilot_pois)
    y_sum    <- 0; y_sq_sum <- 0
    for (i in seq_len(k_eff)) {
      chunk_y  <- y_source[chunk_bounds[[i]]]
      y_sum    <- y_sum    + sum(chunk_y)
      y_sq_sum <- y_sq_sum + sum(chunk_y^2)
    }
    y_mean <- y_sum / n
    y_var  <- max((y_sq_sum - n * y_mean^2) / max(n - 1L, 1L), 0)
    th_cur <- mean(mu_pilot)^2 / max(y_var - y_mean, 1e-4)
    th_cur <- max(th_cur, 0.1)
  }
  
  newton_B <- function(th, B) {
    I_tot <- matrix(0, p, p)
    T_tot <- numeric(p)
    for (i in seq_len(k_eff)) {
      X_i   <- x_cache[[i]]
      Y_i   <- y_source[chunk_bounds[[i]]]
      eta_i <- as.numeric(X_i %*% B)
      mu_i  <- exp(eta_i)
      W_i   <- mu_i^2 / (mu_i + mu_i^2 / th)         
      z_i   <- eta_i + (Y_i - mu_i) / mu_i            
      I_tot <- I_tot + crossprod(X_i * sqrt(W_i))     
      T_tot <- T_tot + as.numeric(t(X_i) %*% (W_i * z_i))  
    }
    I_inv <- solve(I_tot)
    list(B = as.numeric(I_inv %*% T_tot), I_inv = I_inv)
  }

  update_theta <- function(B, th) {
    mu_chunks <- lapply(seq_len(k_eff), function(i) {
      as.numeric(exp(x_cache[[i]] %*% B))
    })
    
    score <- function(t) {
      out <- sum(y_cnts * digamma(y_vals + t)) - n * digamma(t)
      for (i in seq_len(k_eff)) {
        Y_i  <- y_source[chunk_bounds[[i]]]
        mu_i <- mu_chunks[[i]]
        out  <- out +
          sum(log(t / (t + mu_i))) +
          sum((mu_i - Y_i) / (t + mu_i))
      }
      out
    }
    
    lo <- 1e-4
    hi <- max(th, 1)
    for (expand in seq_len(60)) {
      if (score(hi) < 0) break
      hi <- hi * 2
    }
    if (score(hi) >= 0) {
      warning("Could not bracket theta root; returning current theta.")
      return(th)
    }
    uniroot(score, interval = c(lo, hi), tol = 1e-8)$root
  }
  
  converged <- FALSE
  dB  <- Inf
  dth <- Inf
  
  for (iter in seq_len(max_iter)) {
    nb       <- newton_B(th_cur, B_cur)
    B_target <- nb$B
    
    step_size <- 1
    repeat {
      B_try  <- B_cur + step_size * (B_target - B_cur)
      th_try <- tryCatch(update_theta(B_try, th_cur), error = function(e) NA_real_)
      if (all(is.finite(B_try)) && is.finite(th_try)) {
        B_new  <- B_try
        th_new <- th_try
        break
      }
      step_size <- step_size / 2
      if (step_size < 2^-12) {
        warning("Step-halving failed; returning current iterate.")
        B_new  <- B_cur
        th_new <- th_cur
        break
      }
    }
    
    dB  <- max(abs(B_new - B_cur))
    dth <- abs(th_new - th_cur)
    B_cur  <- B_new
    th_cur <- th_new
    
    if (dB < tol && dth < tol) {
      converged <- TRUE
      break
    }
  }
  
  if (!converged) {
    warning(sprintf(
      "drglm did not converge in %d iterations (final dB = %.2e, dtheta = %.2e). Increase max_iter or check model specification.",
      max_iter, dB, dth))
  }
  

  final_nb  <- newton_B(th_cur, B_cur)
  I_inv_out <- final_nb$I_inv
  
  se    <- sqrt(diag(I_inv_out))
  Z     <- B_cur / se
  p_val <- 2 * (1 - pnorm(abs(Z)))
  l_ci  <- B_cur - qnorm(0.975) * se
  u_ci  <- B_cur + qnorm(0.975) * se
  
  coef_table <- data.frame(
    "Estimate"             = B_cur,
    "Incidence Rate Ratio" = exp(B_cur),
    "Std. Error"           = se,
    "z value"              = Z,
    "Pr(>|z|)"             = p_val,
    "95% CI Lower"         = l_ci,
    "95% CI Upper"         = u_ci,
    "95% CI"               = paste0("[", round(l_ci, 6), ", ", round(u_ci, 6), "]"),
    check.names = FALSE,
    row.names   = nms
  )
  
  mu_all <- numeric(n)
  for (i in seq_len(k_eff)) {
    mu_all[chunk_bounds[[i]]] <- exp(as.numeric(x_cache[[i]] %*% B_cur))
  }
  
  cat("Converged         :", converged, "\n")
  cat("Iterations        :", iter, "\n")
  cat("Dispersion (theta):", round(th_cur, 6), "\n\n")
  

  structure(
    list(
      Estimates     = coef_table,         
      coef          = B_cur,               
      theta         = th_cur,              
      fitted.values = mu_all,              
      y             = y_source,            
      df.residual   = n - p,               
      vcov          = I_inv_out,           
      converged     = converged,
      iterations    = iter,
      formula       = formula,
      n             = n,
      p             = p
    ),
    class = "drglm"
  )
}
