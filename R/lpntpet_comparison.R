# Gamma-variate response function:
# h(t) = ((t - tD) / (tP - tD))^alpha * exp(alpha * (1 - (t - tD)/(tP - tD)))
# for t > tD, else 0
#
# This parameterization peaks at t = tP and is commonly used because:
# - onset tD controls delay
# - peak tP controls time-to-peak
# - alpha controls sharpness
# - h(tP)=1, so gamma is interpretable as an amplitude coefficient

# Cumulative trapezoidal integral evaluated at the observed time points.
# The integral at the first observation is defined as zero.
cumtrapz_vec <- function(time, y) {
  stopifnot(length(time) == length(y), length(time) >= 2)
  c(0, cumsum(diff(time) * (y[-length(y)] + y[-1]) / 2))
}

gamma_variate_h <- function(time, tD, tP, alpha) {
  h <- numeric(length(time))
  
  if (tP <= tD) return(h)
  
  idx <- which(time > tD)
  x <- (time[idx] - tD) / (tP - tD)
  
  h[idx] <- (x^alpha) * exp(alpha * (1 - x))
  h[!is.finite(h)] <- 0
  h
}


# Build a sensible catalog of candidate lp-ntPET response functions
# based on delayed gamma-variate curves.
#
# Defaults are chosen to span plausible neurotransmitter responses:
# - onset times from early to late in the scan
# - peaks after onset by a moderate amount
# - several shape parameters
build_lpntpet_basis_catalog <- function(time,
                                        tD_grid = NULL,
                                        tP_grid = NULL,
                                        alpha_grid = c(0.25, 0.5, 1, 2, 4)) {
  tmax <- max(time)
  
  if (is.null(tD_grid)) {
    # candidate onset times across the central part of the scan
    tD_grid <- seq(5, max(5, tmax - 15), length.out = 8)
  }
  
  if (is.null(tP_grid)) {
    # candidate peak times across a broad window
    tP_grid <- seq(10, tmax - 5, length.out = 10)
  }
  
  catalog <- list()
  meta <- list()
  k <- 1
  
  for (tD in tD_grid) {
    for (tP in tP_grid) {
      if (tP <= tD) next
      
      for (alpha in alpha_grid) {
        h <- gamma_variate_h(time, tD = tD, tP = tP, alpha = alpha)
        
        # skip nearly zero basis functions
        if (max(h) < 1e-8) next
        
        catalog[[k]] <- h
        meta[[k]] <- data.frame(
          basis_id = k,
          tD = tD,
          tP = tP,
          alpha = alpha
        )
        k <- k + 1
      }
    }
  }
  
  meta_df <- do.call(rbind, meta)
  
  list(
    basis = catalog,
    meta = meta_df
  )
}



# =========================
# lp-ntPET fitting
# =========================

fit_lpntpet_catalog <- function(time, ct, cr,t_onset, weights = NULL, basis_catalog = NULL) {
  stopifnot(length(time) == length(ct), length(ct) == length(cr))
  
  n <- length(time)
  tmax <- max(time)
  if (is.null(weights)) weights <- rep(1, n)
  
  if (is.null(basis_catalog)) {
    basis_catalog <- build_lpntpet_basis_catalog(time,tD_grid=seq(t_onset, max(5, tmax - 15), length.out = 8),tP_grid=seq(t_onset, tmax - 5, length.out = 10))
    #basis_catalog <- build_lpntpet_basis_catalog(time,tD_grid=seq(t_onset, max(5, tmax - 15), length.out = 20),tP_grid=seq(t_onset, tmax - 5, length.out = 20))
  }
  
  H_list <- basis_catalog$basis
  H_meta <- basis_catalog$meta
  
  if (length(H_list) == 0) {
    stop("No candidate basis functions available.")
  }
  
  # Time integrals common across all bases
  int_cr <- cumtrapz_vec(time, cr)
  int_ct <- cumtrapz_vec(time, ct)
  
  results <- vector("list", length(H_list))
  fitted_list <- vector("list", length(H_list))
  
  for (i in seq_along(H_list)) {
    h_i <- H_list[[i]]
    
    # \int_0^t C_T(\tau) h_i(\tau) d\tau
    int_ct_h <- cumtrapz_vec(time, ct * h_i)
    
    # Regression design according to:
    # C_T(t) = R1 * C_R(t) + k2 * int C_R - k2a * int C_T - gamma * int C_T h_i
    X <- cbind(
      cr,
      int_cr,
      -int_ct,
      -int_ct_h
    )
    
    colnames(X) <- c("R1", "k2", "k2a", "gamma")
    
    fit <- lm.wfit(x = X, y = ct, w = weights)
    beta_hat <- fit$coefficients
    
    fitted <- as.vector(X %*% beta_hat)
    resid <- ct - fitted
    wrss <- sum(weights * resid^2)
    
    p <- ncol(X)
    sigma2 <- wrss / max(1, n - p)
    
    XtWX <- crossprod(X, weights * X)
    vcov_beta <- tryCatch(sigma2 * solve(XtWX),
                          error = function(e) matrix(NA_real_, p, p))
    se_beta <- sqrt(diag(vcov_beta))
    
    results[[i]] <- data.frame(
      basis_id = i,
      R1 = beta_hat[1],
      k2 = beta_hat[2],
      k2a = beta_hat[3],
      gamma = beta_hat[4],
      se_R1 = se_beta[1],
      se_k2 = se_beta[2],
      se_k2a = se_beta[3],
      se_gamma = se_beta[4],
      WRSS = wrss,
      AIC = n * log(wrss / n) + 2 * p,
      BIC = n * log(wrss / n) + log(n) * p
    )
    
    fitted_list[[i]] <- data.frame(
      time = time,
      fitted = fitted,
      residuals = resid,
      h = h_i,
      int_ct_h = int_ct_h
    )
  }
  
  prof <- do.call(rbind, results)
  prof <- cbind(prof, H_meta)
  
  best_idx <- which.min(prof$WRSS)
  best <- prof[best_idx, ]
  best_fit <- fitted_list[[best_idx]]
  
  list(
    estimates = best,
    fitted = best_fit$fitted,
    residuals = best_fit$residuals,
    h_best = best_fit$h,
    int_ct_h_best = best_fit$int_ct_h,
    profile = prof,
    basis_catalog = basis_catalog
  )
}

#Functions for SRTM
trapz_conv_exp <- function(time, cr, k2a) {
  # Computes I(t_i; k2a) = integral_0^{t_i} cr(tau) exp(-k2a (t_i - tau)) d tau
  # using trapezoidal integration on the observed time grid.
  
  n <- length(time)
  I <- numeric(n)
  
  for (i in seq_len(n)) {
    if (i == 1) {
      I[i] <- 0
    } else {
      tau <- time[1:i]
      f <- cr[1:i] * exp(-k2a * (time[i] - tau))
      dt <- diff(tau)
      I[i] <- sum((f[-1] + f[-length(f)]) * dt / 2)
    }
  }
  
  I
}


#Fit SRTM
fit_srtm_profile <- function(time, ct, cr, weights = NULL,
                             k2a_grid = seq(0.001, 1, length.out = 300)) {
  # Profile over k2a and estimate R1, k2 by linear regression for each fixed k2a.
  #
  # Args:
  #   time: vector of times
  #   ct: target TAC
  #   cr: reference TAC
  #   weights: optional regression weights
  #   k2a_grid: candidate values for k2a
  #
  # Returns:
  #   list with best estimates and profiling table
  
  stopifnot(length(time) == length(ct), length(ct) == length(cr))
  
  n <- length(time)
  if (is.null(weights)) weights <- rep(1, n)
  
  results <- vector("list", length(k2a_grid))
  
  for (g in seq_along(k2a_grid)) {
    k2a <- k2a_grid[g]
    
    I_k <- trapz_conv_exp(time, cr, k2a)
    X <- cbind(cr, I_k)
    
    # Weighted least squares using lm.wfit
    fit <- lm.wfit(x = X, y = ct, w = weights)
    beta_hat <- fit$coefficients
    
    # Recover SRTM parameters
    R1_hat <- beta_hat[1]
    k2_hat <- beta_hat[2] + R1_hat * k2a
    
    fitted <- as.vector(X %*% beta_hat)
    resid <- ct - fitted
    wrss <- sum(weights * resid^2)
    
    results[[g]] <- data.frame(
      k2a = k2a,
      R1 = R1_hat,
      k2 = k2_hat,
      beta2 = beta_hat[2],
      WRSS = wrss
    )
  }
  
  prof <- do.call(rbind, results)
  best_idx <- which.min(prof$WRSS)
  best <- prof[best_idx, ]
  
  # Refit at optimum to return fitted values etc.
  I_best <- trapz_conv_exp(time, cr, best$k2a)
  X_best <- cbind(cr, I_best)
  fit_best <- lm.wfit(x = X_best, y = ct, w = weights)
  fitted_best <- as.vector(X_best %*% fit_best$coefficients)
  
  list(
    estimates = best,
    fitted = fitted_best,
    residuals = ct - fitted_best,
    profile = prof,
    I_best = I_best
  )
}

