# LGPE-SRTM with a scalar reference-group level and an anchored
# reference-group time-varying deviation.
library(Matrix)

rbf_kernel <- function(sn, sm, theta0, theta1) {
  theta0 * exp(-(sn - sm)^2 / theta1^2)
}

simulate_tac <- function(t_obs, cr_true, r1, k2, k2a, var_prop) {
  dt_vec <- c(t_obs[1], diff(t_obs))
  m <- length(dt_vec)
  omega_dt <- matrix(0, nrow = m, ncol = m)
  diag(omega_dt) <- 1 / dt_vec
  diag(omega_dt[2:m, 1:(m - 1)]) <- -1 / dt_vec[-1]
  
  ct_true <- solve(diag(k2a) + omega_dt) %*%
    (k2 * diag(1, m) + r1 * omega_dt) %*% cr_true
  
  sigma_t_vec <- var_prop * ct_true
  sigma_r_vec <- var_prop * cr_true
  ct_obs <- ct_true + rnorm(m, sd = sqrt(sigma_t_vec))
  cr_obs <- cr_true + rnorm(m, sd = sqrt(sigma_r_vec))
  
  list(
    ct_true = ct_true,
    cr_true = cr_true,
    ct_obs = ct_obs,
    cr_obs = cr_obs
  )
}


# df <- sim_gamma$df
# predvars <- sim_gamma$predvars
# par_in <- sim_gamma$par_in
# estpat <- sim_gamma$estpat
# control <- sim_gamma$control

# df <- DF_Frick_Amygdala
# predvars <- Predvars
# par_in <- Par_in
# estpat <- Estpat
# control <- Control

# df       = sim$df
# predvars = sim$predvars
# par_in   = sim$par_in
# estpat   = sim$estpat
# control  = sim$control
fit_lgpe_srtm <- function(df, predvars = list(), par_in = list(),
                          estpat = NULL, control = list()) {
  merge_settings <- function(defaults, overrides, argument) {
    keys <- names(overrides)
    if (!is.list(overrides) ||
        (length(overrides) && (is.null(keys) || anyNA(keys) ||
                               any(!nzchar(keys)) || anyDuplicated(keys)))) {
      stop(argument, " must be a list with unique names")
    }
    if (!length(overrides)) return(defaults)
    unknown <- setdiff(keys, names(defaults))
    if (length(unknown)) {
      stop("Unknown ", argument, " setting: ", paste(unknown, collapse = ", "))
    }
    # Replace each top-level value as a whole, including subject-level lists.
    defaults[keys] <- overrides
    defaults
  }
  
  control <- merge_settings(
    list(tol = 1e-4, maxiter = 1000, plot_iterations = FALSE, theta0 = TRUE),
    control, "control"
  )
  if (!is.numeric(control$tol) || length(control$tol) != 1L ||
      !is.finite(control$tol) || control$tol <= 0) {
    stop("control$tol must be a positive finite number")
  }
  if (!is.numeric(control$maxiter) || length(control$maxiter) != 1L ||
      !is.finite(control$maxiter) || control$maxiter < 1 ||
      control$maxiter != floor(control$maxiter)) {
    stop("control$maxiter must be a positive integer")
  }
  if (!is.logical(control$plot_iterations) ||
      length(control$plot_iterations) != 1L || is.na(control$plot_iterations)) {
    stop("control$plot_iterations must be TRUE or FALSE")
  }
  if (length(predvars) == 0L || !is.null(names(predvars))) {
    predvars <- merge_settings(list(R1 = NULL, k2 = NULL, k2a = NULL),
                               predvars, "predvars")
  } else if (!is.list(predvars) || length(predvars) != 3L) {
    stop("Supply predictors as a named list or as a list of length three")
  }
  
  if (!all(vapply(predvars, function(x) {
    is.null(x) || (is.character(x) && !anyNA(x))
  }, logical(1)))) {
    stop("predvars entries must be character vectors of column names")
  }
  predvars <- lapply(predvars, function(x) if (length(x)) x else NULL)
  
  required <- c("id", "time", "ct", "cr", "onset_time", unlist(predvars))
  if (!is.data.frame(df) || !nrow(df) || !all(required %in% names(df))) {
    stop("df must contain observations and columns: ",
         paste(unique(required), collapse = ", "))
  }
  if (anyNA(df$id)) stop("Subject IDs must be present")
  numeric_columns <- setdiff(required, "id")
  if (!all(vapply(df[numeric_columns], function(x) {
    is.numeric(x) && all(is.finite(x))
  }, logical(1)))) {
    stop("Time, TAC, onset and predictor columns must be finite numeric values")
  }
  # Preserve subject order for callers supplying subject-level starts.
  df <- df[order(match(df$id, unique(df$id)), df$time), , drop = FALSE]
  
  
  safe_logdet <- function(A) {
    d <- determinant(A, logarithm = TRUE)
    
    if (d$sign <= 0 || !is.finite(d$modulus)) {
      return(Inf)
    }
    
    as.numeric(d$modulus)
  }
  
  make_gp_cov <- function(grid, theta0, theta1, theta_nu) {
    q <- length(grid)
    K <- matrix(0, nrow = q, ncol = q)
    if (q > 1) {
      for (l1 in 1:(q - 1)) {
        for (l2 in (l1 + 1):q) {
          K[l1, l2] <- K[l2, l1] <-
            rbf_kernel(grid[l1], grid[l2], theta0, theta1)
        }
      }
    }
    diag(K) <- rbf_kernel(0, 0, theta0, theta1) + theta_nu
    K
  }
  
  make_anchored_gp_cov <- function(grid, anchor, theta0, theta1,
                                   theta_nu) {
    q <- length(grid)
    K <- matrix(0, nrow = q, ncol = q)
    for (l1 in seq_len(q)) {
      for (l2 in seq_len(q)) {
        K[l1, l2] <-
          rbf_kernel(grid[l1], grid[l2], theta0, theta1) -
          rbf_kernel(grid[l1], anchor, theta0, theta1) -
          rbf_kernel(anchor, grid[l2], theta0, theta1) +
          rbf_kernel(anchor, anchor, theta0, theta1)
      }
    }
    diag(K) <- diag(K) + theta_nu
    K
  }
  
  p_r1 <- length(predvars[[1]]) + 1
  p_k2 <- length(predvars[[2]]) + 1
  p_k2a <- length(predvars[[3]]) + 1
  
  #p_beta <- length(predvars) #Number of predictors
  ids <- unique(df$id)  #All id:S 
  n <- length(ids)   #Sample size
  t_all <- unique(df$time)   #The common grid (in time)
  m_vec <- c() #The number of unique time points in the sample 
  ts_all_tot <- NULL
  for(j in 1:n){
    df_j <- df[df$id==ids[j],]
    t_j <- df_j$time
    if (length(t_j) < 2L || t_j[1] <= 0 || any(diff(t_j) <= 0)) {
      stop("Subject ", ids[j], " must have at least two distinct positive times")
    }
    if (length(unique(df_j$onset_time)) != 1L) {
      stop("onset_time must be constant within subject ", ids[j])
    }
    onset_time_j <- df_j$onset_time[1]
    ts_j <- t_j - onset_time_j   #The domain variable defining variations in k_2a
    ts_j[which(ts_j<0)] <- 0
    m_vec[j] <- length(t_j)
    ts_all_tot <- c(ts_all_tot,ts_j)
  }
  ts_all <- unique(sort(ts_all_tot))
  n_all <- length(ts_all_tot)
  if (is.null(estpat)) {
    estpat <- list(kernel_par = list(beta_theta0 = matrix(1, 1L, p_k2a)))
  }
  estpat_beta_theta0 <- estpat$kernel_par$beta_theta0
  if (!is.matrix(estpat_beta_theta0) || nrow(estpat_beta_theta0) == 0L ||
      ncol(estpat_beta_theta0) != p_k2a ||
      anyNA(estpat_beta_theta0) ||
      any(!estpat_beta_theta0 %in% c(0, 1)) ||
      any(colSums(estpat_beta_theta0) != 1L) ||
      any(rowSums(estpat_beta_theta0) == 0L)) {
    stop("estpat must assign each k2a term to one kernel amplitude")
  }
  n_restrictions_kernel <- dim(estpat_beta_theta0)[1]
  
  m <- length(ts_all) #Number of unique points in the domain
  if (m < 2) stop("At least two distinct points in ts_all are required")
  m_ref_dev <- m - 1
  
  # k2a fixed effects are parameterised as:
  #   one scalar reference-group level;
  #   m - 1 reference-group deviations, anchored to zero at ts_all[1];
  #   m coefficients for each additional k2a predictor.
  # The total remains p_k2a * m.
  
  p_beta <- p_r1 + p_k2 + p_k2a*m   #The number of betas
  p_alpha <- 2 + m    #The number of alphas
  
  #############################
  #Collect parameter values
  ############################
  par_in <- merge_settings(
    list(
      R1 = rep(1, n),
      k2 = rep(0.5, n),
      k2a = lapply(m_vec, function(frames) rep(0.1, frames)),
      Sigma_alpha_r1 = matrix(0.12^2),
      Sigma_alpha_k2 = matrix(0.10^2),
      Sigma_beta_r1 = diag(1, p_r1),
      Sigma_beta_k2 = diag(1, p_k2),
      Sigma_beta_k2a_intercept = 0.1^2,
      Mu_beta_r1 = matrix(0, p_r1),
      Mu_beta_k2 = matrix(0, p_k2),
      Mu_beta_k2a = matrix(0, m * p_k2a),
      beta_theta0 = rep(1, p_k2a),
      alpha_theta0 = 0.01^2,
      theta1 = 10,
      theta_nu = 1e-9,
      par_beta_theta0 = rep(0.01, n_restrictions_kernel),
      var_prop = 10
    ),
    par_in, "par_in"
  )
  
  #1. Pharmacokinetic parameters (for every subject)
  r1_curr <- par_in$R1
  k2_curr <- par_in$k2
  k2a_curr <- par_in$k2a
  
  #2. Variance parameters
  var_prop_curr <- par_in$var_prop
  sigma_alpha_r1_curr <- par_in$Sigma_alpha_r1
  sigma_alpha_k2_curr <- par_in$Sigma_alpha_k2
  
  #3. Kernel parameters
  beta_theta0_curr <- par_in$beta_theta0  #Can be a vector
  alpha_theta0_curr <- par_in$alpha_theta0   #A scalar
  theta1_curr <- par_in$theta1    #Smoothness parameter
  theta_nu_curr <- par_in$theta_nu
  
  #4a Full Sigma_alpha and Sigma_beta matrices
  sigma_alpha_k2a <- make_gp_cov(ts_all, alpha_theta0_curr,
                                 theta1_curr, theta_nu_curr)
  sigma_alpha_k2a_inv <- solve(sigma_alpha_k2a)
  
  if (is.null(par_in$Sigma_beta_k2a_intercept)) {
    stop("par_in$Sigma_beta_k2a_intercept must specify the prior variance of the scalar k2a reference-group level")
  }
  sigma_beta_k2a_intercept <- as.numeric(par_in$Sigma_beta_k2a_intercept)
  if (length(sigma_beta_k2a_intercept) != 1L ||
      !is.finite(sigma_beta_k2a_intercept) ||
      sigma_beta_k2a_intercept <= 0) {
    stop("par_in$Sigma_beta_k2a_intercept must be one positive finite number")
  }
  
  make_sigma_beta_k2a <- function(beta_theta0, theta1, theta_nu) {
    K <- as.matrix(bdiag(
      sigma_beta_k2a_intercept,
      make_anchored_gp_cov(ts_all[-1], ts_all[1], beta_theta0[1],
                           theta1, theta_nu)
    ))
    if (p_k2a > 1) {
      for (lp in 2:p_k2a) {
        K <- as.matrix(bdiag(
          K,
          make_gp_cov(ts_all, beta_theta0[lp], theta1, theta_nu)
        ))
      }
    }
    K
  }
  
  sigma_beta_k2a <- make_sigma_beta_k2a(beta_theta0_curr,
                                        theta1_curr, theta_nu_curr)
  sigma_beta_k2a_inv <- solve(sigma_beta_k2a)
  
  #sigma_beta <- as.matrix(bdiag(diag(par_in$Sigma_beta_r1),par_in$Sigma_beta_k2,sigma_beta_k2a))
  sigma_beta <- as.matrix(bdiag(par_in$Sigma_beta_r1,par_in$Sigma_beta_k2,sigma_beta_k2a))
  sigma_beta_inv <- solve(sigma_beta)
  
  sigma_alpha <- as.matrix(bdiag(par_in$Sigma_alpha_r1,par_in$Sigma_alpha_k2,sigma_alpha_k2a))
  sigma_alpha_inv <- as.matrix(bdiag(solve(par_in$Sigma_alpha_r1),1/par_in$Sigma_alpha_k2,sigma_alpha_k2a_inv))
  
  #4b. Full Mu_beta vector
  mu_beta <- c(par_in$Mu_beta_r1,par_in$Mu_beta_k2,par_in$Mu_beta_k2a)
  if (length(mu_beta) != p_beta) {
    stop("The combined prior mean must have length p_beta; Mu_beta_k2a must follow the scalar/deviation/group-effect ordering")
  }
  
  
  
  #5. Beta
  beta_curr <- rep(0,p_beta)
  
  #6. Theta0 (beta)  #The number of restricted elements only
  par_beta_theta0_curr <- par_in$par_beta_theta0
  beta_theta0_curr <- c()
  for(k in 1:n_restrictions_kernel){
    index_tmp <- which(estpat_beta_theta0[k,]==1)
    beta_theta0_curr[index_tmp] <- par_beta_theta0_curr[k]
  }
  
  
  
  
  e_list <- list()
  x_list <- list()
  z_list <- list()
  m_vec <- c()
  omega_list <- list()
  y_list <- list()
  x_r1 <- matrix(NA,nrow=n,ncol=p_r1)
  x_k2 <- matrix(NA,nrow=n,ncol=p_k2)
  x_k2a <- matrix(NA,nrow=n,ncol=p_k2a)
  for(j in 1:n){
    df_j <- df[df$id==ids[j],]
    onset_time_j <- df_j$onset_time[1]
    t_j <- df_j$time
    dt_vec <- c(t_j[1],diff(t_j))
    m_j <- length(t_j)
    m_vec[j] <- m_j
    ts_j <- t_j - onset_time_j   #The domain variable defining variations in k_2a
    ts_j[which(ts_j<0)] <- 0
    
    omega_dt <- matrix(0,nrow=m_j,ncol=m_j)
    diag(omega_dt) <- 1/dt_vec
    #diag(omega_dt[2:m_j,1:(m_j-1)]) <- -1/dt_vec[-m_j]
    diag(omega_dt[2:m_j,1:(m_j-1)]) <- -1/dt_vec[-1]
    #omega_dt <- t(omega_dt)
    
    
    omega_list[[j]] <- omega_dt
    
    e_j <- NULL
    for(l in 1:m_j){
      e_j_row_l <- rep(0,m)
      e_j_row_l[which(ts_all %in% ts_j[l])] <- 1
      e_j <- rbind(e_j,e_j_row_l)
    }
    
    
    e_list[[j]] <- e_j
    y_list[[j]] <- omega_dt%*%df_j$ct
    
    #Add the columns in x_j corresponding to R1
    x_list[[j]] <- matrix(0,nrow=m_j,ncol=p_beta)
    z_list[[j]] <- matrix(0,nrow=m_j,ncol=p_alpha)
    
    x_list[[j]][,1] <- omega_dt%*%df_j$cr
    z_list[[j]][,1] <- omega_dt%*%df_j$cr
    if(is.null(predvars[[1]])==FALSE){
      for(k in 2:p_r1){
        x_list[[j]][,k] <- df_j[[predvars[[1]][k-1]]]*omega_dt%*%df_j$cr
      }
    } 
    
    #Add the columns in x_j corresponding to k2
    x_list[[j]][,p_r1+1] <- df_j$cr
    z_list[[j]][,2] <- df_j$cr
    if(is.null(predvars[[2]])==FALSE){
      for(k in 2:p_k2){
        x_list[[j]][,p_r1+k-1] <- df_j[[predvars[[2]][k-1]]]*df_j$cr
      }
    }
    
    #Current column index for x
    lx <- p_r1 + p_k2
    lz <- 2
    
    #Add ct to the x_j matrix
    #x_list[[j]][1:m_j,1:m+lx] <- diag(df_j$ct)%*%e_j
    #z_list[[j]][1:m_j,1:m+lz] <- diag(df_j$ct)%*%e_j
    
    # Scalar reference-group level.
    x_list[[j]][, lx + 1] <- -df_j$ct
    
    # Reference-group time variation, anchored to zero at ts_all[1].
    x_list[[j]][, (lx + 2):(lx + m)] <-
      -diag(df_j$ct) %*% e_j[, -1, drop = FALSE]
    z_list[[j]][, (lz + 1):(lz + m)] <- -diag(df_j$ct) %*% e_j
    
    
    
    if(is.null(predvars[[3]])==FALSE){
      for(k in 2:p_k2a){
        lx_k <- lx + m*(k-1)
        #x_list[[j]][1:m_j,1:m+lx] <- diag(df_j[[predvars[[3]]]][k-1]*df_j$ct)%*%e_j
        x_list[[j]][, (lx_k + 1):(lx_k + m)] <-
          -diag(df_j[[predvars[[3]][k-1]]] * df_j$ct) %*% e_j
      }
    }
    
    #Collect x matrices for r1, k2 and k2a
    x_r1[,1] <- 1
    if(p_r1>1){
      for(j in 1:n){
        df_j <- df[df$id==ids[j],] 
        for(k in 2:p_r1){
          x_r1[j,k] <- df_j[[predvars[[1]][k-1]]][1]
        }
        
        
      }
    }
    
    x_k2[,1] <- 1
    if(p_k2>1){
      for(j in 1:n){
        df_j <- df[df$id==ids[j],] 
        for(k in 2:p_k2){
          x_k2[j,k] <- df_j[[predvars[[2]][k-1]]][1]
        }
        
        
      }
    }
    
    x_k2a[,1] <- 1
    if(p_k2a>1){
      for(j in 1:n){
        df_j <- df[df$id==ids[j],] 
        for(k in 2:p_k2a){
          x_k2a[j,k] <- df_j[[predvars[[3]][k-1]]][1]
        }
        
        
      }
    }
    
  }
  
  
  lx_k <- p_r1 + p_k2
  lz_k <- 2
  woods_index_list <- list()
  woods_index_list[[1]] <- (lx_k + 2):(lx_k + m) #For functional statistics
  k2a_scalar_index <-lx_k+1 #Index of k2a scalar coefficient
  if(is.null(predvars[[3]])==FALSE){
    for(k in 2:p_k2a){
      lx_k <- lx + m*(k-1)
      #x_list[[j]][1:m_j,1:m+lx] <- diag(df_j[[predvars[[3]]]][k-1]*df_j$ct)%*%e_j
      woods_index_list[[k]] <- (lx_k + 1):(lx_k + m)
    }
    
  }
  
  ###################################
  ####Start estimating parameters
  ####################################
  
  y_hat_list <- list()
  
  iter <- 1
  maxabs_dpar <- 1e100
  
  while(iter<control$maxiter & maxabs_dpar>control$tol){
    
    #1. Update beta
    
    v_list <- list()
    v_inv_list <- list()
    xtvx_list <- list()
    xtvy_list <- list()
    sum_xtvx <- matrix(0,nrow=p_beta,ncol=p_beta)
    aa_list <- list()  #for autocorrelation due to ct
    bb_list <- list()  #for autocorrelation due to cr
    sigma_r_list <- list()  #error variance in cr
    sigma_t_list <- list()  #error variance in ct
    sigma_list <- list()  #error variance in the SRTM equation
    sigma_inv_list <- list()
    wr_list <- list()
    wt_list <- list()
    w_list <- list()
    w_inv_list <- list()
    sum_xtvy <- rep(0,p_beta)
    
    for(j in 1:n){
      
      df_j <- df[df$id==ids[j],]
      
      sigma_r_list[[j]] <- var_prop_curr*diag(df_j$wr)
      sigma_t_list[[j]] <- var_prop_curr*diag(df_j$wt)
      
      wr_list[[j]] <- diag(df_j$wr)
      wt_list[[j]] <- diag(df_j$wt)
      
      #aa_list[[j]] <- diag(par_in$k2a[[j]]) + omega_list[[j]]
      aa_list[[j]] <- diag(rep(mean(k2a_curr[[j]]),m_vec[j])) + omega_list[[j]]
      #bb_list[[j]] <- diag(par_in$k2a[[j]]) + par_in$R1[j]*omega_list[[j]]
      bb_list[[j]] <- diag(rep(k2_curr[j],m_vec[j])) + r1_curr[j]*omega_list[[j]]
      #bb_list[[j]] <- r1_curr[j]*omega_list[[j]]
      
      
      
      
      w_list[[j]] <- t(aa_list[[j]])%*%wt_list[[j]]%*%aa_list[[j]] + t(bb_list[[j]])%*%wr_list[[j]]%*%bb_list[[j]]
      w_inv_list[[j]] <- solve(w_list[[j]])
      
      
      
      sigma_list[[j]] <- aa_list[[j]]%*%sigma_t_list[[j]]%*%t(aa_list[[j]]) + bb_list[[j]]%*%sigma_r_list[[j]]%*%t(bb_list[[j]])
      sigma_inv_list[[j]] <- solve(sigma_list[[j]])
      
      v_list[[j]] <-  z_list[[j]]%*%sigma_alpha%*%t(z_list[[j]]) + sigma_list[[j]]
      
      v_inv_list[[j]] <- solve(v_list[[j]])
      
      xtvx_list[[j]] <- t(x_list[[j]])%*%v_inv_list[[j]]%*%x_list[[j]]
      xtvy_list[[j]] <- t(x_list[[j]])%*%v_inv_list[[j]]%*%y_list[[j]]
      
      sum_xtvx <- sum_xtvx + xtvx_list[[j]]
      sum_xtvy <- sum_xtvy + xtvy_list[[j]]
      
    }
    
    # plot(diag(sum_xtvx))
    # plot(diag(sigma_beta_inv))
    # plot(diag(sigma_beta))
    
    #sigma_beta_inv <- solve(as.matrix(sigma_beta))
    lambda <- sigma_beta_inv +  sum_xtvx
    lambda_inv <- solve(lambda)
    beta_new <- lambda_inv%*%(sigma_beta_inv%*%mu_beta + sum_xtvy)    #Add prior mean beta!!!!!!!!!!!!!!!!!!!!
    se_beta <- sqrt(diag(lambda_inv))
    if (isTRUE(control$plot_iterations)) {
      par(mfrow = c(2, 1))
      plot(beta_new[-seq_len(p_r1 + p_k2)],
           main = "k2a fixed effects", ylab = "Estimate")
      plot(se_beta[-seq_len(p_r1 + p_k2)],
           main = "k2a fixed-effect SEs", ylab = "SE")
    }
    se_beta_k2a <- se_beta[-c(1:(p_r1+p_k2))]
    
    #The the difference between new and old beta
    dbeta <- beta_new - beta_curr
    
    beta_curr <- beta_new   #Update current beta
    
    
    beta_r1_curr <- beta_curr[1:p_r1]
    se_beta_r1 <- se_beta[1:p_r1]
    beta_k2_curr <- beta_curr[(p_r1+1):(p_r1+p_k2)]
    se_beta_k2 <- se_beta[(p_r1+1):(p_r1+p_k2)]
    beta_k2a_curr <- beta_curr[-seq_len(p_r1 + p_k2)]
    
    beta_k2a_intercept_curr <- beta_k2a_curr[1]
    beta_k2a_ref_dev_curr <- c(0, beta_k2a_curr[2:m])
    beta_k2a_mat_curr <- matrix(NA_real_, nrow = m, ncol = p_k2a)
    beta_k2a_mat_curr[, 1] <-
      beta_k2a_ref_dev_curr #+       beta_k2a_intercept_curr
    
    if (p_k2a > 1) {
      for (l in 2:p_k2a) {
        idx_l <- (m * (l - 1) + 1):(m * l)
        beta_k2a_mat_curr[, l] <- beta_k2a_curr[idx_l]
      }
    }
    
    
    # Pointwise SEs for the reconstructed functions use the full posterior
    # covariance, including covariance between the scalar and deviations.
    cov_beta_k2a <- lambda_inv[-seq_len(p_r1 + p_k2),
                               -seq_len(p_r1 + p_k2), drop = FALSE]
    se_beta_k2a_mat <- matrix(NA_real_, nrow = m, ncol = p_k2a)
    #ref_map <- cbind(1, rbind(0, diag(m_ref_dev)))
    ref_map <- cbind(0, rbind(0, diag(m_ref_dev)))
    se_beta_k2a_mat[, 1] <- sqrt(diag(
      ref_map %*% cov_beta_k2a[1:m, 1:m, drop = FALSE] %*% t(ref_map)
    ))
    if (p_k2a > 1) {
      for (l in 2:p_k2a) {
        idx_l <- (m * (l - 1) + 1):(m * l)
        se_beta_k2a_mat[, l] <- sqrt(diag(
          cov_beta_k2a[idx_l, idx_l, drop = FALSE]
        ))
      }
    }
    
    
    
    
    
    #2. Update alpha for all subjects
    alpha_new <- matrix(NA,nrow=n,ncol=p_alpha)
    cov_alpha_array <- array(NA,dim=c(p_alpha,p_alpha,n))
    se_alpha <- matrix(NA,nrow=n,ncol=p_alpha)
    for(j in 1:n){
      #alpha_new[j,] <- solve(sigma_alpha_inv + t(z_list[[j]])%*%sigma_inv_list[[j]]%*%z_list[[j]])%*%t(z_list[[j]])%*%sigma_inv_list[[j]]%*%(y_list[[j]] - x_list[[j]]%*%beta_curr)
      cov_alpha_array[,,j] <-solve(sigma_alpha_inv + t(z_list[[j]])%*%sigma_inv_list[[j]]%*%z_list[[j]])
      se_alpha[j,] <- sqrt(diag(cov_alpha_array[,,j]))
      alpha_new[j,] <- cov_alpha_array[,,j]%*%t(z_list[[j]])%*%sigma_inv_list[[j]]%*%(y_list[[j]] - x_list[[j]]%*%beta_curr)
    }
    
    
    alpha_curr <- alpha_new
    alpha_r1_curr <- alpha_curr[,1]
    se_alpha_r1 <- se_alpha[,1]
    alpha_k2_curr <- alpha_curr[,2]
    se_alpha_k2 <- se_alpha[,2]
    alpha_k2a_curr <- alpha_curr[,-c(1:2)]
    se_alpha_k2a <- se_alpha[,-c(1:2)]
    
    #3. Update pharmacokinetic parameters
    
    r1_curr <- x_r1%*%beta_r1_curr + alpha_r1_curr
    se_r1 <- sqrt(x_r1%*%se_beta_r1^2 + se_alpha_r1^2)
    k2_curr <- x_k2%*%beta_k2_curr + alpha_k2_curr
    se_k2 <- sqrt(x_k2%*%se_beta_k2^2 + se_alpha_k2^2)
    k2a_curr <- list()   #This is on the common domain (ts) with m points 
    for(j in 1:n){
      k2a_j <- beta_k2a_mat_curr[,1] + alpha_k2a_curr[j,]
      if(p_k2a>1){
        for(l in 2:p_k2a){
          k2a_j <- k2a_j + x_k2a[j,l]*beta_k2a_mat_curr[,l]
          
        }
      }
      #k2a_curr[[j]] <- e_list[[j]]%*%k2a_j
      k2a_curr[[j]] <- e_list[[j]]%*%k2a_j + beta_k2a_intercept_curr
    }
    
    
    #4. Update noise parameter var_prop
    var_prop_new <- 0
    for(j in 1:n){
      y_j_hat <- x_list[[j]]%*%beta_curr + z_list[[j]]%*%alpha_curr[j,]
      var_prop_new <- var_prop_new + t(y_list[[j]] - y_j_hat)%*%w_inv_list[[j]]%*%(y_list[[j]] - y_j_hat)
      
      y_hat_list[[j]] <- y_j_hat
      
    }
    var_prop_new <- as.numeric(var_prop_new/n_all)
    d_var_prop <- var_prop_new - var_prop_curr
    var_prop_curr <- var_prop_new
    
    #5. Update random effects variances for R1 and k2
    sigma_alpha_r1_new <- sum(alpha_r1_curr^2)/n
    d_sigma_alpha_r1 <- sigma_alpha_r1_new - sigma_alpha_r1_curr
    sigma_alpha_r1_curr <- sigma_alpha_r1_new
    
    sigma_alpha_k2_new <- sum(alpha_k2_curr^2)/n
    d_sigma_alpha_k2 <- sigma_alpha_k2_new - sigma_alpha_k2_curr
    sigma_alpha_k2_curr <- sigma_alpha_k2_new
    
    #6. Update alpha_theta0
    alpha_theta0_new <- mean(diag(cov(alpha_curr[,-c(1:2)])))
    d_alpha_theta0_new <- alpha_theta0_new - alpha_theta0_curr
    alpha_theta0_curr <- alpha_theta0_new
    
    #7. Update kernel parameters
    par_kernel_curr <- list(theta1=theta1_curr,beta_theta0=beta_theta0_curr,alpha_theta0=alpha_theta0_curr,theta_nu=theta_nu_curr)
    
    # par_kernel <- par_kernel_curr
    # sigma_alpha_k2 <- sigma_alpha_k2_curr
    # sigma_alpha_r1 <- sigma_alpha_r1_curr
    logmarglik_fun <- function(par_kernel,sigma_list,mu_beta,y_list,x_list,ts_all,n,m,m_vec,p_r1,p_k2a,sigma_alpha_r1,sigma_alpha_k2){
      
      theta1 <- par_kernel$theta1
      beta_theta0 <- par_kernel$beta_theta0
      alpha_theta0 <- par_kernel$alpha_theta0
      theta_nu <- par_kernel$theta_nu
      
      sigma_alpha_k2a <- make_gp_cov(ts_all, alpha_theta0,
                                     theta1, theta_nu)
      sigma_alpha_k2a_inv <- solve(sigma_alpha_k2a)
      
      sigma_beta_k2a <- make_sigma_beta_k2a(beta_theta0,
                                            theta1, theta_nu)
      sigma_beta_k2a_inv <- solve(sigma_beta_k2a)
      
      sigma_beta <- as.matrix(bdiag(par_in$Sigma_beta_r1,
                                    par_in$Sigma_beta_k2,
                                    sigma_beta_k2a))
      sigma_beta_inv <- solve(sigma_beta)
      
      sigma_alpha <- as.matrix(bdiag(sigma_alpha_r1_curr+theta_nu_curr,sigma_alpha_k2_curr,sigma_alpha_k2a))
      sigma_alpha_inv <- as.matrix(bdiag(1/(sigma_alpha_r1_curr+theta_nu_curr),1/(sigma_alpha_k2_curr+theta_nu_curr),sigma_alpha_k2a_inv))
      
      #create v_list and lambda from, e.g., sigma_alpha and sigma_beta
      v_list <- list()
      v_inv_list <- list()
      lambda <- matrix(0,nrow=p_beta,ncol=p_beta)
      for(j in 1:n){
        v_list[[j]] <- sigma_list[[j]] + z_list[[j]]%*%sigma_alpha%*%t(z_list[[j]])
        v_inv_list[[j]] <- solve(v_list[[j]])
        lambda <- lambda + t(x_list[[j]])%*%v_inv_list[[j]]%*%x_list[[j]]
      }
      lambda <- lambda + sigma_beta_inv
      lambda_inv <- solve(lambda)
      
      sigma_beta_inv_mu <- sigma_beta_inv%*%mu_beta
      
      l1 <- -0.5*log(2*pi)*sum(m_vec)
      
      l2 <- 0
      for(j in 1:n){
        l2 <- l2 + sum(log(eigen(sigma_list[[j]])$values))
      }
      l2 <- -0.5*l2
      
      l3 <- -0.5*sum(log(eigen(sigma_beta)$values))
      
      l4 <- -0.5*sum(log(eigen(lambda)$values))
      
      l5 <- -0.5*n*sum(log(eigen(sigma_alpha)$values))
      
      l6 <- 0
      for(j in 1:n){
        l6 <- l6 + sum(log(eigen(v_inv_list[[j]])$values))
      }
      l6 <- -0.5*l6
      
      
      sum_ytvinvy <- 0
      sum_xtvinvy <- rep(0,p_beta)
      for(j in 1:n){
        sum_ytvinvy <- sum_ytvinvy + t(y_list[[j]])%*%v_inv_list[[j]]%*%y_list[[j]]
        sum_xtvinvy <- sum_xtvinvy + t(x_list[[j]])%*%v_inv_list[[j]]%*%y_list[[j]]
      }
      
      
      mu_sigma_beta_inv_mu <- t(mu_beta)%*%sigma_beta_inv%*%mu_beta
      
      #The last term in l7, the most complicated
      smxvyl <- t(sigma_beta_inv_mu + sum_xtvinvy)%*%lambda_inv%*%(sigma_beta_inv_mu + sum_xtvinvy)
      
      l7 <- -0.5*(sum_ytvinvy + t(mu_beta)%*%sigma_beta_inv%*%mu_beta - smxvyl) 
      
      logmarglik_val <- l1 + l2 + l3 + l4 + l5 + l6 + l7
      return(logmarglik_val)
    }
    
    
    #par_beta_theta0_tmp <- par_beta_theta0_curr
    optfun <- function(par_beta_theta0_tmp){
      beta_theta0_tmp <- c()
      for(k in 1:n_restrictions_kernel){
        index_tmp <- which(estpat_beta_theta0[k,]==1)
        beta_theta0_tmp[index_tmp] <- par_beta_theta0_tmp[k]
      }
      
      par_kernel_tmp <- par_kernel_curr
      par_kernel_tmp$beta_theta0 <- beta_theta0_tmp
      value <- -logmarglik_fun(par_kernel_tmp,sigma_list,mu_beta,y_list,x_list,ts_all,n,m,m_vec,p_r1,p_k2a,sigma_alpha_r1,sigma_alpha_k2) 
      value <- as.numeric(value)
      return(value)
    }
    
    # par_kernel_tmp <- par_kernel_curr
    # par_kernel_tmp$beta_theta0 <- rep(10,2)
    # logmarglik_fun(par_kernel_tmp,sigma_list,mu_beta,y_list,x_list,ts_all,n,m,m_vec,p_r1,p_k2a,sigma_alpha_r1,sigma_alpha_k2)
    # 
    # optfun(0.001)
    # optfun(0.005)
    # optfun(0.008)
    # optfun(0.01)
    # optfun(0.02)
    # optfun(0.03)
    # optfun(0.5)
    # optfun(1)
    # optfun(5)
    
    
    #optobj <- optim(par_beta_theta0_curr,fn=optfun,method="Brent",lower=1e-10,upper=10)
    
    #optfun(par_beta_theta0_curr)
    
    if(control$theta0==TRUE){
      if(length(par_beta_theta0_curr)==1){
        optobj <- optim(par_beta_theta0_curr,fn=optfun,method="Brent",lower=1e-10,upper=10)
      }else{
        optobj <- optim(par_beta_theta0_curr,fn=optfun)
      }
      
      par_beta_theta0_new <- optobj$par
      log_model_evidence <- -optobj$value
    }else{
      par_beta_theta0_new <- par_beta_theta0_curr
      log_model_evidence <- -optfun(par_beta_theta0_curr)
    }
    

    
    # #######Testing
    # if(par_beta_theta0_new<beta_k2a_intercept_curr^2){
    #   par_beta_theta0_new <- beta_k2a_intercept_curr^2
    # }
    # ######END TESTING
    
    
    
    d_par_beta_theta0 <- par_beta_theta0_new - par_beta_theta0_curr
    
    par_beta_theta0_curr <- par_beta_theta0_new
    
    
    #Also update the full beta_theta0 vector
    beta_theta0_tmp <- c()
    for(k in 1:n_restrictions_kernel){
      index_tmp <- which(estpat_beta_theta0[k,]==1)
      beta_theta0_tmp[index_tmp] <- par_beta_theta0_curr[k]
    }
    beta_theta0_curr <- beta_theta0_tmp
    
    #######Testing
    
    par_kernel_curr$beta_theta0 <- beta_theta0_curr
    ###################
    
    
    
    #Upate sigma_alpha and sigma_beta
    sigma_alpha_k2a <- make_gp_cov(ts_all, alpha_theta0_curr,
                                   theta1_curr, theta_nu_curr)
    sigma_alpha_k2a_inv <- solve(sigma_alpha_k2a)
    
    sigma_beta_k2a <- make_sigma_beta_k2a(beta_theta0_curr,
                                          theta1_curr, theta_nu_curr)
    sigma_beta_k2a_inv <- solve(sigma_beta_k2a)
    
    sigma_beta <- as.matrix(bdiag(par_in$Sigma_beta_r1,
                                  par_in$Sigma_beta_k2,
                                  sigma_beta_k2a))
    sigma_beta_inv <- solve(sigma_beta)
    
    sigma_alpha <- as.matrix(bdiag(sigma_alpha_r1_curr+theta_nu_curr,sigma_alpha_k2_curr+theta_nu_curr,sigma_alpha_k2a))
    sigma_alpha_inv <- as.matrix(bdiag(1/(sigma_alpha_r1_curr+theta_nu_curr),1/(sigma_alpha_k2_curr+theta_nu_curr),sigma_alpha_k2a_inv))
    
    
    
    
    d_par <- c(d_var_prop,d_sigma_alpha_r1,d_sigma_alpha_k2,d_par_beta_theta0,dbeta)
    
    
    maxabs_dpar <- max(abs(d_par))
    iter <- iter + 1
  }
  
  
  woods_p_vec <- c()
  woods_df_vec <- c()
  woods_stat_vec <- c()
  woods_scale_vec <- c()
  
  for(k_l in 1:p_k2a){
    idx <- woods_index_list[[k_l]]
    q   <- length(idx)
    
    ## R selects the coefficients belonging to this functional effect
    R <- matrix(0, nrow = q, ncol = p_beta)
    R[cbind(seq_len(q), idx)] <- 1
    
    b <- as.numeric(R %*% beta_curr)
    
    ## ------------------------------------------------------------
    ## 1. Bayesian/posterior covariance used in the Wald statistic
    ## ------------------------------------------------------------
    D <- R %*% lambda_inv %*% t(R)
    D <- (D + t(D)) / 2
    
    ## ------------------------------------------------------------
    ## 2. Frequentist sampling covariance of the penalized estimate
    ##
    ## beta_hat = Lambda^{-1} score
    ## Var(score) = sum_xtvx
    ## ------------------------------------------------------------
    C <- R %*% lambda_inv %*%
      sum_xtvx %*%
      lambda_inv %*% t(R)
    C <- (C + t(C)) / 2
    
    
    ## ------------------------------------------------------------
    ## Stable generalized inverse of D
    ## ------------------------------------------------------------
    ee_D <- eigen(D, symmetric = TRUE)
    
    tol_D <- max(ee_D$values) * 1e-10
    keep_D <- ee_D$values > tol_D
    
    if (!any(keep_D)) {
      
      woods_stat <- 0
      woods_df   <- 0
      woods_p    <- 1
      
    } else {
      
      U <- ee_D$vectors[, keep_D, drop = FALSE]
      d <- ee_D$values[keep_D]
      
      D_mhalf <- U %*%
        diag(1 / sqrt(d), nrow = length(d)) %*%
        t(U)
      
      D_ginv <- U %*%
        diag(1 / d, nrow = length(d)) %*%
        t(U)
      
      ## Observed Wald statistic
      woods_stat <- as.numeric(
        t(b) %*% D_ginv %*% b
      )
      
      
      ## --------------------------------------------------------
      ## Null distribution:
      ##
      ## W ~ sum_j eta_j * ChiSq_1
      ##
      ## eta_j = eigenvalues of
      ## D^{-1/2} C D^{-1/2}
      ## --------------------------------------------------------
      A <- D_mhalf %*% C %*% D_mhalf
      A <- (A + t(A)) / 2
      
      eta <- eigen(A, symmetric = TRUE,
                   only.values = TRUE)$values
      
      ## Remove numerical negative/zero eigenvalues
      tol_eta <- max(abs(eta)) * 1e-10
      eta <- eta[eta > tol_eta]
      
      
      if (length(eta) == 0) {
        
        woods_df    <- 0
        woods_scale <- NA_real_
        woods_p     <- 1
        
      } else {
        
        ## ----------------------------------------------------
        ## Satterthwaite approximation:
        ##
        ## W / scale ~ ChiSq_df
        ## ----------------------------------------------------
        mu_W  <- sum(eta)
        var_W <- 2 * sum(eta^2)
        
        woods_df <-
          2 * mu_W^2 / var_W
        
        woods_scale <-
          var_W / (2 * mu_W)
        
        woods_p <- pchisq(
          woods_stat / woods_scale,
          df = woods_df,
          lower.tail = FALSE
        )
      }
      
      woods_p_vec[k_l] <- woods_p
      woods_df_vec[k_l] <- woods_df
      woods_stat_vec[k_l] <- woods_stat
      woods_scale_vec[k_l] <- woods_scale
      
    }
  }
  
  
  
  
  woods <- data.frame(beta_index=1:p_k2a-1 ,stat=woods_stat_vec,df=woods_df_vec,pval=woods_p_vec,scale=woods_scale_vec)
  
  
  # m_mat <- lambda_inv%*%sum_xtvx
  # 
  # 
  # for(k_l in 1:p_k2a){
  #   p_r_tmp <- length(woods_index_list[[k_l]])
  #   r_tmp <- matrix(0,nrow=p_r_tmp,ncol=p_beta)
  #   r_tmp[1:p_r_tmp,woods_index_list[[k_l]]] <- diag(p_r_tmp)
  #   df_chi2 <- sum(diag(r_tmp%*%m_mat%*%t(r_tmp)))  
  #   beta_tmp <- beta_curr[woods_index_list[[k_l]]]
  #   chi2_val <- t(beta_tmp)%*%solve(r_tmp%*%lambda_inv%*%t(r_tmp))%*%beta_tmp
  #   #chi2_val <- t(beta_tmp)%*%solve(r_tmp%*%cov_beta_k2a%*%t(r_tmp))%*%beta_tmp
  #   
  #   p_val_tmp <- 1-pchisq(chi2_val,df_chi2)
  #   
  #   p_vals_woods[k_l] <- p_val_tmp
  #   df_chi2_woods[k_l] <- df_chi2
  #   chi2_woods[k_l] <- chi2_val
  # }
  
  # woods <- data.frame(beta_index=1:p_k2a-1 ,chi2=chi2_woods,df=df_chi2_woods,pval=p_vals_woods)
  
  
  
  ##########################  
  #Collect results
  ##########################
  res <- list()
  
  #beta CONSIDERING autocorrelation
  res$beta <- beta_curr
  res$se_beta <- se_beta
  res$beta_r1 <- beta_r1_curr
  res$se_beta_r1 <- se_beta_r1
  res$beta_k2 <- beta_k2_curr
  res$se_beta_k2 <- se_beta_k2
  res$beta_k2a <-  beta_k2a_curr
  res$se_beta_k2a <- se_beta_k2a  
  res$beta_k2a_mat <- beta_k2a_mat_curr
  res$se_beta_k2a_mat <- se_beta_k2a_mat
  res$beta_k2a_intercept <- beta_k2a_intercept_curr
  res$se_beta_k2a_intercept <- se_beta_k2a[1]
  res$beta_k2a_ref_dev <- beta_k2a_ref_dev_curr
  res$se_beta_k2a_ref_dev <- c(0, se_beta_k2a[2:m])
  res$cov_beta_k2a <- cov_beta_k2a
  res$idx_beta_k2a_ref_dev <- 2:m
  cov_ref_dev <- cov_beta_k2a[2:m, 2:m, drop = FALSE]
  res$wald_k2a_ref_variation <- as.numeric(
    t(beta_k2a_curr[2:m]) %*% solve(cov_ref_dev,
                                    beta_k2a_curr[2:m])
  )
  res$df_k2a_ref_variation <- m_ref_dev
  res$p_k2a_ref_variation <- pchisq(
    res$wald_k2a_ref_variation,
    df = res$df_k2a_ref_variation,
    lower.tail = FALSE
  )
  
  res$var_prop <- var_prop_curr
  res$sigma_alpha_r1 <- sigma_alpha_r1_curr
  res$sigma_alpha_k2 <- sigma_alpha_k2_curr
  
  res$alpha <- alpha_curr
  res$cov_alpha_array <- cov_alpha_array
  res$se_alpha <- se_alpha
  res$alpha_r1 <- alpha_r1_curr
  res$se_alpha_r1 <- se_alpha_r1
  res$se_alpha_k2 <- se_alpha_k2
  res$se_alpha_k2a <- se_alpha_k2a
  
  res$r1 <- r1_curr
  res$se_r1 <- se_r1
  res$k2 <- k2_curr
  res$se_k2 <- se_k2
  res$k2a <- k2a_curr
  
  res$alpha_theta0 <- alpha_theta0_curr
  res$beta_theta0 <- beta_theta0_curr
  res$theta1 <- theta1_curr
  res$theta_nu <- theta_nu_curr
  
  res$df <- df
  res$n <- n
  res$n_tot <- dim(df)[1]
  res$p_beta <- p_beta
  res$p_alpha <- p_alpha
  
  res$w_list <- w_list
  res$wr_list <- wr_list
  res$wt_list <- wt_list
  
  res$sigma_list <- sigma_list  
  res$sigma_inv_list <- sigma_inv_list
  
  res$sigma_alpha <- sigma_alpha
  res$sigma_beta <- sigma_beta
  res$ts_all <- ts_all
  
  res$log_model_evidence <- log_model_evidence
  res$sum_xtvx <- sum_xtvx
  res$lambda_inv <- lambda_inv
  res$lambda <- lambda
  
  res$y_hat_list <- y_hat_list
  res$omega_list <- omega_list
  res$x_list <- x_list
  res$z_list <- z_list
  res$y_list <- y_list
  
  res$woods <- woods
  
  return(res)
  
  
}
