#' Simulate dynamic PET data under the LGPE-SRTM model
#'
#' Generates synthetic PET time-activity curves (TACs)
#' under the latent Gaussian process extension of the
#' simplified reference tissue model (LGPE-SRTM).
#'
#' @param n Number of subjects.
#' @param t_obs Observation times.
#' @param onset_time Intervention onset time.
#' @param activation_amplitude Amplitude of group activation.
#' @param activation_shape Functional activation shape.
#' @param common_activation
#'   Either "gamma", "logistic", or "none".
#' @param seed Optional random seed.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{df}: simulated PET data
#'     \item \code{truth}: true latent trajectories
#'     \item \code{par_in}: default initialization values
#'     \item \code{predvars}: predictor specification
#'     \item \code{estpat}: estimation restrictions
#'     \item \code{control}: optimization settings
#'   }
#'
#' @export

# n = 50
# t_obs = seq(0.5, 40.5, by = 1)
# onset_time = 20
# activation_amplitude = 0.05
# activation_shape = "gamma"
simulate_lgpe_pet <- function(
    n = 50,
    #t_obs = seq(0.5, 40.5, by = 1),
    t_obs = seq(0.5, 60.5, by = 1),
    #t_obs = seq(0.5, 30.5, by = 1),
    #t_obs = seq(0.5, 60.5, by = 1),
    #t_obs = seq(0.5, 60.5, by = 1),
    #t_obs = seq(0.5, 60.5, by = 1),
    
    #onset_time = 20,
    onset_time = 40,
    #onset_time = 10,
    #onset_time = 10,
    #onset_time = 30,
    
    activation_amplitude = 0.012,
    activation_shape = "gamma",
    common_activation = TRUE,
    seed = NULL
){
  
  ##########################################################
  # Seed
  ##########################################################
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  ##########################################################
  # Predictor structure
  ##########################################################
  
  predvars <- list(
    NULL,   # R1
    NULL,         # k2
    c("group")    # k2a(t)
  )
  
  p_r1  <- length(predvars[[1]]) + 1
  p_k2  <- length(predvars[[2]]) + 1
  p_k2a <- length(predvars[[3]]) + 1
  
  ##########################################################
  # Subject-level variables
  ##########################################################
  
  ids <- 1:n
  
  group <- rbinom(
    n,
    size = 1,
    prob = 0.5
  )
  
  ##########################################################
  # Observation times
  ##########################################################
  
  ts <- t_obs - onset_time
  ts[ts < 0] <- 0
  
  m_j <- length(ts)
  # Shared estimation grid: one point for the entire pre-onset period (ts = 0)
  # plus every distinct post-onset time.
  m <- length(unique(ts))
  
  ##########################################################
  # True scalar parameters
  ##########################################################
  
  #beta_r1_true <- c(1, 0.1)
  beta_r1_true <- 1
  beta_k2_true <- 0.5
  
  sigma_alpha_r1 <- 0.12^2
  sigma_alpha_k2 <- 0.10^2
  
  ##########################################################
  # Subject-specific scalar parameters
  ##########################################################
  
  #x_r1 <- cbind(1, group)
  
  x_r1 <- 1
  
  r1_true <- rnorm(
    n,
    mean = as.vector(x_r1 %*% beta_r1_true),
    sd = sqrt(sigma_alpha_r1)
  )
  
  k2_true <- rnorm(
    n,
    mean = beta_k2_true,
    sd = sqrt(sigma_alpha_k2)
  )
  
  ##########################################################
  # Reference TAC
  ##########################################################
  
  cr_true <- dgamma(
    t_obs,
    shape = 1.5,
    scale = 15
  )
  
  cr_true <- 6000 * cr_true / max(cr_true)
  
  ##########################################################
  # Activation shape
  ##########################################################
  
  if(activation_shape == "gamma"){
    
    act_shape <- dgamma(
      ts,
      #shape = 3,
      #scale = 4
      shape = 4,
      scale = 4
    )
    
    act_shape <- act_shape / max(act_shape)
    
  } else if(activation_shape == "logistic"){
    
    act_shape <- 1 / (1 + exp(-0.5*(ts - 10)))
    act_shape <- act_shape / max(act_shape)
    
  } else if(activation_shape == "none"){
    
    act_shape <- rep(0, m)
    
  } else {
    
    stop("Unknown activation_shape.")
  }
  
  ##########################################################
  # True population-level functional effects
  ##########################################################
  
  beta_k2a_0_true <- rep(0.1, m_j)
  
  beta_k2a_1_true <-
    activation_amplitude * act_shape
  
  ##########################################################
  # GP hyperparameters
  ##########################################################
  
  #alpha_theta0 <- 0.02^2
  alpha_theta0 <- 0.01^2
  #alpha_theta0 <- 0.004^2
  #alpha_theta0 <- 0.002^2
  
  #theta1 <- 40
  #theta1 <- 20
  theta1 <- 10
  theta_nu <- 1e-9
  
  ##########################################################
  # GP covariance matrix
  ##########################################################
  
  sigma_alpha_k2a <- matrix(
    0,
    nrow = m_j,
    ncol = m_j
  )
  
  for(i in 1:m_j){
    
    for(j in 1:m_j){
      
      sigma_alpha_k2a[i, j] <-
        rbf_kernel(
          ts[i],
          ts[j],
          alpha_theta0,
          theta1
        )
    }
  }
  
  diag(sigma_alpha_k2a) <-
    diag(sigma_alpha_k2a) + theta_nu
  
  ##########################################################
  # Subject-specific GP random effects
  ##########################################################
  
  alpha_k2a_true <- mvtnorm::rmvnorm(
    n = n,
    sigma = sigma_alpha_k2a
  )
  
  if(common_activation==FALSE){
    alpha_k2a_true[which(group==0),] <- rep(alpha_k2a_true[which(group==0),1],m_j)  
  }else{
    for(j in 1:n){
      alpha_k2a_true[j,] <- rep(rnorm(1,sd=sqrt(alpha_theta0)),m_j)
    }
  }
  
  
  ##########################################################
  # Subject-specific functional trajectories
  ##########################################################
  
  k2a_true <- vector("list", n)
  
  for(j in 1:n){
    
    k2a_true[[j]] <-
      beta_k2a_0_true +
      group[j] * beta_k2a_1_true +
      alpha_k2a_true[j, ]
  }
  
  ##########################################################
  # Generate TAC observations
  ##########################################################
  
  #var_prop <- 3
  var_prop <- 10
  #var_prop <- 20
  #var_prop <- 30
  
  df_list <- vector("list", n)
  
  for(j in 1:n){
    
    tac_obj <- simulate_tac(
      t_obs,
      cr_true,
      r1_true[j],
      k2_true[j],
      k2a_true[[j]],
      var_prop
    )
    
    df_list[[j]] <- data.frame(
      id = ids[j],
      time = t_obs,
      cr = tac_obj$cr_obs,
      ct = tac_obj$ct_obs,
      wr = tac_obj$cr_obs,
      wt = tac_obj$ct_obs,
      cr_true = tac_obj$cr_true,
      ct_true = tac_obj$ct_true,
      group = group[j],
      onset_time = onset_time
    )
  }
  
  df <- do.call(rbind, df_list)
  
  ##########################################################
  # Estimation restrictions
  ##########################################################
  
  estpat <- list()
  
  estpat$kernel_par <- list()
  
  estpat$kernel_par$beta_theta0 <-
    matrix(
      c(1, 1),
      nrow = 1
    )
  
  ##########################################################
  # Initial values
  ##########################################################
  
  par_in <- list()
  
  ############################
  # Subject-level parameters
  ############################
  
  par_in$R1 <- rep(1, n)
  
  par_in$k2 <- rep(0.5, n)
  
  par_in$k2a <- lapply(
    1:n,
    function(x){
      rep(0.1, m_j)
    }
  )
  
  ############################
  # Random-effect covariance
  ############################
  
  par_in$Sigma_alpha_r1 <-
    matrix(
      sigma_alpha_r1,
      nrow = 1,
      ncol = 1
    )
  
  par_in$Sigma_alpha_k2 <-
    matrix(
      sigma_alpha_k2,
      nrow = 1,
      ncol = 1
    )
  
  ############################
  # Prior covariance matrices
  ############################
  
  par_in$Sigma_beta_r1 <-
    diag(1, p_r1)
  
  par_in$Sigma_beta_k2 <-
    matrix(
      1,
      nrow = 1,
      ncol = 1
    )

  # Prior variance of the scalar reference-group k2a level. This is separate
  # from the GP variance governing deviations from that level.
  par_in$Sigma_beta_k2a_intercept <- 0.1^2
  
  ############################
  # Prior means
  ############################
  
  par_in$Mu_beta_r1 <-
    matrix(
      rep(0, p_r1),
      ncol = 1
    )
  
  par_in$Mu_beta_k2 <-
    matrix(
      0,
      ncol = 1
    )
  
  par_in$Mu_beta_k2a <-
    matrix(
      rep(0, m * p_k2a),
      ncol = 1
    )
  
  ############################
  # Kernel parameters
  ############################
  
  par_in$beta_theta0 <-
    #rep(0.03^2, p_k2a)
  rep(1^2, p_k2a)
  
  par_in$alpha_theta0 <- alpha_theta0
  
  par_in$theta1 <- theta1
  
  par_in$theta_nu <- theta_nu
  
  par_in$par_beta_theta0 <- 0.01
  
  ############################
  # Noise parameter
  ############################
  
  par_in$var_prop <- var_prop
  
  ##########################################################
  # Optimization settings
  ##########################################################
  
  control <- list(
    tol = 1e-4,
    maxiter = 1000,
    plot_iterations = FALSE
  )
  
  ##########################################################
  # Truth object
  ##########################################################
  
  truth <- list(
    beta_k2a_0 = beta_k2a_0_true,
    beta_k2a_1 = beta_k2a_1_true,
    k2a = k2a_true,
    ts = ts
  )
  
  ##########################################################
  # Collect output
  ##########################################################
  
  res <- list()
  
  res$df <- df
  
  res$truth <- truth
  
  res$predvars <- predvars
  
  res$par_in <- par_in
  
  res$estpat <- estpat
  
  res$control <- control
  
  res$beta_r1_true <- beta_r1_true
  res$beta_k2_true <- beta_k2_true
  res$beta_k2a_0_true <- beta_k2a_0_true
  res$r1_true <- r1_true
  res$k2_true <- k2_true
  
  return(res)
}

