library(lgpesrtm)
library(ggplot2)
library(dplyr)
library(patchwork)

N <- 50
Activation_Amplitude <- 0.05

fit_and_extract_k2a <- function(sim, activation_shape_label) {
  
  fit_lgpe <- fit_lgpe_srtm(
    df       = sim$df,
    predvars = sim$predvars,
    par_in   = sim$par_in,
    estpat   = sim$estpat,
    control  = sim$control
  )
  
  Data <- sim$df
  IDs <- unique(Data$id)
  N <- length(IDs)
  time_unique <- sort(unique(Data$time))
  T_onset <- Data$onset_time[1]
  
  lpntpet_fits <- vector("list", N)
  
  for (j in seq_along(IDs)) {
    df_j <- Data[Data$id == IDs[j], ]
    
    lpntpet_fits[[j]] <- fit_lpntpet_catalog(
      time     = df_j$time,
      ct       = df_j$ct,
      cr       = df_j$cr,
      t_onset  = df_j$onset_time[1],
      weights  = 1 / df_j$wt
    )
  }
  
  DF_k2a <- bind_rows(lapply(seq_along(IDs), function(j) {
    
    df_j <- Data[Data$id == IDs[j], ]
    fit_j <- lpntpet_fits[[j]]
    
    k2a_lpntpet <- fit_j$estimates$k2a +
      fit_j$estimates$gamma *
      gamma_variate_h(
        time_unique,
        fit_j$estimates$tD,
        fit_j$estimates$tP,
        fit_j$estimates$alpha
      )
    
    bind_rows(
      data.frame(
        id = IDs[j],
        time = df_j$time,
        value = sim$truth$k2a[[j]],
        method = "True",
        group = df_j$group[1],
        activation_shape = activation_shape_label
      ),
      data.frame(
        id = IDs[j],
        time = df_j$time,
        value = fit_lgpe$k2a[[j]],
        method = "LGPE-SRTM",
        group = df_j$group[1],
        activation_shape = activation_shape_label
      ),
      data.frame(
        id = IDs[j],
        time = df_j$time,
        value = k2a_lpntpet,
        method = "lp-ntPET",
        group = df_j$group[1],
        activation_shape = activation_shape_label
      )
    )
  }))
  
  DF_true_mean <- expand.grid(
    time = time_unique,
    group = c(0, 1),
    method = c("True", "LGPE-SRTM", "lp-ntPET"),
    activation_shape = activation_shape_label
  ) %>%
    mutate(
      true_value = ifelse(
        group == 0,
        sim$truth$beta_k2a_0,
        sim$truth$beta_k2a_0 + sim$truth$beta_k2a_1
      )
    )
  
  list(
    subj = DF_k2a,
    true_mean = DF_true_mean,
    onset = T_onset,
    fit_lgpe = fit_lgpe,
    lpntpet_fits = lpntpet_fits
  )
}





sim_gamma <- simulate_lgpe_pet(
  n = N,
  activation_amplitude = Activation_Amplitude,
  activation_shape = "gamma",
  #seed = 20260534
  seed = 20260566
)



sim_logistic <- simulate_lgpe_pet(
  n = N,
  activation_amplitude = Activation_Amplitude,
  activation_shape = "logistic",
  #seed = 20260534
  seed = 20260565
)

res_gamma <- fit_and_extract_k2a(sim_gamma, "Gamma")
res_logistic <- fit_and_extract_k2a(sim_logistic, "Logistic")

#Collect results for plotting
#################

DF_subj <- bind_rows(res_gamma$subj,res_logistic$subj) %>%
  mutate(
    id = factor(id),
    method = factor(method, levels = c("True", "LGPE-SRTM", "lp-ntPET")),
    group = factor(group, levels = c(0, 1), labels = c("Control", "Activation")),
    activation_shape = factor(activation_shape, levels = c("Gamma", "Logistic"))
  )

DF_mean <- DF_subj %>%
  group_by(activation_shape, method, group, time) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")


#Take the estimated betas for group level 







#Take the betas for the respective methods
DF_mean[DF_mean$activation_shape=="Gamma" & DF_mean$method=="LGPE-SRTM" & DF_mean$group=="Control" & DF_mean$time>res_gamma$fit_lgpe$df$onset_time[1]-1,]$mean_value <- res_gamma$fit_lgpe$beta_k2a_intercept + res_gamma$fit_lgpe$beta_k2a_mat[,1]
DF_mean[DF_mean$activation_shape=="Gamma" & DF_mean$method=="LGPE-SRTM" & DF_mean$group=="Activation" & DF_mean$time>res_gamma$fit_lgpe$df$onset_time[1]-1,]$mean_value <- res_gamma$fit_lgpe$beta_k2a_intercept + res_gamma$fit_lgpe$beta_k2a_mat[,1] + res_gamma$fit_lgpe$beta_k2a_mat[,2]

DF_mean[DF_mean$activation_shape=="Logistic" & DF_mean$method=="LGPE-SRTM" & DF_mean$group=="Control" & DF_mean$time>res_logistic$fit_lgpe$df$onset_time[1]-1,]$mean_value <- res_logistic$fit_lgpe$beta_k2a_intercept + res_logistic$fit_lgpe$beta_k2a_mat[,1]
DF_mean[DF_mean$activation_shape=="Logistic" & DF_mean$method=="LGPE-SRTM" & DF_mean$group=="Activation" & DF_mean$time>res_logistic$fit_lgpe$df$onset_time[1]-1,]$mean_value <- res_logistic$fit_lgpe$beta_k2a_intercept + res_logistic$fit_lgpe$beta_k2a_mat[,1] + res_logistic$fit_lgpe$beta_k2a_mat[,2]




DF_true_mean <- bind_rows(res_gamma$true_mean, res_logistic$true_mean) %>%
  mutate(
    method = factor(method, levels = c("True", "LGPE-SRTM", "lp-ntPET")),
    group = factor(group, levels = c(0, 1), labels = c("Control", "Activation")),
    activation_shape = factor(activation_shape, levels = c("Gamma", "Logistic"))
  )

T_onset <- res_gamma$onset


plot_k2a_shape <- function(shape_name) {
  
  DF_subj_tmp <- DF_subj %>%
    filter(activation_shape == shape_name)
  
  DF_mean_tmp <- DF_mean %>%
    filter(activation_shape == shape_name)
  
  DF_true_tmp <- DF_true_mean %>%
    filter(activation_shape == shape_name)
  
  ggplot() +
    geom_vline(
      xintercept = T_onset,
      linetype = "dotted",
      colour = "grey30",
      linewidth = 0.5
    ) +
    
    geom_line(
      data = DF_subj_tmp,
      aes(x = time, y = value,
          group = id,
          colour = method),
      alpha = 0.7,
      # linewidth = 0.45
      linewidth = 0.2
    ) +
    
    geom_line(
      data = DF_mean_tmp,
      aes(x = time, y = mean_value,
          colour = method),
      linewidth = 1.2
    ) +
    
    geom_line(
      data = DF_true_tmp,
      aes(x = time, y = true_value),
      colour = "black",
      linetype = "33",
      linewidth = 1.25
    ) +
    
    scale_colour_manual(
      values = c(
        "True" = "grey",
        "LGPE-SRTM" = "#F8766D",
        "lp-ntPET" = "#00BFC4"
      )
    ) +
    
    facet_grid(group ~ method) +
    
    coord_cartesian(ylim = c(0.085, 0.1+Activation_Amplitude+0.02), xlim = c(res_gamma$fit_lgpe$df$onset_time[1],max(res_gamma$fit_lgpe$df$time))) +
    
    labs(
      title = paste("Activation shape = ",shape_name),
      x = "Time (min)",
      y = expression(hat(k)[2*a](t))
    ) +
    
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_rect(
        fill = "grey90",
        colour = "grey30"
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}


p_k2a <- plot_k2a_shape("Gamma") / plot_k2a_shape("Logistic")



#Start with method=truth, activation_shape=gamma
beta <- DF_true_mean[DF_true_mean$method=="True" & DF_true_mean$group=="Control" & DF_true_mean$activation_shape=="Gamma" & DF_true_mean$time>T_onset-1,]$true_value
M <- length(beta)
M_j <- length(unique(DF_true_mean$time))
Times <- unique(DF_true_mean$time)
se_beta <- rep(0,length(beta))
#Time <- 
DF_true_gamma_beta0 <- data.frame(activation_shape="Gamma" , method="True" , group="Control", beta=beta,se=se_beta,time=Times[(M_j-M+1):M_j])

beta <- sim_gamma$truth$beta_k2a_1[(M_j-M+1):M_j]
DF_true_gamma_beta1 <- data.frame(activation_shape="Gamma" , method="True" , group="Activation", beta=beta,se=se_beta,time=Times[(M_j-M+1):M_j])

#Second: method=Truth, activation_shape=logistic
beta <- DF_true_mean[DF_true_mean$method=="True" & DF_true_mean$group=="Control" & DF_true_mean$activation_shape=="Logistic" & DF_true_mean$time>T_onset-1,]$true_value
se_beta <- rep(0,M)
DF_true_logistic_beta0 <- data.frame(activation_shape="Logistic" , method="True" , group="Control", beta=beta,se=se_beta,time=Times[(M_j-M+1):M_j])

beta <- sim_logistic$truth$beta_k2a_1[(M_j-M+1):M_j]
se_beta <- rep(0,M)
DF_true_logistic_beta1 <- data.frame(activation_shape="Logistic" , method="True" , group="Activation", beta=beta,se=se_beta,time=Times[(M_j-M+1):M_j])

DF_true <- rbind(DF_true_gamma_beta0,DF_true_gamma_beta1,DF_true_logistic_beta0,DF_true_logistic_beta1)


#third: Method=LGPE-SRTM, activation_shape Gamma
beta <- res_gamma$fit_lgpe$beta_k2a_mat[,1]
se_beta <- res_gamma$fit_lgpe$se_beta_k2a_mat[,1]
DF_lgpe_gamma_beta0 <- data.frame(activation_shape="Gamma",method="LGPE-SRTM",group="Control",beta=beta,se=se_beta,time=Times[(M_j-M+1):M_j])

beta <- res_gamma$fit_lgpe$beta_k2a_mat[,2]
se_beta <- res_gamma$fit_lgpe$se_beta_k2a_mat[,2]
DF_lgpe_gamma_beta1 <- data.frame(activation_shape="Gamma",method="LGPE-SRTM",group="Activation",beta=beta,se=se_beta,time=Times[(M_j-M+1):M_j])

#fourth: Method=LGPE-SRTM, activation_shape Logistic
beta <- res_logistic$fit_lgpe$beta_k2a_mat[,1]
se_beta <- res_logistic$fit_lgpe$se_beta_k2a_mat[,1]
DF_lgpe_logistic_beta0 <- data.frame(activation_shape="Logistic",method="LGPE-SRTM",group="Control",beta=beta,se=se_beta,time=Times[(M_j-M+1):M_j])

beta <- res_logistic$fit_lgpe$beta_k2a_mat[,2]
se_beta <- res_logistic$fit_lgpe$se_beta_k2a_mat[,2]
DF_lgpe_logistic_beta1 <- data.frame(activation_shape="Logistic",method="LGPE-SRTM",group="Activation",beta=beta,se=se_beta,time=Times[(M_j-M+1):M_j])

DF_lgpe <- rbind(DF_lgpe_gamma_beta0,DF_lgpe_gamma_beta1,DF_lgpe_logistic_beta0,DF_lgpe_logistic_beta1)


#Fifth: Method=lp-ntPET, activation_shape=Gamma
DF <- sim_gamma$df  
ID_all <- unique(DF$id)
Group <- c()
k2a_g0_mat <- NULL
k2a_g1_mat <- NULL
for(j in 1:N){
    
  DF_j <- DF[DF$id==ID_all[j],]
  Group[j] <- DF_j$group[1]
  if(Group[j]==0){
    DF_subj <- res_gamma$subj
    k2a_g0_j <- DF_subj[DF_subj$activation_shape=="Gamma" & DF_subj$id==ID_all[j] & DF_subj$method=="lp-ntPET",]$value
    k2a_g0_mat <- rbind(k2a_g0_mat,k2a_g0_j)
  }else{
    DF_subj <- res_gamma$subj
    k2a_g1_j <- DF_subj[DF_subj$activation_shape=="Gamma" & DF_subj$id==ID_all[j] & DF_subj$method=="lp-ntPET",]$value
    k2a_g1_mat <- rbind(k2a_g1_mat,k2a_g1_j)
  }
  
}

n0 <- dim(k2a_g0_mat)[1]
n1 <- dim(k2a_g1_mat)[1]

beta0 <- colMeans(k2a_g0_mat)[(M_j-M+1):M_j]
beta1 <- colMeans(k2a_g1_mat)[(M_j-M+1):M_j] - colMeans(k2a_g0_mat)[(M_j-M+1):M_j]

se_beta0 <- sqrt(diag(cov(k2a_g0_mat[,(M_j-M+1):M_j])))/sqrt(n0)
se_beta1 <- sqrt(diag(cov(k2a_g0_mat[,(M_j-M+1):M_j]))/n0 +diag(cov(k2a_g1_mat[,(M_j-M+1):M_j]))/n1)


DF_lpntpet_gamma_beta0 <- data.frame(activation_shape="Gamma",method="lp-ntPET",group="Control",beta=beta0,se=se_beta0,time=Times[(M_j-M+1):M_j])
DF_lpntpet_gamma_beta1 <- data.frame(activation_shape="Gamma",method="lp-ntPET",group="Activation",beta=beta1,se=se_beta1,time=Times[(M_j-M+1):M_j])


#Sicth: Method=lp-ntPET, activation_shape=Logistic
DF <- sim_logistic$df  
ID_all <- unique(DF$id)
Group <- c()
k2a_g0_mat <- NULL
k2a_g1_mat <- NULL

for(j in 1:N){
  
  DF_j <- DF[DF$id==ID_all[j],]
  Group[j] <- DF_j$group[1]
  if(Group[j]==0){
    DF_subj <- res_logistic$subj
    k2a_g0_j <- DF_subj[DF_subj$activation_shape=="Logistic" & DF_subj$id==ID_all[j] & DF_subj$method=="lp-ntPET",]$value
    k2a_g0_mat <- rbind(k2a_g0_mat,k2a_g0_j)
  }else{
    DF_subj <- res_logistic$subj
    k2a_g1_j <- DF_subj[DF_subj$activation_shape=="Logistic" & DF_subj$id==ID_all[j] & DF_subj$method=="lp-ntPET",]$value
    k2a_g1_mat <- rbind(k2a_g1_mat,k2a_g1_j)
  }
  
}

n0 <- dim(k2a_g0_mat)[1]
n1 <- dim(k2a_g1_mat)[1]

beta0 <- colMeans(k2a_g0_mat)[(M_j-M+1):M_j]
beta1 <- colMeans(k2a_g1_mat)[(M_j-M+1):M_j] - colMeans(k2a_g0_mat)[(M_j-M+1):M_j]

se_beta0 <- sqrt(diag(cov(k2a_g0_mat[,(M_j-M+1):M_j])))/sqrt(n0)
se_beta1 <- sqrt(diag(cov(k2a_g0_mat[,(M_j-M+1):M_j]))/n0 +diag(cov(k2a_g1_mat[,(M_j-M+1):M_j]))/n1)


DF_lpntpet_logistic_beta0 <- data.frame(activation_shape="Logistic",method="lp-ntPET",group="Control",beta=beta0,se=se_beta0,time=Times[(M_j-M+1):M_j])
DF_lpntpet_logistic_beta1 <- data.frame(activation_shape="Logistic",method="lp-ntPET",group="Activation",beta=beta1,se=se_beta1,time=Times[(M_j-M+1):M_j])


DF_lpntpet <- rbind(DF_lpntpet_gamma_beta0,DF_lpntpet_gamma_beta1,DF_lpntpet_logistic_beta0,DF_lpntpet_logistic_beta1)






DF_diff_lpntpet <- data.frame(
  dbeta = DF_lpntpet$beta - DF_true$beta,
  se   = DF_lpntpet$se,
  activation_shape=DF_lpntpet$activation_shape,
  group = DF_lpntpet$group
)

DF_diff_lgpe <- data.frame(
  dbeta = DF_lgpe$beta - DF_true$beta,
  se   = DF_lgpe$se,
  activation_shape=DF_lgpe$activation_shape,
  group = DF_lgpe$group
)




# plot(DF_diff_lpntpet[DF_diff_lpntpet$group=="Activation",]$dbeta,type="l")
# lines(DF_diff_lgpe[DF_diff_lgpe$group=="Activation",]$dbeta,col="red")
# abline(0,0,col="blue")
# 
# 
# plot(DF_diff_lpntpet[DF_diff_lpntpet$group=="Activation",]$se,ylim=c(0,max(DF_diff_lpntpet[DF_diff_lpntpet$group=="Activation",]$se)),type="l")
# lines(DF_diff_lgpe[DF_diff_lgpe$group=="Activation",]$se,col="red")
# abline(0,0,col="blue")


######Plot

DF_plot <- bind_rows(
  DF_lgpe %>% mutate(method = "LGPE-SRTM"),
  DF_lpntpet %>% mutate(method = "lp-ntPET"),
  DF_true %>% mutate(method = "True", se = 0)
) %>%
  mutate(
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se,
    activation_shape = factor(
      activation_shape,
      levels = c("Gamma", "Logistic")
    ),
    method = factor(
      method,
      levels = c("True", "LGPE-SRTM", "lp-ntPET")
    )
  )

DF_plot <- DF_plot %>%
  mutate(
    panel = paste(activation_shape, group, sep = " - ")
  )

method_cols <- c(
  "True"      = "black",
  "LGPE-SRTM" = "#F8766D",
  "lp-ntPET"  = "#00BFC4"
)
###########
#Plot only activation - control
#################
p_activation_beta <- ggplot(
  DF_plot %>% filter(group == "Activation"),
  aes(x = time, y = beta)
) +
  
  # Confidence bands
  geom_ribbon(
    data = DF_plot %>%
      filter(group == "Activation", method != "True"),
    aes(ymin = lower, ymax = upper, fill = method),
    alpha = 0.5
  ) +
  
  # Estimated curves
  geom_line(
    data = DF_plot %>%
      filter(group == "Activation", method != "True"),
    aes(colour = method),
    linewidth = 1
  ) +
  
  # Truth curve on top
  geom_line(
    data = DF_plot %>%
      filter(group == "Activation", method == "True"),
    colour = "black",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  facet_wrap(
    ~ activation_shape,
    scales = "free_y",
    labeller = labeller(
      activation_shape = function(x)
        paste("Activation shape =", tools::toTitleCase(x))
    )
  ) +
  
  scale_colour_manual(values = method_cols[-1]) +
  scale_fill_manual(values = method_cols[-1]) +
  
  labs(
    x = "Time (min)",
    y = expression(beta[1](t)),
    colour = "Method",
    fill = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
################
