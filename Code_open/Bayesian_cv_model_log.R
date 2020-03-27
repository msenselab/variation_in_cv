# Bayesian modeling and simple logarithmic modeling
## ----  Initializing parameters ----
dur = c(0.3,0.49,0.81, 1.33, 2.19, 3.6, 5.92, 9.73, 16)
M_p = exp(mean(log(dur)))
mean_block<- c(exp(mean(log(dur[1:3]))),exp(mean(log(dur[4:6]))),exp(mean(log(dur[7:9]))))

exp_name <- c('Vis/Mix','Vis/Block','Aud/Mix','Aud/Block')
data = msRepr

# assign mean values according to the block number
data$mean <- data$block
data[data$design == 'Mix',]$mean <- M_p
data[data$design == 'Block',]$mean  <- mean_block[data[data$design == 'Block',]$block]

data_m = data %>% 
  group_by(exp,duration) %>%
  summarise(mcv = mean(cv),
            mblock = mean(block))

x<-seq(0.2,16,0.1)
x_block <- as.data.frame( c(seq(0.2,0.82,0.1),
                            seq(1.3,3.7,0.1),
                            seq(5.8,16.1,0.1))
)
x_block <- set_names(x_block, 'm')
x_block$block <- 0
l <-nrow(x_block)
m <-seq(l)

# 2nd parameter in the function, mean values from blocks
for(i in 1:l){
  if(x_block$m[i]< 1) {
    
    m[i] <- mean_block[1]
    x_block$block[i] <- 1
    
  } else if(x_block$m[i] >=1 & x_block$m[i]<4) {
    m[i] <- mean_block[2]
    x_block$block[i] <- 2
  } else {
    m[i] <- mean_block[3] 
    x_block$block[i] <- 3
  }
}

linear_cv <- function (t, par) {
  # Parameters: 1: slope, 2: intercept_short , 3: intercept_medium , 4: intercept_long
  
  slope <- par[1]
  intercept <- par[2]
  y <- intercept + slope*t
  
}

Bayesian_mdur <- function(t, par) {
  
  sp2 <- par[1]^2
  st2 <- par[2]^2
  res <-par[3]
  
  y <- sp2/(sp2+st2)*t+st2/(sp2+st2)*(M_p+res)
  
}

Bayesian_vdur <- function(t, par) {
  
  sp2 <- par[1]^2
  st2 <- par[2]^2
  
  y <- sqrt(sp2*st2/(sp2+st2))
  
}

Bayesian_cv_mean_fit <- function(par, data) {
  
  t <- log(data$duration)
  y <- log(data$reproduction)
  
  var_pred <- Bayesian_vdur(t, par[1:2])
  m_pred <- Bayesian_mdur(t, par[1:3])
  
  a<- -sum(log(dnorm(y, m_pred, var_pred)))
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

Bayesian_cv_sd_fit <- function(par, data, fixed_par) {
  
  t <- log(data$duration)
  y <- data$sdRepr
  
  sig_p <- sqrt(par[2]/(1-fixed_par[1]))
  sig_t <- sqrt(par[2]/fixed_par[1])
  res <- fixed_par[2]
  
  var_pred <- Bayesian_vdur(t, c(sig_p, sig_t))
  m_pred <- Bayesian_mdur(t, c(sig_p, sig_t, res))
  sig_m <- par[1]
  
  pred_sd <- sqrt(predicted_linear_scale_sd(m_pred, var_pred)^2 + sig_m^2)
  
  a <- -sum(log(dlnorm(y, log(pred_sd), 1)))
  
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

Bayesian_cv_full_dist_fit <- function(par, data, fixed_par) {
  
  t <- log(data$duration)
  y <- data$reproduction
  
  sig_p <- sqrt(par[2]/(1-fixed_par[1]))
  sig_t <- sqrt(par[2]/fixed_par[1])
  res <- fixed_par[2]
  
  var_pred <- Bayesian_vdur(t, c(sig_p, sig_t))
  m_pred <- Bayesian_mdur(t, c(sig_p, sig_t, res))
  sig_m <- par[1]
  
  a<- -sum(log(dnormlnorm(y, m_pred, var_pred, sig_m)))
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

# On log scale the model predicts constant SD. This function translates that into an SD on linear
# scale. If I am not mistaken the SD on linear scale should be the SD of the log-normal 
# distribution (not the sigma parameter but the actual SD of the distibution).
predicted_linear_scale_sd <- function(mu, sigma) {
  s <- sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2))
}

linear_cv_fit <- function(par, data) {
  
  t <- data$log_dur
  y_pred <- linear_cv(t, par)
  y <- data$cv
  a <-  sum((y-y_pred)^2)
  
}


# ---- Individual subject fits ----

Bayesian_isub_pre <- data.frame()
bayes_isub_fit_para <- data.frame()
for(i in c(1,3)){
  subs <- unique(filter(vdat, exp==exp_name[i])$sub)
  for(s in subs) {
    data_model = vdat %>% filter(exp == exp_name[i], sub==s) 
    ptm <- proc.time()
    par <- optim(c(0.7,0.3,0), Bayesian_cv_mean_fit, NULL, method = "L-BFGS-B",
                 lower = c(0.01,0.01,-10), upper = c(3,3,10),
                 data_model)
    print(proc.time() - ptm)
    print(paste0(exp_name[i], ': sub=', s))
    print(par$par)
    print(par$value)
    print(par$convergence)
    print(par$message)
    fixed_par = c(par$par[1]^2/(par$par[1]^2+par$par[2]^2), par$par[3])    
    ptm <- proc.time()
    
    data_sd <- data_model %>% group_by(duration) %>% mutate(sdRepr=sd(reproduction)) 
    par <- optim(c(0.2,0.3), Bayesian_cv_sd_fit, NULL, method = "L-BFGS-B",
                 lower = c(0.01,0.01), upper = c(1,3),
                 data=data_sd, fixed_par=fixed_par) 
    
    print(proc.time() - ptm)
    print(paste0(exp_name[i], ': sub=', s)) 
    print(par$par)
    print(par$value)
    print(par$convergence)
    print(par$message)
    
    sig_p <- sqrt(par$par[2]/(1-fixed_par[1]))
    sig_t <- sqrt(par$par[2]/fixed_par[1])
    res <- fixed_par[2]
    
    bayes_isub_fit_para <- rbind2(bayes_isub_fit_para,
                                  data.frame(exp=exp_name[i], sub=s, sig_m=par$par[1], sig_p=sig_p,
                                             sig_t=sig_t, res=fixed_par[2], val=par$value))
    mu = exp(Bayesian_mdur(log(x), c(sig_p, sig_t, res)))
    sigma = sqrt(predicted_linear_scale_sd(Bayesian_mdur(log(x), c(sig_p, sig_t, res)), 
                                           Bayesian_vdur(log(x), c(sig_p, sig_t)))^2 + par$par[1]^2)
    data_temp = data.frame(x=x, y_mdur = mu, 
                           y_vdur = sigma,
                           y_cv=sigma/mu,
                           y_rre = (mu - x)/x, exp = exp_name[i], sub=s)
    Bayesian_isub_pre <- rbind2(Bayesian_isub_pre,data_temp)
    
    print(exp_name[i])
    print(s)

    
    
  }
}
# ---- 1-prior Model - figure_separate ----
# 4 figures from 4 exp using 1st model with 1 prior
# extracting the Bayesian data to integrate into the original figure/Yue
dat1_mix <- data.frame()

for(i in c(1,3)){
  
  Bayesian_sep <- Bayesian_isub_pre %>% filter(exp == exp_name[i]) %>% group_by(x) %>%
    summarize(y_mdur=mean(y_mdur), y_vdur=mean(y_vdur), y_cv=mean(y_cv), y_rre=mean(y_rre))
  #  linear_sep <- linear_pre %>% filter(exp ==exp_name[i])
  datam_sep <- data_m %>% filter(exp == exp_name[i])
  data_model = data %>% filter(exp == exp_name[i]) 
  
  # save the modeling data for using outside this function /Yue
  if(i==1) {  
    dat1_mix <- Bayesian_sep
    dat1_mix$exp = exp_name[i]
  }else{
    Bayesian_sep$exp = exp_name[i]
    dat1_mix = rbind(dat1_mix, Bayesian_sep)
  }
}

## ---- 3-prior model: take 'block' factor into consideration ----
# M_p + bias(block)
# M_p(block)
# each block requires one mean value
# ---- Build up 3-prior Models(Bayesian and linear) ----
Bayesian_mdur_block <- function(t, block, par) {
  
  sp2 <- par[1]^2
  st2 <- par[2]^2
  res <-par[3:5]
  
  y <- sp2/(sp2+st2)*t+st2/(sp2+st2)*(log(mean_block[block])+res[block])
  
}

Bayesian_cv_block <- function(t, block, par) {
  # Parameters: 1: sig_m, 2: cp, 3: c0, 4: res_short, 5: res_medium, 6: res_long
  y <- Bayesian_vdur_block(t,block,par[1:3])/Bayesian_mdur_block(t,block,par[2:6])
}

Bayesian_cv_block_fit <- function(par, data) {
  # Parameters: 1: sig_m, 2: cp, 3: c0, 4: res_short, 5: res_medium, 6: res_long
  block <- data$block
  t <- data$duration
  y <- data$reproduction
  
  var_pred <- Bayesian_vdur_block(t, block, par[1:3])
  m_pred <- Bayesian_mdur_block(t, block, par[2:6])
  
  a<- -sum(log(dnorm(y, m_pred, var_pred)))
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
  
}

Bayesian_cv_mean_fit_block <- function(par, data) {
  
  t <- log(data$duration)
  y <- log(data$reproduction)
  block <- data$block
  
  var_pred <- Bayesian_vdur(t, par[1:2])
  m_pred <- Bayesian_mdur_block(t, block, par[1:5])
  
  a<- -sum(log(dnorm(y, m_pred, var_pred)))
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

Bayesian_cv_sd_fit_block <- function(par, data, fixed_par) {
  
  t <- log(data$duration)
  y <- data$sdRepr
  block <- data$block
  
  sig_p <- sqrt(par[2]/(1-fixed_par[1]))
  sig_t <- sqrt(par[2]/fixed_par[1])
  res <- fixed_par[2:4]
  
  var_pred <- Bayesian_vdur(t, c(sig_p, sig_t))
  m_pred <- Bayesian_mdur_block(t, block, c(sig_p, sig_t, res))
  sig_m <- par[1]
  
  pred_sd <- sqrt(predicted_linear_scale_sd(m_pred, var_pred)^2 + sig_m^2)
  
  a <- -sum(log(dlnorm(y, log(pred_sd), 1)))
  
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

linear_cv_block <- function (t, block, par) {
  # Parameters: 1: slope, 2: intercept_short , 3: intercept_medium , 4: intercept_long
  
  slope <- par[1]
  intercept <- par[2:4]
  y <- intercept[block] + slope*t
  
}

linear_cv_block_fit <- function(par, data) {
  
  block <- data$block
  t <- data$log_dur
  
  y_pred <- linear_cv_block(t, block, par)
  y <- data$cv
  a <-  sum((y-y_pred)^2)
  
}

#  ---- Initializing parameters ---- 
Bayesian_pre_block <- data.frame()
linear_pre_block <- data.frame()

bayes_block_fit_para <- data.frame()
linear_block_fit_para <- data.frame()
aic_cv <- data.frame()

# ---- Model of blocked exp. with variable perceptual representation noise ----

Bayesian_vdur_4 <- function(t, sig_p, sig_t) {
  
  sp2 <- sig_p^2
  st2 <- sig_t^2
  
  y <- sqrt(sp2*st2/(sp2+st2))
  
}

Bayesian_mdur_block_4 <- function(t, block, sig_p, sig_t, res) {
  
  sp2 <- sig_p^2
  st2 <- sig_t^2
  
  y <- sp2/(sp2+st2)*t+st2/(sp2+st2)*(log(mean_block[block])+res[block])
  
}

Bayesian_cv_sd_fit_block_4 <- function(par, data, fixed_par) {
  
  t <- log(data$duration)
  y <- data$sdRepr
  block <- data$block
  sig <- par[2:4]
  
  sig_p <- sqrt(sig[block]/(1-fixed_par[1]))
  sig_t <- sqrt(sig[block]/fixed_par[1])
  res <- fixed_par[2:4]
  
  var_pred <- Bayesian_vdur_4(t, sig_p, sig_t)
  m_pred <- Bayesian_mdur_block_4(t, block, sig_p, sig_t, res)
  sig_m <- par[1]
  
  pred_sd <- sqrt(predicted_linear_scale_sd(m_pred, var_pred)^2 + sig_m^2)
  
  a <- -sum(log(dlnorm(y, log(pred_sd), 1)))
  
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

#  --- Initializing parameters --- 
Bayesian_pre_block_4 <- data.frame()
linear_pre_block <- data.frame()

bayes_block_fit_para_4 <- data.frame()
linear_block_fit_para <- data.frame()
aic_cv <- data.frame()

# ---  Model fitting based on'Block' data  ---
for(i in c(2,4)){
  subs <- unique(filter(vdat, exp==exp_name[i])$sub)
  for(s in subs) {
    data_model = vdat %>% filter(exp == exp_name[i], sub==s) 
    
    ptm <- proc.time()
    par <- optim(c(0.7,0.3,0,0,0), Bayesian_cv_mean_fit_block, NULL, method = "L-BFGS-B", #Yuev4
                 lower = c(0.01,0.01,-10,-10,-10), upper = c(3,3,10,10,10),
                 data_model)
    print(proc.time() - ptm)
    print(paste0(exp_name[i], ': sub=', s))
    print(par$par)
    print(par$value)
    print(par$convergence)
    print(par$message)
    fixed_par = c(par$par[1]^2/(par$par[1]^2+par$par[2]^2), par$par[3], par$par[4], par$par[5])    
    ptm <- proc.time()
    
    data_sd <- data_model %>% group_by(duration) %>% mutate(sdRepr=sd(reproduction)) 
    ptm <- proc.time()
    par <- optim(c(0.2,0.3,0.3,0.3), Bayesian_cv_sd_fit_block_4, NULL, method = "L-BFGS-B",
                 lower = c(0.001,0.001,0.001,0.001), upper = c(3,3,3,3), 
                 data=data_sd, fixed_par=fixed_par) 
    print(proc.time() - ptm)
    print(paste0(exp_name[i], ': sub=', s)) 
    print(par$par)
    print(par$value)
    print(par$convergence)
    print(par$message)
    
    sig_p <- sqrt(par$par[2:4]/(1-fixed_par[1]))
    sig_t <- sqrt(par$par[2:4]/fixed_par[1])
    res <- fixed_par[2:4]
    sig_m <- par$par[1]
    
    bayes_block_fit_para_4 <- rbind2(bayes_block_fit_para_4,
                                     data.frame(exp=exp_name[i], sub=s, sig_m=sig_m, 
                                                sig_p_medium=sig_p[2], sig_p_long=sig_p[3], 
                                                sig_p_short=sig_p[1], sig_t_short=sig_t[1],
                                                sig_t_medium=sig_t[2], sig_t_long=sig_t[3],
                                                res_short=fixed_par[2], 
                                                res_medium=fixed_par[3], 
                                                res_long=fixed_par[4], val=par$value))
    mu = exp(Bayesian_mdur_block_4(log(x_block$m), x_block$block, sig_p[x_block$block], sig_t[x_block$block], res))
    sigma = sqrt(predicted_linear_scale_sd(Bayesian_mdur_block_4(log(x_block$m), x_block$block, sig_p[x_block$block], sig_t[x_block$block], res), 
                                           Bayesian_vdur_4(log(x_block$m), sig_p[x_block$block], sig_t[x_block$block]))^2 + sig_m^2)
    data_temp = data.frame(x=x_block$m, y_mdur = mu, 
                           y_vdur = sigma,
                           y_cv=sigma/mu,
                           y_rre = (mu - x_block$m)/x_block$m, exp = exp_name[i], sub=s)
    Bayesian_pre_block_4 <- rbind2(Bayesian_pre_block_4,data_temp)
    
  }
}

Bayesian_sep_block <- Bayesian_pre_block_4 %>% group_by(exp, x) %>%
  summarize(y_mdur=mean(y_mdur), y_vdur=mean(y_vdur), y_cv=mean(y_cv), y_rre=mean(y_rre))
dat1_block <- data.frame()
dat1_block <- Bayesian_sep_block

for(i in c(2,4)) {
  
  data_model = data %>% filter(exp == exp_name[i])
  Bayesian_sep = Bayesian_sep_block %>% filter(exp == exp_name[i])
  
  mean_data <- group_by(data_model, duration) %>% summarize(mmRRE=mean(mRepErr/duration))
  Fig_rre_model<-
    ggplot(data_model, aes(x=duration, y=mRepErr/duration, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() +  scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10)) +
    geom_point(data=mean_data, aes(x=duration, y=mmRRE, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=mmRRE, group=1), size=1, color="black") +
    geom_line(data=Bayesian_sep, aes(x=x, y=y_rre, group=1), color="blue", size=1) +
    labs(title = exp_name[i])
  print(Fig_rre_model)
  
  mean_data <- group_by(data_model, duration) %>% summarize(mmRepr=mean(mRepr))
  Fig_mdur_model<-
    ggplot(data_model, aes(x=duration, y=mRepr, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() + 
    geom_point(data=mean_data, aes(x=duration, y=mmRepr, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=mmRepr, group=1), size=1, color="black") +
    geom_line(data=Bayesian_sep, aes(x=x, y=y_mdur, group=1), color="blue", size=1) +
    scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
    labs(title = exp_name[i])
  print(Fig_mdur_model)
  
  mean_data <- group_by(data_model, duration) %>% summarize(msdRepr=mean(sdRepr))
  Fig_vdur_model<-
    ggplot(data_model, aes(x=duration, y=sdRepr, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() + 
    geom_point(data=mean_data, aes(x=duration, y=msdRepr, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=msdRepr, group=1), size=1, color="black") +
    geom_line(data=Bayesian_sep, aes(x=x, y=y_vdur, group=1), color="blue", size=1) +
    scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
    labs(title = exp_name[i])
  print(Fig_vdur_model)
  
}

#View(bayesian_block_fit_para)
#View(aic_cv) /Yue: error in integration
# ---- Figure_Block 4 ----
data_block <- data %>% filter(design == 'Block')
data_m_block <- data_m %>% filter(exp == 'Vis/Block'|exp == 'Aud/Block')

Fig_cv_model_block_v4 <- 
  ggplot(data= data_block, aes(x= duration, y = cv)) + 
  scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10))+
  scale_y_continuous(limits = c(0.05,0.75))+
  geom_point(size = 0.4, colour= 'grey')+
  geom_line(data =data_m_block , aes(x= duration, y= mcv), color= 'darkgrey')+
  # data from our model
  geom_line(aes(x=x,y=y_cv),data=Bayesian_sep_block %>%
              filter(exp == 'Vis/Block' | exp =='Aud/Block'),
            color = 'red2',size=0.6 )+
  facet_wrap(~exp, nrow = 1)+
  geom_point(data =data_m_block , aes(x= duration, y= mcv), color= 'grey17',size = 1.6, shape = 17)+
  #annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend-Inf)+
  xlab('Duration(s)') +
  ylab('')+
  theme_minimal()+
  theme(legend.position= "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        axis.ticks = element_line(color = "grey"),
        axis.text.y =element_blank())

Fig_cv_model_block_v4



