# Codes for the Bayesian model fitting and model comparison
# The model assumes that Bayesian integration happens on a logarithmic scale representation
# of the durations, followed by a transformation back to linear scale for the reproduction
# during which some additional decision and motor noise is introduced. The Bayesian integration
# on logarithmic scale on its own would predict constant CV, but the introduction of additional
# duration-independent noise after transforming back to linear scale results in the predicted
# CV being larger for the shortest durations.
#
# For the blocked condition there is additionally a model comparison between 32 different
# models which differ in terms of whether sig_p, sig_s, sig_m and res can differ between short, 
# medium and long duration blocks as well as in whether the bias "res" is added as a shift of 
# the mean of the prior, or later as a bias of the reproduction

## ----- Prepare all necessary data  ----
# This chunk of codes can also be found in the  'Main_Figure.R' with more detailed comments.
# The repetition is included for the code to be run independently.

source('loadPackages.R')
# load data from 4 experiments  
vdat = read.csv("vdata.csv", header = TRUE)
# prepare grand-mean data

msRepr <- vdat %>% group_by(exp, duration,sub) %>% 
  summarise(mRepr = mean(reproduction), sdRepr = sd(reproduction), 
            mPrErr = mean(prod_err,na.rm=TRUE), sdProd = sd(prod_err,na.rm=TRUE), 
            mRepErr = mean(rep_err), sdRepErr = sd(rep_err),
            rre = mean(rep_err/duration))  %>% 
  mutate(cv = sdRepr/mRepr, log_dur = log2(duration)) %>% # cv based on physical duration
  arrange(exp, duration,sub) %>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE)
msRepr$block <- c(rep(1:3,each = 12*3), rep(1:3,each =12*3), rep(1:3, each =11*3), rep(1:3, each =14*3))
msRepr$block <- factor(msRepr$block)


# grand means (collapse subjects)
mRepr <- msRepr %>% group_by(exp, duration) %>% 
  summarise(mRepr = mean(mRepr),  mCV = mean(cv),n = n(),
            seCV = sd(cv)/sqrt(n-1),
            msdProd = mean(sdProd), msdRepr = mean(sdRepr),
            mPrErr = mean(mPrErr), mRE = mean(mRepErr), 
            mRRE = mean(mRepErr/duration),mseRRE = sd(mRepErr/duration)/sqrt(n-1))  %>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE) %>%
  setnames(., 'design','Design')


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

# This function calculates the mean duration according to Bayesian integration with 
# M_p as the empirical mean of the prior, res as a bias of the mean of the subjective prior,
# and sp2 and st2 as the variances of the prior and the likelihood respectively
Bayesian_mdur <- function(t, par) {
  
  sp2 <- par[1]^2
  st2 <- par[2]^2
  res <- par[3]
  
  y <- sp2/(sp2+st2)*t+st2/(sp2+st2)*(log(M_p)+res)
  
}

# This function calculates the standard deviation of the posterior distribution
Bayesian_vdur <- function(sig_p, sig_t) {
  
  sp2 <- sig_p^2
  st2 <- sig_t^2
  
  y <- sqrt(sp2*st2/(sp2+st2))
  
}

# On log scale the model predicts constant SD. This function translates that into an SD on linear
# scale, i.e. to the SD of the log-normal distribution 
# (not the sigma parameter but the actual SD of the distibution).
predicted_linear_scale_sd <- function(mu, sigma) {
  s <- sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2))
}

# Fit a straight line in log-space to get good starting parameters for the full fit
Bayesian_pre_fit <- function(par, data) {

  t <- log(data$duration)
  y <- log(data$reproduction)
  
  sig_p <- par[1]
  sig_t <- par[2]
  res <- par[3]
  
  m_pred <- Bayesian_mdur(t, c(sig_p, sig_t, res))
  v_pred <- Bayesian_vdur(sig_p, sig_t)
  
  a <- -sum(log(dnorm(y, m_pred, v_pred)))
  
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

# Fit the full distribution. Parameters are the same as for the Bayesian_mdur function above,
# plus the additional parameter sig_m which is the stimulus independent component to the 
# variance, which could represent motor noise and decision noise during the reproduction stage
Bayesian_mixed_fit <- function(par, data) {
  
  t <- log(data$duration)
  y <- data$reproduction
  
  sig_p <- par[1]
  sig_t <- par[2]
  res <- par[3]
  
  var_pred <- Bayesian_vdur(sig_p, sig_t)
  m_pred <- Bayesian_mdur(t, c(sig_p, sig_t, res))
  sig_m <- par[4]

  pred_m <- exp(m_pred + var_pred^2/2)
  pred_sd <- sqrt(predicted_linear_scale_sd(m_pred, var_pred)^2 + sig_m^2)
  
  a <- -sum(log(dnorm(y, pred_m, pred_sd)))
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

# ---- Individual subject fits ----

predictions <- data.frame()
# add cv predictions(9 points) for comparison
cv_model_mix <- data.frame()
fit_para <- data.frame()
for(i in c(1,3)){
  subs <- unique(filter(vdat, exp==exp_name[i])$sub)
  for(s in subs) {
    data_model = vdat %>% filter(exp == exp_name[i], sub==s) 
    ptm <- proc.time()
    par <- optim(c(0.6,0.3,0), Bayesian_pre_fit, NULL, method = "L-BFGS-B",
                 lower = c(0.01,0.01,-10), upper = c(3,3,10),
                 data_model, control = c(maxit=10000))
    
    conv1 <- par$convergence
    
    par <- optim(c(par$par,0.1), Bayesian_mixed_fit, NULL, method = "L-BFGS-B",
                 lower = c(0.01,0.01,-10,0), upper = c(3,3,10,3),
                 data_model, control = c(maxit=10000))
    print(proc.time() - ptm)
    print(paste0(exp_name[i], ': sub=', s))
    print(par$par)
    print(par$value)
    print(par$convergence)
    print(par$message)
    
    sig_p <- par$par[1]
    sig_t <- par$par[2]
    w <- sig_p^2/(sig_p^2+sig_t^2)
    res <- par$par[3]
    sig_m <- par$par[4]
    sig_tot <- Bayesian_vdur(sig_p, sig_t)
    conv2 <- par$convergence
    
    fit_para <- rbind2(fit_para, data.frame(exp=exp_name[i], sub=s, sig_p=sig_p,
                                            sig_t=sig_t, sig_m=sig_m, w=w, sig_tot=sig_tot,
                                            res=res, val=par$value, conv1=conv1, conv2=conv2))
    mu = exp(Bayesian_mdur(log(x), c(sig_p, sig_t, res)) + sig_tot^2/2)
    
    sigma = sqrt(predicted_linear_scale_sd(Bayesian_mdur(log(x), c(sig_p, sig_t, res)), 
                                           Bayesian_vdur(sig_p, sig_t))^2 + sig_m^2)
    
    data_temp = data.frame(x=x, y_mdur = mu, 
                           y_vdur = sigma,
                           y_cv=sigma/mu,
                           y_rre = (mu - x)/x, exp = exp_name[i], sub=s)
    
    predictions <- rbind2(predictions,data_temp)

    # extract 9 points of predicted CVs, save in 'cv_model_mix'
    mu_cv = exp(Bayesian_mdur(log(dur), c(sig_p, sig_t, res)) + sig_tot^2/2)
    sigma_cv = sqrt(predicted_linear_scale_sd(Bayesian_mdur(log(dur), c(sig_p, sig_t, res)), 
                                              Bayesian_vdur(sig_p, sig_t))^2 + sig_m^2)    
    # add experimental data into the data frame
    data_model_mean <- data_model %>% group_by(exp, sub, duration)  %>% 
      summarize(rep = mean(reproduction), rep_err = mean(rep_err),
                sd = sd(reproduction), cv = sd/rep)
    
    cv_temp = data.frame(exp = exp_name[i], sub=s,duration = dur, 
                         rep = data_model_mean$rep, rep_err = data_model_mean$rep_err,
                         sd = data_model_mean$sd, cv = data_model_mean$cv,
                         rep_model = mu_cv, rep_err_model = (mu_cv-dur),
                         rre_model = (mu_cv-dur)/dur, sd_model =sigma_cv,
                         cv_model = sigma_cv/mu_cv)
    cv_model_mix = rbind2(cv_model_mix, cv_temp)
  }
}

# ---- 1-prior Model - figure_separate ----
# 4 figures from 4 exp using 1st model with 1 prior
# extracting the Bayesian data to integrate into the original figure/Yue
dat1_mix <- data.frame()

# Plot data + model predictions for visual and auditory conditions
for(i in c(1,3)){
  
  # Extract the part of the predictions associate with one condition
  pred_sep <- predictions %>% filter(exp == exp_name[i]) %>% group_by(x) %>%
    summarize(y_mdur=mean(y_mdur), y_vdur=mean(y_vdur), y_cv=mean(y_cv), y_rre=mean(y_rre))
  
  data_model = data %>% filter(exp == exp_name[i]) 
  
  mean_data <- group_by(data_model, duration) %>% summarize(mcv=mean(cv))
  Fig_cv_model<-
    ggplot(data_model, aes(x=duration, y=cv, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() +  scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10)) +
    geom_point(data=mean_data, aes(x=duration, y=mcv, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=mcv, group=1), size=1, color="black") +
    geom_line(data=pred_sep, aes(x=x, y=y_cv, group=1), color="blue", size=1) +
    labs(title = exp_name[i]) +
    xlab('Duration (s)') +
    ylab('CV')
  print(Fig_cv_model)
  
  Fig_cv_isub_model<-
    ggplot(data_model, aes(x=duration, y=cv, group=sub), color="black") + geom_point() +
    theme_bw() + 
    geom_line(data=filter(predictions, exp== exp_name[i]&x>0.2), aes(x=x, y=y_cv, group=1), color="blue", size=1) +
    scale_x_continuous(trans = "log10", breaks = c(0.3, 1,5,16))+
    facet_wrap(~sub) + 
    labs(title = exp_name[i])+
    xlab('Duration (s)') +
    ylab('CV')
  
  print(Fig_cv_isub_model)
  
  mean_data <- group_by(data_model, duration) %>% summarize(mmRRE=mean(mRepErr/duration))
  Fig_rre_model<-
    ggplot(data_model, aes(x=duration, y=mRepErr/duration, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() +  scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10)) +
    geom_point(data=mean_data, aes(x=duration, y=mmRRE, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=mmRRE, group=1), size=1, color="black") +
    geom_line(data=pred_sep, aes(x=x, y=y_rre, group=1), color="blue", size=1) +
    labs(title = exp_name[i])
  print(Fig_rre_model)
  
  Fig_rre_isub_model<-
    ggplot(data_model, aes(x=duration, y=mRepErr/duration, group=sub), color="black") + geom_point() +
    theme_bw() + 
    geom_line(data=filter(predictions, exp== exp_name[i]&x>0.2), aes(x=x, y=y_rre, group=1), color="blue", size=1) +
   # scale_x_continuous(trans = "log10", breaks = c(0.3, 1,5,16))+
    scale_y_continuous(labels = scales::percent)+
    facet_wrap(~sub) + 
    labs(title = exp_name[i])+
    xlab('Duration (s)') +
    ylab('RRE (%)')
  print(Fig_rre_isub_model)
  
  mean_data <- group_by(data_model, duration) %>% summarize(mmRepr=mean(mRepr))
  Fig_mdur_model<-
    ggplot(data_model, aes(x=duration, y=mRepr, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() + 
    geom_point(data=mean_data, aes(x=duration, y=mmRepr, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=mmRepr, group=1), size=1, color="black") +
    geom_line(data=pred_sep, aes(x=x, y=y_mdur, group=1), color="blue", size=1) +
    labs(title = exp_name[i]) + scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10)) +
    scale_y_continuous(trans = "log10")
  print(Fig_mdur_model)
  
  Fig_mdur_isub_model<-
    ggplot(data_model, aes(x=duration, y=mRepr, group=sub), color="black") + geom_point() +
    theme_bw() + 
    geom_line(data=filter(predictions, exp== exp_name[i]), aes(x=x, y=y_mdur, group=1), color="blue", size=1) +
    facet_wrap(~sub) + 
    labs(title = exp_name[i])
  print(Fig_mdur_isub_model)
  
  mean_data <- group_by(data_model, duration) %>% summarize(msdRepr=mean(sdRepr))
  Fig_vdur_model<-
    ggplot(data_model, aes(x=duration, y=sdRepr, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() + 
    geom_point(data=mean_data, aes(x=duration, y=msdRepr, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=msdRepr, group=1), size=1, color="black") +
    geom_line(data=pred_sep, aes(x=x, y=y_vdur, group=1), color="blue", size=1) +
    scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
    labs(title = exp_name[i])
  print(Fig_vdur_model)
  
  Fig_vdur_isub_model<-
    ggplot(data_model, aes(x=duration, y=sdRepr, group=sub), color="black") + geom_point() +
    theme_bw() + 
    geom_line(data=filter(predictions, exp== exp_name[i]), aes(x=x, y=y_vdur, group=1), color="blue", size=1) +
    facet_wrap(~sub) + 
    labs(title = exp_name[i])
  print(Fig_vdur_isub_model)
  
}


## ---- 3-prior model: take 'block' factor into consideration ----
# M_p + bias(block)
# M_p(block)
# each block requires one mean value
# ---- Build up 3-prior Models(Bayesian and linear) ----

# Same as Bayesian_mdur except that the mean of the prior and the bias to the mean of the 
# subjective prior (res) can differ between blocks
Bayesian_mdur_block <- function(t, block, sig_p, sig_t, res) {
  
  sp2 <- sig_p^2
  st2 <- sig_t^2

  y <- sp2/(sp2+st2)*t+st2/(sp2+st2)*(log(mean_block[block])+res[block])
  
}

# Same as the Bayesian_vdur for the mixed condition model. 
# Repeated here so that the code for the blocked condition model can be run on its own.
Bayesian_vdur <- function(sig_p, sig_t) {
  
  sp2 <- sig_p^2
  st2 <- sig_t^2
  
  y <- sqrt(sp2*st2/(sp2+st2))
  
}

# On log scale the model predicts constant SD. This function translates that into an SD on linear
# scale, i.e. to the SD of the log-normal distribution 
# (not the sigma parameter but the actual SD of the distibution).
# Same as for the mixed condition model. 
# Repeated here so that the code for the blocked condition model can be run on its own.
predicted_linear_scale_sd <- function(mu, sigma) {
  s <- sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2))
}

# Fit a straight line in log-space to get good starting parameters for the full fit
# The "model" parameter specificies which version of the model to run
# Different versions differ in whether sig_p, sig_s, sig_m and res can differ between short, 
# medium and long duration blocks as well as in whether the bias res is added as a shift of 
# the mean of the prior, or later as a bias of the reproduction
Bayesian_pre_fit_block <- function(par, model, data) {
  
  t <- log(data$duration)
  y <- log(data$reproduction)
  block <- data$block
  
  sp_sep <- bitwAnd(model,1) # Separate sig_p for odd-numbered models
  st_sep <- bitwAnd(bitwShiftR(model,1),1) # Separate sig_t 
  res_sep <- bitwAnd(bitwShiftR(model,2),1) # Separate res 
  outside_res <-bitwAnd(bitwShiftR(model,3),1)

  if(sp_sep) {
    sig_p <- par[1:3]
    par <- par[4:length(par)]
  } else {
    sig_p <- rep(par[1], 3)
    par <- par[2:length(par)]
  }
  
  if(st_sep) {
    sig_t <- par[1:3]
    par <- par[4:length(par)]
  } else {
    sig_t <- rep(par[1], 3)
    par <- par[2:length(par)]
  }
  
  if(res_sep) {
    res <- par[1:3]
  } else {
    res <- rep(par[1], 3)
  }
    
  sp <- sig_p[block]
  st <- sig_t[block]
    
  if(outside_res) {
    m_pred <- Bayesian_mdur_block(t, block, sp, st, c(0,0,0)) + res[block]
  } else {
    m_pred <- Bayesian_mdur_block(t, block, sp, st, res)
  }
  v_pred <- Bayesian_vdur(sp, st)
  
  a <- -sum(log(dnorm(y, m_pred, v_pred)))
  
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

# Full model fit for the blocked condition
Bayesian_mixed_fit_block <- function(par, model, data) {
  
  t <- log(data$duration)
  y <- data$reproduction
  block <- data$block
  
  sp_sep <- bitwAnd(model,1) # Separate sig_p for odd-numbered models
  st_sep <- bitwAnd(bitwShiftR(model,1),1) # Separate sig_t 
  res_sep <- bitwAnd(bitwShiftR(model,2),1) # Separate res 
  outside_res <-bitwAnd(bitwShiftR(model,3),1)
  sm_sep <-  bitwAnd(bitwShiftR(model,4),1) # Separate res   
  
  if(sp_sep) {
    sig_p <- par[1:3]
    par <- par[4:length(par)]
  } else {
    sig_p <- rep(par[1], 3)
    par <- par[2:length(par)]
  }
  
  if(st_sep) {
    sig_t <- par[1:3]
    par <- par[4:length(par)]
  } else {
    sig_t <- rep(par[1], 3)
    par <- par[2:length(par)]
  }
  
  if(res_sep) {
    res <- par[1:3]
    par <- par[4:length(par)]
  } else {
    res <- rep(par[1], 3)
    par <- par[2:length(par)]
  }
  
  if(sm_sep) {
    sig_m <- par[1:3]
  } else {
    sig_m <- rep(par[1], 3)
  }
  
  var_pred <- Bayesian_vdur(sig_p[block], sig_t[block])
  if(outside_res) {
    m_pred <- Bayesian_mdur_block(t, block, sig_p[block], sig_t[block], c(0,0,0))
  } else {
    m_pred <- Bayesian_mdur_block(t, block, sig_p[block], sig_t[block], res)
  }
  
  if(outside_res) {
    pred_m <- exp(m_pred + res[block] + var_pred^2/2)
  } else {
    pred_m <- exp(m_pred + var_pred^2/2) 
  }
  
  pred_sd <- sqrt(predicted_linear_scale_sd(m_pred, var_pred)^2 + sig_m[block]^2)
  
  a <- -sum(log(dnorm(y, pred_m, pred_sd)))
  if (is.na(a)) {
    a <- 1e9
  }
  if(a==Inf) {
    a <- 1e9
  }
  return(a)
}

#  ---- Initializing parameters ---- 
predictions_block <- data.frame()
fit_para_block <- data.frame()
cv_model_block <- data.frame()

# ----  Model fitting based on'Block' data  ----
for(i in c(2,4)){
  subs <- unique(filter(vdat, exp==exp_name[i])$sub)
  for(s in subs) {
    data_model = vdat %>% filter(exp == exp_name[i], sub==s) 
    N <- nrow(data_model)
    
    for(model in 0:31) {
      
      sp_sep <- bitwAnd(model,1) # Separate sig_p for odd-numbered models
      st_sep <- bitwAnd(bitwShiftR(model,1),1) # Separate sig_t 
      res_sep <- bitwAnd(bitwShiftR(model,2),1) # Separate res 
      outside_res <-bitwAnd(bitwShiftR(model,3),1)
      sm_sep <- bitwAnd(bitwShiftR(model,4),1) # Separate res   
      
      if(sp_sep) {
        start <- rep(0.6, 3)
        llim <- rep(0.01, 3)
        ulim <- rep(3, 3)
      } else  {
        start <- rep(0.6, 1)
        llim <- rep(0.01, 1)
        ulim <- rep(3, 1)
      }
      
      if(st_sep) {
        start <- c(start, rep(0.3, 3))
        llim <- c(llim, rep(0.01, 3))
        ulim <- c(ulim, rep(3, 3))
      } else  {
        start <- c(start, rep(0.3, 1))
        llim <- c(llim, rep(0.01, 1))
        ulim <- c(ulim, rep(3, 1))
      }
      
      if(res_sep) {
        start <- c(start, rep(0, 3))
        llim <- c(llim, rep(-10, 3))
        ulim <- c(ulim, rep(10, 3))
      } else  {
        start <- c(start, rep(0, 1))
        llim <- c(llim, rep(-10, 1))
        ulim <- c(ulim, rep(10, 1))
      } 
      
      ptm <- proc.time()
      par <- optim(start, Bayesian_pre_fit_block, NULL, method = "L-BFGS-B",
                   lower = llim, upper = ulim,
                   model, data_model, control = c(maxit=10000))
      
      conv1=par$convergence
      
      while(conv1>0) {
        par <- optim(start+runif(length(start),0,0.1), Bayesian_pre_fit_block, NULL, method = "L-BFGS-B",
                     lower = llim, upper = ulim,
                     model, data_model, control = c(maxit=10000))
        conv1=par$convergence
      }
      
      if(sm_sep) {
        start <- c(par$par, rep(0.1, 3))
        llim <- c(llim, rep(0, 3))
        ulim <- c(ulim, rep(3, 3))
      } else  {
        start <- c(par$par, rep(0.1, 1))
        llim <- c(llim, rep(0, 1))
        ulim <- c(ulim, rep(3, 1))
      } 
      npar <- length(start)
      
      par <- optim(start, Bayesian_mixed_fit_block, NULL, method = "L-BFGS-B",
                   lower = llim, upper = ulim,
                   model, data_model, control = c(maxit=10000))
      print(proc.time() - ptm)
      print(paste0(exp_name[i], ': sub=', s))
      print(par$par)
      print(par$value)
      print(par$convergence)
      print(par$message)
    
      conv2=par$convergence
      
      while(conv2>0) {
        par <- optim(start+runif(length(start),0,0.1), Bayesian_mixed_fit_block, NULL, method = "L-BFGS-B",
                     lower = llim, upper = ulim,
                     model, data_model, control = c(maxit=10000))
        conv2=par$convergence
      }
      
      p <- par$par
      if(sp_sep) {
        sig_p <- p[1:3]
        p <- p[4:length(p)]
      } else {
        sig_p <- rep(p[1], 3)
        p <- p[2:length(p)]
      }
      
      if(st_sep) {
        sig_t <- p[1:3]
        p <- p[4:length(p)]
      } else {
        sig_t <- rep(p[1], 3)
        p <- p[2:length(p)]
      }
      
      if(res_sep) {
        res <- p[1:3]
        p <- p[4:length(p)]
      } else {
        res <- rep(p[1], 3)
        p <- p[2:length(p)]
      }
      
      if(sm_sep) {
        sig_m <- p[1:3]
      } else {
        sig_m <- rep(p[1], 3)
      }
      
      w <- sig_p^2/(sig_p^2+sig_t^2)
      sig_tot <- Bayesian_vdur(sig_p, sig_t)
      
      fit_para_block <- rbind2(fit_para_block, data.frame(exp=exp_name[i], sub=s, model=model, sig_p_short=sig_p[1],
                                                          sig_p_medium=sig_p[2], sig_p_long=sig_p[3], sig_t_short=sig_t[1],
                                                          sig_t_medium=sig_t[2], sig_t_long=sig_t[3], sig_m_short=sig_m[1],
                                                          sig_m_medium=sig_m[1], sig_m_long=sig_m[1], res_short=res[1],
                                                          res_medium=res[2], res_long=res[3], w_short=w[1], w_medium=w[2], w_long=w[3], 
                                                          sig_tot_short=sig_tot[1], sig_tot_medium=sig_tot[2], sig_tot_long=sig_tot[3],
                                                          val=par$value, npar=npar, AIC=2*(par$value+npar), BIC=2*par$value+npar*log(N),
                                                          conv1=conv1, conv2=conv2, sp_sep=sp_sep, st_sep=st_sep, sm_sep=sm_sep, 
                                                          res_sep=res_sep, outside_res=outside_res))
      if(outside_res) {
        pred_m = Bayesian_mdur_block(log(x_block$m), x_block$block,
                                     sig_p[x_block$block], sig_t[x_block$block], rep(0, length(x_block$block))) 
        mu = exp(pred_m + res[x_block$block] + sig_tot[x_block$block]^2/2)
        sigma = sqrt(predicted_linear_scale_sd(pred_m, Bayesian_vdur(sig_p[x_block$block], sig_t[x_block$block]))^2 + sig_m[x_block$block]^2)
      } else {
        pred_m = Bayesian_mdur_block(log(x_block$m), x_block$block,
                                     sig_p[x_block$block], sig_t[x_block$block], res)
        mu = exp(pred_m + sig_tot[x_block$block]^2/2)
        sigma = sqrt(predicted_linear_scale_sd(pred_m, Bayesian_vdur(sig_p[x_block$block], sig_t[x_block$block]))^2 + sig_m[x_block$block]^2)
      }
      
      data_temp = data.frame(x=x_block$m, y_mdur = mu, 
                             y_vdur = sigma,
                             y_cv=sigma/mu,
                             y_rre = (mu - x_block$m)/x_block$m, model = model, exp = exp_name[i], sub=s)
      predictions_block <- rbind2(predictions_block,data_temp)
      
      while(model ==13){
        dur_block <- rep(c(1:3),each =3)
        
        # extract 9 points of predicted CVs, save in 'cv_model_block'
        
        pred_m_cv = Bayesian_mdur_block(log(dur), dur_block,
                                        sig_p[dur_block], sig_t[dur_block], rep(0, length(dur_block))) 
        mu_cv = exp(pred_m_cv + res[dur_block] + sig_tot[dur_block]^2/2)
        sigma_cv = sqrt(predicted_linear_scale_sd(pred_m_cv, Bayesian_vdur(sig_p[dur_block], sig_t[dur_block]))^2 + sig_m[dur_block]^2)
        
        
        # add experimental data into the data frame
        data_model_mean <- data_model %>% group_by(exp, sub, duration)  %>% 
          summarize(rep = mean(reproduction), rep_err = mean(rep_err),
                    sd = sd(reproduction), cv = sd/rep)
        
        cv_temp = data.frame(exp = exp_name[i], sub=s,duration = dur, 
                             rep = data_model_mean$rep, rep_err = data_model_mean$rep_err,
                             sd = data_model_mean$sd, cv = data_model_mean$cv,
                             rep_model = mu_cv, rep_err_model = (mu_cv-dur),
                             rre_model = (mu_cv-dur)/dur, sd_model =sigma_cv,
                             cv_model = sigma_cv/mu_cv)
        cv_model_block = rbind2(cv_model_block, cv_temp) 
        model <-0 }
      
      
      
    }
  }
}

pred_sep_block <-  predictions_block %>% group_by(exp, x, model) %>%
  summarize(y_mdur=mean(y_mdur), y_vdur=mean(y_vdur), y_cv=mean(y_cv), y_rre=mean(y_rre))
fit_summary <- group_by(fit_para_block, model) %>% summarize(mBIC=mean(BIC)) %>% arrange(mBIC)
top_mod <- fit_summary[1,]$model

# Plot figures with data and model predictions 
for(i in c(2,4)) {
  
  data_model = data %>% filter(exp == exp_name[i])
  pred_sep = pred_sep_block %>% filter(exp == exp_name[i], model==top_mod)

  mean_data <- group_by(data_model, duration) %>% summarize(mcv=mean(cv))
  Fig_cv_model<-
    ggplot(data_model, aes(x=duration, y=cv, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() +  scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10)) +
    geom_point(data=mean_data, aes(x=duration, y=mcv, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=mcv, group=1), size=1, color="black") +
    geom_line(data=pred_sep, aes(x=x, y=y_cv, group=1), color="blue", size=1) +
    labs(title = exp_name[i])
  print(Fig_cv_model)
  
  Fig_cv_isub_model<-
    ggplot(data_model, aes(x=duration, y=cv, group=sub), color="black") + geom_point() +
    theme_bw() + 
    geom_line(data=filter(predictions_block, exp== exp_name[i], model==top_mod), aes(x=x, y=y_cv, group=1), color="blue", size=1) +
    facet_wrap(~sub) + 
    scale_x_continuous(trans = "log10", breaks = c(0.3, 1,5,16))+
    labs(title = exp_name[i])+
    xlab('Duration (s)') +
    ylab('CV')
  print(Fig_cv_isub_model)

  mean_data <- group_by(data_model, duration) %>% summarize(mmRRE=mean(mRepErr/duration))
  Fig_rre_model<-
    ggplot(data_model, aes(x=duration, y=mRepErr/duration, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() +  scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10)) +
    geom_point(data=mean_data, aes(x=duration, y=mmRRE, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=mmRRE, group=1), size=1, color="black") +
    geom_line(data=pred_sep, aes(x=x, y=y_rre, group=1), color="blue", size=1) +
    labs(title = exp_name[i])
  print(Fig_rre_model)
  
  Fig_rre_isub_model<-
    ggplot(data_model, aes(x=duration, y=mRepErr/duration, group=sub), color="black") + geom_point() +
    theme_bw() + 
    scale_y_continuous(labels = scales::percent)+
    geom_line(data=filter(predictions_block, exp== exp_name[i], model==top_mod), aes(x=x, y=y_rre, group=1), color="blue", size=1) +
    facet_wrap(~sub) + 
    labs(title = exp_name[i])+
    xlab('Duration (s)') +
    ylab('RRE (%)')
  print(Fig_rre_isub_model)
  
  mean_data <- group_by(data_model, duration) %>% summarize(mmRepr=mean(mRepr))
  Fig_mdur_model<-
    ggplot(data_model, aes(x=duration, y=mRepr, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() + 
    geom_point(data=mean_data, aes(x=duration, y=mmRepr, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=mmRepr, group=1), size=1, color="black") +
    geom_line(data=pred_sep, aes(x=x, y=y_mdur, group=1), color="blue", size=1) +
    scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
    labs(title = exp_name[i])
  print(Fig_mdur_model)
  
  Fig_mdur_isub_model<-
    ggplot(data_model, aes(x=duration, y=mRepr, group=sub), color="black") + geom_point() +
    theme_bw() + 
    geom_line(data=filter(predictions_block, exp== exp_name[i], model==top_mod), aes(x=x, y=y_mdur, group=1), color="blue", size=1) +
    facet_wrap(~sub) + 
    labs(title = exp_name[i])
  print(Fig_mdur_isub_model)
  
  mean_data <- group_by(data_model, duration) %>% summarize(msdRepr=mean(sdRepr))
  Fig_vdur_model<-
    ggplot(data_model, aes(x=duration, y=sdRepr, color=sub, group=sub)) + geom_point() +
    geom_line() + theme_bw() + 
    geom_point(data=mean_data, aes(x=duration, y=msdRepr, group=1), size=3, color="black") + 
    geom_line(data=mean_data, aes(x=duration, y=msdRepr, group=1), size=1, color="black") +
    geom_line(data=pred_sep, aes(x=x, y=y_vdur, group=1), color="blue", size=1) +
    scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
    labs(title = exp_name[i])
  print(Fig_vdur_model)
  
  Fig_vdur_isub_model<-
    ggplot(data_model, aes(x=duration, y=sdRepr, group=sub), color="black") + geom_point() +
    theme_bw() + 
    geom_line(data=filter(predictions_block, exp== exp_name[i], model==top_mod), aes(x=x, y=y_vdur, group=1), color="blue", size=1) +
    facet_wrap(~sub) + 
    labs(title = exp_name[i])
  print(Fig_vdur_isub_model)
  
}
