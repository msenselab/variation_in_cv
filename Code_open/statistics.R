## After running 'Main.R' for relevant datasets.
# ---- ttest function against 0 ----
## t-tests function, individual mean slope against 0
t_test <- function(df, u=0) {
  tmp <- t.test(df$slope,mu=u)
  stats <-  c(tmp$parameter,tmp$statistic,tmp$estimate,tmp$p.value)
  names(stats) <- c('df','t','mean','p')
  data.frame (as.list(stats))
}

# ---- ttest function among sessions ----
# group ttest function
ttest_between <- function(df1, df2) {
  tmp <- t.test(df1$regindex,df2$regindex)
  stats <-  c(tmp$parameter,tmp$statistic,tmp$estimate,tmp$p.value)
  names(stats) <- c('df','t','mean1', 'mean2','p')
  data.frame (as.list(stats))
}

# ---- (Figure 4) regression index (RI) from data_mean_slope ----
data_mean_slope$regindex = -data_mean_slope$regindex
d1 <- filter(data_mean_slope, exp == 'Vis/Mix')
d2 <- filter(data_mean_slope, exp == 'Vis/Block')
d3 <- filter(data_mean_slope, exp == 'Aud/Mix')
d4 <- filter(data_mean_slope, exp == 'Aud/Block')

# group regression index ttest
ri12 <- ttest_between(d1,d2)
ri34 <- ttest_between(d3,d4)
ri13 <- ttest_between(d1,d3)
ri24 <- ttest_between(d2,d4)

# ---- (Figure 6) mean CV slope ttest ----
# group ttest function
ttest_cv <- function(df1, df2) {
  tmp <- t.test(df1$slope,df2$slope)
  stats <-  c(tmp$parameter,tmp$statistic,tmp$estimate,tmp$p.value)
  names(stats) <- c('df','t','mean1', 'mean2','p')
  data.frame (as.list(stats))
}
# based on data 'data_cv_slope'
d1_cv <- filter(data_cv_slope, exp == 'Vis/Mix')
d2_cv <- filter(data_cv_slope, exp == 'Vis/Block')
d3_cv <- filter(data_cv_slope, exp == 'Aud/Mix')
d4_cv <- filter(data_cv_slope, exp == 'Aud/Block')

# group regression index ttest
cv12 <- ttest_cv(d1_cv,d2_cv)
cv34 <- ttest_cv(d3_cv,d4_cv)
cv13 <- ttest_cv(d1_cv,d3_cv)
cv24 <- ttest_cv(d2_cv,d4_cv)

#cv  against slope = 0
cv1 <- t_test(d1_cv)
cv2 <- t_test(d2_cv)
cv3 <- t_test(d3_cv)
cv4 <- t_test(d4_cv)

## ---- (Table 1)s tatistics on model parameters ----
# parameters from data:
# Mix: Bayes_isub_fit_pata
# Block: bayes_block_fit_para_4
# mean
data_model_para_mix <- bayes_isub_fit_para %>% group_by(exp) %>%
  summarise(mean_sigm = mean(sig_m),mean_sigp = mean(sig_p),mean_sigt = mean(sig_t),mean_res = mean(res))

data_model_para_block <- bayes_block_fit_para_4 %>% group_by(exp) %>%
  summarise_all(mean, na.rm = TRUE)
# plot mean parameters as tables
library(apaTables)

## ---- ttest ----
## mix condition
t1<- bayes_isub_fit_para%>% filter(exp == 'Vis/Mix')
t2<- bayes_isub_fit_para %>% filter(exp == 'Aud/Mix')
t.test(t1$sig_m,t2$sig_m)
# sig_m:         0.004 **
# sig_t(sig_s) : 0.02 *
# sig_p :     non-sig
# res:        non-sig

## Block condition : only sig_m
# data from: bayes_block_fit_para_r
t11<- bayes_block_fit_para_4%>% filter(exp == 'Vis/Block')
t22<- bayes_block_fit_para_4 %>% filter(exp == 'Aud/Block')
t.test(t11$sig_m,t22$sig_m)

# other 3 parameters in 3 conditions
# (Vis, Aud) * (s,m,l)
t.test(t11$res_short,t22$res_short)
t.test(t11$res_medium,t22$res_medium)
t.test(t11$res_long,t22$res_long)

## between condition
t.test(t1$sig_m,t11$sig_m)
t.test(t2$sig_m,t22$sig_m)

# ANOVA for sig_r(sig_m)
# prepare data 
data_sig_m <- rbind(bayes_isub_fit_para[c(1,2,3)], bayes_block_fit_para_4[c(1,2,3)]) %>%
  separate(.,exp,c('Modality','Condition'),sep='/', remove = FALSE) 
# ANOVA
anova.sig_m <-  aov(sig_m ~ Modality*Condition + Error(sub), data = data_sig_m) %>% tidy(.) 
# using ezanova
library(apaTables)
library(ez)
ezanova_sig_m <- ezANOVA(data = data_sig_m,
                         dv = sig_m,
                         between = .(Modality, Condition),
                         wid = sub ,
                         detailed = TRUE)


# Make APA table and save to png file
#ezanova_sig_m_table <- apa.ezANOVA.table(ezanova_sig_m,
#                                        filename="sig_m_ez_table.png")



