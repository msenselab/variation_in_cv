# Main Codes of figure plots for Figure 2 and 3.
# ----- Install all necessary packages and functions ----
source('loadPackages.R')
# load data from 4 experiments  
vdat = read.csv("vdata.csv", header = TRUE)
## ---- prepare grand-mean data ----
# Grand-mean data from individual participants
# ...were stored in the data.frame 'msRepr'
# Each row of the mean or sd was calculated from 30 repetitions, as each duration
# ...from each participant was tested 30 times.
# The size of the data frame is 441*15, 
# ...in which the row number represents 49(subject number) * 9 (sample duration number);
# ...and each column represents properties that elaberated on in the following:
# The 1st column marked the condition of 4 experiments: 'Vis/Mix', 'Vis/Block', 'Aud/Mix', 'Aud/Block'.
# The 2nd and 3rd columns marked the Modality ('Vis' and 'Aud') 
# ...and the Design ('Mix' and 'Block') of 4 experiments, repectively.
# The 4th column marked the subject number (49 in total).
# The 5th column marked the mean statistical properties collapsed within each subject:
# (column) 6-mean reproduction ; 7-standard devidation of reproduction;
# 8-mean production error; 9-standard deviation of production'
# 10-mean reproduction error; 11-standard deviation of reproduction error;
# 12-relative reproduction error; 13-cofficient of variance;
# The 14th column marked the logarithmic transformation of the sample duraiton.
# The 15th column marked the block number that were in accordance with the temporal range.
# In the 'Mix' conditions (Experiments 2 and 4), all the block numbers were set as 1.

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


# ---- load modeling data ----
source('Bayesian_model_final.R')

# extract and average data from model fittings
# 'predictions' -Mix and 'predictions_block' - Block
predictions_mean <- predictions %>% group_by(exp,x) %>%
  summarize (y_mdur = mean(y_mdur), y_vdur = mean(y_vdur),
             y_rre = mean(y_rre), y_cv = mean(y_cv))%>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE)

predictions_block_mean <- predictions_block %>% filter(model==top_mod) %>% group_by(exp,x) %>%
  summarize (y_mdur = mean(y_mdur), y_vdur = mean(y_vdur),
             y_rre = mean(y_rre), y_cv = mean(y_cv))%>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE)

# ---- Figure 2: Mean RREs ---- ----
# Relative Reproduction Errors(RREs)
# Configuration of subplot titles
rre.labs <- as_labeller( c('Aud'='Audition', 'Vis'='Vision'))

## Figure 2
#Relative Reproduction Errors (RREs) (grey and black dots) and predicted RREs (dark and light red lines) 
# by the two-stage model, as a function of the sample interval. The left panel depicts the auditory sessions, 
# the right panel the visual sessions. Black hollow dots represent estimates from ‘blocked’ conditions and 
# ...grey solid dots represent ‘mixed’ conditions. 

Figure2 <-
  ggplot(mRepr, aes(x=duration, y = mRRE))+
  geom_point(aes(alpha = Design, shape = Design),size = 2) +
  scale_alpha_manual(values = c(0.4,1))+
  geom_errorbar(aes(ymin = mRRE-mseRRE, ymax = mRRE+mseRRE), width = 0.1) + 
  scale_shape_manual(values = c(4,20))+ 
  # Modeling data
  # Mix data
  geom_line(aes(x=x,y=y_rre),  color = "black",data = predictions_mean,linetype = "dashed")+ 
  facet_wrap(~modality,nrow=1,labeller=rre.labs)+
  #scale_color_manual(values = myPairs[c(2,4)], guide = FALSE)+
  # Block data
  geom_line(aes(x=x,y=y_rre),data=predictions_block_mean %>% filter(x>0.2&x<0.9), size =0.5, alpha = 0.7)+ 
  geom_line(aes(x=x,y=y_rre),data=predictions_block_mean %>% filter(x>1&x<3.7), size =0.5, alpha = 0.7)+ 
  geom_line(aes(x=x,y=y_rre),data=predictions_block_mean %>% filter(x>5),  size =0.5, alpha = 0.7)+ 
  scale_x_continuous(trans = "log10", limits = c(0.25,30),breaks = c(0.3, 1,5,16))+
  geom_abline(aes(intercept = 0, slope = 0),linetype = 9)+
  ylab("Relative Reproduction Error (RRE) (%)")+
  xlab("Durations (s)")+
  scale_y_continuous(labels = scales::percent, limits = c(-0.6,1),breaks = c(-0.5,0,0.5)) + 
  facet_wrap(~modality,nrow=1,labeller=rre.labs)+
  theme_minimal()+
  theme(legend.position = c(0.15,0.15),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        axis.ticks = element_line(color = "grey")) 

# ---- Figure 3 : CVs ----
# configurations of the subplot titles
cv.labs <- as_labeller( c('Aud'='Audition', 'Vis'='Vision'))

# Figure 3 A
figure3_cv <-
  ggplot(data= mRepr, aes(x= duration, y = mCV))+  
  # mean points
  geom_point(aes(alpha = Design, shape = Design), size = 2)+
  scale_alpha_manual(values = c(0.4,1))+
  
  geom_errorbar(aes(ymin = mCV-seCV, ymax = mCV+seCV, 
                    alpha= Design), width = 0.1) + 
  scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10))+
  scale_y_continuous(limits = c(0.05,0.65))+
  scale_shape_manual(values = c(4,20))+
  facet_rep_wrap(~modality, nrow = 1,labeller=cv.labs)+
  #Byesian model lines
  # Mix data
  geom_line(aes(x=x,y=y_cv),data=predictions_mean %>% filter(x>0.2), linetype = "dashed")+ 
  # scale_color_manual(values = myPairs[c(2,4)], guide=FALSE)+
  # block data
  geom_line(aes(x=x,y=y_cv),data=predictions_block_mean %>% filter(x>0.2&x<0.9),  size =0.5, alpha = 0.6)+ 
  geom_line(aes(x=x,y=y_cv),data=predictions_block_mean %>% filter(x>1&x<3.7), size =0.5, alpha = 0.6)+ 
  geom_line(aes(x=x,y=y_cv),data=predictions_block_mean %>% filter(x>5),size =0.5, alpha = 0.6)+ 
  xlab('Duartions (s)') +
  ylab('Coefficient of Variance (CV)')+
  theme_minimal()+
  theme(legend.position = c(0.15,0.85),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        axis.ticks = element_line(color = "grey")
  )

# Figure 3 B
# prepare dataset
data_cv_slope <- slope_estimation(msRepr,'cv ~ log_dur')
data_cv_mslope <- data_cv_slope %>% 
  group_by(exp) %>% 
  summarise(mslope = mean(slope),
            seslope = sd(slope)/sqrt(n())) %>%
  setnames(., "exp","Experiment")


# Figure plot
Fig_meanslopes <- ggplot(data_cv_mslope, aes(x =Experiment, y = mslope)) +
  
  geom_bar( stat = "identity",width = 0.6, fill = "black") +
  geom_text(aes(x= 'Vis/Mix', y = 0.005,label = "**")) +
  
  geom_errorbar(aes(ymin = mslope-seslope,ymax = 0),width = 0.4)+
  scale_y_continuous(labels = scales::percent_format(accuracy =0.01))+  
  theme_minimal()+
  
  theme(legend.position = "null",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #  axis.text.x=element_blank(),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        axis.ticks = element_line(color = "grey"))+
  xlab('Conditions') +
  ylim(c(-0.05,0.025))+
  ylab('CV Slopes ')

Fig_meanslopes_mid <- plot_grid(NULL,Fig_meanslopes, NULL, nrow =3,rel_heights = c(2,5,1))

# Figure 3
# panels A: Mean CVs (colored dots) as a function of the sample duration, 
# separately for the auditory (left panel) and visual (right panel) sessions. 
# The smaller colored dots along the vertical line centered on each sample duration represent the CVs of individual participants. 
# The colored lines represent model fittings using the proposed Bayesian-Estimator model. 
# Red dots (lines) indicate the ‘block(ed)’ condition and cyan dots (lines) the ‘mix(ed)’ condition. 
# panel B: Mean Slopes of the CV as a linear function of logarithmic duration, for the four experiments. 
# Error bars indicate one standard error. 
# The slope was obtained by estimating parameter b of the linear function CV = a + b(Duration)

Figure3 <- plot_grid(figure3_cv, Fig_meanslopes_mid, nrow =1,rel_widths = c(2,1), labels = c('A','B'))


# ---- Figure 4 : model v.s. estimates ----
## prepare data from modeling
# 'Mix' conditions
fit_para %>% group_by (exp) %>% 
  summarize (mean_sig_p = mean(sig_p), sd_sig_p = sd(sig_p),
             mean_sig_t = mean(sig_t), sd_sig_t = sd(sig_t),
             mean_sig_m = mean(sig_m), sd_sig_m = sd(sig_m),
             mean_w = mean(w), sd_w = sd(w),
             mean_res1 = mean(res), sd_res1 = sd(res),
             mean_res2 = mean(res2), sd_res2 = sd(res2)) -> para_mix_mean
# 'Block' conditions
fit_para_block %>% filter(model == top_mod) -> fit_para_block_top
fit_para_block_top %>% group_by (exp) %>% 
  summarize (mean_sig_p1 = mean(sig_p_short), sd_sig_p1 = sd(sig_p_short),
             mean_sig_p2 = mean(sig_p_medium), sd_sig_p2 = sd(sig_p_medium),
             mean_sig_p3 = mean(sig_p_long), sd_sig_p3 = sd(sig_p_long),
             mean_sig_t = mean(sig_t_short), sd_sig_t1 = sd(sig_t_short),
             mean_sig_r = mean(sig_m_short), sd_sig_r1 = sd(sig_m_short),
             mean_w1 = mean(w_short), sd_w1 = sd(w_short),
             mean_w2 = mean(w_medium), sd_w2 = sd(w_medium),
             mean_w3 = mean(w_long), sd_w3 = sd(w_long),
             mean_res1 = mean(res_short), sd_res1 = sd(res_short),
             mean_res21 = mean(res2_short), sd_res21 = sd(res2_short),
             mean_res22 = mean(res2_medium), sd_res22 = sd(res2_medium),
             mean_res23 = mean(res2_long), sd_res23 = sd(res2_long),
             mean_aic = mean(AIC), mean_bic = mean(BIC)) -> para_block_mean

data_model_exp <- rbind(cv_model_mix, cv_model_block)
data_model_exp$rre <- data_model_exp$rep_err/data_model_exp$duration

data_model_exp_mean <- data_model_exp %>% group_by(exp,duration) %>%
  summarize(n = n(), cv_se = sd(cv)/sqrt(n-1),cv = mean(cv), 
            cv_model_se = sd(cv_model)/sqrt(n-1),cv_model = mean(cv_model),
            rep_se = sd(rep)/sqrt(n-1),rep = mean(rep), 
            rep_model_se = sd(rep_model)/sqrt(n-1), rep_model = mean(rep_model),
            sd_se = sd(sd)/sqrt(n-1),sd = mean(sd), 
            sd_model_se = sd(sd_model)/sqrt(n-1), sd_model = mean(sd_model),
            rep_err_se = sd(rep_err)/sqrt(n-1),rep_err = mean(rep_err), 
            rep_err_model_se = sd(rep_err_model)/sqrt(n-1), rep_err_model = mean(rep_err_model),
            rre_se = sd(rre)/sqrt(n-1),rre = mean(rre),
            rre_model_se = sd(rre_model)/sqrt(n-1),rre_model = mean(rre_model)
  )

Figure4 <- ggplot()+
  geom_point(data=data_model_exp, aes(x=rre, y = rre_model),color = "grey")+
  # add mean points
  geom_point(data = data_model_exp_mean, aes(x= rre, y= rre_model),
             size = 2, color = "Black", shape = 20)+
  geom_errorbar(data = data_model_exp_mean,
                aes(x= rre, ymin=rre_model - rre_model_se, 
                    ymax= rre_model+ rre_model_se),
                width= 0.05, color ="black")+
  geom_errorbarh(data = data_model_exp_mean, 
                 aes(y= rre_model, xmin=rre-rre_se, 
                     xmax= rre + rre_se),
                 height= 0.01,color = "black")+   
  theme_minimal()+
  theme(legend.position = "NULL",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        axis.ticks = element_line(color = "grey")
  )+
  xlab('RREs (%)')+
  ylab('Predicted RREs (%)')+
  scale_colour_manual(values = myPairs) +
  scale_x_continuous(labels = scales::percent, limits = c(-1,1.5)) + 
  scale_y_continuous(labels = scales::percent, limits = c(-1,1.5)) + 
  # add diagonal lines
  geom_abline(slope = 1, linetype = 'dashed')+
  facet_wrap(~exp)


