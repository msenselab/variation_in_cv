# Main Codes for figure plots from Figs 3 to 6
## ----- load rawdata and install necessary packages ----
source('loadPackages.R')
# change the colour pairs used in the bar figures
myPairs <- brewer.pal(9, "Paired")[c(6,2,5,1)]
# load raw data (after filtering out outliers)
vdat = read.csv("vdata.csv", header = TRUE)
## ---- prepare mean data ----
# average means subject-wise 
msRepr <- vdat %>% group_by(exp, duration,sub) %>% 
  summarise(mRepr = mean(reproduction), sdRepr = sd(reproduction), 
            mPrErr = mean(prod_err,na.rm=TRUE), sdProd = sd(prod_err,na.rm=TRUE), 
            mRepErr = mean(rep_err), sdRepErr = sd(rep_err),
            rre = mean(rep_err/duration))  %>% 
  mutate(cv = sdRepr/mRepr, log_dur = log2(duration)) %>% # cv based on physical duration
  arrange(exp, duration,sub) %>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE)
# msRepr$block <- rep(1:3)
msRepr$block <- c(rep(1:3,each = 12*3), rep(1:3,each =12*3), rep(1:3, each =11*3), rep(1:3, each =14*3))
msRepr$block <- factor(msRepr$block)

# grand means (collapse subjects)
mRepr <- msRepr %>% group_by(exp, duration) %>% 
  summarise(mRepr = mean(mRepr),  mCV = mean(cv),n = n(),
            msdProd = mean(sdProd), msdRepr = mean(sdRepr),
            mPrErr = mean(mPrErr), mRE = mean(mRepErr), 
            mRRE = mean(mRepErr/duration),mseRRE = sd(mRepErr/duration)/sqrt(n-1))  %>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE)


# ---- load modeling data ----
source('Bayesian_cv_model_log.R')

# ---- Figure 3 Mean RREs ---- ----
# Relative Reproduction Errors(RREs)
library("lemon")
fig_rre <- ggplot(mRepr, aes(x=duration, y = mRRE, colour= design))+
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mRRE-mseRRE, ymax = mRRE+mseRRE), width = 0.1) + 
  # individual points
  #geom_point(data = msRepr, aes(x=duration, y = rre, colour = design), size = 0.5) +
  scale_x_continuous(trans = "log10", limits = c(0.25,30),breaks = c(0.3, 1,5,10,16))+
  geom_abline(aes(intercept = 0, slope = 0),linetype = 9)+
  # add modeling data
  geom_line(data = dat1_mix %>% separate(.,exp,c('modality','design'),sep='/', remove = FALSE) , aes(x = x, y = y_rre),color = "cyan", size = 0.35)+
  geom_line(data = dat1_block %>% separate(.,exp,c('modality','design'),sep='/', remove = FALSE), aes(x = x, y = y_rre),color = "red",size = 0.35)+
  ylab("Relative Reproduction Error (RRE) (%)")+
  xlab("Duration (s)")+
  scale_y_continuous(labels = scales::percent, limits = c(-0.6,1.6),breaks = c(-0.5,0,0.5,1)) + 
  facet_wrap(~modality,nrow=1)+
  theme_minimal()+
  theme(legend.position = c(0.1,0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        # axis.text.y = element_blank(),
        axis.ticks = element_line(color = "grey"))+
  labs(colour = "Condition") 
 #coord_flip()

# ---- Figure 4: Regression Index(TI) ----
# regression index
data_mean_slope <- slope_estimation(msRepr,'mRepr ~ duration')
data_mean_slope$regindex <- 1-data_mean_slope$slope

fig_ri <- data_mean_slope %>% 
  group_by(exp) %>% 
  summarise(mslope = mean(slope),
            mregindex = mean(regindex),
            seregindex = sd(regindex)/sqrt(n())) %>% 
  ggplot(., aes(x= exp, y = mregindex, fill = exp),
         draw_key = draw_key_polygon) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.6) +
  geom_errorbar(aes(ymin = mregindex-seregindex, ymax = mregindex+seregindex,width = 0.6, colour = exp))+
  theme_classic() +
  legend_pos('null') +
  scale_colour_manual(values = myPairs) +
  scale_fill_manual(values = myPairs) +
  geom_signif(comparisons = list(c("Vis/Mix","Vis/Block"), size = 0.6),
              annotations = "***",
              y_position = 0.5) +  
  geom_signif(comparisons = list(c("Aud/Mix","Aud/Block")),
              annotations = "*",
              y_position = 0.45) +
  ylim(0,0.65)+
  ylab('Regression Index') +
  xlab('Condition')
# ---- Figure 5: CVs ----
data = msRepr
data_m = data %>% 
  group_by(exp,duration) %>%
  summarise(mcv = mean(cv),
            mblock = mean(block))%>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE) 

dat1_mix <- dat1_mix %>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE)

dat1_block <- dat1_block%>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE)

Fig_cv_model <-
  ggplot(data= data, aes(x= duration, y = cv, colour= design))+  
  scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10))+
  scale_y_continuous(limits = c(0.05,0.75))+
  #invididual points
  geom_point(size = 0.2,alpha = 0.6)+
  #Byesian model lines
  geom_line(aes(x=x,y=y_cv),data=dat1_mix %>% filter(x>0.2),size = 0.35)+ 
  geom_line(aes(x=x,y=y_cv),data=dat1_block %>% filter(x>0.2),size = 0.35)+ 
  # mean points
  #geom_line(data =data_m , aes(x= duration, y= mcv), color= design,linetype = design)+
  geom_point(data =data_m, aes(x= duration, y= mcv), size = 1.2) +
  # legend_pos(c(0.8,0.8))+
  facet_rep_wrap(~modality, nrow = 1)+
  xlab("Duration (s)") +
  ylab('Coefficient of Variance')+
  theme_minimal()+
  theme(legend.position = c(0.1,0.8),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        axis.ticks = element_line(color = "grey")
  )


# ---- Figure 6: CV slopes ----
data_cv_slope <- slope_estimation(msRepr,'cv ~ log_dur')
data_cv_mslope <- data_cv_slope %>% 
  group_by(exp) %>% 
  summarise(mslope = mean(slope),
            seslope = sd(slope)/sqrt(n()))

Fig_meanslopes <- ggplot( data_cv_mslope, aes(x =exp, y = mslope, fill = exp)) +
  geom_bar( stat = "identity",width = 0.6) +
  geom_errorbar(aes(ymin = mslope-seslope,ymax = 0,width = 0.6, colour = exp))+
  geom_signif(comparisons = list(c("Vis/Mix","Vis/Block")),
              annotations = "*",
              y_position = 0.01) +  
  geom_signif(comparisons = list(c("Vis/Mix","Aud/Mix")),
              annotations = "*",
              y_position = 0.015) +
  theme_classic() +
  legend_pos("null") +
  scale_colour_manual(values = myPairs) +
  scale_fill_manual(values = myPairs) +
  xlab('Condition') +
 # ylim(c(0,0.057))+
  ylab('CV Slopes ')


# ---- (New) Figure 2 ----
fig2_new <- plot_grid(fig_rre, fig_ri, nrow =1,rel_widths = c(2,1), labels = c('A','B'))
# ---- (New) Figure 3 ----
fig3_new <- plot_grid(Fig_cv_model, Fig_meanslopes, nrow =1,rel_widths = c(2,1), labels = c('A','B'))
