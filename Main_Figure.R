# Main Codes of figure plots from Figs 3 to 6
## ----- install necessary packages ----
source('Codes/loadPackages.R')
# load data from 4 experiments  
vdat = read.csv("Code_open/vdata.csv", header = TRUE)
## ---- prepare grand-mean data ----

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
source('Bayesian_model_final_local.R')

# extract and average data from model fittings
# 'predictions' -Mix and 'predictions_block' - Block
predictions_mean <- predictions %>% group_by(exp,x) %>%
  summarize (y_mdur = mean(y_mdur), y_vdur = mean(y_vdur),
             y_rre = mean(y_rre), y_cv = mean(y_cv))%>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE)

predictions_block_mean <- predictions_block %>% group_by(exp,x) %>%
  summarize (y_mdur = mean(y_mdur), y_vdur = mean(y_vdur),
             y_rre = mean(y_rre), y_cv = mean(y_cv))%>%
  separate(.,exp,c('modality','design'),sep='/', remove = FALSE)

# ---- Figure 2 Mean RREs ---- ----
# Relative Reproduction Errors(RREs)
# Configuration of subplot titles

rre.labs <- as_labeller( c('Aud'='Audition', 'Vis'='Vision'))
# Figure 2
Figure2_rre <-
  ggplot(mRepr, aes(x=duration, y = mRRE,color = Design))+
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mRRE-mseRRE, ymax = mRRE+mseRRE), width = 0.1) + 
  scale_color_manual(values = c("grey43","black"))+
  scale_shape_manual(values = c(1,16))+ 
  #Modeling data
  geom_line(aes(x=x,y=y_rre), data = predictions_mean,color = "#E31A1C")+ 
  facet_wrap(~modality,nrow=1,labeller=rre.labs)+
  
  # 3 block data
  geom_line(aes(x=x,y=y_rre),data=predictions_block_mean %>% filter(x>0.2&x<0.9),color = "#FB9A99",  size =0.5)+ 
  geom_line(aes(x=x,y=y_rre),data=predictions_block_mean %>% filter(x>1&x<3.7), color = "#FB9A99", size =0.5)+ 
  geom_line(aes(x=x,y=y_rre),data=predictions_block_mean %>% filter(x>5), color = "#FB9A99",  size =0.5)+ 
  
  scale_x_continuous(trans = "log10", limits = c(0.25,30),breaks = c(0.3, 1,5,16))+
  geom_abline(aes(intercept = 0, slope = 0),linetype = 9)+
  ylab("Relative Reproduction Error (RRE) (%)")+
  xlab("Durations (s)")+
  scale_y_continuous(labels = scales::percent, limits = c(-0.6,1),breaks = c(-0.5,0,0.5)) + 
  facet_wrap(~modality,nrow=1,labeller=rre.labs)+
  theme_minimal()+
  theme(legend.position = c(0.15,0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        axis.ticks = element_line(color = "grey")) 


# ---- Figure 3 (A+B): CVs ----
# configurations of the subplot titles
cv.labs <- as_labeller( c('Aud'='Audition', 'Vis'='Vision'))

# Figure 3
Figure3_cv <-
  ggplot(data= mRepr, aes(x= duration, y = mCV, colour= Design))+  
  #mean points
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = mCV-seCV, ymax = mCV+seCV, 
                    colour= Design), width = 0.1) + 
  scale_x_continuous(trans = "log10",breaks = c(0.3,1,2,5,10))+
  scale_y_continuous(limits = c(0.05,0.65))+
  scale_color_manual(values = c("grey43","black"))+
  scale_shape_manual(values = c(1,16))+
  facet_rep_wrap(~modality, nrow = 1,labeller=cv.labs)+
  #Byesian model lines
  geom_line(aes(x=x,y=y_cv),data=predictions_mean %>% filter(x>0.2), color = "#E31A1C")+ 
  # 3 block data
  geom_line(aes(x=x,y=y_cv),data=predictions_block_mean %>% filter(x>0.2&x<0.9), color = "#FB9A99",  size =0.5)+ 
  geom_line(aes(x=x,y=y_cv),data=predictions_block_mean %>% filter(x>1&x<3.7), color = "#FB9A99", size =0.5)+ 
  geom_line(aes(x=x,y=y_cv),data=predictions_block_mean %>% filter(x>5), color = "#FB9A99",  size =0.5)+ 
  
  # legend_pos(c(0.8,0.8))+
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

# ---- Figure: CV slopes ----

data_cv_slope <- slope_estimation(msRepr,'cv ~ log_dur')
data_cv_mslope <- data_cv_slope %>% 
  group_by(exp) %>% 
  summarise(mslope = mean(slope),
            seslope = sd(slope)/sqrt(n())) %>%
  setnames(., "exp","Experiment")

myPairs <- brewer.pal(9, "Paired")[c(6,2,5,1)]

Fig_meanslopes <- ggplot(data_cv_mslope, aes(x =Experiment, y = mslope)) +
  
  geom_bar( stat = "identity",width = 0.6, aes(color = Experiment,fill = Experiment)) +
  geom_errorbar(aes(ymin = mslope-seslope,ymax = 0,color = Experiment),width = 0.4)+
  
  scale_color_manual(values = myPairs)+
  scale_fill_manual(values = myPairs)+
  
  geom_signif(comparisons = list(c("Vis/Mix","Vis/Block")),
              annotations = "*",
              y_position = 0.01) +  
  geom_signif(comparisons = list(c("Vis/Mix","Aud/Mix")),
              annotations = "*",
              y_position = 0.02) +
  theme_minimal()+
  theme(legend.position = "null",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        axis.ticks = element_line(color = "grey"))+
  xlab('Conditions') +
  ylim(c(-0.05,0.025))+
  ylab('CV Slopes ')

Fig_meanslopes_mid <- plot_grid(NULL,Fig_meanslopes, NULL, nrow =3,rel_widths = c(1,10,1))

# ---- (New) Figure 3 ----
 fig3_new <- plot_grid(Figure3_cv, Fig_meanslopes_mid, nrow =1,rel_widths = c(2,0.7), labels = c('A','B'))

