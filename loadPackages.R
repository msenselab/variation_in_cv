# Install necessary packages
usePackage <- function(pk){
  for (p in pk) {
    if (!is.element(p, installed.packages()[,1])) install.packages(p, dependencies =TRUE)
    library(p, character.only = TRUE)
  }
}
usePackage(c('knitr','R.matlab','data.table','xtable','dplyr','tidyr',
             'grid','gridExtra','broom','ggplot2','nlme','ez','RColorBrewer',
             'scales','ggsignif','cowplot','lemon','tidyverse','ggsignif','apaTables'))
library(RColorBrewer)
library(tidyverse)
library(ggsignif)
library(lemon)


options(digits=3)

theme_set(theme_bw(base_size = 16))

# Functions: control legends etc. theme
legend_pos <-function(l) {
  theme( plot.title = element_text(face="bold",size = 12),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = l,
         legend.title = element_blank(),
         legend.text = element_text(size=12),
         legend.key = element_blank() # switch off the rectangle around symbols
         #axis.line.x = element_line(color="black", size = 0.2),
         #axis.line.y = element_line(color="black", size = 0.2)
  )
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# function for slope estimation from each subject
slope_estimation <- function(df,lm_formula) {
  df %>% ungroup %>%
    group_by(exp, sub) %>% 
    do(tidy(lm(as.formula(lm_formula), data=., na.action=na.exclude))) %>% 
    filter(., term != '(Intercept)') %>%
    #    select(., estimate) %>%
    rename(., slope = estimate) %>%
    separate(.,exp,c('modality','design'),sep='/', remove = FALSE)
}
