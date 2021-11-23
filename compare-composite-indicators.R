# Compare composite indicators: 
# model based versus weighted sum

#CJ Brown 
# 2021-11-10
#

rm(list = ls())
library(ggplot2)
library(dplyr)
library(DataGLMRepeat)
library(cowplot)
library(rethinking)
library(tidyr)
library(patchwork)

source("model-functions.R")

load("../models/2020-12-11_model-fit-no-inun-logflow.rda")

#
# Weighted sum aggregate indicator 
#

convert_01 <- function(x){
  (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

z <- cbind(convert_01(dat2$ln_cpue),
 convert_01(dat2$NDVI),
 convert_01(dat2$pasture_bio))

Z <- rowMeans(z, na.rm = TRUE)


# 
# Model based ecological condition 
# 
post <- extract.samples(fitm1) %>% data.frame()
getsamps <- function(post, grep_string, datin){
  ib <- grep(grep_string, names(post))
  post_expected_bio <- as.matrix(post[,ib])
  
  dout <- t(apply((post_expected_bio), 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75,0.975))) %>%
    data.frame() %>% bind_cols(datin)
  
  dout
}

nu_est <- getsamps(post, "nu[.]", dat2)
nu_est$comp_ind <- Z
nu_mean <- getsamps(post, "nu_hat[.]", dat2)

gnu <- ggplot(nu_est) + 
  aes(x = Year, color = NULL) +
  geom_hline(yintercept = 0) + 
  geom_line(aes(y = X50.)) + 
  geom_ribbon(aes(ymin = X25., ymax = X75.), alpha = 0.35, color = NA, fill = "seagreen") + 
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5, color = NA, fill = "seagreen") + 
  theme_classic() + 
  ylab("Ecological condition \n indicator") 
gnu


# ------------ 
# Compare condition indicators 
# ------------ 
years_select <- c(1991, 2013, 2011)
nu_est <- mutate(nu_est,year2 = 
                     ifelse(Year %in% years_select,
                            Year,
                            NA))
nu_est$years_missing <- is.na(nu_est$NDVI)


g_covar <- ggplot(nu_est) + 
  aes(x = comp_ind,y = X50., color = NULL) +
  geom_point(aes(col = years_missing)) +
  geom_linerange(aes(ymin = X25., ymax = X75.), alpha = 0.35) + 
  # geom_linerange(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5) + 
  geom_label(aes(label = year2),
             nudge_x = 0.07, nudge_y = 0.25, 
             check_overlap = TRUE,
             size = 2) +
  theme_classic() + 
  stat_smooth() + 
  ylab("Model-based \n indicator") +
  xlab("Additive indicator") +
  scale_colour_manual(values = c("black", 'red')) +
  theme(legend.position = "none")
g_covar

gnu <- ggplot(nu_est) + 
  aes(x = Year,y = X50., color = NULL) +
  geom_point(aes(col = years_missing))+
  geom_line()+
  geom_linerange(aes(ymin = X25., ymax = X75.), alpha = 0.35) + 
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5) + 
  theme_classic() + 
  ylab("Model-based \n indicator") +
  xlab("Year")+
  scale_colour_manual(values = c("black", 'red')) +
  theme(legend.position = "none")
gnu

gcomp <- ggplot(nu_est) + 
  aes(x = Year,y = comp_ind, color = NULL) +
  geom_point(aes(col = years_missing))+
  geom_line() + 
  theme_classic() + 
  xlab("Year") +
  ylab("Additive indicator")+
  scale_colour_manual(values = c("black", 'red')) +
  theme(legend.position = "none")
gcomp


cor(nu_est$X50., nu_est$comp_ind)

pw <- g_covar / gnu / gcomp

pw2 <- pw  + plot_annotation(tag_levels = 'A') + 
  plot_layout() & 
  theme(plot.tag = element_text(size = 16)) 

ggsave(pw2, file = "../figures/2020-06-15_composite-indicator-plot.png",
       width = 3, height = 6)


