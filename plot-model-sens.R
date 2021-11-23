
# ----------------- #
# Plotting model fits from sensitivity analyses 
# ----------------- #
# CJ Brown 2021-11-10
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(rstan)
library(rethinking)
library(patchwork)

files <- c("2021-07-07_model-fit-no-inun-logflow-Binit20perc-q1.rda",
           "2021-07-07_model-fit-no-inun-logflow-Binit50perc-q1.rda",
           "2021-07-07_model-fit-no-inun-logflow-Binit80perc-q1.rda",
           "2020-12-11_model-fit-no-inun-logflow.rda",
           "2021-07-07_model-fit-no-inun-logflow-Binit50perc.rda",
           "2021-07-07_model-fit-no-inun-logflow-Binit80perc.rda"
           )

#Simulations matching above files
mods <- data.frame(Binit = rep(c(0.2, 0.5, 0.8), 2), 
                   qchange = rep(c(1, 0), each = 3))


#
#Function to extract samples for temporal params
#

getsamps <- function(post, grep_string, datin){
  ib <- grep(grep_string, names(post))
  post_expected_bio <- as.matrix(post[,ib])
  
  dout <- t(apply((post_expected_bio), 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75,0.975))) %>%
    data.frame() %>% bind_cols(datin)
  
  dout
}

brsq <- function(post,grep_string, grep_string_sd){
  #from: https://statmodeling.stat.columbia.edu/2017/12/21/r-squared-bayesian-regression-models/
  ib <- grep(grep_string, names(post))
  pred <- as.matrix(post[,ib])
  isigma <- grep(grep_string_sd, names(post))
  i <- 1:nrow(post)
  sigma <- as.matrix(post[,isigma])
  cvals <- lapply(i, function(x) signif(var(pred[x,])/(var(pred[x,]) + (sigma[x]^2)),2))
  quantile(unlist(cvals), probs = c(0.025, 0.25, 0.5, 0.75,0.975))
}

# brsq(post,"lnCPUE_hat[.]", "sigma_cpue")

txtplot <- theme(axis.text = element_text(size = 14),
                 axis.title = element_text(size = 16))
txtplot2 <- theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14),
                  axis.text.x = element_text(size = 10))

# ------------ 
# Extract results from each  model run 
# ------------ 
rout <- NULL
pout <- NULL 

for (i in 1:length(files)){
  load(paste0("../models/",files[i]))
  post <- extract.samples(fitm1) %>% data.frame()
  
  #Bayesian correlations
  r_cpue <- brsq(post,"lnCPUE_hat[.]", "sigma_cpue")
  r_ndvi <- brsq(post,"ndvi_hat[.]", "sigma_ndvi")
  r_pasture <- brsq(post,"pasture_hat[.]", "sigma_pasture")
  
  rsqout <- data.frame(Variable = c("lnCPUE", "NDVI", "Pasture biomass"), 
                       rbind(r_cpue, r_ndvi, r_pasture)) %>%
    cbind(mods[i,])
  
  rout <- rbind(rout, rsqout)
  
  qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  paramquants <- rbind(beta_nu = quantile(post$beta_nu, probs = qs),
                       beta_ndvi = quantile(post$beta_ndvi, probs = qs), 
                       beta_pasture = quantile(post$beta_pasture, probs = qs),
                       beta_u = quantile(post$beta_u, probs = qs)
  ) %>%
    data.frame() %>%
    mutate(param = row.names(.)) %>%
    cbind(mods[i,])
  
  pout <- rbind(pout, paramquants)

}

#
# Plots 
#
pd <- position_dodge(width=0.1)
grsq <- ggplot(rout) + 
  aes(x = Binit, y = X50., color = factor(qchange)) + 
  geom_point(position = pd) + 
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.),
                 position = pd) + 
  facet_wrap(~Variable) + 
  xlab("Initial biomass (prop. of K)") + 
  ylab(bquote('Bayesian '~R^2)) + 
  labs(color = "Catchability \n increase \n p.a. (%)") +
  txtplot2 +
  theme_classic()

gparams <- ggplot(pout) + 
  aes(x = Binit, y = X50., color = factor(qchange)) + 
  geom_point(position = pd) + 
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.),
                 position = pd) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  facet_wrap(~param, scales = "free") + 
  xlab("Initial biomass (prop. of K)") + 
  ylab("Value") + 
  labs(color = "Catchability \n increase \n p.a. (%)") +
  txtplot +
  theme_classic()
gparams

ggsave(gparams, file = "../figures/2021-08-04_parameter-sensitivity.png")
ggsave(grsq, 
       file = "../figures/2021-08-04_parameter-sensitivity-rsq.png",
       width = 5, height = 2)
