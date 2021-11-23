# Estimate reference condition and emergence time
# Runs for each of the sensitivty analyses 
#CJ Brown 
# 2021-11-10
#

rm(list = ls())
library(ggplot2)
library(dplyr)
library(forecast)
library(DataGLMRepeat)
library(cowplot)
library(rethinking)
library(tidyr)
library(patchwork)

source("model-functions.R")

models <- list.files("../models/", pattern = "2021-08-03_model")
nmodels <- length(models)

td_dat <- NULL 

for (i in 1:nmodels){
  load(paste0("../models/",models[i]))
  post <- extract.samples(fitm1) %>% data.frame()
  post$imcmc <- 1:nrow(post)
  nmcmc <- nrow(post)
  
  ntrends <- 10 #trends for simulations
  trends <- seq(-0.3, 0, length.out = ntrends)
  
  nyears <- 30 #max years for simulations. 
  nsims <- nmcmc #simulations for flow baselines, up to nmcmc
  
  qbrel <- quantile(post$Brel.26)
  init_fraction1 <- 0.5 #just BMSY
  #qbrel[3] #median of most recent year
  
  effort_rel <- 1 #relative to EMSY
  
  # --------------------------
  # Historical flow variation and model sims
  # -------------------------- 
  
  #
  # Fit ARMA 
  #
  
  acf(datstan$flow)
  m1 <- Arima(dat2$Streamflow_std, c(6, 0, 3))
  
  #
  # Equilibrium effort 
  #
  EMSY <- quantile(post$EMSY, 0.5)
  eseries <- rep(EMSY * effort_rel, nyears)
  
  #
  # Check at equilibrium 
  #
  
  flowseed <- 30
  flowsim <- simulate(m1, nsim = nyears,
                      seed = flowseed, future = FALSE) %>%
    as.numeric()
  flowsim <- (flowsim - mean(flowsim))/sd(flowsim)
  
  flowsim <- flowsim + (0:(nyears - 1)*(-0.1))

  #
  # Simulate historical flow series and from model 
  #
  
  flowseeds <- 1:nsims
  datsim <- expand.grid(flowseed = flowseeds,
                        trend = trends)
  
  
  #note below code always uses the same mcmc sample for the 
  # same flow seed
  #residuals are random
  datout <- datsim %>%
    group_by(flowseed, trend) %>%
    DataGLMRepeat::with_groups({
      flowsim <- simulate(m1, nsim = nyears,
                          seed = flowseed, future = FALSE) %>%
        as.numeric()
      flowsim <- flowsim + (0:(nyears - 1)*trend)
      
      x <- simmod(flowseed, post, nyears,
                  init_fraction1, eseries,
                  flowsim, ((nmcmc*2):(nmcmc*3))[flowseed])
      
      x$trend <- trend
      x
    })
  
  dout2 <- do.call("rbind", datout)
  
  #unstandardize
  
  dout2$ndvi_mean <- (dout2$ndvi_mean*sd(dat2$NDVI, na.rm = TRUE)) + mean(dat2$NDVI, na.rm = TRUE)
  dout2$pasture_mean <- (dout2$pasture_mean*sd(dat2$pasture_bio, na.rm = TRUE)) + mean(dat2$pasture_bio, na.rm = TRUE)
  
  #
  # Plot indicator reference points 
  #
  
  #summarize quantiles for each year
  dsum <- dout2 %>% group_by(n, trend) %>%
    summarize(across(lnCPUE:pasture_mean, ~quantile(.x, 0.05, na.rm = TRUE),
                     .names = "{.col}-q05"),
              across(lnCPUE:pasture_mean, ~quantile(.x, 0.2, na.rm = TRUE),
                     .names = "{.col}-q20"),
              across(lnCPUE:pasture_mean, ~quantile(.x, 0.4, na.rm = TRUE),
                     .names = "{.col}-q40"),
              across(lnCPUE:pasture_mean, ~quantile(.x, 0.5, na.rm = TRUE),
                     .names = "{.col}-q50"),
              across(lnCPUE:pasture_mean, ~quantile(.x, 0.6, na.rm = TRUE),
                     .names = "{.col}-q60"),
              across(lnCPUE:pasture_mean, ~quantile(.x, 0.8, na.rm = TRUE),
                     .names = "{.col}-q80"),
              across(lnCPUE:pasture_mean, ~quantile(.x, 0.95, na.rm = TRUE),
                     .names = "{.col}-q95")) %>%
    pivot_longer(3:51, names_to = "Variable", 
                 values_to = "Val") %>%
    separate(Variable, into = c("Var", "Quant"), sep = "-") %>%
    pivot_wider(names_from = "Quant",
                values_from = "Val")
  
  dloop <- expand.grid(Var = c("lnCPUE_mean", "ndvi_mean", "pasture_mean"),
                       trends = trends[trends<0], 
                       qs = c("q05", "q20", "q40"))
  
  dloop$td <- lapply(1:nrow(dloop), function(x) 
    get_td(dsum, dloop$qs[x],dloop$trends[x], 
           dloop$Var[x])) %>%
    as.numeric()
  dloop$model <- models[i]
  td_dat <- c(td_dat, list(dloop))
              
}

td_dat <- bind_rows(td_dat)
td_dat$Quantile <- factor(td_dat$qs, labels = c("5%", "20%", "40%"))
td_dat$Var <- factor(td_dat$Var, labels = c("CPUE", "NDVI", "Pasture"))
td_dat$`Initial biomass (%)` <-  with(td_dat, 
                      as.numeric(substr(model,start = (nchar(models[1])-12),
                                        stop = (nchar(models[1])-11))))
td_dat$qtrend <- paste0(with(td_dat, 
                      as.numeric(substr(model,start = (nchar(models[1])-4),
                                        stop = (nchar(models[1])-4)))), "% p.a.")

# ------------ 
# Plots 
# ------------ 
td_dat_cpue <- filter(td_dat, Var == "CPUE")
g1 <- ggplot(td_dat_cpue) + 
  aes(x = trends, y = td, color = `Initial biomass (%)` ,
      group = model) + 
  geom_line() + 
  facet_grid(qtrend~Quantile) + 
  xlab("Trend in flow (multiples /yr)") + 
  ylab("Emergence time (years)") +
  theme_classic()

ggsave("../figures/2021-08-06_emergence-times-sensitivity.png", 
       g1)
