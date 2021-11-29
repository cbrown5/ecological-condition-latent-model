# Fit model to simulated data 
#CJ Brown 
# 2021-11-10
#
# Simulate dataseries with different amounts of missing
# data 

library(ggplot2)
library(dplyr)
library(forecast)
library(DataGLMRepeat) #not on cran, available from github/cbrown5
library(cowplot)
library(rethinking)
library(tidyr)
library(rstan)

source("model-functions.R")

load("../models/2021-08-03_model-fit-no-inun-logflow-Binit20perc-q0.rda")
post <- extract.samples(fitm1) %>% data.frame()
postmed <- data.frame(t(apply(post, 2, median)))
post_mean <- data.frame(t(apply(post, 2, mean)))

N1 <- 26
init_fraction1 <- 0.2
nsims <- 30
r_est <- 0.3

params <- c("beta_u",
            "beta_nu",
            "beta_ndvi",
            "beta_pasture",
            "a_ndvi",
            "a_pasture",
            "lnr",
            "lnK",
            "q",
            "sigma_ndvi",
            "sigma_pasture",
            "sigma_cpue" #CPUE SD
)

# ----------------- #
# Missing years of data
# ----------------- #

yrs_missing <- c(1, 6, 11, 16) #first year, so 11 = 10 years missing, 1 =
# yrs_missing <- c(1, 6) #first year, so 11 = 10 years missing, 1 =
# no years missing 

# --------------------------
# Historical flow variation and model sims
# -------------------------- 

#
# Fit ARMA 
#
nyears <- nrow(dat2)
acf(datstan$flow)
m1 <- Arima(dat2$Streamflow_std, c(6, 0, 3))
AIC(m1)
summary(m1)
acf(resid(m1))
plot(simulate(m1, nyears, seed = 42, future = FALSE))
points(simulate(m1, nyears, seed = 42, future = FALSE))
lines(dat2$Streamflow_std, col = "red")

#
# Effort series 
#
eseries <- dat2$Effort/365
 

# ------------ 
# Simulate historical flow series and from model 
# ------------ 

flowseeds <- 1:nsims
datsim <- expand.grid(flowseed = flowseeds,
                      yrs_missing = yrs_missing)

compiled_model <- stan_model("indicator-model.stan")

datout <- datsim %>%
  group_by(flowseed, yrs_missing) %>%
  DataGLMRepeat::with_groups({
    #
    # Simulate data 
    #
    
    flowsim <- simulate(m1, nsim = nyears,
                        seed = flowseed, future = FALSE) %>%
      as.numeric()
  
    x <- simmod(1, postmed, nyears, 
                init_fraction1, eseries,
                flowsim, 1)
    
    #Setup data and centre NDVI and pasture variables
    datstan2 <- within(datstan, {
      lnCPUE = x$lnCPUE
      catches = x$catches
      pasture = x$pasture - mean(x$pasture)
      Npasture = length(x$pasture)
      i_pasture = 1:length(x$pasture)
      flow = flowsim
      
      #missing data for NDVI
      
      Nndvi = nyears - yrs_missing + 1
      i_ndvi = yrs_missing:nyears
      ndvi = x$ndvi[yrs_missing:nyears] - mean(x$ndvi, na.rm = TRUE)
    })
    
    
    options(mc.cores = 3) #parallel::detectCores())

    fitm2 <- sampling(compiled_model, data = datstan2,
                  iter=3000, chains=3, thin = 5, 
                  init = list(
                    list(lnr = log(r_est) , lnq = log(0.05), sigma_cpue = 0.1,
                         sigma_u =0.01, lnK = log(3000), 
                         beta_u = 0.01,
                         beta_nu = 1.5,
                         beta_ndvi = 0.3,
                         beta_pasture = 0.05
                    ),
                    list(lnr = log(r_est)*1.2 , lnq = log(0.1), sigma_cpue = 0.15,
                         sigma_u = 0.1, lnK = log(5000),
                         beta_u = 0.05,
                         beta_nu = 2,
                         beta_ndvi = 0.2,
                         beta_pasture = 0.15
                    ),
                    list(lnr = log(r_est)*0.6 , lnq = log(0.02), sigma_cpue = 0.05,
                         sigma_u = 0.2, lnK = log(2000),
                         beta_u = 0.1,
                         beta_nu = 1,
                         beta_ndvi = 0.45,
                         beta_pasture = 0.3
                    )))
    
    
    
    #
    # Summary
    #
    s <- summary(fitm2)
    rn <- row.names(s$summary)
    i <- !grepl("\\[", rn)
    sout <- data.frame(signif(s$summary[i,c(1, 4, 6, 8, 9, 10)], 3))
    sout$isim <- flowseed
    sout$yrs_missing <- yrs_missing
    
    rhats <- data.frame(param = row.names(sout), Rhat = sout$Rhat)
    
    postsim <- extract.samples(fitm2) %>% data.frame()
    ps2 <- postsim[,names(postsim) %in% params] %>%
      pivot_longer(1:length(params), names_to = "param",
                   values_to = "val") %>%
      group_by(param) %>%
      summarize(med = median(val),
                mn = mean(val), 
                upr = quantile(val, 0.975),
                lwr = quantile(val, 0.025),
                upr50 = quantile(val, 0.75),
                lwr50 = quantile(val, 0.25))
    
    p2 <- postmed %>% pivot_longer(1:ncol(postmed), names_to = "param",
                                   values_to = "truth")
    psall <- ps2 %>% left_join(p2) %>% left_join(rhats)
    psall$isim <- flowseed
    psall$yrs_missing <- yrs_missing
  
    #
    # R2 values 
    #
    fitcors <- data.frame(ndvi = NA,
                          pasture = NA,
                          lncpue = NA,
                          nu = NA)
    i <- grepl("ndvi_pred", rn)
    fitcors$ndvi <- cor(x$ndvi,data.frame(s$summary[i,])$X50.)^2
    i <- grepl("pasture_pred", rn)
    fitcors$pasture <- cor(x$pasture,data.frame(s$summary[i,])$X50.)^2
    i <- grepl("lnCPUE_hat", rn)
    fitcors$lncpue = cor(datstan2$lnCPUE,data.frame(s$summary[i,])$X50.)^2
    i <- grepl("nu\\[", rn)
    fitcors$nu <- cor(x$nu,data.frame(s$summary[i,])$X50.)^2
    fitcors$sim <- flowseed
    fitcors$yrs_missing <- yrs_missing
    fitout <- list(din = datstan2, #sfit = sout,
                             paramsum = psall,
                   fitcors = fitcors)
    fitout$yrs_missing <- yrs_missing
    fitout
  })

save(datout,params, file = "../models/2021-11-24_param-sims.rda")


