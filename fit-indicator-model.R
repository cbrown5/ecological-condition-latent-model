# Fit the multivariate indicator model to data
# including state-space model for barramundi CPUE
#
#plots some exploratory plots of results
#
#CJ Brown 2021-11-10
#
# Note this script was rerun for different settings of: 
# bf_guess (initial biomass) (base setting = 0.2)
# catchability_increase (base setting = 1)
# See Figures S2-4 for results. 

rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(rstan)
library(rethinking)

dat <- read.csv("Data for Barramundi model.csv")
dat$catch <- dat$Barramundicatch
dat$days <- dat$Effort
n <- nrow(dat)

#calculates params for lognormal based on mean and SD
logNormalParams <- function(mB, sigma_B){
  a <-2*log(mB) - 0.5 * log(sigma_B^2 + mB^2) 
  b <- sqrt(-2*log(mB)+log(sigma_B^2 + mB^2))
  return(list(lmean = a, lsigma = b))
}

s2 <- function(x){
  (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

#
# Prep data 
#
catchability_increase <- 1 #1.01 = 1% per year

#Calculate what a 1SD change in log(flow) means
sd(log(dat$Streamflow), na.rm = TRUE)
exp(sd(log(dat$Streamflow), na.rm = TRUE))

dat2 <- dat %>%
  mutate(Streamflow = log(Streamflow_wetseason)) %>%
  filter(!is.na(Streamflow)) %>%
  mutate(#days_std = days/365,
          days_std = (days/365)*catchability_increase^(Year-1990),
         ln_cpue = log(catch/days_std),
         Streamflow_std = (Streamflow - mean(Streamflow, na.rm = TRUE))/sd(Streamflow, na.rm = TRUE),
         ndvi_std = s2(NDVI),
         inun_std = s2(inundation),
         pasture_std = s2(pasture_bio)
         ) 

cor(dat2, use = "pairwise.complete.obs")

#Data summaries 
datsum <- dat2 %>% 
  summarize(across(inundation:ln_cpue, list(mn = mean, sd = sd), na.rm = TRUE))

# ----------------- #
# Prior params 
# ----------------- #

#
# B1 and B0 
#

#Guess B1 and B0
r_est <- 0.3

#sensitivty analysis
bf_guess <- 0.2 #fraction biomass is of Bzero in 1989
b1_fract <- 10*(bf_guess/0.2) # Multiple B1 is of catch in 1989, relative to 0.2, which
# was used in the main analysis

B0_CV <- 1 # CV for lognormal on B0 - prior 

B1_guess <- dat2$catch[1]*b1_fract
B0_guess <- B1_guess*(1/bf_guess) #assume B1 is XX% of B0
B0params <- logNormalParams(B0_guess, B0_guess * B0_CV)
plot(density(exp(rnorm(10000, B0params$lmean, B0params$lsigma))))
quantile(exp(rnorm(10000, B0params$lmean, B0params$lsigma)), c(0.1, 0.9))

#Guess q 
# C = BEq
#CPUE = Bq
exp(dat2$ln_cpue)/B1_guess
median(exp(dat2$ln_cpue)/B1_guess)

rparams <- logNormalParams(r_est, r_est*0.3)
# plot(density(exp(rnorm(10000, rparams$lmean, rparams$lsigma))))

#
# Obs and process errors
#

#Choose parameters so you get desired multiplier on CPUE/biomass

#Prior quantiles on SD
qexp(c(0.05,0.95), 2)
qexp(c(0.05,0.95), 14)


#Process 
qlwr <- 0.05
quant <- 0.9
rgaus <- rnorm(50000, mean = 0, sd = rexp(50000, 14))
quantile(exp(rgaus), c(qlwr, quant+qlwr)) 
# ~ 1.15

#Observation
rgaus <- rnorm(50000, mean = 0, sd = rexp(50000, 2))
quantile(exp(rgaus), c(qlwr, quant+qlwr)) 
# ~ 1.5

# ----------------- #
# Model 
# ----------------- #

#
#Setup data 
#

datstan <- with(dat2, {
  x = list(N = nrow(dat2),
       flow = Streamflow_std, 
       #stock model 
       lnCPUE = ln_cpue,
       catches = catch,
       logK_mean = B0params$lmean,
       logK_sd = B0params$lsigma,
       init_fraction = bf_guess,
       logr_mean = rparams$lmean,
       logr_sd = rparams$lsigma,

       # N for indicators
       Nndvi = sum(!is.na(ndvi_std)),
       Npasture = sum(!is.na(pasture_std)),
       
       #Indices for years
       i_ndvi = which(!is.na(ndvi_std)),
       i_pasture = which(!is.na(pasture_std))
  )
  #Indicator data 
  x$ndvi = ndvi_std[x$i_ndvi]
  x$pasture = pasture_std[x$i_pasture]
  x
})

#
# Run model 
#
options(mc.cores = 3) #parallel::detectCores())
              
fitm1 <- stan(file = "indicator-model.stan", data = datstan,
               iter=5000, chains=3, thin = 5,
               init = list(
                 list(lnr = log(r_est) , lnq = log(0.05), sigma_cpue = 0.1,
                      sigma_u =0.05, lnK = log(3000), beta_u = 0.05,
                      beta_nu = 1
                 ),
                 list(lnr = log(r_est)*1.2 , lnq = log(0.1), sigma_cpue = 0.15,
                      sigma_u = 0.1, lnK = log(5000),
                      beta_u = 0.1,
                      beta_nu = 2
                 ),
                 list(lnr = log(r_est)*0.6 , lnq = log(0.02), sigma_cpue = 0.05,
                      sigma_u = 0.2, lnK = log(2000),
                      beta_u = 0.15,
                      beta_nu = 0.1
                 )
               ),
               control = list(max_treedepth = 12, adapt_delta = 0.8)
)

save(dat2, datsum, datstan, 
     fitm1, file = paste0("../models/2021-08-03_model-fit-no-inun-logflow-Binit",bf_guess*100,
                          "perc-q",(catchability_increase-1)*100,".rda"))


# ----------------- #
# Model checks 
# ----------------- #

#choose model to follow up on
# shinystan::launch_shinystan(fitm1)

{s <- summary(fitm1)
  rn <- row.names(s$summary)
  i <- !grepl("\\[", rn)
  signif(s$summary[i,c(1, 4, 6, 8, 9, 10)], 3)}

post <- extract.samples(fitm1) %>% data.frame()

par(mfrow = c(1,2))
plot(density(post$sigma_cpue))
plot(density(rexp(10000, 2)))
quantile(post$sigma_cpue, c(0.05, 0.95))

plot(density(post$sigma_u))
plot(density(rexp(10000, 14)))
quantile(post$sigma_u, c(0.05, 0.5, 0.95))

plot(post$beta_u)
plot(density(post$beta_u))

plot(density(exp(rnorm(10000, rparams$lmean, rparams$lsigma))))
lines(density(post$r))

plot(density(exp(rnorm(10000, B0params$lmean, B0params$lsigma))))
lines(density(post$K))

plot(density(post$q))

plot(post$q, post$K)

plot(post$sigma_u, post$sigma_cpue)
plot(post$sigma_cpue, post$lp__)

par(mfrow = c(2,2))
plot(density(post$beta_ndvi))
plot(density(post$beta_pasture))
plot(density(post$beta_nu))

par(mfrow = c(1,1))
# ----------------- #
# Model predictions
# ----------------- #
{txtplot2 <- theme(axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12))



getsamps <- function(post, grep_string, datin){
  ib <- grep(grep_string, names(post))
  post_expected_bio <- as.matrix(post[,ib])
  
  dout <- t(apply((post_expected_bio), 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75,0.975))) %>%
    data.frame() %>% bind_cols(datin)
  
  dout
}
}

{bio_est <- getsamps(post, "\\bBrel[.]", dat2)

gbio <- ggplot(bio_est) + 
  geom_hline(yintercept = 0.5, color = "grey20")+ 
  geom_hline(yintercept = 1, color = "grey20")+ 
  aes(x = Year, color = NULL) +
  # geom_line(aes(y = 0.4)) + 
  geom_ribbon(aes(ymin = X25., ymax = X75.), alpha = 0.35, fill = "blue") + 
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5, fill = "blue") + 
  geom_line(aes(y = X50.)) + 
  theme_classic()  + 
  xlab("Year") + 
  ylab("Relative biomass") + 
  txtplot2

cpue_est <- getsamps(post, "lnCPUE_hat[.]", dat2)

gcpue <- ggplot(cpue_est) + 
  aes(x = Year, y = exp(ln_cpue), color = NULL) +
  geom_point() + 
  geom_line(aes(y = exp(X50.))) + 
  geom_ribbon(aes(ymin = exp(X2.5.), ymax = exp(X97.5.)), alpha = 0.5, color = NA,
              fill = "red") + 
  theme_classic() + 
  ylab("CPUE") + 
  xlab("") + 
  txtplot2
}

cowplot::plot_grid(gbio, gcpue)

acf(cpue_est$X50. - cpue_est$ln_cpue)  

#
# Plot indicators 
#


ind_est <- list(getsamps(post, "ndvi_pred[.]", dat2),
 getsamps(post, "pasture_pred[.]", dat2))
thisind <- c("ndvi_std", "pasture_std")

par(mfrow = c(1,2))
ginds <- NULL
for (i in 1:length(ind_est)){
  gtemp <- ggplot(ind_est[[i]]) + 
    aes_string(x = "Year", y = thisind[i], color = NULL) +
    geom_point() + 
    geom_line(aes(y = (X50.))) + 
    geom_ribbon(aes(ymin = (X2.5.), ymax = (X97.5.)), alpha = 0.5, color = NA,
                fill = "red") + 
    theme_classic() + 
    ylab(thisind[i]) + 
    xlab("") + 
    txtplot2

  ginds <- c(ginds, list(gtemp))
  acf(cpue_est$X50.[!is.na(cpue_est[,thisind[i]])] - cpue_est[!is.na(cpue_est[,thisind[i]]),thisind[i]],
      main = thisind[i])  
}
ginds <- c(list(gcpue), ginds)
cowplot::plot_grid(plotlist = ginds)

  
  