
# ----------------- #
# Plotting model fits 
# ----------------- #
# CJ Brown 2021-11-10
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(rstan)
library(patchwork)

#version with flow logged
load("../models/2021-08-03_model-fit-no-inun-logflow-Binit20perc-q0.rda")


 post <- rstan::extract(fitm1) %>% data.frame()

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
txtplot2 <- theme(axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12))

#
#CIs
#

qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
quantile(post$beta_nu, probs = qs)
quantile(post$beta_ndvi, probs = qs)
quantile(post$beta_pasture, probs = qs)
quantile(post$beta_u, probs = qs)

sum(post$beta_nu>0)/nrow(post)
sum(post$beta_ndvi>0)/nrow(post)
sum(post$beta_pasture>0)/nrow(post)
sum(post$beta_u>0)/nrow(post)

sum((post$beta_ndvi*post$beta_nu)>0)/nrow(post)
sum((post$beta_pasture * post$beta_nu)>0)/nrow(post)
sum((post$beta_u * post$beta_nu)>0)/nrow(post)


# ----------------- #
# Biomass 
# ----------------- #

bio_est <- getsamps(post, "\\bBrel[.]", dat2)

gbio <- ggplot(bio_est) + 
  geom_hline(yintercept = 0.5, color = "grey20")+ 
  geom_hline(yintercept = 1, color = "grey20")+ 
  aes(x = Year, color = NULL) +
  # geom_line(aes(y = 0.4)) + 
  geom_ribbon(aes(ymin = X25., ymax = X75.), alpha = 0.35, fill = "blue") + 
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5, fill = "blue") + 
  geom_line(aes(y = X50.)) + 
  theme_classic()  + 
  ylim(0, 0.7) + 
  xlab("Year") + 
  ylab("Relative biomass") + 
  txtplot
gbio

# ----------------- #
# Ecological condition 
# ----------------- #

nu_est <- getsamps(post, "nu[.]", dat2)
nu_mean <- getsamps(post, "nu_hat[.]", dat2)

gnu <- ggplot(nu_est) + 
  aes(x = Year, color = NULL) +
  geom_hline(yintercept = 0) + 
  geom_line(aes(y = X50.)) + 
  geom_ribbon(aes(ymin = X25., ymax = X75.), alpha = 0.35, color = NA, fill = "seagreen") + 
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5, color = NA, fill = "seagreen") + 
  theme_classic() + 
  ylab("Ecological condition \n index") + 
  txtplot
gnu

gnuflow <- ggplot() + 
  geom_hline(yintercept = 0) + 
  geom_line(data = nu_mean, aes(x = Streamflow, y = X50.), size = 1.5, color = "tomato") + 
  # geom_ribbon(data =nu_mean, aes(x = Streamflow, ymin = X2.5., ymax = X97.5.), color= NA, fill = "red",
  # alpha = 0.5) + 
  
  geom_linerange(data = nu_est, aes(x = Streamflow, ymin = X2.5., ymax = X97.5.), 
                 color = "grey40") + 
  geom_linerange(data = nu_est, aes(x = Streamflow, ymin = X25., ymax = X75.), 
                 size = 1.5,
                 color = "grey40") + 
  geom_point(data = nu_est, aes(x = Streamflow, y = X50.)) + 
  theme_classic() + 
  ylab("Ecological condition \n indicator") + 
  txtplot
gnuflow

# ----------------- #
# Data plots 
# ----------------- #

#Bayesian correlations
r_cpue <- brsq(post,"lnCPUE_hat[.]", "sigma_cpue")
r_ndvi <- brsq(post,"ndvi_hat[.]", "sigma_ndvi")
r_pasture <- brsq(post,"pasture_hat[.]", "sigma_pasture")

rsqout <- data.frame(Variable = c("lnCPUE", "NDVI", "Pasture biomass"), 
                     rbind(r_cpue, r_ndvi, r_pasture))
write.csv(rsqout, "../figures/2021-08-03_model_rsq-logflow.csv", row.names = FALSE)

#
# CPUE 
#

cpue_est <- getsamps(post, "lnCPUE_hat[.]", dat2)

gcpue <- ggplot(cpue_est) + 
  aes(x = Year, y = exp(ln_cpue), color = NULL) +
  geom_point() + 
  geom_line(aes(y = exp(X50.))) + 
  geom_ribbon(aes(ymin = exp(X2.5.), ymax = exp(X97.5.)), alpha = 0.5, color = NA,
              fill = "red") + 
  theme_classic() + 
  ylab("CPUE \n (tonnes/yr)") + 
  xlab("") + 
  txtplot2
gcpue

#
# NDVI
#

ndvi_est <- getsamps(post, "ndvi_pred[.]", dat2)

gndvi <- ggplot(ndvi_est) + 
  aes(x = Year, y = ndvi_std, color = NULL) +
  geom_point() + 
  geom_line(aes(y = X50.)) + 
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5, color = NA,
              fill = "red") + 
  theme_classic() + 
  ylab("NDVI \n (standardized)") + 
  xlab("") + 
  txtplot2
gndvi



#
# pasture
#

pasture_est <- getsamps(post, "pasture_pred[.]", dat2)

gpasture <- ggplot(pasture_est) + 
  aes(x = Year, y = pasture_std, color = NULL) +
  geom_point() + 
  geom_line(aes(y = X50.)) + 
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5, color = NA,
              fill = "red") + 
  theme_classic() + 
  ylab("Pasture biomass \n (standardized)") + 
  xlab("") + 
  txtplot2
gpasture

#
#Process errors
#
#not of interest? 

u_est <- getsamps(post, "^u[.]", dat2)

gmu <- ggplot(u_est) + 
  geom_line(aes(x = Year, y = X50.)) + 
  geom_ribbon(aes(x = Year, ymin = X2.5., ymax = X97.5.), alpha = 0.5, color = NA,
              fill = "red") + 
  theme_classic() + 
  ylab("Barramundi process errors") + 
  txtplot
gmu

# ----------------- #
# Flow 
# ----------------- #

gflow <- ggplot(dat2) + 
  aes(x = Year, color = NULL) +
  geom_line(aes(y = exp(Streamflow)), size = 2, color = "darkblue") + 
  theme_classic() + 
  ylab(expression(Flow~(m^{3}/s))) + 
  xlab("") + 
  txtplot2
gflow

# ----------------- #
# Plot compilation
# ----------------- #

pw <- (gcpue / gndvi /gpasture/gflow) | (gnu  / gbio)

pw2 <- pw  + plot_annotation(tag_levels = 'A') + 
  plot_layout(ncol = 2, widths = c(1,2)) & 
  theme(plot.tag = element_text(size = 16)) 

ggsave(pw2, file = "../figures/2022-06-18_data-fits-logflow.png",
       width = 10, height = 8)

# gflow
# gcpue
# gndvi
# ginun
# 
# gbio
# gnu
# gnuflow
# gmu
# 
# ----------------- #
# Checking residuals for AC
# ----------------- #

lmax <- 8 #lags to check

 png(filename = "../figures/2021-08-03_ACF-logflow.png")
par(mfrow = c(2,2))
acf(nu_est$X50. - nu_mean$X50., main = "Ecosystem condition",
    lag.max = lmax)
acf(cpue_est$X50. - dat2$ln_cpue, main = "log(CPUE)", 
    lag.max = lmax)
acf(ndvi_est$X50.[datstan$i_ndvi] - dat2$ndvi_std[datstan$i_ndvi], main = "NDVI", 
    lag.max = lmax)
acf(pasture_est$X50.[datstan$i_pasture] - dat2$pasture_std[datstan$i_pasture], main = "Pasture biomass", lag.max = lmax)
 dev.off()  
