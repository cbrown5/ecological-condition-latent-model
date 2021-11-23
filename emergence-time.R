# Estimate reference condition and emergence time
#CJ Brown 
# 2021-02-15
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

load("../models/2021-08-03_model-fit-no-inun-logflow-Binit20perc-q0.rda")
#file not provided on github due to file size

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

effort_rel <- c(1, 1.5) #relative to EMSY

# --------------------------
# Historical flow variation and model sims
# -------------------------- 


#
# Fit ARMA 
#

acf(datstan$flow)
m1 <- Arima(dat2$Streamflow_std, c(6, 0, 3))
AIC(m1)
summary(m1)
acf(resid(m1))
plot(simulate(m1, nyears, seed = 42, future = FALSE))
points(simulate(m1, nyears, seed = 42, future = FALSE))
lines(dat2$Streamflow_std, col = "red")
y <- simulate(m1, nyears, seed = 42, future = FALSE) +
  1:nyears*1
plot(y)

#
# Equilibrium effort 
#
EMSY <- quantile(post$EMSY, 0.5)

#
# Check at equilibrium 
#

flowseed <- 30
flowsim <- simulate(m1, nsim = nyears,
                    seed = flowseed, future = FALSE) %>%
  as.numeric()
flowsim <- (flowsim - mean(flowsim))/sd(flowsim)
 
flowsim <- flowsim + (0:(nyears - 1)*(-0.1))
 plot(flowsim)
 lines(dat2$Streamflow_std)

 eseries <- rep(EMSY * effort_rel[2], nyears)

x <- simmod(100, post, nyears, 
            init_fraction1, eseries,
            flowsim, 100)
datplot <- data.frame(n = x$n, lnCPUE = x$lnCPUE, 
                      ln_cpue_data = c(dat2$ln_cpue, rep(NA, 3)))
ggplot(datplot) + 
  aes(x = n) + 
  geom_line(aes(y = lnCPUE)) + 
  geom_line(aes(y = ln_cpue_data), col = "red") 

#
# Simulate historical flow series and from model 
#

flowseeds <- 1:nsims
datsim <- expand.grid(flowseed = flowseeds,
                      trend = trends,
                      effort_rel = effort_rel)


#note below code always uses the same mcmc sample for the 
# same flow seed
#residuals are random

datout <- datsim %>%
  group_by(flowseed, trend, effort_rel) %>%
  DataGLMRepeat::with_groups({
    flowsim <- simulate(m1, nsim = nyears,
                        seed = flowseed, future = FALSE) %>%
      as.numeric()
    flowsim <- flowsim + (0:(nyears - 1)*trend)
    eseries <- rep(EMSY * effort_rel, nyears)
    x <- simmod(flowseed, post, nyears,
           init_fraction1, eseries,
           flowsim, ((nmcmc*2):(nmcmc*3))[flowseed])

    x$trend <- trend
    x$effort_rel <- effort_rel
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
dsum <- dout2 %>% group_by(n, trend, effort_rel) %>%
  summarize(across(lnCPUE:nu, ~quantile(.x, 0.05, na.rm = TRUE),
                   .names = "{.col}-q05"),
            across(lnCPUE:nu, ~quantile(.x, 0.2, na.rm = TRUE),
                   .names = "{.col}-q20"),
            across(lnCPUE:nu, ~quantile(.x, 0.4, na.rm = TRUE),
                   .names = "{.col}-q40"),
            across(lnCPUE:nu, ~quantile(.x, 0.5, na.rm = TRUE),
                   .names = "{.col}-q50"),
            across(lnCPUE:nu, ~quantile(.x, 0.6, na.rm = TRUE),
                   .names = "{.col}-q60"),
            across(lnCPUE:nu, ~quantile(.x, 0.8, na.rm = TRUE),
                   .names = "{.col}-q80"),
            across(lnCPUE:nu, ~quantile(.x, 0.95, na.rm = TRUE),
                   .names = "{.col}-q95")) %>%
  pivot_longer(4:53, names_to = "Variable",
               values_to = "Val") %>%
  separate(Variable, into = c("Var", "Quant"), sep = "-") %>%
  pivot_wider(names_from = "Quant",
              values_from = "Val")
  

get_td(dsum, "q20", trends[1], "lnCPUE_mean")
get_td(dsum, "q20", trends[1], "ndvi_mean")
get_td(dsum, "q20", trends[1], "pasture_mean")
get_td(dsum, "q20", trends[1], "SP")
get_td(dsum, "q20", trends[1], "nu")

dloop <- expand.grid(Var = c("lnCPUE_mean", "ndvi_mean", "pasture_mean", "nu"),
                      trends = trends[trends<0], 
                     effort_rel = effort_rel,
                      qs = c("q05", "q20", "q40"))

dloop$td <- lapply(1:nrow(dloop), function(x) 
  get_td_effort(dsum, dloop$qs[x],dloop$trends[x], 
         dloop$Var[x], dloop$effort_rel[x])) %>%
  as.numeric()

td_dat <- dloop
td_dat$Quantile <- as.numeric(substr(td_dat$qs, start = 2, stop= 5))

# ------------ 
# Plots 
# ------------ 
dsum2 <- filter(dsum, Var != "SP") %>%
  mutate(Trend = exp(trend))

dsum3 <- dsum2 %>%
  filter(Var %in% c("lnCPUE_mean", "ndvi_mean", "pasture_mean", "nu")) 
dsum3$Var <- factor(dsum3$Var, labels = c("ln(CPUE)", "NDVI", 
                                          "Latent condition", 
"Pasture"))

exp(filter(dsum3, n == 30 & Var == "ln(CPUE)" & trend == 0)$q80)

#units for fishing are tonnes of fish per year of fishing (since I 
# divide days by 365)
dsum3_temp <- dsum3 %>% 
  filter(trend %in% trends[(length(trends)-4):(length(trends))]) %>%
  filter(effort_rel == effort_rel[1])

dsum3_tempa <- filter(dsum3_temp, Var == "ln(CPUE)")
g1a <-  ggplot(dsum3_tempa) + 
  aes(x = n, y = exp(q50), color = Trend, group = Trend) + 
  geom_line() + 
  geom_ribbon(data = filter(dsum3_tempa, trend == 0), 
              aes(ymin = exp(q05), ymax = exp(q95)), 
              alpha = 0.25,
              color = NA) + 
  geom_ribbon(data = filter(dsum3_tempa, trend == 0),
              aes(ymin = exp(q40), ymax = exp(q60)), 
              alpha = 0.25,
              color = NA) + 
  geom_ribbon(data = filter(dsum3_tempa, trend == 0),
              aes(ymin = exp(q20), ymax = exp(q80)), 
              alpha = 0.25,
              color = NA) + 
  theme_classic() + 
  guides(color = FALSE) + 
  xlab("Year") + 
  ylab("Catch per unit effort \n (Tonnes per fishing year)")
g1a

dsum3_tempb <- filter(dsum3_temp, Var == "NDVI")
g1b <-  ggplot(dsum3_tempb) + 
  aes(x = n, y = q50, color = Trend, group = Trend) + 
  geom_line() + 
  geom_ribbon(data = filter(dsum3_tempb, trend == 0), 
              aes(ymin = q05, ymax = q95), 
              alpha = 0.25,
              color = NA) + 
  geom_ribbon(data = filter(dsum3_tempb, trend == 0),
              aes(ymin = q40, ymax = q60), 
              alpha = 0.25,
              color = NA) + 
  geom_ribbon(data = filter(dsum3_tempb, trend == 0),
              aes(ymin = q20, ymax = q80), 
              alpha = 0.25,
              color = NA) + 
  theme_classic() + 
  guides(color = FALSE) + 
  xlab("Year") + 
  ylab("NDVI")
g1b


dsum3_tempc <- filter(dsum3_temp, Var == "Pasture")
g1c <-  ggplot(dsum3_tempc) + 
  aes(x = n, y = q50, color = Trend, group = Trend) + 
  geom_line() + 
  geom_ribbon(data = filter(dsum3_tempc, trend == 0), 
              aes(ymin = q05, ymax = q95), 
              alpha = 0.25,
              color = NA) + 
  geom_ribbon(data = filter(dsum3_tempc, trend == 0),
              aes(ymin = q40, ymax = q60), 
              alpha = 0.25,
              color = NA) + 
  geom_ribbon(data = filter(dsum3_tempc, trend == 0),
              aes(ymin = q20, ymax = q80), 
              alpha = 0.25,
              color = NA) + 
  theme_classic() + 
  labs(color = "Trend \n (multiples /year)") +
  xlab("Year") + 
  ylab("Pasture biomass")
g1c

g1 <- g1a + g1b + g1c + 
  plot_annotation(tag_levels = 'A')



# Latent condition emergence plot, not included in paper
# (its just the average of the other indicators)
dsum3_tempd <- filter(dsum3_temp, Var == "Latent condition")
g1d <-  ggplot(dsum3_tempd) + 
  aes(x = n, y = q50, color = Trend, group = Trend) + 
  geom_line() + 
  geom_ribbon(data = filter(dsum3_tempd, trend == 0), 
              aes(ymin = q05, ymax = q95), 
              alpha = 0.25,
              color = NA) + 
  geom_ribbon(data = filter(dsum3_tempd, trend == 0),
              aes(ymin = q40, ymax = q60), 
              alpha = 0.25,
              color = NA) + 
  geom_ribbon(data = filter(dsum3_tempd, trend == 0),
              aes(ymin = q20, ymax = q80), 
              alpha = 0.25,
              color = NA) + 
  theme_classic() + 
  labs(color = "Trend \n (multiples /year)") +
  xlab("Year") + 
  ylab("Latent condition")
g1d


#
# Emergence time
#
td_dat$Quantile <- td_dat$qs
td_dat$Quantile <- factor(td_dat$qs, labels = c("5%", "20%", "40%"))
td_dat$Var <- factor(td_dat$Var, labels = c("CPUE", "NDVI", "Pasture",
                                            "Latent condition"))

#Filter out latent condition (its just the average
# emergence time of the other indicators)
g2 <- td_dat %>%
  filter(effort_rel == effort_rel[1] & Var != "Latent condition") %>%
  ggplot() + 
  aes(x = exp(trends), y = td, color = Var, 
      group = Var, linetype = Var) + 
  geom_line() + 
  facet_wrap(~Quantile) + 
  theme_classic() + 
  scale_x_log10() + 
  scale_color_manual(values = c("lightblue","steelblue", "darkblue", "black")) + 
  xlab("Trend in flow (multiples /yr)") + 
  ylab("Emergence time (years)") +
  labs(color = "Indicator", linetype = "Indicator")
g2


#
# Comparison with different effort rates 
#


g3 <- td_dat %>%
  filter(Var == "CPUE") %>%
  ggplot() + 
  aes(x = exp(trends), y = td, color = factor(effort_rel)) + 
  geom_line() + 
  facet_wrap(~Quantile) + 
  theme_classic() + 
  scale_x_log10() + 
  scale_color_manual(values = c("steelblue", "darkblue", "black")) + 
  xlab("Trend in flow (multiples /yr)") + 
  ylab("Emergence time (years)") +
  labs(color = "Relative effort")
g3

#
# Save plots 
#

ggsave("../figures/2021-08-05_emergence_plot.png", g1,
       width = 8, height = 3)
ggsave(g2, file = "../figures/2021-08-05_emergence-times.png",
       width = 6, height = 3)

ggsave(g3, file = "../figures/2021-08-05_emergence-times-increased-effort.png",
       width = 6, height = 3)






