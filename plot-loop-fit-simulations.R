# Plots for fits of model to simulated data 
#CJ Brown 
# 2020-11-10
#

library(ggplot2)
library(dplyr)
library(forecast)
library(DataGLMRepeat)
library(cowplot)
library(rethinking)
library(tidyr)
library(rstan)

load("../models/2021-11-24_param-sims.rda")

load("../models/2021-08-03_model-fit-no-inun-logflow-Binit20perc-q0.rda")
nsims <- length(datout)

#
# Extract and plot parameter coverage 
#

 params <- lapply(1:nsims, function(x) {
   datout[[x]]$paramsum})
 
params <- do.call("rbind", params) %>%
  mutate(covered = (truth >= lwr) & (truth <= upr),
         covered50 = (truth >= lwr50) & (truth <= upr50),
         error = med-truth)

#
#Coverage table 
#

paramdf <- params %>%
  filter(Rhat < 1.05) %>% 
  group_by(param, yrs_missing) %>%
  summarize(n = n(),
            coverage = sum(covered)/n,
            coverage50 = sum(covered50)/n) %>%
  filter(!is.na(coverage))
paramdf$yrs_missing <- paramdf$yrs_missing - 1
write.csv(data.frame(paramdf), "figures/coverage-stats.csv", 
          row.names = FALSE)

#
# R2 values 
#

R2 <- lapply(1:nsims, function(x) {
  datout[[x]]$fitcors}) %>% 
  do.call("rbind", .) %>%
  pivot_longer(ndvi:nu, names_to = "Variable",
               values_to = "R2") %>%
  group_by(yrs_missing, Variable) %>%
  summarize(meanR2 = mean(R2),
            sdR2 = sd(R2))

gp1 <- ggplot(R2) + 
  aes(x = yrs_missing, y = meanR2) + 
  geom_point(color = "grey", size = 2) + 
  geom_linerange(aes(ymin = meanR2 - sdR2, ymax = meanR2 + sdR2),
                 size = 1, color = "grey") + 
  facet_wrap(~Variable) +
  theme_classic() + 
  ylab(bquote(~R^2)) + 
  xlab("Years of missing data") + 
  ylim(0,1)
gp1

ggsave(gp1, file = "../figures/2021-11-29_R2-plot-sensitivity.png",
       width = 7, height = 3.5)

