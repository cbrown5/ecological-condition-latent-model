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

load("../models/2021-06-14_param-sims.rda")

load("../models/2020-12-11_model-fit-no-inun-logflow.rda")
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

#Coverage table 
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

p2 <- params %>%
  filter(param %in% c("beta_u", "beta_pasture", "beta_ndvi", "beta_nu")) %>%
  filter(Rhat < 1.05) %>% 
  group_by(param, yrs_missing) %>% 
  arrange(mn, .by_group = TRUE) %>%
  mutate(id = seq_len(n())) %>%
  ungroup()
p2$yrs_missing2 <- p2$yrs_missing-1
ptruth <- p2 %>% select(param, truth) %>% distinct()

gp1 <- ggplot(p2) + 
  aes(x = id, y = mn) + 
  geom_point(color = "grey", size = 0.9) + 
  geom_linerange(aes(ymin = lwr50, ymax = upr50),size = 1, color = "grey") + 
  geom_linerange(aes(ymin = lwr, ymax = upr), color = "grey") + 
  facet_grid(param ~ yrs_missing2, scales = "free") + 
  geom_hline(data = ptruth, aes(yintercept = truth), linetype = 2) + 
  geom_hline(aes(yintercept = 0)) + 
  theme_classic() + 
  ylab("Parameter estimate") + 
  xlab("Rank order")
gp1

ggsave(gp1, file = "../figures/2021-06-11_coverage-plot.png",
       width = 7, height = 3.5)

