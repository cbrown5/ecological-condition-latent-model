# Estimate effect sizes
#CJ Brown 
# 2021-11-10
#

# Calculate relative change in each indicator, for a given change in flow or condition
# Changes in condition are 1SD
# changes in flow are equivalent to a doubling of flow
# (streamflow changes by up to 15 times from its lowest value)
# For NDVI and pasture biomass, these changes are additive
# For Barramundi changes are multiplicative 



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
post <- extract.samples(fitm1) %>% data.frame()
post$imcmc <- 1:nrow(post)
probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)

n <- 100
dfnu <- expand.grid(imcmc = 1:nrow(post), 
                    nu =seq(-2, 2, length.out = n))
post2 <- select(post, imcmc, a_ndvi, beta_ndvi,
                a_pasture, beta_pasture) %>%
  inner_join(dfnu) %>%
  mutate(ndvi = exp(a_ndvi + beta_ndvi*nu),
         pasture = exp(a_pasture + beta_pasture*nu))

postsum <- post2 %>%
  group_by(nu) %>%
  summarize(ndvi_mu = quantile(ndvi, 0.5))


ggplot(postsum) + 
  aes(x = nu, y = ndvi_mu) + 
  geom_line()

#
#A one SD change in conditions results in a change of indicators of:
#

#NDVI: 
quantile(post$beta_ndvi, probs = probs)
quantile(post$beta_ndvi*post$beta_nu, probs = probs)
#log pasture
quantile(post$beta_pasture, probs = probs)
quantile(post$beta_pasture*post$beta_nu, probs = probs)

#Relative change in surplus production for a multiple change in flow
#(algebra is in my blue notebook)
relSPFlow <- function(beta_u, beta_nu, flow_mult){
  flow_mult*exp(beta_u * beta_nu)
}

#relative change in SP for an increase in u
relSPu <- function(beta_u, u_diff){
  exp(beta_u * u_diff)
}

# Quantile function for across
quants <- list(
  quant025 = ~quantile(.x, probs = 0.025)
)

quantdf <- post %>% select(beta_u, beta_nu, 
                           beta_pasture, beta_ndvi) %>%
  mutate(SPflow = relSPFlow(beta_u, beta_nu, 2),
         SPu = relSPu(beta_u, 1),
         pasture_flow = beta_pasture * beta_nu*log(2),
         ndvi_flow = beta_ndvi * beta_nu*log(2)) %>%
  select(beta_pasture, beta_ndvi,
         pasture_flow, ndvi_flow,
         SPflow,
         SPu) %>%
  pivot_longer(beta_pasture:SPu, names_to = "Var", values_to = "val") %>%
  group_by(Var) %>%
  summarize(median = quantile(val, 0.5),
            lwr025 = quantile(val, 0.025),
            upr975 = quantile(val, 0.975),
            lwr25 = quantile(val, 0.25),
            upr75 = quantile(val, 0.75)
            ) %>%
  ungroup() %>%
  mutate(influ = if_else(grepl("flow", Var), "Flow", "Condition"),
    Var2 = case_when(
    grepl("ndvi", Var) ~ "NDVI \n difference",
    grepl("pasture", Var) ~ "Pasture \n difference",
    grepl("SP", Var) ~ "Surplus prodn \n multiple"
  ))

g1 <- ggplot(filter(quantdf, influ == "Flow")) + 
  aes(x = Var2, y = median) + 
  geom_point(size = 4, alpha = 0.8, shape = 3) +
  geom_linerange(aes(ymin = lwr25, ymax = upr75), 
                 size = 3, alpha = 0.85) + 
  geom_linerange(aes(ymin = lwr025, ymax = upr975),
                 alpha = 0.6, size = 1.25) + 
  theme_classic() +
  geom_hline(yintercept =0) +
  ylab("Sensitivity \n to a doubling of flow") + 
  xlab("")
  
g2 <- ggplot(filter(quantdf, influ == "Condition")) + 
  aes(x = Var2, y = median) + 
  geom_point(size = 4, alpha = 0.8, shape = 3) +
  geom_linerange(aes(ymin = lwr25, ymax = upr75), 
                 size = 3, alpha = 0.85) + 
  geom_linerange(aes(ymin = lwr025, ymax = upr975),
                 alpha = 0.6, size = 1.25) + 
  theme_classic() +
  geom_hline(yintercept =0) +
  ylab("Sensitivity to a \n 1 S.D. increase in condition") + 
  xlab("") 

pw <- g1 + g2


pw2 <- pw  + plot_annotation(tag_levels = 'A') + 
  plot_layout(ncol = 2, widths = c(1,1)) & 
  theme(plot.tag = element_text(size = 14)) 

ggsave(pw2, file = "../figures/2021-08-06_sensitivity.png",
       width = 6, height = 3)



