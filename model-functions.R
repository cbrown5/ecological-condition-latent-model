#Functions for the model 
# CJ Brown 2021-10-11


#Function simulate from posterior vals
simmod <- function(i, post, N, init_fraction, effort, flow, iseed){
  with(post[i,], {
    set.seed(iseed)
    K <- exp(lnK)
    nu_mean <- beta_nu * flow #flow relationship
    nu <- rnorm(N, nu_mean, 1)
    r <- exp(lnr)
    
    # Latent regressions
    pasture_mean <- beta_pasture * nu + a_pasture
    ndvi_mean <- beta_ndvi * nu + a_ndvi
    u <- rnorm(N, beta_u * nu, sigma_u)
    
    #Biomass model
    B <- catches <- numeric(N)
    B[1] <- K * init_fraction
    catches[1] <- effort[1] * q * B[1]
    
    for (t in 2:N){
      Btemp <- (B[t-1] + r*B[t-1]*(1-(B[t-1]/K)) - catches[t-1])*exp(u[t])
      
      if (Btemp < 0.001) 
        Btemp <- 0.001
      
      B[t] <- Btemp
      catches[t] <- effort[t]*q*B[t]
    }
    
    lnCPUE_mean <- log(q) + log(B) 
    
    #Surplus production
    SP <- c(diff(B) + catches[1:(N-1)], NA)
    #
    # Observation models
    #
    ndvi <- rnorm(N,ndvi_mean, sigma_ndvi)
    pasture <- rnorm(N,pasture_mean, sigma_pasture)
    
    
    # sampling 
    lnCPUE <- rnorm(N, lnCPUE_mean, sigma_cpue)
    return(data.frame(n = 1:N,
                      iter = i,
                      lnCPUE = lnCPUE, 
                      ndvi = ndvi,
                      SP = SP, 
                      pasture = pasture,
                      lnCPUE_mean = lnCPUE_mean,
                      ndvi_mean = ndvi_mean,
                      pasture_mean = pasture_mean,
                      nu = nu,
                      catches = catches,
                      flow = flow))
  })
  
}

#
#Simulate model with random flow series
#
simflow <- function(i, post, N, init_fraction, effort,
                    fmod, year_int, #year flow loss starts
                    flow_loss, #absolute flow loss
                    iseedflow, #seed for flow simulation
                    iseedsim,#seed for model sims
                    flowsd,#to standardize 
                    flowmn
){
  
  f1 <- stats::simulate(fmod, seed = iseedflow)*flow_loss
  f1 <- (f1[1:N] - flowmn)/flowsd
  f1[year_int:N] <- f1[year_int:N]
  simmod(i, post, N1, init_fraction1, effort, f1, iseedsim)
}


#Function to get emergence times for a given variable,
#trend and quantile. Designed to work with the 'dsum' dataframe
get_td <- function(x, this_quant, this_trend, this_var){
  x1 <- dplyr::filter(x, trend == 0 & Var == this_var)
  x2 <- dplyr::filter(x, trend == this_trend & Var == this_var)
  if (this_trend > 0 ){
    i <- which(x2$q50 > x1[,this_quant])
  } else {
    i <- which(x2$q50 < x1[,this_quant])
  }
  min(x2$n[i])
}

#as above, but with different efforts 
get_td_effort <- function(x, this_quant, this_trend, this_var, this_effort){
  x1 <- dplyr::filter(x, trend == 0 & Var == this_var & effort_rel == this_effort)
  x2 <- dplyr::filter(x, trend == this_trend & Var == this_var & effort_rel == this_effort)
  if (this_trend > 0 ){
    i <- which(x2$q50 > x1[,this_quant])
  } else {
    i <- which(x2$q50 < x1[,this_quant])
  }
  min(x2$n[i])
}


