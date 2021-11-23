//Indicator model 
// CJ Brown 2021-11-10

data{
  //all variables 
  int<lower=1> N; 
  vector[N] flow; 
  
  //stock model 
  vector[N] lnCPUE; 
  vector<lower=0>[N] catches; 
  
  //Stock param settings
  real logK_mean;
  real<lower=0> logK_sd;
  real logr_mean;
  real<lower=0> logr_sd;
  real<lower=0> init_fraction;
  
  //
  //other indicators 
  //
  
  //Sample sizes
  int<lower = 1> Nndvi; 
  int<lower = 1> Npasture;
  
  //Data
  vector[Nndvi] ndvi;
  vector[Npasture] pasture;
  
  //Indicator for years
  int<lower = 1> i_ndvi [Nndvi];
  int<lower = 1> i_pasture [Npasture];
  
}

parameters{
  //LV
  real<lower=0> beta_nu;
  vector[N] nu_raw;
  
  //stock model 
  real lnK;
  real lnr;
  real<lower = 0.001> q;
  real beta_u;
  real<lower=0> sigma_u;
  vector[N] u_raw; 
  real<lower=0> sigma_cpue;
  
  //Other indicators
  real beta_ndvi; //regression on latent 
  real beta_pasture;
  real a_ndvi; //intercepts
  real a_pasture;
  
  real<lower=0> sigma_ndvi;
  real<lower=0> sigma_pasture;
  
}

transformed parameters{
 real<lower=0> K;
 real<lower=0> r;
 real Btemp;
 vector<lower=0>[N] B; 
 vector[N] u;
 vector[N] lnCPUE_hat;
 vector[N] nu_hat;
 vector[N] nu;
 vector[Nndvi] ndvi_hat;
 vector[Npasture] pasture_hat;
 
 
 //
 // Latent variables
 //
 nu_hat = beta_nu * flow; //flow relationship
 nu = nu_hat + nu_raw; //implies nu ~ normal(nu_hat, 1)
 
 //
 //Biomass model
 //
 K = exp(lnK);
 r = exp(lnr);
 u = beta_u * nu + u_raw*sigma_u; //implies u ~ normal(beta_u*nu, sigma_u)

 B[1] = K * init_fraction;
 
 for (t in 2:N){
      Btemp = (B[t-1] + r*B[t-1]*(1-(B[t-1]/K)) - catches[t-1]) * exp(u[t]);
  
   if (Btemp < 0.001){
      B[t] = 0.001;
   } else {
    B[t] = Btemp;
   }
 }
 
  lnCPUE_hat = log(q) + log(B);
 
 //
 // Latent regressions
 //
 ndvi_hat = beta_ndvi * nu[i_ndvi] + a_ndvi;
 pasture_hat = beta_pasture * nu[i_pasture] + a_pasture;

}

model{ 

  //
  //Ecological condition 
  //
  nu_raw ~ std_normal();
  beta_nu ~ exponential(1); //informed prior to prevent overfitting/inseparability 

  //
  // Stock params 
  //
  lnK ~ normal(logK_mean, logK_sd);
  lnr ~ normal(logr_mean, logr_sd);
  q ~ uniform(0.001, 0.15); //v2
  beta_u ~ normal(0, 0.2);
  
  // process error
  sigma_u ~ exponential(14);
  u_raw ~ std_normal();
  
  //Indicator params 
  a_ndvi ~ normal(0,10);
  a_pasture ~ normal(0,10);
  
  beta_ndvi ~ normal(0,1);
  beta_pasture ~ normal(0,1);
  
  
  //Observation errors
  sigma_cpue ~  exponential(2);
  sigma_ndvi ~  exponential(0.1);
  sigma_pasture ~  exponential(0.1);
  
  // Observations
  lnCPUE ~ normal(lnCPUE_hat, sigma_cpue);
  ndvi ~ normal(ndvi_hat, sigma_ndvi);
  pasture ~ normal(pasture_hat, sigma_pasture);
}

generated quantities{
  real<lower = 0> MSY;
  real<lower = 0> FMSY;
  real<lower = 0> EMSY;
  real<lower = 0> BMSY;
  vector<lower=0>[N] Brel; 
  vector[N] ll_cpue;
  vector[Nndvi] ll_ndvi;
  vector[Npasture] ll_pasture;
  vector[N] ndvi_pred; 
  vector[N] pasture_pred; 
  

  MSY = exp(lnr)*K/4;
  EMSY =  exp(lnr)/(2*q);
  FMSY =  exp(lnr)/2;
  BMSY = K/2;
  Brel = B/K;
  
  ndvi_pred = beta_ndvi * nu + a_ndvi;
  pasture_pred = beta_pasture * nu + a_pasture;
  
  
   // Log liks 
  for (n in 1:N){
    ll_cpue[n] = normal_lpdf( lnCPUE[n]|lnCPUE_hat[n], sigma_cpue);
  }
  
  for (n in 1:Nndvi){
    ll_ndvi[n] = normal_lpdf( ndvi[n]|ndvi_hat[n], sigma_ndvi);
  }

  for (n in 1:Npasture){
    ll_pasture[n] = normal_lpdf( pasture[n]|pasture_hat[n], sigma_pasture);
  }
  
}

  

