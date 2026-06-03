 functions {
   real partial_sum(array[] int ind, int start, int end, 
                   vector Y,
                   vector logdilution,
                   vector loginjection_volume,
                   array[] int metabolite_number,
                   vector r_metabolite,
                   vector r_dilution,
                   vector r_injection,
                   vector r_dilution_injection_interaction,
                   vector r_log_sigma,
                   vector r_sigma_dilution,
                   vector r_sigma_injection,
                   vector r_sigma_dilution_injection_interaction,
                   vector r_nu) {
                     
    real lp = 0;

    for (n in start : end) {
    real mu = r_metabolite[metabolite_number[n]]+
    r_dilution[metabolite_number[n]]*logdilution[n]+
    r_injection[metabolite_number[n]]*loginjection_volume[n]+
    r_dilution_injection_interaction[metabolite_number[n]]*loginjection_volume[n]*logdilution[n];
    real sigma = exp(r_log_sigma[metabolite_number[n]]+
    r_sigma_dilution[metabolite_number[n]]*logdilution[n]+
    r_sigma_injection[metabolite_number[n]]*loginjection_volume[n]+
    r_sigma_dilution_injection_interaction[metabolite_number[n]]*loginjection_volume[n]*logdilution[n]);
    real nu = r_nu[metabolite_number[n]];
    lp = lp + student_t_lpdf(Y[n] | nu, mu, sigma);
      
    }
    return lp;
  }
  
 
  real prop_gt(vector x, real c) {
    real count = 0;
    for (k in 1:num_elements(x))
      count += x[k] > c;
    return count / num_elements(x);
  }

}

data {
  int<lower=1> N;     // total number of observations
  vector[N] Y;        // response variable
  vector[N] dilution_factor;   // logdilutions
  vector[N] injection_volume;   // logdilutions
  int<lower=1> number_metabolites ; // total number of metabolits
  array[N] int metabolite_number;   // metabolite indeces
  int prior_only;     // should the likelihood be ignored?
}

transformed data {
  int grainsize = 1;
  array[N] int ind = rep_array(1, N);
  vector[N] logdilution = log(dilution_factor/0.5);   // logdilutions
  vector[N] loginjection_volume = log(injection_volume/5);   // logdilutions
  
  vector[10] dilution_factor_sim = [0.2500000, 0.3333333, 0.4166667, 0.5000000, 0.5833333, 0.6666667, 0.7500000, 0.8333333, 0.9166667, 1.0000000]';
  vector[10] injection_volume_sim = [2.500000, 3.333333,  4.166667,  5.000000,  5.833333,  6.666667,  7.500000,  8.333333,  9.166667, 10.000000]';
  vector[10] logdilution_factor_sim = log(dilution_factor_sim/0.5);
  vector[10] loginjection_volume_sim = log(injection_volume_sim/5);
  matrix[10,10] Xd;
  matrix[10,10] Yi;

  for (i in 1:10)
    for (j in 1:10) {
      Xd[i,j] = logdilution_factor_sim[i];
      Yi[i,j] = loginjection_volume_sim[j];
    }
  
}

parameters {
  real mean_metabolite;  // regression coefficients
  real mean_dilution; 
  real mean_injection;
  real mean_dilution_injection_interaction;
  real mean_sigma_dilution; 
  real mean_sigma_injection;
  real mean_sigma_dilution_injection_interaction;
  real mean_log_sigma;  
  real<lower=0> sd_metabolite;  // group-level standard deviations
  real<lower=0> sd_dilution;  
  real<lower=0> sd_injection; 
  real<lower=0> sd_dilution_injection_interaction;
  real<lower=0> sd_log_sigma; 
  real<lower=0> sd_sigma_dilution;  
  real<lower=0> sd_sigma_injection; 
  real<lower=0> sd_sigma_dilution_injection_interaction;
  vector[number_metabolites] r_metabolite;  // actual group-level effects
  vector[number_metabolites] r_dilution;
  vector[number_metabolites] r_injection;
  vector[number_metabolites] r_dilution_injection_interaction;
  vector[number_metabolites] r_log_sigma;
  vector[number_metabolites] r_sigma_dilution;
  vector[number_metabolites] r_sigma_injection;
  vector[number_metabolites] r_sigma_dilution_injection_interaction;
  vector<lower=3>[number_metabolites] r_nu;  // degrees of freedom or shape
}

transformed parameters {

  real lprior = 0;  // prior contributions to the log posterior

  lprior += normal_lpdf(mean_metabolite | 11, 2); // mean 11.12239, scale 1.734539
  lprior += normal_lpdf(mean_dilution | 1, 1);    // about half the scale of the data
  lprior += normal_lpdf(mean_injection | 1, 1);
  lprior += normal_lpdf(mean_dilution_injection_interaction | 0, 0.25); // 1/4 of main
  
  lprior += normal_lpdf(mean_log_sigma | -1.609438, 0.5); //~20% CV:  ~5%-40%
  lprior += normal_lpdf(mean_sigma_dilution | 0, 0.25);  // close to zero
  lprior += normal_lpdf(mean_sigma_injection | 0, 0.25); // close to zero
  lprior += normal_lpdf(mean_sigma_dilution_injection_interaction | 0, 0.0625); // close to zero

  lprior += normal_lpdf(sd_metabolite | 0, 2) - normal_lccdf(0 | 0, 2);
  lprior += normal_lpdf(sd_dilution | 0, 1) -  normal_lccdf(0 | 0, 1);
  lprior += normal_lpdf(sd_injection | 0, 1) -  normal_lccdf(0 | 0, 1);
  lprior += normal_lpdf(sd_dilution_injection_interaction | 0, 0.25) -  normal_lccdf(0 | 0, 0.25);
  
  lprior += normal_lpdf(sd_log_sigma | 0, 0.5) -  normal_lccdf(0 | 0, 0.5);
  lprior += normal_lpdf(sd_sigma_dilution | 0, 0.25) -  normal_lccdf(0 | 0, 0.25);
  lprior += normal_lpdf(sd_sigma_injection | 0, 0.25) -  normal_lccdf(0 | 0, 0.25);
  lprior += normal_lpdf(sd_sigma_dilution_injection_interaction | 0, 0.0625) -  normal_lccdf(0 | 0, 0.0625);
  
  lprior += gamma_lpdf(r_nu | 2, 0.1) - gamma_lccdf(3 | 2, 0.1);
}
model {
 
// likelihood 
  if (!prior_only) {

   target += normal_lpdf(r_metabolite | mean_metabolite, sd_metabolite);
   target += normal_lpdf(r_dilution | mean_dilution, sd_dilution);
   target += normal_lpdf(r_injection | mean_injection, sd_injection);
   target += normal_lpdf(r_dilution_injection_interaction | mean_dilution_injection_interaction, sd_dilution_injection_interaction);
   target += normal_lpdf(r_sigma_dilution | mean_sigma_dilution, sd_sigma_dilution);
   target += normal_lpdf(r_sigma_injection | mean_sigma_injection, sd_sigma_injection);
   target += normal_lpdf(r_sigma_dilution_injection_interaction | mean_sigma_dilution_injection_interaction, sd_sigma_dilution_injection_interaction);
   target += normal_lpdf(r_log_sigma | mean_log_sigma, sd_log_sigma);
  // target += student_t_lpdf(Y | nu, mu, sigma);
  
    target += reduce_sum(partial_sum, ind, grainsize, Y, logdilution, loginjection_volume,
                         metabolite_number, r_metabolite, r_dilution, r_injection, r_dilution_injection_interaction, r_log_sigma, r_sigma_dilution, r_sigma_injection, r_sigma_dilution_injection_interaction, r_nu);
   }
  // priors
  target += lprior;
}
generated quantities {
  array[10,10,2] real SVR;
  array[10,10,2] real pSVR;
 for (i in 1:10) {    
   for (j in 1:10) {
    vector[number_metabolites] svr_temp;
     for (n in 1:number_metabolites) {
       real slope = r_dilution[n]+r_injection[n]+r_dilution_injection_interaction[n]*(Yi[i,j]+Xd[i,j]);
       real sd = exp(r_log_sigma[n]+r_sigma_dilution[n]*Xd[i,j]+r_sigma_injection[n]*Yi[i,j]+
       r_sigma_dilution_injection_interaction[n]*Yi[i,j]*Xd[i,j])*sqrt(r_nu[n]/(r_nu[n]-2));
       svr_temp[n] = slope/2/sd;
       }
    SVR[i,j,1] = mean(svr_temp[1:184]);
    SVR[i,j,2] = mean(svr_temp[185:number_metabolites]);
    pSVR[i,j,1] = prop_gt(svr_temp[1:184], 2/log(2));
    pSVR[i,j,2] = prop_gt(svr_temp[185:number_metabolites], 2/log(2));
    }}
    
array[3] int idx = {1,4,10};

array[3,3] vector[number_metabolites] svrm;
 for (i in 1:3 ) {    
   for (j in 1:3) {
     for (n in 1:number_metabolites) {
       real slope = r_dilution[n]+r_injection[n]+r_dilution_injection_interaction[n]*(Yi[idx[i],idx[j]]+Xd[idx[i],idx[j]]);
       real sd = exp(r_log_sigma[n]+r_sigma_dilution[n]*Xd[idx[i],j]+r_sigma_injection[n]*Yi[idx[i],idx[j]]+
       r_sigma_dilution_injection_interaction[n]*Yi[idx[i],idx[j]]*Xd[idx[i],idx[j]])*sqrt(r_nu[n]/(r_nu[n]-2));
       svrm[i,j,n] = slope/2/sd;
       }
    }}
}

