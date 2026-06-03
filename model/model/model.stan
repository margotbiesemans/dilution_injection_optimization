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
  
  lprior += normal_lpdf(mean_log_sigma | -1.609438, 0.5); //~20% CV:  ~5%-50%
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
generated quantities {}

