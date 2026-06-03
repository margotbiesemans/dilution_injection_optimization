 functions {
   
  real lkj_corr_cholesky_point_lower_tri2_lpdf(matrix cor_L, real point_mu_lower, real point_scale_lower) {
    real lpdf = lkj_corr_cholesky_lpdf(cor_L | 1);
    int d = rows(cor_L);
    matrix[d,d] cor = multiply_lower_tri_self_transpose(cor_L);
    lpdf += normal_lpdf(cor[2,1]  | point_mu_lower, point_scale_lower);
    return(lpdf);
 }
  
  real partial_sum(array[] int ind, int start, int end, 
                   vector Y,
                   vector logdilution,
                   vector loginjection_volume,
                   array[] int metabolite_number,
                   array[] int chrom_number,
                   array[] int  posnegdata,
                   vector r_metabolite,
                   vector r_dilution,
                   vector r_injection,
                   vector r_dilution_injection_interaction,
                   vector r_log_sigma,
                   vector r_sigma_dilution,
                   vector r_sigma_injection,
                   vector r_sigma_dilution_injection_interaction,
                   vector r_nu,
                   matrix r_chrom,
                   vector sigma_shared)  {
                     
    real lp = 0;

    for (n in start : end) {
    real mu = r_metabolite[metabolite_number[n]]+
    r_dilution[metabolite_number[n]]*logdilution[n]+
    r_injection[metabolite_number[n]]*loginjection_volume[n]+
    r_dilution_injection_interaction[metabolite_number[n]]*loginjection_volume[n]*logdilution[n]+
    sigma_shared[posnegdata[n]]*r_chrom[chrom_number[n],posnegdata[n]];
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
  
  int<lower=1> number_chroms ; // total number of metabolits
  array[N] int chrom_number;   // metabolite indeces
  array[N] int posnegdata;   // metabolite indeces
  vector[number_metabolites] posneg;
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
  
  real dmean_metabolite;  // regression coefficients
  real dmean_dilution; 
  real dmean_injection;
  real dmean_dilution_injection_interaction;
  real dmean_sigma_dilution; 
  real dmean_sigma_injection;
  real dmean_sigma_dilution_injection_interaction;
  real dmean_log_sigma;  
  
  real<lower=0> sd_metabolite;  // group-level standard deviations
  real<lower=0> sd_dilution;  
  real<lower=0> sd_injection; 
  real<lower=0> sd_dilution_injection_interaction;
  real<lower=0> sd_log_sigma; 
  real<lower=0> sd_sigma_dilution;  
  real<lower=0> sd_sigma_injection; 
  real<lower=0> sd_sigma_dilution_injection_interaction;
 
 cholesky_factor_corr[2] L_ms;
 cholesky_factor_corr[2] L_di;
 
  matrix[number_metabolites, 2] r_ms;  // actual group-level effects

  sum_to_zero_vector[number_metabolites] z_di1, z_di2;  // actual group-level effects
  sum_to_zero_vector[number_metabolites] z_dilution_injection_interaction;
  vector[number_metabolites] r_sigma_dilution;
  vector[number_metabolites] r_sigma_injection;
  vector[number_metabolites] r_sigma_dilution_injection_interaction;
  
  sum_to_zero_vector[number_chroms] z_chrom1, z_chrom2;
  vector<lower=0>[2] sigma_shared;   // strength of common variation

  vector<lower=3>[number_metabolites] r_nu;  // degrees of freedom or shape
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  vector[number_metabolites] r_dilution_injection_interaction = mean_dilution_injection_interaction *(1-posneg)+
                                                                dmean_dilution_injection_interaction*posneg + 
                                                                sd_dilution_injection_interaction*z_dilution_injection_interaction;
  matrix[number_metabolites, 2] r_di;
  
  for (i in 1:number_metabolites) {
  {
  vector[2] mu;
  vector[2] z;
  mu[1] = mean_dilution *(1-posneg[i])  + dmean_dilution  * posneg[i];
  mu[2] = mean_injection*(1-posneg[i]) + dmean_injection * posneg[i];
  z[1] = z_di1[i];
  z[2] = z_di2[i];
  r_di[i,] = (mu + diag_pre_multiply([sd_dilution, sd_injection], L_di) * z)';
  }
  } 
  
  lprior += normal_lpdf(mean_metabolite | 11, 2); // mean 11.12239, scale 1.734539
  lprior += normal_lpdf(mean_dilution | 1, 1);    // about half the scale of the data
  lprior += normal_lpdf(mean_injection | 1, 1);
  lprior += normal_lpdf(mean_dilution_injection_interaction | 0, 0.25); // 1/4 of main
  
  lprior += normal_lpdf(mean_log_sigma | -1.609438, 0.5); //~20% CV:  ~5%-50%
  lprior += normal_lpdf(mean_sigma_dilution | 0, 0.25);  // close to zero
  lprior += normal_lpdf(mean_sigma_injection | 0, 0.25); // close to zero
  lprior += normal_lpdf(mean_sigma_dilution_injection_interaction | 0, 0.0625); // close to zero

  lprior += normal_lpdf(dmean_metabolite | 11, 2);
  lprior += normal_lpdf(dmean_dilution | 1, 1);    
  lprior += normal_lpdf(dmean_injection | 1, 1);
  lprior += normal_lpdf(dmean_dilution_injection_interaction | 0, 0.25);

  lprior += normal_lpdf(dmean_log_sigma | -1.609438, 0.5); 
  lprior += normal_lpdf(dmean_sigma_dilution | 0, 0.25);  
  lprior += normal_lpdf(dmean_sigma_injection | 0, 0.25); 
  lprior += normal_lpdf(dmean_sigma_dilution_injection_interaction | 0, 0.0625); 

  lprior += normal_lpdf(sd_metabolite | 0, 2) - normal_lccdf(0 | 0, 2);
  lprior += normal_lpdf(sd_dilution | 0, 1) -  normal_lccdf(0 | 0, 1);
  lprior += normal_lpdf(sd_injection | 0, 1) -  normal_lccdf(0 | 0, 1);
  lprior += normal_lpdf(sd_dilution_injection_interaction | 0, 0.25) -  normal_lccdf(0 | 0, 0.25);
  
  lprior += normal_lpdf(sd_log_sigma | 0, 0.5) -  normal_lccdf(0 | 0, 0.5);
  lprior += normal_lpdf(sd_sigma_dilution | 0, 0.25) -  normal_lccdf(0 | 0, 0.25);
  lprior += normal_lpdf(sd_sigma_injection | 0, 0.25) -  normal_lccdf(0 | 0, 0.25);
  lprior += normal_lpdf(sd_sigma_dilution_injection_interaction | 0, 0.0625) -  normal_lccdf(0 | 0, 0.0625);

  lprior +=normal_lpdf(sigma_shared | 0, 1) - normal_lccdf(0 | 0, 1);
  
}
model {
 
// likelihood 
  if (!prior_only) {


   for (i in 1:number_metabolites) {
   target +=multi_normal_cholesky_lpdf(r_ms[i,]' | [mean_metabolite*(1-posneg[i])+dmean_metabolite*posneg[i],
                                                    mean_log_sigma*(1-posneg[i])+dmean_log_sigma*posneg[i]], diag_pre_multiply([sd_metabolite, sd_log_sigma], L_ms));
   // target +=multi_normal_cholesky_lpdf(r_di[i,]' | [mean_dilution+dmean_dilution*posneg[i],
   //                                                  mean_injection+dmean_injection*posneg[i]], diag_pre_multiply([sd_dilution, sd_injection], L_di));
 }      
      
   target += lkj_corr_cholesky_lpdf(L_ms | 2);
   target += lkj_corr_cholesky_point_lower_tri2_lpdf(L_di | 0.9, 0.05);
                                               
   target += normal_lpdf(z_dilution_injection_interaction | 0 , sd_dilution_injection_interaction);
   target += normal_lpdf(r_sigma_dilution | mean_sigma_dilution*(1-posneg)+dmean_sigma_dilution*posneg, sd_sigma_dilution);
   target += normal_lpdf(r_sigma_injection | mean_sigma_injection*(1-posneg)+dmean_sigma_injection*posneg, sd_sigma_injection);
   target += normal_lpdf(r_sigma_dilution_injection_interaction | mean_sigma_dilution_injection_interaction*(1-posneg)+dmean_sigma_dilution_injection_interaction*posneg, sd_sigma_dilution_injection_interaction);
 
   target += gamma_lpdf(r_nu | 2, 0.1)- gamma_lccdf(3 | 2, 0.1);
  // target += student_t_lpdf(Y | nu, mu, sigma);
  
   target += normal_lpdf(z_chrom1 | 0, 1);
   target += normal_lpdf(z_chrom2 | 0, 1);
   
   target += normal_lpdf(z_di1 | 0, 1);
   target += normal_lpdf(z_di2 | 0, 1);

   target += reduce_sum(partial_sum, ind, grainsize, Y, logdilution, loginjection_volume,
                         metabolite_number, chrom_number, posnegdata, r_ms[1:number_metabolites,1], r_di[1:number_metabolites,1], r_di[1:number_metabolites,2], r_dilution_injection_interaction, r_ms[1:number_metabolites,2], r_sigma_dilution, r_sigma_injection, r_sigma_dilution_injection_interaction, r_nu, append_col(z_chrom1, z_chrom2),sigma_shared);
   

   }
  // priors
  target += lprior;
}
generated quantities {
  
  corr_matrix[2] rho_ms;
  corr_matrix[2] rho_di;
  vector[number_metabolites] r_dilution = r_di[1:number_metabolites,1];
  vector[number_metabolites] r_injection = r_di[1:number_metabolites,2];
  vector[number_metabolites] r_metabolite = r_ms[1:number_metabolites,1];
  vector[number_metabolites] r_log_sigma = r_ms[1:number_metabolites,2];
   
  rho_ms = multiply_lower_tri_self_transpose(L_ms);
  rho_di = multiply_lower_tri_self_transpose(L_di);
}

