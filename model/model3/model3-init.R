# Create Stan initial values
#
# This function must return something that can be passed to the `init` argument
#   of `cmdstanr::sample()`. There are several options; see `?cmdstanr::sample`
#   for details.
#
# `.data` represents the list returned from `make_standata()` for this model.
#   This is provided in case any of your initial values are dependent on some
#   aspect of the data (e.g. the number of rows).
#
# `.args` represents the list of attached arguments that will be passed through to
#   cmdstanr::sample(). This is provided in case any of your initial values are
#   dependent on any of these arguments (e.g. the number of chains).
#
# Note: you _don't_ need to pass anything to either of these arguments, you only
#   use it within the function. `bbr` will pass in the correct objects when it calls
#   `make_init()` under the hood.
#
make_init <- function(.data, .args) {
  
  temp = data.frame(metabolite_number=.data$metabolite_number, Y = .data$Y)
  r_metabolite_init <- temp%>%group_by(metabolite_number)%>%
    summarize(mY = median(Y)) %>%
    select(mY)%>% 
    pull()
  
  function(){
    list(mean_metabolite = rnorm(1, 10, 0.25),
         mean_dilution  = rnorm(1, 1, 0.25),
         mean_injection  = rnorm(1, 1, 0.25),
         mean_dilution_injection_interaction  = rnorm(1, 0, 0.25),
         mean_sigma_dilution  = rnorm(1, 0, 0.25),
         mean_sigma_injection  = rnorm(1, 0, 0.25),
         mean_sigma_dilution_injection_interaction  = rnorm(1, 0, 0.25),
         mean_log_sigma = rnorm(1, -1, 0.25),
         dmean_metabolite = rnorm(1, 0, 0.25),
         dmean_dilution  = rnorm(1, 0, 0.25),
         dmean_injection  = rnorm(1, 0, 0.25),
         dmean_dilution_injection_interaction  = rnorm(1, 0, 0.25),
         dmean_sigma_dilution  = rnorm(1, 0, 0.25),
         dmean_sigma_injection  = rnorm(1, 0, 0.25),
         dmean_sigma_dilution_injection_interaction  = rnorm(1, 0, 0.25),
         dmean_log_sigma = rnorm(1, 0, 0.25),
         L_ms = t(chol(matrix(c(1,0.5,0.5,1), nrow=2))),
         L_di = t(chol(matrix(c(1,-0.5,-0.5,1), nrow=2))),
         sd_metabolite = exp(rnorm(1, log(0.5), 0.25)),
         sd_dilution = exp(rnorm(1, log(0.5), 0.25)),
         sd_injection = exp(rnorm(1, log(0.5), 0.25)),
         sd_dilution_injection_interaction = exp(rnorm(1, log(0.5), 0.25)),
         sd_sigma_dilution = exp(rnorm(1, log(0.5), 0.25)),
         sd_sigma_injection = exp(rnorm(1, log(0.5), 0.25)),
         sd_sigma_dilution_injection_interaction = exp(rnorm(1, log(0.5), 0.25)),
         sd_log_sigma = exp(rnorm(1, log(0.5), 0.25)),
         r_log_sigma = rep(-1,.data$number_metabolites),
         r_nu = rgamma(.data$number_metabolites,2,0.1) + 3.1,
         r_chrom = matrix(rep(0,.data$number_chroms*2), ncol =2),
         sigma_shared = exp(rnorm(2, log(0.5), 0.25)),
         r_ms = matrix(c(r_metabolite_init, rep(-1,.data$number_metabolites)),ncol = 2),
         r_di = matrix(c(rep(1,.data$number_metabolites), rep(1,.data$number_metabolites)),ncol = 2),
         r_injection = rep(1,.data$number_metabolites),
         r_dilution_injection_interaction = rep(0,.data$number_metabolites),
         r_sigma_dilution = rep(0,.data$number_metabolites),
         r_sigma_injection = rep(0,.data$number_metabolites),
         r_sigma_dilution_injection_interaction = rep(0,.data$number_metabolites))
  }
}
