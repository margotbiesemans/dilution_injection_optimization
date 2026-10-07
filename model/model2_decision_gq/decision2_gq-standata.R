# Create Stan data
#
# This function must return the list that will be passed to `data` argument
#   of `cmdstanr::sample()`
#
# The `.dir` argument represents the absolute path to the directory containing
#   this file. This is useful for building file paths to the input files you will
#   load. Note: you _don't_ need to pass anything to this argument, you only use
#   it within the function. `bbr` will pass in the correct path when it calls
#   `make_standata()` under the hood.
make_standata <- function(.dir) {
  
  xdata <- readr::read_csv(here::here("data", "derived", "data-posneg.csv")) 
  metabolite_number = match(xdata$metabolite_number, unique(xdata$metabolite_number))
  posneg = if_else(unique(metabolite_number) <= 184, 1, 0)
  data <- with(xdata,
               list(N = nrow(xdata),
                    Y = log_intensity,
                    dilution_factor = dilution_factor,
                    injection_volume = injection_volume,
                    posneg = posneg,
                    number_metabolites=length(unique(xdata$metabolite_number)),
                    metabolite_number =  match(xdata$metabolite_number, unique(xdata$metabolite_number)),
                    prior_only=0))
}