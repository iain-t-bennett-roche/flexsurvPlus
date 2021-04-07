
rm(list = ls())
require(flexsurvPlus)

adtte <- sim_adtte(seed = 2020, # for reproducibility
                   rho = 0.6 # defines a correlation on underlying survival times
) 

# subset OS data for reference arm and rename
OS_data <- adtte %>%
  filter(PARAMCD=="OS", ARMCD == "A") %>%
  transmute(USUBJID,
            OS_days = AVAL,
            OS_event = 1- CNSR
  )


set.seed(2358)

# define a number of samples to generate. For illustration 100 are taken
# In practice a larger number is preferable.
n.sim <- 100

PSM_bootstraps_OS <- boot(
  statistic = bootPSM, # bootstrap function
  R=n.sim, # number of bootstrap samples
  data=OS_data, # the dataset
  time_var="OS_days",  # the time variable
  event_var="OS_event", # the event variable coded as 1 for event
  model.type="One arm", # again for speed only a single model is specified
  distr = "gengamma", 
  int_name = "A" # needed for meta data even with a single arm model
)

# This is error message, suspect needs some checking in the run_one_arm that coefs are only extracted for converged models

#Fitting one arm models
#Fitting model gengamma
#Fitting one arm models
#Fitting model gengamma
#Fitting one arm models
#Fitting model gengamma
#Fitting one arm models
#A warning occured in gengamma:
#  simpleWarning in flexsurv::flexsurvreg(formula = model.formula, data = data, dist = dist): Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite. 

#Fitting model gengamma
#Error in t.default(unlist(coef)) : argument is not a matrix

