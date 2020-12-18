

#library(flexsurvPlus)
#local version for testing
devtools::install()

library(flexsurvPlus)
library(tibble)
library(dplyr)
library(survival)
library(survminer)
library(tidyr)

# Read in ADaM data and rename variables of interest

adtte <- read.csv(system.file("extdata", "adtte.csv", package = "flexsurvPlus", mustWork = TRUE))

# subset OS data and rename
OS_data <- adtte %>%
  filter(PARAMCD=="OS") %>%
  mutate(ARM = as.factor(ARM)) %>%
  rename(OS_days = AVAL,
         OS_event = EVNT,
         OS_cens = CNSR
  ) %>%
  select(-PARAMCD, -PARAM)

# subset PFS data and rename
PFS_data <- adtte %>%
  filter(PARAMCD=="PFS") %>%
  rename(PFS_days = AVAL,
         PFS_event = EVNT,
         PFS_cens = CNSR
  ) %>%
  select(-PARAMCD, -PARAM)

analysis_data <- left_join(OS_data, PFS_data)
analysis_data$ARM <- relevel(as.factor(analysis_data$ARM), ref="B")

# Fitting the models

#These functions have been used to estimate 3 types of model:
  
# Treatment-effect models (TEM) with a common shape parameter for both treatments
# Independent shape models with different shape parameters for each treatment
# Separate models for each treatment - 2 models in total


## Fit treatment-effect models (TEM)

#Fit the models for seven standard distributions
models_TEM_OS <- fit_models(model.formula=Surv(OS_days, OS_event) ~ ARM,
                            data=analysis_data,
                            distr = c('exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma'))
models_TEM_OS

#get parameter estimates and model fit statistics
params_TEM_OS <- get_params(models=models_TEM_OS)
params_TEM_OS

# Meta data needed for labelling purposes
parout_TEM_OS <- data.frame(
  Study_name = "Study ABC",
  Treatment = "All",
  Datacut = "Final",
  Population = "ITT",
  Endpoint = "OS",
  Endpoint_def = "NA",
  Model = "Common shape - Treatment",
  SampleID = "Mean estimate")

psm_TEM_OS <- runPSM(model.formula = Surv(OS_days, OS_event) ~ ARM,
                     data = analysis_data,
                     metadata = parout_TEM_OS,
                     dist = c('exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma'),
                     model.type = 'Common shape',
                     toi_name="A", 
                     ref_name="B")
psm_TEM_OS

## Seperate models for each treatment

#split data set into individual arms
analysis_data_TrtA <- filter(analysis_data, ARM=="A")
analysis_data_TrtB <- filter(analysis_data, ARM=="B")

# Fit models to Trt A 

# Meta data needed for labelling purposes
parout_TrtA_OS <- data.frame(
  Study_name = "Study ABC",
  Treatment = "All",
  Datacut = "Final",
  Population = "ITT",
  Endpoint = "OS",
  Endpoint_def = "NA",
  Model = "Independent",
  SampleID = "Mean estimate")

psm_armA_OS <- runPSM(model.formula = Surv(OS_days, OS_event) ~ 1,
                      data = analysis_data_TrtA,
                      metadata = parout_TrtA_OS,
                      dist = c('exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma'),
                      model.type = 'Separate',
                      toi_name="A")

# Fit models to Trt B

# Meta data needed for labeling purposes
parout_TrtB_OS <- data.frame(
  Study_name = "Study ABC",
  Treatment = "All",
  Datacut = "Final",
  Population = "ITT",
  Endpoint = "OS",
  Endpoint_def = "NA",
  Model = "Independent",
  SampleID = "Mean estimate")

psm_armB_OS <- runPSM(model.formula = Surv(OS_days, OS_event) ~ 1,
                      data = analysis_data_TrtB,
                      metadata = parout_TrtB_OS,
                      dist = c('exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma'),
                      model.type = 'Separate',
                      toi_name="B")

# Combine the 2 independent models
psm_indep_OS <- full_join(psm_armA_OS$parameters, psm_armB_OS$parameters)

## Independant shape models 

# These models are similar to the separate model for each treatment as they allow all parameters to vary by treatment, however it is performed using 1 model for 2 treatments, rather than 2 models as above. They have been fit using the <tt>runPSM</tt> function.


# Independent shape models ------------------------------------------------
#This only applies for models with 2 or more parameters,i.e. not exponential

#Enter metadata
metadata_OS_ind <- data.frame(
  Study_name = "Study ABC",
  Treatment = "All",
  Datacut = "Final",
  Population = "ITT",
  Endpoint = "OS",
  Endpoint_def = "NA",
  Model = "Independent shape",
  SampleID = "Mean estimate"
)

# All distributions that have a parameter named shape - weibull, gompertz, llogis, gamma
psm_OS_ind_shape <- runPSM(model.formula=Surv(OS_days, OS_event) ~ ARM + shape(ARM),
                           data=analysis_data,
                           dist = c('weibull', 'gompertz', 'llogis',  'gamma'),
                           model.type = 'Independent shape',
                           metadata = metadata_OS_ind,
                           toi_name="A", 
                           ref_name="B")

#For the log normal distribution the shape parameter is named sdlog which
#requires a different formula
psm_OS_ind_lnorm <- runPSM(model.formula=Surv(OS_days, OS_event) ~ ARM + sdlog(ARM),
                           data=analysis_data,
                           dist = 'lnorm',
                           model.type = 'Independent shape',
                           metadata = metadata_OS_ind,
                           toi_name="A", 
                           ref_name="B")

#Generalised gamma model has two shape parameters in the flexsurv
#parameterisation, sigma and Q.
#Allowing both sigma and Q to be different between
#treatment arms is equivalent to fitting separate models to each arm
#Allowing only one of sigma or Q to be different between treatments is possible but hard to interpret
#Not implemented in the underlying functions
psm_OS_ind_gengamma <- runPSM(model.formula=Surv(OS_days, OS_event) ~ ARM + sigma(ARM) + Q(ARM),
                              data=analysis_data,
                              dist = 'gengamma',
                              model.type = 'Independent shape',
                              metadata = metadata_OS_ind,
                              toi_name="A", 
                              ref_name="B")


## Output parameters to excel

# Once all the models have been fit, they can be combined to output to excel such that there is 1 row per model.

parameters_all <- full_join(psm_OS_ind_shape$parameters, psm_OS_ind_lnorm$parameters) %>%
  full_join(psm_OS_ind_gengamma$parameters) %>%
  full_join(psm_TEM_OS$parameters) %>%
  full_join(psm_indep_OS) %>%
  #reorder columns
  select(colnames(metadata_OS_ind),
         starts_with("exp"), starts_with("weibull"), starts_with("gompertz"), starts_with("lnorm"),
         starts_with("llogis"), starts_with("gengamma"), starts_with("gamma")) %>%
  arrange(Model)


# check differences between models

library(reshape2)
mdf <- parameters_all %>%
  melt() %>%
  transmute(Model, variable, value) %>%
  arrange(variable, Model) %>%
  dcast(variable~Model)

mdf %>%
  transmute(variable, `Independent shape`, Independent, diff = abs(`Independent shape`-Independent)) %>%
  filter(diff > 0.001)
