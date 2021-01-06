rm(list = ls())

#library(flexsurvPlus)
#local version for testing
devtools::install()

library(flexsurvPlus)
library(tibble)
library(dplyr)
library(survival)
library(survminer)
library(tidyr)
library(boot)

adtte <- sim_adtte(seed = 2020, rho = 0.6) # simulate data with a medium correlation between PFS & OS

# subset OS data and rename
OS_data <- adtte %>%
  filter(PARAMCD=="OS") %>%
  transmute(USUBJID,
            ARMCD,
            OS_days = AVAL,
            OS_event = 1- CNSR
  )

# subset PFS data and rename
PFS_data <- adtte %>%
  filter(PARAMCD=="PFS") %>%
  transmute(USUBJID,
            ARMCD,
            PFS_days = AVAL,
            PFS_event = 1- CNSR
  )

analysis_data <- left_join(OS_data, PFS_data, by = c("USUBJID", "ARMCD"))


# fit basic models

n.sim <- 100


PSM_bootstraps_PFS <- boot(
  statistic = bootPSM, # bootstrap function
  R=n.sim, # number of bootstrap samples
  data=analysis_data,
  time_var="PFS_days",
  event_var="PFS_event",
  model.type=c("Common shape", "Separate"),
  distr = c('weibull',
            'gamma'),
  strata_var = "ARMCD",
  int_name="A",
  ref_name = "B"
)

PSM_boot_data_PFS <- bootPSMtidy(PSM_bootstraps_PFS)

PSM_bootstraps_PFS$seed


PSM_bootstraps_PFS$t0
PSM_bootstraps_PFS$t

PSM_bootstraps_PFS$call$int_name

PSM_boot_data_PFS <- as.data.frame(PSM_bootstraps_PFS$t)
colnames(PSM_boot_data_PFS) <- names(PSM_bootstraps_PFS$t0)
head(PSM_boot_data_PFS)

PSM_boot_data_PFS

# add some meta data
PSM_boot_data_PFS


# Read in ADaM data and rename variables of interest


adtte <- sim_adtte(seed = 2020, rho = 0.6)

# subset OS data and rename
OS_data <- adtte %>%
  filter(PARAMCD=="OS") %>%
  transmute(USUBJID,
            ARMCD,
            OS_days = AVAL,
            OS_event = 1- CNSR
  )

data=OS_data
time_var="OS_days"
event_var="OS_event"
model.type= c("Common shape", "Independent shape", "Separate")
distr = c('exp',
                 'weibull',
                 'gompertz',
                 'lnorm',
                 'llogis',
                 'gengamma',
                 'gamma')
strata_var = "ARMCD"
int_name="A"
ref_name = "B"

i = 1:500
data_boot <- data[i,]



bootPSM <- function(data, i, ...){
  boot_data = data[i,]
  output <- runPSM(data = boot_data, ...)
  return(output$parameters_vector)
}







a <- runPSM(data=OS_data,
                     time_var="OS_days",
                     event_var="OS_event",
                     model.type= c("Common shape", "Independent shape"),
                     distr = c('exp'),
                     strata_var = "ARMCD",
                     int_name="A",
                     ref_name = "B")

data_standard=Format_data(data, time_var, event_var, strata_var, int_name, ref_name)
dist = "exp"
model.formula=Surv(Time, Event==1) ~ ARM
flexsurv::flexsurvreg(formula=model.formula, data=data_standard, dist="exp")

data_standard$ARM
levels(data_standard$ARM)

models <- fit_models(model.formula=model.formula, distr = "exp", data=data_standard)

flexsurvPlus:::fit_models(data = dat, distr = "exp", model.formula = model.formula)
  
b <- runPSM(data_boot, time_var, event_var, 
                 model.type = "Separate", distr,
                 strata_var, int_name, ref_name)
flexsurv.dists
a$parameters
psm_OS_all

i=1:500

adata <- tibble(x=1:10,y=x^2)

fx <- function(data, i, opt1 = "tibble"){
  
  data_boot = data[i,]
  
  stats <- summarise(data_boot, x = sum(x), y = sum(y))
  
  if (opt1 == "tibble"){
    rc <- stats 
  } else {
    rc <- c(stats$x, stats$y)
    names(rc) <- names(stats)
  } 
  
  return(rc)  
}

fx(adata)
fx(adata, opt1 = "v")
seed = 1234
a <- boot(data = adata, statistic = fx, R = 5)

set.seed(1234)
b <- boot(data = adata, statistic = fx, R = 5, opt1 = "v")
set.seed(1234)
c <- boot(data = adata, statistic = fx, R = 5, opt1 = "v")

boot.array(b)
boot.array(c)


b$t0
b$t
names(b)
b$data
b$stype
selected_models <- c("weibull.comshp","weibull.indshp","gamma.comshp")
boot.array(b, indices = TRUE)


# Create data set for plots
curvefits_PFS <- get_curvefits(models = psm_PFS_all$models[selected_models],
                               time = seq(from=0, to = 500, by = 10))

# format data ready for the plot
curvefits_df_PFS <- bind_rows(curvefits_PFS$curvefits, .id = "Dist") %>%
  pivot_wider(names_from  = c(ARM), values_from  = c(est))

analysis_data_TrtA <- filter(analysis_data, ARMCD=="A")
analysis_data_TrtB <- filter(analysis_data, ARMCD=="B")

# Get KM estimates
km.est.PFS.TrtA <- survfit(Surv(PFS_days, PFS_event) ~ 1 , data = analysis_data_TrtA, conf.type = 'plain')
km.est.PFS.TrtB <- survfit(Surv(PFS_days, PFS_event) ~ 1 , data = analysis_data_TrtB, conf.type = 'plain')


# Plot survival curves for arm A
plot.PFS.A <- ggsurvplot(km.est.PFS.TrtA,
                         combine = TRUE,
                         censor = FALSE,
                         risk.table = TRUE,
                         conf.int = FALSE,
                         break.x.by = 100,
                         xlim = c(0, 500),
                         xlab = "Time (days)",
                         size = 0.72,
                         linetype = c(6,6,6,6, 1, 6, 6,6,6),
                         title  =  'Kaplan-Meier of OS with standard distribution overlays',
                         legend.title = '',
                         legend.labs = c("KM"),
                         surv.median.line = 'hv',
                         palette = c(rgb(0,0,0,max=30)),
                         risk.table.y.text.col = F
)


plot.PFS.A$plot <- plot.PFS.A$plot +
  geom_line(aes(x = time, y = Int,
                colour = Dist,
                linetype = Dist),
            size = 0.72,
            data = curvefits_df_PFS) +
  scale_color_manual(values = c(blue, pink, orange, red))

flexsurvPlus:::validate_standard_data(data=analysis_data,
                                      time_var="OS_days",
                                      event_var="OS_event",
                                      strata_var = "ARMCD",
                                      int_name="A",
                                      ref_name = "B")


km.est.OS <- survfit(Surv(OS_days, OS_event) ~ ARMCD , data = analysis_data, conf.type = 'plain')
km.est.PFS <- survfit(Surv(PFS_days, PFS_event) ~ ARMCD , data = analysis_data, conf.type = 'plain')

#KM_list <- list(OS = km.est.OS, PFS = km.est.PFS)

KM_plot_OS <- ggsurvplot(km.est.OS, risk.table = TRUE, data = analysis_data,
                         break.time.by = 100,
                         conf.int = FALSE,
                         censor=FALSE,
                         legend.title = '',
                         xlab = paste0('Overall survival (Days)'),
                         size = 0.72,
                         xlim = c(0, 700))
KM_plot_OS

KM_plot_PFS <- ggsurvplot(km.est.PFS, risk.table = TRUE, data = analysis_data,
                          break.time.by = 100,
                          conf.int = FALSE,
                          censor=FALSE,
                          legend.title = '',
                          xlab = paste0('Progression-free survival (Days)'),
                          size = 0.72,
                          xlim = c(0, 700))
KM_plot_PFS

# Fitting the models

psm_OS <- runPSM(data=analysis_data,
                          time_var="OS_days",
                          event_var="OS_event",
                          model.type=c("Separate","Common shape"),
                          distr = c('exp',
                                    'weibull',
                                    'gompertz',
                                    'lnorm',
                                    'llogis',
                                    'gengamma',
                                    'gamma'),
                          strata_var = "ARM",
                          int_name="A",
                          ref_name = "B")


get_curvefits(psm_OS$models, time = c(0,1,2))

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













###############################################################

devtools::install()

rm(list = ls())

# simulate data with a medium correlation between PFS & OS on the patient level
adtte <- sim_adtte(seed = 2020, rho = 0.6) 

# subset OS data and rename
OS_data <- adtte %>%
  filter(PARAMCD=="OS") %>%
  transmute(USUBJID,
            ARMCD,
            OS_days = AVAL,
            OS_event = 1- CNSR
  )

# subset PFS data and rename
PFS_data <- adtte %>%
  filter(PARAMCD=="PFS") %>%
  transmute(USUBJID,
            ARMCD,
            PFS_days = AVAL,
            PFS_event = 1- CNSR
  )

analysis_data <- left_join(OS_data, PFS_data, by = c("USUBJID", "ARMCD"))




# Bootstrapping the endpoints

## Ignoring correlation between endpoints


set.seed(2358)
n.sim <- 100

PSM_bootstraps_PFS <- boot(
  statistic = bootPSM, # bootstrap function
  R=n.sim, # number of bootstrap samples
  data=analysis_data,
  time_var="PFS_days",
  event_var="PFS_event",
  model.type=c("Common shape"),
  distr = c('weibull'),
  strata_var = "ARMCD",
  int_name="B",
  ref_name = "A"
)

PSM_bootstraps_OS <- boot(
  statistic = bootPSM, # bootstrap function
  R=n.sim, # number of bootstrap samples
  data=analysis_data,
  time_var="OS_days",
  event_var="OS_event",
  model.type=c("Separate"),
  distr = c('gamma'),
  strata_var = "ARMCD",
  int_name="B",
  ref_name = "A"
)

# we can see these are using different bootstrap samples by looking at the bootstrap indexes

index_PFS <- boot.array(PSM_bootstraps_PFS, indices = TRUE)
index_OS  <- boot.array(PSM_bootstraps_OS, indices = TRUE)

all(index_OS == index_PFS)




n.sim <- 100

set.seed(2020)

PSM_bootstraps_PFScor <- boot(
  statistic = bootPSM, # bootstrap function
  R=n.sim, # number of bootstrap samples
  data=analysis_data,
  time_var="PFS_days",
  event_var="PFS_event",
  model.type=c("Common shape"),
  distr = c('weibull'),
  strata_var = "ARMCD",
  int_name="B",
  ref_name = "A"
)

set.seed(2020)

PSM_bootstraps_OScor <- boot(
  statistic = bootPSM, # bootstrap function
  R=n.sim, # number of bootstrap samples
  data=analysis_data,
  time_var="OS_days",
  event_var="OS_event",
  model.type=c("Separate"),
  distr = c('gamma'),
  strata_var = "ARMCD",
  int_name= "B",
  ref_name = "A"
)


# we can check these are using the same bootstrap samples by comparing the indexes selected

index_PFScor <- boot.array(PSM_bootstraps_PFScor, indices = TRUE)
index_OScor  <- boot.array(PSM_bootstraps_OScor, indices = TRUE)

all(index_OScor == index_PFScor)



# To illustrate the impact we will look at plotting the difference between the extrapolated mean PFS and mean OS. To do this we need
# to calculate the mean PFS and mean OS for the reference and intervention arm for each sample

os_means <- bind_rows(
  mutate(bootPSMtidy(PSM_bootstraps_OS), Bootstrap = "Uncorrelated"), 
  mutate(bootPSMtidy(PSM_bootstraps_OScor), Bootstrap = "Correlated")
  ) %>%
  transmute(SampleID, Bootstrap, 
    os_mean_time_ref = flexsurv::mean_gamma(shape = gamma.shape.ref, rate = gamma.rate.ref),
    os_mean_time_int = flexsurv::mean_gamma(shape = gamma.shape.int, rate = gamma.rate.int),
    os_mean_delta = os_mean_time_int - os_mean_time_ref
  ) 

pfs_means <- bind_rows(
  mutate(bootPSMtidy(PSM_bootstraps_PFS), Bootstrap = "Uncorrelated"), 
  mutate(bootPSMtidy(PSM_bootstraps_PFScor), Bootstrap = "Correlated")
) %>%
  transmute(SampleID, Bootstrap, 
            pfs_mean_time_ref = flexsurv::mean_weibull(shape = weibull.shape.ref, scale = weibull.scale.ref),
            pfs_mean_time_int = flexsurv::mean_weibull(shape = weibull.shape.int, scale = weibull.scale.int),
            pfs_mean_delta = pfs_mean_time_int - pfs_mean_time_ref
  ) 


mean_durations <- pfs_means %>%
  left_join(os_means, by = c("SampleID", "Bootstrap"))

# we now transform the data so we have PFS and OS data in one row

mean_durations %>% 
  ggplot(aes(x = pfs_mean_delta, y = os_mean_delta, color = Bootstrap)) +
  theme_bw() +
  geom_point() +
  stat_ellipse() +
  coord_cartesian(xlim = c(100, 300), ylim = c(0, 2500)) +
  xlab("Difference in mean (weibull) PFS (B vs A)") +
  ylab("Difference in mean (gamma) OS (B vs A)") +
  ggtitle("Each point is a bootstrap sample")

  


expfir <- flexsurvreg(Surv(OS_days, event = OS_event == 1) ~ARMCD, dist ="exp", data = analysis_data )
expfir$coefficients
mean_exp(expfi$parameters$exp.rate.ref)
mean_exp(expfi$parameters$exp.rate.int)
mean_exp(-7.09)
mean_exp(-7.09-0.56)

survfit(Surv(OS_days, event = OS_event == 1) ~ARMCD, data = analysis_data) %>%
  plot()

get_curvefits(expfi$models, time = 200)$curvefits

flexsurv:::summary.flexsurvreg() 
summary(expfir, type = "mean")
exp(coef(expfir))
mean_exp
plot(expfir)

require(flexsurv)
mean_weibullPH(0.086, 7.5)
mean_weibull(0.086, 7.5)

coefs.flexsurvreg

flexsurv.dists()


flexsurv::mean_weibull(shape, scale = 1)

flexsurv::meam



