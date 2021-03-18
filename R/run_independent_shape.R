#' Run a complete parametric survival analysis for an independent shape model
#'
#' Fits \code{\link{flexsurv}} models using \code{\link{flexsurvreg}}
#' containing a covariate for treatment on shape and scale parameter.
#'
#' @param data A data frame containing individual patient data for the relevant
#'   time to event outcomes. This is passed to the \code{data} argument of the
#'   \code{\link{fit_models}} function
#' @param time_var Name of time variable in 'data'. Variable must be numerical and >0.
#' @param  event_var Name of event variable in 'data'. Variable must be
#'   numerical and contain 1's to indicate an event and 0 to indicate a censor.
#' @param  strata_var Name of stratification variable in "data". This is usually
#'   the treatment variable and must be categorical.
#' @param int_name Character to indicate the name of the treatment of interest,
#'   must be a level of the "strata_var" column in "data", used for labeling
#'   the parameters.
#'  @param ref_name Character to indicate the name of the reference treatment,
#'    must be a level of the "strata_var" column in "data", used for labeling
#'    the parameters.
#' @param distr A vector string of distributions, see dist argument in
#'   \code{\link{flexsurvreg}}. This is passed to the \code{distr} argument of
#'   the \code{\link{fit_models}} function. Default is all available distributions (see below).
#' @details Possible distributions include:
#' \itemize{
#'   \item Exponential ('exp')
#'   \item Weibull ('weibull')
#'   \item Gompertz ('gompertz')
#'   \item Log-normal ('lnorm')
#'   \item Log-logistic ('llogis')
#'   \item Generalized gamma ('gengamma')
#'   \item Gamma ('gamma')
#'   \item Generalised F ('genf')
#'   }
#'
#'   The model fit is in the form Surv(Time, Event==1) ~ ARM + shape(ARM). The
#'   scale parameter is derived directly from the model for the reference
#'   category, however for the intervention arm, this is calculated as reference
#'   scale + treatment effect (scale). The shape parameter is derived directly
#'   from the model for the reference category, however for the intervention
#'   arm, this is calculated as reference shape + treatment effect (shape).
#' @return A list containing 'models' (output from \code{\link{fit_models}}), 'model_summary' (output from\code{\link{get_params}}) and
#'   'parameters', a data frame containing the coefficients of each flexsurv model.
#' \itemize{
#'   \item 'models' is a list of flexsurv objects for each distribution specified
#'   \item 'model_summary' is a tibble object containing the fitted model objects, the parameter
#'   estimates (\code{\link{coef}}),  \code{\link{AIC}} and \code{\link{BIC}}
#'   from flexsurv objects.
#'   \item 'parameters' is a data frame with with one row which contains the coefficients for all of the flexsurv models specified.
#'    The column names are in the format 'distribution.parameter.TreatmentName', for example, weibull.shape.Intervention refers to the shape parameter
#'     of the weibull distribution for the intervention treatment and 'gengamma.mu.Reference' refers to the mu parameter
#'     of the generalised gamma distribution for the reference treatment. Columns with 'TE' at the end are the treatment effect coefficients
#'      (applicable to the scale and shape parameters for independent shape models).}
#'
#' @export
run_independent_shape <- function(data,
                   time_var, event_var,
                   distr = c('exp',
                             'weibull',
                             'gompertz',
                             'lnorm',
                             'llogis',
                             'gengamma',
                             'gamma', 
                             'genf'),
                   strata_var,
                   int_name, ref_name){

  #test that only valid distributions have been provided
  #This is also tested within fit_models. Consider eliminating here to avoid redundancy
  allowed_dist <- c('exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma', 'genf')
  assertthat::assert_that(
    all(distr %in% allowed_dist),
    msg = "Only the following distributions are supported: 'exp', weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma', 'genf' "
  )

  # standardise variable names
  data_standard=Format_data(data, time_var, event_var, strata_var, int_name, ref_name)

  # Model formulas
  model.formula.int = Surv(Time, Event==1) ~ ARM
  model.formula.shape = Surv(Time, Event==1) ~ ARM + shape(ARM)
  model.formula.sdlog = Surv(Time, Event==1) ~ ARM + sdlog(ARM)
  model.formula.sigma_Q = Surv(Time, Event==1) ~ ARM + sigma(ARM) + Q(ARM)
  model.formula.sigma_Q_P = Surv(Time, Event==1) ~ ARM + sigma(ARM) + Q(ARM) + P(ARM)
  
  models <- list()
  params <- tibble()
  params_out <- tibble(.rows = 1)
  
  message("Fitting independant shape models")
  
  if('exp' %in% distr){
    models.exp <- fit_models(model.formula=model.formula.int, distr = "exp", data=data_standard)
  }
  
  if('weibull' %in% distr){
    models.weib <- fit_models(model.formula=model.formula.shape, distr = "weibull", data=data_standard)
  }
  
  if('gompertz' %in% distr){
    models.gomp <- fit_models(model.formula=model.formula.shape, distr = "gompertz", data=data_standard)
  }
  
  if('llogis' %in% distr){
    models.llogis <- fit_models(model.formula=model.formula.shape, distr = "llogis", data=data_standard)
  }
  
  if('gamma' %in% distr){
    models.gamma <- fit_models(model.formula=model.formula.shape, distr = "gamma", data=data_standard)
  }
  
  if('lnorm' %in% distr){
    models.lnorm <- fit_models(model.formula=model.formula.sdlog, distr = "lnorm", data=data_standard)
  }
  
  if('gengamma' %in% distr){
    models.gengamma <- fit_models(model.formula=model.formula.sigma_Q, distr = "gengamma", data=data_standard)
  }
    
  if('genf' %in% distr){
    models.genf <- fit_models(model.formula=model.formula.sigma_Q_P, distr = "genf", data=data_standard)
  }
  
  #Combine models
  models <- c(models.exp, models.weib, models.gomp, models.llogis, models.gamma, models.lnorm, models.gengamma, models.genf)  
  
  #get parameter estimates and model fit statistics
  params <- get_params(models=models)
  
  #Extract parameter estimates
  param_out <- t(unlist(params$coef)) %>% as.data.frame()
  
  if('exp' %in% params$name){
    param_out <- param_out %>%
      dplyr::mutate(
        exp.rate.int = exp.rate + exp.ARMInt,
        exp.rate.ref = exp.rate,
        exp.rate.TE = exp.ARMInt) %>%
      dplyr::select(-exp.rate, -exp.ARMInt)
    
    
  }
  
  if('weibull' %in% params$name){
    param_out <- param_out %>%
      dplyr::mutate(
        weibull.scale.int = weibull.scale + weibull.ARMInt,
        weibull.scale.ref = weibull.scale,
        weibull.shape.int = weibull.shape,
        weibull.shape.ref = weibull.shape,
        weibull.scale.TE = weibull.ARMInt,
        weibull.shape.TE = `weibull.shape(ARMInt)`) %>%
      select(-weibull.scale, -weibull.shape, -weibull.ARMInt, -`weibull.shape(ARMInt)`)
    
  }
  
  if('gompertz' %in% params$name){
    param_out <- param_out %>%
      dplyr::mutate(
        gompertz.rate.int = gompertz.rate + gompertz.ARMInt,
        gompertz.rate.ref = gompertz.rate,
        gompertz.shape.int = gompertz.shape,
        gompertz.shape.ref = gompertz.shape,
        gompertz.rate.TE = gompertz.ARMInt,
        gompertz.shape.TE = `gompertz.shape(ARMInt)`) %>%
      select(-gompertz.rate, -gompertz.shape, -gompertz.ARMInt, -`gompertz.shape(ARMInt)`)
    
    
  }
  
  if('llogis' %in% params$name){
    param_out <- param_out %>%
      dplyr::mutate(
        llogis.scale.int = llogis.scale + llogis.ARMInt,
        llogis.scale.ref = llogis.scale,
        llogis.shape.int = llogis.shape,
        llogis.shape.ref = llogis.shape,
        llogis.scale.TE = llogis.ARMInt,
        llogis.shape.TE = `llogis.shape(ARMInt)`) %>%
      select(-llogis.scale, -llogis.shape, -llogis.ARMInt, -`llogis.shape(ARMInt)`)
    
  }
  
  if('gamma' %in% params$name){
    param_out <- param_out %>%
      dplyr::mutate(
        gamma.rate.int = gamma.rate + gamma.ARMInt,
        gamma.rate.ref = gamma.rate,
        gamma.shape.int = gamma.shape,
        gamma.shape.ref = gamma.shape,
        gamma.rate.TE = gamma.ARMInt,
        gamma.shape.TE = `gamma.shape(ARMInt)`) %>%
      select(-gamma.rate, -gamma.shape, -gamma.ARMInt, -`gamma.shape(ARMInt)`)
    
  }
  
  if('lnorm' %in% params$name){
    param_out <- param_out %>%
      dplyr::mutate(
        lnorm.meanlog.int = lnorm.meanlog + lnorm.ARMInt,
        lnorm.meanlog.ref = lnorm.meanlog,
        lnorm.sdlog.int = lnorm.sdlog,
        lnorm.sdlog.ref = lnorm.sdlog,
        lnorm.meanlog.TE = lnorm.ARMInt,
        lnorm.sdlog.TE = `lnorm.sdlog(ARMInt)`) %>%
      select(-lnorm.meanlog, -lnorm.sdlog, -lnorm.ARMInt, -`lnorm.sdlog(ARMInt)`)
    
  }
  
  if('gengamma' %in% params$name){
    param_out <- param_out %>%
      dplyr::mutate(
        gengamma.mu.int = gengamma.mu + gengamma.ARMInt,
        gengamma.mu.ref = gengamma.mu,
        gengamma.sigma.int = gengamma.sigma,
        gengamma.sigma.ref = gengamma.sigma,
        gengamma.Q.int = gengamma.Q,
        gengamma.Q.ref = gengamma.Q,
        gengamma.mu.TE = gengamma.ARMInt,
        gengamma.sigma.TE = `gengamma.sigma(ARMInt)`,
        gengamma.Q.TE = `gengamma.Q(ARMInt)`) %>%
      select(-gengamma.mu, -gengamma.sigma, -gengamma.Q, -gengamma.ARMInt, -`gengamma.sigma(ARMInt)`, -`gengamma.Q(ARMInt)`)
    
  }
  
  if('genf' %in% params$name){
    param_out <- param_out %>%
      dplyr::mutate(
        genf.mu.int = genf.mu + genf.ARMInt,
        genf.mu.ref = genf.mu,
        genf.sigma.int = genf.sigma,
        genf.sigma.ref = genf.sigma,
        genf.Q.int = genf.Q,
        genf.Q.ref = genf.Q,
        genf.P.int = genf.P,
        genf.P.ref = genf.P,
        genf.mu.TE = genf.ARMInt,
        genf.sigma.TE = `genf.sigma(ARMInt)`,
        genf.Q.TE = `genf.Q(ARMInt)`,
        genf.P.TE = `genf.P(ARMInt)`) %>%
      select(-genf.mu, -genf.sigma, -genf.Q, -genf.P, -genf.ARMInt, -`genf.sigma(ARMInt)`, -`genf.Q(ARMInt)`, -`genf.P(ARMInt)`)
  }
  
  # rename models so can bind with others without conflicts
  models.out <- models
  names(models.out) <- paste0("indshp.", names(models.out))
  
  params.out <- params
  params.out$name <- paste0("indshp.", params.out$name)
  


  #######################################################
  # prepare parameter outputs
  # this function exponentiates values that coef returns on the log scale
  # e.g. weibull shape and scale
  # this further simplifies other function use
  param_final <- post_process_param_out(param_out)
  
  # as a data frame with metadata 
  param_df <- param_final %>%
    dplyr::mutate(Model="Independent shape", Intervention_name=int_name, Reference_name=ref_name)
  
  params_all <- param_df 
  
  if('exp' %in% distr){
    params_all$exp.rate.int <- ifelse("exp.rate.int" %in% names(params_all), params_all$exp.rate.int, NA) 
    params_all$exp.rate.ref <- ifelse("exp.rate.ref" %in% names(params_all), params_all$exp.rate.ref, NA) 
    params_all$exp.rate.TE <- ifelse("exp.rate.TE" %in% names(params_all), params_all$exp.rate.TE, NA) 
  }
  
  if('weibull' %in% distr){
    params_all$weibull.shape.int <- ifelse("weibull.shape.int" %in% names(params_all), params_all$weibull.shape.int, NA) 
    params_all$weibull.scale.int <- ifelse("weibull.scale.int" %in% names(params_all), params_all$weibull.scale.int, NA) 
    params_all$weibull.shape.ref <- ifelse("weibull.shape.ref" %in% names(params_all), params_all$weibull.shape.ref, NA) 
    params_all$weibull.scale.ref <- ifelse("weibull.scale.ref" %in% names(params_all), params_all$weibull.scale.ref, NA) 
    params_all$weibull.scale.TE <- ifelse("weibull.scale.TE" %in% names(params_all), params_all$weibull.scale.TE, NA) 
    params_all$weibull.shape.TE <- ifelse("weibull.shape.TE" %in% names(params_all), params_all$weibull.shape.TE, NA) 
  }
  
  if('gompertz' %in% distr){
    params_all$gompertz.shape.int <- ifelse("gompertz.shape.int" %in% names(params_all), params_all$gompertz.shape.int, NA) 
    params_all$gompertz.rate.int <- ifelse("gompertz.rate.int" %in% names(params_all), params_all$gompertz.rate.int, NA) 
    params_all$gompertz.shape.ref <- ifelse("gompertz.shape.ref" %in% names(params_all), params_all$gompertz.shape.ref, NA) 
    params_all$gompertz.rate.ref <- ifelse("gompertz.rate.ref" %in% names(params_all), params_all$gompertz.rate.ref, NA) 
    params_all$gompertz.rate.TE <- ifelse("gompertz.rate.TE" %in% names(params_all), params_all$gompertz.rate.TE, NA) 
    params_all$gompertz.shape.TE <- ifelse("gompertz.shape.TE" %in% names(params_all), params_all$gompertz.shape.TE, NA) 
    
  }
  
  if('llogis' %in% distr){
    params_all$llogis.shape.int <- ifelse("llogis.shape.int" %in% names(params_all), params_all$llogis.shape.int, NA) 
    params_all$llogis.scale.int <- ifelse("llogis.scale.int" %in% names(params_all), params_all$llogis.scale.int, NA) 
    params_all$llogis.shape.ref <- ifelse("llogis.shape.ref" %in% names(params_all), params_all$llogis.shape.ref, NA) 
    params_all$llogis.scale.ref <- ifelse("llogis.scale.ref" %in% names(params_all), params_all$llogis.scale.ref, NA) 
    params_all$llogis.scale.TE <- ifelse("llogis.scale.TE" %in% names(params_all), params_all$llogis.scale.TE, NA) 
    params_all$llogis.shape.TE <- ifelse("llogis.shape.TE" %in% names(params_all), params_all$llogis.shape.TE, NA) 
  }
  
  if('gamma' %in% distr){
    params_all$gamma.shape.int <- ifelse("gamma.shape.int" %in% names(params_all), params_all$gamma.shape.int, NA) 
    params_all$gamma.rate.int <- ifelse("gamma.rate.int" %in% names(params_all), params_all$gamma.rate.int, NA) 
    params_all$gamma.shape.ref <- ifelse("gamma.shape.ref" %in% names(params_all), params_all$gamma.shape.ref, NA) 
    params_all$gamma.rate.ref <- ifelse("gamma.rate.ref" %in% names(params_all), params_all$gamma.rate.ref, NA) 
    params_all$gamma.rate.TE <- ifelse("gamma.rate.TE" %in% names(params_all), params_all$gamma.rate.TE, NA) 
    params_all$gamma.shape.TE <- ifelse("gamma.shape.TE" %in% names(params_all), params_all$gamma.shape.TE, NA) 
  }
  
  if('lnorm' %in% distr){
    params_all$lnorm.meanlog.int <- ifelse("lnorm.meanlog.int" %in% names(params_all), params_all$lnorm.meanlog.int, NA) 
    params_all$lnorm.sdlog.int <- ifelse("lnorm.sdlog.int" %in% names(params_all), params_all$lnorm.sdlog.int, NA) 
    params_all$lnorm.meanlog.ref <- ifelse("lnorm.meanlog.ref" %in% names(params_all), params_all$lnorm.meanlog.ref, NA) 
    params_all$lnorm.sdlog.ref <- ifelse("lnorm.sdlog.ref" %in% names(params_all), params_all$lnorm.sdlog.ref, NA) 
    params_all$lnorm.meanlog.TE <- ifelse("lnorm.meanlog.TE" %in% names(params_all), params_all$lnorm.meanlog.TE, NA) 
    params_all$lnorm.sdlog.TE <- ifelse("lnorm.sdlog.TE" %in% names(params_all), params_all$lnorm.sdlog.TE, NA) 
    
  }
  
  if('gengamma' %in% distr){
    params_all$gengamma.mu.int <- ifelse("gengamma.mu.int" %in% names(params_all), params_all$gengamma.mu.int, NA) 
    params_all$gengamma.sigma.int <- ifelse("gengamma.sigma.int" %in% names(params_all), params_all$gengamma.sigma.int, NA) 
    params_all$gengamma.Q.int <- ifelse("gengamma.Q.int" %in% names(params_all), params_all$gengamma.Q.int, NA) 
    params_all$gengamma.mu.ref <- ifelse("gengamma.mu.ref" %in% names(params_all), params_all$gengamma.mu.ref, NA) 
    params_all$gengamma.sigma.ref <- ifelse("gengamma.sigma.ref" %in% names(params_all), params_all$gengamma.sigma.ref, NA) 
    params_all$gengamma.Q.ref <- ifelse("gengamma.Q.ref" %in% names(params_all), params_all$gengamma.Q.ref, NA) 
    params_all$gengamma.mu.TE <- ifelse("gengamma.mu.TE" %in% names(params_all), params_all$gengamma.mu.TE, NA) 
    params_all$gengamma.sigma.TE <- ifelse("gengamma.sigma.TE" %in% names(params_all), params_all$gengamma.sigma.TE, NA) 
    params_all$gengamma.Q.TE <- ifelse("gengamma.Q.TE" %in% names(params_all), params_all$gengamma.Q.TE, NA) 
  }
  
  if('genf' %in% distr){
    params_all$genf.mu.int <- ifelse("genf.mu.int" %in% names(params_all), params_all$genf.mu.int, as.numeric(NA)) 
    params_all$genf.sigma.int <- ifelse("genf.sigma.int" %in% names(params_all), params_all$genf.sigma.int, as.numeric(NA)) 
    params_all$genf.Q.int <- ifelse("genf.Q.int" %in% names(params_all), params_all$genf.Q.int, as.numeric(NA)) 
    params_all$genf.P.int <- ifelse("genf.P.int" %in% names(params_all), params_all$genf.P.int, as.numeric(NA)) 
    params_all$genf.mu.ref <- ifelse("genf.mu.ref" %in% names(params_all), params_all$genf.mu.ref, as.numeric(NA)) 
    params_all$genf.sigma.ref <- ifelse("genf.sigma.ref" %in% names(params_all), params_all$genf.sigma.ref, as.numeric(NA)) 
    params_all$genf.Q.ref <- ifelse("genf.Q.ref" %in% names(params_all), params_all$genf.Q.ref, as.numeric(NA)) 
    params_all$genf.P.ref <- ifelse("genf.P.ref" %in% names(params_all), params_all$genf.P.ref, as.numeric(NA)) 
    params_all$genf.mu.TE <- ifelse("genf.mu.TE" %in% names(params_all), params_all$genf.mu.TE, as.numeric(NA)) 
    params_all$genf.sigma.TE <- ifelse("genf.sigma.TE" %in% names(params_all), params_all$genf.sigma.TE, as.numeric(NA)) 
    params_all$genf.Q.TE <- ifelse("genf.Q.TE" %in% names(params_all), params_all$genf.Q.TE, as.numeric(NA)) 
    params_all$genf.P.TE <- ifelse("genf.P.TE" %in% names(params_all), params_all$genf.P.TE, as.numeric(NA)) 
    
  }
  
  params_all <- params_all %>%
    select(-Model, -Intervention_name, -Reference_name)
  
  
  # as a vector version with just numerics - needed for bootstrapping
  paramV <- as.numeric(params_all)
  names(paramV) <- paste0("indshp.", colnames(params_all))
  
  #######################################################
  
  #collect and return output
  output <- list(
    models = models.out,
    model_summary = params.out,
    parameters = param_df,
    parameters_vector = paramV
  )
  return(output)
  
}
