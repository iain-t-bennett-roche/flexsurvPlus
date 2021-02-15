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
  
  if('exp' %in% distr){
    models.exp <- fit_models(model.formula=model.formula.int, distr = "exp", data=data_standard)
    params.exp <- get_params(models=models.exp)
    param_out.exp <- t(unlist(params.exp$coef)) %>% as.data.frame() %>%
      dplyr::mutate(
        exp.rate.int = exp.rate + exp.ARMInt,
        exp.rate.ref = exp.rate,
        exp.rate.TE = exp.ARMInt) %>%
      dplyr::select(-exp.rate, -exp.ARMInt)
    
    # append this model to output 
    models$indshp.exp <- models.exp$exp
    params_out <- dplyr::bind_cols(params_out, param_out.exp)
    params.exp$name <- "indshp.exp"
    params <- dplyr::bind_rows(params, params.exp)
    
  }

  if('weibull' %in% distr){
    models.weib <- fit_models(model.formula=model.formula.shape, distr = "weibull", data=data_standard)
    params.weib <- get_params(models=models.weib)
    param_out.weib <- t(unlist(params.weib$coef)) %>% as.data.frame() %>%
      dplyr::mutate(
        weibull.scale.int = weibull.scale + weibull.ARMInt,
        weibull.scale.ref = weibull.scale,
        weibull.shape.int = weibull.shape + `weibull.shape(ARMInt)`,
        weibull.shape.ref = weibull.shape,
        weibull.scale.TE = weibull.ARMInt,
        weibull.shape.TE = `weibull.shape(ARMInt)`) %>%
      select(-weibull.scale, -weibull.shape, -weibull.ARMInt, -`weibull.shape(ARMInt)`)
    
    # append this model to output 
    models$indshp.weibull <- models.weib$weibull
    params_out <- dplyr::bind_cols(params_out, param_out.weib)
    params.weib$name <- "indshp.weibull"
    params <- dplyr::bind_rows(params, params.weib)
  }

  if('gompertz' %in% distr){
    models.gomp <- fit_models(model.formula=model.formula.shape, distr = "gompertz", data=data_standard)
    params.gomp <- get_params(models=models.gomp)
    param_out.gomp <- t(unlist(params.gomp$coef)) %>% as.data.frame() %>%
      dplyr::mutate(
        gompertz.rate.int = gompertz.rate + gompertz.ARMInt,
        gompertz.rate.ref = gompertz.rate,
        gompertz.shape.int = gompertz.shape + `gompertz.shape(ARMInt)`,
        gompertz.shape.ref = gompertz.shape,
        gompertz.rate.TE = gompertz.ARMInt,
        gompertz.shape.TE = `gompertz.shape(ARMInt)`) %>%
      select(-gompertz.rate, -gompertz.shape, -gompertz.ARMInt,-`gompertz.shape(ARMInt)`)
    
    # append this model to output 
    models$indshp.gompertz <- models.gomp$gompertz
    params_out <- dplyr::bind_cols(params_out, param_out.gomp)
    params.gomp$name <- "indshp.gompertz"
    params <- dplyr::bind_rows(params, params.gomp)
  }

  if('llogis' %in% distr){
    models.llogis <- fit_models(model.formula=model.formula.shape, distr = "llogis", data=data_standard)
    params.llogis <- get_params(models=models.llogis)
    param_out.llogis <- t(unlist(params.llogis$coef)) %>% as.data.frame() %>%
      dplyr::mutate(
        llogis.scale.int = llogis.scale + llogis.ARMInt,
        llogis.scale.ref = llogis.scale,
        llogis.shape.int = llogis.shape + `llogis.shape(ARMInt)`,
        llogis.shape.ref = llogis.shape,
        llogis.scale.TE = llogis.ARMInt,
        llogis.shape.TE = `llogis.shape(ARMInt)`) %>%
      select(-llogis.scale, -llogis.shape, -llogis.ARMInt, -`llogis.shape(ARMInt)`)

    # append this model to output 
    models$indshp.llogis <- models.llogis$llogis
    params_out <- dplyr::bind_cols(params_out, param_out.llogis)
    params.llogis$name <- "indshp.llogis"
    params <- dplyr::bind_rows(params, params.llogis)
  }

  if('gamma' %in% distr){
    models.gamma <- fit_models(model.formula=model.formula.shape, distr = "gamma", data=data_standard)
    params.gamma <- get_params(models=models.gamma)
    param_out.gamma <- t(unlist(params.gamma$coef)) %>% as.data.frame() %>%
      dplyr::mutate(
        gamma.rate.int = gamma.rate + gamma.ARMInt,
        gamma.rate.ref = gamma.rate,
        gamma.shape.int = gamma.shape + `gamma.shape(ARMInt)`,
        gamma.shape.ref = gamma.shape,
        gamma.rate.TE = gamma.ARMInt,
        gamma.shape.TE = `gamma.shape(ARMInt)`) %>%
      select(-gamma.rate, -gamma.shape, -gamma.ARMInt, -`gamma.shape(ARMInt)`)
    
    # append this model to output 
    models$indshp.gamma <- models.gamma$gamma
    params_out <- dplyr::bind_cols(params_out, param_out.gamma)
    params.gamma$name <- "indshp.gamma"
    params <- dplyr::bind_rows(params, params.gamma)
  }

  if('lnorm' %in% distr){
    models.lnorm <- fit_models(model.formula=model.formula.sdlog, distr = "lnorm", data=data_standard)
    params.lnorm <- get_params(models=models.lnorm)
    param_out.lnorm <- t(unlist(params.lnorm$coef)) %>% as.data.frame() %>%
      dplyr::mutate(
        lnorm.meanlog.int = lnorm.meanlog + lnorm.ARMInt,
        lnorm.meanlog.ref = lnorm.meanlog,
        lnorm.sdlog.int = lnorm.sdlog + `lnorm.sdlog(ARMInt)`,
        lnorm.sdlog.ref = lnorm.sdlog,
        lnorm.meanlog.TE = lnorm.ARMInt,
        lnorm.sdlog.TE = `lnorm.sdlog(ARMInt)`) %>%
      select(-lnorm.meanlog, -lnorm.sdlog, -lnorm.ARMInt, -`lnorm.sdlog(ARMInt)`)

    # append this model to output 
    models$indshp.lnorm <- models.lnorm$lnorm
    params_out <- dplyr::bind_cols(params_out, param_out.lnorm)
    params.lnorm$name <- "indshp.lnorm"
    params <- dplyr::bind_rows(params, params.lnorm)
  }

  if('gengamma' %in% distr){
    models.gengamma <- fit_models(model.formula=model.formula.sigma_Q, distr = "gengamma", data=data_standard)
    params.gengamma <- get_params(models=models.gengamma)
    param_out.gengamma <- t(unlist(params.gengamma$coef)) %>% as.data.frame() %>%
      dplyr::mutate(
        gengamma.mu.int = gengamma.mu + gengamma.ARMInt,
        gengamma.mu.ref = gengamma.mu,
        gengamma.sigma.int = gengamma.sigma + `gengamma.sigma(ARMInt)`,
        gengamma.sigma.ref = gengamma.sigma,
        gengamma.Q.int = gengamma.Q + `gengamma.Q(ARMInt)`,
        gengamma.Q.ref = gengamma.Q,
        gengamma.mu.TE = gengamma.ARMInt,
        gengamma.sigma.TE = `gengamma.sigma(ARMInt)`,
        gengamma.Q.TE = `gengamma.Q(ARMInt)`) %>%
      select(-gengamma.mu, -gengamma.sigma, -gengamma.Q, -gengamma.ARMInt, -`gengamma.sigma(ARMInt)`, -`gengamma.Q(ARMInt)`)

    # append this model to output 
    models$indshp.gengamma <- models.gengamma$gengamma
    params_out <- dplyr::bind_cols(params_out, param_out.gengamma)
    params.gengamma$name <- "indshp.gengamma"
    params <- dplyr::bind_rows(params, params.gengamma)
  }
  
  if('genf' %in% distr){
    models.genf <- fit_models(model.formula=model.formula.sigma_Q_P, distr = "genf", data=data_standard)
    params.genf <- get_params(models=models.genf)
    param_out.genf <- t(unlist(params.genf$coef)) %>% as.data.frame() %>%
      dplyr::mutate(
        genf.mu.int = genf.mu + genf.ARMInt,
        genf.mu.ref = genf.mu,
        genf.sigma.int = genf.sigma + `genf.sigma(ARMInt)`,
        genf.sigma.ref = genf.sigma,
        genf.Q.int = genf.Q + `genf.Q(ARMInt)`,
        genf.Q.ref = genf.Q,
        genf.P.int = genf.P + `genf.P(ARMInt)`,
        genf.P.ref = genf.P,
        genf.mu.TE = genf.ARMInt,
        genf.sigma.TE = `genf.sigma(ARMInt)`,
        genf.Q.TE = `genf.Q(ARMInt)`,
        genf.P.TE = `genf.P(ARMInt)`) %>%
      select(-genf.mu, -genf.sigma, -genf.Q, -genf.P, -genf.ARMInt, -`genf.sigma(ARMInt)`, -`genf.Q(ARMInt)`, -`genf.P(ARMInt)`)
    
    # append this model to output 
    models$indshp.genf <- models.genf$genf
    params_out <- dplyr::bind_cols(params_out, param_out.genf)
    params.genf$name <- "indshp.genf"
    params <- dplyr::bind_rows(params, params.genf)
  }

  #######################################################
  # prepare parameter outputs
  # this function exponentiates values that coef returns on the log scale
  # e.g. weibull shape and scale
  # this further simplifies other function use
  param_final <- post_process_param_out(params_out)
  
  # as a data frame with metadata 
  param_df <- param_final %>%
    dplyr::mutate(Model="Independent shape", Intervention_name=int_name, Reference_name=ref_name)
  
  # as a vector version with just numerics - needed for bootstrapping
  paramV <- as.numeric(param_final)
  names(paramV) <- paste0("indshp.", colnames(param_final))
  
  #######################################################
  

  # prepare parameter outputs
  
  output <- list(
    models = models,
    model_summary = params,
    parameters = param_df,
    parameters_vector = paramV
  )
  return(output)
}
