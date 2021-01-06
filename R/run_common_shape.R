#' Run a complete parametric survival analysis for a common shape model
#'
#' Fits \code{\link{flexsurv}} models using \code{\link{flexsurvreg}} containing a covariate for treatment on scale parameter only.
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
#'   the \code{\link{fit_models}} function
#' @details Possible distributions include:
#' \itemize{
#'   \item Exponential ('exp')
#'   \item Weibull ('weibull')
#'   \item Gompertz ('gompertz')
#'   \item Log-normal ('lnorm')
#'   \item Log-logistic ('llogis')
#'   \item Generalized gamma ('gengamma')
#'   \item Gamma ('gamma')
#'   }
#'   The model fit is in the form Surv(Time, Event==1) ~ ARM.
#'   The shape parameter is the same for each treatment, and derived directly from the model (no additional manipulation is required).
#'   The scale parameter is derived directly from the model for the reference category, however for the intervention arm, this is calculated as reference shape + treatment effect (shape).
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
#'      (applicable to the scale parameter only for the common shape model).
#'   }
#' @export
run_common_shape <- function(data,
                   time_var, event_var,
                   distr = c('exp',
                             'weibull',
                             'gompertz',
                             'lnorm',
                             'llogis',
                             'gengamma',
                             'gamma'),
                   strata_var,
                   int_name, ref_name){


  #test that only valid distributions have been provided
  #This is also tested within fit_models. Consider eliminating here to avoid redundancy
  allowed_dist <- c('exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma')
  assertthat::assert_that(
    all(distr %in% allowed_dist),
    msg = "Only the following distributions are supported: 'exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma' "
  )

  # standardise variable names
  data_standard=Format_data(data, time_var, event_var, strata_var, int_name, ref_name)
  model.formula=Surv(Time, Event==1) ~ ARM

  #Fit the models for seven standard distributions
  models <- fit_models(model.formula=model.formula, distr = distr, data=data_standard)
  
  #get parameter estimates and model fit statistics
  params <- get_params(models=models)

  #Extract parameter estimates
  param_out <- t(unlist(params$coef)) %>% as.data.frame()

  if('exp' %in% distr){
    param_out <- param_out %>%
      dplyr::mutate(
        exp.rate.int = exp.rate + exp.ARMInt,
        exp.rate.ref = exp.rate,
        exp.rate.TE = exp.ARMInt) %>%
      dplyr::select(-exp.rate, -exp.ARMInt)

      
  }

  if('weibull' %in% distr){
    param_out <- param_out %>%
      dplyr::mutate(
        weibull.scale.int = weibull.scale + weibull.ARMInt,
        weibull.scale.ref = weibull.scale,
        weibull.shape.int = weibull.shape,
        weibull.shape.ref = weibull.shape,
        weibull.scale.TE = weibull.ARMInt) %>%
      select(-weibull.scale, -weibull.shape, -weibull.ARMInt)
    
  }

  if('gompertz' %in% distr){
    param_out <- param_out %>%
      dplyr::mutate(
        gompertz.rate.int = gompertz.rate + gompertz.ARMInt,
        gompertz.rate.ref = gompertz.rate,
        gompertz.shape.int = gompertz.shape,
        gompertz.shape.ref = gompertz.shape,
        gompertz.rate.TE = gompertz.ARMInt) %>%
      select(-gompertz.rate, -gompertz.shape, -gompertz.ARMInt)
  
    
  }

  if('llogis' %in% distr){
    param_out <- param_out %>%
      dplyr::mutate(
        llogis.scale.int = llogis.scale + llogis.ARMInt,
        llogis.scale.ref = llogis.scale,
        llogis.shape.int = llogis.shape,
        llogis.shape.ref = llogis.shape,
        llogis.scale.TE = llogis.ARMInt) %>%
      select(-llogis.scale, -llogis.shape, -llogis.ARMInt)

  }

  if('gamma' %in% distr){
    param_out <- param_out %>%
      dplyr::mutate(
        gamma.rate.int = gamma.rate + gamma.ARMInt,
        gamma.rate.ref = gamma.rate,
        gamma.shape.int = gamma.shape,
        gamma.shape.ref = gamma.shape,
        gamma.rate.TE = gamma.ARMInt) %>%
      select(-gamma.rate, -gamma.shape, -gamma.ARMInt)
    
  }

  if('lnorm' %in% distr){
    param_out <- param_out %>%
      dplyr::mutate(
        lnorm.meanlog.int = lnorm.meanlog + lnorm.ARMInt,
        lnorm.meanlog.ref = lnorm.meanlog,
        lnorm.sdlog.int = lnorm.sdlog,
        lnorm.sdlog.ref = lnorm.sdlog,
        lnorm.meanlog.TE = lnorm.ARMInt) %>%
      select(-lnorm.meanlog, -lnorm.sdlog, -lnorm.ARMInt)

  }

  if('gengamma' %in% distr){
    param_out <- param_out %>%
      dplyr::mutate(
        gengamma.mu.int = gengamma.mu + gengamma.ARMInt,
        gengamma.mu.ref = gengamma.mu,
        gengamma.sigma.int = gengamma.sigma,
        gengamma.sigma.ref = gengamma.sigma,
        gengamma.Q.int = gengamma.Q,
        gengamma.Q.ref = gengamma.Q,
        gengamma.mu.TE = gengamma.ARMInt) %>%
      select(-gengamma.mu, -gengamma.sigma, -gengamma.Q, -gengamma.ARMInt)

  }

  # rename models so can bind with others without conflicts
  models.out <- models
  names(models.out) <- paste0("comshp.", names(models.out))
  
  params.out <- params
  params.out$name <- paste0("comshp.", params.out$name)
  
  #######################################################
  # prepare parameter outputs
  # this function exponentiates values that coef returns on the log scale
  # e.g. weibull shape and scale
  # this further simplifies other function use
  param_final <- post_process_param_out(param_out)
  
  # as a data frame with metadata 
  param_df <- param_final %>%
    dplyr::mutate(Model="Common shape", Intervention_name=int_name, Reference_name=ref_name)

  # as a vector version with just numerics - needed for bootstrapping
  paramV <- as.numeric(param_final)
  names(paramV) <- paste0("comshp.", colnames(param_final))
  
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
