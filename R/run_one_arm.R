#' Run a complete parametric survival analysis for one-arm
#'
#' Fits a single \code{\link{flexsurv}} models using \code{\link{flexsurvreg}}.
#'
#' @param data A data frame containing individual patient data for the relevant
#'   time to event outcomes. This is passed to the \code{data} argument of the
#'   \code{\link{fit_models}} function
#' @param time_var Name of time variable in 'data'. Variable must be numerical and >0.
#' @param  event_var Name of event variable in 'data'. Variable must be
#'   numerical and contain 1's to indicate an event and 0 to indicate a censor.
#' @param int_name Character to indicate the name of the treatment of interest,
#'   must be a level of the "strata_var" column in "data", used for labeling
#'   the parameters.
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
#'   The model fit is in the form Surv(Time, Event==1) ~ 1 and is fit to the entire data (no strata). The parameters for each
#'   treatment, are derived directly from the model (no additional manipulation
#'   is required).
#' @return A list containing 'models' (output from \code{\link{fit_models}}), 'model_summary' (output from\code{\link{get_params}}) and
#'   'parameters', a data frame containing the coefficients of each flexsurv model.
#' \itemize{
#'   \item 'models' is a list of flexsurv objects for each distribution specified
#'   \item 'model_summary' is a tibble object containing the fitted model objects, the parameter
#'   estimates (\code{\link{coef}}),  \code{\link{AIC}} and \code{\link{BIC}}
#'   from flexsurv objects.
#'   \item 'parameters' is a data frame with with one row which contains the coefficients for all of the flexsurv models specified.
#'    The column names are in the format 'distribution.parameter.TreatmentName', for example, weibull.shape.Intervention refers to the shape parameter
#'     of the weibull distribution for the intervention treatment.}
#'
#' @export
run_one_arm <- function(data,
                   time_var, event_var,
                   distr = c('exp',
                             'weibull',
                             'gompertz',
                             'lnorm',
                             'llogis',
                             'gengamma',
                             'gamma', 
                             'genf'),
                      int_name){


  #test that only valid distributions have been provided
  #This is also tested within fit_models. Consider eliminating here to avoid redundancy
  allowed_dist <- c('exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma', 'genf')
  assertthat::assert_that(
    all(distr %in% allowed_dist),
    msg = "Only the following distributions are supported: 'exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma', 'genf' "
  )

  # standardise variable names
  data_standard=Format_data_onearm(data, time_var, event_var, int_name)
  model.formula.one.arm=Surv(Time, Event==1) ~ 1

  #Fit the models for seven standard distributions
  models.int <- fit_models(model.formula=model.formula.one.arm, distr = distr, data=data_standard$dat.int)

  #get parameter estimates and model fit statistics
  params.int <- get_params(models=models.int)

  #Extract parameter estimates
  param_out.int <- t(unlist(params.int$coef))

  # Rename the parameter from the
  # exponential model to be consistent with output from other models
  suppressWarnings(colnames(param_out.int)[colnames(param_out.int) == 'exp'] <- "exp.rate")
  suppressWarnings(colnames(param_out.int) <- paste0(colnames(param_out.int),".int"))


  # rename for output
  names(models.int) <- paste0("onearm.int.", names(models.int))

  models <- c(models.int)
  
  params.int$name <-  paste0("onearm.int.", params.int$name)

  params <- dplyr::bind_rows(params.int)
  
  param_out <- cbind(param_out.int)  %>%
    as.data.frame() 
  
  #######################################################
  # prepare parameter outputs
  # this function exponentiates values that coef returns on the log scale
  # e.g. weibull shape and scale
  # this further simplifies other function use
  param_final <- post_process_param_out(param_out)
  
  # as a data frame with metadata 
  param_df <- param_final %>%
    dplyr::mutate(Model="One arm", Intervention_name=int_name)
  
  # as a vector version with just numerics - needed for bootstrapping
  paramV <- as.numeric(param_final)
  names(paramV) <- paste0("onearm.", colnames(param_final))
  
  #######################################################
  #collect and return output
  output <- list(
    models = models,
    model_summary = params,
    parameters = param_df,
    parameters_vector = paramV
  )

  return(output)
}
