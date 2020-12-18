#' Run a complete parametric survival analysis for one model with multiple distributions
#'
#' @param data A data frame containing individual patient data for the relevant
#'   time to event outcomes. This is passed to the \code{data} argument of the
#'   \code{\link{fit_models}} function
#' @param time_var Name of time variable in 'data'. Variable must be numerical and >0.
#' @param  event_var Name of event variable in 'data'. Variable must be
#'   numerical and contain 1's to indicate an event and 0 to indicate a censor.
#' @param model.type Character vector indicating the type of model formula
#'   provided. Permitted values are
#' \itemize{
#'   \item 'Common shape' a model with a single covariate for the effect of
#'   treatment on the scale parameter of the model
#'   \item 'Independent shape' a model with a single covariate for treatment
#'   that affects both the scale and shape parameters of the model
#'   \item 'Separate' a model with no covariates typically fitted separately to
#'   data from each treatment group in a study
#'  }
#' @param  strata_var Name of stratification variable in "data". This is usually
#'   the treatment variable and must be categorical.
#' @param int_name Character to indicate the name of the treatment of interest,
#'   must be a level of the "strata_var" column in "data", used for labelling
#'   the parameters.
#' @param ref_name Character to indicate the name of the reference treatment,
#'    must be a level of the "strata_var" column in "data", used for labelling
#'    the parameters.
#' @param distr A vector string of distributions, see dist argument in
#'   \code{\link{flexsurvreg}}. This is passed to the \code{distr} argument of
#'   the \code{\link{fit_models}} function
#'
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
#'
#' @return A list containing 'models' (output from \code{\link{fit_models}}), 'model_summary' (output from \code{\link{get_params}}) and
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
#'      (applicable to the scale and shape parameters for independent shape models, applicable to the scale parameter only for the common shape
#'      model and not applicable for the separate model).}
#'
#' @export
runPSM <- function(data,
                   time_var, event_var,
                   model.type,
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


  #test that a legitimate value for model type has been provided
  assertthat::assert_that(
    all(model.type %in% c('Common shape', 'Independent shape', 'Separate')),
    msg = "Only the following model types are supported are supported: 'Common shape', 'Independent shape', 'Separate' "
  )



  #For models with common shape, calculate the location
  #parameter for the treatment arm for each distribution

  if(model.type == 'Common shape'){
    output <- run_common_shape(data, time_var, event_var,distr,strata_var, int_name, ref_name)
  }

  # For separate models split the data by treatment and for 2 separate models
  # for each distribution

  if(model.type == 'Separate'){
    output <- run_separate(data, time_var, event_var,distr,strata_var, int_name, ref_name)
  }

  #For models with independent shape calculate the scale and shape parameters for the
  #treatment arm for each distribution
  if(model.type == 'Independent shape'){
    output <- run_independent_shape(data, time_var, event_var,distr,strata_var, int_name, ref_name)
  }

return(output)
}
