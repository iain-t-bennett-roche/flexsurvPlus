#' Run a complete parametric survival analysis for two endpoints with multiple distributions
#'
#' Run a complete parametric survival analysis for two endpoints with multiple
#' distributions. For use when performing bootstrapping on two correlated
#' endpoints such as PFS and OS.
#'
#' @param data A data frame containing individual patient data for the relevant
#'   time to event outcomes. This is passed to the \code{data} argument of the
#'   \code{\link{fit_models}} function
#' @param endpoint1 Character to indicate the name of the first endpoint, used for labelling.
#' @param endpoint2 Character to indicate the name of the second endpoint, used for labelling.
#' @param time_var Name of time variable in 'data' of first endpoint. Variable must be numerical and >0.
#' @param  event_var Name of event variable in 'data' of first endpoint.
#'   Variable must be numerical and contain 1's to indicate an event and 0 to
#'   indicate a censor.
#' @param time_var2 Name of time variable in 'data' of second endpoint. Variable must be numerical and >0.
#' @param  event_var2 Name of event variable in 'data' of second endpoint.
#'   Variable must be numerical and contain 1's to indicate an event and 0 to
#'   indicate a censor.
#' @param model.type Character vector indicating the type of model formula provided. Permitted values are
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
#'   must be a level of the "strata_var" column in "data", used for labelling
#'   the parameters.
#' @param distr A vector string of distributions, see dist argument in
#'   \code{\link{flexsurvreg}}. This is passed to the \code{distr} argument of
#'   the \code{\link{fit_models}} function
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @details  This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing parametric survival modelling using
#'   {\link{flexsurv}} and returning the parameters of the survival distributions.
#'    This is used as the 'statistic' argument in
#'   the boot function.
#'
#' #' Possible distributions include:
#' \itemize{
#'   \item Exponential ('exp')
#'   \item Weibull ('weibull')
#'   \item Gompertz ('gompertz')
#'   \item Log-normal ('lnorm')
#'   \item Log-logistic ('llogis')
#'   \item Generalized gamma ('gengamma')
#'   \item Gamma ('gamma')
#'   }
#' @return The 'parameters' object from the \code{\link{runPSM}} function.
#'
#'    'parameters' is a data frame with with one row which contains the coefficients for all of the flexsurv models specified.
#'    The column names are in the format 'distribution.parameter.TreatmentName', for example, weibull.shape.Intervention refers to the shape parameter
#'     of the weibull distribution for the intervention treatment and 'gengamma.mu.Reference' refers to the mu parameter
#'     of the generalised gamma distribution for the reference treatment. Columns with 'TE' at the end are the treatment effect coefficients
#'      (applicable to the scale and shape parameters for independent shape models, applicable to the scale parameter only for the common shape
#'      model and not applicable for the separate model).
#'
#' @export
bootPSM <- function(data,
                   endpoint1, endpoint2,
                   time_var, event_var,
                   time_var2 = NULL, event_var2 = NULL,
                   model.type,
                   distr = c('exp',
                             'weibull',
                             'gompertz',
                             'lnorm',
                             'llogis',
                             'gengamma',
                             'gamma'),
                   strata_var,
                   int_name, ref_name,
                   i){

  data_boot <- data[i,]

  # call run_PSM (this does validity checks on inputs so no need to duplicate checks here)
  # check if have 1 or more time vars included
  tryCatch({
    
    if (!is.null(time_var2)){
      # 2 endpoints - assume correlated
      
      output1 <- run_PSM(data_boot, time_var = time_var, event_var = event_var, 
                         distr = distr, strata_var = strata_var, int_name = int_name, ref_name  = ref_name, model.type = model.type)
      
      
      output2 <- run_PSM(data_boot, time_var = time_var2, event_var = event_var2, 
                         distr = distr, strata_var = strata_var, int_name = int_name, ref_name  = ref_name, model.type = model.type)
      
      # remove labels
      output1$parameters <- output1$parameters %>%
        select(-Model, -`Intervention name`, -`Reference name`)
      
      output2$parameters <- output2$parameters %>%
        select(-Model, -`Intervention name`, -`Reference name`)
      
      colnames(output1$parameters) <- paste0(colnames(output1$parameters), ".", endpoint1)
      colnames(output2$parameters) <- paste0(colnames(output2$parameters), ".", endpoint2)
      
      
      all_params <- cbind(output1$parameters, output2$parameters)
      all_params_out <- data.matrix(all_params)
      
    } else{
      # only 1 endpoints so format is matching standard return format
      output <- runPSM(data_boot, time_var = time_var, event_var = event_var, 
                       distr = distr, strata_var = strata_var, int_name = int_name, ref_name  = ref_name, model.type = model.type)
      
      params 
    }
  
  
   
    return(all_params_out)
  })
}
