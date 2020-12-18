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
                   distr = c('weibull',
                             'gompertz',
                             'lnorm',
                             'llogis',
                             'gengamma',
                             'gamma'),
                   strata_var,
                   int_name, ref_name){

  #test that only valid distributions have been provided
  #This is also tested within fit_models. Consider eliminating here to avoid redundancy
  allowed_dist <- c('weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma')
  assertthat::assert_that(
    all(distr %in% allowed_dist),
    msg = "Only the following distributions are supported: 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma' "
  )

  # standardise variable names
  data_standard=Format_data(data, time_var, event_var, strata_var, int_name, ref_name)

  # Model formulas
  model.formula.shape = Surv(Time, Event==1) ~ ARM + shape(ARM)
  model.formula.sdlog = Surv(Time, Event==1) ~ ARM + sdlog(ARM)
  model.formula.sigma_Q = Surv(Time, Event==1) ~ ARM + sigma(ARM) + Q(ARM)



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

  }

  #Re-organise columns
  param_out <- cbind(param_out.weib, param_out.llogis, param_out.gomp, param_out.gamma, param_out.lnorm, param_out.gengamma)


  colnames(param_out) <- gsub("int", int_name, colnames(param_out))
  colnames(param_out) <- gsub("ref", ref_name, colnames(param_out))



  #Re-organise columns
  param_final <- param_out %>%
    #groups parameters by distribution in the order given in the dist argument
    dplyr::mutate(Model="Independent shape", `Intervention name`=int_name, `Reference name`=ref_name)

  models <- list(weibull=models.weib$weibull, llogis=models.llogis$llogis, gompertz=models.gomp$gompertz, gamma=models.gamma$gamma, lnorm=models.lnorm$lnorm, gengamma=models.gengamma$gengamma)

  params <- bind_rows(params.weib, params.llogis, params.gomp, params.gamma, params.lnorm, params.gengamma)

  #collect and return output
  output <- list(
    models = models,
    model_summary = params,
    parameters = param_final
  )
  return(output)
}
