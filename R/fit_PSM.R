
#' Fit flexsurv models
#'
#' Fits single or multiple \code{\link{flexsurv}} models using \code{\link{flexsurvreg}}.
#'
#' @param model.formula A survival model formula in the general form shown
#'   below. Note that variable names must match the corresponding columns in the
#'   data. This is passed to the \code{formula} argument of the
#'   \code{\link{fit_models}} function
#'  \itemize{
#'   \item Surv(Time, Event==1) ~ ARM is a model with a single covariate for the
#'   effect of treatment
#'   \item Surv(Time, Event==1) ~ 1 is a model with no covariates typically
#'   fitted to data from a single treatment group
#'  }
#' @param data A data frame containing individual patient data for the relevant
#'   time to event outcomes. This is passed to the \code{data} argument of the
#'   \code{\link{fit_models}} function
#' @param dist A vector string of distributions, see dist argument in
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
#'   \item Generalised F ('genf')
#'   }

#'
#' @return A list containing flexsurv objects.
#' @seealso \code{\link{flexsurvreg}}
#'
#' @export
fit_models <- function(model.formula,
                       distr = c('exp',
                                 'weibull',
                                 'gompertz',
                                 'lnorm',
                                 'llogis',
                                 'gengamma',
                                 'gamma',
                                 'genf'),
                       data) {


  #test that only valid distributions have been provided
  allowed_dist <- c('exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma', 'genf')
  assertthat::assert_that(
    all(distr %in% allowed_dist),
    msg = "Only the following distributions are supported: 'exp', 'weibull', 'gompertz', 'lnorm', 'llogis', 'gengamma', 'gamma', 'genf'"
  )

  runFLEX <- function(dist
                      ) {
    model <- flexsurv::flexsurvreg(formula=model.formula, data=data, dist=dist)
  }

  #list of flexsurv objects
  output <- lapply(distr, function(x) runFLEX(x))
  names(output) <- distr
  output
}

#' Get parameter estimates from \code{\link{fit_models}} objects
#'
#' Manipulates \code{\link{fit_models}} objects to get the parameter estimates and AIC and BIC values.
#'
#' @param models Object from \code{\link{fit_models}}
#'
#' @return A tibble object containing the fitted model objects, the parameter
#'   estimates (\code{\link{coef}}),  \code{\link{AIC}} and \code{\link{BIC}}
#'   from flexsurv objects.
#'
#' @seealso \code{\link{fit_models}} \code{\link{flexsurvreg}}
#'
#' @export
get_params <- function(models) {

  #test inputs before proceeding
  input.class <- sapply(models, class)
  assertthat::assert_that(
    all(input.class == "flexsurvreg"),
    msg = "get_params expects a list of 'flexsurvreg' objects as input. At least one of your inputs is not a flexsurvreg object"
  )

  output <-   tibble::enframe(models) %>%
    dplyr::mutate(
      coef= lapply(models, coef), #get coefficient estimates for each model
      AIC = sapply(models, AIC), #get AIC
      BIC = sapply(models, BIC) #get BIC
    )
}


