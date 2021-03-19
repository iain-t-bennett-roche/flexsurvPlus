
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
  
  
  
  runFLEX <- function(dist){
    tryCatch(model <- flexsurv::flexsurvreg(formula=model.formula, data=data, dist=dist),
             error = function(e){
               message("An error occurred in ",dist,":\n", e)
               return("error")
             },
             warning = function(w){
               message("A warning occured in ",dist,":\n",  w)
               return("warning")
             },
             finally = {
               message("Fitting model ", dist)
               
             })
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
#' @return A data frame containing the \code{\link{AIC}}, \code{\link{BIC}}
#'   from flexsurv objects and a status column of whether the model converged or produced an error or warning.
#'
#' @seealso \code{\link{fit_models}} \code{\link{flexsurvreg}}
#'
#' @export
get_params <- function(models) {
  
  # Filter on flexsurv models
  flexsurvreg.test <- sapply(models, function(x) class(x)=="flexsurvreg")
  models.flexsurv  <- models[flexsurvreg.test]
  
  #models that produced an error or warning
  models.nonconverge <- models[flexsurvreg.test==FALSE]
  
  #test inputs before proceeding
  input.class <- sapply(models.flexsurv, class)
  assertthat::assert_that(
    all(input.class == "flexsurvreg"),
    msg = "get_params expects a list of 'flexsurvreg' objects as input. At least one of your inputs is not a flexsurvreg object"
  )
  
  output <-   tibble::enframe(models.flexsurv) %>%
    dplyr::mutate(
      Dist = name,
      AIC = sapply(models.flexsurv, AIC), #get AIC
      BIC = sapply(models.flexsurv, BIC), #get BIC
      Status = "Converged"
    ) %>%
    dplyr::select(-value, -name) %>% 
    as.data.frame()
  
  if(length(models.nonconverge)>0){
  non_converged_output  <- data.frame(
    Dist = names(models.nonconverge),
    AIC = NA,
    BIC = NA,
    Status = unlist(models.nonconverge))
  rownames(non_converged_output) <- c()
  
  output <- rbind(output, non_converged_output)
  }
  
  output
  
}

