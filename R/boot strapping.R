#' Run a complete parametric survival analysis 
#'
#' This function is a wrapper for the \code{\link{runPSM}} function intended to be run with the \code{\link{boot}} package. 
#' It enables running a complete parametric survival analysis for use when performing bootstrapping to explore uncertainty.
#' By re-using random seeds for each bootstrap sample it is possible to maintain correlations across multiple endpoints.
#'
#' @param data A data frame containing individual patient data for the relevant
#'   time to event outcomes. This is passed to the \code{data} argument of the
#'   \code{\link{runPSM}} function
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param ... Additional parameters as used by \code{\link{runPSM}}
#' @details  This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing parametric survival modelling using
#'   {\link{flexsurv}} and returning the parameters of the survival distributions.
#'    This is used as the 'statistic' argument in the boot function.
#' @return The 'parameters_vector' object from the \code{\link{runPSM}} function.
#'
#'    'parameters_vector' is a vector which contains the coefficients for all of the flexsurv models specified.
#'    The column names are in the format 'modeltype.distribution.parameter.TreatmentName', for example, comshp.weibull.shape.Int refers to the shape parameter
#'     of the common shape weibull distribution for the intervention treatment and 'indshp.gengamma.mu.ref' refers to the mu parameter
#'     of the independent shape generalised gamma distribution for the reference treatment. Columns with 'TE' at the end are the treatment effect coefficients
#'      (applicable to the scale and shape parameters for independent shape models, applicable to the scale parameter only for the common shape
#'      model and not applicable for the separate model).
#'
#' @export
bootPSM <- function(data,
                    i,
                   ...){

  data_boot <- data[i,]
  # call run_PSM (this does validity checks on inputs so no need to duplicate checks here)
  tryCatch({
    output <- runPSM(data = data_boot, ...)
    return(output$parameters_vector)
  })
}

#' Process bootstrap samples
#'
#' This function does some minimal data manipulation to make the object returned by 
#' bootPSM more useful and to match the \code{parameters} data frame that the \code{\link{runPSM}} function returns
#'
#' @param X The object created by running \code{\link{boot}} with \code{\link{runPSM}}
#' @details  This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing parametric survival modelling using
#'   {\link{flexsurv}} and returning the parameters of the survival distributions.
#'    This is used as the 'statistic' argument in the boot function.
#' @return A data frame similar to the 'parameters' object from the \code{\link{runPSM}} function.
#'
#' @export
bootPSMtidy <- function(X){
  
  assertthat::assert_that(
    X$call$statistic == "bootPSM", 
    msg = "This function only works with bootstrap samples created using boot with statistic = bootPSM"
    )
  
  
  # extract the samples, include the main estimates and add names
  samples <- as.data.frame(rbind(X$t0, X$t))
  
  # extract the model types and transpose - can be com, ind or sep
  
  params <- tibble()
  
  # any common shape models?
  if (any(grepl("comshp.",names(samples), fixed = TRUE))){
    
    com.samples <- samples[,grepl("comshp.",names(samples), fixed = TRUE)]
  
    # rename  
    names(com.samples) <- gsub("comshp.", "", names(com.samples))  
    
    # add metadata
    com.samples <- com.samples %>%
      mutate(SampleID = 0:X$R,
             Model = "Common shape", 
             Reference_name = X$call$ref_name,
             Intervention_name = X$call$int_name
             )
    
    params <- dplyr::bind_rows(params, com.samples)
  }
  
  # any independent shape models?
  if (any(grepl("indshp.",names(samples), fixed = TRUE))){
    
    ind.samples <- samples[,grepl("indshp.",names(samples), fixed = TRUE)]
    
    # rename  
    names(ind.samples) <- gsub("indshp.", "", names(ind.samples))  
    
    # add metadata
    ind.samples <- ind.samples %>%
      mutate(SampleID = 0:X$R,
             Model = "Independent shape", 
             Reference_name = X$call$ref_name,
             Intervention_name = X$call$int_name
             )
    
    params <- dplyr::bind_rows(params, ind.samples)
  }
  
  # any separate models?
  if (any(grepl("sep.",names(samples), fixed = TRUE))){
    
    sep.samples <- samples[,grepl("sep.",names(samples), fixed = TRUE)]
    
    # rename  
    names(sep.samples) <- gsub("sep.", "", names(sep.samples))  
    
    # add metadata
    sep.samples <- sep.samples %>%
      mutate(SampleID = 0:X$R,
             Model = "Seperate", 
             Reference_name = X$call$ref_name,
             Intervention_name = X$call$int_name
             )
    
    params <- dplyr::bind_rows(params, sep.samples)
  }
  
  # any one arm models?
  if (any(grepl("onearm.",names(samples), fixed = TRUE))){
    
    onearm.samples <- samples[,grepl("onearm.",names(samples), fixed = TRUE)]
    
    # rename  
    names(onearm.samples) <- gsub("onearm.", "", names(onearm.samples))  
    
    # add metadata
    onearm.samples <- onearm.samples %>%
      mutate(SampleID = 0:X$R,
             Model = "One arm", 
             Reference_name = X$call$ref_name,
             Intervention_name = X$call$int_name
      )
    
    params <- dplyr::bind_rows(params, onearm.samples)
  }
  
  # return the processed dataframe and reorder
  metadata.cols <- which(names(params) %in% c("Model", "Reference_name", "Intervention_name", "SampleID"))
  data.cols <- which(!names(params) %in% c("Model", "Reference_name", "Intervention_name", "SampleID"))
  
  # add metadata and stack
  params_df <- params[,c(metadata.cols, data.cols)] %>%
    dplyr::arrange(SampleID, Model)
  
  return(params_df)
}

