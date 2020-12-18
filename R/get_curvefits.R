#' Get curve fits from \code{\link{fit_models}} objects
#'
#' Manipulates \code{\link{fit_models}} objects to get the curve fits.
#'
#' @param models Object from \code{\link{fit_models}}
#' @param time list of times
#'
#' @return A data frame containing time and survival from flexsurv objects.
#'
#' @seealso \code{\link{fit_models}} \code{\link{flexsurvreg}}
#'
#' @export
get_curvefits <- function(models, time) {
  output <-   tibble::enframe(models) %>%
    dplyr::mutate(
      curvefits= purrr::map(value, summary, type="survival",B=0,t=time, ci = FALSE, tidy = T)
    )
  names(output$curvefits) <- names(models)
    output
}
