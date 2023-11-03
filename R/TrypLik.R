#' Calculate Tryptase Variant Likelihoods
#'
#' @param wt wildtype
#' @param fs frameshift
#' @param areads alpha counts
#' @param breads beta counts
#' @param dreads delta counts
#' @param pop population to use (0=EUR, 1=AFR, 2=EAS, 3=SAS) as either a number
#'   or character string
#'
#' @return a `data.frame`
#'
#' @export
TrypLik <- function(wt, fs, areads, breads, dreads, pop = 0) {

  if (length(wt) > 1) {
    inputs <- wt
    if (length(inputs) == 9) {
      wt <- inputs[5]
      fs <- inputs[6]
      areads <- inputs[7]
      breads <- inputs[8]
      dreads <- inputs[9]
    } else if (length(inputs) == 5) {
      wt <- inputs[1]
      fs <- inputs[2]
      areads <- inputs[3]
      breads <- inputs[4]
      dreads <- inputs[5]
    } else {
      stop("Unrecognised input format")
    }
  }

  pop <- select_population(pop)
  .Call("C_tryplik",
        as.integer(wt), as.integer(fs),
        as.integer(areads), as.integer(breads), as.integer(dreads),
        as.integer(pop))
}

#' Select a genetic-ancestry population
#'
#' @param pop population to use (0=EUR, 1=AFR, 2=EAS, 3=SAS) as either a number
#'   or character string
#' @importFrom stats setNames
#' @export
select_population <- function(pop = NULL) {
  supported_pops <- c("EUR", "AFR", "EAS", "SAS")
  stopifnot("pop must be provided" = !is.null(pop))
  if (is.numeric(pop)) {
    stopifnot("pop must be between 0 and 3" = pop >= 0 && pop <= length(supported_pops) - 1L)
    popname <- supported_pops[pop+1L]
  } else {
    stopifnot("pop must be one of `EUR, AFR, EAS, SAS`" = pop %in% supported_pops)
    popname <- pop
    pop <- which(popname == supported_pops)
  }
  stats::setNames(pop, popname)
}
