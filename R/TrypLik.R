#' Calculate Tryptase Variant Likelihoods
#'
#' @param wt wildtype
#' @param fs frameshift
#' @param areads alpha counts
#' @param breads beta counts
#' @param dreads delta counts
#' @param pop population to use (1=, 2=)
#'
#' @return a `data.frame`
#'
#' @export
TrypLik <- function(wt, fs, areads, breads, dreads, pop = 0) {
  .Call("C_tryplik",
        as.integer(wt), as.integer(fs),
        as.integer(areads), as.integer(breads), as.integer(dreads),
        as.integer(pop))
}
