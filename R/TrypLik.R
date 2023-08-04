#' TrypLik
#'
#' @export
TrypLik <- function(wt, fs, areads, breads, dreads, pop = 0) {
  .Call("C_tryplik",
        as.integer(wt), as.integer(fs),
        as.integer(areads), as.integer(breads), as.integer(dreads),
        as.integer(pop))
}
