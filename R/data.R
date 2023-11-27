#' 1000 Genomes Tryptase Data
#'
#' A subset of 503 samples from the 1000 Genomes Project, processed
#' for tryptase allele counts. See the background vignette of `TrypLik`
#' for additional details.
#'
#' @format ## `genomes`
#' A data frame with 503 rows and 3 columns:
#' \describe{
#'   \item{2A/D}{double ratio of alpha to delta alleles}
#'   \item{2B/D}{double ratio of beta to delta alleles}
#'   \item{Geno Group}{Genotype classification group}
#'   ...
#' }
#' @usage data(genomes)
#' @source <https://www.internationalgenome.org/1000-genomes-summary/>
#' @examples
#' data(genomes)
#' head(genomes)
#'
"genomes"

#' Retrieve the 1000 Genomes data, processed for tryptase allele counts
#'
#' See the background vignette of `TrypLik`
#' for additional details.
#'
#' @returns A `data.frame` of 2A/D, 2B/D, and the inferred genotype group for the
#'   503 relevant samples
#' @export
#' @importFrom utils data
#'
#' @examples
#' head(tryptase_1KG())
tryptase_1KG <- function() {
    genomes <- NULL
    utils::data("genomes", package = "TrypLik", envir = environment())
    genomes
}
