% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{genomes}
\alias{genomes}
\title{1000 Genomes Tryptase Data}
\format{
\subsection{\code{genomes}}{

A data frame with 503 rows and 3 columns:
\describe{
\item{2A/D}{double ratio of alpha to delta alleles}
\item{2B/D}{double ratio of beta to delta alleles}
\item{Geno Group}{Genotype classification group}
...
}
}
}
\source{
\url{https://www.internationalgenome.org/1000-genomes-summary/}
}
\usage{
data(genomes)
}
\description{
A subset of 503 samples from the 1000 Genomes Project, processed
for tryptase allele counts. See the background vignette of \code{TrypLik}
for additional details.
}
\examples{
data(genomes)
head(genomes)

}
\keyword{datasets}
