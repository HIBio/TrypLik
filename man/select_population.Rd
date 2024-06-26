% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrypLik.R
\name{select_population}
\alias{select_population}
\title{Select a genetic-ancestry population}
\usage{
select_population(pop = NULL)
}
\arguments{
\item{pop}{population to use (0=EUR, 1=AFR, 2=EAS, 3=SAS) as either a number
or character string}
}
\value{
a named integer mapping the population name to the index
}
\description{
Select a genetic-ancestry population
}
\examples{
# by index
select_population(pop = 0)

# by name
select_population(pop = "EUR")

}
