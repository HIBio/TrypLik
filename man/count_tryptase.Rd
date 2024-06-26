% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_tryptase.R
\name{count_tryptase}
\alias{count_tryptase}
\title{Count Tryptase Alleles in a .sam File}
\usage{
count_tryptase(sam_file)
}
\arguments{
\item{sam_file}{path to the temp.sam file generated by the pre-processing}
}
\value{
A vector of counts suitable for input to \code{TrpLik}. See Details.
}
\description{
Count Tryptase Alleles in a .sam File
}
\details{
The first four of values are the numbers of reads mapping uniquely to ⍺,
βI, βII or βIII alleles in the specific 83 bp region where
all four alleles can be distinguished from each other.

The next two are the numbers of reads that map uniquely to β ‘wildtype’
(i.e. not frameshifted) and β frameshifted sequences.

The final three are the numbers of reads that can be uniquely mapped to
⍺, β, or δ alleles respectively.
}
\examples{
\dontrun{
count_tryptase("temp.sam")
}
}
