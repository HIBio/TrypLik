
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TrypLik

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/HIBio/TrypLik)](https://github.com/HIBio/TrypLik/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/HIBio/TrypLik)](https://github.com/HIBio/TrypLik/pulls)
<!-- badges: end -->

The goal of `TrypLik` is to support tryptase copy number estimation from
sequence data.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `TrypLik` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("TrypLik")
```

And the development version from
[GitHub](https://github.com/HIBio/TrypLik) with:

``` r
BiocManager::install("HIBio/TrypLik")
```

## Example

Provided the sequencing data has been processed (see the background
vignette included in this package, and accompanying shell scripts in
`inst/`) to produce the reads mapped to specific tryptase alleles, these
can be provided to the `TrypLik` function

``` r
library("TrypLik")

TrypLik(27, 0, 382, 356, 357)
#>   Alpha_count Beta_count Beta_FS_count Posterior_likelihood
#> 1           2          2             0              0.99985
#> 2           2          3             0              0.00002
#> 3           3          2             0              0.00011
#> 4           3          3             0              0.00002
```

This function is an R wrapper for the C function, producing a native
`data.frame`.

The C code is provided as a standalone source file in `inst/tryplik.c`
in which the R-related code has been commented out. The R wrapper uses
the file `src/tryplik.c` which only differs in the additional R
components.

Once compiled, the standalone version can be used in a command-line
pipeline; following the processing outlined in the background vignette,
one could use

``` bash
samtools markdup -r positionsort.sam temp.sam
./in_tryptasenew2 >> tryptase.out

cat tryptase.out | awk '{printf "./Tryplik %d %d %d %d %d | sort -nk2 | tail -n 1\n", $6, $7, $8, $9, $10}'
```

## Citation

Below is the citation output from using `citation('TrypLik')` in R.
Please run this yourself to check for any updates on how to cite
**TrypLik**.

``` r
print(citation('TrypLik'), bibtex = TRUE)
#> Warning in citation("TrypLik"): no date field in DESCRIPTION file of package
#> 'TrypLik'
#> Warning in citation("TrypLik"): could not determine year for 'TrypLik' from
#> package DESCRIPTION file
#> 
#> To cite package 'TrypLik' in publications use:
#> 
#>   Carroll J, Wall J (????). _TrypLik: Calculates Tryptase Variant
#>   Likelihoods_. R package version 0.1.0.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {TrypLik: Calculates Tryptase Variant Likelihoods},
#>     author = {Jonathan Carroll and Jeff Wall},
#>     note = {R package version 0.1.0},
#>   }
```

Please note that the `TrypLik` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `TrypLik` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.16/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://HIBio.github.io/TrypLik) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.16/biocthis)*.
