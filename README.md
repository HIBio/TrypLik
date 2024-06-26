
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
[CRAN](http://cran.r-project.org/).

<!-- Then install `TrypLik` from [Bioconductor](http://bioconductor.org/) using the following code: -->
<!-- ```{r 'install', eval = FALSE} -->
<!-- if (!requireNamespace("BiocManager", quietly = TRUE)) { -->
<!--     install.packages("BiocManager") -->
<!-- } -->
<!-- BiocManager::install("TrypLik") -->
<!-- ``` -->

Install the development version from
[GitHub](https://github.com/HIBio/TrypLik) with:

``` r
BiocManager::install("HIBio/TrypLik")
```

or, if not using the Bioconductor ecosystem:

``` r
remotes::install_github("HIBio/TrypLik")
```

## Preprocessing

Provided the sequencing data has been processed (see the background
vignette included in this package, and accompanying shell scripts in
`inst/`) to produce the reads mapped to specific tryptase alleles, these
can be provided to the `TrypLik` function

``` r
library("TrypLik")

TrypLik(26, 0, 347, 316, 304)
#>   Alpha_count Beta_count Beta_FS_count Posterior_likelihood
#> 1           2          2             0              0.99890
#> 2           2          3             0              0.00008
#> 3           3          2             0              0.00082
#> 4           3          3             0              0.00020
```

This function is an R wrapper for the C function, producing a native
`data.frame`.

The C code is provided as a standalone source file in
[`inst/tryplik.c`](inst/tryplik.c) in which the R-related code has been
commented out. The R wrapper uses the file
[`src/tryplik.c`](src/tryplik.c) which only differs in the additional R
components.

Once compiled, the standalone version can be used in a command-line
pipeline; following the processing outlined in the background vignette,
one could use

``` bash
[...]
samtools markdup -r positionsort.sam temp.sam
./count_tryptase >> tryptase.out

cat tryptase.out | awk '{printf "./Tryplik %d %d %d %d %d | sort -nk2 | tail -n 1\n", $6, $7, $8, $9, $10}'
```

The processing can be entirely performed from R, provided that both
`samtools` and `bwa` are installed. To preprocess the file
`HG00100.final.cram` from 1000 Genomes, producing a `.sam` file suitable
for counting

``` r
sam_file <- preprocess("HG00100.final.cram")
# Processing... this may take some time
# [bwa_index] Pack FASTA... 0.00 sec
# [bwa_index] Construct BWT for the packed sequence...
# [bwa_index] 0.00 seconds elapse.
# [bwa_index] Update BWT... 0.00 sec
# [bwa_index] Pack forward-only FASTA... 0.00 sec
# [bwa_index] Construct SA from BWT and Occ... 0.00 sec
# [main] Version: 0.7.17-r1188
# [main] CMD: bwa index consensus.fa
# [main] Real time: 0.004 sec; CPU: 0.004 sec
# [M::bam2fq_mainloop] discarded 228 singletons
# [M::bam2fq_mainloop] processed 27228 reads
# [M::bwa_idx_load_from_disk] read 0 ALT contigs
# [M::process] read 27000 sequences (4050000 bp)...
# [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (0, 1190, 0, 0)
# [M::mem_pestat] skip orientation FF as there are not enough pairs
# [M::mem_pestat] analyzing insert size distribution for orientation FR...
# [M::mem_pestat] (25, 50, 75) percentile: (348, 429, 530)
# [M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 894)
# [M::mem_pestat] mean and std.dev: (422.49, 140.29)
# [M::mem_pestat] low and high boundaries for proper pairs: (1, 1076)
# [M::mem_pestat] skip orientation RF as there are not enough pairs
# [M::mem_pestat] skip orientation RR as there are not enough pairs
# [M::mem_process_seqs] Processed 27000 reads in 1.101 CPU sec, 1.064 real sec
# [main] Version: 0.7.17-r1188
# [main] CMD: bwa mem consensus.fa tempR1.fastq tempR2.fastq
# [main] Real time: 1.097 sec; CPU: 1.133 sec
# result file temp.sam is here: temp.sam

sam_file
# [1] "temp.sam"

counts <- count_tryptase(sam_file)
counts
# [1]  24   0  22   0  26   0 347 316 304

TrypLik(counts)
#   Alpha_count Beta_count Beta_FS_count Posterior_likelihood
# 1           2          2             0              0.99890
# 2           2          3             0              0.00008
# 3           3          2             0              0.00082
# 4           3          3             0              0.00020
```

(the `TrypLik` function takes either a single vector of counts or values
provided to specific arguments).

## Citation

Below is the citation output from using `citation('TrypLik')` in R.
Please run this yourself to check for any updates on how to cite
**TrypLik**.

``` r
print(citation("TrypLik"), bibtex = TRUE)
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

- Package development is possible thanks to
  *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.16/BiocCheck)*.
- The [documentation website](http://HIBio.github.io/TrypLik) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

<!-- For more details, check the `dev` directory. -->

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.16/biocthis)*.
