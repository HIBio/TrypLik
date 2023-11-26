# Files included with TrypLik

The following files are included in this package's `inst/` directory:

* [consensus.fa](consensus.fa): Manually compiled consensus file based on the [supplemental data][consensus] of [Maun et al. 2019][maun2019]
* [count_tryptase.sh](count_tryptase.sh): bash script to count tryptase alleles in the pre-processed `temp.sam` file. An R version of this is available as `count_tryptase()`
* [evaluate.sh](evaluate.sh): bash script to run a compiled `TrypLik` C program with inputs taken from `tryptase.out`, produced via pre-processing
* [preprocess_single.sh](preprocess_single.sh): bash script to pre-process a single `.cram` file
* [preprocess.sh](preprocess.sh): bash script to pre-process several `.cram` files

[consensus]: https://www.cell.com/cms/10.1016/j.cell.2019.09.009/attachment/3a53f18d-da9f-4021-929c-0e6231b9d7ad/mmc2.pdf
[maun2019]: https://www.cell.com/cell/fulltext/S0092-8674(19)31015-3
