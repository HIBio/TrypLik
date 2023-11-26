## code to prepare `genomes` dataset goes here

genomes <- read.csv(system.file("extdata", "1KGP_EUR_figure_raw.csv", package = "TrypLik"))
colnames(genomes) <- c("2A/D", "2B/D", "Geno Group")
lvls <- c(
    "0α:4β" = 0,
    "1α:3β" = 1,
    "2α:2β" = 2,
    "everything else" = 3
)
genomes$`Geno Group` <- names(lvls[match(genomes$`Geno Group`, lvls)])

usethis::use_data(genomes, overwrite = TRUE)
