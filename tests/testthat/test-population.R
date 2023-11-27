test_that("population selection via index", {
    expect_identical(select_population(pop = 0), c(EUR = 0L))
    expect_identical(select_population(pop = 1), c(AFR = 1L))
    expect_identical(select_population(pop = 2), c(EAS = 2L))
    expect_identical(select_population(pop = 3), c(SAS = 3L))
})

test_that("population selection via name", {
    expect_identical(select_population(pop = "EUR"), c(EUR = 0L))
    expect_identical(select_population(pop = "AFR"), c(AFR = 1L))
    expect_identical(select_population(pop = "EAS"), c(EAS = 2L))
    expect_identical(select_population(pop = "SAS"), c(SAS = 3L))
})

test_that("invalid inputs fail", {
    expect_error(select_population(pop = "FOO"), regexp = "pop must be one of")
    expect_error(select_population(pop = "1"), regexp = "pop must be one of")
    expect_error(select_population(pop = 9), regexp = "pop must be between")
})
