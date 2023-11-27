test_that("invalid input fails", {
    expect_error(TrypLik(0), regexp = "argument.*missing")
    expect_error(TrypLik(26, 0, 347, 316, 304, pop = 8), regexp = "pop must be between")
    expect_error(TrypLik(26, 0, 347, 316, pop = 0), regexp = "argument.*missing")
    expect_error(TrypLik(c(26, 0, 347, 316, 304, 999), pop = 0), regexp = "Unrecognised input format")
    expect_error(TrypLik(c(26, 0, 347), pop = 0), regexp = "Unrecognised input format")
})

test_that("correct usage works silently", {
    expect_silent(TrypLik(26, 0, 347, 316, 304))
    expect_silent(TrypLik(c(26, 0, 347, 316, 304)))
    expect_silent(TrypLik(26, 0, 347, 316, 304, pop = 0))
    expect_silent(TrypLik(c(26, 0, 347, 316, 304), pop = 0))
    expect_silent(TrypLik(26, 0, 347, 316, 304, pop = 1))
    expect_silent(TrypLik(c(26, 0, 347, 316, 304), pop = 1))
})

test_that("default argument works", {
    expect_identical(TrypLik(26, 0, 347, 316, 304), TrypLik(26, 0, 347, 316, 304, pop = 0))
})

test_that("appropriate structure is returned", {
    expect_s3_class(TrypLik(26, 0, 347, 316, 304), "data.frame")
    expect_s3_class(TrypLik(c(26, 0, 347, 316, 304)), "data.frame")

    expect_identical(colnames(TrypLik(26, 0, 347, 316, 304)), c("Alpha_count", "Beta_count", "Beta_FS_count", "Posterior_likelihood"))
})

test_that("preprocess fails gracefully when cram file is missing", {
  expect_error(preprocess("foo.cram"))
})
