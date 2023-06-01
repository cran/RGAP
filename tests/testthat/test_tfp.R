library(RGAP)

data("gap")
tsList <- amecoData2input(gap[["Italy"]], alpha = 0.65)
model <- TFPmodel(tsl = tsList, cycle = "RAR2")

test_that("TFP model", {

  expect_s3_class(model, "TFPmodel")
  expect_equal(is.TFPmodel(model, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = model, cran = FALSE)

})

test_that("TFP MLE fit", {

  skip_on_cran()
  f <- fit(model = model)

  expect_s3_class(f, "TFPfit")
  expect_equal(RGAP:::is.TFPfit(f, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = f, cran = FALSE)

})

test_that("TFP bayesian fit", {

  skip_on_cran()
  f <- fit(model = model)
  set.seed(5)
  fBayes <- fit(model = model, method = "bayesian", R = 1000, thin = 2, MLEfit = f)

  expect_s3_class(fBayes, "TFPfit")
  expect_equal(RGAP:::is.TFPfit(fBayes, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = fBayes, cran = FALSE)

})
