library(RGAP)

data("gap")
tsList <- amecoData2input(gap[["France"]], alpha = 0.65)
exoType <- initializeExo(varNames = "ws", D = 2, L = 0)
model <- NAWRUmodel(tsl = tsList, exoType = exoType)
parRestr <- initializeRestr(model = model, type = "hp")

test_that("NAWRU model", {

  expect_s3_class(model, "NAWRUmodel")
  expect_equal(is.NAWRUmodel(model, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = model, cran = FALSE)

})

test_that("NAWRU MLE fit", {

  skip_on_cran()
  f <-  fit(model, parRestr = parRestr)

  expect_s3_class(f, "NAWRUfit")
  expect_equal(RGAP:::is.NAWRUfit(f, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = f, cran = FALSE)

})

test_that("NAWRU bayesian fit", {

  skip_on_cran()
  f <-  fit(model, parRestr = parRestr)
  set.seed(5)
  fBayes <- fit(model = model, method = "bayesian", R = 10000, thin = 10, MLEfit = f, burnin = 3000)

  expect_s3_class(fBayes, "NAWRUfit")
  expect_equal(RGAP:::is.NAWRUfit(fBayes, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = fBayes, cran = FALSE)

})
