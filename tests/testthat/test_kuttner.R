library(RGAP)

data("gap")
tsList <- as.list(gap[["Netherlands"]][,c("cpih","gdp")])
tsList$infl <- diff(tsList$cpih)
model <- KuttnerModel(tsl = tsList)
parRestr <- initializeRestr(model, type = "hp")


test_that("Kuttner model", {

  expect_s3_class(model, "KuttnerModel")
  expect_equal(is.KuttnerModel(model, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = model, cran = FALSE)

})

test_that("Kuttner MLE fit", {

  skip_on_cran()
  f <- fit(model, parRestr = parRestr)

  expect_s3_class(f, "KuttnerFit")
  expect_equal(RGAP:::is.KuttnerFit(f, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = f, cran = FALSE)

})


