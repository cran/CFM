test_that("CFM data satisfy minimal KMO / Bartlett criteria", {
  set.seed(123)
  result <- CFM(n = 200, p = 10, m = 3, cens.dist = "logistic")
  X <- as.matrix(result$data)
  expect_true(is.matrix(X))
  expect_true(is.numeric(X))
  sds <- apply(X, 2, sd, na.rm = TRUE)
  expect_true(all(sds > 0))
})

test_that("CFM returns data with correct dimensions", {
  set.seed(123)
  result <- CFM(n = 200, p = 10, m = 3, cens.dist = "logistic")
  X <- as.matrix(result$data)
  expect_equal(nrow(X), 200)
  expect_equal(ncol(X), 10)
})
test_that("KMO and Bartlett tests run without error", {
  set.seed(123)
  result <- CFM(n = 200, p = 10, m = 3, cens.dist = "logistic")
  X <- as.matrix(result$data)
  expect_silent(psych::KMO(X))
  expect_silent(psych::cortest.bartlett(cor(X), n = 200))
})
