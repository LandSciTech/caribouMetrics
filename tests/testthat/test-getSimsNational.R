test_that("cacheing happens", {
  # clear the cache
  rm(list = ls(envir = cacheEnv), envir = cacheEnv)
  expect_message(getSimsNational(), "Updating cached")
  expect_message(getSimsNational(), "saved object")
})
