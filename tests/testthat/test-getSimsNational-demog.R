test_that("cacheing happens", {
  # clear the cache
  rm(list = ls(envir = cacheEnv), envir = cacheEnv)
  expect_message(expect_warning(getSimsNational()), "Updating cached")
  expect_message(getSimsNational(), "saved object")
})
