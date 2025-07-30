test_that("cacheing happens", {
  # clear the cache
  rm(list = ls(envir = cacheEnv), envir = cacheEnv)
  expect_message(expect_warning(getSimsInitial()), "Updating cached")
  expect_message(getSimsInitial(), "saved object")
})
