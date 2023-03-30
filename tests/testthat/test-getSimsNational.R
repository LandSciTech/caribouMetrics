test_that("cacheing happens", {
  expect_message(getSimsNational(), "Updating cached")
  expect_message(getSimsNational(), "saved object")
})
