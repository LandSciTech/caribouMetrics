test_that("cacheing happens", {
  expect_message(getSimsNational(), "will be saved")
  expect_message(getSimsNational(), "saved object")
})
