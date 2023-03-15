test_that("cacheing happens", {
  expect_warning(expect_message(getSimsNational(), "will be saved"), "expected survival")
  expect_message(getSimsNational(), "saved object")
})
