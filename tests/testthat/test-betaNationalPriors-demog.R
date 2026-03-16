test_that("defaults work", {
  expect_is(betaNationalPriors(), "list")
})

test_that("results are as expected", {
  # when mod is higher result is higher
  defs <- eval(formals(betaNationalPriors)) %>% as.list() %>% lapply(eval)
  defs1 <- rapply(defs, function(x){x+0.5}, classes = "numeric", how = "replace")
  
  defPriors <- betaNationalPriors(defs)
  def1Priors <- betaNationalPriors(defs1)
  purrr::map2(defPriors, def1Priors, ~(.y >= .x)) %>% unlist() %>% all() %>% 
    expect_true()
})

test_that("bbouNationalPriors works",{
  bbou_pr <- bbouNationalPriors(5, 6, month = "both")
  
  expect_named(bbou_pr, c("priors_recruitment", "priors_survival"))
  expect_named(bbou_pr$priors_survival, paste0(c("", "", "monthly_", "monthly_"),
                                               c("b0_mu", "b0_sd")))
})
