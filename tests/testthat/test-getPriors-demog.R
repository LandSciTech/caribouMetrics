test_that("defaults work", {
  expect_is(getPriors(), "list")
})

test_that("results are as expected", {
  # when mod is higher result is higher
  defs <- eval(formals(getPriors)) %>% as.list() %>% lapply(eval)
  defs1 <- rapply(defs, function(x){x+0.5}, classes = "numeric", how = "replace")
  
  defPriors <- getPriors(defs)
  def1Priors <- getPriors(defs1)
  purrr::map2(defPriors, def1Priors, ~(.y >= .x)) %>% unlist() %>% all() %>% 
    expect_true()
})

test_that("getBBNationalInformativePriors works", {
  lowA <- getBBNationalInformativePriors(Anthro = 1:6, fire_excl_anthro = 1:6*5/10)
  
  highA <- getBBNationalInformativePriors(Anthro = 91:96, fire_excl_anthro = 1:6*5/10)
  
  # survival is lower when anthro is high
  expect_gt(lowA$priors_survival["b0_mu"], highA$priors_survival["b0_mu"])
  
  # annual survival is lower than monthly
  annual <- getBBNationalInformativePriors(Anthro = 1:6, fire_excl_anthro = 1:6*5/10, 
                                           month = FALSE)
  
  expect_gt(lowA$priors_survival["b0_mu"], annual$priors_survival["b0_mu"])
  
})