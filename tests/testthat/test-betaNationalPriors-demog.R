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

test_that("bbouNationalPriors works", {
  lowA <- bbouNationalPriors(Anthro = 1:6, fire_excl_anthro = 1:6*5/10)
  
  highA <- bbouNationalPriors(Anthro = 91:96, fire_excl_anthro = 1:6*5/10)
  
  # survival is lower when anthro is high
  expect_gt(lowA$priors_survival["b0_mu"], highA$priors_survival["b0_mu"])
  
  # annual survival is lower than monthly
  annual <- bbouNationalPriors(Anthro = 1:6, fire_excl_anthro = 1:6*5/10, 
                                           month = FALSE)
  
  expect_gt(lowA$priors_survival["b0_mu"], annual$priors_survival["b0_mu"])
  
})
