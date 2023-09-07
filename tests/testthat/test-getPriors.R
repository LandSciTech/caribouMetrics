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