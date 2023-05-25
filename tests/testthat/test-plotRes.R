test_that("simplest plot works", {
  # default input data
  mod <- caribouBayesianIPM(Nchains = 1, Niter = 100, Nburn = 10,
                            Nthin = 2)

  
  mod_tab <- getOutputTables(mod, simNational = getSimsNational(),
                             getKSDists = FALSE)

  plotRes(mod_tab, parameter = "Recruitment") %>% 
    expect_s3_class("ggplot")
  
  # also works with out obs or sim
  mod_tab$obs.all <- NULL
  mod_tab$sim.all <- NULL
  
  plotRes(mod_tab, parameter = "Recruitment") %>% 
    expect_s3_class("ggplot")
})
