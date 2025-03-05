# default input data
mod <- caribouBayesianPM(niters=100)


mod_tab <- suppressWarnings(getOutputTables(mod, simInitial = getSimsInitial()))

param_nms <- c(
  "Adult female survival", "Recruitment", "Population growth rate", 
  "Female population size",
  "Adjusted recruitment"
)

test_that("simplest plot works", {
  plotRes(mod_tab, parameter = "Recruitment") %>% 
    expect_s3_class("ggplot")
  
  # also works with out obs or sim
  mod_tab2 <- mod_tab
  mod_tab2$obs.all <- NULL
  mod_tab2$sim.all <- NULL
  
  plotRes(mod_tab2, parameter = "Recruitment") %>% 
    expect_s3_class("ggplot")
})

test_that("all parameters work",{
  plotRes(mod_tab, param_nms) %>% expect_type("list")
})

