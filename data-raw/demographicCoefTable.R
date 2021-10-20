# get file from extdata
fl <- system.file("extdata/populationGrowthTable.csv", package = "caribouMetrics")
popGrowthTableJohnsonECCC <- read.csv(fl)
usethis::use_data(popGrowthTableJohnsonECCC, overwrite = TRUE)
