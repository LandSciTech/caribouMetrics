# get file from extdata
fl <- system.file("extdata/FNLC_Lookup_table.csv", package = "caribouMetrics")
fnlcToResType <- read.csv(fl)
usethis::use_data(fnlcToResType, overwrite = TRUE)
# have to run threshTable script to add to the internal data