## code to prepare `resTypeCode` dataset goes here
resTypeCode <- data.frame(ResourceType = c("CON", "DEC", "DTN", "LGOP", "LGTP",
                                      "LGW", "MIX", "ST", "other"), 
                          code = 1:9, stringsAsFactors = FALSE)

usethis::use_data(resTypeCode, overwrite = TRUE)
