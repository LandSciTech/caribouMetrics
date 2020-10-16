## code to prepare `rfuToResType` 

# Convert LSL to csv look up table Forest Group Classification by Forest Unit
# (OLT - Northwest Only) copied from LSL lines 1991-2040 (same as lines
# 2067-2116) and " escaped using Find and replace 


rfuToResType <- "If .FU(i) = \"PRDOM\" Then
									.class2(i) = 2	'CON
								End If
								If .FU(i) = \"PWDOM\" Then
									.class2(i) = 2	'CON
								End If
								If .FU(i) = \"SBMX1\" Then
									.class2(i) = 2	'CON
								End If
								If .FU(i) = \"SBLOW\" Then
									.class2(i) = 2	'CON
								End If
								If .FU(i) = \"PJDOM\" Then
									.class2(i) = 2	'CON
								End If
								If .FU(i) = \"SBDOM\" Then
									.class2(i) = 2	'CON
								End If
								If .FU(i) = \"OCLOW\" Then
									.class2(i) = 2	'CON
								End If
								If .FU(i) = \"BFDOM\" Then
									.class2(i) = 2	'CON
								End If
								If .FU(i) = \"PJMX1\" Then
									.class2(i) = 2	'CON
								End If
								If .FU(i) = \"CONMX\" Then
									.class2(i) = 3	'MIX
								End If
								If .FU(i) = \"HRDMX\" Then
									.class2(i) = 3	'MIX
								End If
								If .FU(i) = \"PRWMX\" Then
									.class2(i) = 4	'MIX
								End If
								If .FU(i) = \"HRDOM\" Then
									.class2(i) = 4	'DEC
								End If
								If .FU(i) = \"PODOM\" Then
									.class2(i) = 4	'DEC
								End If
								If .FU(i) = \"BWDOM\" Then
									.class2(i) = 4	'DEC
								End If
								If .FU(i) = \"PODOM\" Then
									.class2(i) = 4	'DEC
								End If
								If .FU(i) = \"OTHHD\" Then
									.class2(i) = 4	'DEC" %>% 
  stringr::str_squish() %>% 
  stringr::str_extract_all("[[:upper:]]{3,}\\d?", simplify = TRUE) %>% 
  matrix(ncol = 2, byrow = TRUE) %>% data.frame(stringsAsFactors = FALSE) %>% 
  setNames(c("RegionalForestUnit", "ResourceType")) %>% 
  distinct()

usethis::use_data(rfuToResType, overwrite = TRUE)
