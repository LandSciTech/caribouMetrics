## code to prepare `plcToResType` dataset goes here

# Transcribed from lines 2308 - 2319 and 
# 3378-3393 of the LSL script Create_Caribou_RSPF_Grid.lsl

# plcToResType <- read.csv("inputTables/BorealCaribouPLCtoResType.csv", 
#                          stringsAsFactors = FALSE)

# could also construct from lsl directly
# lines 2308 - 2319 
plcToResType <- ".OWAT_PLC(j) = g.class.1(j)
				.TWAT_PLC(j) = g.class.2(j)
				.ST_PLC(j) = g.class.10(j)
				.DEC_PLC(j) = g.class.11(j)
				.MIX_PLC(j) = g.class.12(j)
				.CON_PLC(j) = g.class.13(j)
				.CSWP_PLC(j)= g.class.19(j)
				.OFEN_PLC(j) = g.class.20(j)
				.TFEN_PLC(j) = g.class.21(j)
				.OBOG_PLC(j) = g.class.22(j)
				.TBOG_PLC(j) = g.class.23(j)" %>% 
  stringr::str_squish() %>% stringr::str_remove_all("_PLC") %>% 
  stringr::str_extract_all("[[:upper:]]{2,}|\\d{1,2}", simplify = TRUE) %>% 
  matrix(ncol = 2, byrow = TRUE) %>% as.data.frame(stringsAsFactors = FALSE) %>% 
  setNames(c("ResourceType", "PLCCode")) %>% 
  mutate(ResourceType = case_when(ResourceType %in% c("OFEN", "OBOG") ~ "LGOP",
                                  ResourceType %in% c("CSWP", "TBOG", "TFEN") ~ "LGTP",
                                  ResourceType %in% c("OWAT", "TWAT") ~ "LGW",
                                  TRUE ~ ResourceType),
         PLCCode = as.numeric(PLCCode)) %>% 
  bind_rows(data.frame(ResourceType = "other", 
                       PLCCode = which(!1:29 %in% .$PLCCode), 
                       stringsAsFactors = FALSE)) 


# lines 3378
".LGW_PLC(j) = (.OWAT_PLC(j) + .TWAT_PLC(j)) / gsarea
						.OWAT_PLC(j) = .OWAT_PLC(j) / gsarea
						.TWAT_PLC(j) = .TWAT_PLC(j) / gsarea
						.LGOP_PLC(j) = (.OFEN_PLC(j) + .OBOG_PLC(j)) / gsarea
						.LGTP_PLC(j) = (.CSWP_PLC(j) + .TBOG_PLC(j) + .TFEN_PLC(j)) / gsarea
						.CSWP_PLC(j) = .CSWP_PLC(j) / gsarea
						.OFEN_PLC(j) = .OFEN_PLC(j) / gsarea 
						.TFEN_PLC(j) = .TFEN_PLC(j) / gsarea 
						.OBOG_PLC(j) = .OBOG_PLC(j) / gsarea 
						.TBOG_PLC(j) = .TBOG_PLC(j) / gsarea 
						.ST_PLC(j) = .ST_PLC(j) / gsarea
						.DEC_PLC(j) = .DEC_PLC(j) / gsarea 
						.MIX_PLC(j) = .MIX_PLC(j) / gsarea 
						.CON_PLC(j) = .CON_PLC(j) / gsarea"

usethis::use_data(plcToResType, overwrite = TRUE)
