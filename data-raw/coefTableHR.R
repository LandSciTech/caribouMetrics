## code to prepare `coefTableHR` dataset goes here
library(tidyverse)

# Make table from LSL script and record window size #===========================
#_S5 = 5000 and _S6 = 10000

# copied from lines 4962 - 6266 of Caribou_Range_Specific_RSPF.lsl
# edited for Pagwachuan Fall because DEC was missing but was in tables in report/paper
inText <- {"#region Model code for Berens Caribou Range
					If RR1 = 1 Then
						If gs1.S_area(j) >= Threshold Then	
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-64.3391 * gs.DEC_S6(j)) + _
							(-17.7959 * gs.LGTP_S6(j)) + _
							(-8.87683 * gs.MIX_S6(j)) + _
							(-0.0231416 * gs.TDENLF_S6(j)) + _
							(0.117523 * gs.ESK_S6(j)) + _
							(1.19193 * gs.DTN_S6(j)) + _
							(2.45984 * gs.ST_S6(j)) + _
							(3.90644 * gs.LGW_S6(j)) + _
							(8.04834 * gs.CON_S6(j)) + _
							(48.407 * gs.LGOP_S6(j)) + (-4.11371)
							
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > BESp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-80.9814 * gs.DEC_S6(j)) + _
							(-15.1026 * gs.MIX_S6(j)) + _
							(-6.73918 * gs.LGTP_S6(j)) + _
							(-1.2562 * gs.ST_S6(j)) + _
							(-0.611864 * gs.DTN_S6(j)) + _
							(-0.463687 * gs.ESK_S6(j)) + _
							(-0.04311 * gs.TDENLF_S6(j)) + _
							(4.07661 * gs.LGW_S6(j)) + _
							(5.85445 * gs.CON_S6(j)) + _
							(45.3663 * gs.LGOP_S6(j)) + (-2.39214)

							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > BESu Then
								.cf_SU_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _
							(-16.814 * gs.DEC_S6(j)) + _
							(-5.48341 * gs.ST_S6(j)) + _
							(-4.2601 * gs.DTN_S6(j)) + _
							(-4.10281 * gs.MIX_S6(j)) + _
							(-0.115161 * gs.TDENLF_S6(j)) + _
							(0.428106 * gs.ESK_S6(j)) + _
							(1.23197 * gs.CON_S6(j)) + _
							(1.48339 * gs.LGW_S6(j)) + _
							(8.05471 * gs.LGTP_S6(j)) + _
							(19.703 * gs.LGOP_S6(j)) + (-0.932354)
							 							 
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > BEFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _	
							(-16.0938 * gs.ST_S6(j)) + _
							(-14.7214 * gs.DTN_S6(j)) + _
							(-11.4149 * gs.DEC_S6(j)) + _
							(-11.1377 * gs.LGW_S6(j)) + _
							(-8.63904 * gs.CON_S6(j)) + _
							(-8.54187 * gs.MIX_S6(j)) + _
							(-3.41591 * gs.ESK_S6(j)) + _
							(-0.898634 * gs.TDENLF_S6(j)) + _
							(1.37538 * gs.LGOP_S6(j)) + _
							(1.38735 * gs.LGTP_S6(j)) + (8.68769)

							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > BEWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If
					#endregion
					
					#region Model code for Sydney Caribou Range
					If RR2 = 1 Then
						If gs2.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-64.3391 * gs.DEC_S6(j)) + _
							(-17.7959 * gs.LGTP_S6(j)) + _
							(-8.87683 * gs.MIX_S6(j)) + _
							(-0.0231416 * gs.TDENLF_S6(j)) + _
							(0.117523 * gs.ESK_S6(j)) + _
							(1.19193 * gs.DTN_S6(j)) + _
							(2.45984 * gs.ST_S6(j)) + _
							(3.90644 * gs.LGW_S6(j)) + _
							(8.04834 * gs.CON_S6(j)) + _
							(48.407 * gs.LGOP_S6(j)) + (-4.11371)
							
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > SYSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-80.9814 * gs.DEC_S6(j)) + _
							(-15.1026 * gs.MIX_S6(j)) + _
							(-6.73918 * gs.LGTP_S6(j)) + _
							(-1.2562 * gs.ST_S6(j)) + _
							(-0.611864 * gs.DTN_S6(j)) + _
							(-0.463687 * gs.ESK_S6(j)) + _
							(-0.04311 * gs.TDENLF_S6(j)) + _
							(4.07661 * gs.LGW_S6(j)) + _
							(5.85445 * gs.CON_S6(j)) + _
							(45.3663 * gs.LGOP_S6(j)) + (-2.39214)
															 
							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > SYSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-16.814 * gs.DEC_S6(j)) + _
							(-5.48341 * gs.ST_S6(j)) + _
							(-4.2601 * gs.DTN_S6(j)) + _
							(-4.10281 * gs.MIX_S6(j)) + _
							(-0.115161 * gs.TDENLF_S6(j)) + _
							(0.428106 * gs.ESK_S6(j)) + _
							(1.23197 * gs.CON_S6(j)) + _
							(1.48339 * gs.LGW_S6(j)) + _
							(8.05471 * gs.LGTP_S6(j)) + _
							(19.703 * gs.LGOP_S6(j)) + (-0.932354)
							 
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > SYFa Then
								.cf_FA_Use(j) = 1
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-16.0938 * gs.ST_S6(j)) + _
							(-14.7214 * gs.DTN_S6(j)) + _
							(-11.4149 * gs.DEC_S6(j)) + _
							(-11.1377 * gs.LGW_S6(j)) + _
							(-8.63904 * gs.CON_S6(j)) + _
							(-8.54187 * gs.MIX_S6(j)) + _
							(-3.41591 * gs.ESK_S6(j)) + _
							(-0.898634 * gs.TDENLF_S6(j)) + _
							(1.37538 * gs.LGOP_S6(j)) + _
							(1.38735 * gs.LGTP_S6(j)) + (8.68769)
							 
							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > SYWi Then
								.cf_WI_Use(j) = 1
							End If
							#endregion
						End If
					End If
					#endregion

					#region Model code for Churchill Caribou Range
					If RR3 = 1 Then
						If gs3.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-13.8214 * gs.DEC_S5(j)) + _
							(-1.78556 * gs.MIX_S5(j)) + _
							(-0.55255 * gs.ESK_S5(j)) + _
							(0.0160082 * gs.TDENLF_S5(j)) + _
							(2.32037 * gs.DTN_S5(j)) + _
							(4.02012 * gs.ST_S5(j)) + _
							(4.42165 * gs.LGW_S5(j)) + _
							(5.34232 * gs.CON_S5(j)) + _
							(7.17763 * gs.LGTP_S5(j)) + _
							(16.1389 * gs.LGOP_S5(j)) + (-5.01081)
							 
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > CHSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-41.8592 * gs.DEC_S5(j)) + _
							(-0.454738 * gs.ESK_S5(j)) + _
							(-0.0350838 * gs.TDENLF_S5(j)) + _
							(0.530534 * gs.ST_S5(j)) + _
							(0.802784 * gs.MIX_S5(j)) + _
							(2.15437 * gs.DTN_S5(j)) + _
							(2.77072 * gs.LGW_S5(j)) + _
							(4.08734 * gs.CON_S5(j)) + _
							(5.32027 * gs.LGTP_S5(j)) + _
							(44.061 * gs.LGOP_S5(j)) + (-3.95866)
 							 							
							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > CHSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-56.4116 * gs.DEC_S5(j)) + _
							(-5.23716 * gs.MIX_S5(j)) + _
							(-0.356225 * gs.DTN_S5(j)) + _
							(-0.128353 * gs.TDENLF_S5(j)) + _
							(2.86515 * gs.LGTP_S5(j)) + _
							(4.40582 * gs.CON_S5(j)) + (-2.15043)
							 							 						 
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > CHFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-31.845 * gs.DEC_S5(j)) + _
							(-10.6597 * gs.MIX_S5(j)) + _
							(-0.708967 * gs.DTN_S5(j)) + _
							(-0.275126 * gs.ESK_S5(j)) + _
							(0.0332388 * gs.TDENLF_S5(j)) + _
							(0.346165 * gs.LGTP_S5(j)) + _
							(1.80009 * gs.LGW_S5(j)) + _
							(6.18843 * gs.CON_S5(j)) + _
							(6.26025 * gs.ST_S5(j)) + _
							(32.0258 * gs.LGOP_S5(j)) + (-3.99216)

							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > CHWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If
					#endregion

					#region Model code for Brightsand Caribou Range
					If RR4 = 1 Then
						If gs4.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-13.8214 * gs.DEC_S5(j)) + _
							(-1.78556 * gs.MIX_S5(j)) + _
							(-0.55255 * gs.ESK_S5(j)) + _
							(0.0160082 * gs.TDENLF_S5(j)) + _
							(2.32037 * gs.DTN_S5(j)) + _
							(4.02012 * gs.ST_S5(j)) + _
							(4.42165 * gs.LGW_S5(j)) + _
							(5.34232 * gs.CON_S5(j)) + _
							(7.17763 * gs.LGTP_S5(j)) + _
							(16.1389 * gs.LGOP_S5(j)) + (-5.01081)
							
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > BSSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-41.8592 * gs.DEC_S5(j)) + _
							(-0.454738 * gs.ESK_S5(j)) + _
							(-0.0350838 * gs.TDENLF_S5(j)) + _
							(0.530534 * gs.ST_S5(j)) + _
							(0.802784 * gs.MIX_S5(j)) + _
							(2.15437 * gs.DTN_S5(j)) + _
							(2.77072 * gs.LGW_S5(j)) + _
							(4.08734 * gs.CON_S5(j)) + _
							(5.32027 * gs.LGTP_S5(j)) + _
							(44.061 * gs.LGOP_S5(j)) + (-3.95866)
							 							
							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > BSSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-56.4116 * gs.DEC_S5(j)) + _
							(-5.23716 * gs.MIX_S5(j)) + _
							(-0.356225 * gs.DTN_S5(j)) + _
							(-0.128353 * gs.TDENLF_S5(j)) + _
							(2.86515 * gs.LGTP_S5(j)) + _
							(4.40582 * gs.CON_S5(j)) + (-2.15043)
							 							 
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > BSFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-34.5627 * gs.DEC_S5(j)) + _
							(-11.2342 * gs.MIX_S5(j)) + _
							(-1.57714 * gs.DTN_S5(j)) + _
							(-0.421174 * gs.LGTP_S5(j)) + _
							(-0.282713 * gs.ESK_S5(j)) + _
							(0.00354971 * gs.TDENLF_S5(j)) + _
							(0.819689 * gs.LGW_S5(j)) + _
							(5.24557 * gs.ST_S5(j)) + _
							(5.33736 * gs.CON_S5(j)) + _
							(30.3627 * gs.LGOP_S5(j)) + (-3.1063)			
														 
							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > BSWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						Else
							.cf_SP_Use(j) = MISSING
						End If
					End If
					#endregion
					
					#region Model code for Nipigon Caribou Range
					If RR5 = 1 Then
						If gs5.S_area(j) >= Threshold Then
							#region Code for Spring Use Model							
							logitSP_Use = _ 
							(-44.8817 * gs.DEC_S5(j)) + _
							(-1.61879 * gs.ST_S5(j)) + _
							(-0.0591871 * gs.TDENLF_S5(j)) + _
							(0.445199 * gs.ESK_S5(j)) + _
							(0.63021 * gs.CON_S5(j)) + _
							(0.775091 * gs.DTN_S5(j)) + _
							(2.33635 * gs.LGW_S5(j)) + _
							(7.39955 * gs.MIX_S5(j)) + _
							(7.82448 * gs.LGTP_S5(j)) + _
							(20.2548 * gs.LGOP_S5(j)) + (-3.19835)
							 
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > NPSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-24.0817 * gs.DEC_S5(j)) + _
							(-0.902747 * gs.CON_S5(j)) + _
							(-0.120702 * gs.LGW_S5(j)) + _
							(-0.108422 * gs.TDENLF_S5(j)) + _
							(0.724132 * gs.ESK_S5(j)) + _
							(1.1134 * gs.DTN_S5(j)) + _
							(2.76126 * gs.ST_S5(j)) + _
							(2.87828 * gs.MIX_S5(j)) + _
							(4.37304 * gs.LGTP_S5(j)) + _
							(26.2326 * gs.LGOP_S5(j)) + (-2.46477)

							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > NPSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-38.5131 * gs.DEC_S5(j)) + _
							(-0.0657281 * gs.TDENLF_S5(j)) + _
							(0.351276 * gs.ESK_S5(j)) + _
							(0.61502 * gs.CON_S5(j)) + _
							(1.10202 * gs.ST_S5(j)) + _
							(1.2464 * gs.LGW_S5(j)) + _
							(1.44648 * gs.DTN_S5(j)) + _
							(4.35131 * gs.MIX_S5(j)) + _
							(7.98351 * gs.LGTP_S5(j)) + _
							(27.1499 * gs.LGOP_S5(j)) + (-2.81546)
							
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > NPFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-2.84384 * gs.DTN_S5(j)) + _
							(-2.10533 * gs.DEC_S5(j)) + _
							(-0.969895 * gs.LGW_S5(j)) + _
							(-0.428307 * gs.ESK_S5(j)) + _
							(-0.107917 * gs.TDENLF_S5(j)) + _
							(0.960811 * gs.LGTP_S5(j)) + _
							(1.16347 * gs.CON_S5(j)) + _
							(2.07816 * gs.MIX_S5(j)) + _
							(3.26042 * gs.ST_S5(j)) + _
							(6.35273 * gs.LGOP_S5(j)) + (-1.81318)
							 
							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > NPWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If				
					#endregion
					
					#region Model code for Pagwachuan Caribou Range
					If RR6 = 1 Then
						If gs6.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-42.0006 * gs.DTN_S6(j)) + _
							(-37.6389 * gs.DEC_S6(j)) + _
							(-4.58361 * gs.MIX_S6(j)) + _
							(-4.5621 * gs.ST_S6(j)) + _
							(-2.7476 * gs.LGTP_S6(j)) + _
							(-0.35184 * gs.TDENLF_S6(j)) + _
							(1.40204 * gs.LGW_S6(j)) + _
							(1.6818 * gs.CON_S6(j)) + _
							(3.55455 * gs.ESK_S6(j)) + _
							(4.39236 * gs.LGOP_S6(j)) + (-0.15753) 
							  								 
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > PWSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-104.387 * gs.DEC_S6(j)) + _
							(-41.3792 * gs.DTN_S6(j)) + _
							(-2.01557 * gs.ST_S6(j)) + _
							(-1.83914 * gs.MIX_S6(j)) + _
							(-0.123711 * gs.TDENLF_S6(j)) + _
							(0.360641 * gs.LGTP_S6(j)) + _
							(2.24184 * gs.ESK_S6(j)) + _
							(5.30014 * gs.CON_S6(j)) + _
							(7.37853 * gs.LGW_S6(j)) + _
							(7.53761 * gs.LGOP_S6(j)) + (-3.73633)
							 
							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > PWSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _
              (-39.0287 * gs.DEC_S6(j)) + _
							(-34.2672 * gs.DTN_S6(j)) + _
							(-14.1546 * gs.MIX_S6(j)) + _
							(-13.7814 * gs.ST_S6(j)) + _
							(-6.47953 * gs.CON_S6(j)) + _
							(-6.26648 * gs.LGTP_S6(j)) + _
							(-6.11699 * gs.LGW_S6(j)) + _
							(-3.48516 * gs.LGOP_S6(j)) + _
							(-0.738744 * gs.TDENLF_S6(j)) + _
							(2.86622 * gs.ESK_S6(j)) + (6.50218)

							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > PWFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _ 
							(-3.48497 * gs.DEC_S6(j)) + _
							(-0.59024 * gs.TDENLF_S6(j)) + _
							(3.51487 * gs.MIX_S6(j)) + _
							(5.2854 * gs.ESK_S6(j)) + _
							(9.68799 * gs.CON_S6(j)) + _
							(10.5729 * gs.LGOP_S6(j)) + _
							(10.9739 * gs.LGTP_S6(j)) + _
							(15.4934 * gs.ST_S6(j)) + _
							(16.6945 * gs.LGW_S6(j)) + _
							(23.8352 * gs.DTN_S6(j)) + (-11.7408)
							 							 							 
							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > PWWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If				
					#endregion
					
					#region Model code for Kesagami Caribou Range
					If RR7 = 1 Then
						If gs7.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-158.952 * gs.DEC_S6(j)) + _
							(-11.6382 * gs.MIX_S6(j)) + _
							(0.199959 * gs.TDENLF_S6(j)) + _
							(0.623459 * gs.LGTP_S6(j)) + _
							(1.7108 * gs.ESK_S6(j)) + _
							(7.48504 * gs.CON_S6(j)) + _
							(9.0546 * gs.DTN_S6(j)) + _
							(10.4209 * gs.LGW_S6(j)) + _
							(12.4208 * gs.LGOP_S6(j)) + _
							(17.5233 * gs.ST_S6(j)) + (-7.62303)
							 
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > KSSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-101.209 * gs.DEC_S6(j)) + _
							(-15.7408 * gs.MIX_S6(j)) + _
							(-3.4 * gs.LGTP_S6(j)) + _
							(-0.603098 * gs.DTN_S6(j)) + _
							(0.0610013 * gs.TDENLF_S6(j)) + _
							(0.625172 * gs.ESK_S6(j)) + _
							(2.57208 * gs.CON_S6(j)) + _
							(6.69949 * gs.LGW_S6(j)) + _
							(7.49653 * gs.LGOP_S6(j)) + _
							(14.6237 * gs.ST_S6(j)) + (-3.16522)
							 
							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > KSSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-67.3413 * gs.DEC_S6(j)) + _
							(-26.0146 * gs.MIX_S6(j)) + _
							(-7.84446 * gs.DTN_S6(j)) + _
							(-4.4653 * gs.LGOP_S6(j)) + _
							(-1.06279 * gs.LGTP_S6(j)) + _
							(-0.0620399 * gs.TDENLF_S6(j)) + _
							(2.24476 * gs.ESK_S6(j)) + _
							(2.88131 * gs.CON_S6(j)) + _
							(8.13394 * gs.LGW_S6(j)) + _
							(14.278 * gs.ST_S6(j)) + (-2.68833)
							
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > KSFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _ 
							(-195.376 * gs.DEC_S6(j)) + _
							(-11.2609 * gs.MIX_S6(j)) + _
							(-5.08062 * gs.DTN_S6(j)) + _
							(-1.281 * gs.LGTP_S6(j)) + _
							(-0.929046 * gs.ESK_S6(j)) + _
							(0.171649 * gs.TDENLF_S6(j)) + _
							(3.61196 * gs.LGW_S6(j)) + _
							(6.54676 * gs.LGOP_S6(j)) + _
							(7.00059 * gs.CON_S6(j)) + _
							(14.7874 * gs.ST_S6(j)) + (-5.66009)
							 
							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > KSWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If				
					#endregion
					
					#region Model code for Lake Superior Coast Caribou Range
					If RR8 = 1 Then
						If gs8.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _
							0.0
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > LSSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							0.0
							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > LSSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							0.0
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > LSFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							0.0
							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > LSWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If
					#endregion
									
					#region Model code for Discontinuous Distribution Caribou Range
					If RR9 = 1 Then
						If gs9.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _
							0.0
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > DDSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							0.0
							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > DDSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							0.0
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > DDFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model							
							logitWI_Use = _
							0.0
							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > DDWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If
					#endregion
					
					#region Model code for Swan Caribou Range
					If RR10 = 1 Then
						If gs10.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-44.2195 * gs.DEC_S6(j)) + _
							(-9.06067 * gs.MIX_S6(j)) + _
							(-4.24696 * gs.LGOP_S6(j)) + _
							(-0.441428 * gs.TDENLF_S6(j)) + _
							(-0.127303 * gs.ESK_S6(j)) + _
							(0.0817173 * gs.DTN_S6(j)) + _
							(2.1867 * gs.ST_S6(j)) + _
							(3.07141 * gs.CON_S6(j)) + _
							(6.74145 * gs.LGW_S6(j)) + _
							(13.2148 * gs.LGTP_S6(j)) + (-4.40276)
							
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > SWSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-19.569 * gs.DEC_S6(j)) + _
							(-2.41043 * gs.MIX_S6(j)) + _
							(-0.782122 * gs.ESK_S6(j)) + _
							(0.00946266 * gs.TDENLF_S6(j)) + _
							(4.26762 * gs.CON_S6(j)) + _
							(4.38107 * gs.ST_S6(j)) + _
							(5.93138 * gs.DTN_S6(j)) + _
							(6.26119 * gs.LGW_S6(j)) + _
							(6.91973 * gs.LGTP_S6(j)) + _
							(14.7378 * gs.LGOP_S6(j)) + (-6.50782)

							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > SWSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-23.3775 * gs.DEC_S6(j)) + _
							(-17.6436 * gs.MIX_S6(j)) + _
							(-10.1348 * gs.LGOP_S6(j)) + _
							(-9.46093 * gs.DTN_S6(j)) + _
							(-3.5667 * gs.CON_S6(j)) + _
							(-2.96147 * gs.ST_S6(j)) + _
							(-0.363172 * gs.TDENLF_S6(j)) + _
							(0.21628 * gs.ESK_S6(j)) + _
							(4.35684 * gs.LGW_S6(j)) + _
							(6.72984 * gs.LGTP_S6(j)) + (1.48626)
							 
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > SWFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-50.0395 * gs.DEC_S6(j)) + _
							(-0.376008 * gs.TDENLF_S6(j)) + _
							(-0.246234 * gs.ESK_S6(j)) + _
							(-0.109764 * gs.DTN_S6(j)) + _
							(1.38047 * gs.MIX_S6(j)) + _
							(1.72236 * gs.LGW_S6(j)) + _
							(3.13675 * gs.CON_S6(j)) + _
							(3.75243 * gs.LGTP_S6(j)) + _
							(4.67817 * gs.ST_S6(j)) + _
							(17.4244 * gs.LGOP_S6(j)) + (-4.48403)

							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > SWWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If				
					#endregion
					
					#region Model code for Spirit Caribou Range
					If RR11 = 1 Then
						If gs11.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-44.2195 * gs.DEC_S6(j)) + _
							(-9.06067 * gs.MIX_S6(j)) + _
							(-4.24696 * gs.LGOP_S6(j)) + _
							(-0.441428 * gs.TDENLF_S6(j)) + _
							(-0.127303 * gs.ESK_S6(j)) + _
							(0.0817173 * gs.DTN_S6(j)) + _
							(2.1867 * gs.ST_S6(j)) + _
							(3.07141 * gs.CON_S6(j)) + _
							(6.74145 * gs.LGW_S6(j)) + _
							(13.2148 * gs.LGTP_S6(j)) + (-4.40276)
							
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > SPSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-19.569 * gs.DEC_S6(j)) + _
							(-2.41043 * gs.MIX_S6(j)) + _
							(-0.782122 * gs.ESK_S6(j)) + _
							(0.00946266 * gs.TDENLF_S6(j)) + _
							(4.26762 * gs.CON_S6(j)) + _
							(4.38107 * gs.ST_S6(j)) + _
							(5.93138 * gs.DTN_S6(j)) + _
							(6.26119 * gs.LGW_S6(j)) + _
							(6.91973 * gs.LGTP_S6(j)) + _
							(14.7378 * gs.LGOP_S6(j)) + (-6.50782)

							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > SPSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-23.3775 * gs.DEC_S6(j)) + _
							(-17.6436 * gs.MIX_S6(j)) + _
							(-10.1348 * gs.LGOP_S6(j)) + _
							(-9.46093 * gs.DTN_S6(j)) + _
							(-3.5667 * gs.CON_S6(j)) + _
							(-2.96147 * gs.ST_S6(j)) + _
							(-0.363172 * gs.TDENLF_S6(j)) + _
							(0.21628 * gs.ESK_S6(j)) + _
							(4.35684 * gs.LGW_S6(j)) + _
							(6.72984 * gs.LGTP_S6(j)) + (1.48626)
							 
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > SPFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-50.0395 * gs.DEC_S6(j)) + _
							(-0.376008 * gs.TDENLF_S6(j)) + _
							(-0.246234 * gs.ESK_S6(j)) + _
							(-0.109764 * gs.DTN_S6(j)) + _
							(1.38047 * gs.MIX_S6(j)) + _
							(1.72236 * gs.LGW_S6(j)) + _
							(3.13675 * gs.CON_S6(j)) + _
							(3.75243 * gs.LGTP_S6(j)) + _
							(4.67817 * gs.ST_S6(j)) + _
							(17.4244 * gs.LGOP_S6(j)) + (-4.48403)

							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > SPWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If				
					#endregion
					
					#region Model code for Kinloch Caribou Range
					If RR12 = 1 Then
						If gs12.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-44.2195 * gs.DEC_S6(j)) + _
							(-9.06067 * gs.MIX_S6(j)) + _
							(-4.24696 * gs.LGOP_S6(j)) + _
							(-0.441428 * gs.TDENLF_S6(j)) + _
							(-0.127303 * gs.ESK_S6(j)) + _
							(0.0817173 * gs.DTN_S6(j)) + _
							(2.1867 * gs.ST_S6(j)) + _
							(3.07141 * gs.CON_S6(j)) + _
							(6.74145 * gs.LGW_S6(j)) + _
							(13.2148 * gs.LGTP_S6(j)) + (-4.40276)
							
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > KLSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-2.22557 * gs.LGOP_S6(j)) + _
							(-2.13973 * gs.DEC_S6(j)) + _
							(-1.08059 * gs.MIX_S6(j)) + _
							(-0.811039 * gs.ST_S6(j)) + _
							(-0.329841 * gs.TDENLF_S6(j)) + _
							(0.12532 * gs.ESK_S6(j)) + _
							(0.146291 * gs.DTN_S6(j)) + _
							(2.03843 * gs.CON_S6(j)) + _
							(7.24798 * gs.LGW_S6(j)) + _
							(14.3959 * gs.LGTP_S6(j)) + (-4.44855)
							
							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > KLSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-23.3775 * gs.DEC_S6(j)) + _
							(-17.6436 * gs.MIX_S6(j)) + _
							(-10.1348 * gs.LGOP_S6(j)) + _
							(-9.46093 * gs.DTN_S6(j)) + _
							(-3.5667 * gs.CON_S6(j)) + _
							(-2.96147 * gs.ST_S6(j)) + _
							(-0.363172 * gs.TDENLF_S6(j)) + _
							(0.21628 * gs.ESK_S6(j)) + _
							(4.35684 * gs.LGW_S6(j)) + _
							(6.72984 * gs.LGTP_S6(j)) + (1.48626)
							 
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > KLFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-50.0395 * gs.DEC_S6(j)) + _
							(-0.376008 * gs.TDENLF_S6(j)) + _
							(-0.246234 * gs.ESK_S6(j)) + _
							(-0.109764 * gs.DTN_S6(j)) + _
							(1.38047 * gs.MIX_S6(j)) + _
							(1.72236 * gs.LGW_S6(j)) + _
							(3.13675 * gs.CON_S6(j)) + _
							(3.75243 * gs.LGTP_S6(j)) + _
							(4.67817 * gs.ST_S6(j)) + _
							(17.4244 * gs.LGOP_S6(j)) + (-4.48403)

							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > KLWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If				
					#endregion
					
					#region Model code for Ozhiski Caribou Range
					If RR13 = 1 Then
						If gs13.S_area(j) >= Threshold Then	
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-44.2195 * gs.DEC_S6(j)) + _
							(-9.06067 * gs.MIX_S6(j)) + _
							(-4.24696 * gs.LGOP_S6(j)) + _
							(-0.441428 * gs.TDENLF_S6(j)) + _
							(-0.127303 * gs.ESK_S6(j)) + _
							(0.0817173 * gs.DTN_S6(j)) + _
							(2.1867 * gs.ST_S6(j)) + _
							(3.07141 * gs.CON_S6(j)) + _
							(6.74145 * gs.LGW_S6(j)) + _
							(13.2148 * gs.LGTP_S6(j)) + (-4.40276)
							
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > OZSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-19.569 * gs.DEC_S6(j)) + _
							(-2.41043 * gs.MIX_S6(j)) + _
							(-0.782122 * gs.ESK_S6(j)) + _
							(0.00946266 * gs.TDENLF_S6(j)) + _
							(4.26762 * gs.CON_S6(j)) + _
							(4.38107 * gs.ST_S6(j)) + _
							(5.93138 * gs.DTN_S6(j)) + _
							(6.26119 * gs.LGW_S6(j)) + _
							(6.91973 * gs.LGTP_S6(j)) + _
							(14.7378 * gs.LGOP_S6(j)) + (-6.50782)

							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > OZSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-23.3775 * gs.DEC_S6(j)) + _
							(-17.6436 * gs.MIX_S6(j)) + _
							(-10.1348 * gs.LGOP_S6(j)) + _
							(-9.46093 * gs.DTN_S6(j)) + _
							(-3.5667 * gs.CON_S6(j)) + _
							(-2.96147 * gs.ST_S6(j)) + _
							(-0.363172 * gs.TDENLF_S6(j)) + _
							(0.21628 * gs.ESK_S6(j)) + _
							(4.35684 * gs.LGW_S6(j)) + _
							(6.72984 * gs.LGTP_S6(j)) + (1.48626)
							 
							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > OZFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-50.0395 * gs.DEC_S6(j)) + _
							(-0.376008 * gs.TDENLF_S6(j)) + _
							(-0.246234 * gs.ESK_S6(j)) + _
							(-0.109764 * gs.DTN_S6(j)) + _
							(1.38047 * gs.MIX_S6(j)) + _
							(1.72236 * gs.LGW_S6(j)) + _
							(3.13675 * gs.CON_S6(j)) + _
							(3.75243 * gs.LGTP_S6(j)) + _
							(4.67817 * gs.ST_S6(j)) + _
							(17.4244 * gs.LGOP_S6(j)) + (-4.48403)

							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > OZWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If				
					#endregion
					
					#region Model code for Missisa Caribou Range
					If RR14 = 1 Then
						If gs14.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-47.0453 * gs.TDENLF_S6(j)) + _
							(-7.00107 * gs.DTN_S6(j)) + _
							(-3.20956 * gs.ST_S6(j)) + _
							(-2.00532 * gs.LGTP_S6(j)) + _
							(-1.26999 * gs.CON_S6(j)) + _
							(-0.252612 * gs.LGOP_S6(j)) + _
							(0.690524 * gs.LGW_S6(j)) + _
							(1.65949 * gs.LGMD_S6(j)) + (-0.223601)
							 
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > MSSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-54.014 * gs.TDENLF_S6(j)) + _
							(-9.0874 * gs.LGMD_S6(j)) + _
							(-1.09398 * gs.DTN_S6(j)) + _
							(6.03112 * gs.ST_S6(j)) + _
							(7.91777 * gs.LGOP_S6(j)) + _
							(9.73417 * gs.LGW_S6(j)) + _
							(12.1113 * gs.CON_S6(j)) + _
							(12.3221 * gs.LGTP_S6(j)) + (-10.9329)

							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > MSSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-11.7367 * gs.LGMD_S6(j)) + _
							(-3.85024 * gs.DTN_S6(j)) + _
							(-1.1584 * gs.TDENLF_S6(j)) + _
							(1.31976 * gs.LGTP_S6(j)) + _
							(3.02311 * gs.CON_S6(j)) + _
							(4.23221 * gs.LGW_S6(j)) + _
							(4.23322 * gs.ST_S6(j)) + _
							(4.61501 * gs.LGOP_S6(j)) + (-5.98728)

							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > MSFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-32.4414 * gs.TDENLF_S6(j)) + _
							(-4.96709 * gs.LGMD_S6(j)) + _
							(-2.54278 * gs.LGTP_S6(j)) + _
							(-1.55981 * gs.LGW_S6(j)) + _
							(1.38101 * gs.CON_S6(j)) + _
							(4.07806 * gs.DTN_S6(j)) + _
							(10.5872 * gs.ST_S6(j)) + _
							(10.8757 * gs.LGOP_S6(j)) + (-5.28692)
							 									 
							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > MSWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End if
					End If
					#endregion
					
					#region Model code for James Bay Caribou Range
					If RR15 = 1 Then
						If gs15.S_area(j) >= Threshold Then
							#region Code for Spring Use Model
							logitSP_Use = _ 
							(-8.12364 * gs.LGMD_S6(j)) + _
							(-6.76679 * gs.LGOP_S6(j)) + _
							(-5.14577 * gs.CON_S6(j)) + _
							(-0.871053 * gs.LGW_S6(j)) + _
							(-0.0389879 * gs.TDENLF_S6(j)) + _
							(0.783623 * gs.DTN_S6(j)) + _
							(2.6777 * gs.ST_S6(j)) + _
							(3.41592 * gs.LGTP_S6(j)) + (-1.79609)
							 
							 odds_SP = 2.718281828 ^ logitSP_Use
							.pSP_Use(j) = odds_SP / (1.0 + odds_SP)
							.cf_SP_Use(j) = 0 
							 
							If .pSP_Use (j) > JBSp Then
								.cf_SP_Use(j) = 1 
							End If
							#endregion
							
							#region Code for Summer Use Model
							logitSU_Use = _ 
							(-14.0371 * gs.LGMD_S6(j)) + _
							(-6.2329 * gs.LGOP_S6(j)) + _
							(-4.58141 * gs.CON_S6(j)) + _
							(-1.01115 * gs.LGW_S6(j)) + _
							(-0.917044 * gs.DTN_S6(j)) + _
							(-0.0964143 * gs.TDENLF_S6(j)) + _
							(1.2406 * gs.ST_S6(j)) + _
							(2.80308 * gs.LGTP_S6(j)) + (-1.09174)
							 
							odds_SU = 2.718281828 ^ logitSU_Use
							.pSU_Use(j) = odds_SU / (1.0 + odds_SU)
							.cf_SU_Use(j) = 0 
							 
							If .pSU_Use (j) > JBSu Then
								.cf_SU_Use(j) = 1 
							End If 
							#endregion
							
							#region Code for Fall Use Model
							logitFA_Use = _ 
							(-18.8983 * gs.DTN_S6(j)) + _
							(-15.4255 * gs.LGMD_S6(j)) + _
							(-15.0418 * gs.CON_S6(j)) + _
							(-11.1888 * gs.LGW_S6(j)) + _
							(-9.69171 * gs.TDENLF_S6(j)) + _
							(-9.0392 * gs.LGOP_S6(j)) + _
							(-6.80688 * gs.LGTP_S6(j)) + _
							(3.07579 * gs.ST_S6(j)) + (6.65858)

							odds_FA = 2.718281828 ^ logitFA_Use
							.pFA_Use(j) = odds_FA / (1.0 + odds_FA)
							.cf_FA_Use(j) = 0 
							 
							If .pFA_Use (j) > JBFa Then
								.cf_FA_Use(j) = 1 
							End If
							#endregion
								
							#region Code for Winter Use Model
							logitWI_Use = _
							(-4.8709 * gs.DTN_S6(j)) + _
							(-4.66829 * gs.CON_S6(j)) + _
							(-3.62043 * gs.LGW_S6(j)) + _
							(-0.851921 * gs.LGMD_S6(j)) + _
							(0.0889229 * gs.TDENLF_S6(j)) + _
							(0.559626 * gs.LGTP_S6(j)) + _
							(3.65392 * gs.LGOP_S6(j)) + _
							(11.3358 * gs.ST_S6(j)) + (-3.18369)
							 
							 
							odds_WI = 2.718281828 ^ logitWI_Use
							.pWI_Use(j) = odds_WI / (1.0 + odds_WI)
							.cf_WI_Use(j) = 0 
							 
							If .pWI_Use (j) > JBWi Then
								.cf_WI_Use(j) = 1 
							End If
							#endregion
						End If
					End If			
					#endregion"}

extractCoefTableLSLSeason <- function(text){
  data.frame(
    Season = text %>%
      str_extract("(?<=for )[[:upper:]].*(?= Use Model)"),
    Variable = c(text %>% 
                   str_extract_all("(?<=gs.)[[:upper:]]{2,6}(?=_S\\d\\(j\\))") %>%
                   .[which(!is.na(.))] %>% unlist(), "CONST"),
    Coefficient = c(text %>% 
                      str_extract_all("(?<=\\()(\\d|[[:punct:]])*(?= \\* gs)") %>%
                      .[which(!is.na(.))] %>% unlist(), text %>% 
                      str_extract("(?<=\\+ \\()(\\d|[[:punct:]])*(?=\\))")) %>% as.numeric(), 
    WinArea = text %>% 
      str_extract("(?<=gs.[[:upper:]]{2,6}_)S\\d(?=\\(j\\))"), 
    stringsAsFactors = FALSE
  )
}
extractCoefTableLSL <- function(text){
  rng <- text[1] %>%
    str_extract("(?<=code for )[[:upper:]].*(?= Caribou Range)")
  text <- text[2:5]
  map_dfr(text, extractCoefTableLSLSeason) %>%
    mutate(Range = rng, 
           WinArea = ifelse(WinArea == "S6", 10000, 5000))
}

coefTableHR <- inText %>% str_squish() %>%  str_split("#region Model") %>% 
  unlist() %>% 
  str_split("#region Code") %>% 
  .[c(2:8,11:16)] %>% #  discontinuous/Lake Superior not implemented in LSL
  map_dfr(extractCoefTableLSL)


usethis::use_data(coefTableHR, overwrite = TRUE)
