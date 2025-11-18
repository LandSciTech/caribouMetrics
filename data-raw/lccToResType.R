## code to prepare `lccToResType` dataset goes here
# not currently in use but keep this in case we want to use FRI data again
CON_ls <- expand.grid(PLCCode=c(1,6,7),ResourceType="CON")
DEC_ls <- expand.grid(PLCCode=c(2,11,12),ResourceType="DEC")
DTN_ls <- expand.grid(PLCCode=c(34,35),ResourceType="DTN")
LGOP_ls <- expand.grid(PLCCode=c(19,31,32),ResourceType="LGOP")
LGTP_ls <- expand.grid(PLCCode=c(8,10),ResourceType="LGTP")
LGW_ls <- expand.grid(PLCCode=c(37,38),ResourceType="LGW")
MIX_ls <- expand.grid(PLCCode=c(3,4,5,13,14,15),ResourceType="MIX")
ST_ls <- expand.grid(PLCCode=c(9,20),ResourceType="ST")
OTHER_ls <- expand.grid(PLCCode=c(16,17,18,21,22,23,24,25,26,27,28,29,30,33,36,39),ResourceType="other")

lccToResType <- rbind(CON_ls,DEC_ls,DTN_ls,LGOP_ls,LGTP_ls,LGW_ls,MIX_ls,ST_ls,OTHER_ls)
lccToResType <- lccToResType[order(lccToResType$PLCCode),]

usethis::use_data(lccToResType, overwrite = TRUE)
