
rm(list=ls())
require(plyr)
require(dplyr)
require(plantecophys)
source("functions/functions.R")

##read aci function 
source("functions/bala_c4test_remkoedits_v1.R")

## read the csv file
df<- read.csv("rawdata/maize25_aci.csv")


##rename the variable as defined in functions
df$PPFD<- df$PARi
df$ALEAF<- df$Photo


##starting guess value for function
vpmax_guess<- 70
vcmax_guess<-30
gbs_guess<- 0.003



f1<- fitc4(df)
f1
plot(f1)


