rm(list=ls())
require(dplyr)
#source("functions/bala_c4test_20150501_edits.R") ## this file contain all the function created for C4 fitting
source("functions/load packages.R")
## read dataframe
df<- read.csv("rawdata/c4test.csv")




df$ALEAF<- df$Photo

vpmax_guess<- 70
vcmax_guess<-30
jmax_guess<- 400

fit<- fitc4(df)
coef(fit)


fitj<- fitc4j(df)
coef(fitj)
