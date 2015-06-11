rm(list=ls())
require(dplyr)
source("functions/bala_c4test.R") ## this file contain all the function created for C4 fitting
source("functions/load packages.R")
## read dataframe
df<- read.csv("rawdata/aci_c4.csv")

df<- df[, c("code","growth","Photo","Cond", "Ci","Tleaf", "PARi","CO2R", "RH_S")]

df$spp<- lab_spp(df)

df$treat<- df$growth

df$PPFD<- df$PARi
df$ALEAF<- df$Photo

df$spp<- as.factor(df$spp)

c4<- droplevels(df%>%
                  filter(!spp== "P. bisulcatum ")%>%
                  filter(!spp== "P. milliodes "))

c3<- droplevels(df%>% filter(spp== "P. bisulcatum "  | spp== "P. milliodes "))

c425<- c4%>% filter(Tleaf<30)
c435<-c4%>% filter(Tleaf>30)

c325 <- c3%>% filter(Tleaf<30)
c335<-c3%>% filter(Tleaf>30)

vpmax_guess<- 70
vcmax_guess<-30



## need to subset dat to remove cilliaris cool code e1.2. work out on it latter


mz<-  filter(c425, code== "m1.1")
fitmz<- fitc4(mz)
summary(fitmz)

plot(fitmz)



fit<- dlply(c425, .(code),function(x) fitc4(x))

const <- ldply(fit, function(x) coef(x))



## In 35C aci, no prblem in the data collected from warm grown plants
## need to subset dat to remove cilliaris  code e1.2, and coloratum f1.1




n<- filter(c435, !spp=="P. coloratum ")
n1<- filter(n, !spp== "C. ciliaris ")
unique(n1$spp)

fit_w<- dlply(n1, .(code),function(x) fitc4(x))
const <- ldply(fit_w, function(x) coef(x))
const


## fit aci for C3 grasses
fit_c3<- dlply(c3, .(code),function(x) fitaci(x))
const <- ldply(fit_w, function(x) coef(x))





plotBy

require(plotBy)


