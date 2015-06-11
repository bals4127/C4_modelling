#Get Remko's R package
library(plantecophys)

####################### FUNCTIONS
#functions I copied from FARAO because not sure how to reference them as part of package
.Rgas <- function()8.314
Tk <- function(x)x+273.15

# Arrhenius
arrh <- function(Tleaf, Ea){
  exp((Ea * (Tk(Tleaf) - 298.15)) / (298.15 * .Rgas() * Tk(Tleaf))) 
}
#Gamma Star
TGammaStar <- function(Tleaf, 
                       Egamma=37830.0, 
                       value25=42.75){  
  
  value25*arrh(Tleaf,Egamma)
}
#Quadratic equation
QUADM <- function(A,B,C){
  
  (- B - sqrt(B^2 - 4*A*C)) / (2*A)
  
}

#Generic rectangular hyperbola - for fitting to C4 A-Ci curve
recthyp <- function(Amax=50, Km=100, Ci=NONE) {
  Amax*Ci/(Ci+Km)
}

# The quadratic solution to the optimal gs model for a generic function 
# of form Assim = Vmax*(Ci-GStar)/(Ci+Km)
# lambda is defined as: maximise A - lambda E
# A is is umol CO2 m-2 s-1, gs and E are in mol H2O m-2 s-1 
# lambda is passed in mmol CO2 mol-1 H2O (it should be around about 1 in these units)
fquad <- function(lambda=2, Assim=30, VPD=1, Ca=400, Gstar=42, Km=100) {
  
  Press = 100 # KPa
  L = 1.6 * VPD / Press * lambda*1E3  #units umol CO2 mol-1 H2O
  a = Km + Gstar - L
  b = 2 * L * Gstar - 2 * Ca * (Gstar + Km)
  c = Ca * Ca * (Gstar+Km) - L*Ca*(Gstar+Km)+L*Gstar*Km
  Ci = QUADM(a,b,c)
  gs = 1.6 * Assim / (Ca - Ci)  # DO NOT FORGET THE 1.6!!!
  return(gs)
  
}

# The BBOpti model
fBBOpti <- function(g1=4, Assim=NONE, VPD = NONE, Ca = NONE) {
  
  1.6 * (1.0 + g1/sqrt(VPD)) * Assim/Ca
       
}

# The BB Leuning model
fBBL <- function(a1=1, D0 = 1.5, Assim=NONE, VPD = NONE, Ca = NONE) {
  
  a1 * Assim/Ca / (1+VPD/D0)
  
}

# wrappers for FARAO C4 functions, so they only return A and Gs, to which we are fitting
# This allows us to fit numerical solution
fAC4 <- function(Vcmax, Vpmax25, VPD, PPFD, Ci, Tair)
                      AciC4(Vcmax=Vcmax, VPMAX25=Vpmax25,
                      VPD=VPD, PPFD=PPFD, Ci=Ci, Tleaf=Tair, Q10F=0.0693)$ALEAF
fC4num <- function(lambda, Vcmax, Vpmax, VPD, PPFD, Ca, Tair)
                      FARAO(lambda=lambda/1000, #is in mol CO2 mol-1 H2O
                      Vcmax=Vcmax, VPMAX25=Vpmax, Q10F=0.0693,
                      VPD=VPD, PPFD=PPFD, Ca=Ca, Tair=Tair, photo="BOTH", C4=TRUE)$GS
################ END FUNCTIONS

##### READ IN DATA
# Reading in gas exchange data - then do one species at a time
alldata <- read.csv("averaged_data.csv")
#speciesdata <- subset(alldata, Species=="PC")
#speciesdata <- subset(alldata, Species=="FB")
#speciesdata <- subset(alldata, Species=="PM")
#speciesdata <- subset(alldata, Species=="PD")
#speciesdata <- subset(alldata, Species=="FP")
#speciesdata <- subset(alldata, Species=="PB")
speciesdata <- subset(alldata, Species=="SL")

speciesdata <- subset(speciesdata, VpdL < 2.3) # needed for Flaverias because the A-Ci at high VPD is very low
attach(speciesdata)

##### FIT A-CI CURVES
#fit curves to A-Ci curve - firstly rectangular hyperbola, then full model
fitAciC4rect <- nls(Photo ~ recthyp(Amax, Km, Ci), 
                start=list(Amax=40, Km=100),trace=TRUE)
fitAciC4num <- nls(Photo ~ fAC4(Vcmax, Vpmax25, VpdL, PARi, Ci, Tair),
                   start=list(Vcmax=40, Vpmax25 = 60),trace=TRUE,
                   control=(nls.control(warnOnly=TRUE)))
coef(summary(fitAciC4rect))
coef(summary(fitAciC4num))

#Visualise A-Ci curves
fittedArect <-recthyp(coef(fitAciC4rect)[1],coef(fitAciC4rect)[2],Ci)
Vcmax25=coef(fitAciC4num)[1]
Vpmax25=coef(fitAciC4num)[2]
fittedAnum <- fAC4(Vcmax25,Vpmax25,VpdL,PARi,Ci,Tair)
plot(Ci, Photo, col=VpdL*2)
points(Ci, fittedArect,type="l")
points(Ci, fittedAnum, col="orange")
legend("bottomright", c("VPD=1","VPD=1.5","VPD=2","VPD=2.5","hyperbola","C4 model"),
       lwd=2, col=c(2,3,4,5,"black","orange"))

###### OPTIMAL STOMATAL FITS
# fit gs to COND using different versions of optimal calculation

#for C4, first use quadratic solution with Km fitted directly to A-Ci curve
Amrect <- coef(fitAciC4rect)[1]
Kmrect <- coef(fitAciC4rect)[2]
fitC4 <- nls(Cond ~ fquad(lambda,Assim=Photo,VPD=VpdL,Ca=CO2S,Km=Kmrect,Gstar=0),
             start=list(lambda=0.3),trace=TRUE)
f <- function(x){sum((fquad(x, Assim=Photo,VPD=VpdL,Ca=CO2S,Km=Kmrect,Gstar=0)-Cond)^2)}
minlam <- optimize(f,c(0.05,2),tol=0.0001)
minlam

#next try quadratic solution for C3 Jmax limited
Gstar=mean(TGammaStar(Tleaf))
fitC3 <- nls(Cond ~ fquad(lambda, Assim=Photo, VPD=VpdL, Ca=CO2S, Km=2*Gstar,Gstar=Gstar),
             start=list(lambda=0.5), trace=TRUE)
f <- function(x){sum((fquad(x, Assim=Photo,VPD=VpdL,Ca=CO2S,Km=2*Gstar,Gstar=Gstar)-Cond)^2)}
minlam <- optimize(f,c(0.05,2),tol=0.0001)
minlam
g1lam <- sqrt(3*Gstar*101/1.6/1000/coef(fitC3)[1])
g1lam

#Next use BBOpti and calculate lambda from g1 - omitting as this seems to give wrong values? 
fitBBOpti <- nls(Cond ~ fBBOpti(g1, Photo, VpdL, CO2S),
                 start=list(g1=4), trace = TRUE)
fitBBL <- nls(Cond ~ fBBL(a1, D0, Photo, VpdL, CO2S),
                 start=list(a1=4, D0=1.5), trace = TRUE)
g1 <-coef(fitBBOpti)[1]
lambdag1 <- 3*Gstar*101/1.6/(g1*g1)/1000
lambdag1

#fit to numerical function - gives warning
fitC4num <- nls(Cond ~ fC4num(lambda, Vcmax=Vcmax25, Vpmax=Vpmax25, VPD=VpdL, PPFD=PARi, Ca=CO2S, Tair=Tair), 
                start=list(lambda=0.447),trace=TRUE,
                control=(nls.control(warnOnly=TRUE,minFactor=1/8192)))
f <- function(x){sum((fC4num(x, Vcmax=Vcmax25, Vpmax=Vpmax25, VPD=VpdL, PPFD=PARi, Ca=CO2S, Tair=Tair)-Cond)^2)}
minlam <- optimize(f,c(0.05,2),tol=0.0001)
minlam

#Calculate predicted values
BBM <- Photo / CO2S / sqrt(VpdL)
lambdaC3 <- coef(fitC3)[1]
predC3gs <- fquad(lambda=lambdaC3,Assim=Photo,VPD=VpdL,Ca=CO2S,Gstar,2*Gstar)
lambdaC4 <- coef(fitC4)[1]
predC4gs <- fquad(lambda=lambdaC4,Assim=Photo,VPD=VpdL,Ca=CO2S,0.0,Kmrect)
lambdaC4num <- coef(fitC4num)[1]
predC4gsnum <- fC4num(lambdaC4num, Vcmax25, Vpmax25, VpdL, PARi, CO2S, 25)
predC4l3gs <- fquad(lambda=lambdaC4,Assim=Photo,VPD=VpdL,Ca=CO2S,0.0,Kmrect)
g1 <- coef(fitBBOpti)[1]
predBBOpti <- fBBOpti(g1,Photo,VpdL,CO2S)

####### VISUALISE fits to data
#vs measured gs
plot(Cond,predC4gsnum,col=VpdL*2)
points(Cond,predC4gs,col="purple")
points(Cond,predC3gs,col="black")
points(Cond,predBBOpti,col="pink")
#vs VPD
plot(CO2S,Cond,col=VpdL*2)
points(CO2S,predC4gsnum,col="black")
points(CO2S,predC4gs,col="purple")
#vs BBM
plot(BBM,Cond,col=VpdL*2)
points(BBM,predC3gs,col="black")
points(BBM,predC4gs,col="purple")
points(BBM,predBBOpti,col="orange")
points(BBM,predC4gsnum,col="pink")
#vs each other
plot(predC4gs,predC4gsnum,col=VpdL*2)
points(predC3gs,predC4gsnum,col="orange")
points(predC3gs,predC4l3gs,col="red")

##### STATS 
#Compare models with data and with other models - gives R2 values
summary(lm(predC4gsnum~predC4gs-1))
summary(lm(predC4gs~predC3gs))
summary(lm(predC4gsnum~Cond))
summary(lm(predC4gs~Cond))
summary(lm(predC3gs~Cond))
summary(lm(predBBOpti~Cond))


#STUFF - attempting to visualise A - lambda E
cis <- seq (10,250,length=20)
Amod <- AciC4(Vcmax=Vcmax25,VPMAX25=Vpmax25,Ci=cis,Tleaf=26.4,Q10F=0.0693)
Aquad <- recthyp(Amax=Amrect,Km=Kmrect,Ci=cis)
plot(cis,Amod$ALEAF)
points(cis,Aquad,col="red")
Emod <- 1.6*Amod$ALEAF/(550.0-cis)*2.5/101  # in mol m-2 s-1
Equad <- 1.6*Aquad/(550.0-cis)*2.5/101
points(cis,Amod$ALEAF-3*Emod*1E3,col="blue")
points(cis,Aquad-lambdaC4num*Equad*1E3,col="orange")
points(cis,Aquad-lambdaC4*Equad*1E3,col="green")

lamval <- seq(0.1,1,length=30)
ssq <- lapply(lamval, function(x){sum((fC4num(x, Vcmax=Vcmax25, Vpmax=Vpmax25, VPD=VpdL, PPFD=PARi, Ca=CO2S, Tair=Tair)-Cond)^2)})
plot(lamval,ssq)
f <- function(x){sum((fC4num(x, Vcmax=Vcmax25, Vpmax=Vpmax25, VPD=VpdL, PPFD=PARi, Ca=CO2S, Tair=Tair)-Cond)^2)}
minlam <- optimize(f,c(0.05,2),tol=0.0001)