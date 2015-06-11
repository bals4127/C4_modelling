rm(list=ls())
require(plyr)
source("functions/functions.R")
require(plantecophys)
### ACi model for Vpmax and Vc max only
c4_acij<- function (Ci, PPFD = 1500, Tleaf = 25, VPMAX25 = 120, Vcmax = 60, JMAX25 = 400,
                    Vpr = 80, alpha = 0, gbs = 0.003, O2 = 210,x = 0.4, THETA = 0.7, 
                    RD0 = 1, RTEMP = 25, TBELOW = 0,DAYRESP = 1, FRM = 0.5, ...) 
{
  Vcmax <- Vcmax
  Vpmax <- VPMAX25
  Jmax <- JMAX25
  TK <- Tleaf + 273.15
  low_gammastar <- 0.00019305
  Kc <- 650 
  Kp <- 80 
  Ko <- 450 
  K <- Kc * (1 + O2/Ko)
  Rd<-1
  Rm <- FRM * Rd
  Vp <- pmin(Ci * Vpmax/(Ci + Kp), Vpr)
  a.c <- 1 - (alpha * Kc)/(0.047 * Ko)
  b.c <- -((Vp - Rm + gbs * Ci) + (Vcmax - Rd) + gbs * K +alpha /0.047 * (low_gammastar * Vcmax +Rd * Kc/Ko))
  c.c <- (Vcmax - Rd) * (Vp - Rm + gbs * Ci) - (Vcmax * gbs *low_gammastar * O2 + Rd * gbs * K)
  
  A.enzyme <- (-b.c - sqrt(b.c^2 - 4 * a.c * c.c))/(2 * a.c)
  
  Qp2 <- PPFD *0.85* (1 - 0.15)/2
  
  J <- (1/(2 * THETA)) * (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 -4 * THETA * Qp2 * Jmax))
  
  a.j <- 1 - 7 * low_gammastar * alpha/(3 * 0.047)
  
  b.j <- -((x * J/2 - Rm + gbs * Ci) + ((1 - x) * J/3 - Rd) + gbs * (7 * low_gammastar * O2/3) + alpha * low_gammastar/0.047 * 
             ((1 - x) * J/3 +7* Rd/3))
  c.j <- ((x * J/2 - Rm + gbs * Ci) * ((1 - x) * J/3 - Rd) -gbs * low_gammastar * O2 * ((1 - x) * J/3 - 7 * Rd/3))
  
  lc<-0.9999
  
  A.light <- (-b.j - sqrt(b.j^2 - 4 * lc*a.j * c.j))/(2 * lc*a.j)
  
  An <- pmin(A.enzyme, A.light)
  
  Ac <- A.enzyme
  
  Aj <- A.light
  
  
  shape2 <- 0.999
  Ad <- (Ac + Aj - sqrt((Ac + Aj)^2 - 4 * shape2 * Ac * Aj))/(2 * shape2) - Rd
  Ac <- Ac - Rd
  Aj <- Aj - Rd
  
  df <- data.frame (ALEAF = Ad, Ci = Ci,An = An, Ac = Ac, Aj = Aj, Vp = Vp, 
                    Rd = Rd, Tleaf = Tleaf, PPFD = PPFD)
  return(df)
  
}


## fiting the model

fitc4t<- function(dat, varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "PPFD"))
{  
  m <- as.list(match.call())
  a <- as.list(formals(fitc4t))
  f <- names(formals(c4_acij))
  
  extrapars <- setdiff(names(m), c(names(a), ""))
  
  for (i in seq_along(extrapars)) {
    if (extrapars[i] %in% f) {
      val <- eval(m[[extrapars[i]]])
      formals(c4_acij)[extrapars[i]] <- val
    }
    else {
      warning("Parameter ", extrapars[i], " not recognized.")
    }
  }
  
  
  acifun_wrap <- function(Ci, ...) {
    r <- c4_acij(Ci = Ci, ...)
    r$ALEAF
  }
  
  aciSSj <- function(VPMAX25, Vcmax, JMAX25) {
    Photo_mod <- acifun_wrap(dat$Ci, PPFD = dat$PPFD,VPMAX25=VPMAX25,Vcmax = Vcmax, JMAX25  = JMAX25, Rd = Rd, Tleaf = dat$Tleaf)
    SS <- sum((dat$ALEAF - Photo_mod)^2)
    return(SS)
  }
  d <- 0.3
  n <- 20
  gg <- expand.grid(VPMAX25 = seq(vpmax_guess * (1 - d),vpmax_guess* (1 + d), length = n),
                    Vcmax = seq(vcmax_guess * (1 - d),vcmax_guess* (1 + d), length = n),
                    JMAX25 = seq(jmax_guess * (1 - d),jmax_guess* (1 + d), length = n))
  
  m<- with(gg, mapply(aciSSj,VPMAX25=VPMAX25, Vcmax = Vcmax,JMAX25=JMAX25))
  
  ii <- which.min(m)
  vpmax_guess <- gg$VPMAX25[ii]
  vcmax_guess<- gg$Vcmax[ii]
  jmax_guess<- gg$JMAX25[ii]
  
  nlsfit <- nlsLM(ALEAF ~ acifun_wrap(Ci, PPFD = PPFD, VPMAX25= VPMAX25,Vcmax = Vcmax, JMAX25=JMAX25), 
                  data = dat, control = nls.control(maxiter = 500, minFactor = 1/10000), 
                  start = list(VPMAX25= vpmax_guess,Vcmax = vcmax_guess, JMAX25= jmax_guess))
  
  p <- coef(nlsfit)
  acirun <- c4_acij(Ci = dat$Ci,  VPMAX25 = p[[1]],  Vcmax  = p[[2]], JMAX25 = p[[3]],PPFD = dat$PPFD, Tleaf = dat$Tleaf)
  acirun$Ameas <- dat$ALEAF
  acirun$Ci <- dat$Ci
  names(acirun)[names(acirun) == "ALEAF"] <- "Amodel"
  avars <- match(c("Ci", "Ameas", "Amodel"), names(acirun))
  acirun <- acirun[, c(avars, setdiff(1:ncol(acirun), avars))]
  
  l <- list()
  l$df <- acirun[order(acirun$Ci), ]
  l$pars <- summary(nlsfit)$coefficients[, 1:2]
  l$nlsfit <- nlsfit
  formals(c4_acij)$Tleaf<- mean(dat$Tleaf)
  formals(c4_acij)$PPFD <- mean(dat$PPFD)
  formals(c4_acij)$Vpmax <- l$pars[1]
  formals(c4_acij)$Vcmax <- l$pars[2]
  formals(c4_acij)$JMAX25 <- l$pars[3]
  l$Vpmax_guess <- vpmax_guess
  l$Vcmax_guess <- vcmax_guess
  class(l) <- "acifit"
  return(l)
}


