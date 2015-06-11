

# source("functions/functions.R")
### ACi model for Vpmax and Vc max only
c4_aci<- function (Ci, PPFD = 1500, Tleaf = 25, VPMAX25 = 120, Vcmax = 60, 
                   Vpr = 80, alpha = 0, gbs = 0.003, O2 = 210, JMAX25 = 400, 
                   x = 0.4, THETA = 0.7, RD0 = 1, RTEMP = 25, TBELOW = 0,
                   low_gammastar = 0.00019305, Kc = 650,Kp = 80, Ko= 450,
                   DAYRESP = 1, FRM = 0.5, Rd= 1, ...) 
{
  Vcmax <- Vcmax
  Vpmax <- VPMAX25
  Jmax <- JMAX25
  TK <- Tleaf + 273.15
  K <- Kc * (1 + O2/Ko)
  Rd<- Rd
  Rm <- FRM * Rd
  Vp <- pmin(Ci * Vpmax/(Ci + Kp), Vpr)
  a.c <- 1 - (alpha * Kc)/(0.047 * Ko)
  b.c <- -((Vp - Rm + gbs * Ci) + (Vcmax - Rd) + gbs * K + 
             alpha * low_gammastar/0.047 * (low_gammastar * Vcmax + 
                                              Rd * Kc/Ko))
  c.c <- (Vcmax - Rd) * (Vp - Rm + gbs * Ci) - (Vcmax * gbs * 
                                                  low_gammastar * O2 + Rd * gbs * K)
  
  A.enzyme <- (-b.c - sqrt(b.c^2 - 4 * a.c * c.c))/(2 * a.c)
  
  Qp2 <- PPFD * (1 - 0.15)/2
  
  J <- (1/(2 * THETA)) * (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 
                                              4 * THETA * Qp2 * Jmax))
  
  a.j <- 1 - 7 * low_gammastar * alpha/(3 * 0.047)
  
  b.j <- -((x * J/2 - Rm + gbs * Ci) + ((1 - x) * J/3 - Rd) + 
             gbs * (7 * low_gammastar * O2/3) + alpha * low_gammastar/0.047 * 
             ((1 - x) * J/3 + Rd))
  c.j <- ((x * J/2 - Rm + gbs * Ci) * ((1 - x) * J/3 - Rd) - 
            gbs * low_gammastar * O2 * ((1 - x) * J/3 - 7 * Rd/3))
  
  
  A.light <- (-b.j - sqrt(b.j^2 - 4 * a.j * c.j))/(2 * a.j)
  
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

fitc4 <- function(dat, varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", 
        Ci = "Ci", PPFD = "PPFD", Rd= "Rd"),algorithm="default", ...)
{  
  m <- as.list(match.call())
   a <- as.list(formals(fitc4))
   f <- names(formals(c4_aci))
   extrapars <- setdiff(names(m), c(names(a), ""))
   for (i in seq_along(extrapars)) {
     if (extrapars[i] %in% f) {
       val <- eval(m[[extrapars[i]]])
       formals(c4_aci)[extrapars[i]] <- val
     }
     else {
       warning("Parameter ", extrapars[i], " not recognized.")
     }
   }
   
  acifun_wrap <- function(Ci, ...) {
    r <- c4_aci(Ci = Ci,...)
    r$ALEAF
  }
  
  aciSS <- function(VPMAX25, Vcmax) {
    Photo_mod <- acifun_wrap(dat$Ci, PPFD = dat$PPFD,VPMAX25=VPMAX25,Vcmax = Vcmax, JMAX25  = 400, Rd = dat$Rd, Tleaf = dat$Tleaf)
    SS <- sum((dat$ALEAF - Photo_mod)^2)
    return(SS)
  }
  
  aciSS(vpmax_guess, vcmax_guess)
  
  d <- 0.3
  n <- 20
  
  gg <- expand.grid(VPMAX25 = seq(vpmax_guess * (1 - d),vpmax_guess* (1 + d), length = n),Vcmax = seq(vcmax_guess * (1 - d),vcmax_guess* (1 + d), length = n))
  ## 
  m <- with(gg, mapply(aciSS,VPMAX25=VPMAX25, Vcmax = Vcmax))
  ii <- which.min(m)
  
  vpmax_guess <- gg$VPMAX25[ii]
  vcmax_guess<- gg$Vcmax[ii]
  
  
  nlsfit <- nls(ALEAF ~ acifun_wrap(Ci, PPFD = PPFD, VPMAX25= VPMAX25,Vcmax = Vcmax), 
                data = dat, control = nls.control(maxiter = 500, minFactor = 1/10000), 
                start = list(VPMAX25= vpmax_guess,Vcmax = vcmax_guess))
  
  p <- coef(nlsfit)
  acirun <- c4_aci(Ci = dat$Ci,  VPMAX25 = p[[1]],  Vcmax  = p[[2]], PPFD = dat$PPFD, Tleaf = dat$Tleaf, Rd= dat$Rd)
  
  acirun$Ameas <- dat$ALEAF
  acirun$Ci <- dat$Ci
  names(acirun)[names(acirun) == "ALEAF"] <- "Amodel"
  avars <- match(c("Ci", "Ameas", "Amodel"), names(acirun))
  acirun <- acirun[, c(avars, setdiff(1:ncol(acirun), avars))]
  
  l <- list()
  l$df <- acirun[order(acirun$Ci), ]
  l$pars <- summary(nlsfit)$coefficients[, 1:2]
  l$nlsfit <- nlsfit
  formals(c4_aci)$Tleaf<- mean(dat$Tleaf)
  formals(c4_aci)$PPFD <- mean(dat$PPFD)
  formals(c4_aci)$VPMAX25 <- l$pars[1]
  formals(c4_aci)$Vcmax <- l$pars[2]
  l$c4_aci <- c4_aci
  l$Vpmax_guess <- vpmax_guess
  l$Vcmax_guess <- vcmax_guess
  class(l) <- "c4acifit"
  return(l)
}



print.c4acifit <- function (x, ...){
  cat("Result of C4 fitaci.\n\n")
  cat("Data and predictions:\n")
  print(x$df)
  cat("\nEstimated parameters:\n")
  print(x$pars)

    cat("\n\n")
  cat("Parameter settings:\n")
  fm <- formals(x$c4_aci)
  pars <- c("Vpr","THETA","Kc")
  fm <- fm[pars]
  cat(paste0(names(fm), " = ", unlist(fm), "\n"))
}




plot.c4acifit  <- function (x, what = c("data", "model"), xlim = NULL, ylim = NULL, 
          whichA = c("Ac", "Vp", "Amin"), add = FALSE, pch = 19, addzeroline = TRUE, 
          addlegend = !add, lwd = c(1, 2), 
          ...) 
{
  if (is.null(ylim)) 
    ylim <- with(x$df, c(min(Ameas), 1.1 * max(Ameas)))
  if (is.null(xlim)) 
    xlim <- with(x$df, c(0, max(Ci)))
  if (length(lwd) == 1) 
    lwd <- c(lwd, lwd)
  Ci <- with(x$df, seq(min(Ci), max(Ci), length = 101))
  pred <- x$c4_aci(Ci = Ci)
  if (!add) {
    with(x$df, plot(Ci, Ameas, type = "n", ylim = ylim, xlim = xlim, 
                    xlab = expression(italic(C)[i] ~ ~(ppm)), ylab = expression(italic(A)[net] ~ 
                                                                                  ~(mu * mol ~ m^-2 ~ s^-1)), ...))
  }
  if ("data" %in% what) 
    with(x$df, points(Ci, Ameas, pch = pch, ...))
  if ("model" %in% what) {
    if ("Vp" %in% whichA) 
      with(pred, points(Ci, Vp - Rd, type = "l", col = "blue", 
                        lwd = lwd[1]))
    if ("Ac" %in% whichA) 
      with(pred, points(Ci, Ac - Rd, type = "l", col = "red", 
                        lwd = lwd[1]))
    if ("Amin" %in% whichA) 
      with(pred, points(Ci, ALEAF, type = "l", col = "black", 
                        lwd = lwd[2]))
  }

  if (addzeroline) 
    abline(h = 0, lty = 3)
  if (addlegend) {
    legend("bottomright", c(expression(italic(A)[c]), expression(italic(V)[p]), 
                            "Limiting rate"), lty = 1, lwd = c(1, 1, 2), col = c("red", 
                                                                                 "blue", "black"))
  }
}




