> fitaci
function (dat, varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", 
                               Ci = "Ci", PPFD = "PARi"), Tcorrect = TRUE, quiet = FALSE, 
          startValgrid = TRUE, algorithm = "default", debug = FALSE, 
          ...) 
{
  m <- as.list(match.call())
  a <- as.list(formals(fitaci))
  f <- names(formals(Photosyn))
  extrapars <- setdiff(names(m), c(names(a), ""))
  for (i in seq_along(extrapars)) {
    if (extrapars[i] %in% f) {
      val <- eval(m[[extrapars[i]]])
      formals(Photosyn)[extrapars[i]] <- val
    }
    else {
      warning("Parameter ", extrapars[i], " not recognized.")
    }
  }
  photpars <- formals(Photosyn)
  removevars <- c("whichA")
  photpars <- photpars[-which(names(photpars) %in% removevars)]
  if (!varnames$PPFD %in% names(dat)) {
    dat$PPFD <- 1800
    if (!quiet) 
      warning("PARi not in dataset; assumed PARi = 1800.")
  }
  else dat$PPFD <- dat[, varnames$PPFD]
  if (!varnames$Tleaf %in% names(dat)) {
    dat$Tleaf <- 25
    if (!quiet) 
      warning("Tleaf not in dataset; assumed Tleaf = 25.")
  }
  else dat$Tleaf <- dat[, varnames$Tleaf]
  dat$Ci <- dat[, varnames$Ci]
  dat$ALEAF <- dat[, varnames$ALEAF]
  TcorrectVJ <- Tcorrect
  acifun_wrap <- function(Ci, ...) {
    r <- Photosyn(Ci = Ci, Tcorrect = TcorrectVJ, ...)
    r$ALEAF
  }
  Rd_guess <- 1.5
  maxCi <- max(dat$Ci)
  mi <- which.max(dat$Ci)
  maxPhoto <- dat$ALEAF[mi]
  Tl <- dat$Tleaf[mi]
  gammastar <- TGammaStar(Tl)
  VJ <- (maxPhoto + Rd_guess)/((maxCi - gammastar)/(maxCi + 
                                                      2 * gammastar))
  Jmax_guess <- VJ * 4
  if (Tcorrect) {
    Teffect <- TJmax(Tl, EaJ = 39676.89, delsJ = 641.3615, 
                     EdVJ = 2e+05)
    Jmax_guess <- Jmax_guess/Teffect
  }
  dato <- dat[dat$Ci < 150 & dat$Ci > 60 & dat$ALEAF > 0, ]
  if (nrow(dato) > 0) {
    Km <- TKm(dato$Tleaf)
    gammastar <- TGammaStar(dato$Tleaf)
    vcmax <- with(dato, (ALEAF + Rd_guess)/((Ci - gammastar)/(Ci + 
                                                                Km)))
    Vcmax_guess <- median(vcmax)
  }
  else {
    Vcmax_guess <- Jmax_guess/1.8
  }
  if (Tcorrect) {
    Teffect <- TVcmax(Tl, EaV = 82620.87, delsC = 645.1013, 
                      EdVC = 0)
    Vcmax_guess <- Vcmax_guess/Teffect
  }
  if (startValgrid) {
    aciSS <- function(Vcmax, Jmax, Rd) {
      Photo_mod <- acifun_wrap(dat$Ci, PPFD = dat$PPFD, 
                               Vcmax = Vcmax, Jmax = Jmax, Rd = Rd, Tleaf = dat$Tleaf)
      SS <- sum((dat$ALEAF - Photo_mod)^2)
      return(SS)
    }
    d <- 0.3
    n <- 20
    gg <- expand.grid(Vcmax = seq(Vcmax_guess * (1 - d), 
                                  Vcmax_guess * (1 + d), length = n), Rd = seq(Rd_guess * 
                                                                                 (1 - d), Rd_guess * (1 + d), length = n))
    m <- with(gg, mapply(aciSS, Vcmax = Vcmax, Jmax = Jmax_guess, 
                         Rd = Rd))
    ii <- which.min(m)
    Vcmax_guess <- gg$Vcmax[ii]
    Rd_guess <- gg$Rd[ii]
  }
  if (debug) 
    browser()
  nlsfit <- nls(ALEAF ~ acifun_wrap(Ci, PPFD = PPFD, Vcmax = Vcmax, 
                                    Jmax = Jmax, Rd = Rd, Tleaf = Tleaf), algorithm = algorithm, 
                data = dat, control = nls.control(maxiter = 500, minFactor = 1/10000), 
                start = list(Vcmax = Vcmax_guess, Jmax = Jmax_guess, 
                             Rd = Rd_guess))
  p <- coef(nlsfit)
  acirun <- Photosyn(Ci = dat$Ci, Vcmax = p[[1]], Jmax = p[[2]], 
                     Rd = p[[3]], PPFD = dat$PPFD, Tleaf = dat$Tleaf, Tcorrect = Tcorrect)
  acirun$Ameas <- dat$ALEAF
  acirun$ELEAF <- NULL
  acirun$GS <- NULL
  acirun$Ca <- NULL
  names(acirun)[names(acirun) == "ALEAF"] <- "Amodel"
  avars <- match(c("Ci", "Ameas", "Amodel"), names(acirun))
  acirun <- acirun[, c(avars, setdiff(1:ncol(acirun), avars))]
  l <- list()
  l$df <- acirun[order(acirun$Ci), ]
  l$pars <- summary(nlsfit)$coefficients[, 1:2]
  l$nlsfit <- nlsfit
  l$Tcorrect <- Tcorrect
  formals(Photosyn)$Tleaf <- mean(dat$Tleaf)
  formals(Photosyn)$PPFD <- mean(dat$PPFD)
  formals(Photosyn)$Vcmax <- l$pars[1]
  formals(Photosyn)$Jmax <- l$pars[2]
  formals(Photosyn)$Rd <- l$pars[3]
  formals(Photosyn)$Tcorrect <- Tcorrect
  l$Photosyn <- Photosyn
  l$Ci_transition <- findCiTransition(l$Photosyn)
  l$Vcmax_guess <- Vcmax_guess
  l$Jmax_guess <- Jmax_guess
  l$Rd_guess <- Rd_guess
  class(l) <- "acifit"
  return(l)
}