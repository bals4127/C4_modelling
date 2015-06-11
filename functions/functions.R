lab_mkcat<- function(dfr)
{
 dfr$mkat[dfr$species == "P. bisulcatum "] <- 2.6
 dfr$mkat[dfr$species == "P. milliodes "] <- 2.2
 dfr$mkat[dfr$species == "P. monticola "] <- 5.3
 dfr$mkat[dfr$species == "C. ciliaris "] <- 6
 dfr$mkat[dfr$species == "P. coloratum "] <- 3.4
 dfr$mkat[dfr$species == "E. curvula "] <- 5.3
 dfr$mkat[dfr$species == "P. maximum "] <- 5.3
 dfr$mkat[dfr$species == "C. gayana "] <- 5.7
 dfr$mkat[dfr$species == "S. bicolor "] <- 5.8
 dfr$mkat[dfr$species == "Z. mays "] <- 5.5
 return(dfr$mkat)
}

## function to add treatment column is dataset
lab_treat<- function (dfr){
  dfr$treat<- gsub("[:a-z:]", " ",dfr$code) # remove all alphabets from code
  dfr$treat<- as.numeric(dfr$treat)         # conver in numeric vector
  dfr$treat<- round(dfr$treat, 0)           # round up into single number
  dfr$treat[dfr$treat > 2.9999]<- "warm"    # name treatment equal or bigger than three as a warm
  dfr$treat[dfr$treat < 3]<- "cool"         # name as cool
  return(dfr$treat)                        # return vector of label
  }


### function to get species label
lab_spp<- function (dfr){
  dfr$test4<- gsub("[:1-9:]", " ",dfr$code)           # remove numbers from code
  dfr$sppa<- gsub("^a .","P. bisulcatum",dfr$test4)   #give lable based on the alphabet coding 
  dfr$sppc<- gsub("^c .","P. milliodes",dfr$sppa)
  dfr$sppd<- gsub("^d .","P. monticola",dfr$sppc)
  dfr$sppe<- gsub("^e .","C. ciliaris",dfr$sppd)
  dfr$sppf<- gsub("^f .","P. coloratum",dfr$sppe)
  dfr$sppg<- gsub("^g .","E. curvula",dfr$sppf)
  dfr$spph<- gsub("^h .","P. maximum",dfr$sppg)
  dfr$sppi<- gsub("^i .","C. gayana",dfr$spph)
  dfr$sppj<- gsub("^j .","S. bicolor",dfr$sppi)
  dfr$sppm<- gsub("^m .","Z. mays",dfr$sppj)
  return(dfr$sppm)                             # return the vector of species label
}




se<- function(x, ...) {
  se <- sd(x, ...)/sqrt(sum(!is.na(x)))
  return(se)
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
    # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
    # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
    # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(sum(!is.na(datac$N)))
  
   # Calculate standard error of the mean
    # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
    return(datac)
}


bar <- function(dv, factors, dataframe, percentage=FALSE, errbar=!percentage, half.errbar=TRUE, conf.level=.95, 
                xlab=NULL, ylab=NULL, main=NULL, names.arg=NULL, bar.col="black", whisker=.015,args.errbar=NULL,
                legend=TRUE, legend.text=NULL, args.legend=NULL,legend.border=FALSE, box=TRUE, args.yaxis=NULL, 
                mar=c(5,4,3,2),...){
  axes=!percentage
  dv.name<-substitute(dv)
  if(length(dv.name)>1) stop("'dv' only takes one variable")
  dv.name<-as.character(dv.name)
  dv<-dataframe[[dv.name]]
  fnames<-substitute(factors)
  if(length(fnames)==1){
    factors<-as.character(fnames)
    nf<-1
  }else{
    factors<-as.character(fnames[-1L])
    nf<-length(factors)
  }
  if(nf>2) stop("This function accepts no more than 2 factors \n",
                "\t-i.e., it only plots one-way or two-way designs.")
  if(percentage & errbar){
    warning("percentage=TRUE; error bars were not plotted")
    errbar<-FALSE
  }
  if(!percentage) xbars<-tapply(dv, dataframe[,factors], mean, na.rm=TRUE)
  else {
    xbars<-tapply(dv, list(interaction(dataframe[,factors], lex.order=TRUE)), mean, na.rm=TRUE)
    if(sum(na.omit(dv)!=0&na.omit(dv)!=1)>0) 
      stop("Data points in 'dv' need to be 0 or 1 in order to set 'percentage' to TRUE")
    xbars<-rbind(xbars, 1-xbars)*100
  }
  if(errbar){
    se<-tapply(dv, dataframe[,factors], sd, na.rm=TRUE)/sqrt(tapply(dv, dataframe[,factors], length))
    lo.bar<-xbars-se
    hi.bar<-xbars+se 
  }
  extras<-list(...)
  if(legend & !percentage){
    if(is.null(legend.text))
      legend.text<-sort(unique(dataframe[[factors[1]]]))
    args.legend.temp<-list(x="topright", bty=if(!legend.border)"n" else "o",
                           inset=c(0,0))
    if(is.list(args.legend))
      args.legend<-modifyList(args.legend.temp, args.legend)
    else 
      args.legend<-args.legend.temp
  } else if(legend & percentage){
    if(is.null(legend.text)) 
      legend.text<-c("1", "0")
    args.legend.temp<-list(x="topright", bty=if(!legend.border)"n" else "o",
                           inset=c(0,0))
    if(is.list(args.legend))
      args.legend<-modifyList(args.legend.temp, args.legend)
    else 
      args.legend<-args.legend.temp
  } else if(!legend){
    args.legend<-NULL
    legend.text<-NULL
  }
  if(errbar && legend && !percentage) ymax<-max(hi.bar)+max(hi.bar)/20
  else if(errbar && legend && percentage) ymax<-115
  else if(errbar && !legend) ymax <- max(xbars)
  else if(!errbar && legend && percentage) ymax<-110  
  else if(!errbar) ymax<-max(xbars) + max(xbars)/20
  if(!percentage){
    args.barplot<-list(beside=TRUE, height=xbars, ylim=c(0, ymax), main=main, names.arg=names.arg,
                       col=hcl(h=seq(0,270, 270/(length(unique(dataframe[[factors[1]]]))))[-length(unique(dataframe[[factors[1]]]))]),
                       legend.text=legend.text, args.legend=args.legend, xpd=TRUE,
                       xlab=if(is.null(xlab)) factors[length(factors)] else xlab,
                       ylab=if(is.null(ylab)) dv.name else ylab, axes=axes)
  }else{
    args.barplot<-list(beside=TRUE, height=xbars, ylim=c(0, ymax),  main=main, names.arg=names.arg,
                       col=hcl(h=seq(0,270, 270/(length(unique(dataframe[[factors[1]]]))))[-length(unique(dataframe[[factors[1]]]))]),
                       legend.text=legend.text, args.legend=args.legend, xpd=TRUE,
                       xlab=if(is.null(xlab)) " "[length(factors)] else xlab,
                       ylab=if(is.null(ylab)) "percentage" else ylab, axes=axes)                        
  }
  args.barplot<-modifyList(args.barplot, extras)
  errbars = function(xvals, cilo, cihi, whisker, nc, args.errbar = NULL, half.errbar=TRUE) {
    if(half.errbar){
      cilo<-(cihi+cilo)/2
    }
    fixedArgs.bar = list(matlines, x=list(xvals), 
                         y=lapply(split(as.data.frame(t(do.call("rbind", 
                                                                list(cihi, cilo)))),1:nc),matrix, 
                                  nrow=2, byrow=T))
    allArgs.bar = c(fixedArgs.bar, args.errbar)
    whisker.len = whisker*(par("usr")[2] - par("usr")[1])/2
    whiskers = rbind((xvals - whisker.len)[1,],
                     (xvals + whisker.len)[1,])
    fixedArgs.lo = list(matlines, x=list(whiskers),   
                        y=lapply(split(as.data.frame(t(do.call("rbind", 
                                                               list(cilo, cilo)))), 1:nc), matrix, nrow=2, byrow=T))
    allArgs.bar.lo = c(fixedArgs.lo, args.errbar)
    fixedArgs.hi = list(matlines, x=list(whiskers), 
                        y=lapply(split(as.data.frame(t(do.call("rbind", 
                                                               list(cihi, cihi)))), 1:nc), matrix, nrow=2, byrow=T))
    allArgs.bar.hi = c(fixedArgs.hi, args.errbar)  
    invisible(do.call(mapply, allArgs.bar))
    if(!half.errbar) invisible(do.call(mapply, allArgs.bar.lo))
    invisible(do.call(mapply, allArgs.bar.hi))
  }
  par(mar=mar)
  errloc<-as.vector(do.call(barplot, args.barplot))
  if(errbar){
    errloc<-rbind(errloc, errloc)
    lo.bar<-matrix(as.vector(lo.bar))
    hi.bar<-matrix(as.vector(hi.bar))
    args.errbar.temp<-list(col=bar.col, lty=1)
    args.errbar<-if(is.null(args.errbar)|!is.list(args.errbar)) 
      args.errbar.temp
    else if(is.list(args.errbar)) 
      modifyList(args.errbar.temp, args.errbar)
    errbars(errloc, cilo=lo.bar, cihi=hi.bar, nc=1, whisker=whisker, 
            args.errbar=args.errbar, half.errbar=half.errbar)
  }
  if(box) box()
  if(percentage){
    args.yaxis.temp<-list(at=seq(0,100, 20), las=1)
    args.yaxis<-if(!is.list(args.yaxis)) args.yaxis.temp else modifyList(args.yaxis.temp, args.yaxis)
    do.call(axis, c(side=2, args.yaxis))
  }
}


lab_sub<- function(dfr)
{
  dfr$mkat[dfr$species == "P. bisulcatum "] <- "C3"
  dfr$mkat[dfr$species == "P. milliodes "] <-"C3-C4"
  dfr$mkat[dfr$species == "P. monticola "] <- "PEP-CK"
  dfr$mkat[dfr$species == "C. ciliaris "] <- "NADP-ME"
  dfr$mkat[dfr$species == "P. coloratum "] <- "NAD-ME"
  dfr$mkat[dfr$species == "E. curvula "] <- "NAD-ME"
  dfr$mkat[dfr$species == "P. maximum "] <- "PEP-CK"
  dfr$mkat[dfr$species == "C. gayana "] <- "PEP-CK"
  dfr$mkat[dfr$species == "S. bicolor "] <- "NADP-ME"
  dfr$mkat[dfr$species == "Z. mays "] <- "NADP-ME"
  return(dfr$mkat)
}

###Tidy anova data output
tidy.anova <- function(x, ...) {
  nn <- c("df", "sumsq", "meansq", "statistic", "p.value")
  ret <- fix_data_frame(x, nn[1:ncol(x)])
  
  if ("term" %in% colnames(ret)) {
    # if rows had names, strip whitespace in them
    ret <- ret %>% mutate(term = stringr::str_trim(term))
  }
  ret
}

tidy.Anova <- function(x, ...) {
  nn <- c( "sumsq", "df","F value" , "p.value")
  ret <- fix_data_frame(x, nn[1:ncol(x)])
  
  if ("term" %in% colnames(ret)) {
    # if rows had names, strip whitespace in them
    ret <- ret %>% mutate(term = stringr::str_trim(term))
  }
  ret
}


tidy.Anova.lmer <- function(x, ...) {
  nn <- c( "Chisq", "df","p.value")
  ret <- fix_data_frame(x, nn[1:ncol(x)])
  
  if ("term" %in% colnames(ret)) {
    # if rows had names, strip whitespace in them
    ret <- ret %>% mutate(term = stringr::str_trim(term))
  }
  ret
}

##tidy ANOVA list

tidy.aovlist <- function(x, ...) {
  # must filter out Intercept stratum since it has no dimensions
  if (names(x)[1L] == "(Intercept)") {
    x <- x[-1L]
  }
  
  ret <- plyr::ldply(x, tidy, .id = "stratum")
  # get rid of leading and trailing whitespace in term and stratum columns
  ret <- ret %>% mutate(term = stringr::str_trim(term),
                        stratum = stringr::str_trim(stratum))
  ret
}



sma_const<- function (mod)
{
  c<-mod$groupsummary
  d<- round(c[, c(2,3,5,8)], 2)
  dat<- cbind(c[,1], d)
  return(dat)
}
