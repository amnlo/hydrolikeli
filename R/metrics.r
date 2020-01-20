calc.reliability <- function(Qsim,Qobs){# as proposed by Evin et al. (2014)
  N <- length(Qobs)
  if(any(Qobs<0, na.rm=TRUE)) stop("Qobs was smaller than 0")
  if(dim(Qsim)[2]!=N) stop("dimensions of Qsim and Qobs don't match")
  ec <- apply(Qsim, 2, ecdf)
  omega <- numeric(N)
  for(i in 1:N){
    omega[i] <- ec[[i]](Qobs[i])
  }
  empi <- ecdf(omega)
  devi <- numeric(N)
  for(i in 1:N){
    devi[i] <- ec[[i]](Qobs[i]) - empi(ec[[i]](Qobs[i]))
  }
  reli <- 1-2/N*sum(abs(devi),na.rm=TRUE)
  return(reli)
}
calc.spread <- function(Qsim,Qobs){# as proposed by McInerney et al. (2017)
  if(any(Qobs<0)) stop("Qobs was smaller than 0")
  if(dim(Qsim)[2]!=length(Qobs)) stop("dimensions of Qsim and Qobs don't match")
  spread <- sum(apply(Qsim, 2, sd))/sum(Qobs)
  return(spread)
}
calc.flashiness <- function(x){ ## according to Baker et al. 2004: A new flashiness index: characteristics and applications to midwestern rivers and streams. American Water Association.
  fi <- sum(abs(x[2:length(x)]-x[1:(length(x)-1)]))/sum(x)
  return(fi)
}
calc.dQdt <- function(x, standardize=TRUE){
  if(standardize){stnd=0.5*(x[2:length(x)]+x[1:(length(x)-1)])}else{stnd=1}
  dQdt <- c(NA,diff(x)/stnd)
  return(dQdt)
}
calc.nse <- function(pred,obs){
  if(length(obs)!=length(pred)) stop("non-equal lengths")
  nse <- 1 - sum((obs-pred)^2)/sum((obs-mean(obs))^2)
  return(nse)
}
calc.crps <- function(pred,obs){
  if(length(obs)!=ncol(pred)) stop("non-equal lengths or wrong array dimensions")
  res <- crps_sample(y=c(obs), dat=t(as.matrix(pred)))
  return(res)
}
calc.metrics <- function(sudriv, dat=NA, xlim=NA, file.out=NA, vars=NA, ...){
  ind.sel.tot <- select.ind(sudriv,xlim=xlim,ind.sel=NA)
  ##if(xlim=="calib") ind.sel.tot <- ind.sel.tot[ind.sel.tot %in% sudriv$layout$calib]#ATTENTION: This is necessary since 'select.ind' does not fully consider layout$calib, but just consideres the range() of layout$calib. Changing the function select.ind() would be a major operation, since it is used often.
  if(is.na(vars[1])) vars <- unique(sudriv$layout$layout[ind.sel.tot,"var"]) #ATTENTION: this $layout is hard-coded here (and below), but if xlim="pred", $pred.layout would be the right choice (confidence 60%)
  cat("vars: ",vars,"\n")
  metrics <- matrix(ncol = length(vars), nrow = 13)
  colnames(metrics) <- vars
  rownames(metrics) <- c("reli","spread","spread.parunc","nse.de","nse.med","nse.sd","sferr.det","sferr.med","sferr.sd","flash.det","flash.med","flash.sd","flash.obs")
  for(var.curr in vars){ #loop over the variables in the selected time (xlim)
    ind.sel = ind.sel.tot[sudriv$layout$layout[ind.sel.tot,"var"]==var.curr]
    pl <- subset(sudriv$layout$pred.layout, var==var.curr)
    pr <- sudriv$predicted$sample[,sudriv$layout$pred.layout$var==var.curr]
    if(!is.null(sudriv$predicted$sample.parunc)){
      pru <- sudriv$predicted$sample.parunc[,sudriv$layout$pred.layout$var==var.curr]
    }else{
      pru <- NULL
    }
    ly <- sudriv$layout$layout[ind.sel,]
    ##predobs <- match(paste(pl$var, pl$time), paste(ly$var, ly$time))
    ##obspred <- match(paste(ly$var, ly$time), paste(pl$var, pl$time))
    ##Qsim = sudriv$predicted$sample[,predobs[!is.na(predobs)]]
    lmp <- ifelse(all(is.na(sudriv$layout$lump[ind.sel])), FALSE, TRUE)
    Qsim <- t(apply(pr, 1, function(y,x,xout) approx(x=x,y=y,xout=xout)$y, x=pl$time, xout=ly$time))
    if(!is.null(pru)){
      Qsim.prunc <- t(apply(pru, 1, function(y,x,xout) approx(x=x,y=y,xout=xout)$y, x=pl$time, xout=ly$time))
    }else{
      Qsim.prunc <- NULL
    }
    if(lmp) Qsim <- t(apply(Qsim, 1, function(x) as.numeric(tapply(x, sudriv$layout$lump[ind.sel], mean))))
    if(lmp & !is.null(Qsim.prunc)) Qsim.prunc <- t(apply(Qsim.prunc, 1, function(x) as.numeric(tapply(x, sudriv$layout$lump[ind.sel], mean))))
    cat("Qsim: ", dim(Qsim), "\n")
    print(any(is.na(Qsim)))
    obs = sudriv$observations[ind.sel]
    if(lmp) obs <- as.numeric(tapply(obs, sudriv$layout$lump[ind.sel], mean))
    cat("obs: ", length(obs), "\n")
    flash = c(apply(Qsim,1,calc.flashiness))
    flash.obs = calc.flashiness(obs)
    det = sudriv$predicted$det[1,sudriv$layout$pred.layout$var==var.curr]
    det = approx(x=pl$time, y=det, xout=ly$time)$y
    if(lmp) det <- as.numeric(tapply(det, sudriv$layout$lump[ind.sel], mean))
    cat("det: ", length(det), "\n")
    flash.det = calc.flashiness(det)
    nse.det = 1-sum((det-obs)^2)/sum((obs - mean(obs))^2)
    nse = c(apply(Qsim,1,calc.nse,obs=obs))
    ##crps = c(calc.crps(Qsim,obs=obs))
    ##crps = crps/sudriv$layout$timestep.fac
    strmflw.tot = c(apply(Qsim,1,sum))
    strmflw.err = (sum(obs)-strmflw.tot)/sum(obs)*100
    strmflw.err.det = (sum(obs)-sum(det))/sum(obs)*100
    reliability = calc.reliability(Qsim=Qsim,Qobs=obs)
    likpars = ifelse(as.logical(sudriv$likelihood$tran), exp(sudriv$likelihood$parameters), sudriv$likelihood$parameters)
    names(likpars) = names(sudriv$likelihood$parameters)
    spread = calc.spread(Qsim=Qsim,Qobs=obs)
    if(!is.null(Qsim.prunc)){
      spread.prunc = calc.spread(Qsim=Qsim.prunc,Qobs=obs)
    }else{
      spread.prunc = NULL
    }
    metrics[,var.curr] <- c(reliability,spread,ifelse(is.null(spread.prunc),NA,spread.prunc),nse.det,median(nse),sd(nse),strmflw.err.det,median(strmflw.err),sd(strmflw.err),flash.det,median(flash),sd(flash),flash.obs)
  }
  if(is.na(file.out)){
    sudriv$predicted$metrics <- metrics
    return(sudriv)
  }else{
    metrics <- round(metrics,digits=3)
    write.table(metrics, file=file.out, quote=FALSE)
  }
}
