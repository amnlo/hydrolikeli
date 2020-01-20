pred.stats <- function(list.sudriv, auto=NA, time.recess=NA, mu=NA, rep.mu.times=NA, biased=FALSE, tme.orig="1985-01-01", brn.in=0, n.sample=500){
  ## extracts some internal states of the likelihood calculation
  sudriv <- list.sudriv[[1]]
  time <- sudriv$layout$layout$time[c(sudriv$layout$calib, sudriv$layout$pred)] ## ATTENTION: what if we want to sample for calib or valid alone?
  time <- as.POSIXlt(x=tme.orig)+time*60*60*sudriv$layout$timestep.fac*ifelse(sudriv$layout$time.units=="days",24,1)
  dat <- data.frame(x=rep(time, times=length(list.sudriv)), det=NA, quant=NA, white.noise=NA, unifvalue=NA, innovation=NA, obs=NA, lp.num.pred=NA, case = NA, taus=NA)
  i <- 1
  for(sudriv in list.sudriv){
    ## if(!biased){
    ##     sudriv$layout$pred <- c(sudriv$layout$calib , sudriv$layout$pred)
    ##     sudriv$layout$calib <- sudriv$layout$pred
    ## }else{
    ##     sudriv$layout$pred <- sudriv$layout$calib
    ## }
    det    <- c(sampling_wrapper(sudriv, sample.par=FALSE, n.sample=1, sample.likeli=FALSE))
    ## sample.calibpred <- sampling_wrapper(sudriv, brn.in=brn.in, sample.par=TRUE, n.sample=n.sample, biased=biased, sample.likeli=TRUE)
    ##sudriv$layout$calib <- sudriv$layout$pred
    x0 <- c(sudriv$model$parameters[as.logical(sudriv$model$par.fit)], sudriv$likelihood$parameters[as.logical(sudriv$likelihood$par.fit)])
    if(biased) x0 <- c(x0, mu)
    sudriv$layout$calib <- c(sudriv$layout$calib, sudriv$layout$pred)
    ll <-logposterior(x0=x0, sudriv=sudriv, prior=FALSE, verbose=FALSE)
    lp.num.pred <- logposterior(x0=x0, sudriv=sudriv, prior=TRUE, verbose=TRUE)
    dat[((i-1)*length(time)+1):(i*length(time)),] <- data.frame(x=time, det=det, quant=ll$quant, white.noise=ll$white.noise, unifvalue=ll$unifvalue, innovation=ll$innovation, obs=sudriv$observations[c(sudriv$layout$calib)], lp.num.pred=lp.num.pred, case=names(list.sudriv)[i], taus=ifelse(rep(is.null(ll$taus),length(time)), NA, ll$taus))
    ##obs   <- data.frame(x=time, value=sudriv$observations[sudriv$layout$pred])
    i <- i + 1
  }
  return(dat)
}
