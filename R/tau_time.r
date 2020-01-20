tau.time <- function(dat, sudriv, obspred="obspred", calc.S=FALSE){
  ## get the time-course of the time-dependent tau
  par.likeli<- ifelse(as.logical(sudriv$likelihood$tran), exp(sudriv$likelihood$parameters), sudriv$likelihood$parameters)
  names(par.likeli) <- names(sudriv$likelihood$parameters)
  a <- par.likeli["C1Wv_Qstream_a_lik"]
  b <- par.likeli["C1Wv_Qstream_b_lik"]
  c <- par.likeli["C1Wv_Qstream_c_lik"]
  d <- par.likeli["C1Wv_Qstream_d_lik"]
  sd0 <- par.likeli["C1Wv_Qstream_sd01_lik"]
  Q0 <- par.likeli["C1Wv_Qstream_Q01_lik"]
  ttP<- par.likeli["C1Wv_Qstream_ttP_lik"]
  tkP<- par.likeli["C1Wv_Qstream_tkP_lik"]
  tkQ<- par.likeli["C1Wv_Qstream_tkQ_lik"]
  l<- par.likeli["C1Wv_Qstream_l_lik"]
  m<- par.likeli["C1Wv_Qstream_m_lik"]
  tau.max<- par.likeli["C1Wv_Qstream_taumax_lik"]
  tau.min<- par.likeli["C1Wv_Qstream_taumin_lik"]
  pacf.tmp <- numeric(nrow(dat))
  sudriv$layout$calib <- sudriv$layout$calib[!is.na(sudriv$observations[sudriv$layout$calib])]
  sudriv$layout$pred <- sudriv$layout$pred[!is.na(sudriv$observations[sudriv$layout$pred])]
  dt.vec <- diff(sudriv$layout$layout$time[c(sudriv$layout$calib, sudriv$layout$pred)])
  n <- length(dat$det)
  if(grepl("obs", obspred)){
    resid <- (dat$obs - dat$det)/(a*sd0*(dat$det/Q0)^c+b + d*rollmean(c(0,pmax(diff(dat$det),0)),k=2,fill=0))
    f.lik <- function(taulog, x1, x2, dt){sum(dnorm(x2, mean=x1*exp(-dt/exp(taulog)), sd=sqrt(1-exp(-2*dt/exp(taulog))), log=TRUE))}
    ##mat.tau <- matrix(NA, nrow=nrow(dat), ncol=5)
    vec.tau <- numeric(nrow(dat))
    for(i in 3:(nrow(dat)-2)){
      ##pacf.tmp[i] <- optimize(f=f.lik, interval=c(-3,6.2146), maximum=TRUE, x1=resid[i-1], x2=resid[i], dt=dt.vec[i-1])$maximum
      ##mat.tau[row(mat.tau)==(col(mat.tau)+(i-1))] <- optimize(f=f.lik, interval=c(-3,6.9), maximum=TRUE, x1=resid[(i-1):(i+3)], x2=resid[i:(i+4)], dt=dt.vec[(i-1):(i+3)])$maximum
      vec.tau[i] <- optimize(f=f.lik, interval=c(-3,6.9), maximum=TRUE, x1=resid[(i-3):(i+1)], x2=resid[(i-2):(i+2)], dt=dt.vec[(i-3):(i+1)])$maximum
    }
    ##dat$tau.obs <- rowMeans(exp(mat.tau))
    dat$tau.obs <- exp(vec.tau)
    dat$tau.obs.rel <- dat$tau.obs*dat$det/(a*sd0*(dat$det/Q0)^c + b+d*rollmean(c(0,pmax(diff(dat$det),0)),k=2,fill=0))
  }
  if(grepl("pred", obspred)){
    library("deSolve")
    parms <- c(k=1e-1, k2=0.5)
    p <- sudriv$input$inputobs[sudriv$layout$layout$time[c(sudriv$layout$calib, sudriv$layout$pred)], c(1,2)] ## ATTENTION: we assume that the input and the layou time are the same, and that no precipitation was observed in the missing intervals
    print(sum(p[,2]!=0))
    p[,2] <- pmax(rollmean(p[,2], k=5, fill=0), 0)
    pt1 <- as.numeric(p[,1]) - 0.5/24
    pt1 <- rep(pt1[2:(length(pt1))], each=2)
    pt1 <- c(p[1,1] - 0.5/24, pt1, p[nrow(p),1]+0.5/24)
    pt2 <- c(p[1,2], rep(p[2:nrow(p),2], each=2), p[nrow(p),2])
    precip=cbind(pt1, pt2)
    rhs.S <- function(t, y, parms, precip){
      P <- approx(x=precip[,1], y=precip[,2], xout=t)$y
      return(list(P-parms["k"]*y))
    }
    S0 <- as.numeric(mean(p[,2])/(parms["k"]))
    print(S0)
    tme <- proc.time()
    ## if(calc.S){
    ##     S <- ode(y=S0, times=as.numeric(p[,1]), func=rhs.S, parms=parms, precip=precip, method="euler", rtol=1e-2, atol=1e-2)[,2]
    ##     dat$S <- S
    ## }
    S <- rep(NA,nrow(p))
    S[1] <- S0
    for(i in 2:length(S)){
      S[i] <- S[i-1]*exp(-tkP*dt.vec[i-1]) + (p[i-1,2]+p[i,2])/2##*(1-exp(-parms["k"]*(dt.vec[i-1]))) ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
    }
    dat$S <- S
    cat("time solver: ", proc.time() - tme, "\n")
    gamma <- 1
    RC <- 0.6
    ##tau.calc.tmp <- tau.min + (tau.max-tau.min)*exp(-gamma*p[,2]/(parms["k"]*dat$S)-1/gamma*dat$S*parms["k"]/mean(p[,2])/RC)
    dQpos <- pmax(c(0,diff(dat$det)/dt.vec/(dat$det[1:(n-1)] + Q0*a*(b/(1+b))^c)),0)
    ## tau.calc.tmp <- tau.min + (tau.max-tau.min)*exp(-gamma*p[,2]/dat$S-dat$det/Q0*tkQ)
    tau.calc.tmp <- tau.min + (tau.max-tau.min)*exp(-tkP*dQpos^l-tkQ*(dat$det/Q0)^m)
    print(summary(tau.calc.tmp))
    dat$tau.calc <- NA
    dat$tau.calc <- tau.calc.tmp
  }
  return(dat)
}
