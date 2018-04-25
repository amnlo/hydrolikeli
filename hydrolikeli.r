## ===================================================================================
## Likelihood Evaluation
## ===================================================================================
## These are functions that evaluate the density of the likelihood (and the prior), e.g. for inference

logposterior <- function(x0, sudriv, prior, apprx=FALSE, verbose=TRUE, biased=FALSE, rep.mu.times=NA, time.recess=NA, auto=NA){
    flp <- sudriv$likelihood$par.fit
    fmp <- sudriv$model$par.fit
    l.fit <- sum(c(fmp,flp))
    if(biased){
        ##overdim <- length(x0) - l.fit
        ##if(overdim %% 2 != 0) stop("the rest of x0 must be of even length")
        mu <- x0[(l.fit+1):length(x0)]
        ##shift <- x0[(l.fit+overdim/2+1):length(x0)]
        x0 <- x0[1:l.fit]
    }else{
        mu <- 1
        rep.mu.times <- length(sudriv$layout$calib)
    }
    if(all(is.na(time.recess))){
        time.recess <- rep(-1, length(sudriv$layout$calib))
    }
    if(all(is.na(auto))){
        auto <- rep(FALSE, length(sudriv$layout$calib))
    }
    ## =======================================================
    ## update the likelihood parameters with the ones from x0
    ## make sure they are within the bounds
    par.lik.fit <- x0[(length(x0)-sum(flp)+1):length(x0)]
    if(!sudriv$settings$OPT.bounded){
        lower_bound <- sudriv$likelihood$lb[as.logical(flp)]
        upper_bound <- sudriv$likelihood$ub[as.logical(flp)]
        par.lik.fit <- constrain_parameters(par.lik.fit, lower_bound, upper_bound)
    }
    ## update likelihood parameters
    sudriv$likelihood$parameters[which(flp != 0)] <- par.lik.fit

    ## =======================================================
    ## update model parameters with the ones from x0
    par.mod.fit <- x0[1:(length(x0)-sum(flp))]
    ## make sure they are within the bounds
    if(!sudriv$settings$OPT.bounded){
        lower_bound <- sudriv$model$args$parLo[as.logical(fmp)]
        upper_bound <- sudriv$model$args$parHi[as.logical(fmp)]
        par.mod.fit <- constrain_parameters(par.mod.fit, lower_bound, upper_bound)
    }
    ## update model parameters
    sudriv$model$parameters[which(fmp != 0)] <- par.mod.fit

    ## =======================================================
    ## prepare arguments for the likelihood
    likeli.args           <- list()
    likeli.args$par.model <- sudriv$model$parameters
    likeli.args$run.model <- run.model
    likeli.args$layout    <- sudriv$layout
    likeli.args$y.obs     <- sudriv$observations
    likeli.args$P         <- sudriv$input$P.roll[sudriv$layout$calib]
    likeli.args$par.likeli<- ifelse(as.logical(sudriv$likelihood$tran), exp(sudriv$likelihood$parameters), sudriv$likelihood$parameters)
    names(likeli.args$par.likeli) <- names(sudriv$likelihood$parameters)
    likeli.args$mu        <- mu
    likeli.args$rep.mu.times <- rep.mu.times
    likeli.args$time.recess <- time.recess
    likeli.args$auto      <- auto
    likeli.args$apprx     <- apprx
    likeli.args$sudriv    <- sudriv
    likeli.args$verbose   <- verbose
    f.likeli <- sudriv$likelihood$f.likeli

    ## =======================================================
    ## calculate logprior
    if(prior){
        pri.m <- sudriv$model$prior
        pri.m$distdef <- pri.m$distdef[as.logical(fmp)]
        pri.l <- sudriv$likelihood$prior
        pri.l$distdef <- pri.l$distdef[as.logical(flp)]
        pri.mu <- sudriv$likelihood$prior
        pri.mu$distdef <- rep(list(c("lognormal", "1", "0.3")), length(mu))
	if(sum(fmp)>0){
	        args.pdf.model       <- c(list(z=as.numeric(sudriv$model$parameters[as.logical(fmp)])), pri.m)
	        logpri.modelpar      <- do.call(calcpdf_mv, args.pdf.model)
                ## cat("logpri.modelpar: ", logpri.modelpar, "\n")
	}else{
		logpri.modelpar <- 0
	}
	if(sum(flp)>0){
	        args.pdf.likeli      <- c(list(z=as.numeric(sudriv$likelihood$parameters[as.logical(flp)])), pri.l)
		logpri.likelipar     <- do.call(calcpdf_mv, args.pdf.likeli)
                ## cat("logpri.likelipar: ", logpri.likelipar, "\n")
	}else{
		logpri.likelipar <- 0
	}
        if(biased){
            args.pdf.mu      <- c(list(z=as.numeric(mu)), pri.mu)
            logpri.likelipar <- logpri.likelipar + do.call(calcpdf_mv, args.pdf.mu)
        }
    }else{
        logpri.likelipar <- 0
        logpri.modelpar  <- 0
    }
    ## =======================================================
    ## calculate loglikelihood
    if(is.finite(logpri.modelpar) & is.finite(logpri.likelipar)){
        tme <- proc.time()
        loglikeli <- do.call(f.likeli, likeli.args)
        ## cat("loglik: ", loglikeli, "\n")
        if(!verbose){return(loglikeli)}
    }else{
        return(-Inf)
    }
    ## =======================================================
    ## calculate logposterior
    logpost <- loglikeli + logpri.likelipar + logpri.modelpar
    return(logpost)
}

logposterior.neg <- function(x0, sudriv, prior, apprx=FALSE){
    ## =======================================================
    ## update the likelihood parameters with the ones from x0
    flp <- sudriv$likelihood$par.fit
    l.fit.lik <- sum(flp)
    par.lik.fit <- x0[(length(x0)-l.fit.lik+1):length(x0)]
    ## make sure they are within the bounds
    if(!sudriv$settings$OPT.bounded){
        lower_bound <- sudriv$likelihood$lb[as.logical(flp)]
        upper_bound <- sudriv$likelihood$ub[as.logical(flp)]
        par.lik.fit <- constrain_parameters(par.lik.fit, lower_bound, upper_bound)
    }
    ## update likelihood parameters
    sudriv$likelihood$parameters[which(flp != 0)] <- par.lik.fit

    ## =======================================================
    ## update model parameters with the ones from x0
    fmp <- sudriv$model$par.fit
    par.mod.fit <- x0[1:(length(x0)-l.fit.lik)]
    ## make sure they are within the bounds
    if(!sudriv$settings$OPT.bounded){
        lower_bound <- sudriv$model$args$parLo[as.logical(fmp)]
        upper_bound <- sudriv$model$args$parHi[as.logical(fmp)]
        par.mod.fit <- constrain_parameters(par.mod.fit, lower_bound, upper_bound)
    }
    ## update model parameters
    sudriv$model$parameters[which(fmp != 0)] <- par.mod.fit

    ## =======================================================
    ## prepare arguments for the likelihood
    likeli.args           <- list()
    likeli.args$par.model <- sudriv$model$parameters
    likeli.args$run.model <- run.model
    likeli.args$layout    <- sudriv$layout
    likeli.args$y.obs     <- sudriv$observations
    likeli.args$par.likeli<- ifelse(as.logical(sudriv$likelihood$tran), exp(sudriv$likelihood$parameters), sudriv$likelihood$parameters)
    names(likeli.args$par.likeli) <- names(sudriv$likelihood$parameters)
    likeli.args$apprx     <- apprx
    likeli.args$sudriv    <- sudriv
    f.likeli <- sudriv$likelihood$f.likeli

    ## =======================================================
    ## calculate logprior
    if(prior){
        pri.m <- sudriv$model$prior
        pri.m$distdef <- pri.m$distdef[as.logical(fmp)]
        pri.l <- sudriv$likelihood$prior
        pri.l$distdef <- pri.l$distdef[as.logical(flp)]
	if(sum(fmp)>0){
	        args.pdf.model       <- c(list(z=as.numeric(sudriv$model$parameters)[as.logical(fmp)]), pri.m)
	        logpri.modelpar      <- do.call(calcpdf_mv, args.pdf.model)
	}else{
		logpri.modelpar <- 0
	}
	if(sum(flp)>0){
	        args.pdf.likeli      <- c(list(z=as.numeric(sudriv$likelihood$parameters[as.logical(flp)])), pri.l)
		logpri.likelipar     <- do.call(calcpdf_mv, args.pdf.likeli)
	}else{
		logpri.likelipar <- 0
	}
    }else{
        logpri.likelipar <- 0
        logpri.modelpar  <- 0
    }

    ## =======================================================
    ## calculate loglikelihood
    if(is.finite(logpri.modelpar) & is.finite(logpri.likelipar)){
        loglikeli <- do.call(f.likeli, likeli.args)
    }else{
        loglikeli <- 0
    }
    ## =======================================================
    ## calculate logposterior
    logpost <- loglikeli + logpri.likelipar + logpri.modelpar
    return(-1*logpost)
}

LogLikelihoodHydrology_mod <- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n)
        scale <- sd[1]/sd.t
        if ( Q.obs[1] > 0 ){
            log.p[1] <- dt((Q.obs[1]-Q.sim[1])/scale,df=df,log=TRUE) - log(scale)}
        else{
            log.p[1] <- pt(-1*Q.sim[1]/scale,df=df,log=TRUE)}
        for ( i in 2:n ){
            if ( Q.sim[i]<Q.sim[i-1] & Q.obs[i-1]>0 ){
                delta <- (t.obs[i]-t.obs[i-1])/tau
                eta.1 <- qnorm(pt((Q.obs[i-1]-Q.sim[i-1])/scale,df=df),mean=0,sd=1)
                scale <- sd[i]/sd.t
                if ( Q.obs[i] > 0 ){
                    eta.2 <- qnorm(pt((Q.obs[i]-Q.sim[i])/scale,df=df),mean=0,sd=1)
                    log.p[i] <- dt((Q.obs[i]-Q.sim[i])/scale,df=df,log=TRUE) - log(scale) +
                        dnorm(eta.2,mean=eta.1*exp(-delta),sd=sqrt(1-exp(-2*delta)),log=TRUE) -
                        dnorm(eta.2,mean=0,sd=1,log=TRUE)
                }else{
                    eta.2 <- qnorm(pt(-1*Q.sim[i]/scale,df=df),mean=0,sd=1)
                    log.p[i] <- pnorm(eta.2,mean=eta.1*exp(-delta),sd=sqrt(1-exp(-2*delta)),log=TRUE)}}
             else{
                 scale <- sd[i]/sd.t
                 if ( Q.obs[i] > 0 ){
                     log.p[i] <- dt((Q.obs[i]-Q.sim[i])/scale,df=df,log=TRUE) - log(scale)}
                 else{
                     log.p[i] <- pt(-1*Q.sim[i]/scale,df=df,log=TRUE)}}
        }
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    ss <- sum(loglikeli)
    return(ss)
}

LogLikelihoodHydrology_mod_fast <- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n-1)
        scale.b <- sd[1:(n-1)]/sd.t
        scale.f <- sd[2:n]/sd.t
        scale.ini <- sd[1]/sd.t
        if ( Q.obs[1] > 0 ){
            log.p.ini <- dt((Q.obs[1]-Q.sim[1])/scale.ini,df=df,log=TRUE) - log(scale.ini)
        }else{
            log.p.ini <- pt(Q.sim[1]/scale.ini,df=df,log=TRUE)
        }
        auto  <- Q.sim[2:n] < Q.sim[1:(n-1)]
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        delta <- (t.obs[2:n] - t.obs[1:(n-1)])/tau
        dq.f  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f
        dq.b  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b
        eta.1 <- qnorm(pt(dq.b,df=df),mean=0,sd=1)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        eta.2.aq <- qnorm(pt(dq.f[ind.auto.qob],df=df),mean=0,sd=1)
        eta.2.aq2<- qnorm(pt(-1*Q.sim[2:n][ind.auto.qob2]/scale.f[ind.auto.qob2],df=df),mean=0,sd=1)
        log.p[ind.auto.qob] <- dt(dq.f[ind.auto.qob],df=df,log=TRUE) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq,mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]),sd=sqrt(1-exp(-2*delta[ind.auto.qob])),log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2,mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])),log=TRUE)
        log.p[!(auto & qob)] <- dt(dq.f[!(auto&qob)],df=df,log=TRUE) - log(scale.f[!(auto&qob)])
        log.p[!(auto & qob) & !qob2] <- pt(-1*Q.sim[2:n][!(auto & qob) & !qob2]/scale.f[!(auto & qob) & !qob2],df=df,log=TRUE)
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    ss <- sum(loglikeli)
    return(ss)
}

LogLikelihoodHydrology_mod_skewt <- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf ){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n)
        scale <- sd[1]/sd.t
        if ( Q.obs[1] > 0 ){
            log.p[1] <- mydskt((Q.obs[1]-Q.sim[1])/scale,df=df,gamma=gamma,log=TRUE) - log(scale)}
        else{
            log.p[1] <- log(pskt(-1*Q.sim[1]/scale,df=df,gamma=gamma))}
        for ( i in 2:n ){
            if ( Q.sim[i]<Q.sim[i-1] & Q.obs[i-1]>0 ){
                delta <- (t.obs[i]-t.obs[i-1])/tau
                eta.1 <- qnorm(pskt((Q.obs[i-1]-Q.sim[i-1])/scale,df=df,gamma=gamma),mean=0,sd=1)
                scale <- sd[i]/sd.t
                if ( Q.obs[i] > 0 ){
                    eta.2 <- qnorm(pskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma),mean=0,sd=1)
                    log.p[i] <- mydskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma,log=TRUE) - log(scale) +
                        dnorm(eta.2,mean=eta.1*exp(-delta),sd=sqrt(1-exp(-2*delta)),log=TRUE) -
                        dnorm(eta.2,mean=0,sd=1,log=TRUE)
                }else{
                    eta.2 <- qnorm(pskt(-1*Q.sim[i]/scale,df=df,gamma=gamma),mean=0,sd=1)
                    log.p[i] <- pnorm(eta.2,mean=eta.1*exp(-delta),sd=sqrt(1-exp(-2*delta)),log=TRUE)}}
             else{
                 scale <- sd[i]/sd.t
                 if ( Q.obs[i] > 0 ){
                     log.p[i] <- mydskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma,log=TRUE) - log(scale)}
                 else{
                     log.p[i] <- log(pskt(-1*Q.sim[i]/scale,df=df,gamma=gamma))}}
        }
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    ss <- sum(loglikeli)
    return(ss)
}

LogLikelihoodHydrology_la <- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n)
        scale <- sd[1]/sd.t
        if ( Q.obs[1] > 0 ){
            log.p[1] <- dt((Q.obs[1]-Q.sim[1])/scale,df=df,log=TRUE) - log(scale)}
        else{
            log.p[1] <- pt(Q.sim[1]/scale,df=df,log=TRUE)}
        for ( i in 2:n ){
            if (Q.obs[i-1]>0 ){
                delta <- (t.obs[i]-t.obs[i-1])/tau
                eta.1 <- qnorm(pt((Q.obs[i-1]-Q.sim[i-1])/scale,df=df),mean=0,sd=1)
                if(is.infinite(eta.1)) return(-Inf)
                scale <- sd[i]/sd.t
                if ( Q.obs[i] > 0 ){
                    eta.2 <- qnorm(pt((Q.obs[i]-Q.sim[i])/scale,df=df),mean=0,sd=1)
                    a <- dnorm(eta.2,mean=eta.1*exp(-delta),sd=sqrt(1-exp(-2*delta)),log=TRUE)
                    log.p[i] <- dt((Q.obs[i]-Q.sim[i])/scale,df=df,log=TRUE) - log(scale) +
                        a - dnorm(eta.2,mean=0,sd=1,log=TRUE)
                }else{
                    eta.2 <- qnorm(pt(-1*Q.sim[i]/scale,df=df),mean=0,sd=1)
                    log.p[i] <- pnorm(eta.2,mean=eta.1*exp(-delta),sd=sqrt(1-exp(-2*delta)),log=TRUE)}}
             else{
                 scale <- sd[i]/sd.t
                 if ( Q.obs[i] > 0 ){
                     log.p[i] <- dt((Q.obs[i]-Q.sim[i])/scale,df=df,log=TRUE) - log(scale)}
                 else{
                     log.p[i] <- pt(-1*Q.sim[i]/scale,df=df,log=TRUE)}}
        }
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    ss <- sum(loglikeli)
    return(ss)
}

LogLikelihoodHydrology_la_fast <- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n-1)
        scale.b <- sd[1:(n-1)]/sd.t
        scale.f <- sd[2:n]/sd.t
        scale.ini <- sd[1]/sd.t
        if ( Q.obs[1] > 0 ){
            log.p.ini <- dt((Q.obs[1]-Q.sim[1])/scale.ini,df=df,log=TRUE) - log(scale.ini)
        }else{
            log.p.ini <- pt(Q.sim[1]/scale.ini,df=df,log=TRUE)
        }
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        delta <- (t.obs[2:n] - t.obs[1:(n-1)])/tau
        dq.f  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f
        dq.b  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b
        eta.1 <- qnorm(pt(dq.b,df=df),mean=0,sd=1)
        ind.auto.qob <- qob & qob2
        ind.auto.qob2<- qob & !qob2
        eta.2.aq <- qnorm(pt(dq.f[ind.auto.qob],df=df),mean=0,sd=1)
        eta.2.aq2<- qnorm(pt(-1*Q.sim[2:n][ind.auto.qob2]/scale.f[ind.auto.qob2],df=df),mean=0,sd=1)
        log.p[ind.auto.qob] <- dt(dq.f[ind.auto.qob],df=df,log=TRUE) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq,mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]),sd=sqrt(1-exp(-2*delta[ind.auto.qob])),log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2,mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])),log=TRUE)
        log.p[!qob] <- dt(dq.f[!qob],df=df,log=TRUE) - log(scale.f[!qob])
        log.p[!qob & !qob2] <- pt(-1*Q.sim[2:n][!qob & !qob2]/scale.f[!qob & !qob2],df=df,log=TRUE)
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    ss <- sum(loglikeli)
    return(ss)
}

LogLikelihoodHydrology_la_fast_skewt<- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf ){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n-1)
        scale.b <- sd[1:(n-1)]/sd.t
        scale.f <- sd[2:n]/sd.t
        scale.ini <- sd[1]/sd.t
        if ( Q.obs[1] > 0 ){
            log.p.ini <- log(mydskt((Q.obs[1]-Q.sim[1])/scale.ini,df=df,gamma=gamma)) - log(scale.ini)
        }else{
            log.p.ini <- log(pskt(Q.sim[1]/scale.ini,df=df,gamma=gamma))
        }
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        delta <- (t.obs[2:n] - t.obs[1:(n-1)])/tau
        dq.f  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f
        dq.b  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b
        eta.1 <- qnorm(pskt(dq.b,df=df,gamma=gamma),mean=0,sd=1)
        ind.auto.qob <- qob & qob2
        ind.auto.qob2<- qob & !qob2
        eta.2.aq <- qnorm(pskt(dq.f[ind.auto.qob],df=df,gamma=gamma),mean=0,sd=1)
        eta.2.aq2<- qnorm(pskt(-1*Q.sim[2:n][ind.auto.qob2]/scale.f[ind.auto.qob2],df=df,gamma=gamma),mean=0,sd=1)
        log.p[ind.auto.qob] <- log(mydskt(dq.f[ind.auto.qob],df=df,gamma=gamma)) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq,mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]),sd=sqrt(1-exp(-2*delta[ind.auto.qob])),log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2,mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])),log=TRUE)
        log.p[!qob] <- log(mydskt(dq.f[!qob],df=df,gamma=gamma)) - log(scale.f[!qob])
        log.p[!qob & !qob2] <- log(pskt(-1*Q.sim[2:n][!qob & !qob2]/scale.f[!qob & !qob2],df=df,gamma=gamma))
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    ss <- sum(loglikeli)
    return(ss)
}


LogLikelihoodHydrology_la2_fast <- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")] + 1 ## ATTENTION: note the +1 here
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n-1)
        scale.b <- sd[1:(n-1)]/sd.t
        scale.f <- sd[2:n]/sd.t
        scale.ini <- sd[1]/sd.t
        if ( Q.obs[1] > 0 ){
            log.p.ini <- dt((Q.obs[1]-Q.sim[1])/scale.ini,df=df,log=TRUE) - log(scale.ini)
        }else{
            log.p.ini <- pt(Q.sim[1]/scale.ini,df=df,log=TRUE)
        }
        auto  <- ((1/(t.obs[2:n]-t.obs[1:(n-1)])) * Q.sim[2:n]/Q.sim[1:(n-1)]) < d
        auto[Q.sim[2:n]==0 & Q.sim[1:(n-1)]==0] <- FALSE
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        delta <- (t.obs[2:n] - t.obs[1:(n-1)])/tau
        dq.f  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f
        dq.b  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b
        eta.1 <- qnorm(pt(dq.b,df=df),mean=0,sd=1)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        eta.2.aq <- qnorm(pt(dq.f[ind.auto.qob],df=df),mean=0,sd=1)
        eta.2.aq2<- qnorm(pt(-1*Q.sim[2:n][ind.auto.qob2]/scale.f[ind.auto.qob2],df=df),mean=0,sd=1)
        log.p[ind.auto.qob] <- dt(dq.f[ind.auto.qob],df=df,log=TRUE) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq,mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]),sd=sqrt(1-exp(-2*delta[ind.auto.qob])),log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2,mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])),log=TRUE)
        log.p[!(auto & qob)] <- dt(dq.f[!(auto&qob)],df=df,log=TRUE) - log(scale.f[!(auto&qob)])
        log.p[!(auto & qob) & !qob2] <- pt(-1*Q.sim[2:n][!(auto & qob) & !qob2]/scale.f[!(auto & qob) & !qob2],df=df,log=TRUE)
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    ss <- sum(loglikeli)
    return(ss)
}

LogLikelihoodHydrology_la2_fast_skewt <- function(par.model, run.model, layout, y.obs, par.likeli, apprx, verbose, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")] + 1 ## ATTENTION: note the +1 here
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf ){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n-1)
        scale.b <- sd[1:(n-1)]/sd.t
        scale.f <- sd[2:n]/sd.t
        scale.ini <- sd[1]/sd.t
        if ( Q.obs[1] > 0 ){
            log.p.ini <- log(mydskt((Q.obs[1]-Q.sim[1])/scale.ini,df=df,gamma=gamma)) - log(scale.ini)
        }else{
            log.p.ini <- log(pskt(Q.sim[1]/scale.ini,df=df,gamma=gamma))
        }
        auto  <- ((1/(t.obs[2:n]-t.obs[1:(n-1)])) * Q.sim[2:n]/Q.sim[1:(n-1)]) < d
        auto[Q.sim[2:n]==0 & Q.sim[1:(n-1)]==0] <- FALSE
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        delta <- (t.obs[2:n] - t.obs[1:(n-1)])/tau
        dq.f  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f
        dq.b  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b
        eta.1 <- qnorm(pskt(dq.b,df=df,gamma=gamma),mean=0,sd=1)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        eta.2.aq <- qnorm(pskt(dq.f[ind.auto.qob],df=df,gamma=gamma),mean=0,sd=1)
        eta.2.aq2<- qnorm(pskt(-1*Q.sim[2:n][ind.auto.qob2]/scale.f[ind.auto.qob2],df=df,gamma=gamma),mean=0,sd=1)
        log.p[ind.auto.qob] <- log(mydskt(dq.f[ind.auto.qob],df=df,gamma=gamma)) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq,mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]),sd=sqrt(1-exp(-2*delta[ind.auto.qob])),log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2,mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])),log=TRUE)
        log.p[!(auto & qob)] <- log(mydskt(dq.f[!(auto&qob)],df=df,gamma=gamma)) - log(scale.f[!(auto&qob)])
        log.p[!(auto & qob) & !qob2] <- log(pskt(-1*Q.sim[2:n][!(auto & qob) & !qob2]/scale.f[!(auto & qob) & !qob2],df=df,gamma=gamma))
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){return(list(smm=sum(loglikeli), loglik=loglikeli, eta=c(eta.1, eta.2.aq[length(eta.2.aq)])))}
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la3_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, apprx, verbose, ...){ ## la3 uses the precipitation as a criterion to decide if autocorrelation is broken or not
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        P.var <- P[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf | gamma != 1){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            ##myfun2 <- function(x,df,scale){2*x/scale*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            ##m3 <- integrate(f=myfun2, lower=0, upper=Inf, df=df, scale=scale)$value
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            ##Ey <- m3 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*(Q.sim/Q0)^c + b
        log.p <- rep(NA,n-2)
        scale.b <- sd[2:(n-1)]/sd.t
        scale.f <- sd[3:n]/sd.t
        scale.ini <- sd[1:2]/sd.t
        log.p.ini <- c(0,0)
        for(i in 1:2){
            if ( Q.obs[i] > 0 ){
                log.p.ini[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale.ini[i],df=df,gamma=gamma)) - log(scale.ini[i])
            }else{
                log.p.ini[i] <- log(pskt(-Q.sim[i]/scale.ini[i],df=df,gamma=gamma))
            }
        }
        auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        auto[Q.obs[3:n]==0 & Q.obs[2:(n-1)]==0 & Q.obs[1:(n-2)]==0] <- FALSE
        auto[Q.sim[3:n]==0 & Q.sim[2:(n-1)]==0 & Q.sim[1:(n-2)]==0] <- FALSE
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[2:(n-1)] > 0
        qob2  <- Q.obs[3:n] > 0
        delta <- (t.obs[3:n] - t.obs[2:(n-1)])/tau
        dq.f  <- (Q.obs[3:n]-Q.sim[3:n])/scale.f
        dq.b  <- (Q.obs[2:(n-1)]-Q.sim[2:(n-1)])/scale.b
        eta.1 <- qnorm(pskt(dq.b,df=df,gamma=gamma),mean=0,sd=1)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        eta.2.aq <- qnorm(pskt(dq.f[ind.auto.qob],df=df,gamma=gamma),mean=0,sd=1)
        eta.2.aq2<- qnorm(pskt(-1*Q.sim[3:n][ind.auto.qob2]/scale.f[ind.auto.qob2],df=df,gamma=gamma),mean=0,sd=1)
        log.p[ind.auto.qob] <- log(mydskt(dq.f[ind.auto.qob],df=df,gamma=gamma)) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq,mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]),sd=sqrt(1-exp(-2*delta[ind.auto.qob])),log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2,mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])),log=TRUE)
        log.p[!(auto & qob)] <- log(mydskt(dq.f[!(auto&qob)],df=df,gamma=gamma)) - log(scale.f[!(auto&qob)])
        log.p[!(auto & qob) & !qob2] <- log(pskt(-1*Q.sim[3:n][!(auto & qob) & !qob2]/scale.f[!(auto & qob) & !qob2],df=df,gamma=gamma))
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){return(list(smm=sum(loglikeli), loglik=loglikeli, eta=c(NA, eta.1, eta.2.aq[length(eta.2.aq)]), auto=c(FALSE,FALSE,auto)))}
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la4_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, apprx, verbose, ...){ ## la4 tries to keep at least some part of the mass balance by setting the mean of the skewed t distribution to the observations.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        P.var <- P[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        myfun1 <- function(x,df){2*x*dt(x,df)}
        m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
        m2 <- df/(df-2)
        Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
        Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
        sd.t <- sqrt(Ex2 - Ex^2)
        sd <- a*sd0*(Q.sim/Q0)^c + b
        ## scale <- sd/sd.t
        ## m3 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value*scale
        ## Ey <- m3 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        log.p <- rep(NA,n-2)
        scale.b <- sd[2:(n-1)]/sd.t
        scale.f <- sd[3:n]/sd.t
        scale.ini <- sd[1:2]/sd.t
        log.p.ini <- c(0,0)
        for(i in 1:2){
            if ( Q.obs[i] > 0 ){
                log.p.ini[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale.ini[i],df=df,gamma=gamma)) - log(scale.ini[i])
            }else{
                log.p.ini[i] <- log(pskt(-Q.sim[i]/scale.ini[i],df=df,gamma=gamma))
            }
        }
        auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        auto[Q.obs[3:n]==0 & Q.obs[2:(n-1)]==0 & Q.obs[1:(n-2)]==0] <- FALSE
        auto[Q.sim[3:n]==0 & Q.sim[2:(n-1)]==0 & Q.sim[1:(n-2)]==0] <- FALSE
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[2:(n-1)] > 0
        qob2  <- Q.obs[3:n] > 0
        delta <- (t.obs[3:n] - t.obs[2:(n-1)])/tau
        dq.f  <- (Q.obs[3:n]-Q.sim[3:n])/scale.f+Ex
        dq.b  <- (Q.obs[2:(n-1)]-Q.sim[2:(n-1)])/scale.b+Ex
        eta.1 <- qnorm(pskt(dq.b,df=df,gamma=gamma),mean=0,sd=1)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        eta.2.aq <- qnorm(pskt(dq.f[ind.auto.qob],df=df,gamma=gamma),mean=0,sd=1)
        eta.2.aq2<- qnorm(pskt((-1*Q.sim[3:n][ind.auto.qob2])/scale.f[ind.auto.qob2]+Ex,df=df,gamma=gamma),mean=0,sd=1)
        log.p[ind.auto.qob] <- log(mydskt(dq.f[ind.auto.qob],df=df,gamma=gamma)) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq,mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]),sd=sqrt(1-exp(-2*delta[ind.auto.qob])),log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2,mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])),log=TRUE)
        log.p[!(auto & qob)] <- log(mydskt(dq.f[!(auto&qob)],df=df,gamma=gamma)) - log(scale.f[!(auto&qob)])
        log.p[!(auto & qob) & !qob2] <- log(pskt((-1*Q.sim[3:n][!(auto & qob) & !qob2])/scale.f[!(auto & qob) & !qob2]+Ex,df=df,gamma=gamma))
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){return(list(smm=sum(loglikeli), loglik=loglikeli, eta=c(NA, eta.1, eta.2.aq[length(eta.2.aq)]), auto=c(FALSE,FALSE,auto)))}
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la5_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, apprx, verbose, ...){ ## la5 tries the cumulative stochastic process for eta.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        sigma.slpe <- par.likeli[paste(var.curr, "_sigmaslpe_lik", sep = "")]
        sigma.min  <- par.likeli[paste(var.curr, "_sigmamin_lik", sep = "")]
        m          <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        P.var <- P[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf | gamma != 1){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            ##myfun2 <- function(x,df,scale){2*x/scale*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            ##m3 <- integrate(f=myfun2, lower=0, upper=Inf, df=df, scale=scale)$value
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            ##Ey <- m3 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*(Q.sim/Q0)^c + b
        log.p <- rep(NA,n)
        a <- log.p
        b <- log.p
        r <- log.p
        scale.ini <- sd[1:2]/sd.t
        log.p.ini <- c(0,0)
        eta <- rep(NA,n)
        state <- rep(NA,n)
        x     <- rep(NA,n)
        for(i in 1:2){
            if ( Q.obs[i] > 0 ){
                eta[i] <- qnorm(pskt((Q.obs[i]-Q.sim[i])/scale.ini[i],df=df,gamma=gamma),mean=0,sd=1)
                log.p.ini[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale.ini[i],df=df,gamma=gamma)) - log(scale.ini[i])
            }else{
                eta[i] <- qnorm(pskt(-1*Q.sim[i]/scale.ini[i],df=df,gamma=gamma),mean=0,sd=1)
                log.p.ini[i] <- log(pskt(-Q.sim[i]/scale.ini[i],df=df,gamma=gamma))
            }
        }
        log.p[1:2] <- log.p.ini
        auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        auto[Q.obs[3:n]==0 & Q.obs[2:(n-1)]==0 & Q.obs[1:(n-2)]==0] <- FALSE
        auto[Q.sim[3:n]==0 & Q.sim[2:(n-1)]==0 & Q.sim[1:(n-2)]==0] <- FALSE
        scale <- sd[1]/sd.t
        state[2] <- 0
        x[2] <- 0
        for ( i in 3:n ){
            state[i] <- state[i-1]+P.var[i-1]-psi*state[i-1]
            sigma.slpe.tmp <- sigma.slpe*state[i] + sigma.min
            if (auto[i-2]){
                end <- i-1
                x[i] <- x[i-1] + m*eta[i-1]
                delta <- (t.obs[i]-t.obs[i-1])/tau
                scale <- sd[i]/sd.t
                if ( Q.obs[i] > 0 ){
                    eta[i] <- qnorm(pskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma),mean=0,sd=1)
                    a[i] <- dnorm(eta[i],mean=eta[i-1]-x[i],sd=sigma.slpe.tmp,log=TRUE)
                    b[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma))
                    r[i] <- dnorm(eta[i],mean=0,sd=1,log=TRUE)
                    log.p[i] <- b[i] - log(scale) + a[i] - r[i]
                }else{
                    eta[i] <- qnorm(pskt(-1*Q.sim[i]/scale,df=df,gamma=gamma),mean=0,sd=1)
                    log.p[i] <- pnorm(eta[i],mean=eta[i-1]-x,sd=sigma.slpe.tmp,log=TRUE)
                }
            }else{
                x[i] <- 0
                scale <- sd[i]/sd.t
                 if ( Q.obs[i] > 0 ){
                     eta[i] <- qnorm(pskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma),mean=0,sd=1)
                     log.p[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma)) - log(scale)}
                 else{
                     eta[i] <- qnorm(pskt(-1*Q.sim[i]/scale,df=df,gamma=gamma),mean=0,sd=1)
                     log.p[i] <- log(pskt(-1*Q.sim[i]/scale,df=df,gamma=gamma))
                 }
             }
        }
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){return(list(smm=sum(loglikeli), loglik=loglikeli, eta=eta, a=a, b=b, r=r, auto=c(FALSE,FALSE,auto)))}
    return(sum(loglikeli))
}

LogLikelihoodHydrology_pr2_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, apprx, verbose, ...){ ## pr2 uses a OU process with mean != 0, inferred from each event separately
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        P.var <- P[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf | gamma != 1){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            ##myfun2 <- function(x,df,scale){2*x/scale*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            ##m3 <- integrate(f=myfun2, lower=0, upper=Inf, df=df, scale=scale)$value
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            ##Ey <- m3 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*(Q.sim/Q0)^c + b
        log.p <- rep(NA,n)
        a <- log.p
        b <- log.p
        r <- log.p
        scale.ini <- sd[1:2]/sd.t
        log.p.ini <- c(0,0)
        eta <- rep(NA,n)
        for(i in 1:2){
            if ( Q.obs[i] > 0 ){
                eta[i] <- qnorm(pskt((Q.obs[i]-Q.sim[i])/scale.ini[i],df=df,gamma=gamma),mean=0,sd=1)
                log.p.ini[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale.ini[i],df=df,gamma=gamma)) - log(scale.ini[i])
            }else{
                eta[i] <- qnorm(pskt(-1*Q.sim[i]/scale.ini[i],df=df,gamma=gamma),mean=0,sd=1)
                log.p.ini[i] <- log(pskt(-Q.sim[i]/scale.ini[i],df=df,gamma=gamma))
            }
        }
        log.p[1:2] <- log.p.ini
        auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        auto[Q.obs[3:n]==0 & Q.obs[2:(n-1)]==0 & Q.obs[1:(n-2)]==0] <- FALSE
        auto[Q.sim[3:n]==0 & Q.sim[2:(n-1)]==0 & Q.sim[1:(n-2)]==0] <- FALSE
        scale <- sd[1]/sd.t
        j <- ifelse(auto[2], 1, 0)
        for ( i in 3:n ){
            if (auto[i-2]){
                if(i > 3){
                    if(!auto[i-3]) j <- j + 1
                }
                delta <- (t.obs[i]-t.obs[i-1])/tau
                scale <- sd[i]/sd.t
                if ( Q.obs[i] > 0 ){
                    eta[i] <- qnorm(pskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma),mean=0,sd=1)
                    a[i] <- dnorm(eta[i], mean=mu[j]+(eta[i-1]-mu[j])*exp(-delta), sd=sqrt(1-exp(-2*delta)),log=TRUE)
                    b[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma))
                    r[i] <- dnorm(eta[i],mean=0,sd=1,log=TRUE)
                    log.p[i] <- b[i] - log(scale) + a[i] - r[i]
                }else{
                    eta[i] <- qnorm(pskt(-1*Q.sim[i]/scale,df=df,gamma=gamma),mean=0,sd=1)
                    log.p[i] <- pnorm(eta[i],mean=mu[j]+(eta[i-1]-mu[j])*exp(-delta),sd=sqrt(1-exp(-2*delta)),log=TRUE)
                }
            }else{
                scale <- sd[i]/sd.t
                 if ( Q.obs[i] > 0 ){
                     eta[i] <- qnorm(pskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma),mean=0,sd=1)
                     log.p[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale,df=df,gamma=gamma)) - log(scale)}
                 else{
                     eta[i] <- qnorm(pskt(-1*Q.sim[i]/scale,df=df,gamma=gamma),mean=0,sd=1)
                     log.p[i] <- log(pskt(-1*Q.sim[i]/scale,df=df,gamma=gamma))
                 }
             }
        }
        if(j != length(mu)){
            cat("length mu: ", length(mu), "\n")
            cat("j: ", j, "\n")
            stop("mu seems to be longer than expected...")
        }
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){return(list(smm=sum(loglikeli), loglik=loglikeli, eta=eta, a=a, b=b, r=r, auto=c(FALSE,FALSE,auto)))}
    return(sum(loglikeli))
}

LogLikelihoodHydrology_pr2_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, apprx, verbose, ...){ ## pr2 uses a OU process with mean != 0, inferred from each event separately
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        P.var <- P[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf | gamma != 1){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            ##myfun2 <- function(x,df,scale){2*x/scale*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            ##m3 <- integrate(f=myfun2, lower=0, upper=Inf, df=df, scale=scale)$value
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            ##Ey <- m3 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*(Q.sim/Q0)^c + b
        log.p <- rep(NA,n-2)
        scale.b <- sd[2:(n-1)]/sd.t
        scale.f <- sd[3:n]/sd.t
        scale.ini <- sd[1:2]/sd.t
        log.p.ini <- c(0,0)
        for(i in 1:2){
            if ( Q.obs[i] > 0 ){
                log.p.ini[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale.ini[i],df=df,gamma=gamma)) - log(scale.ini[i])
            }else{
                log.p.ini[i] <- log(pskt(-Q.sim[i]/scale.ini[i],df=df,gamma=gamma))
            }
        }
        auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        auto[Q.obs[3:n]==0 & Q.obs[2:(n-1)]==0 & Q.obs[1:(n-2)]==0] <- FALSE
        auto[Q.sim[3:n]==0 & Q.sim[2:(n-1)]==0 & Q.sim[1:(n-2)]==0] <- FALSE
        mu.tmp <- rep(mu, times=rep.mu.times)
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[2:(n-1)] > 0
        qob2  <- Q.obs[3:n] > 0
        delta <- (t.obs[3:n] - t.obs[2:(n-1)])/tau
        dq.f  <- (Q.obs[3:n]-Q.sim[3:n])/scale.f
        dq.b  <- (Q.obs[2:(n-1)]-Q.sim[2:(n-1)])/scale.b
        eta.1 <- qnorm(pskt(dq.b,df=df,gamma=gamma),mean=0,sd=1)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        eta.2.aq <- qnorm(pskt(dq.f[ind.auto.qob],df=df,gamma=gamma),mean=0,sd=1)
        eta.2.aq2<- qnorm(pskt(-1*Q.sim[3:n][ind.auto.qob2]/scale.f[ind.auto.qob2],df=df,gamma=gamma),mean=0,sd=1)
        log.p[ind.auto.qob] <- log(mydskt(dq.f[ind.auto.qob],df=df,gamma=gamma)) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq, mean=mu.tmp[ind.auto.qob]+(eta.1[ind.auto.qob]-mu.tmp[ind.auto.qob])*exp(-delta[ind.auto.qob]),sd=sqrt(1-exp(-2*delta[ind.auto.qob])),log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2, mean=mu.tmp[ind.auto.qob2]+(eta.1[ind.auto.qob2]-mu.tmp[ind.auto.qob2])*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])),log=TRUE)
        log.p[!(auto & qob)] <- log(mydskt(dq.f[!(auto&qob)],df=df,gamma=gamma)) - log(scale.f[!(auto&qob)])
        log.p[!(auto & qob) & !qob2] <- log(pskt(-1*Q.sim[3:n][!(auto & qob) & !qob2]/scale.f[!(auto & qob) & !qob2],df=df,gamma=gamma))
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){return(list(smm=sum(loglikeli), loglik=loglikeli, eta=c(NA, eta.1, eta.2.aq[length(eta.2.aq)]), auto=c(FALSE,FALSE,auto), mu = c(mu.tmp[1], mu.tmp[1], mu.tmp)))}
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la6_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, apprx, verbose, ...){ ## la6 assumes that the quantiles can be described by an damped harmonic oscillator + OU-deviations
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        freq <- 2*pi/par.likeli[paste(var.curr, "_wavlen_lik", sep = "")]
        damp   <- par.likeli[paste(var.curr, "_damp_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        P.var <- P[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf | gamma != 1){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            ##myfun2 <- function(x,df,scale){2*x/scale*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            ##m3 <- integrate(f=myfun2, lower=0, upper=Inf, df=df, scale=scale)$value
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            ##Ey <- m3 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*(Q.sim/Q0)^c + b
        log.p <- rep(NA,n-2)
        scale.b <- sd[2:(n-1)]/sd.t
        scale.f <- sd[3:n]/sd.t
        delta <- (t.obs[3:n] - t.obs[2:(n-1)])/tau
        scale.ini <- sd[1:2]/sd.t
        log.p.ini <- c(0,0)
        for(i in 1:2){
            if ( Q.obs[i] > 0 ){
                log.p.ini[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale.ini[i],df=df,gamma=gamma)) - log(scale.ini[i])
            }else{
                log.p.ini[i] <- log(pskt(-Q.sim[i]/scale.ini[i],df=df,gamma=gamma))
            }
        }
        auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        auto[Q.obs[3:n]==0 & Q.obs[2:(n-1)]==0 & Q.obs[1:(n-2)]==0] <- FALSE
        auto[Q.sim[3:n]==0 & Q.sim[2:(n-1)]==0 & Q.sim[1:(n-2)]==0] <- FALSE
        mu.tmp <- rep(mu, times=rep.mu.times)
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[2:(n-1)] > 0
        qob2  <- Q.obs[3:n] > 0
        dq.f  <- (Q.obs[3:n]-Q.sim[3:n])/scale.f
        dq.b  <- (Q.obs[2:(n-1)]-Q.sim[2:(n-1)])/scale.b
        eta.1 <- qnorm(pskt(dq.b,df=df,gamma=gamma),mean=0,sd=1) - mu.tmp*exp(-damp*time.recess)*cos(freq*time.recess)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        quant.2.aq <- qnorm(pskt(dq.f[ind.auto.qob],df=df,gamma=gamma),mean=0,sd=1)
        eta.2.aq <-  quant.2.aq - mu.tmp[ind.auto.qob]*exp(-damp*time.recess[ind.auto.qob])*cos(freq*time.recess[ind.auto.qob])
        eta.2.aq2<- qnorm(pskt(-1*Q.sim[3:n][ind.auto.qob2]/scale.f[ind.auto.qob2],df=df,gamma=gamma),mean=0,sd=1) - mu.tmp[ind.auto.qob2]*exp(-damp*time.recess[ind.auto.qob2])*cos(freq*time.recess[ind.auto.qob2])

        log.p[ind.auto.qob] <- log(mydskt(dq.f[ind.auto.qob],df=df,gamma=gamma)) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq, mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]), sd=sqrt(1-exp(-2*delta[ind.auto.qob])), log=TRUE) - dnorm(quant.2.aq,mean=0,sd=1,log=TRUE)

        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2, mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])), log=TRUE)

        log.p[!(auto & qob)] <- log(mydskt(dq.f[!(auto&qob)],df=df,gamma=gamma)) - log(scale.f[!(auto&qob)])
        log.p[!(auto & qob) & !qob2] <- log(pskt(-1*Q.sim[3:n][!(auto & qob) & !qob2]/scale.f[!(auto & qob) & !qob2],df=df,gamma=gamma))
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){
        quant <- rep(NA, length(loglikeli)-2)
        quant[ind.auto.qob] <- quant.2.aq
        quant <- c(0,0,quant)
        return(list(smm=sum(loglikeli), loglik=loglikeli, eta=c(NA, eta.1, eta.2.aq[length(eta.2.aq)]), quant=quant, auto=c(FALSE,FALSE,auto), mu = c(mu.tmp[1], mu.tmp[1], mu.tmp)))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la7_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## la7 puts a deterministic exponential decay with a certain shift on top of the deterministic hydrological model. The quantiles are assumed to be OU-deviations from that combination.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        freq <- 2*pi/par.likeli[paste(var.curr, "_wavlen_lik", sep = "")]
        damp   <- par.likeli[paste(var.curr, "_damp_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        P.var <- P[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
##        auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        auto[Q.obs==0] <- FALSE
        auto <- auto[3:n]
        ##mu.tmp <- rep(mu*Q.sim[3:n][which(time.recess == 0)], times=rep.mu.times)
        mu.tmp <- rep(mu, times=rep.mu.times)
        ##shift.tmp <- rep(shift, times=rep.mu.times)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf | gamma != 1){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            ##myfun2 <- function(x,df,scale){2*x/scale*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            ##m3 <- integrate(f=myfun2, lower=0, upper=Inf, df=df, scale=scale)$value
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            ##Ey <- m3 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
##        sd <- a*sd0*(Q.sim/Q0)^c + b
        ##exp.det <- (c(mu.tmp[1],mu.tmp) - c(shift.tmp[1], shift.tmp))*exp(-damp*c(0,time.recess)) + c(shift.tmp[1], shift.tmp)
        exp.det <- mu.tmp*Q.sim
        exp.det[!is.finite(time.recess)] <- 0
        Q.sim <- Q.sim + exp.det
        ##Q.sim <- pmax(Q.sim, 0)
        if(any(Q.sim<=0)) return(-Inf)
        sd <- a*sd0*(Q.sim/Q0)^c + b
##        if(any(!is.finite(sd))){"sd is strange"); stop()}
        log.p <- rep(NA,n-2)
        scale.b <- sd[2:(n-1)]/sd.t
        scale.f <- sd[3:n]/sd.t
        delta <- (t.obs[3:n] - t.obs[2:(n-1)])/tau
        scale.ini <- sd[1:2]/sd.t
        log.p.ini <- c(0,0)
        for(i in 1:2){
            if ( Q.obs[i] > 0 ){
                log.p.ini[i] <- log(mydskt((Q.obs[i]-Q.sim[i])/scale.ini[i],df=df,gamma=gamma)) - log(scale.ini[i])
            }else{
                log.p.ini[i] <- log(pskt(-Q.sim[i]/scale.ini[i],df=df,gamma=gamma))
            }
        }
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[2:(n-1)] > 0
        qob2  <- Q.obs[3:n] > 0
        dq.f  <- (Q.obs[3:n]-Q.sim[3:n])/scale.f
        dq.b  <- (Q.obs[2:(n-1)]-Q.sim[2:(n-1)])/scale.b
        eta.1 <- qnorm(pskt(dq.b,df=df,gamma=gamma),mean=0,sd=1)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        eta.2.aq <- qnorm(pskt(dq.f[ind.auto.qob],df=df,gamma=gamma),mean=0,sd=1)
        ##eta.2.aq <-  quant.2.aq - mu.tmp[ind.auto.qob]*exp(-damp*time.recess[ind.auto.qob])*cos(freq*time.recess[ind.auto.qob])
        eta.2.aq2<- qnorm(pskt(-1*Q.sim[3:n][ind.auto.qob2]/scale.f[ind.auto.qob2],df=df,gamma=gamma),mean=0,sd=1)

        log.p[ind.auto.qob] <- log(mydskt(dq.f[ind.auto.qob],df=df,gamma=gamma)) - log(scale.f[ind.auto.qob]) +
            dnorm(eta.2.aq, mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]), sd=sqrt(1-exp(-2*delta[ind.auto.qob])), log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)

        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2, mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])), log=TRUE)

        log.p[!(auto & qob)] <- log(mydskt(dq.f[!(auto&qob)],df=df,gamma=gamma)) - log(scale.f[!(auto&qob)])
        log.p[!(auto & qob) & !qob2] <- log(pskt(-1*Q.sim[3:n][!(auto & qob) & !qob2]/scale.f[!(auto & qob) & !qob2],df=df,gamma=gamma))
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){
        return(list(smm=sum(loglikeli), loglik=loglikeli, eta=c(NA, eta.1, eta.2.aq[length(eta.2.aq)]), exp.det=exp.det, auto=c(FALSE,FALSE,auto), mu = mu.tmp))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la8_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## la8 assumes that the quantiles follow an integrated OU-process, where the mean of the OU-process is depending on the quantile. The quantile has another stochastic component depending on the rainfall.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        sd01<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tauslpe_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        mn.rev <- par.likeli[paste(var.curr, "_mnrev_lik", sep = "")]
        sigmaslpe <- par.likeli[paste(var.curr, "_sigmaslpe_lik", sep = "")]
        sigmaquant     <- par.likeli[paste(var.curr, "_sdquant_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[2:(n-1)]/sd.t1
        scale.b2 <- sd2[2:(n-1)]/sd.t2
        scale.f1 <- sd1[3:n]/sd.t1
        scale.f2 <- sd2[3:n]/sd.t2
        delta1 <- (t.obs[3:n] - t.obs[2:(n-1)])/tau1
        delta2 <- (t.obs[3:n] - t.obs[2:(n-1)])/tau2
        scale.ini <- ifelse(auto[1:2], sd2[1:2]/sd.t2, sd1[1:2]/sd.t1)
        log.p.ini <- c(0,0)
        for(i in 1:2){
            df.tmp <- ifelse(auto[i], df2, df1)
            gamma.tmp <- ifelse(auto[i], gamma2, gamma1)
            if ( Q.obs[i] > 0 ){
                log.p.ini[i] <- mydskt((Q.obs[i]-Q.sim[i])/scale.ini[i]+Exb2,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini[i])
            }else{
                log.p.ini[i] <- mypskt(-Q.sim[i]/scale.ini[i]+Exb2,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        }
        auto <- auto[3:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[2:(n-1)] > 0
        qob2  <- Q.obs[3:n] > 0
        dq.f1  <- (Q.obs[3:n]-Q.sim[3:n])/scale.f1+Exb1
        dq.b1  <- (Q.obs[2:(n-1)]-Q.sim[2:(n-1)])/scale.b1+Exb1
        dq.f2  <- (Q.obs[3:n]-Q.sim[3:n])/scale.f2+Exb2
        dq.b2  <- (Q.obs[2:(n-1)]-Q.sim[2:(n-1)])/scale.b2+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        mu <- -mn.rev*quant.b2
        eta.b <- c(0,diff(quant.b2))
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        eta.f <- c(eta.b[2], diff(quant.f2))

        a <- dnorm(eta.f[ind.auto.qob], mean=mu[ind.auto.qob] + (eta.b[ind.auto.qob] - mu[ind.auto.qob])*exp(-delta2[ind.auto.qob]), sd=sigmaslpe*sqrt(1-exp(-2*delta2[ind.auto.qob])), log=TRUE)
        b <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sigmaquant*sqrt(1-exp(-2*delta1[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + a - dnorm(eta.f[ind.auto.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.f[ind.auto.qob2], mean=eta.b[ind.auto.qob2]*exp(-delta2[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta2[ind.auto.qob2])), log=TRUE)

        log.p[ind.auto2.qob] <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob]) + b - dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-delta1[ind.auto2.qob2]),sd=sqrt(1-exp(-2*delta1[ind.auto2.qob2])), log=TRUE) ## develop: this is not correct yet...
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        return(list(smm=sum(loglikeli), loglik=loglikeli, eta=c(NA, eta.b, eta.f[length(eta.f)]), quant=c(NA, ifelse(auto, quant.b2, quant.b1), quant.f2[length(quant.f2)]), auto=c(FALSE,FALSE,auto), a=a, b=b))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la9b_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## reparametrized version of la9, where tau is normalized by Qdet/sd
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        p1  <- par.likeli[paste(var.curr, "_p_lik", sep = "")]
        ## p2  <- par.likeli[paste(var.curr, "_p2_lik", sep = "")]
        p2 <- p1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ##tau1 <--1/log(par.likeli[paste(var.curr, "_ar1_lik", sep = "")])
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        ## gamma2 <- gamma1
        tau.max <- 20
        tau.min <- 0
        k       <- 0.1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        Q.sim <- pmax(Q.sim , 0)
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 ##+ pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd1 <- rollmean(sd1, k=2, fill=sd1[length(sd1)])
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd2 <- rollmean(sd2, k=2, fill=sd2[length(sd2)])
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        S <- rep(NA,n)
        RC <- 0.6
        S0 <- mean(P)*RC/(k)
        S[1] <- S0
        for(i in 2:n){
             S[i] <- S[i-1]*exp(-k*(t.obs[i]-t.obs[i-1])) + (P[i-1]+P[i])/2 ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
        }
        alpha <- 1
        S <- pmax(S,0)
        taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-1/alpha*S*k/mean(P)/RC)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-1/alpha*sqrt(S/S0))
        gmms <- 1/taus
        delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ## mns <- quant.b1[1]*exp(-delta1)
        ## sds <- sqrt(1-exp(-2*delta1))

        a <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])), log=TRUE)
        b <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + a - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        ## ind2 <- auto & qob2 & time.recess[2:n]==0
        ## log.p[ind2] <- mydskt(dq.f2[ind2],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind2])
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-delta2[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta2[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f2[ind.auto.qob3],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        log.p[ind.auto2.qob] <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob]) + b - dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        ## ind1 <- !auto & qob2 & time.recess[2:n]==-Inf
        ## log.p[ind1] <- mydskt(dq.f1[ind1],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind1])
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-delta1[ind.auto2.qob2]),sd=sqrt(1-exp(-2*delta1[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- mydskt(dq.f1[ind.auto2.qob3],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob3])
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])))
        unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1)
        unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2)
        unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1,gamma=gamma1)
        unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), auto=c(FALSE,auto), a=a, b=b, sd=sd1))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la9c_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## reparametrized version of la9, where tau is normalized by Qdet/sd
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        p1  <- par.likeli[paste(var.curr, "_p_lik", sep = "")]
        ## p2  <- par.likeli[paste(var.curr, "_p2_lik", sep = "")]
        p2 <- p1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ##tau1 <--1/log(par.likeli[paste(var.curr, "_ar1_lik", sep = "")])
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ttP<- par.likeli[paste(var.curr, "_ttP_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ##tau.max <- 20
        ##tau.min <- 0
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                ## Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1) ## if we use this line, sd.t1 is not the "real" SD of the final distribution DQ, since we do not consider the influence of gamma on SD in the correct way.
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
                m1 <- 1
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                ## Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1) ## if we use this line, sd.t2 is not the "real" SD of the final distribution DQ, since we do not consider the influence of gamma on SD in the correct way.
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
                m2 <- 1
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        Q.sim <- pmax(Q.sim , 0)
        sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ##sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd1 <- rollmean(sd1, k=2, fill=sd1[length(sd1)])
        ##sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd2 <- rollmean(sd2, k=2, fill=sd2[length(sd2)])
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        S <- rep(NA,n)
        RC <- 0.6
        S0 <- mean(P)*RC/(tkP)
        S[1] <- S0
        ## for(i in 2:n){
        ##      S[i] <- S[i-1]*exp(-tkP*(t.obs[i]-t.obs[i-1])) + (P[i-1]+P[i])/2 ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
        ## }
        alpha <- 1
        S <- pmax(S,0)
        taus <- numeric(n)
        ## dQsim <- c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)]))
        ## dQpos <- pmax(rollmean(dQsim, k=3, fill=0), 0)
        taus <- tau.min + (tau.max-tau.min)*exp(-pmax(P-Plim,0)/Q.sim*tkP)
        taus <- rollmean(taus, k=3,fill=0)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-(pmax(P-ttP,0))/Q.sim*tkP-Q.sim/Q01*tkQ)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-1/alpha*sqrt(S/S0))
        gmms <- 1/taus
        dt <- (t.obs[2:n]-t.obs[1:(n-1)])
        sgmm <- gmms[2:n] + gmms[1:(n-1)]
        ## delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        ## delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1##+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1##+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2##+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2##+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ## mns <- quant.b1[1]*exp(-delta1)
        ## sds <- sqrt(1-exp(-2*delta1))

        densEta <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-sgmm[ind.auto.qob]/2*dt[ind.auto.qob]), sd=sqrt(1-exp(-sgmm[ind.auto.qob]*dt[ind.auto.qob])), log=TRUE)
        densEtb <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-sgmm[ind.auto2.qob]/2*dt[ind.auto2.qob]), sd=sqrt(1-exp(-sgmm[ind.auto2.qob]*dt[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + densEta - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        ## ind2 <- auto & qob2 & time.recess[2:n]==0
        ## log.p[ind2] <- mydskt(dq.f2[ind2],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind2])
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-sgmm[ind.auto.qob2]/2*dt[ind.auto.qob2]),sd=sqrt(1-exp(-sgmm[ind.auto.qob2]*dt[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f2[ind.auto.qob3],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        c <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,log=TRUE)  - log(scale.f1[ind.auto2.qob])
        d <- dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto2.qob] <- c + densEtb - d
        ## ind1 <- !auto & qob2 & time.recess[2:n]==-Inf
        ## log.p[ind1] <- mydskt(dq.f1[ind1],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind1])
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-sgmm[ind.auto2.qob2]/2*dt[ind.auto2.qob2]),sd=sqrt(1-exp(-sgmm[ind.auto2.qob2]*dt[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- mydskt(dq.f1[ind.auto2.qob3],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob3])
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(any(loglikeli==Inf,na.rm=T)){warning("infinite loglikeli encountered. Returning -Inf"); return(-Inf)}
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-sgmm[ind.auto2.qob]/2*dt[ind.auto2.qob]), sd=sqrt(1-exp(-sgmm[ind.auto2.qob]*dt[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-sgmm[ind.auto.qob]/2*dt[ind.auto.qob]), sd=sqrt(1-exp(-sgmm[ind.auto.qob]*dt[ind.auto.qob])))
        unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1)
        unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2)
        unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1,gamma=gamma1)
        unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), auto=c(FALSE,auto), a=densEta, b=densEtb, c=c, d=d, sd=sd1))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la9d_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## reparametrized version of la9, where tau is normalized by Qdet/sd
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        p1  <- par.likeli[paste(var.curr, "_p_lik", sep = "")]
        ## p2  <- par.likeli[paste(var.curr, "_p2_lik", sep = "")]
        p2 <- p1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ##tau1 <--1/log(par.likeli[paste(var.curr, "_ar1_lik", sep = "")])
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ttP<- par.likeli[paste(var.curr, "_ttP_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ##tau.max <- 20
        ##tau.min <- 0
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                ## Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1) ## if we use this line, sd.t1 is not the "real" SD of the final distribution DQ, since we do not consider the influence of gamma on SD in the correct way.
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
                m1 <- 1
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                ## Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1) ## if we use this line, sd.t2 is not the "real" SD of the final distribution DQ, since we do not consider the influence of gamma on SD in the correct way.
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
                m2 <- 1
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        Q.sim <- pmax(Q.sim , 0)
        ## dQsim <- pmax(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])),0)
        ## dQpos1 <- pmax(rollmean(dQsim, k=3, fill=0, align="left"), 0)
        S <- numeric(length(Q.sim))
        S[1] <- mean(P)*0.6/psi
        P <- pmax(P,0)
        for(i in 2:n){
            S[i] <- S[i-1]*exp(-psi*(t.obs[i]-t.obs[i-1])) + P[i] ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
            ## S[i] <- S[i-1]/(1+S[i-1]*psi*(t.obs[i]-t.obs[i-1])) + P[i] ## in case alpha = 2
        }
        S <- pmax(S,0)
        sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1 + ifelse(P>0,d*P^e,0)##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2 + ifelse(P>0,d*P^e,0)##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ##sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd1 <- rollmean(sd1, k=2, fill=sd1[length(sd1)])
        ##sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd2 <- rollmean(sd2, k=2, fill=sd2[length(sd2)])
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        ## dQpos2 <- abs(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])/(0.5*Q.sim[1:(n-1)] + 0.5*Q.sim[2:n] + psi)))
        ## dQpos2 <- rollmean(dQpos2, k=3, fill=0, align="left")
        P.rel <- ifelse(P==0, 0, P/S/psi)
        taus <- tau.min + (tau.max-tau.min)*exp(-tkP*P.rel^l - tkQ*(S*psi)^m)
        ## taus <- rollmean(taus, k=3,fill=0)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-(pmax(P-ttP,0))/Q.sim*tkP-Q.sim/Q01*tkQ)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-1/alpha*sqrt(S/S0))
        gmms <- 1/taus
        dt <- (t.obs[2:n]-t.obs[1:(n-1)])
        ## sgmm <- gmms[2:n] + gmms[1:(n-1)]
        gmms <- gmms[2:n]
        ## delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        ## delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1##+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1##+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2##+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2##+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ## mns <- quant.b1[1]*exp(-delta1)
        ## sds <- sqrt(1-exp(-2*delta1))

        densEta <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])), log=TRUE)
        densEtb <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + densEta - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        ## ind2 <- auto & qob2 & time.recess[2:n]==0
        ## log.p[ind2] <- mydskt(dq.f2[ind2],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind2])
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-gmms[ind.auto.qob2]*dt[ind.auto.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto.qob2]*dt[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f2[ind.auto.qob3],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        c <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,log=TRUE)  - log(scale.f1[ind.auto2.qob])
        d <- dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto2.qob] <- c + densEtb - d
        ## ind1 <- !auto & qob2 & time.recess[2:n]==-Inf
        ## log.p[ind1] <- mydskt(dq.f1[ind1],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind1])
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-gmms[ind.auto2.qob2]*dt[ind.auto2.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto2.qob2]*dt[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- mydskt(dq.f1[ind.auto2.qob3],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob3])
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(any(loglikeli==Inf,na.rm=T)){warning("infinite loglikeli encountered. Returning -Inf"); return(-Inf)}
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])))
        unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1)
        unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2)
        unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1,gamma=gamma1)
        unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), auto=c(FALSE,auto), a=densEta, b=densEtb, sd=sd1, taus=taus, Qs=S*psi))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la9e_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## reparametrized version of la9, where tau is normalized by Qdet/sd
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd01 <- Q01/5 ## in case Q01 is fitted
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ##tau1 <--1/log(par.likeli[paste(var.curr, "_ar1_lik", sep = "")])
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ##tau.max <- 20
        ##tau.min <- 0
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        ## P <- pmax(P,0)
        ## gamma1 <- pmin(gamma1+tkQ*P^m,5)
        if ( any(df1 < Inf) | any(gamma1 != 1) | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(any(df1<Inf)){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
                m1 <- 1
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
                m2 <- 1
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
            ## Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1) ## for further use, this is the mean of the unscaled (original) skewed t-distribution
            ## Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            Exb1 <- 0 ## ATTENTION: this means that Qdet is at the mode, not the mean...
            Exb2 <- 0
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        Q.sim <- pmax(Q.sim , 0)
        ## dQsim <- pmax(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])),0)
        ## dQpos1 <- pmax(rollmean(dQsim, k=3, fill=0, align="left"), 0)

        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + Q01*b1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + Q02*b2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*Q.sim[Q.sim<Q01]/Q01 + b1*Q01
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01/c1*(c1-1+(Q.sim[Q.sim>=Q01]/Q01)^c1)
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*Q.sim[Q.sim<Q02]/Q02 + b2*Q02
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02/c2*(c2-1+(Q.sim[Q.sim>=Q02]/Q02)^c2)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*(Q.sim[Q.sim<Q01]/Q01)^c1 + b1
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01*(1-c1^2+c1^2*(Q.sim[Q.sim>=Q01]/Q01)^(1/c1))
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*(Q.sim[Q.sim<Q02]/Q02)^c2 + b2
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02*(1-c2^2+c2^2*(Q.sim[Q.sim>=Q02]/Q02)^(1/c2))

        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        ## gamma1 <- gamma1[2:n]
        ## dQpos2 <- abs(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])/(0.5*Q.sim[1:(n-1)] + 0.5*Q.sim[2:n] + psi)))
        ## dQpos2 <- rollmean(dQpos2, k=3, fill=0, align="left")
        P <- pmax(P,0)
        taus <- tau.min + (tau.max-tau.min)*exp(-tkP*P^l - tkQ*Q.sim^m)
        ## taus[((1:length(taus)) %% 20) == 0] <- 0
        ## taus <- rollmean(taus, k=3,fill=0)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-(pmax(P-ttP,0))/Q.sim*tkP-Q.sim/Q01*tkQ)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-1/alpha*sqrt(S/S0))
        gmms <- 1/taus
        dt <- (t.obs[2:n]-t.obs[1:(n-1)])
        gmms <- gmms[2:n]
        ## delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        ## delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2[1], gamma1[1])
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ## mns <- quant.b1[1]*exp(-delta1)
        ## sds <- sqrt(1-exp(-2*delta1))

        densEta <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])), log=TRUE)
        densEtb <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + densEta - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        ## ind2 <- auto & qob2 & time.recess[2:n]==0
        ## log.p[ind2] <- mydskt(dq.f2[ind2],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind2])
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-gmms[ind.auto.qob2]*dt[ind.auto.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto.qob2]*dt[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f2[ind.auto.qob3],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        c <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,log=TRUE)  - log(scale.f1[ind.auto2.qob])
        d <- dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto2.qob] <- c + densEtb - d
        ## ind1 <- !auto & qob2 & time.recess[2:n]==-Inf
        ## log.p[ind1] <- mydskt(dq.f1[ind1],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind1])
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-gmms[ind.auto2.qob2]*dt[ind.auto2.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto2.qob2]*dt[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- mydskt(dq.f1[ind.auto2.qob3],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob3])
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(any(loglikeli==Inf,na.rm=T)){warning("infinite loglikeli encountered. Returning -Inf"); return(-Inf)}
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        white.noise1 <- (quant.f1 - quant.b1*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        white.noise2 <- (quant.f2 - quant.b2*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])))
        unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1)
        unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2)
        unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1,gamma=gamma1[length(gamma1)])
        unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), white.noise=c(ifelse(auto, white.noise2, white.noise1), ifelse(auto[length(auto)], white.noise2[length(white.noise2)], white.noise1[length(white.noise1)])), auto=c(FALSE,auto), a=densEta, b=densEtb, sd=sd1, taus=taus))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la9esimp_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## simplified version of la9(e), where we first rescale, and then skew the distribution DQ.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd01 <- Q01/5 ## in case Q01 is fitted
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ##tau1 <--1/log(par.likeli[paste(var.curr, "_ar1_lik", sep = "")])
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ##tau.max <- 20
        ##tau.min <- 0
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        ## P <- pmax(P,0)
        ## gamma1 <- pmin(gamma1+tkQ*P^m,5)
        if ( any(df1 < Inf) | any(gamma1 != 1) | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(any(df1<Inf)){
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
            }else{
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
            }else{
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            chunk1 <- sqrt(Ex1 - Exb1^2)
            chunk2 <- sqrt(Ex2 - Exb2^2)
            Exb1 <- 0 ## ATTENTION: this means that Qdet is at the mode. If we remove these two lines, Qdet is at the mean.
            Exb2 <- 0
        }else{
            chunk1 <- 1
            chunk2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        Q.sim <- pmax(Q.sim , 0)
        ## dQsim <- pmax(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])),0)
        ## dQpos1 <- pmax(rollmean(dQsim, k=3, fill=0, align="left"), 0)

        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + Q01*b1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + Q02*b2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*Q.sim[Q.sim<Q01]/Q01 + b1*Q01
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01/c1*(c1-1+(Q.sim[Q.sim>=Q01]/Q01)^c1)
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*Q.sim[Q.sim<Q02]/Q02 + b2*Q02
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02/c2*(c2-1+(Q.sim[Q.sim>=Q02]/Q02)^c2)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*(Q.sim[Q.sim<Q01]/Q01)^c1 + b1
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01*(1-c1^2+c1^2*(Q.sim[Q.sim>=Q01]/Q01)^(1/c1))
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*(Q.sim[Q.sim<Q02]/Q02)^c2 + b2
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02*(1-c2^2+c2^2*(Q.sim[Q.sim>=Q02]/Q02)^(1/c2))

        log.p <- rep(NA,n-2)
        sdt.b1 <- sd1[1:(n-1)]/chunk1
        sdt.b2 <- sd2[1:(n-1)]/chunk2
        sdt.f1 <- sd1[2:n]/chunk1
        sdt.f2 <- sd2[2:n]/chunk2
        ## gamma1 <- gamma1[2:n]
        ## dQpos2 <- abs(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])/(0.5*Q.sim[1:(n-1)] + 0.5*Q.sim[2:n] + psi)))
        ## dQpos2 <- rollmean(dQpos2, k=3, fill=0, align="left")
        P <- pmax(P-Plim,0)
        dt <- (t.obs[2:n]-t.obs[1:(n-1)])
        taus <- tau.min + (tau.max-tau.min)*exp(-tkP*P^l - tkQ*Q.sim^m)
        ## taus <- tau.max*tkP^l/(P^l+tkP^l)
        ## taus[((1:length(taus)) %% 20) == 0] <- 0
        ## taus <- rollmean(taus, k=3,fill=0)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-(pmax(P-ttP,0))/Q.sim*tkP-Q.sim/Q01*tkQ)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-1/alpha*sqrt(S/S0))
        gmms <- 1/taus
        gmms <- gmms[2:n]
        ## delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        ## delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        sdt.ini <- ifelse(auto[1], sd2[1]/chunk2, sd1[1]/chunk1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2[1], gamma1[1])
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1]+Exb1*sdt.ini),df=df.tmp,gamma=gamma.tmp,sigm=sdt.ini,log=TRUE)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]+Exb1*sdt.ini,df=df.tmp,gamma=gamma.tmp,sigm=sdt.ini,log.p=TRUE)
            }
        auto <- auto[2:n]
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- Q.obs[2:n]-Q.sim[2:n]+Exb1*sdt.f1
        dq.b1  <- Q.obs[1:(n-1)]-Q.sim[1:(n-1)]+Exb1*sdt.b1
        dq.f2  <- Q.obs[2:n]-Q.sim[2:n]+Exb2*sdt.f2
        dq.b2  <- Q.obs[1:(n-1)]-Q.sim[1:(n-1)]+Exb2*sdt.b2
        lower0 <- dq.b1<0
        higher0<- dq.b1>=0
        quant.b1 <- numeric(length(dq.b1))
        quant.b1[lower0]  <- qnorm(mypskt(dq.b1[lower0],df=df1,gamma=gamma1,sigm=sdt.b1[lower0],log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.b1[higher0] <- qnorm(mypskt(dq.b1[higher0],df=df1,gamma=gamma1,sigm=sdt.b1[higher0],low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        lower0 <- dq.b2<0
        higher0<- dq.b2>=0
        quant.b2 <- numeric(length(dq.b2))
        quant.b2[lower0] <- qnorm(mypskt(dq.b2[lower0],df=df2,gamma=gamma2,sigm=sdt.b2[lower0],log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.b2[higher0] <- qnorm(mypskt(dq.b2[higher0],df=df2,gamma=gamma2,sigm=sdt.b2[higher0],low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        lower0 <- dq.f1<0
        higher0<- dq.f1>=0
        quant.f1 <- numeric(length(dq.f1))
        quant.f1[lower0] <- qnorm(mypskt(dq.f1[lower0],df=df1,gamma=gamma1,sigm=sdt.f1[lower0],log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.f1[higher0] <- qnorm(mypskt(dq.f1[higher0],df=df1,gamma=gamma1,sigm=sdt.f1[higher0],low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        lower0 <- dq.f2<0
        higher0<- dq.f2>=0
        quant.f2 <- numeric(length(dq.f2))
        quant.f2[lower0] <- qnorm(mypskt(dq.f2[lower0],df=df2,gamma=gamma2,sigm=sdt.f2[lower0],log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.f2[higher0] <- qnorm(mypskt(dq.f2[higher0],df=df2,gamma=gamma2,sigm=sdt.f2[higher0],low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)

        densEta <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])), log=TRUE)
        densEtb <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,sigm=sdt.f2[ind.auto.qob],log=TRUE) + densEta - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-gmms[ind.auto.qob2]*dt[ind.auto.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto.qob2]*dt[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f2[ind.auto.qob3],df=df2,gamma=gamma2,sigm=sdt.f2[ind.auto.qob3],log=TRUE)
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        c <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,sigm=sdt.f1[ind.auto2.qob],log=TRUE)
        d <- dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto2.qob] <- c + densEtb - d
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-gmms[ind.auto2.qob2]*dt[ind.auto2.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto2.qob2]*dt[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- mydskt(dq.f1[ind.auto2.qob3],df=df1,gamma=gamma1,sigm=sdt.f1[ind.auto2.qob3],log=TRUE)
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(any(loglikeli==Inf,na.rm=T)){warning("infinite loglikeli encountered. Returning -Inf"); return(-Inf)}
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        white.noise1 <- (quant.f1 - quant.b1*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        white.noise2 <- (quant.f2 - quant.b2*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])))
        unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1,sigm=sdt.b1)
        unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2,sigm=sdt.b2)
        unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1,gamma=gamma1[length(gamma1)],sigm=sdt.f1[length(sdt.f1)])
        unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2,sigm=sdt.f2[length(sdt.f2)])
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), white.noise=c(ifelse(auto, white.noise2, white.noise1), ifelse(auto[length(auto)], white.noise2[length(white.noise2)], white.noise1[length(white.noise1)])), auto=c(FALSE,auto), a=densEta, b=densEtb, sd=sd1, taus=taus))
    }
    return(sum(loglikeli))
}


LogLikelihoodHydrology_la9f_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## reparametrized version of la9, where df=f(t)
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        sd01 <- Q01/5 ## in case Q01 is fitted
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ##tau1 <--1/log(par.likeli[paste(var.curr, "_ar1_lik", sep = "")])
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ##tau.max <- 20
        ##tau.min <- 0
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        P <- pmax(P,0)
        Q.sim <- pmax(Q.sim,0)
        df1 <- pmin(pmax(df1 - ifelse(Q.sim==0,Inf,tkQ*P^m),3),200)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        ## P <- pmax(P,0)
        ## gamma1 <- pmin(gamma1+tkQ*P^m,5)
        if ( any(df1 < Inf) | any(gamma1 != 1) | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(any(df1<Inf)){
                m1 <- df1/(df1-2)
                ## Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1) ## if we use this line, sd.t1 is not the "real" SD of the final distribution DQ, since we do not consider the influence of gamma on SD in the correct way.
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
                m1 <- 1
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                ## Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1) ## if we use this line, sd.t2 is not the "real" SD of the final distribution DQ, since we do not consider the influence of gamma on SD in the correct way.
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
                m2 <- 1
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        Q.sim <- pmax(Q.sim , 0)
        ## dQsim <- pmax(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])),0)
        ## dQpos1 <- pmax(rollmean(dQsim, k=3, fill=0, align="left"), 0)

        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + Q01*b1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + Q02*b2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*Q.sim[Q.sim<Q01]/Q01 + b1*Q01
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01/c1*(c1-1+(Q.sim[Q.sim>=Q01]/Q01)^c1)
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*Q.sim[Q.sim<Q02]/Q02 + b2*Q02
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02/c2*(c2-1+(Q.sim[Q.sim>=Q02]/Q02)^c2)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*(Q.sim[Q.sim<Q01]/Q01)^c1 + b1
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01*(1-c1^2+c1^2*(Q.sim[Q.sim>=Q01]/Q01)^(1/c1))
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*(Q.sim[Q.sim<Q02]/Q02)^c2 + b2
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02*(1-c2^2+c2^2*(Q.sim[Q.sim>=Q02]/Q02)^(1/c2))

        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1[1:(n-1)]
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1[2:n]
        scale.f2 <- sd2[2:n]/sd.t2
        ## gamma1 <- gamma1[2:n]
        ## dQpos2 <- abs(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])/(0.5*Q.sim[1:(n-1)] + 0.5*Q.sim[2:n] + psi)))
        ## dQpos2 <- rollmean(dQpos2, k=3, fill=0, align="left")
        P <- pmax(P,0)
        taus <- tau.min + (tau.max-tau.min)*exp(-tkP*P^l - tkQ*Q.sim^m)
        ## taus[((1:length(taus)) %% 20) == 0] <- 0
        ## taus <- rollmean(taus, k=3,fill=0)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-(pmax(P-ttP,0))/Q.sim*tkP-Q.sim/Q01*tkQ)
        ##taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-1/alpha*sqrt(S/S0))
        gmms <- 1/taus
        dt <- (t.obs[2:n]-t.obs[1:(n-1)])
        gmms <- gmms[2:n]
        ## delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        ## delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])*gmms[1:(n-1)]
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1[1])
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2[1], df1[1])
            gamma.tmp <- ifelse(auto[1], gamma2[1], gamma1[1])
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        df1 <- df1[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1##+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1##+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2##+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2##+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ## mns <- quant.b1[1]*exp(-delta1)
        ## sds <- sqrt(1-exp(-2*delta1))

        densEta <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])), log=TRUE)
        densEtb <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + densEta - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        ## ind2 <- auto & qob2 & time.recess[2:n]==0
        ## log.p[ind2] <- mydskt(dq.f2[ind2],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind2])
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-gmms[ind.auto.qob2]*dt[ind.auto.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto.qob2]*dt[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f2[ind.auto.qob3],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        c <- mydskt(dq.f1[ind.auto2.qob],df=df1[ind.auto2.qob],gamma=gamma1,log=TRUE)  - log(scale.f1[ind.auto2.qob])
        d <- dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto2.qob] <- c + densEtb - d
        ## ind1 <- !auto & qob2 & time.recess[2:n]==-Inf
        ## log.p[ind1] <- mydskt(dq.f1[ind1],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind1])
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-gmms[ind.auto2.qob2]*dt[ind.auto2.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto2.qob2]*dt[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- mydskt(dq.f1[ind.auto2.qob3],df=df1[ind.auto2.qob3],gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob3])
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(any(loglikeli==Inf,na.rm=T)){warning("infinite loglikeli encountered. Returning -Inf"); return(-Inf)}
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        white.noise1 <- (quant.f1 - quant.b1*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        white.noise2 <- (quant.f2 - quant.b2*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])))
        unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1)
        unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2)
        unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1[length(df1)],gamma=gamma1[length(gamma1)])
        unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), white.noise=c(ifelse(auto, white.noise2, white.noise1), ifelse(auto[length(auto)], white.noise2[length(white.noise2)], white.noise1[length(white.noise1)])), auto=c(FALSE,auto), a=densEta, b=densEtb, sd=sd1, taus=taus))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la9_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## la9 assumes that the quantiles follow an integrated OU-process, where the mean of the OU-process is depending on the quantile. The quantile has another stochastic component depending on the rainfall.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ##tau1 <--1/log(par.likeli[paste(var.curr, "_ar1_lik", sep = "")])
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        ## gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau1
        delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau2
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)

        a <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])), log=TRUE)
        b <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + a - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        ## ind2 <- auto & qob2 & time.recess[2:n]==0
        ## log.p[ind2] <- mydskt(dq.f2[ind2],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind2])
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-delta2[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta2[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f2[ind.auto.qob3],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        log.p[ind.auto2.qob] <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob]) + b - dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        ## ind1 <- !auto & qob2 & time.recess[2:n]==-Inf
        ## log.p[ind1] <- mydskt(dq.f1[ind1],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind1])
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-delta1[ind.auto2.qob2]),sd=sqrt(1-exp(-2*delta1[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- mydskt(dq.f1[ind.auto2.qob3],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob3])
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])))
        unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1)
        unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2)
        unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1,gamma=gamma1)
        unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), auto=c(FALSE,auto), a=a, b=b))
    }
    return(sum(loglikeli))
}


LogLikelihoodHydrology_la92_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## la92 assumes that tau has two different values in two different arbitrary periods delineated by <auto>.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        ## Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        ## sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        ## a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        ## b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ## tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        tau2 <- tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        ## mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( any(df1 < Inf) | any(gamma1 != 1) | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(any(df1<Inf)){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
                m1 <- 1
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
                m2 <- 1
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
            ## Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1) ## for further use, this is the mean of the unscaled (original) skewed t-distribution
            ## Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            Exb1 <- 0 ## ATTENTION: this means that Qdet is at the mode, not the mean...
            Exb2 <- 0
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        ## mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        ## Q.sim <- mu.tmp*Q.sim
        Q.sim <- pmax(Q.sim, 0)
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + Q01*b1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + Q02*b2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ## sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 + d*pmax(P-Plim,0)^e
        ## sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2 + d*pmax(P-Plim,0)^e
        ## sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ## sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau1
        delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau2
        ##auto <- ifelse(P>Plim,FALSE,TRUE) ## ATTENTION: auto of the input is overwritten here.
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                ## log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                ## log.p.ini <- mypskt(-Q.sim[1]/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)

        a <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])), log=TRUE)
        b <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + a - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-delta2[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta2[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f2[ind.auto.qob3],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)
        ## ind2 <- auto & qob2 & time.recess[2:n]==0
        ## log.p[ind2] <- mydskt(dq.f2[ind2],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind2])

        log.p[ind.auto2.qob] <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob]) + b - dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-delta1[ind.auto2.qob2]),sd=sqrt(1-exp(-2*delta1[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- mydskt(dq.f1[ind.auto2.qob3],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob3])
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)
        ## ind1 <- P[2:n] > 1.5*Q.sim[2:n]##!auto & qob2 & time.recess[2:n]==-Inf
        ## log.p[ind1] <- mydskt(dq.f1[ind1],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind1])

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        ## innovation0 <- pnorm(quant.f1[ind1], mean=0, sd=1)
        white.noise1 <- (quant.f1 - quant.b1*exp(-delta1))/sqrt(1-exp(-2*delta1))
        white.noise2 <- (quant.f2 - quant.b2*exp(-delta2))/sqrt(1-exp(-2*delta2))
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])))
        innovation=ifelse(auto,innovation2,innovation1)
        ## innovation[ind1] <- innovation0
        innovation <- c(NA,innovation)
        unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1)
        unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2)
        unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1,gamma=gamma1)
        unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=innovation, auto=c(FALSE,auto), a=a, b=b, sd=sd1, white.noise=c(ifelse(auto, white.noise2, white.noise1), ifelse(auto[length(auto)], white.noise2[length(white.noise2)], white.noise1[length(white.noise1)])), taus=tau1))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la93_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## la93 extends the idea of la92 by assuming a time-dependent tau!
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        k      <- par.likeli[paste(var.curr, "_k_lik", sep = "")]
        ## gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau.max) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ## approximate or solve some integrals for the auxiliary reservoir and the time-dependent tau:
        S <- numeric(length(P))
        S0 <- mean(P)/(3*k)
        int <- numeric(length(P))
        int2 <- numeric(length(P))
        dt <- t.obs[2:n] - t.obs[1:(n-1)] ## ATTENTION: we assume that the input times are equal to the layout times!!!
        int[1] <- 0
        int2[1] <- 0
        S[1] <- 3*S0
        f.S <- function(dt, St0, P, k) St0*exp(-k*dt)+P*dt
        f.tau <- function(t, t0, tau.max, P, k, S0, St0){
            S.tmp <- St0*exp(-k*(t-t0))+P*(t-t0)
            res <- pmin(1/tau.max*exp((P*(t-t0)/(k*S.tmp)+(S.tmp/S0))), 10)
            return(res)
        }
        f.tau2 <- function(t, t0, tau.max, P, k, S0, St0, I0, subd=100){
            S.tmp <- St0*exp(-k*(t-t0))+P*(t-t0) ## ATTENTION the precipitation is distributed evenly over [t0,t]
            tt <- pmin(1/tau.max*exp((P*(t-t0)/(k*S.tmp)+(S.tmp/S0))), 10)
            res <- tt*exp(2*(I0+integrate(f.tau, lower=0, upper=t, t0=t0, tau.max=tau.max, P=P, k=k, S0=S0, St0=St0, subdivisions=subd, rel.tol=rel.tol)$value))
            return(res)
        }
        subdivisions <- 100
        rel.tol <- 0.01
        for(i in 2:length(S)){
            S[i] <- integrate(f.S, lower=0, upper=dt[i-1], St0=S[i-1], k=k, P=P[i], subdivisions=subdivisions, rel.tol=rel.tol)$value ## ATTENTION: this causes problems if P[i] != 0 and dt[i-1] is larger than the time step for P
            int[i] <- integrate(f.tau, lower=0, upper=dt[i-1], t0=0, tau.max=tau.max, P=P[i], k=k, S0=S0, St0=S[i-1], subdivisions=subdivisions, rel.tol=rel.tol)$value ## ATTENTION: see comment 1 line above
            int2[i] <- integrate(f.tau2, lower=0, upper=dt[i-1], t0=0, tau.max=tau.max, P=P[i], k=k, S0=S0, St0=S[i-1], I0=int[i-1], subd=subdivisions, subdivisions=subdivisions, rel.tol=rel.tol)$value
        }
        int <- int[2:n]
        int2 <- int2[2:n]
        int.cumu <- cumsum(int)
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1+Exb1
        dq.b  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1+Exb1
        quant.b <- qnorm(mypskt(dq.b,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        quant.f <- qnorm(mypskt(dq.f,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        mns <- quant.b*exp(-int)
        print(summary(mns))
        sds <- sqrt(2*exp(-2*int.cumu)*int2)
        print(summary(sds))
        a <- dnorm(quant.f[ind.auto.qob], mean=mns[ind.auto.qob], sd=sds[ind.auto.qob], log=TRUE)
        log.p[ind.auto.qob] <- mydskt(dq.f[ind.auto.qob],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto.qob]) + a - dnorm(quant.f[ind.auto.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(quant.f[ind.auto.qob2], mean=mns[ind.auto.qob2],sd=sds[ind.auto.qob2], log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f[ind.auto.qob3],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        innovation <- pnorm(quant.f[ind.auto.qob], mean=quant.b[ind.auto.qob]*exp(-int[ind.auto.qob]), sd=sqrt(2*exp(-2*int.cumu[ind.auto.qob])*int2[ind.auto.qob]))
        innovation <- c(NA,innovation)
        unifvalue <- mypskt(dq.b,df=df1,gamma=gamma1)
        unifvalue.f <- mypskt(dq.f[length(dq.f)],df=df1,gamma=gamma1)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(unifvalue, unifvalue.f[length(unifvalue.f)]), quant= c(quant.b, quant.f[length(quant.f)]), innovation=innovation, auto=c(FALSE,auto), a=a))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la94_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## la93 extends the idea of la92 by assuming a time-dependent tau!
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        k      <- par.likeli[paste(var.curr, "_k_lik", sep = "")]
        ## gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau.max) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ## approximate or solve some integrals for the auxiliary reservoir and the time-dependent tau:
        S <- numeric(length(P))
        S0 <- mean(P)/(3*k)
        int <- numeric(length(P))
        int2 <- numeric(length(P))
        dt <- t.obs[2:n] - t.obs[1:(n-1)] ## ATTENTION: we assume that the input times are equal to the layout times!!!
        int[1] <- 0
        int2[1] <- 0
        S[1] <- 3*S0
        f.S <- function(dt, St0, P, k) St0*exp(-k*dt)+P*dt
        f.tau <- function(t, t0, tau.max, P, k, S0, St0){
            S.tmp <- St0*exp(-k*(t-t0))+P*(t-t0)
            res <- pmin(1/tau.max*exp((P*(t-t0)/(k*S.tmp)+(S.tmp/S0))), 10)
            return(res)
        }
        f.tau2 <- function(t, t0, tau.max, P, k, S0, St0, I0, subd=100){
            S.tmp <- St0*exp(-k*(t-t0))+P*(t-t0) ## ATTENTION the precipitation is distributed evenly over [t0,t]
            tt <- pmin(1/tau.max*exp((P*(t-t0)/(k*S.tmp)+(S.tmp/S0))), 10)
            res <- tt*exp(2*(I0+integrate(f.tau, lower=0, upper=t, t0=t0, tau.max=tau.max, P=P, k=k, S0=S0, St0=St0, subdivisions=subd, rel.tol=rel.tol)$value))
            return(res)
        }
        subdivisions <- 100
        rel.tol <- 0.01
        for(i in 2:length(S)){
            S[i] <- integrate(f.S, lower=0, upper=dt[i-1], St0=S[i-1], k=k, P=P[i], subdivisions=subdivisions, rel.tol=rel.tol)$value ## ATTENTION: this causes problems if P[i] != 0 and dt[i-1] is larger than the time step for P
            int[i] <- integrate(f.tau, lower=0, upper=dt[i-1], t0=0, tau.max=tau.max, P=P[i], k=k, S0=S0, St0=S[i-1], subdivisions=subdivisions, rel.tol=rel.tol)$value ## ATTENTION: see comment 1 line above
            int2[i] <- integrate(f.tau2, lower=0, upper=dt[i-1], t0=0, tau.max=tau.max, P=P[i], k=k, S0=S0, St0=S[i-1], I0=int[i-1], subd=subdivisions, subdivisions=subdivisions, rel.tol=rel.tol)$value
        }
        int <- int[2:n]
        int2 <- int2[2:n]
        int.cumu <- cumsum(int)
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1+Exb1
        dq.b  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1+Exb1
        quant.b <- qnorm(mypskt(dq.b,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        quant.f <- qnorm(mypskt(dq.f,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        mns <- quant.b*exp(-int)
        print(summary(mns))
        sds <- sqrt(2*exp(-2*int.cumu)*int2)
        print(summary(sds))
        a <- dnorm(quant.f[ind.auto.qob], mean=mns[ind.auto.qob], sd=sds[ind.auto.qob], log=TRUE)
        log.p[ind.auto.qob] <- mydskt(dq.f[ind.auto.qob],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto.qob]) + a - dnorm(quant.f[ind.auto.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(quant.f[ind.auto.qob2], mean=mns[ind.auto.qob2],sd=sds[ind.auto.qob2], log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f[ind.auto.qob3],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        innovation <- pnorm(quant.f[ind.auto.qob], mean=quant.b[ind.auto.qob]*exp(-int[ind.auto.qob]), sd=sqrt(2*exp(-2*int.cumu[ind.auto.qob])*int2[ind.auto.qob]))
        innovation <- c(NA,innovation)
        unifvalue <- mypskt(dq.b,df=df1,gamma=gamma1)
        unifvalue.f <- mypskt(dq.f[length(dq.f)],df=df1,gamma=gamma1)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(unifvalue, unifvalue.f[length(unifvalue.f)]), quant= c(quant.b, quant.f[length(quant.f)]), innovation=innovation, auto=c(FALSE,auto), a=a))
    }
    return(sum(loglikeli))
}


LogLikelihoodHydrology_la9_compound <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## la9 assumes that the quantiles follow an integrated OU-process, where the mean of the OU-process is depending on the quantile. The quantile has another stochastic component depending on the rainfall.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        ## Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        ## sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        ## a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        ## b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ## tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        tau2 <- tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        pout1 <- par.likeli[paste(var.curr, "_p_lik", sep = "")]
        ## pout2 <- par.likeli[paste(var.curr, "_p2_lik", sep = "")]
        pout2 <- pout1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if(df1<Inf){
            m1 <- df1/(df1-2)
            ## Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1) ## if we use this line, sd.t1 is not the "real" SD of the final distribution DQ, since we do not consider the influence of gamma on SD in the correct way.
            Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
            m1 <- 1
        }else{
            m1 <- 1
            Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
        }
        if(df2<Inf){
            m2 <- df2/(df2-2)
            ## Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1) ## if we use this line, sd.t2 is not the "real" SD of the final distribution DQ, since we do not consider the influence of gamma on SD in the correct way.
            Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
            m2 <- 1
        }else{
            m2 <- 1
            Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
        }
        Ez1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
        Ez2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
        dfout <- 3
        Ek1 <- dfout/(dfout-2) ## dfout are the degrees of freedom of the outer t-distribution used to describe the outliers
        Ek2 <- dfout/(dfout-2)
        Ex1 <- (1-pout1)*Ez1 + pout1*Ek1
        Ex2 <- (1-pout2)*Ez2 + pout2*Ek2
        sd.t1 <- sqrt(Ex1 - ((1-pout1)*Exb1)^2)
        sd.t2 <- sqrt(Ex2 - ((1-pout2)*Exb2)^2)
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        Q.sim <- pmax(Q.sim, 0)
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + Q01*b1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + Q02*b2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau1
        delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau2
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        pout.ini <- ifelse(auto[1], pout2, pout1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- dcompound((Q.obs[1]-Q.sim[1])/scale.ini,df=df.tmp,df2=dfout,gamma=gamma.tmp,p=pout.ini,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- pcompound(-Q.sim[1]/scale.ini,df=df.tmp,df2=dfout,gamma=gamma.tmp,p=pout.ini,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1##+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1##+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2##+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2##+Exb2
        quant.b1 <- qnorm(pcompound(dq.b1,df=df1,df2=dfout,gamma=gamma1,p=pout1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(pcompound(dq.b2,df=df2,df2=dfout,gamma=gamma2,p=pout2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        quant.f1 <- qnorm(pcompound(dq.f1,df=df1,df2=dfout,gamma=gamma1,p=pout1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(pcompound(dq.f2,df=df2,df2=dfout,gamma=gamma2,p=pout2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)

        a <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])), log=TRUE)
        b <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- dcompound(dq.f2[ind.auto.qob],df=df2,df2=dfout,gamma=gamma2,p=pout2,log=TRUE) - log(scale.f2[ind.auto.qob]) + a - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-delta2[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta2[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- dcompound(dq.f2[ind.auto.qob3],df=df2,df2=dfout,gamma=gamma2,p=pout2,log=TRUE) - log(scale.f2[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        log.p[ind.auto2.qob] <- dcompound(dq.f1[ind.auto2.qob],df=df1,df2=dfout,gamma=gamma1,p=pout1,log=TRUE) - log(scale.f1[ind.auto2.qob]) + b - dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-delta1[ind.auto2.qob2]),sd=sqrt(1-exp(-2*delta1[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- dcompound(dq.f1[ind.auto2.qob3],df=df1,df2=dfout,gamma=gamma1,p=pout1,log=TRUE) - log(scale.f1[ind.auto2.qob3])
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        white.noise1 <- (quant.f1 - quant.b1*exp(-delta1))/sqrt(1-exp(-2*delta1))
        white.noise2 <- (quant.f2 - quant.b2*exp(-delta2))/sqrt(1-exp(-2*delta2))
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])))
        innovation=ifelse(auto,innovation2,innovation1)
        innovation <- c(NA,innovation)
        unifvalue1 <- pcompound(dq.b1,df=df1,df2=dfout,gamma=gamma1,p=pout1)
        unifvalue2 <- pcompound(dq.b2,df=df2,df2=dfout,gamma=gamma2,p=pout2)
        unifvalue.f1 <- pcompound(dq.f1[length(dq.f1)],df=df1,df2=dfout,gamma=gamma1,p=pout1)
        unifvalue.f2 <- pcompound(dq.f2[length(dq.f2)],df=df2,df2=dfout,gamma=gamma2,p=pout2)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=innovation, auto=c(FALSE,auto), a=a, b=b, sd=sd1, white.noise=c(ifelse(auto, white.noise2, white.noise1), ifelse(auto[length(auto)], white.noise2[length(white.noise2)], white.noise1[length(white.noise1)]))))
    }
    return(sum(loglikeli))
}


LogLikelihoodHydrology_la9_mode_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## la9 assumes that the quantiles follow an integrated OU-process, where the mean of the OU-process is depending on the quantile. The quantile has another stochastic component depending on the rainfall.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ##c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ##tau1 <--1/log(par.likeli[paste(var.curr, "_ar1_lik", sep = "")])
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ##df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            ##Exb1 <- 0
            ##Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau1
        delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau2
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1##+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1##+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2##+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2##+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        ind.auto.qob3<- auto & !qob & qob2
        ind.auto2.qob3<- !auto & !qob & qob2
        ind.auto.qob4<- auto & !qob & !qob2
        ind.auto2.qob4<- !auto & !qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)

        a <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])), log=TRUE)
        b <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + a - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        ## ind2 <- auto & qob2 & time.recess[3:n]==0
        ## log.p[ind2] <- mydskt(dq.f2[ind2],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind2])
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-delta2[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta2[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- mydskt(dq.f2[ind.auto.qob3],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob3])
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        log.p[ind.auto2.qob] <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob]) + b - dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        ## ind1 <- !auto & qob2 & time.recess[3:n]==-Inf
        ## log.p[ind1] <- mydskt(dq.f1[ind1],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind1])
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-delta1[ind.auto2.qob2]),sd=sqrt(1-exp(-2*delta1[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- mydskt(dq.f1[ind.auto2.qob3],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob3])
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    if(apprx){
        loglikeli[is.nan(loglikeli)] <- -Inf
    }else{
        if(any(is.nan(loglikeli))){
            nprob <- sum(is.nan(loglikeli))
            if(nprob<50)cat("problem quantiles at: ", which(is.nan(loglikeli)), "\n")
            stop(paste("numerical problems: ", nprob, " infinite quantile(s)", sep=""))
        }
    }
    if(!verbose){
        innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])))
        innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])))
        unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1)
        unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2)
        unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1,gamma=gamma1)
        unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2)
        return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), auto=c(FALSE,auto), a=a, b=b))
    }
    return(sum(loglikeli))
}


LogLikelihoodHydrology_la10_fast_skewt <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){ ## la10 assumes that the quantiles during the events follow an OU-process and the quantiles in the recessions and low-flows follow an OU-process with correlated (OU) noise.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    ##y.obs[na.y.obs] <- 0 ## the skewt distribution cannot handle NAs, so we have to replace them by an arbitrary value. The loglikeli should be set 0 for the cases of na.y.obs later...
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        ##Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        ##sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        ##a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        ##b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ##c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau1 <--1/log(par.likeli[paste(var.curr, "_ar1_lik", sep = "")])
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ##df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ##gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            m1 <- df1/(df1-2)
            m2 <- df2/(df2-2)
            Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau1) | is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        log.p <- rep(NA,n-2)
        scale.b1 <- sd1[1:(n-1)]/sd.t1
        scale.b2 <- sd2[1:(n-1)]/sd.t2
        scale.f1 <- sd1[2:n]/sd.t1
        scale.f2 <- sd2[2:n]/sd.t2
        delta1 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau1
        delta2 <- (t.obs[2:n] - t.obs[1:(n-1)])/tau2
        scale.ini <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1])/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log=TRUE) - log(scale.ini)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]/scale.ini+Exb2,df=df.tmp,gamma=gamma.tmp,log.p=TRUE)
            }
        auto <- auto[2:n]
        ## cat("d: ", d, "\n")
        ## x <- 1:(n-1)
        ## sec <- 2000:5000
        ## plot(x[sec], Q.sim[2:n][sec], type="l")
        ## points(x[sec][!auto[sec]], Q.sim[2:n][sec][!auto[sec]], pch=19, col="red")
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f1+Exb1
        dq.b1  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b1+Exb1
        dq.f2  <- (Q.obs[2:n]-Q.sim[2:n])/scale.f2+Exb2
        dq.b2  <- (Q.obs[1:(n-1)]-Q.sim[1:(n-1)])/scale.b2+Exb2
        quant.b1 <- qnorm(mypskt(dq.b1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.b2 <- qnorm(mypskt(dq.b2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.auto.qob <- auto & qob & qob2
        ind.auto2.qob <- !auto &qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        ind.auto2.qob2<- !auto & qob & !qob2
        quant.f1 <- qnorm(mypskt(dq.f1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        quant.f2 <- qnorm(mypskt(dq.f2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)

        a <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-delta2[ind.auto.qob]), sd=sqrt(1-exp(-2*delta2[ind.auto.qob])), log=TRUE)
        b <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-delta1[ind.auto2.qob]), sd=sqrt(1-exp(-2*delta1[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- mydskt(dq.f2[ind.auto.qob],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind.auto.qob]) + a - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        ## ind1 <- auto & qob2 & time.recess[3:n]==0
        ## log.p[ind1] <- mydskt(dq.f2[ind1],df=df2,gamma=gamma2,log=TRUE) - log(scale.f2[ind1])
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-delta2[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta2[ind.auto.qob2])), log=TRUE)

        log.p[ind.auto2.qob] <- mydskt(dq.f1[ind.auto2.qob],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind.auto2.qob]) + b - dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        ## ind2 <- !auto & qob2 & time.recess[3:n]==-Inf
        ## log.p[ind2] <- mydskt(dq.f1[ind2],df=df1,gamma=gamma1,log=TRUE) - log(scale.f1[ind2])
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-delta1[ind.auto2.qob2]),sd=sqrt(1-exp(-2*delta1[ind.auto2.qob2])), log=TRUE)
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- -723
    }
    if(!verbose){
        return(list(smm=sum(loglikeli), loglik=loglikeli, quant=c(ifelse(auto, quant.b2, quant.b1), quant.f2[length(quant.f2)]), auto=c(FALSE,auto), a=a, b=b))
    }
    return(sum(loglikeli))
}



LogLikelihoodHydrology_mod_n<- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n)
        if ( Q.obs[1] > 0 ){
            log.p[1] <- dnorm(Q.obs[1],mean=Q.sim[1],sd=sd[1],log=TRUE)
        }else{
            log.p[1] <- pnorm(0,mean=Q.sim[1],sd=sd[1],log=TRUE)
        }
        for ( i in 2:n ){
            if ( Q.sim[i]<Q.sim[i-1] & Q.obs[i-1]>0 ){
                delta <- (t.obs[i]-t.obs[i-1])/tau
                eta.1 <- qnorm(pnorm(Q.obs[i-1],mean=Q.sim[i-1],sd=sd[i-1]),mean=0,sd=1)
                if ( Q.obs[i] > 0 ){
                    eta.2 <- qnorm(pnorm(Q.obs[i],mean=Q.sim[i],sd=sd[i]),mean=0,sd=1)
                    log.p[i] <- dnorm(Q.obs[i],mean=Q.sim[i],sd=sd[i],log=TRUE) +
                        dnorm(eta.2,mean=eta.1*exp(-delta),sd=sqrt(1-exp(-2*delta)),log=TRUE) -
                        dnorm(eta.2,mean=0,sd=1,log=TRUE)
                }else{
                    eta.2 <- qnorm(pnorm(-1*Q.sim[i],mean=0,sd=sd[i]), mean=0, sd=1)
                    log.p[i] <- pnorm(eta.2,mean=eta.1*exp(-delta),sd=sqrt(1-exp(-2*delta)),log=TRUE)
                }
            }else{
                if ( Q.obs[i] > 0 ){
                    log.p[i] <- dnorm(Q.obs[i],mean=Q.sim[i],sd=sd[i],log=TRUE)
                }else{
                    log.p[i] <- pnorm(0,mean=Q.sim[i],sd=sd[i],log=TRUE)
                }
            }
        }
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    ss <- sum(loglikeli)
    return(ss)
}

LogLikelihoodHydrology_mod_n_fast <- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n-1)
        sd.b <- sd[1:(n-1)]
        sd.f <- sd[2:n]
        sd.ini <- sd[1]
        if ( Q.obs[1] > 0 ){
            log.p.ini <- dnorm(Q.obs[1]-Q.sim[1],sd=sd.ini,log=TRUE)
        }else{
            log.p.ini <- pnorm(0,mean=Q.sim[1],sd=sd.ini,log=TRUE)
        }
        auto  <- Q.sim[2:n] < Q.sim[1:(n-1)]
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        delta <- (t.obs[2:n] - t.obs[1:(n-1)])/tau
        dq.f  <- Q.obs[2:n]-Q.sim[2:n]
        dq.b  <- Q.obs[1:(n-1)]-Q.sim[1:(n-1)]
        eta.1 <- qnorm(pnorm(dq.b,sd=sd.b),mean=0,sd=1)
        ind.auto.qob <- auto & qob & qob2
        ind.auto.qob2<- auto & qob & !qob2
        eta.2.aq <- qnorm(pnorm(dq.f[ind.auto.qob],sd=sd.f[ind.auto.qob]),mean=0,sd=1)
        eta.2.aq2<- qnorm(pnorm(-1*Q.sim[2:n][ind.auto.qob2],sd=sd.f[ind.auto.qob2]),mean=0,sd=1)
        log.p[ind.auto.qob] <- dnorm(dq.f[ind.auto.qob],sd=sd.f[ind.auto.qob],log=TRUE) + dnorm(eta.2.aq,mean=eta.1[ind.auto.qob]*exp(-delta[ind.auto.qob]),sd=sqrt(1-exp(-2*delta[ind.auto.qob])),log=TRUE) - dnorm(eta.2.aq,mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(eta.2.aq2,mean=eta.1[ind.auto.qob2]*exp(-delta[ind.auto.qob2]),sd=sqrt(1-exp(-2*delta[ind.auto.qob2])),log=TRUE)
        log.p[!(auto & qob)] <- dnorm(dq.f[!(auto&qob)],sd=sd.f[!(auto&qob)],log=TRUE)
        log.p[!(auto & qob) & !qob2] <- pnorm(-1*Q.sim[2:n][!(auto & qob) & !qob2],sd=sd.f[!(auto & qob) & !qob2],log=TRUE)
        log.p <- c(log.p.ini, log.p)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    ss <- sum(loglikeli)
    return(ss)
}


LogLikelihoodHydrology_noautocorr <- function(par.model, run.model, layout, y.obs, par.likeli, apprx, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if (is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n)
        scale <- sd[1]/sd.t
        if ( Q.obs[1] > 0 ){
            log.p[1] <- dt((Q.obs[1]-Q.sim[1])/scale,df=df,log=TRUE) - log(scale)}
        else{
            log.p[1] <- pt(Q.sim[1]/scale,df=df,log=TRUE)}
        for ( i in 2:n ){
            scale <- sd[i]/sd.t
            if ( Q.obs[i] > 0 ){
                log.p[i] <- dt((Q.obs[i]-Q.sim[i])/scale,df=df,log=TRUE) - log(scale)
            }else{
                log.p[i] <- pt(-1*Q.sim[i]/scale,df=df,log=TRUE)
            }
        }
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_wls <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if (is.na(a))
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        log.p <- rep(NA,n)
        log.p <- dnorm((Q.sim-Q.obs)/sd, mean=0, sd=1, log=TRUE) - log(sd)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){
        return(list(smm=sum(loglikeli), loglik=loglikeli, auto=c(FALSE,FALSE,auto)))
    }
    return(sum(loglikeli))
}

LogLikelihoodHydrology_AR1 <- function(par.model, run.model, layout, y.obs, P, par.likeli, mu, rep.mu.times, time.recess, auto, apprx, verbose, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    ##na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ar1<- par.likeli[paste(var.curr, "_ar1_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if (is.na(a))
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*(Q.sim/Q0)^c+b
        log.p <- rep(NA,n)
        res <- (Q.obs-Q.sim)/sd
        log.p[2:n] <- dnorm(res[2:n], mean=res[1:(n-1)]*ar1, sd=sqrt(1-ar1^2), log=TRUE) - log(sd[2:n])
        log.p[1] <- dnorm(res[1], mean=0, sd=1, log=TRUE) - log(sd[1])
        loglikeli[ind.var] <- log.p
    }
    ##loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    if(apprx){
        loglikeli[is.infinite(loglikeli)] <- min(loglikeli[!is.infinite(loglikeli)])
    }
    if(!verbose){
        return(list(smm=sum(loglikeli), loglik=loglikeli, quant=res))
    }
    return(sum(loglikeli))
}


LogLikelihoodHydrology_iidnorm <- function(par.model, run.model, layout, y.obs, par.likeli, ...){
    L = layout$layout
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    na.y.obs <- is.na(y.obs)
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    for(var.curr in vars){
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if (is.na(a))
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a
        log.p <- rep(NA,n)
        log.p <- dnorm(Q.sim-Q.obs, mean=0, sd=sd, log=TRUE)
        loglikeli[ind.var] <- log.p
    }
    loglikeli[na.y.obs] <- 0
    loglikeli[is.nan(loglikeli)] <- -Inf
    return(sum(loglikeli))
}

## =======================================================================================================
## Likelihood Samplers
## =======================================================================================================
## These are functions that sample from the likelihood (e.g. for prediction)

sampling_wrapper <- function(sudriv, brn.in=0, sample.par=TRUE, n.sample=1, biased=FALSE, sample.likeli=TRUE, rep.mu.times=NA, time.recess=NA, auto=NA, mu=NA, eta=NA){
    ## sample from a population (sample) of parameters
    if(all(is.na(time.recess))){
        time.recess <- rep(-1, nrow(su$layout$layout))
    }
    if(all(is.na(auto))){
        auto <- rep(FALSE, nrow(su$layout$layout))
    }
    if(sample.par){
        if(is.null(sudriv$parameter.sample)){
            warning("Sudriv object does not contain parameter sample to draw from. Drawing from the prior ...")
            ## Develop: draw parameter samples from prior distribution...
        }else{ ## Draw parameter sample from existing sample (representing e.g. posterior)
            ndim <- length(dim(sudriv$parameter.sample))
            if(brn.in >= nrow(sudriv$parameter.sample)) stop("brn.in is longer than chain ...")
            if(ndim==3) s <- sudriv$parameter.sample[(brn.in+1):nrow(sudriv$parameter.sample),,]
            if(ndim==2) s <- sudriv$parameter.sample[(brn.in+1):nrow(sudriv$parameter.sample),]
            ind.chosen <- sample(x=1:nrow(s), size=n.sample, replace=TRUE)
            if(ndim==3) wlk.chosen <- sample(x=1:(dim(s)[3]), size=n.sample, replace=TRUE)
        }
    }
    likeli.sample <- matrix(nrow=n.sample, ncol=length(sudriv$layout$calib)+length(sudriv$layout$pred))
    rtruncnorm <- function (n, a = -Inf, b = Inf, mean = 0, sd = 1){
    if (length(n) > 1)
        n <- length(n)
    if (length(n) > 1)
        n <- length(n)
    else if (!is.numeric(n))
        stop("non-numeric argument n.")
    .Call("do_rtruncnorm", as.integer(n), a, b, mean, sd, PACKAGE="truncnorm")}

    for(i in 1:n.sample){
        if(sample.par){
            if(ndim==3) x0 <- s[ind.chosen[i],,wlk.chosen[i]]
            if(ndim==2) x0 <- s[ind.chosen[i],]
            flp <- sudriv$likelihood$par.fit
            fmp <- sudriv$model$par.fit
            l.fit.lik <- sum(flp)
            l.fit <- sum(c(fmp,flp))
            if(biased){
                x0 <- x0[1:l.fit]
            }
            ## =======================================================
            ## update the likelihood parameters with the ones from x0
            par.lik.fit <- x0[(length(x0)-sum(flp)+1):length(x0)]
            ## make sure they are within the bounds
            if(!sudriv$settings$OPT.bounded){
                lower_bound <- sudriv$likelihood$lb[as.logical(flp)]
                upper_bound <- sudriv$likelihood$ub[as.logical(flp)]
                par.lik.fit <- constrain_parameters(par.lik.fit, lower_bound, upper_bound)
            }
            ## update likelihood parameters
            sudriv$likelihood$parameters[which(flp != 0)] <- par.lik.fit

            ## =======================================================
            ## update model parameters with the ones from x0
            par.mod.fit <- x0[1:(length(x0)-l.fit.lik)]
            ## make sure they are within the bounds
            if(!sudriv$settings$OPT.bounded){
                lower_bound <- sudriv$model$args$parLo[as.logical(fmp)]
                upper_bound <- sudriv$model$args$parHi[as.logical(fmp)]
                par.mod.fit <- constrain_parameters(par.mod.fit, lower_bound, upper_bound)
            }
            ## update model parameters
            sudriv$model$parameters[which(fmp != 0)] <- par.mod.fit
        }
        if(sample.likeli){
            ## =======================================================
            ## prepare arguments for the likelihood sampler
            if(biased){
                if(is.na(mu[1])){
                    mu <- rtruncnorm(length(rep.mu.times), a=0.001, b=Inf, mean=1, sd=0.34)
                }
            }else{
                mu <- rep(1, length(rep.mu.times))
            }
            likeli.args           <- list()
            likeli.args$par.model <- sudriv$model$parameters
            likeli.args$run.model <- run.model
            likeli.args$layout    <- sudriv$layout
            likeli.args$P         <- sudriv$input$P.roll##[sudriv$layout$pred]
            likeli.args$par.likeli<- ifelse(as.logical(sudriv$likelihood$tran), exp(sudriv$likelihood$parameters), sudriv$likelihood$parameters)
            likeli.args$mu <- mu
            likeli.args$time.recess <- time.recess
            likeli.args$rep.mu.times <- rep.mu.times
            likeli.args$auto <- auto
            names(likeli.args$par.likeli) <- names(sudriv$likelihood$parameters)
            likeli.args$sudriv    <- sudriv
            f.sample <- sudriv$likelihood$f.sample
            ## =======================================================
            ## sample from the likelihood
            likeli.sample[i,] <- do.call(f.sample, likeli.args)
        }else{## in this case, we just run the deterministic model (propagate parameter uncertainty only)
            par <- sudriv$model$parameters
            ## if(is.null(sudriv$layout$pred)){
            L <- sudriv$layout
            L$layout <- L$layout[c(sudriv$layout$calib,sudriv$layout$pred),]
            ## }else{
            ##     L <- sudriv$layout
            ##     L$layout    <- L$layout[L$pred,]
            ## }
            likeli.sample[i,] <- as.numeric(run.model(par=par, layout=L, sudriv=sudriv)$incld.lmpd)
        }
    }
    return(likeli.sample)
}

LogLikelihoodHydrology_mod_sample <- function(par.model, run.model, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        samp.var <- Q.sim
        scale <- sd[1]/sd.t
        samp.var[1] <- max(samp.var[1] + rt(n=1, df=df) * scale, 0)
        for ( i in 2:n ){
            if ( Q.sim[i]<Q.sim[i-1] & samp.var[i-1]>0 ){
                delta <- (t.mod[i]-t.mod[i-1])/tau
                eta.1 <- qnorm(pt((samp.var[i-1]-Q.sim[i-1])/scale,df=df),mean=0,sd=1)
                eta.2 <- rnorm(n=1, mean=eta.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + qt(pnorm(eta.2), df=df) * scale, 0)
            }else{
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + rt(n=1, df=df) * scale, 0)
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la_sample <- function(par.model, run.model, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        samp.var <- Q.sim
        scale <- sd[1]/sd.t
        samp.var[1] <- max(samp.var[1] + rt(n=1, df=df) * scale, 0)
        for ( i in 2:n ){
            if (samp.var[i-1]>0 ){
                delta <- (t.mod[i]-t.mod[i-1])/tau
                eta.1 <- qnorm(pt((samp.var[i-1]-Q.sim[i-1])/scale,df=df),mean=0,sd=1)
                eta.2 <- rnorm(n=1, mean=eta.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + qt(pnorm(eta.2), df=df) * scale, 0)
            }else{
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + rt(n=1, df=df) * scale, 0)
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la_skewt_sample <- function(par.model, run.model, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf ){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        samp.var <- Q.sim
        scale <- sd[1]/sd.t
        samp.var[1] <- max(samp.var[1] + rt(n=1, df=df) * scale, 0)
        for ( i in 2:n ){
            if (samp.var[i-1]>0 ){
                delta <- (t.mod[i]-t.mod[i-1])/tau
                eta.1 <- qnorm(pskt((samp.var[i-1]-Q.sim[i-1])/scale,df=df,gamma=gamma),mean=0,sd=1)
                eta.2 <- rnorm(n=1, mean=eta.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + qskt(pnorm(eta.2), df=df,gamma=gamma) * scale, 0)
            }else{
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + rskt(n=1, df=df,gamma=gamma) * scale, 0)
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la2_sample <- function(par.model, run.model, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")] + 1 ## ATTENTION: note the +1 here
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        samp.var <- Q.sim
        scale <- sd[1]/sd.t
        samp.var[1] <- max(samp.var[1] + rt(n=1, df=df) * scale, 0)
        for ( i in 2:n ){
            if(Q.sim[i]== 0 & Q.sim[i-1]==0){
                dec <- FALSE
            }else{
                dec <- ((1/(t.mod[i]-t.mod[i-1])) *(Q.sim[i]/Q.sim[i-1])) < d
            }
            if (dec  & samp.var[i-1]>0 ){
                delta <- (t.mod[i]-t.mod[i-1])/tau
                eta.1 <- qnorm(pt((samp.var[i-1]-Q.sim[i-1])/scale,df=df),mean=0,sd=1)
                eta.2 <- rnorm(n=1, mean=eta.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + qt(pnorm(eta.2), df=df) * scale, 0)
            }else{
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + rt(n=1, df=df) * scale, 0)
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la2_skewt_sample <- function(par.model, run.model, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")] + 1 ## ATTENTION: note the +1 here
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf ){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        samp.var <- Q.sim
        scale <- sd[1]/sd.t
        samp.var[1] <- max(samp.var[1] + rt(n=1, df=df) * scale, 0)
        for ( i in 2:n ){
            if(Q.sim[i]== 0 & Q.sim[i-1]==0){
                dec <- FALSE
            }else{
                dec <- ((1/(t.mod[i]-t.mod[i-1])) *(Q.sim[i]/Q.sim[i-1])) < d
            }
            if (dec  & samp.var[i-1]>0 ){
                delta <- (t.mod[i]-t.mod[i-1])/tau
                eta.1 <- qnorm(pskt((samp.var[i-1]-Q.sim[i-1])/scale,df=df,gamma=gamma),mean=0,sd=1)
                eta.2 <- rnorm(n=1, mean=eta.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + qskt(pnorm(eta.2), df=df,gamma=gamma) * scale, 0)
            }else{
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + rskt(n=1, df=df,gamma=gamma) * scale, 0)
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la3_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")] ## ATTENTION: note the +1 here
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf ){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        sd <- a*sd0*(Q.sim/Q0)^c + b
        samp.var <- Q.sim
        scale <- sd[1:2]/sd.t
        samp.var[1:2] <- pmax(samp.var[1:2] + rskt(n=2, df=df,gamma=gamma) * scale, 0)
        for ( i in 3:n ){
            if(Q.sim[i]== 0 & Q.sim[i-1]==0 & Q.sim[i-2]==0){
                dec <- FALSE
            }else{
                dec <- P.var[i]<= d & P.var[i-1]<=d & P.var[i-2]<=d
            }
            if (dec  & samp.var[i-1]>0 ){
                delta <- (t.mod[i]-t.mod[i-1])/tau
                eta.1 <- qnorm(pskt((samp.var[i-1]-Q.sim[i-1])/scale[length(scale)],df=df,gamma=gamma),mean=0,sd=1)
                eta.2 <- rnorm(n=1, mean=eta.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + qskt(pnorm(eta.2), df=df,gamma=gamma) * scale, 0)
            }else{
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + rskt(n=1, df=df,gamma=gamma) * scale, 0)
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la5_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        sigma.slpe <- par.likeli[paste(var.curr, "_sigmaslpe_lik", sep = "")]
        sigma.min  <- par.likeli[paste(var.curr, "_sigmamin_lik", sep = "")]
        m          <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf | gamma != 1){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        sd <- a*sd0*(Q.sim/Q0)^c + b
        samp.var <- Q.sim
        scale <- sd[1:2]/sd.t
        samp.var[1:2] <- pmax(samp.var[1:2] + rskt(n=2, df=df, gamma=gamma) * scale, 0)
        eta   <- rep(NA,n)
        state <- rep(NA,n)
        x <- rep(NA,n)
        x[2] <- 0
        eta[2] <- rnorm(1, mean=0, sd=1)
        state[2] <- 0
        for ( i in 3:n ){
            if(Q.sim[i]== 0 & Q.sim[i-1]==0 & Q.sim[i-2]==0){
                dec <- FALSE
            }else{
                dec <- P.var[i]<= d & P.var[i-1]<=d & P.var[i-2]<=d
            }
            state[i] <- state[i-1] + P.var[i-1] - psi*state[i-1]
            sigma.slpe.tmp <- sigma.slpe*state[i] + sigma.min
            if (dec  & samp.var[i-1]>0 ){
                x[i] <- x[i-1] + m*eta[i-1]
                eta[i] <- rnorm(n=1, mean=eta[i-1] - x[i], sd=sigma.slpe.tmp)
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + qskt(pnorm(eta[i]), df=df,gamma=gamma) * scale, 0)
            }else{
                scale <- sd[i]/sd.t
                eta[i] <- rnorm(1, mean=0, sd=1)
                samp.var[i] <- max(samp.var[i] + qskt(pnorm(eta[i]), df=df,gamma=gamma) * scale, 0)
                x[i] <- 0
            }
        }
        samp[ind.var] <- samp.var
    }
    return(samp)
}

LogLikelihoodHydrology_la7_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess,...){## la7 puts a deterministic exponential decay with a certain shift on top of the deterministic hydrological model. The quantiles are assumed to be OU-deviations from that combination.
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")] ## ATTENTION: note the +1 here
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        damp   <- par.likeli[paste(var.curr, "_damp_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        mu.tmp <- rep(mu, times=rep.mu.times)
        ## calculate actual value of likelihood
        sd.t <- 1
        if ( df < Inf ){## calculate sd of skewed t distribution
            myfun1 <- function(x,df){2*x*dt(x,df)}
            m1 <- integrate(f=myfun1, lower=0, upper=Inf, df=df)$value
            m2 <- df/(df-2)
            Ex <-  m1 * (gamma^2 - 1/(gamma^2)) / (gamma + 1/gamma)
            Ex2 <- m2 * (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            sd.t <- sqrt(Ex2 - Ex^2)
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        sd <- a*sd0*(Q.sim/Q0)^c + b
        exp.det <- c(mu.tmp[1],mu.tmp)*sd*exp(-damp*c(0,time.recess))
        exp.det[is.na(exp.det)] <- 0
        Q.sim[2:n] <- pmax(Q.sim[2:n] + exp.det, 0)
        sd <- a*sd0*(Q.sim/Q0)^c + b
        samp.var <- Q.sim
        scale <- sd[1:2]/sd.t
        samp.var[1:2] <- pmax(samp.var[1:2] + rskt(n=2, df=df,gamma=gamma) * scale, 0)
        for ( i in 3:n ){
            if(Q.sim[i-1]==0){
                dec <- FALSE
            }else{
                dec <- P.var[i]<= d & P.var[i-1]<=d & P.var[i-2]<=d
            }
            if (dec  & samp.var[i-1]>0 ){
                delta <- (t.mod[i]-t.mod[i-1])/tau
                eta.1 <- qnorm(pskt((samp.var[i-1]-Q.sim[i-1])/scale[length(scale)],df=df,gamma=gamma),mean=0,sd=1)
                eta.2 <- rnorm(n=1, mean=eta.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + qskt(pnorm(eta.2), df=df,gamma=gamma) * scale, 0)
            }else{
                scale <- sd[i]/sd.t
                samp.var[i] <- max(samp.var[i] + rskt(n=1, df=df,gamma=gamma) * scale, 0)
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la8_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        sd01<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        mn.rev <- par.likeli[paste(var.curr, "_mnrev_lik", sep = "")]
        sigmaslpe <- par.likeli[paste(var.curr, "_sigmaslpe_lik", sep = "")]
        sd.min     <- par.likeli[paste(var.curr, "_sdmin_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 2
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 2
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd <- a*sd0*(Q.sim/Q0)^c + b
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1:2], sd2[1:2]/sd.t2, sd1[1:2]/sd.t1)
        delta <- (t.mod[2]-t.mod[1])/ifelse(auto[2], tau2, tau1)
        quant.1 <- rnorm(1)
        quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
        eta.1 <- quant.2 - quant.1
        df.tmp <- ifelse(auto[1:2], df2, df1)
        gamma.tmp <- ifelse(auto[1:2], gamma2, gamma1)
        Exb.tmp <- ifelse(auto[1:2], Exb2, Exb1)
        samp.var[1:2] <- pmax(samp.var[1:2] + (qskt(pnorm(c(quant.1,quant.2)), df=df.tmp,gamma=gamma.tmp)-Exb.tmp)*scale, 0)
        for ( i in 3:n ){
            if (samp.var[i-1]>0 & auto[i]==auto[i-1]){
                tau <- ifelse(auto[i], tau2, tau1)
                delta <- (t.mod[i]-t.mod[i-1])/tau
                if(!auto[i]){
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])/scale[length(scale)]+Exb1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                    scale <- sd1[i]/sd.t1
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1)-Exb1) * scale, 0)
                }else{
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])/scale[length(scale)]+Exb2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    mu <- -mn.rev*eta.1
                    eta.2 <- rnorm(n=1, mean=mu+(eta.1-mu)*exp(-delta), sd=sigmaslpe*sqrt(1-exp(-2*delta)))
                    quant.2 <- quant.1 + eta.2
                    scale <- sd2[i]/sd.t2
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2)-Exb2) * scale, 0)
                    eta.1 <- eta.2
                }
            }else{
                if(!auto[i]){
                    scale <- sd1[i]/sd.t1
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1)-Exb1) * scale, 0)
                }else{
                    eta.1 <- rnorm(1, 0, sigmaslpe)
                    quant.2 <- quant.1 + eta.1
                    scale <- sd2[i]/sd.t2
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2)-Exb2) * scale, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la9_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        ## gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd <- a*sd0*(Q.sim/Q0)^c + b
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        delta <- (t.mod[2]-t.mod[1])/ifelse(auto[2], tau2, tau1)
        quant.1 <- rnorm(1)
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        Exb.tmp <- ifelse(auto[1], Exb2, Exb1)
        samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp)-Exb.tmp)*scale, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                tau <- ifelse(auto[i], tau2, tau1)
                delta <- (t.mod[i]-t.mod[i-1])/tau
                if(!auto[i]){
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1]+Exb1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    scale <- sd2[i]/sd.t2
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1]+Exb2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
                if(is.infinite(samp.var[i])){print(quant.1);print(quant.2);print(auto[i]);print(sd2[i]);stop("infinites in prediction")}
            }else{
                if(!auto[i]){
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la9b_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        ## gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- pmax(Q.sim , 0)
        Q.sim <- mu.tmp*Q.sim
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 ##+ pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd1 <- rollmean(sd1, k=2, fill=sd1[length(sd1)])
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd2 <- rollmean(sd2, k=2, fill=sd2[length(sd2)])
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        S <- rep(NA,n)
        tau.max <- 20
        tau.min <- 0
        k <- 0.1
        RC <- 0.6
        alpha <- 1
        S0 <- mean(P)*RC/k
        S[1] <- S0
        for(i in 2:n){
            S[i] <- S[i-1]*exp(-k*(t.mod[i]-t.mod[i-1])) + (P[i-1]+P[i])/2 ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
        }
        taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-1/alpha*S/S0)
        gmms <- 1/taus
        delta1 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        delta2 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        quant.1 <- rnorm(1)
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        Exb.tmp <- ifelse(auto[1], Exb2, Exb1)
        samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp)-Exb.tmp)*scale, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                if(!auto[i]){
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1]+Exb1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta1[i-1]), sd=sqrt(1-exp(-2*delta1[i-1])))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    scale <- sd2[i]/sd.t2
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1]+Exb2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta2[i-1]), sd=sqrt(1-exp(-2*delta2[i-1])))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
                if(is.infinite(samp.var[i])){print(quant.1);print(quant.2);print(auto[i]);print(sd2[i]);stop("infinites in prediction")}
            }else{
                if(!auto[i]){
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la9c_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ttP<- par.likeli[paste(var.curr, "_ttP_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- pmax(Q.sim , 0)
        Q.sim <- mu.tmp*Q.sim
        sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ## sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd1 <- rollmean(sd1, k=2, fill=sd1[length(sd1)])
        ## sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd2 <- rollmean(sd2, k=2, fill=sd2[length(sd2)])
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        S <- rep(NA,n)
        ##tau.min <- 0
        RC <- 0.6
        alpha <- 1
        S0 <- mean(P)*RC/tkP
        S[1] <- S0
        ## for(i in 2:n){
        ##      S[i] <- S[i-1]*exp(-tkP*(t.mod[i]-t.mod[i-1])) + (P[i-1]+P[i])/2 ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
        ## }
        ## taus <- tau.min + (tau.max-tau.min)*exp(-P/S-Q.sim/Q01*tkQ)
        taus <- numeric(n)
        taus <- tau.min + (tau.max-tau.min)*exp(-pmax(P-Plim,0)/Q.sim*tkP)
        taus <- rollmean(taus, k=3,fill=0)
        ## taus[1] <- tau.min + (tau.max-tau.min)*exp(-(pmax(P[1]-Plim,0)/Q.sim[1])*tkP-Q.sim[1]/Q01*tkQ)
        ## taus[2:n] <- tau.min + (tau.max-tau.min)*exp(-(pmax(0.5*(P[2:n]+P[1:(n-1)])-Plim,0)/0.5/(Q.sim[2:n]+Q.sim[1:(n-1)]))*tkP-0.5*(Q.sim[2:n]+Q.sim[1:(n-1)])/Q01*tkQ)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-pmax(P-ttP,0)/Q.sim*tkP-Q.sim/Q01*tkQ)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-Q.sim/Q01/3)
        gmms <- 1/taus
        dt <- (t.mod[2:n]-t.mod[1:(n-1)])
        sgmm <- gmms[2:n] + gmms[1:(n-1)]
        ## delta1 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        ## delta2 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        quant.1 <- rnorm(1)
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        Exb.tmp <- ifelse(auto[1], Exb2, Exb1)
        ## samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp)-Exb.tmp)*scale, 0)
        samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp))*scale, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                if(!auto[i]){
                    ## quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1]+Exb1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1],df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-sgmm[i-1]/2*dt[i-1]), sd=sqrt(1-exp(-sgmm[i-1]*dt[i-1])))
                    ## samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1))*sd1[i]/sd.t1, 0)
                }else{
                    scale <- sd2[i]/sd.t2
                    ## quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1]+Exb2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1],df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-sgmm[i-i]/2*dt[i-1]), sd=sqrt(1-exp(-sgmm[i-1]*dt[i-1])))
                    ## samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2))*sd2[i]/sd.t2, 0)
                }
                if(is.infinite(samp.var[i])){print(quant.1);print(quant.2);print(auto[i]);print(sd2[i]);stop("infinites in prediction")}
            }else{
                if(!auto[i]){
                    ## samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1))*sd1[i]/sd.t1, 0)
                }else{
                    ## samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2))*sd2[i]/sd.t2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la9d_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ttP<- par.likeli[paste(var.curr, "_ttP_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        Q.sim <- pmax(Q.sim , 0)
        ## dQsim <- pmax(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])),0)
        ## dQpos1 <- pmax(rollmean(dQsim, k=3, fill=0, align="left"), 0)
        S <- numeric(length(Q.sim))
        S[1] <- mean(P)*0.6/psi
        P <- pmax(P,0)
        for(i in 2:n){
            S[i] <- S[i-1]*exp(-psi*(t.mod[i]-t.mod[i-1])) + P[i] ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
            ## S[i] <- S[i-1]/(1+S[i-1]*psi*(t.obs[i]-t.obs[i-1])) + P[i] ## in case alpha = 2
        }
        S <- pmax(S,0)
        sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1 + ifelse(P>0,d*P^e,0)##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2 + ifelse(P>0,d*P^e,0)##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ## sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd1 <- rollmean(sd1, k=2, fill=sd1[length(sd1)])
        ## sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd2 <- rollmean(sd2, k=2, fill=sd2[length(sd2)])
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        ## for(i in 2:n){
        ##      S[i] <- S[i-1]*exp(-tkP*(t.mod[i]-t.mod[i-1])) + (P[i-1]+P[i])/2 ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
        ## }
        ## taus <- tau.min + (tau.max-tau.min)*exp(-P/S-Q.sim/Q01*tkQ)
        P.rel <- ifelse(P==0, 0, P/S/psi)
        taus <- tau.min + (tau.max-tau.min)*exp(-tkP*P.rel^l - tkQ*(S*psi)^m)
        ## taus <- numeric(n)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-tkP*(pmax(P,0)/(Q.sim+Q01*b1))^l - tkQ*(pmax(P,0))^m)
        ## taus[1] <- tau.min + (tau.max-tau.min)*exp(-(pmax(P[1]-Plim,0)/Q.sim[1])*tkP-Q.sim[1]/Q01*tkQ)
        ## taus[2:n] <- tau.min + (tau.max-tau.min)*exp(-(pmax(0.5*(P[2:n]+P[1:(n-1)])-Plim,0)/0.5/(Q.sim[2:n]+Q.sim[1:(n-1)]))*tkP-0.5*(Q.sim[2:n]+Q.sim[1:(n-1)])/Q01*tkQ)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-pmax(P-ttP,0)/Q.sim*tkP-Q.sim/Q01*tkQ)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-Q.sim/Q01/3)
        gmms <- 1/taus
        dt <- (t.mod[2:n]-t.mod[1:(n-1)])
        gmms <- gmms[2:n]
        ## delta1 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        ## delta2 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        quant.1 <- rnorm(1)
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        Exb.tmp <- ifelse(auto[1], Exb2, Exb1)
        ## samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp)-Exb.tmp)*scale, 0)
        samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp))*scale, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                if(!auto[i]){
                    ## quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1]+Exb1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1],df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-gmms[i-1]*dt[i-1]), sd=sqrt(1-exp(-2*gmms[i-1]*dt[i-1])))
                    ## samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1))*sd1[i]/sd.t1, 0)
                }else{
                    scale <- sd2[i]/sd.t2
                    ## quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1]+Exb2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1],df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-gmms[i-i]*dt[i-1]), sd=sqrt(1-exp(-2*gmms[i-1]*dt[i-1])))
                    ## samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2))*sd2[i]/sd.t2, 0)
                }
                if(is.infinite(samp.var[i])){print(quant.1);print(quant.2);print(auto[i]);print(sd2[i]);stop("infinites in prediction")}
            }else{
                if(!auto[i]){
                    ## samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1))*sd1[i]/sd.t1, 0)
                }else{
                    ## samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2))*sd2[i]/sd.t2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la9e_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd01 <- Q01/5 ## in case Q01 is fitted
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ttP<- par.likeli[paste(var.curr, "_ttP_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
                m1 <- 1
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
                m2 <- 1
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
            ## Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1) ## for further use, this is the mean of the unscaled (original) skewed t-distribution
            ## Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            Exb1 <- 0 ## ATTENTION: this means that Qdet is at the mode, not the mean...
            Exb2 <- 0
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        Q.sim <- pmax(Q.sim , 0)
        ## dQsim <- pmax(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])),0)
        ## dQpos1 <- pmax(rollmean(dQsim, k=3, fill=0, align="left"), 0)
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1 + ifelse(P>0,d*P^e,0)##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2 + ifelse(P>0,d*P^e,0)##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + Q01*b1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + Q02*b2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*Q.sim[Q.sim<Q01]/Q01 + b1*Q01
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01/c1*(c1-1+(Q.sim[Q.sim>=Q01]/Q01)^c1)
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*Q.sim[Q.sim<Q02]/Q02 + b2*Q02
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02/c2*(c2-1+(Q.sim[Q.sim>=Q02]/Q02)^c2)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*(Q.sim[Q.sim<Q01]/Q01)^c1 + b1
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01*(1-c1^2+c1^2*(Q.sim[Q.sim>=Q01]/Q01)^(1/c1))
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*(Q.sim[Q.sim<Q02]/Q02)^c2 + b2
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02*(1-c2^2+c2^2*(Q.sim[Q.sim>=Q02]/Q02)^(1/c2))

        ## sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd1 <- rollmean(sd1, k=2, fill=sd1[length(sd1)])
        ## sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd2 <- rollmean(sd2, k=2, fill=sd2[length(sd2)])
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        ## for(i in 2:n){
        ##      S[i] <- S[i-1]*exp(-tkP*(t.mod[i]-t.mod[i-1])) + (P[i-1]+P[i])/2 ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
        ## }
        ## taus <- tau.min + (tau.max-tau.min)*exp(-P/S-Q.sim/Q01*tkQ)
        P <- pmax(P,0)
        taus <- tau.min + (tau.max-tau.min)*exp(-tkP*P^l - tkQ*Q.sim^m)
        ## taus[((1:length(taus)) %% 20) == 0] <- 0
        ## taus <- numeric(n)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-tkP*(pmax(P,0)/(Q.sim+Q01*b1))^l - tkQ*(pmax(P,0))^m)
        ## taus[1] <- tau.min + (tau.max-tau.min)*exp(-(pmax(P[1]-Plim,0)/Q.sim[1])*tkP-Q.sim[1]/Q01*tkQ)
        ## taus[2:n] <- tau.min + (tau.max-tau.min)*exp(-(pmax(0.5*(P[2:n]+P[1:(n-1)])-Plim,0)/0.5/(Q.sim[2:n]+Q.sim[1:(n-1)]))*tkP-0.5*(Q.sim[2:n]+Q.sim[1:(n-1)])/Q01*tkQ)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-pmax(P-ttP,0)/Q.sim*tkP-Q.sim/Q01*tkQ)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-Q.sim/Q01/3)
        gmms <- 1/taus
        dt <- (t.mod[2:n]-t.mod[1:(n-1)])
        gmms <- gmms[2:n]
        ## delta1 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        ## delta2 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        quant.1 <- rnorm(1)
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        Exb.tmp <- ifelse(auto[1], Exb2, Exb1)
        ## samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp)-Exb.tmp)*scale, 0)
        samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp))*scale, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                if(!auto[i]){
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1]+Exb1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-gmms[i-1]*dt[i-1]), sd=sqrt(1-exp(-2*gmms[i-1]*dt[i-1])))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    scale <- sd2[i]/sd.t2
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1]+Exb2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-gmms[i-i]*dt[i-1]), sd=sqrt(1-exp(-2*gmms[i-1]*dt[i-1])))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
                if(is.infinite(samp.var[i])){print(quant.1);print(quant.2);print(auto[i]);print(sd2[i]);stop("infinites in prediction")}
            }else{
                if(!auto[i]){
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_la9esimp_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    options(warn=2)
    ## if(is.null(layout$pred)){
    layout$layout <- layout$layout[c(layout$calib,layout$pred),]
    L <- layout$layout
    P <- P[c(layout$calib,layout$pred)]
    ## }else{
    ##     L             <- layout$layout[layout$pred,]
    ##     layout$layout <- L
    ## }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd01 <- Q01/5 ## in case Q01 is fitted
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ttP<- par.likeli[paste(var.curr, "_ttP_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( any(df1 < Inf) | any(gamma1 != 1) | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(any(df1<Inf)){
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
            }else{
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
            }else{
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            chunk1 <- sqrt(Ex1 - Exb1^2)
            chunk2 <- sqrt(Ex2 - Exb2^2)
            Exb1 <- 0 ## ATTENTION: this means that Qdet is at the mode. If we remove these two lines, Qdet is at the mean.
            Exb2 <- 0
        }else{
            chunk1 <- 1
            chunk2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        Q.sim <- pmax(Q.sim , 0)
        ## dQsim <- pmax(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])),0)
        ## dQpos1 <- pmax(rollmean(dQsim, k=3, fill=0, align="left"), 0)
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1 + ifelse(P>0,d*P^e,0)##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2 + ifelse(P>0,d*P^e,0)##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + Q01*b1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + Q02*b2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*Q.sim[Q.sim<Q01]/Q01 + b1*Q01
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01/c1*(c1-1+(Q.sim[Q.sim>=Q01]/Q01)^c1)
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*Q.sim[Q.sim<Q02]/Q02 + b2*Q02
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02/c2*(c2-1+(Q.sim[Q.sim>=Q02]/Q02)^c2)

        ## sd1 <- numeric(length(Q.sim))
        ## sd1[Q.sim<Q01] <-  a1*sd01*(Q.sim[Q.sim<Q01]/Q01)^c1 + b1
        ## sd1[Q.sim>=Q01] <- b1*Q01 + a1*sd01*(1-c1^2+c1^2*(Q.sim[Q.sim>=Q01]/Q01)^(1/c1))
        ## sd2 <- numeric(length(Q.sim))
        ## sd2[Q.sim<Q02] <-  a2*sd02*(Q.sim[Q.sim<Q02]/Q02)^c2 + b2
        ## sd2[Q.sim>=Q02] <- b2*Q02 + a2*sd02*(1-c2^2+c2^2*(Q.sim[Q.sim>=Q02]/Q02)^(1/c2))

        ## sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 + pmax(0.5*rollmean(c(0,diff(Q.sim)),k=2,fill=0), 0)
        ##sd1 <- rollmean(sd1, k=2, fill=sd1[length(sd1)])
        ## sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd2 <- rollmean(sd2, k=2, fill=sd2[length(sd2)])
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        sdt <- ifelse(auto, sd2/chunk2, sd1/chunk1)
        ## for(i in 2:n){
        ##      S[i] <- S[i-1]*exp(-tkP*(t.mod[i]-t.mod[i-1])) + (P[i-1]+P[i])/2 ## ATTENTION the times of P have to be at the layout times in order for this to work, but usualy they are at the input times...
        ## }
        ## taus <- tau.min + (tau.max-tau.min)*exp(-P/S-Q.sim/Q01*tkQ)
        P <- pmax(P-Plim,0)
        taus <- tau.min + (tau.max-tau.min)*exp(-tkP*P^l - tkQ*Q.sim^m)
        ## taus[((1:length(taus)) %% 20) == 0] <- 0
        ## taus <- numeric(n)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-tkP*(pmax(P,0)/(Q.sim+Q01*b1))^l - tkQ*(pmax(P,0))^m)
        ## taus[1] <- tau.min + (tau.max-tau.min)*exp(-(pmax(P[1]-Plim,0)/Q.sim[1])*tkP-Q.sim[1]/Q01*tkQ)
        ## taus[2:n] <- tau.min + (tau.max-tau.min)*exp(-(pmax(0.5*(P[2:n]+P[1:(n-1)])-Plim,0)/0.5/(Q.sim[2:n]+Q.sim[1:(n-1)]))*tkP-0.5*(Q.sim[2:n]+Q.sim[1:(n-1)])/Q01*tkQ)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-pmax(P-ttP,0)/Q.sim*tkP-Q.sim/Q01*tkQ)
        ## taus <- tau.min + (tau.max-tau.min)*exp(-alpha*P/(k*S)-Q.sim/Q01/3)
        gmms <- 1/taus
        dt <- (t.mod[2:n]-t.mod[1:(n-1)])
        gmms <- gmms[2:n]
        ## delta1 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        ## delta2 <- (t.mod[2:n] - t.mod[1:(n-1)])*gmms[1:(n-1)]
        quant.1 <- rnorm(1)
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        Exb.tmp <- ifelse(auto[1],Exb2[1],Exb1[1])
        chunk.tmp <- ifelse(auto[1],chunk2[1],chunk1[1])
        samp.var[1] <- pmax(samp.var[1] + myqskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp,sigm=sdt[1])-sdt[1]*Exb.tmp, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                if(!auto[i]){
                    et.old <- samp.var[i-1]-Q.sim[i-1]+Exb1*sdt[i-1]
                    if(et.old>=0){
                        quant.1 <- qnorm(mypskt(et.old,df=df1,gamma=gamma1,sigm=sdt[i-1],low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    }else{
                        quant.1 <- qnorm(mypskt(et.old,df=df1,gamma=gamma1,sigm=sdt[i-1],log.p=TRUE),mean=0,sd=1,log.p=TRUE)
                    }
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-gmms[i-1]*dt[i-1]), sd=sqrt(1-exp(-2*gmms[i-1]*dt[i-1])))
                    samp.var[i] <- max(samp.var[i] + myqskt(pnorm(quant.2), df=df1,gamma=gamma1,sigm=sdt[i])-sdt[i]*Exb1, 0)
                }else{
                    et.old <- samp.var[i-1]-Q.sim[i-1]+Exb2*sdt[i-1]
                    if(et.old>=0){
                        quant.1 <- qnorm(mypskt(et.old,df=df2,gamma=gamma2,sigm=sdt[i-1],low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    }else{
                        quant.1 <- qnorm(mypskt(et.old,df=df2,gamma=gamma2,sigm=sdt[i-1],log.p=TRUE),mean=0,sd=1,log.p=TRUE)
                    }
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-gmms[i-i]*dt[i-1]), sd=sqrt(1-exp(-2*gmms[i-1]*dt[i-1])))
                    samp.var[i] <- max(samp.var[i] + myqskt(pnorm(quant.2), df=df2,gamma=gamma2,sigm=sdt[i])-sdt[i]*Exb2, 0)
                }
                if(is.infinite(samp.var[i])){print(quant.1);print(quant.2);print(auto[i]);print(sd2[i]);stop("infinites in prediction")}
            }else{
                if(!auto[i]){
                    samp.var[i] <- max(samp.var[i] + myqskt(pnorm(rnorm(1)),df=df1,gamma=gamma1,sigm=sdt[i])-sdt[i]*Exb1, 0)
                }else{
                    samp.var[i] <- max(samp.var[i] + myqskt(pnorm(rnorm(1)),df=df2,gamma=gamma2,sigm=sdt[i])-sdt[i]*Exb2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    options(warn=0)
    return(samp)
}


LogLikelihoodHydrology_la92_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        ## Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd01_lik", sep = "")]
        ## sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        ## a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        ## b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        ## tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        tau2 <- tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        ## mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( any(df1 < Inf) | any(gamma1 != 1) | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(any(df1<Inf)){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
                m1 <- 1
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
                m2 <- 1
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
            ## Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1) ## for further use, this is the mean of the unscaled (original) skewed t-distribution
            ## Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            Exb1 <- 0 ## ATTENTION: this means that the Qdet is at the mode, not the mean...
            Exb2 <- 0
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        ## mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        ## Q.sim <- mu.tmp*Q.sim
        Q.sim <- pmax(Q.sim, 0)
        ## sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ## sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + Q01*b1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + Q02*b2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        ## sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1 + d*pmax(P-Plim,0)^e
        ## sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2 + d*pmax(P-Plim,0)^e
        ##sd <- a*sd0*(Q.sim/Q0)^c + b
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        delta <- (t.mod[2]-t.mod[1])/ifelse(auto[2], tau2, tau1)
        quant.1 <- rnorm(1)
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        Exb.tmp <- ifelse(auto[1], Exb2, Exb1)
        samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp)-Exb.tmp)*scale, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                tau <- ifelse(auto[i], tau2, tau1)
                delta <- (t.mod[i]-t.mod[i-1])/tau
                if(!auto[i]){
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1]+Exb1,df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    scale <- sd2[i]/sd.t2
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1]+Exb2,df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
                if(is.infinite(samp.var[i])){print(quant.1);print(quant.2);print(auto[i]);print(sd2[i]);stop("infinites in prediction")}
            }else{
                if(!auto[i]){
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}


LogLikelihoodHydrology_la9_compound_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ## c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        ## gamma2 <- gamma1
        pout1 <- par.likeli[paste(var.curr, "_p_lik", sep = "")]
        pout2 <- par.likeli[paste(var.curr, "_p2_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if(df1<Inf){
            m1 <- df1/(df1-2)
            Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
        }else{
            m1 <- 1
            Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
        }
        if(df2<Inf){
            m2 <- df2/(df2-2)
            Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
        }else{
            m2 <- 1
            Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
        }
        Ez1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
        Ez2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
        dfout <- 5
        Ek1 <- dfout/(dfout-2) ## 2.5 are the degrees of freedom of the outer t-distribution used to describe the outliers
        Ek2 <- dfout/(dfout-2)
        Ex1 <- (1-pout1)*Ez1 + pout1*Ek1
        Ex2 <- (1-pout2)*Ez2 + pout2*Ek2
        sd.t1 <- sqrt(Ex1 - ((1-pout1)*Exb1)^2)
        sd.t2 <- sqrt(Ex2 - ((1-pout2)*Exb2)^2)
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        sd1 <- a1*sd01*((Q.sim/Q01+b1)/(1+b1))^c1
        sd2 <- a2*sd02*((Q.sim/Q02+b2)/(1+b2))^c2
        ##sd <- a*sd0*(Q.sim/Q0)^c + b
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        delta <- (t.mod[2]-t.mod[1])/ifelse(auto[2], tau2, tau1)
        quant.1 <- rnorm(1)
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        pout.tmp <- ifelse(auto[1], pout2, pout1)
        Exb.tmp <- ifelse(auto[1], Exb2, Exb1)
        samp.var[1] <- pmax(samp.var[1] + (qcompound(pnorm(quant.1), df=df.tmp,df2=dfout,gamma=gamma.tmp,p=pout.tmp)-Exb.tmp)*scale, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                tau <- ifelse(auto[i], tau2, tau1)
                delta <- (t.mod[i]-t.mod[i-1])/tau
                if(!auto[i]){
                    quant.1 <- qnorm(pcompound((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1]+Exb1,df=df1,df2=dfout,gamma=gamma1,p=pout1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                    samp.var[i] <- max(samp.var[i] + (qcompound(pnorm(quant.2), df=df1,df2=dfout,gamma=gamma1,p=pout1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    scale <- sd2[i]/sd.t2
                    quant.1 <- qnorm(pcompound((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1]+Exb2,df=df2,df2=dfout,gamma=gamma2,p=pout2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                    samp.var[i] <- max(samp.var[i] + (qcompound(pnorm(quant.2), df=df2,df2=dfout,gamma=gamma2,p=pout2)-Exb2)*sd2[i]/sd.t2, 0)
                }
            }else{
                if(!auto[i]){
                    samp.var[i] <- max(samp.var[i] + (rcompound(n=1, df=df1,df2=dfout,gamma=gamma1,p=pout1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    samp.var[i] <- max(samp.var[i] + (rcompound(n=1, df=df2,df2=dfout,gamma=gamma2,p=pout2)-Exb2)*sd2[i]/sd.t2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}


LogLikelihoodHydrology_la9_skewt_sample_synth<- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, eta, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        ##Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        ##sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        ##a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        ##b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ##c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ##df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ##gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            Exb1 <- 0
            Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd <- a*sd0*(Q.sim/Q0)^c + b
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        delta <- (t.mod[2]-t.mod[1])/ifelse(auto[2], tau2, tau1)
        quant.1 <- eta[1]
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        Exb.tmp <- ifelse(auto[1], Exb2, Exb1)
        samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp)-Exb.tmp)*scale, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                tau <- ifelse(auto[i], tau2, tau1)
                delta <- (t.mod[i]-t.mod[i-1])/tau
                if(!auto[i]){
                    quant.2 <- eta[i]
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    scale <- sd2[i]/sd.t2
                    quant.2 <- eta[i]
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
            }else{
                if(!auto[i]){
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1)-Exb1)*sd1[i]/sd.t1, 0)
                }else{
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2)-Exb2)*sd2[i]/sd.t2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}


LogLikelihoodHydrology_la9_mode_skewt_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        sd01<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        sd02<- par.likeli[paste(var.curr, "_sd02_lik", sep = "")]
        ## sd02 <- sd01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ##c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        c2 <- c1
        tau1<- par.likeli[paste(var.curr, "_tau1_lik", sep = "")]
        tau2<- par.likeli[paste(var.curr, "_tau2_lik", sep = "")]
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ##df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        ind.var <- L[,1]==var.curr
        mu.tmp <- rep(mu, times=rep.mu.times)
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1){## calculate sd of skewed t distribution
            if(df1<Inf){
                m1 <- df1/(df1-2)
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*df1/(df1-1)
            }else{
                m1 <- 1
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                m2 <- df2/(df2-2)
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*df2/(df2-1)
            }else{
                m2 <- 1
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            Ex1 <- m1 * (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            Ex2 <- m2 * (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            sd.t1 <- sqrt(Ex1 - Exb1^2)
            sd.t2 <- sqrt(Ex2 - Exb2^2)
        }else{
            sd.t1 <- 1
            sd.t2 <- 1
            ##Exb1 <- 0
            ##Exb2 <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        ##if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        ##{ cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate the likelihood based on a t-distribution
        mu.tmp[time.recess<0 | is.infinite(time.recess)] <- 1
        Q.sim <- mu.tmp*Q.sim
        sd1 <- a1*sd01*(Q.sim/Q01)^c1 + b1
        sd2 <- a2*sd02*(Q.sim/Q02)^c2 + b2
        ##sd <- a*sd0*(Q.sim/Q0)^c + b
        samp.var <- Q.sim
        ##scale <- sd[1:2]/sd.t
        scale <- ifelse(auto[1], sd2[1]/sd.t2, sd1[1]/sd.t1)
        delta <- (t.mod[2]-t.mod[1])/ifelse(auto[2], tau2, tau1)
        quant.1 <- rnorm(1)
        df.tmp <- ifelse(auto[1], df2, df1)
        gamma.tmp <- ifelse(auto[1], gamma2, gamma1)
        ##Exb.tmp <- ifelse(auto[1], Exb2, Exb1)
        samp.var[1] <- pmax(samp.var[1] + (qskt(pnorm(quant.1), df=df.tmp,gamma=gamma.tmp))*scale, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                tau <- ifelse(auto[i], tau2, tau1)
                delta <- (t.mod[i]-t.mod[i-1])/tau
                if(!auto[i]){
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t1/sd1[i-1],df=df1,gamma=gamma1,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df1,gamma=gamma1))*sd1[i]/sd.t1, 0)
                }else{
                    scale <- sd2[i]/sd.t2
                    quant.1 <- qnorm(mypskt((samp.var[i-1]-Q.sim[i-1])*sd.t2/sd2[i-1],df=df2,gamma=gamma2,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                    samp.var[i] <- max(samp.var[i] + (qskt(pnorm(quant.2), df=df2,gamma=gamma2))*sd2[i]/sd.t2, 0)
                }
            }else{
                if(!auto[i]){
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df1,gamma=gamma1))*sd1[i]/sd.t1, 0)
                }else{
                    samp.var[i] <- max(samp.var[i] + (rskt(n=1, df=df2,gamma=gamma2))*sd2[i]/sd.t2, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}

LogLikelihoodHydrology_mod_n_sample <- function(par.model, run.model, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        tau<- par.likeli[paste(var.curr, "_tau_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( is.na(tau) | is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        samp.var <- Q.sim
        samp.var[1] <- pmax(rnorm(n=1, mean=Q.sim[1], sd=sd[1]), 0)
        for ( i in 2:n ){
            if ( Q.sim[i]<Q.sim[i-1] & samp.var[i-1]>0 ){
                delta <- (t.mod[i]-t.mod[i-1])/tau
                eta.1 <- qnorm(pnorm((samp.var[i-1]-Q.sim[i-1]),mean=0,sd=sd[i-1]),mean=0,sd=1)
                eta.2 <- rnorm(n=1, mean=eta.1*exp(-delta), sd=sqrt(1-exp(-2*delta)))
                samp.var[i] <- pmax(qnorm(pnorm(eta.2), mean=Q.sim[i], sd=sd[i]), 0)
            }else{
                samp.var[i] <- pmax(rnorm(n=1, mean=Q.sim[i], sd=sd[i]), 0)
            }
        samp[ind.var] <- samp.var
        }
    }
    return(samp)
}


LogLikelihoodHydrology_noautocorr_sample <- function(par.model, run.model, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        sd.t <- 1; if ( df < Inf ) sd.t <- sqrt(df/(df-2))
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if (is.na(a) | is.na(b) | is.na(c) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        samp.var <- Q.sim
        scale <- sd/sd.t
        samp.var <- samp.var + rt(n=n, df=df) * scale
        samp.var <- pmax(samp.var, 0)
        samp[ind.var] <- samp.var
    }
    return(samp)
}

LogLikelihoodHydrology_wls_sample <- function(par.model, run.model, layout, par.likeli, ...){
    if(is.null(layout$pred)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    sample <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if (is.na(a))
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*((Q.sim/Q0+b)/(1+b))^c
        sample[ind.var] <- rnorm(n, mean=Q.sim, sd=sd)
    }
    return(sample)
}

LogLikelihoodHydrology_AR1_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, ...){
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$pred,]
        layout$layout <- L
    }
    y.mod <- as.numeric(run.model(par=par.model, layout=layout, ...)$incld.lmpd)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    loglikeli <- numeric(length=nrow(L))
    samp <- rep(NA, nrow(L))
    for(var.curr in vars){
        Q0 <- par.likeli[paste(var.curr, "_Q0_lik", sep = "")]
        sd0<- par.likeli[paste(var.curr, "_sd0_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        ar1<- par.likeli[paste(var.curr, "_ar1_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        t.obs <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if (is.na(a))
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        sd <- a*sd0*(Q.sim/Q0)^c+b
        res <- rep(NA,length(Q.sim))
        res[1] <- rnorm(1, mean=0, sd=1)
        for(i in 2:length(res)){
            res[i] <- rnorm(1, mean=res[i-1]*ar1, sd=sqrt(1-ar1^2))
        }
        samp.var <- Q.sim + res*sd
        samp[ind.var] <- samp.var
    }
    return(samp)
}

## Auxiliary functions
mydskt <- function (x, df=Inf, gamma=1, sigm=1, log=FALSE){
    result <- rep(NA, length(x))
    if(length(sigm)==1){
        fac <- ifelse(df<Inf,sqrt(df/(df-2)),1)/sigm
        if(log){
            result[x < 0] <- log(2/(gamma + 1/gamma)*fac) + dt(gamma * fac*x[x < 0],df=df,log=TRUE)
            result[x >= 0] <-log(2/(gamma + 1/gamma)*fac) + dt(fac*x[x >= 0]/gamma,df=df,log=TRUE)
        }else{
            result[x < 0] <- 2/(gamma + 1/gamma)*fac * dt(gamma * fac*x[x < 0],df=df,log=FALSE)
            result[x >= 0] <- 2/(gamma + 1/gamma)*fac * dt(fac*x[x >= 0]/gamma,df=df,log=FALSE)
        }
        return(result)
    }else{
        if(length(x)!=length(sigm))stop("different length in mydskt")
        fac1 <- ifelse(df<Inf,sqrt(df/(df-2)),1)/(sigm[x<0])
        fac2 <- ifelse(df<Inf,sqrt(df/(df-2)),1)/(sigm[x>=0])
        if(log){
            result[x < 0] <- log(2/(gamma + 1/gamma)*fac1) + dt(gamma * fac1*x[x < 0],df=df,log=TRUE)
            result[x >= 0] <-log(2/(gamma + 1/gamma)*fac2) + dt(fac2*x[x >= 0]/gamma,df=df,log=TRUE)
        }else{
            result[x < 0] <- 2/(gamma + 1/gamma)*fac1 * dt(gamma * fac1*x[x < 0],df=df,log=FALSE)
            result[x >= 0] <- 2/(gamma + 1/gamma)*fac2 * dt(fac2*x[x >= 0]/gamma,df=df,log=FALSE)
        }
        return(result)
    }
}

mypskt <- function (x, df, gamma = 1, sigm=1, low=TRUE, log.p=FALSE){
    result <- rep(NA, length(x))
    if(length(sigm)==1){
        fac <- ifelse(df<Inf,sqrt(df/(df-2)),1)/sigm
        if(low & !log.p){
            result[x < 0] <- 2/(gamma^2 + 1) * pt(gamma * fac*x[x < 0], df, lower.tail=low, log.p=log.p)
            result[x >= 0] <- 1/(gamma^2 + 1) + 2/(1 + (1/gamma^2)) * (pt(fac*x[x >= 0]/gamma, df, lower.tail=low, log.p=log.p) - 1/2)
        }
        if(low & log.p){
            result[x < 0] <- log(2/(gamma^2+1)) + pt(fac*x[x<0]*gamma, df=df, log.p=TRUE)
            result[x >=0 ] <- log(1/(gamma^2 + 1) + 2/(1 + (1/gamma^2)) * (pt(fac*x[x >= 0]/gamma, df=df) - 1/2))
        }
        if(!low & !log.p){
            result[x < 0] <- 1 - 2/(gamma^2+1)*pt(fac*x[x<0]*gamma, df=df)
            result[x >=0 ] <- 1 - 1/(gamma^2 + 1) - 2/(1 + (1/gamma^2)) * (pt(fac*x[x >= 0]/gamma, df=df) - 1/2)
        }
        if(!low & log.p){
            result[x < 0] <- log(1 - 2/(gamma^2+1)*pt(fac*x[x<0]*gamma, df=df))
            result[x >= 0] <- log(2/((1/gamma^2)+1)) + pt(fac*x[x>=0]/gamma, df=df, lower.tail=FALSE, log.p=TRUE)
        }
        return(result)
    }else{
        fac1 <- ifelse(df<Inf,sqrt(df/(df-2)),1)/(sigm[x<0])
        fac2 <- ifelse(df<Inf,sqrt(df/(df-2)),1)/(sigm[x>=0])
        if(low & !log.p){
            result[x < 0] <- 2/(gamma^2 + 1) * pt(gamma * fac1*x[x < 0], df, lower.tail=low, log.p=log.p)
            result[x >= 0] <- 1/(gamma^2 + 1) + 2/(1 + (1/gamma^2)) * (pt(fac2*x[x >= 0]/gamma, df, lower.tail=low, log.p=log.p) - 1/2)
        }
        if(low & log.p){
            result[x < 0] <- log(2/(gamma^2+1)) + pt(fac1*x[x<0]*gamma, df=df, log.p=TRUE)
            result[x >=0 ] <- log(1/(gamma^2 + 1) + 2/(1 + (1/gamma^2)) * (pt(fac2*x[x >= 0]/gamma, df=df) - 1/2))
        }
        if(!low & !log.p){
            result[x < 0] <- 1 - 2/(gamma^2+1)*pt(fac1*x[x<0]*gamma, df=df)
            result[x >=0 ] <- 1 - 1/(gamma^2 + 1) - 2/(1 + (1/gamma^2)) * (pt(fac2*x[x >= 0]/gamma, df=df) - 1/2)
        }
        if(!low & log.p){
            result[x < 0] <- log(1 - 2/(gamma^2+1)*pt(fac1*x[x<0]*gamma, df=df))
            result[x >= 0] <- log(2/((1/gamma^2)+1)) + pt(fac2*x[x>=0]/gamma, df=df, lower.tail=FALSE, log.p=TRUE)
        }
        return(result)
    }
}

myqskt <- function (x, df, gamma = 1, sigm=1){
    result <- rep(NA, length(x))
    ind1 <- x<mypskt(0,df=df,gamma=gamma,sigm=sigm)
    ind2 <- x>=mypskt(0,df=df,gamma=gamma,sigm=sigm)
    if(length(sigm)==1){
        result[ind1]  <- sigm*ifelse(df<Inf,sqrt(df-2)/sqrt(df),1)/gamma * qt(x[ind1]*(1+gamma^2)/2, df=df)
        result[ind2] <- sigm*gamma*ifelse(df<Inf,sqrt(df-2)/sqrt(df),1) * qt(0.5*(x[ind2]-1/(1+gamma^2))*(1+1/gamma^2) + 0.5, df=df)
        return(result)
    }else{
        sigm1 <- sigm[ind1]
        sigm2 <- sigm[ind2]
        result[ind1]  <- sigm1*ifelse(df<Inf,sqrt(df-2)/sqrt(df),1)/gamma * qt(x[ind1]*(1+gamma^2)/2, df=df)
        result[ind2] <- sigm2*gamma*ifelse(df<Inf,sqrt(df-2)/sqrt(df),1) * qt(0.5*(x[ind2]-1/(1+gamma^2))*(1+1/gamma^2) + 0.5, df=df)
        return(result)
    }
}

dcompound <- function(x, df=Inf, df2=5, gamma=1, p=1e-3, log=FALSE){
    result <- (1-p)*mydskt(x, df=df, gamma=gamma, log=FALSE) + p*mydskt(x, df=df2, gamma=1, log=FALSE)
    if(log) result <- log(result)
    return(result)
}
pcompound <- function(x, df=Inf, df2=5, gamma=1, p=1e-3, low=TRUE, log.p=FALSE){
    result <- (1-p)*mypskt(x, df=df, gamma=gamma, low=low, log.p=FALSE) + p*mypskt(x, df=df2, gamma=1, low=low, log.p=FALSE)
    if(log.p) result <- log(result)
    return(result)
}
rcompound <- function(n, df=Inf, df2=5, gamma=1, p=1e-3){
    k <- runif(n)
    result <- (1-p)*qcompound(k, df=df, df2=df2, gamma=gamma, p=p) + p*qcompound(k, df=df2, df2=df2, gamma=1, p=p)
    return(result)
}
qcompound <- function(prob, df=Inf, df2=5, gamma=1, p=1e-3, low=TRUE, log.p=FALSE){
    if(log.p)  prob <- exp(prob)
    if(!low) prob <- 1-prob
    result <- (1-p)*qskt(prob, df=df, gamma=gamma) + p*qskt(prob, df=df2, gamma=1)
    return(result)
}
