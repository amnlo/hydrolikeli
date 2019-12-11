## ===================================================================================
## Likelihood Evaluation
## ===================================================================================
## These are functions that evaluate the density of the likelihood (and the prior), e.g. for inference

wrap.loglik <- function(param, logposterior, sudriv, scaleshift=NA, mnprm=NA){
    ## This is a wrapper for the logposterior function, to connect it to the time-dependent parameter framework of Peter Reichert.
    ind.timedep <- unlist(lapply(param, length))>1
    any.timedep <- FALSE
    if(sum(ind.timedep)>0){
        any.timedep <- TRUE
        if(!all(names(param)[ind.timedep] == names(sudriv$model$parameters)[as.logical(sudriv$model$timedep$pTimedep)])) stop("order of timdependent parameters is wrong.")
        parmat <- do.call(cbind, param[ind.timedep]) # make a matrix from timedep. parameters in list
        parmat <- parmat[,((1:ncol(parmat)) %% 2) == 0, drop=FALSE] # keep only the values (every second column, order has to agree)
        parmat <- as.matrix(parmat)
        colnames(parmat) <- NULL
        if(!all(is.na(mnprm))){ ## shift mean of some timedep-parameters for which the mean was re-parameterized as a constant parameter
          ind.mnprm <- names(param) %in% mnprm
          if(!all(mnprm %in% names(param))) stop("some parameters of ", mnprm, " not found")
          mnprm <- unlist(param[ind.mnprm])
          for(mn.curr in names(mnprm)){
            ind.parmat <- match(gsub("_fmean","",mn.curr), names(sudriv$model$parameters)[sudriv$model$timedep$pTimedep])
            if(as.logical(sudriv$model$args$parTran[names(sudriv$model$parameters)==gsub("_fmean","",mn.curr)]) | (mn.curr %in% rownames(scaleshift))){
              parmat[,ind.parmat] <- parmat[,ind.parmat] + mnprm[mn.curr] # if the parameter is log-transformed
            }else{
              parmat[,ind.parmat] <- parmat[,ind.parmat] * mnprm[mn.curr] # if it is not log-transformed
            }
          }
          param[ind.mnprm] <- NULL
          ind.timedep <- unlist(lapply(param, length))>1
        }
        if(!all(is.na(scaleshift))){ ## back-transform parameter with sigmoid transformation
          if(ncol(scaleshift)!=2) stop("dimension of scaleshift is not right")
          for(i in 1:ncol(parmat)){
            td.curr <- names(sudriv$model$parameters[sudriv$model$timedep$pTimedep])[i]
            if(td.curr %in% rownames(scaleshift)) parmat[,i] <- sigm.trans(parmat[,i], scale=scaleshift[td.curr,1], shift=scaleshift[td.curr,2])
          }
        }        
        ## force time course within bounds after addition or multiplication with fmean parameter
        lo <- sudriv$model$args$parLo[sudriv$model$timedep$pTimedep]
        hi <- sudriv$model$args$parHi[sudriv$model$timedep$pTimedep]
        for(i in 1:ncol(parmat)){
          parmat[,i] <- pmin(pmax(parmat[,i], lo[i]), hi[i])
        }
        sudriv$model$timedep$par <- parmat
    }
    prs.su.wtd <- c(sudriv$model$parameters[as.logical(sudriv$model$par.fit)], sudriv$likelihood$parameters[as.logical(sudriv$likelihood$par.fit)])
    prs.su     <- prs.su.wtd[!(names(prs.su.wtd) %in% names(param)[ind.timedep])]
    x0 <- numeric(length=sum(c(sudriv$model$par.fit, sudriv$likelihood$par.fit))) ## create the vector of time-constant parameters that is fitted
    prs <- unlist(param[!ind.timedep])
    if(length(prs.su)!=length(prs)) stop("number of supplied constant parameters does not agree with number of fitted parameters")
    prs <- prs[match(names(prs.su), names(prs))]
    if(any(is.na(prs))) stop("sudriv object contains parameters that are designated as fitted but are not supplied by param")
    ind.timedep.su <- names(prs.su.wtd) %in% names(param)[ind.timedep]
    x0[!ind.timedep.su] <- prs
    if(any.timedep) x0[ind.timedep.su]  <- as.numeric(parmat[1,]) ## the value of x0 for the time-dependent parameter should not matter, since it is taken from su$model$timedep$par
    lik <- logposterior(x0=x0, sudriv=sudriv, prior=FALSE) # calculate the log-likelihood
    return(lik)
}
sigm.trans <- function(x, scale=1, shift=0){
    scale/(1+exp(-x)) + shift
}
sigm.trans.inv <- function(x, scale=1, shift=0){
    -log(scale/(x-shift)-1)
}
logposterior <- function(x0, sudriv, prior, mode=TRUE, apprx=FALSE, verbose=TRUE, auto=NA, weight.equally=FALSE){
    flp <- sudriv$likelihood$par.fit
    fmp <- sudriv$model$par.fit
    l.fit <- sum(c(fmp,flp))
    if(length(x0) != l.fit) stop("length of x0 and par.fit in sudriv does not agree")
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
    par.likeli <- sudriv$likelihood$parameters
    par.likeli[which(flp != 0)] <- par.lik.fit

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
    likeli.args$run.model <- run.model
    likeli.args$layout    <- sudriv$layout
    likeli.args$y.obs     <- sudriv$observations
    likeli.args$P         <- sudriv$input$P.roll[sudriv$layout$calib]
    likeli.args$par.likeli<- ifelse(as.logical(sudriv$likelihood$tran), exp(par.likeli), par.likeli)
    names(likeli.args$par.likeli) <- names(sudriv$likelihood$parameters)
    likeli.args$auto      <- auto
    likeli.args$mode      <- mode
    likeli.args$apprx     <- apprx
    likeli.args$sudriv    <- sudriv
    likeli.args$verbose   <- verbose
    likeli.args$weight.equally <- weight.equally
    f.likeli <- sudriv$likelihood$f.likeli

    ## =======================================================
    ## calculate logprior
    if(prior){
        pri.m <- sudriv$model$prior
        pri.m$distdef <- pri.m$distdef[as.logical(fmp)]
        pri.l <- sudriv$likelihood$prior
        pri.l$distdef <- pri.l$distdef[as.logical(flp)]
        pri.hyper <- sudriv$hyperparameters$prior
	if(sum(fmp)>0){
	        args.pdf.model       <- c(list(z=as.numeric(sudriv$model$parameters[as.logical(fmp)])), pri.m)
	        logpri.modelpar      <- do.call(calcpdf_mv, args.pdf.model)
                ## cat("logpri.modelpar: ", logpri.modelpar, "\n")
	}else{
		logpri.modelpar <- 0
	}
	if(sum(flp)>0){
	        args.pdf.likeli      <- c(list(z=as.numeric(par.lik.fit)), pri.l)
		logpri.likelipar     <- do.call(calcpdf_mv, args.pdf.likeli)
                ## cat("logpri.likelipar: ", logpri.likelipar, "\n")
	}else{
		logpri.likelipar <- 0
	}
        if(!is.null(sudriv$hyperparameters)){
            z <- numeric(length=length(sudriv$hyperparameters)-1)
            for(i in 1:length(z)){
                allpr <- c(sudriv$model$parameters, sudriv$likelihood$parameters)
                names(allpr) <- gsub("%", "_", names(allpr))
                z[i] <- with(as.list(allpr), expr=eval(parse(text=sudriv$hyperparameters[[i]]$formula)))
            }
            args.pdf.hyper       <- c(list(z=z), pri.hyper)
            logpri.hyperpar      <- do.call(calcpdf_mv, args.pdf.hyper)
        }else{
            logpri.hyperpar <- 0
        }
    }else{
        logpri.likelipar <- 0
        logpri.modelpar  <- 0
        logpri.hyperpar  <- 0
    }
    ## =======================================================
    ## calculate loglikelihood
    if(is.finite(logpri.modelpar) & is.finite(logpri.likelipar) & is.finite(logpri.hyperpar)){
        tme <- proc.time()
        loglikeli <- do.call(f.likeli, likeli.args)
        ## cat("loglik: ", loglikeli, "\n")
        if(!verbose){return(loglikeli)}
    }else{
        return(-Inf)
    }
    ## =======================================================
    ## calculate logposterior
    logpost <- loglikeli + logpri.likelipar + logpri.modelpar + logpri.hyperpar
    return(logpost)
}

LogLikelihoodHydrology_la9esimp_fast_skewt <- function(run.model, layout, y.obs, P, par.likeli, auto, mode, apprx, verbose, weight.equally, ...){ ## simplified version of la9(e), where we first rescale, and then skew the distribution DQ.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        layout$lump   <- layout$lump[layout$calib]
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(layout=layout, ...)$incld.lmpd)
    if(any(!is.numeric(y.mod))) stop("y.mod contains non-numeric values")
    loglikeli <- numeric(length=nrow(L))
    vars <- unique(L[,1])
    P.all <- P # copies, since they are overwritten later on
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
	if(!any(grepl("GLOB_Mult_.*_a_lik", names(par.likeli)))){
		mult_a <- 1
	}else{
		mult_a <- ifelse(grepl("Wv",var.curr), par.likeli["GLOB_Mult_Q_a_lik"], 1)*ifelse(grepl("Tc",var.curr), par.likeli["GLOB_Mult_T_a_lik"], 1)
	}
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]*par.likeli["GLOB_Mult_a_lik"]*mult_a
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]*par.likeli["GLOB_Mult_b_lik"]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]*par.likeli["GLOB_Mult_taumin_lik"]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]*par.likeli["GLOB_Mult_taumax_lik"]*ifelse(grepl("Wv",var.curr), par.likeli["GLOB_Mult_Q_taumax_lik"], 1)*ifelse(grepl("Tc",var.curr), par.likeli["GLOB_Mult_T_taumax_lik"], 1)
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        P     <- P.all[ind.var]
        n <- sum(ind.var)
        if ( any(df < Inf) | any(gamma != 1) | df<Inf | gamma!=1){## calculate sd of skewed t distribution
            if(any(df<Inf)){
                Exb <- 2*(gamma^2-(1/(gamma^2)))/(gamma+1/gamma)*gamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2))*sqrt(df*(df-2))/(df-1)
            }else{
                Exb <- 2*(gamma^2-(1/(gamma^2)))/(gamma+1/gamma)*1/sqrt(2*pi)
            }
            Ex <- (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            chunk <- sqrt(Ex - Exb^2)
            if(mode){
                Exb <- 0 ## ATTENTION: this means that Qdet is at the mode. If we remove these two lines, Qdet is at the mean.
            }
        }else{
            chunk <- 1
            Exb <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }

        Q.sim <- pmax(Q.sim , 0)
        sd <- a*Q01*(Q.sim/Q01)^c + Q01*b + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)

        log.p <- rep(NA,n-2)
        sdt.b <- sd[1:(n-1)]/chunk
        sdt.f <- sd[2:n]/chunk
        P <- pmax(P-Plim,0)
        dt <- (t.obs[2:n]-t.obs[1:(n-1)])
        if(any(dt<=0)) stop(paste("dt was zero or negative for variable", var.curr))
        taus <- tau.min + (tau.max-tau.min)*exp(-tkP*P^l - tkQ*Q.sim^m)
        gmms <- 1/taus
        gmms <- gmms[2:n]
        sdt.ini <- sd[1]/chunk
        log.p.ini <- c(0,0)
            if ( Q.obs[1] > 0 ){
                log.p.ini <- mydskt((Q.obs[1]-Q.sim[1]+Exb*sdt.ini),df=df[1],gamma=gamma[1],sigm=sdt.ini,log=TRUE)
            }else{
                log.p.ini <- mypskt(-Q.sim[1]+Exb*sdt.ini,df=df[1],gamma=gamma[1],sigm=sdt.ini,log.p=TRUE)
            }
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f  <- Q.obs[2:n]-Q.sim[2:n]+Exb*sdt.f
        dq.b  <- Q.obs[1:(n-1)]-Q.sim[1:(n-1)]+Exb*sdt.b
        if(any(!is.finite(c(dq.b[1],dq.f)))){warning("encountered non-finite dq.b/f");return(-Inf)}
        lower0 <- dq.b<0
        higher0<- dq.b>=0
        quant.b <- numeric(length(dq.b))
        quant.b[lower0]  <- qnorm(mypskt(dq.b[lower0],df=df,gamma=gamma,sigm=sdt.b[lower0],log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.b[higher0] <- qnorm(mypskt(dq.b[higher0],df=df,gamma=gamma,sigm=sdt.b[higher0],low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        ind.qob <- qob & qob2
        ind.qob2<- qob & !qob2
        ind.qob3<- !qob & qob2
        ind.qob4<- !qob & !qob2
        lower0 <- dq.f<0
        higher0<- dq.f>=0
        quant.f <- numeric(length(dq.f))
        quant.f[lower0] <- qnorm(mypskt(dq.f[lower0],df=df,gamma=gamma,sigm=sdt.f[lower0],log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.f[higher0] <- qnorm(mypskt(dq.f[higher0],df=df,gamma=gamma,sigm=sdt.f[higher0],low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)

        densEt <- dnorm(quant.f[ind.qob], mean=quant.b[ind.qob]*exp(-gmms[ind.qob]*dt[ind.qob]), sd=sqrt(1-exp(-2*gmms[ind.qob]*dt[ind.qob])), log=TRUE)

        c <- mydskt(dq.f[ind.qob],df=df,gamma=gamma,sigm=sdt.f[ind.qob],log=TRUE)
        d <- dnorm(quant.f[ind.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.qob] <- c + densEt - d
        log.p[ind.qob2] <- pnorm(quant.f[ind.qob2], mean=quant.b[ind.qob2]*exp(-gmms[ind.qob2]*dt[ind.qob2]),sd=sqrt(1-exp(-2*gmms[ind.qob2]*dt[ind.qob2])), log=TRUE)
        log.p[ind.qob3] <- mydskt(dq.f[ind.qob3],df=df,gamma=gamma,sigm=sdt.f[ind.qob3],log=TRUE)
        log.p[ind.qob4] <- pnorm(quant.f[ind.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        if(weight.equally) log.p <- log.p/length(log.p)
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
        ## white.noise1 <- (quant.f1 - quant.b1*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        ## white.noise2 <- (quant.f2 - quant.b2*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        ## innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])))
        ## innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])))
        ## unifvalue1 <- mypskt(dq.b1,df=df1,gamma=gamma1,sigm=sdt.b1)
        ## unifvalue2 <- mypskt(dq.b2,df=df2,gamma=gamma2,sigm=sdt.b2)
        ## unifvalue.f1 <- mypskt(dq.f1[length(dq.f1)],df=df1,gamma=gamma1[length(gamma1)],sigm=sdt.f1[length(sdt.f1)])
        ## unifvalue.f2 <- mypskt(dq.f2[length(dq.f2)],df=df2,gamma=gamma2,sigm=sdt.f2[length(sdt.f2)])
        ## return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), white.noise=c(ifelse(auto, white.noise2, white.noise1), ifelse(auto[length(auto)], white.noise2[length(white.noise2)], white.noise1[length(white.noise1)])), auto=c(FALSE,auto), a=densEta, b=densEtb, sd=sd1, taus=taus))
        return(list(smm=sum(loglikeli), loglik=loglikeli))
        warning("!verbose is not fully implemented yet for multiple variables...")
        return(NA)
    }
    ##print(warnings())
    return(sum(loglikeli))
}

LogLikelihoodHydrology_la9compound_fast_skewt <- function(run.model, layout, y.obs, P, par.likeli, auto, mode, apprx, verbose, weight.equally, ...){ ## simplified version of la9(e), where we first rescale, and then skew the distribution DQ.
    if(is.null(layout$calib)){
        L <- layout$layout
    }else{
        L             <- layout$layout[layout$calib,]
        layout$layout <- L
        layout$lump   <- layout$lump[layout$calib]
        y.obs         <- y.obs[layout$calib]
    }
    y.mod <- as.numeric(run.model(layout=layout, ...)$incld.lmpd)
    if(any(!is.numeric(y.mod))) stop("y.mod contains non-numeric values")
    loglikeli <- numeric(length=nrow(L))
    vars <- unique(L[,1])
    P.all <- P # copies, since they are overwritten later on
    auto.all <- auto
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        Q02 <- par.likeli[paste(var.curr, "_Q02_lik", sep = "")]
        ## Q02 <- Q01
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]*par.likeli["GLOB_Mult_a_lik"]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]*par.likeli["GLOB_Mult_a_lik"]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]*par.likeli["GLOB_Mult_b_lik"]
        b2  <- par.likeli[paste(var.curr, "_b2_lik", sep = "")]
        ## b2 <- b1
        c1  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        c2  <- par.likeli[paste(var.curr, "_c2_lik", sep = "")]
        ## c2 <- c1
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]*par.likeli["GLOB_Mult_taumin_lik"]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]*par.likeli["GLOB_Mult_taumax_lik"]*ifelse(grepl("Wv",var.curr), par.likeli["GLOB_Mult_Q_taumax_lik"], 1)*ifelse(grepl("Tc",var.curr), par.likeli["GLOB_Mult_T_taumax_lik"], 1)
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        df.out <- par.likeli[paste(var.curr, "_dfout_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        gamma.out <- par.likeli[paste(var.curr, "_gammaout_lik", sep = "")]
        pout      <- par.likeli[paste(var.curr, "_pout_lik", sep = "")]
        ##tau.max <- 20
        ##tau.min <- 0
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        Q.obs <- y.obs[ind.var]
        t.obs <- L[ind.var,2]
        auto  <- auto.all[ind.var]
        P     <- P.all[ind.var]
        n <- sum(ind.var)
        ##ma.p <- ma.p[ind.var][3:n]
        ##auto  <- P.var[3:n] <= d & P.var[2:(n-1)] <=d & P.var[1:(n-2)] <= d
        ## calculate actual value of likelihood
        ## P <- pmax(P,0)
        ## gamma1 <- pmin(gamma1+tkQ*P^m,5)
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1 | gamma.out!=1 | df.out<Inf){## calculate sd of skewed t distribution
            if(df1<Inf){
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
            }else{
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
            }else{
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            if(df.out<Inf){
                Exb.out <- 2*(gamma.out^2-(1/(gamma.out^2)))/(gamma.out+1/gamma.out)*gamma((df.out+1)/2)/(sqrt(df.out*pi)*gamma(df.out/2))*sqrt(df.out*(df.out-2))/(df.out-1)
            }else{
                Exb.out <- 2*(gamma.out^2-(1/(gamma.out^2)))/(gamma.out+1/gamma.out)*1/sqrt(2*pi)
            }
            chunk1 <- (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            chunk2 <- (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            chunk.out <- (gamma.out^3 + 1/(gamma.out^3)) / (gamma.out + 1/gamma.out)
            Ex.comp <- (1-pout)*Exb1 + pout*Exb.out
            Varx1 <- sqrt(chunk1 - Exb1^2)
            Varx2 <- sqrt(chunk2 - Exb2^2)
            Varx.out <- sqrt(chunk.out - Exb.out^2)
            Var.comp <- (1-pout)*(Varx1 + Exb1^2) + pout*(Varx.out + Exb.out^2) - Ex.comp^2
            if(mode){
                Ex.comp <- 0 ## ATTENTION: this means that Qdet is at the mode. If we remove this line, Qdet is at the mean.
            }
        }else{
            Var.comp <- 1
            Ex.comp <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }
        if ( n != length(Q.obs) ) { cat("*** length of Q.obs not equal to length of t\n"); return(NA) }
        if ( is.na(a1) | is.na(b1) | is.na(c1) )
        { cat("*** provide named parameters a, b, c and tau\n"); return(NA) }
        ## evaluate truncated normal-based log likelihood:
        Q.sim <- pmax(Q.sim , 0)
        ## dQsim <- pmax(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])),0)
        ## dQpos1 <- pmax(rollmean(dQsim, k=3, fill=0, align="left"), 0)

        sd1 <- a1*Q01*(Q.sim/Q01)^c1 + Q01*b1 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        sd2 <- a2*Q02*(Q.sim/Q02)^c2 + Q02*b2 + d*pmax(P-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)

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
        sdt.b1 <- sd1[1:(n-1)]/Var.comp
        sdt.b2 <- sd2[1:(n-1)]/Var.comp
        sdt.f1 <- sd1[2:n]/Var.comp
        sdt.f2 <- sd2[2:n]/Var.comp
        ## gamma1 <- gamma1[2:n]
        ## dQpos2 <- abs(c(0,diff(Q.sim)/(t.obs[2:n]-t.obs[1:(n-1)])/(0.5*Q.sim[1:(n-1)] + 0.5*Q.sim[2:n] + psi)))
        ## dQpos2 <- rollmean(dQpos2, k=3, fill=0, align="left")
        P <- pmax(P-Plim,0)
        dt <- (t.obs[2:n]-t.obs[1:(n-1)])
        if(any(dt<=0)) stop(paste("dt was zero or negative for variable", var.curr))
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
        sdt.ini <- ifelse(auto[1], sd2[1]/Var.comp, sd1[1]/Var.comp)
        log.p.ini <- c(0,0)
            df.tmp <- ifelse(auto[1], df2, df1)
            gamma.tmp <- ifelse(auto[1], gamma2[1], gamma1[1])
            if ( Q.obs[1] > 0 ){
                log.p.ini <- dcompound((Q.obs[1]-Q.sim[1]+Ex.comp*sdt.ini),df=df.tmp,df.out=df.out,gamma=gamma.tmp,gamma.out=gamma.out,sigm=sdt.ini,p=pout,log=TRUE)
            }else{
                log.p.ini <- pcompound(-Q.sim[1]+Ex.comp*sdt.ini,df=df.tmp,df.out=df.out,gamma=gamma.tmp,gamma.out=gamma.out,sigm=sdt.ini,p=pout,log.p=TRUE)
            }
        auto <- auto[2:n]
        qob   <- Q.obs[1:(n-1)] > 0
        qob2  <- Q.obs[2:n] > 0
        dq.f1  <- Q.obs[2:n]-Q.sim[2:n]+Ex.comp*sdt.f1
        dq.b1  <- Q.obs[1:(n-1)]-Q.sim[1:(n-1)]+Ex.comp*sdt.b1
        dq.f2  <- Q.obs[2:n]-Q.sim[2:n]+Ex.comp*sdt.f2
        dq.b2  <- Q.obs[1:(n-1)]-Q.sim[1:(n-1)]+Ex.comp*sdt.b2
        lower0 <- dq.b1<0
        higher0<- dq.b1>=0
        quant.b1 <- numeric(length(dq.b1))
        quant.b1[lower0]  <- qnorm(pcompound(dq.b1[lower0],df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt.b1[lower0],p=pout,log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.b1[higher0] <- qnorm(pcompound(dq.b1[higher0],df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt.b1[higher0],p=pout,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        lower0 <- dq.b2<0
        higher0<- dq.b2>=0
        quant.b2 <- numeric(length(dq.b2))
        quant.b2[lower0] <- qnorm(pcompound(dq.b2[lower0],df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt.b2[lower0],p=pout,log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.b2[higher0] <- qnorm(pcompound(dq.b2[higher0],df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt.b2[higher0],p=pout,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
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
        quant.f1[lower0] <- qnorm(pcompound(dq.f1[lower0],df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt.f1[lower0],p=pout,log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.f1[higher0] <- qnorm(pcompound(dq.f1[higher0],df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt.f1[higher0],p=pout,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
        lower0 <- dq.f2<0
        higher0<- dq.f2>=0
        quant.f2 <- numeric(length(dq.f2))
        quant.f2[lower0] <- qnorm(pcompound(dq.f2[lower0],df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt.f2[lower0],p=pout,log.p=TRUE),mean=0,sd=1,log.p=TRUE)
        quant.f2[higher0] <- qnorm(pcompound(dq.f2[higher0],df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt.f2[higher0],p=pout,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)

        densEta <- dnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])), log=TRUE)
        densEtb <- dnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])), log=TRUE)

        log.p[ind.auto.qob] <- dcompound(dq.f2[ind.auto.qob],df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt.f2[ind.auto.qob],p=pout,log=TRUE) + densEta - dnorm(quant.f2[ind.auto.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto.qob2] <- pnorm(quant.f2[ind.auto.qob2], mean=quant.b2[ind.auto.qob2]*exp(-gmms[ind.auto.qob2]*dt[ind.auto.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto.qob2]*dt[ind.auto.qob2])), log=TRUE)
        log.p[ind.auto.qob3] <- dcompound(dq.f2[ind.auto.qob3],df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt.f2[ind.auto.qob3],p=pout,log=TRUE)
        log.p[ind.auto.qob4] <- pnorm(quant.f2[ind.auto.qob4], mean=0,sd=1, log=TRUE)

        c <- dcompound(dq.f1[ind.auto2.qob],df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt.f1[ind.auto2.qob],p=pout,log=TRUE)
        d <- dnorm(quant.f1[ind.auto2.qob],mean=0,sd=1,log=TRUE)
        log.p[ind.auto2.qob] <- c + densEtb - d
        log.p[ind.auto2.qob2] <- pnorm(quant.f1[ind.auto2.qob2], mean=quant.b1[ind.auto2.qob2]*exp(-gmms[ind.auto2.qob2]*dt[ind.auto2.qob2]),sd=sqrt(1-exp(-2*gmms[ind.auto2.qob2]*dt[ind.auto2.qob2])), log=TRUE)
        log.p[ind.auto2.qob3] <- dcompound(dq.f1[ind.auto2.qob3],df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt.f1[ind.auto2.qob3],p=pout,log=TRUE)
        log.p[ind.auto2.qob4] <- pnorm(quant.f1[ind.auto2.qob4], mean=0,sd=1, log=TRUE)

        log.p <- c(log.p.ini, log.p)
        if(weight.equally) log.p <- log.p/length(log.p)
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
        ## white.noise1 <- (quant.f1 - quant.b1*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        ## white.noise2 <- (quant.f2 - quant.b2*exp(-gmms*dt))/sqrt(1-exp(-2*gmms*dt))
        ## innovation1 <- pnorm(quant.f1[ind.auto2.qob], mean=quant.b1[ind.auto2.qob]*exp(-gmms[ind.auto2.qob]*dt[ind.auto2.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto2.qob]*dt[ind.auto2.qob])))
        ## innovation2 <- pnorm(quant.f2[ind.auto.qob], mean=quant.b2[ind.auto.qob]*exp(-gmms[ind.auto.qob]*dt[ind.auto.qob]), sd=sqrt(1-exp(-2*gmms[ind.auto.qob]*dt[ind.auto.qob])))
        ## unifvalue1 <- pcompound(dq.b1,df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt.b1,p=pout)
        ## unifvalue2 <- pcompound(dq.b2,df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt.b2,p=pout)
        ## unifvalue.f1 <- pcompound(dq.f1[length(dq.f1)],df=df1,df.out=df.out,gamma=gamma1[length(gamma1)],gamma.out=gamma.out,sigm=sdt.f1[length(sdt.f1)],p=pout)
        ## unifvalue.f2 <- pcompound(dq.f2[length(dq.f2)],df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt.f2[length(sdt.f2)],p=pout)
        ## return(list(smm=sum(loglikeli), loglik=loglikeli, unifvalue=c(ifelse(auto, unifvalue2, unifvalue1), ifelse(auto[length(auto)], unifvalue.f2, unifvalue.f1)), quant=c(ifelse(auto, quant.b2, quant.b1), ifelse(auto[length(auto)], quant.f2[length(quant.f2)], quant.f1[length(quant.f1)])), innovation=c(NA, ifelse(auto,innovation2,innovation1)), white.noise=c(ifelse(auto, white.noise2, white.noise1), ifelse(auto[length(auto)], white.noise2[length(white.noise2)], white.noise1[length(white.noise1)])), auto=c(FALSE,auto), a=densEta, b=densEtb, sd=sd1, taus=taus))
        return(list(smm=sum(loglikeli), loglik=loglikeli))
        warning("!verbose is not fully implemented yet for multiple variables...")
    }
    return(sum(loglikeli))
}

## =======================================================================================================
## Likelihood Samplers
## =======================================================================================================
## These are functions that sample from the likelihood (e.g. for prediction)

sampling_wrapper <- function(sudriv, brn.in=0, sample.par=TRUE, n.sample=1, sample.likeli=TRUE, auto=NA, mode=TRUE, eta=NA){
    ## sample from a population (sample) of parameters
    if(all(is.na(auto))){
        auto <- rep(FALSE, nrow(su$layout$pred.layout))
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
    likeli.sample <- matrix(nrow=n.sample, ncol=nrow(sudriv$layout$pred.layout))
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
            likeli.args           <- list()
            likeli.args$run.model <- run.model
            likeli.args$layout    <- sudriv$layout
            likeli.args$P         <- sudriv$input$P.roll.pred##[sudriv$layout$pred]
            likeli.args$par.likeli<- ifelse(as.logical(sudriv$likelihood$tran), exp(sudriv$likelihood$parameters), sudriv$likelihood$parameters)
            names(likeli.args$par.likeli) <- names(sudriv$likelihood$parameters)
            likeli.args$auto <- auto
            likeli.args$mode <- mode
            likeli.args$sudriv    <- sudriv
            likeli.args$lump <- FALSE
            f.sample <- sudriv$likelihood$f.sample
            ## =======================================================
            ## sample from the likelihood
            likeli.sample[i,] <- do.call(f.sample, likeli.args)
        }else{## in this case, we just run the deterministic model (propagate parameter uncertainty only)
            par <- sudriv$model$parameters
            L <- sudriv$layout
            L$layout <- L$pred.layout
            likeli.sample[i,] <- as.numeric(run.model(layout=L, sudriv=sudriv, lump=FALSE)$original)
        }
        cat(i," / ", n.sample, "\n")
    }
    return(likeli.sample)
}
sampling_wrapper_timedep <- function(sudriv, brn.in=0, sample.par=TRUE, n.sample=1, sample.likeli=TRUE, auto=NA, mode=TRUE, eta=NA, scaleshift=NA, mnprm=NA){
  ## sample from a population (sample) of parameters for the superfelx driver with the timedependent parameters
  if(all(is.na(auto))){
    auto <- rep(FALSE, nrow(su$layout$pred.layout))
  }
  if(sample.par){
    if(is.null(sudriv$parameter.sample.const)){
      warning("Sudriv object does not contain parameter sample to draw from. Drawing from the prior ...")
      ## Develop: draw parameter samples from prior distribution...
    }else{ ## Draw parameter sample from existing sample (representing e.g. posterior)
      if(brn.in >= nrow(sudriv$parameter.sample.const)) stop("brn.in is longer than chain ...")
      s.cnst <- sudriv$parameter.sample.const[(brn.in+1):nrow(sudriv$parameter.sample.const),]
      for(i in names(sudriv$parameter.sample.timedep)){
        sudriv$parameter.sample.timedep[[i]][-1,] <- sudriv$parameter.sample.timedep[[i]][-1,][(brn.in+1):nrow(sudriv$parameter.sample.const),]
      }
      ind.chosen <- sample(x=1:nrow(s.cnst), size=n.sample, replace=TRUE)
    }
  }
  likeli.sample <- matrix(nrow=n.sample, ncol=nrow(sudriv$layout$pred.layout))
  rtruncnorm <- function (n, a = -Inf, b = Inf, mean = 0, sd = 1){
    if (length(n) > 1)
      n <- length(n)
    if (length(n) > 1)
      n <- length(n)
    else if (!is.numeric(n))
      stop("non-numeric argument n.")
    .Call("do_rtruncnorm", as.integer(n), a, b, mean, sd, PACKAGE="truncnorm")}
  
  any.timedep <- FALSE
  if(length(sudriv$parameter.sample.timedep)>0){
    any.timedep <- TRUE
    if(!all(is.na(scaleshift))){ ## back-transform parameter with sigmoid transformation
      if(nrow(scaleshift)!=length(sudriv$parameter.sample.timedep) | ncol(scaleshift)!=2) stop("dimension of scaleshift is not right")
      for(i in 1:length(sudriv$parameter.sample.timedep)){
        sudriv$parameter.sample.timedep[[i]] <- sigm.trans(sudriv$parameter.sample.timedep[[i]], scale=scaleshift[i,1], shift=scaleshift[i,2])
      }
    }
    if(!all(is.na(mnprm))){ ## shift mean of some timedep-parameters for which the mean was re-parameterized as a constant parameter
      td <- names(sudriv$parameter.sample.timedep)
      ind.mnprm <- which(td %in% mnprm)
      if(!all(mnprm %in% td)) stop("some parameters of ", mnprm, " not found")
      #if(sum(ind.mnprm)!=1) stop(paste0("too many or too few parameters of ", mnprm, " found"))
      for(i in 1:length(ind.mnprm)){
        mn <- sudriv$parameter.sample.const[,paste0(mnprm[i],"_fmean")]
        sudriv$parameter.sample.timedep[[ind.mnprm[i]]][-1,] <- sudriv$parameter.sample.timedep[[ind.mnprm[i]]][-1,] * mn
      }
    }
  }
  
  flp <- sudriv$likelihood$par.fit
  fmp <- sudriv$model$par.fit
  l.fit.lik <- sum(flp)
  l.fit <- sum(c(fmp,flp))
  ind.timedep <- names(c(sudriv$model$parameters[as.logical(fmp)],sudriv$likelihood$parameters[as.logical(flp)])) %in% names(sudriv$parameter.sample.timedep)
  for(i in 1:n.sample){
    if(sample.par){
      x0 <- rep(NA, l.fit)
      ## create the vector of time-constant parameters
      x0.cnst <- s.cnst[ind.chosen[i],]
      x0[!ind.timedep] <- x0.cnst
      ## create matrix of time dependent parameters
      x0.tmdp <- do.call(cbind, lapply(sudriv$parameter.sample.timedep, function(x,ind.chosen,i) x[ind.chosen[i],], ind.chosen=ind.chosen, i=i))
      sudriv$model$timedep$par <- x0.tmdp
      if(any.timedep) x0[ind.timedep]  <- as.numeric(x0.tmdp[1,]) ## the value of x0 for the time-dependent parameter should not matter, since it is taken from su$model$timedep$par

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
      likeli.args           <- list()
      likeli.args$run.model <- run.model
      likeli.args$layout    <- sudriv$layout
      likeli.args$P         <- sudriv$input$P.roll.pred##[sudriv$layout$pred]
      likeli.args$par.likeli<- ifelse(as.logical(sudriv$likelihood$tran), exp(sudriv$likelihood$parameters), sudriv$likelihood$parameters)
      names(likeli.args$par.likeli) <- names(sudriv$likelihood$parameters)
      likeli.args$auto <- auto
      likeli.args$mode <- mode
      likeli.args$sudriv    <- sudriv
      likeli.args$lump <- FALSE
      f.sample <- sudriv$likelihood$f.sample
      ## =======================================================
      ## sample from the likelihood
      likeli.sample[i,] <- do.call(f.sample, likeli.args)
    }else{## in this case, we just run the deterministic model (propagate parameter uncertainty only)
      par <- sudriv$model$parameters
      L <- sudriv$layout
      L$layout <- L$pred.layout
      likeli.sample[i,] <- as.numeric(run.model(layout=L, sudriv=sudriv, lump=FALSE)$original)
    }
    cat(i," / ", n.sample, "\n")
  }
  return(likeli.sample)
}
LogLikelihoodHydrology_la9esimp_skewt_sample <- function(run.model, P, layout, par.likeli, auto, mode, ...){
    options(warn=2)
    layout$layout <- layout$pred.layout
    L <- layout$layout
    y.mod <- as.numeric(run.model(layout=layout, ...)$original)
    if(any(is.na(y.mod))) stop("y.mod contains nas")
    vars <- unique(L[,1])
    samp <- numeric(length=nrow(L))
    for(var.curr in vars){
        Q01 <- par.likeli[paste(var.curr, "_Q01_lik", sep = "")]
        a  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]*par.likeli["GLOB_Mult_a_lik"]
        b  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]*par.likeli["GLOB_Mult_b_lik"]
        c  <- par.likeli[paste(var.curr, "_c_lik", sep = "")]
        d  <- par.likeli[paste(var.curr, "_d_lik", sep = "")]
        e  <- par.likeli[paste(var.curr, "_e_lik", sep = "")]
        Plim  <- par.likeli[paste(var.curr, "_Plim_lik", sep = "")]
        ttP<- par.likeli[paste(var.curr, "_ttP_lik", sep = "")]
        tkP<- par.likeli[paste(var.curr, "_tkP_lik", sep = "")]
        tkQ<- par.likeli[paste(var.curr, "_tkQ_lik", sep = "")]
        l  <- par.likeli[paste(var.curr, "_l_lik", sep = "")]
        m  <- par.likeli[paste(var.curr, "_m_lik", sep = "")]
        psi  <- par.likeli[paste(var.curr, "_psi_lik", sep = "")]
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]*par.likeli["GLOB_Mult_taumin_lik"]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]*par.likeli["GLOB_Mult_taumax_lik"]*ifelse(grepl("Wv",var.curr), par.likeli["GLOB_Mult_Q_taumax_lik"], 1)*ifelse(grepl("Tc",var.curr), par.likeli["GLOB_Mult_T_taumax_lik"], 1)
        df <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( any(df < Inf) | any(gamma != 1) | df<Inf | gamma!=1){## calculate sd of skewed t distribution
            if(any(df<Inf)){
                Exb <- 2*(gamma^2-(1/(gamma^2)))/(gamma+1/gamma)*gamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2))*sqrt(df*(df-2))/(df-1)
            }else{
                Exb <- 2*(gamma^2-(1/(gamma^2)))/(gamma+1/gamma)*1/sqrt(2*pi)
            }
            Ex <- (gamma^3 + 1/(gamma^3)) / (gamma + 1/gamma)
            chunk <- sqrt(Ex - Exb^2)
            if(mode){
                Exb <- 0 ## ATTENTION: this means that Qdet is at the mode. If we remove these two lines, Qdet is at the mean.
            }
        }else{
            chunk <- 1
            Exb <- 0
        }
        ## consistency checks:
        if ( n != length(Q.sim) ) { cat("*** length of Q.sim not equal to length of t\n"); return(NA) }

        Q.sim <- pmax(Q.sim , 0)
        sd <- a*Q01*(Q.sim/Q01)^c + Q01*b + d*pmax(P.var-Plim,0)^e##rollmean(c(0,pmax(diff(Q.sim),0)),k=2,fill=0)
        samp.var <- Q.sim
        sdt <- ifelse(auto, sd2/chunk2, sd/chunk)
        P.var <- pmax(P.var-Plim,0)
        taus <- tau.min + (tau.max-tau.min)*exp(-tkP*P.var^l - tkQ*Q.sim^m)
        gmms <- 1/taus
        dt <- (t.mod[2:n]-t.mod[1:(n-1)])
        gmms <- gmms[2:n]
        quant.1 <- rnorm(1)
        samp.var[1] <- pmax(samp.var[1] + myqskt(pnorm(quant.1), df=df,gamma=gamma,sigm=sdt[1])-sdt[1]*Exb[1], 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                et.old <- samp.var[i-1]-Q.sim[i-1]+Exb*sdt[i-1]
                if(et.old>=0){
                    quant.1 <- qnorm(mypskt(et.old,df=df,gamma=gamma,sigm=sdt[i-1],low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                }else{
                    quant.1 <- qnorm(mypskt(et.old,df=df,gamma=gamma,sigm=sdt[i-1],log.p=TRUE),mean=0,sd=1,log.p=TRUE)
                }
                quant.2 <- rnorm(n=1, mean=quant.1*exp(-gmms[i-1]*dt[i-1]), sd=sqrt(1-exp(-2*gmms[i-1]*dt[i-1])))
                samp.var[i] <- max(samp.var[i] + myqskt(pnorm(quant.2), df=df,gamma=gamma,sigm=sdt[i])-sdt[i]*Exb, 0)
                if(is.infinite(samp.var[i])){print(quant.1);print(quant.2);stop("infinites in prediction")}
            }else{
                samp.var[i] <- max(samp.var[i] + myqskt(pnorm(rnorm(1)),df=df,gamma=gamma,sigm=sdt[i])-sdt[i]*Exb, 0)
            }
        samp[ind.var] <- samp.var
        }
    }
    options(warn=0)
    return(samp)
}

LogLikelihoodHydrology_la9compound_sample <- function(par.model, run.model, P, layout, par.likeli, mu, rep.mu.times, time.recess, auto, tau.Qdet, ...){
    options(warn=2)
    layout$layout <- layout$layout[layout$pred,]
    L <- layout$layout
    P <- P[layout$pred]
    y.mod <- as.numeric(run.model(layout=layout, ...)$incld.lmpd)
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
        a1  <- par.likeli[paste(var.curr, "_a_lik", sep = "")]*par.likeli["GLOB_Mult_a_lik"]
        a2  <- par.likeli[paste(var.curr, "_a2_lik", sep = "")]*par.likeli["GLOB_Mult_a_lik"]
        ## a2 <- a1
        b1  <- par.likeli[paste(var.curr, "_b_lik", sep = "")]*par.likeli["GLOB_Mult_b_lik"]
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
        tau.min <- par.likeli[paste(var.curr, "_taumin_lik", sep = "")]*par.likeli["GLOB_Mult_taumin_lik"]
        tau.max <- par.likeli[paste(var.curr, "_taumax_lik", sep = "")]*par.likeli["GLOB_Mult_taumax_lik"]*ifelse(grepl("Wv",var.curr), par.likeli["GLOB_Mult_Q_taumax_lik"], 1)*ifelse(grepl("Tc",var.curr), par.likeli["GLOB_Mult_T_taumax_lik"], 1)
        ## tau2 <- 1.48*tau1
        df1 <- par.likeli[paste(var.curr, "_df_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        ## df2 <- par.likeli[paste(var.curr, "_df2_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        df2 <- df1
        df.out <- par.likeli[paste(var.curr, "_dfout_lik", sep = "")] + 2 ## ATTENTION: note the +2 here...
        gamma1 <- par.likeli[paste(var.curr, "_gamma_lik", sep = "")]
        ## gamma2 <- par.likeli[paste(var.curr, "_gamma2_lik", sep = "")]
        gamma2 <- gamma1
        gamma.out <- par.likeli[paste(var.curr, "_gammaout_lik", sep = "")]
        pout      <- par.likeli[paste(var.curr, "_pout_lik", sep = "")]
        ind.var <- L[,1]==var.curr
        Q.sim <- y.mod[ind.var]
        P.var <- P[ind.var]
        t.mod <- L[ind.var,2]
        n <- sum(ind.var)
        ## calculate actual value of likelihood
        if ( df1 < Inf | gamma1 != 1 | df2<Inf | gamma2!=1 | gamma.out!=1 | df.out<Inf){## calculate sd of skewed t distribution
            if(df1<Inf){
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*gamma((df1+1)/2)/(sqrt(df1*pi)*gamma(df1/2))*sqrt(df1*(df1-2))/(df1-1)
            }else{
                Exb1 <- 2*(gamma1^2-(1/(gamma1^2)))/(gamma1+1/gamma1)*1/sqrt(2*pi)
            }
            if(df2<Inf){
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*gamma((df2+1)/2)/(sqrt(df2*pi)*gamma(df2/2))*sqrt(df2*(df2-2))/(df2-1)
            }else{
                Exb2 <- 2*(gamma2^2-(1/(gamma2^2)))/(gamma2+1/gamma2)*1/sqrt(2*pi)
            }
            if(df.out<Inf){
                Exb.out <- 2*(gamma.out^2-(1/(gamma.out^2)))/(gamma.out+1/gamma.out)*gamma((df.out+1)/2)/(sqrt(df.out*pi)*gamma(df.out/2))*sqrt(df.out*(df.out-2))/(df.out-1)
            }else{
                Exb.out <- 2*(gamma.out^2-(1/(gamma.out^2)))/(gamma.out+1/gamma.out)*1/sqrt(2*pi)
            }
            chunk1 <- (gamma1^3 + 1/(gamma1^3)) / (gamma1 + 1/gamma1)
            chunk2 <- (gamma2^3 + 1/(gamma2^3)) / (gamma2 + 1/gamma2)
            chunk.out <- (gamma.out^3 + 1/(gamma.out^3)) / (gamma.out + 1/gamma.out)
            Ex.comp <- (1-pout)*Exb1 + pout*Exb.out
            Varx1 <- sqrt(chunk1 - Exb1^2)
            Varx2 <- sqrt(chunk2 - Exb2^2)
            Varx.out <- sqrt(chunk.out - Exb.out^2)
            Var.comp <- (1-pout)*(Varx1 + Exb1^2) + pout*(Varx.out + Exb.out^2) - Ex.comp^2
            Ex.comp <- 0 ## ATTENTION: this means that Qdet is at the mode. If we remove this line, Qdet is at the mean.
        }else{
            Var.comp <- 1
            Ex.comp <- 0
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
        sdt <- ifelse(auto, sd2/Var.comp, sd1/Var.comp)
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
        samp.var[1] <- pmax(samp.var[1] + qcompound(pnorm(quant.1,lower.tail=quant.1<0), df=df.tmp,df.out=df.out,gamma=gamma.tmp,gamma.out=gamma.out,sigm=sdt[1],p=pout, low=quant.1<0)-sdt[1]*Ex.comp, 0)
        for ( i in 2:n ){
            if(samp.var[i-1]>0){
                if(!auto[i]){
                    et.old <- samp.var[i-1]-Q.sim[i-1]+Ex.comp*sdt[i-1]
                    if(et.old>=0){
                        quant.1 <- qnorm(pcompound(et.old,df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt[i-1],p=pout,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    }else{
                        quant.1 <- qnorm(pcompound(et.old,df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt[i-1],p=pout,log.p=TRUE),mean=0,sd=1,log.p=TRUE)
                    }
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-gmms[i-1]*dt[i-1]), sd=sqrt(1-exp(-2*gmms[i-1]*dt[i-1])))
                    samp.var[i] <- max(samp.var[i] + qcompound(pnorm(quant.2,lower.tail=quant.2<0), df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt[i],p=pout,low=quant.2<0)-sdt[i]*Ex.comp, 0)
                }else{
                    et.old <- samp.var[i-1]-Q.sim[i-1]+Ex.comp*sdt[i-1]
                    if(et.old>=0){
                        quant.1 <- qnorm(pcompound(et.old,df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt[i-1],p=pout,low=FALSE,log.p=TRUE),mean=0,sd=1,low=FALSE,log.p=TRUE)
                    }else{
                        quant.1 <- qnorm(pcompound(et.old,df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt[i-1],p=pout,log.p=TRUE),mean=0,sd=1,log.p=TRUE)
                    }
                    quant.2 <- rnorm(n=1, mean=quant.1*exp(-gmms[i-i]*dt[i-1]), sd=sqrt(1-exp(-2*gmms[i-1]*dt[i-1])))
                    samp.var[i] <- max(samp.var[i] + qcompound(pnorm(quant.2,lower.tail=quant.2<0), df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt[i],p=pout,low=quant.2<0)-sdt[i]*Ex.comp, 0)
                }
                if(is.infinite(samp.var[i])){print(quant.1);print(quant.2);print(auto[i]);print(sd2[i]);stop("infinites in prediction")}
            }else{
                if(!auto[i]){
                    z <- rnorm(1)
                    samp.var[i] <- max(samp.var[i] + qcompound(pnorm(z,lower.tail=z<0),df=df1,df.out=df.out,gamma=gamma1,gamma.out=gamma.out,sigm=sdt[i],p=pout,low=z<0)-sdt[i]*Ex.comp, 0)
                }else{
                    z <- rnorm(1)
                    samp.var[i] <- max(samp.var[i] + qcompound(pnorm(z,lower.tail=z<0),df=df2,df.out=df.out,gamma=gamma2,gamma.out=gamma.out,sigm=sdt[i],p=pout,low=z<0)-sdt[i]*Ex.comp, 0)
                }
            }
        samp[ind.var] <- samp.var
        }
    }
    options(warn=0)
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

dcompound <- function(x, df=Inf, df.out=5, gamma=1, gamma.out=1, sigm=1, p=1e-3, log=FALSE){
    result <- (1-p)*mydskt(x, df=df, gamma=gamma, sigm=sigm, log=FALSE) + p*mydskt(x, df=df.out, gamma=gamma.out, sigm=sigm, log=FALSE)
    if(log) result <- log(result)
    return(result)
}
pcompound <- function(x, df=Inf, df.out=5, gamma=1, gamma.out=1, sigm=1, p=1e-3, low=TRUE, log.p=FALSE){
    result <- (1-p)*mypskt(x, df=df, gamma=gamma, sigm=sigm, low=low, log.p=FALSE) + p*mypskt(x, df=df.out, gamma=gamma.out, sigm=sigm, low=low, log.p=FALSE)
    if(log.p) result <- log(result)
    return(result)
}
rcompound <- function(n, df=Inf, df.out=5, gamma=1, gamma.out=1, sigm=1, p=1e-3){
    k <- runif(n)
    result <- (1-p)*qcompound(k, df=df, df.out=df.out, gamma=gamma, gamma.out=gamma.out, sigm=sigm, p=p) + p*qcompound(k, df=df, df.out=df.out, gamma=1, gamma.out=gamma.out, sigm=sigm, p=p)
    return(result)
}
qcompound <- function(prob, df=Inf, df.out=5, gamma=1, gamma.out=1, sigm=1, p=1e-3, low=TRUE, log.p=FALSE){
    ##print(prob)
    ##print(low)
    lb <- myqskt(ifelse(low, prob, 1-prob), df=df.out, gamma=gamma.out, sigm=sigm)
    ub <- myqskt(ifelse(low, prob, 1-prob), df=df, gamma=gamma, sigm=sigm)
    ##print(sort(c(lb,ub)))
    minim <- function(x, prob, df, df.out, gamma, gamma.out, sigm, p, low, log.p){prob - pcompound(x, df=df, df.out=df.out, gamma=gamma, gamma.out=gamma.out, sigm=sigm, p=p, low=low, log.p=log.p)}
    result <- uniroot(f=minim, interval=sort(c(lb,ub)), extendInt=ifelse(low, "downX", "upX"), prob=prob, df=df, df.out=df.out, gamma=gamma, gamma.out=gamma.out, sigm=sigm, p=p, low=low, log.p=log.p)$root
    ##(1-p)*myqskt(prob, df=df, gamma=gamma,sigm=sigm) + p*myqskt(prob, df=df.out, gamma=gamma.out,sigm=sigm)
    return(result)
}
