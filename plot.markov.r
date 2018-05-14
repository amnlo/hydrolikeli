s.m.mcmc.wrapper <- function(log.posterior, max.iter, sudriv, ...){
    done.iter <- 0
    iter.curr <- 5000
    while(done.iter < max.iter){
        if(done.iter != 0) init.range <- redef.init.range(sudriv)
        result.s.m = s.m.mcmc(log.posterior, max.iter=iter.curr, sudriv=sudriv, ...)
        s <- result.s.m$samples
        sudriv$parameter.sample <- aperm(s, perm=c(2,3,1))
        colnames(sudriv$parameter.sample) <- c(names(sudriv$model$parameters)[as.logical(sudriv$model$par.fit)], names(sudriv$likelihood$parameters)[as.logical(sudriv$likelihood$par.fit)])##, paste("mu", 1:(1*n.mu), sep=""))
        sudriv$posterior.sample <- t(result.s.m$log.p)
        done.iter <- done.iter + iter.curr
    }
    return(sudriv)
}
redef.init.range <- function(sudriv){
    logpost <- quantile(sudriv$posterior.sample[nrow(sudriv$posterior.sample),], 0.9)
    s <- remove.chains(sudriv, brn.in=nrow(sudriv$posterior.sample)-2, logpost=logpost)$sample
    init.range <- t(apply(s[nrow(s),,], 1, quantile, probs=c(0.1,0.9)))
    return(init.range)
}
adapt.prior <- function(sudriv){##adapt prior based on the timestep we are using
    fac  <- ifelse(sudriv$layout$time.units=="days", sudriv$layout$timestep.fac*24, sudriv$layout$timestep.fac)
    adpt <- grepl("K_Q", names(sudriv$model$parameters)) & sudriv$model$par.fit
    for(i in which(adpt)){
        len.parpri <- length(sudriv$model$prior$distdef[[i]])
        if(as.logical(sudriv$model$args$parTran[i])){
            sudriv$model$prior$distdef[[i]][(2:len.parpri)[-2]] <- as.character(as.numeric(sudriv$model$prior$distdef[[i]][(2:len.parpri)[-2]])+log(fac))
        }else{
            sudriv$model$prior$distdef[[i]][2:len.parpri] <- as.character(as.numeric(sudriv$model$prior$distdef[[i]][2:len.parpri])*fac)
        }
    }
    adpt <- grepl("tau", names(sudriv$likelihood$parameters)) & sudriv$likelihood$par.fit
    for(i in which(adpt)){
        len.parpri <- length(sudriv$likelihood$prior$distdef[[i]])
        if(as.logical(sudriv$likelihood$tran[i])){
            sudriv$likelihood$prior$distdef[[i]][(2:len.parpri)[-2]] <- as.character(as.numeric(sudriv$likelihood$prior$distdef[[i]][(2:len.parpri)[-2]])-log(fac))
        }else{
            sudriv$likelihood$prior$distdef[[i]][2:len.parpri] <- as.character(as.numeric(sudriv$likelihood$prior$distdef[[i]][2:len.parpri])/fac)
        }
    }
    return(sudriv)
}
select.maxlikpars <- function(sudriv){
    ## select parameters with highest posterior probability
    ndim <- length(dim(sudriv$parameter.sample))
    if(ndim==3){
        m <- which(sudriv$posterior.sample == max(sudriv$posterior.sample), arr.ind=TRUE)[1,]
        m.par <- sudriv$parameter.sample[m[1],1:sum(c(sudriv$model$par.fit, sudriv$likelihood$par.fit)),m[2]]
    }
    if(ndim==2){
        m <- which(sudriv$posterior.sample == max(sudriv$posterior.sample))[1]
        m.par <- sudriv$parameter.sample[m,1:sum(c(sudriv$model$par.fit, sudriv$likelihood$par.fit))]
    }
    sudriv$model$parameters[as.logical(sudriv$model$par.fit)] <- m.par[1:sum(sudriv$model$par.fit)]
    sudriv$likelihood$parameters[as.logical(sudriv$likelihood$par.fit)] <- m.par[(sum(sudriv$model$par.fit)+1):length(m.par)]
    return(sudriv)
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
interp.precip <- function(time.p, precip, layout, avr.prev=TRUE){## ATTENTION: behaviour is unclear when precipitation contains missings.
    if(avr.prev){
        p.intp <- approx(x=time.p, y=precip, xout=layout$time)$y
    }else{
        t.tmp1 <- time.p - 0.01*abs(diff(time.p[1:2]))
        t.tmp2 <- time.p + 0.01*abs(diff(time.p[1:2]))
        times.p.double <- c(time.p[1]+diff(time.p[1:2]), c(t(cbind(t.tmp1,t.tmp2))))
        p.p.double     <- c(rep(precip, each=2), precip[length(precip)])
        p.intp <- approx(x=times.p.double, y=p.p.double, xout=layout$time)$y
    }
    return(p.intp)
}

plot.markov.hist <- function(sudriv, brn.in = 0, pridef = NULL, v.line=NULL, lower.logpost=NA){
    ## Visualizes marginal parameter distributions of Markov Chains
    par.names <- c(names(sudriv$model$parameters)[as.logical(sudriv$model$par.fit)], names(sudriv$likelihood$parameters)[as.logical(sudriv$likelihood$par.fit)])
    par.trans <- c(sudriv$model$args$parTran[as.logical(sudriv$model$par.fit)], sudriv$likelihood$tran[as.logical(sudriv$likelihood$par.fit)])

    ## remove burn-in
    ndim <- length(dim(sudriv$parameter.sample))
    if(ndim==3){s <- sudriv$parameter.sample[,par.names,]; w.names <- paste("w", 1:(dim(s)[ndim]), sep="")}
    if(ndim==2) s <- sudriv$parameter.sample[,par.names]
    post <- sudriv$posterior.sample
    if(ndim==3){
        s.brn <- s[(brn.in+1):(dim(s)[1]),,]
        post.brn <- post[(brn.in+1):nrow(post),]
        if(!is.na(lower.logpost)){
            rm.ind <- unique(which(post.brn < lower.logpost, arr.ind=TRUE)[,2])
            if(length(rm.ind)>0){
                post.brn <- post.brn[,-rm.ind]
                s.brn <- s.brn[,,-rm.ind]
                cat("removed ", length(rm.ind), " walkers\n")
            }else{warning("no chains were removed...")}
        }
        post.brn <- t(post.brn)
    }
    if(ndim==2){
        s.brn <- s[(brn.in+1):(dim(s)[1]),]
        post.brn <- post[(brn.in+1):length(post)]
    }
    ## create data frame for ggplot-object
    a <- s.brn
    value <- array(a, dim=c(prod(dim(a)), 1))
##    x    <- rep(1:(dim(a)[1]), times=prod(dim(a)[2:length(dim(a))]))
    param    <- rep(par.names, each = dim(a)[1], times = ifelse(ndim==3, dim(a)[3], 1))
    walker    <- ifelse(ndim==3, rep(w.names, each = dim(a)[1]*dim(a)[2]), NA)
    a.re <- data.frame(value=value, param=param, walker=walker, y=NA, pri=0)
    ## plot only values between certain quantiles
    ##if(length(probs) != 2){warning("only length 2 probs allowed"); return(NA)}
    ## quants <- tapply(a.re$value, a.re$param, quantile, probs=probs)
    ## for(pcurr in names(quants)){
    ##     a.re <- subset(a.re, param!=pcurr | (value > quants[[pcurr]][1] & value < quants[[pcurr]][2]))
    ## }
    ## back-transform parameters to original scale
    ind.trans <- a.re$param %in% c(names(sudriv$model$parameters)[as.logical(sudriv$model$args$parTran)], names(sudriv$likelihood$parameters)[as.logical(sudriv$likelihood$tran)])
    a.re[ind.trans,"value"] <- exp(a.re[ind.trans,"value"])

    ## insert prior into data frame to be plotted:
    if(!is.null(pridef)){
        l.pri <- 1000
        a.pri <- data.frame(value=rep(NA, l.pri*length(par.names)), param=NA, walker=NA, y=NA, pri=1)
        j <- 1
        for(par.curr in par.names){
            if(which(par.curr==par.names) %in% which(as.logical(par.trans)) & pridef[[par.curr]][1] %in% c("normal", "Normal", "norm", "Norm", "normaltrunc")){
                pridef[[par.curr]][1] <- "lognormal"
                m <- as.numeric(pridef[[par.curr]][2])
                s <- as.numeric(pridef[[par.curr]][3])
                pridef[[par.curr]][2] <- exp(m + s^2/2)
                pridef[[par.curr]][3] <- as.numeric(pridef[[par.curr]][2])*sqrt(exp(s^2)-1)
            }
            ##g.obj <- g.obj + stat_function(data=subset(a.re, param==par.curr), fun = calcpdf, args=list(distpar=pridef[[par.curr]], log=FALSE))
            mu <- as.numeric(pridef[[par.curr]][2])
            sd <- as.numeric(pridef[[par.curr]][3])
            if(is.na(sd)) sd <- mu
            pri.x <- seq(mu - 3*sd, mu + 3*sd, length.out=l.pri)
            rang  <- range(subset(a.re, param==par.curr)$value)
            if(rang[1]>pri.x[1] & rang[2]<pri.x[l.pri]){
                pri.x <- seq(rang[1], rang[2], length.out=l.pri)
            }
            pri.dens <- calcpdf(pri.x, distpar=pridef[[par.curr]], log=FALSE)
            a.pri[(l.pri*(j-1)+1):(l.pri*j),] <- data.frame(value=pri.x, param=par.curr, walker=NA, y=pri.dens, pri=1)
            j <- j + 1
        }
        a.re <- rbind(a.re, a.pri)
    }
    if(!is.null(v.line)){
        vl <- v.line[levels(as.factor(a.re$param))]
        vline.dat <- data.frame(param=levels(as.factor(a.re$param)), vl=vl)
    }

    ## actual plotting
    g.obj <- ggplot(mapping=aes(x=value)) + geom_density(data=subset(a.re, pri==0), fill="blue", alpha=0.1) + geom_line(mapping=aes(y=y), data=subset(a.re, pri==1)) + facet_wrap("param", nrow=floor(sqrt(dim(a)[2])), scales="free")
    ##geom_line(data=subset(a.re, pri==1)) +
    ## g.obj <- g.obj + scale_x_continuous(trans="log2")
    if(!is.null(v.line)){
        g.obj <- g.obj + geom_vline(aes(xintercept=vl), data=vline.dat)
    }
    plot(g.obj)
}

remove.chains <- function(sudriv, brn.in=0, logpost=NA, lower=TRUE){
    par.names <- c(names(sudriv$model$parameters)[as.logical(sudriv$model$par.fit)], names(sudriv$likelihood$parameters)[as.logical(sudriv$likelihood$par.fit)])
    s <- sudriv$parameter.sample[(brn.in+1):(dim(sudriv$parameter.sample)[1]),par.names,]
    post <- sudriv$posterior.sample[(brn.in+1):nrow(sudriv$posterior.sample),]
    if(!is.na(logpost)){
        if(lower){
            rm.ind <- unique(which(post < logpost, arr.ind=TRUE)[,2])
        }else{
            rm.ind <- unique(which(post > logpost, arr.ind=TRUE)[,2])
        }
        if(length(rm.ind)>0){
            post <- post[,-rm.ind]
            s <- s[,,-rm.ind]
        }else{warning("no chains were removed...")}
    }
    return(list(sample=s, post=post))
}

plot.markov.chain <- function(sudriv, brn.in = 0, thin=1, lower.logpost=NA){
    ## Visualizes marginal parameter distributions of Markov Chains
    ndim <- length(dim(sudriv$parameter.sample))
    par.names <- c(names(sudriv$model$parameters)[as.logical(sudriv$model$par.fit)], names(sudriv$likelihood$parameters)[as.logical(sudriv$likelihood$par.fit)])
    if(ndim==3){
        rm.chains <- remove.chains(sudriv, brn.in=brn.in, logpost=lower.logpost)
        s <- rm.chains$sample
        s <- s[(1:nrow(s))%%thin==0,,]
        post <- rm.chains$post
        post <- post[(1:nrow(post))%%thin==0,]
        w.names <- paste("w", 1:(dim(s)[3]), sep="")
    }
    if(ndim==2){
        s <- sudriv$parameter.sample[(brn.in+1):(dim(sudriv$parameter.sample)[1]),par.names]
        s <- s[(1:nrow(s))%%thin==0,]
        post <- sudriv$posterior.sample[(brn.in+1):length(sudriv$posterior.sample)]
        post <- post[(1:length(post))%%thin==0]
        w.names <- 1
    }
    ## create data frame for ggplot-object
                                        #a <- aperm(s, c(2,3,1))
    a <- s
    value <- array(a, dim=c(prod(dim(a)), 1))
    x    <- rep(1:(dim(a)[1]), times=prod(dim(a)[2:length(dim(a))]))
    param    <- rep(par.names, each = dim(a)[1], times = ifelse(ndim==3, dim(a)[3], 1))
    walker <- NA
    if(ndim==3) walker    <- rep(w.names, each = dim(a)[1]*dim(a)[2])
    a.re <- data.frame(x=x, value=value, param=param, walker=walker)
    post.chain <- data.frame(x=1:(dim(a)[1]), value=c(post), param="log.post", walker=rep(w.names, each=dim(a)[1]))
    a.re <- rbind(a.re, post.chain)
    ## back-transform parameters to original scale
    ind.trans <- a.re$param %in% c(names(sudriv$model$parameters)[as.logical(sudriv$model$args$parTran)], names(sudriv$likelihood$parameters)[as.logical(sudriv$likelihood$tran)])
    a.re[ind.trans,"value"] <- exp(a.re[ind.trans,"value"])


    ## actual plotting
    g.obj <- ggplot(data=a.re, mapping=aes(x=x,y=value, color=walker)) + geom_line() + facet_wrap("param", nrow=floor(sqrt(dim(a)[2])), scales="free") + theme(legend.position="none")
    plot(g.obj)
}

plot.cor <- function(sudriv, brn.in=0, thin=1, lower.logpost=NA){
    par.names <- c(names(sudriv$model$parameters)[as.logical(sudriv$model$par.fit)], names(sudriv$likelihood$parameters)[as.logical(sudriv$likelihood$par.fit)])
    ndim <- length(dim(sudriv$parameter.sample))
    if(ndim==3){
        s <- sudriv$parameter.sample[(brn.in+1):(dim(sudriv$parameter.sample)[1]),par.names,]
        s <- s[(1:nrow(s))%%thin==0,,]
        post <- sudriv$posterior.sample[(brn.in+1):nrow(sudriv$posterior.sample),]
        post <- post[(1:nrow(post))%%thin==0,]
        if(!is.na(lower.logpost)){
            rm.ind <- unique(which(post < lower.logpost, arr.ind=TRUE)[,2])
            if(length(rm.ind)>0){
                post <- post[,-rm.ind]
                s <- s[,,-rm.ind]
            }else{warning("no chains were removed...")}
        }
        df <- as.data.frame(matrix(c(aperm(s, perm=c(1,3,2))), prod(dim(s)[c(1,3)]), dim(s)[2]))
    }
    if(ndim==2){
        s <- sudriv$parameter.sample[(brn.in+1):(dim(sudriv$parameter.sample)[1]),par.names]
        s <- s[(1:nrow(s))%%thin==0,]
        df <- as.data.frame(s)
    }
    ## back-transform parameters to original scale
    ind.trans <- par.names %in% c(names(sudriv$model$parameters)[as.logical(sudriv$model$args$parTran)], names(sudriv$likelihood$parameters)[as.logical(sudriv$likelihood$tran)])
    for(i in which(ind.trans)){df[,i] <- exp(df[,i])}
    ## scale model parameters to the desired time format
    if(!is.null(sudriv$model$par.time)){
        for(i in par.names){
            if(i %in% names(sudriv$model$parameters)[sudriv$model$par.time==1]) df[,which(par.names==i)] <- df[,which(par.names==i)]*sudriv$layout$timestep.fac
            if(i %in% names(sudriv$model$parameters)[sudriv$model$par.time==-1]) df[,which(par.names==i)] <- df[,which(par.names==i)]/sudriv$layout$timestep.fac
        }
    }
    ## scale likelihood parameters to the desired time format
    if(!is.null(sudriv$likelihood$par.time)){
        for(i in par.names){
            if(i %in% names(sudriv$likelihood$parameters)[sudriv$likelihood$par.time==1]) df[,which(par.names==i)] <- df[,which(par.names==i)]*sudriv$layout$timestep.fac
            if(i %in% names(sudriv$likelihood$parameters)[sudriv$likelihood$par.time==-1]) df[,which(par.names==i)] <- df[,which(par.names==i)]/sudriv$layout$timestep.fac
        }
    }
    colnames(df) <- gsub("%", "", colnames(s))
    colnames(df) <- gsub("_lik", "", colnames(df))
    colnames(df) <- gsub("C1Wv_Qstream_", "", colnames(df))
    colnames(df) <- gsub("U1W", "", colnames(df))
    labels <- colnames(df)
    labels <- gsub("Cmlt_E", "C[E]", labels)
    labels <- gsub("Smax_UR", "S[max]", labels)
    labels <- gsub("K_Qb_UR", "k[u]", labels)
    labels <- gsub("K_Qq_FR", "k[f]", labels)
    labels <- gsub("taumax", "tau[max]", labels)
    labels <- substr(labels, start=nchar(labels)-10, stop=nchar(labels))
    myBreaks <- function(x){
        brks <- signif(as.numeric(quantile(x, probs=c(0.2,0.8))), 2)
        if(brks[1]==brks[2] | brks[1]<min(x) | brks[2]>quantile(x, 0.95)){
            brks <- signif(as.numeric(quantile(x, probs=c(0.2,0.8))), 3)
        }
        if(any(brks>=1000)){
            brks <- signif(as.numeric(quantile(x, probs=c(0.3,0.7))), 3)
        }
        return(brks)
    }
    myBreaks.diag <- function(x){
        brks <- signif(as.numeric(quantile(x, probs=c(0.2,0.5,0.8))), 2)
        if(brks[1]==brks[2] | brks[1]<min(x) | brks[2]>max(x)){
            brks <- signif(as.numeric(quantile(x, probs=c(0.2,0.5,0.8))), 3)
        }
        return(brks)
    }
    myfun <- function(data, mapping, ...){
        ggplot(data = data, mapping = mapping) + geom_point(size=0.2) + scale_y_continuous(breaks=myBreaks) + scale_x_continuous(breaks=myBreaks)
    }
    myfun.diag <- function(data, mapping, ...){
        ggplot(data = data, mapping = mapping)+ geom_density() + scale_x_continuous(breaks=myBreaks.diag)
    }
    ##    pm <- ggpairs(df, lower=list(continuous=wrap("points", size=0.2))) + theme_bw(base_size=20)
    pm <- ggpairs(df, lower=list(continuous=myfun), diag=list(continuous=myfun.diag), columnLabels=labels, labeller="label_parsed") + theme_bw(base_size=17)
    pm
}

plot.predictions <- function(list.su, probs=NA, n.samp=0, rand=TRUE, xlim=NA, ylim=NA, tme.orig="1000-01-01", lp.num.pred=NA,plt=TRUE){
    ## create data frame for ggplot-object
    sudriv <- list.su[[1]]
    sudriv$predicted$det <- sudriv$predicted$det/sudriv$layout$timestep.fac
    sudriv$predicted$sample <- sudriv$predicted$sample/sudriv$layout$timestep.fac
    sudriv$observations <- sudriv$observations/sudriv$layout$timestep.fac
    n.case <- length(list.su)
    ind.sel <- select.ind(list.su[[1]], xlim=xlim, ind.sel=NA)
    if(sum(ind.sel)==0){warning("no time period selected"); return(NA)}
    time <- sudriv$layout$layout$time[c(sudriv$layout$calib,sudriv$layout$pred)][ind.sel]
    ## time <- as.POSIXlt(x=tme.orig)+time*60*60
    time <- as.POSIXlt(x=tme.orig)+time*60*60*sudriv$layout$timestep.fac*ifelse(sudriv$layout$time.units=="days",24,1)
    if(n.samp > 0){## plot actual realisations
        if(n.case>1) stop("plotting realizations not implemented for multiple models")
        if(rand){
            sudriv$predicted$sample <- sudriv$predicted$sample[sample(1:nrow(sudriv$predicted$sample),n.samp),ind.sel,drop=FALSE]
        }else{
            sudriv$predicted$sample <- sudriv$predicted$sample[1:min(n.samp, nrow(sudriv$predicted$sample)),ind.sel,drop=FALSE]
        }
        dms <- dim(sudriv$predicted$sample)
        preds <- array(t(sudriv$predicted$sample), dim=c(prod(dms), 1))
        preds <- data.frame(x=rep(time, dms[1]), value=c(preds), simu=rep(paste("sim", 1:(dms[1]), sep=""), each = dms[2]))
        obs   <- data.frame(x=time, value=sudriv$observations[c(sudriv$layout$calib,sudriv$layout$pred)][ind.sel], simu="obs")
        dt <- sudriv$predicted$det[ind.sel]
        det <-   data.frame(x=time, value = c(dt))
        preds <- rbind(preds, obs)
        ## actual plotting
        g.obj <- ggplot(data=subset(preds, simu != "obs"), mapping=aes(x=x,y=value)) + geom_line(mapping=aes(color=simu), data=subset(preds, simu != "obs"), size=1.05) + geom_point(mapping=aes(color=simu), data=subset(preds, simu=="obs")) + geom_line(mapping=aes(color=simu), data=subset(preds, simu=="obs")) + geom_line(data=det, color="black", size=1.3)+ theme_bw(base_size=24)
    }else{## plot uncertainty bands
        if(!is.na(probs[1])){
            if(n.case>1) {warning("plotting uncertainty bands not implemented for multiple models. Set probs=NA in order not to plot uncertainty bands.");return()}
            ss <- sudriv$predicted$sample[,ind.sel]
            quants <- apply(ss, 2, quantile, probs=probs)
        }else{
            quants <- rbind(NA,NA)
        }
        dt <- sudriv$predicted$det[ind.sel]
        if(n.case>1){for(i in 2:n.case){dt <- c(dt,list.su[[i]]$predicted$det[ind.sel]/list.su[[i]]$layout$timestep.fac)}}
        pred <- data.frame(x=rep(time,n.case), value = as.numeric(dt), model=rep(names(list.su),each=length(time)), lower=c(quants[1,]), upper=c(quants[2,]))
        obs   <- data.frame(x=time, value=sudriv$observations[c(sudriv$layout$calib,sudriv$layout$pred)][ind.sel], model="observed")
        obs$lower <- obs$value
        obs$upper <- obs$value
        obspred <- rbind(pred, obs)
        g.obj <- ggplot(data=obspred, mapping=aes(x=x, y=value, linetype=model, color=model)) + geom_line(size=0.8) + labs(linetype="", color="", x="", y=ifelse(sudriv$layout$time.units=="hours", "Streamflow [mm/h]", "Streamflow [mm/d]")) + theme(axis.text=element_text(size=24)) + theme_bw(base_size=24) ##+ scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
        if(!is.na(probs[1])){
            outside <- obs$value > c(apply(sudriv$predicted$sample[,ind.sel], 2, quantile, probs=probs)[2,]) | obs$value < c(apply(sudriv$predicted$sample[,ind.sel], 2, quantile, probs=probs)[1,])
            frc <- round(1 - sum(outside)/length(outside), 2)
            mbe <- (sum(obs$value)-sum(dt))/sum(obs$value)*100
            nse <- 1-sum((dt-obs$value)^2)/sum((obs$value - mean(obs$value))^2)
            g.obj <- g.obj + geom_ribbon(mapping=aes(ymin=lower, ymax=upper), data=obspred, alpha=0.3, linetype=0) +labs(caption=paste("MBE: ", round(mbe), "%, NSE: ", round(nse,2), ", Logpost calib: ", round(lp.num.pred[1]), ", Frac. in bounds: ", frc, sep="")) ##+ annotate("text", x=Inf, y=Inf, label=paste("Fraction inside limits:", frc), hjust=1, vjust=1)
        }
    }
    if(!is.na(ylim[1])) g.obj <- g.obj + coord_cartesian(ylim=ylim)
    if(plt){
        plot(g.obj)
    }else{
        return(g.obj)
    }
}

pred.stats <- function(list.sudriv, auto=NA, time.recess=NA, mu=NA, rep.mu.times=NA, biased=FALSE, tme.orig="1985-01-01", brn.in=0, n.sample=500){
    ## create data frame for ggplot-object
    ## plot uncertainty bands
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
        det    <- c(sampling_wrapper(sudriv, sample.par=FALSE, n.sample=1, biased=biased, sample.likeli=FALSE, rep.mu.times=rep.mu.times, auto=auto, time.recess=time.recess, mu=mu))
        ## sample.calibpred <- sampling_wrapper(sudriv, brn.in=brn.in, sample.par=TRUE, n.sample=n.sample, biased=biased, sample.likeli=TRUE)
        ##sudriv$layout$calib <- sudriv$layout$pred
        x0 <- c(sudriv$model$parameters[as.logical(sudriv$model$par.fit)], sudriv$likelihood$parameters[as.logical(sudriv$likelihood$par.fit)])
        if(biased) x0 <- c(x0, mu)
        sudriv$layout$calib <- c(sudriv$layout$calib, sudriv$layout$pred)
        ll <-logposterior(x0=x0, sudriv=sudriv, prior=FALSE, verbose=FALSE, biased=biased, rep.mu.times=rep.mu.times)
        lp.num.pred <- logposterior(x0=x0, sudriv=sudriv, prior=TRUE, verbose=TRUE, biased=biased, rep.mu.times=rep.mu.times)
        dat[((i-1)*length(time)+1):(i*length(time)),] <- data.frame(x=time, det=det, quant=ll$quant, white.noise=ll$white.noise, unifvalue=ll$unifvalue, innovation=ll$innovation, obs=sudriv$observations[c(sudriv$layout$calib)], lp.num.pred=lp.num.pred, case=names(list.sudriv)[i], taus=ifelse(rep(is.null(ll$taus),length(time)), NA, ll$taus))
        ##obs   <- data.frame(x=time, value=sudriv$observations[sudriv$layout$pred])
        i <- i + 1
    }
    return(dat)
}
select.ind <- function(sudriv, xlim, ind.sel){
## create data frame for ggplot-object
    if(is.na(xlim[1])) xlim <- c(-Inf, Inf)
    if(xlim[1]=="pred") xlim <- c(sudriv$layout$layout$time[sudriv$layout$pred[1]],sudriv$layout$layout$time[sudriv$layout$pred[length(sudriv$layout$pred)]])
    if(xlim[1]=="calib") xlim <- c(sudriv$layout$layout$time[sudriv$layout$calib[1]],sudriv$layout$layout$time[sudriv$layout$calib[length(sudriv$layout$calib)]])
    if(is.na(ind.sel[1])){
        ind.sel <- which(sudriv$layout$layout$time[c(sudriv$layout$calib, sudriv$layout$pred)] >= xlim[1] & sudriv$layout$layout$time[c(sudriv$layout$calib, sudriv$layout$pred)] <= xlim[2])
    }else{
        ind.sel <- ind.sel[sudriv$layout$layout$time[c(sudriv$layout$calib, sudriv$layout$pred)][ind.sel] >= xlim[1] & sudriv$layout$layout$time[c(sudriv$layout$calib, sudriv$layout$pred)][ind.sel] <= xlim[2]]
    }
    return(ind.sel)
}
plot.ts.quantiles <- function(dat, sudriv, xlim=NA, ind.sel=NA, precip=FALSE, plim=0){
    ind.sel <- select.ind(sudriv=sudriv, xlim=xlim, ind.sel=ind.sel)
    n.case <- length(unique(dat$case))
    if(n.case > 1){
        for(i in 2:(n.case)){
            ind.sel <- c(ind.sel, ind.sel + (i-1)*sum(dat$case==dat$case[1]))
        }
    }
    ## mbe.calib <- (sum(dat$det[1:length(sudriv$layout$calib)])/sum(sudriv$observations[sudriv$layout$calib])-1)*100
    dat.ind.sel <- dat[ind.sel,]
    if(precip){
        dat.ind.sel[,"precip"] <- dat.ind.sel[,"precip"] > plim
        n.p <- which(colnames(dat)=="precip")
        dat.rect <- data.frame(from=dat.ind.sel[1,"x"], to=dat.ind.sel[2,"x"])
        for(i in 2:nrow(dat.ind.sel)){
            if(dat.ind.sel[i,n.p] & ! dat.ind.sel[i-1,n.p]) dat.rect <- rbind(dat.rect,data.frame(from=dat.ind.sel[i,"x"], to=dat.ind.sel[i,"x"]))
            if(!dat.ind.sel[i,n.p] &  dat.ind.sel[i-1,n.p]) dat.rect[nrow(dat.rect),2] <- dat.ind.sel[i-1,"x"]
        }
    }
    g.obj1 <- ggplot() + geom_point(mapping=aes(x=x,y=quant, shape=case), data=dat.ind.sel) + geom_line(mapping=aes(x=x,y=quant, linetype=case), data = dat.ind.sel)
    if(precip) g.obj1 <- g.obj1 + geom_rect(data=dat.rect, aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf), alpha=0.3)
    g.obj1 <- g.obj1 + labs(x="", y=expression(eta), shape="Error Model", linetype="Error Model")+ theme_bw(base_size=24) + theme(axis.text=element_text(size=24), legend.position="none") ##+ scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
    plot(g.obj1)
}
plot.ts.tau <- function(dat, sudriv, xlim=NA, ind.sel=NA){
    ind.sel <- select.ind(sudriv=sudriv, xlim=xlim, ind.sel=ind.sel)
    n.case <- length(unique(dat$case))
    if(n.case > 1){
        for(i in 2:(n.case)){
            ind.sel <- c(ind.sel, ind.sel + (i-1)*sum(dat$case==dat$case[1]))
        }
    }
    unt <- ifelse(su$layout$time.units=="hours", "h", "d")
    g.obj1 <- ggplot(data=dat[ind.sel,], mapping=aes(x=x, y=taus*su$layout$timestep.fac, shape=case, linetype=case)) + geom_point() + geom_line() + labs(x="", y=bquote(tau~"["*.(unt)*"]"), shape="Error Model", linetype="Error Model") + theme(axis.text=element_text(size=24)) + theme_bw(base_size=24) ##+ scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
    plot(g.obj1)
}
plot.dens.quantiles <- function(dat, sudriv, xlim=NA, ind.sel=NA, distr="norm"){
    ind.sel <- select.ind(sudriv=sudriv, xlim=xlim, ind.sel=ind.sel)
    n.case <- length(unique(dat$case))
    if(n.case > 1){
        for(i in 2:(n.case)){
            ind.sel <- c(ind.sel, ind.sel + (i-1)*sum(dat$case==dat$case[1]))
        }
    }
    if(distr=="norm"){
        g.obj1 <- ggplot(data=dat[ind.sel,], mapping=aes(sample=quant, shape=case, col=case)) + stat_qq() + stat_qq(geom="line") + geom_abline(slope=1, intercept=0) + expand_limits(y=-4) + theme_bw(base_size=24) + theme(axis.text=element_text(size=24), legend.text=element_text(size=21)) + labs(shape="Model", col="Model",  x="Normal theoretical quantiles", y=expression(eta))
        plot(g.obj1)
    }else{
        g.obj1 <- ggplot(data=dat[ind.sel,], mapping=aes(sample=unifvalue, linetype=case, col=case)) + stat_qq(distribution=stats::qunif, geom="line", size=1.5) + geom_abline(slope=1, intercept=0) + scale_x_continuous(breaks=seq(0,1,by=0.1), expand=c(0.001,0)) + scale_y_continuous(breaks=seq(0,1,by=0.1), expand=c(0.001,0))+ expand_limits(y=0) + theme_bw(base_size=24) + theme(axis.text=element_text(size=24), legend.text=element_text(size=21)) + labs(linetype="Model", col="Model", x="Uniform theoretical quantiles", y=expression(F[D[Q]](Q[obs]))) + expand_limits(x=c(0,1), y=c(0,1))
        print(summary(dat$unifvalue))
        g.obj2 <- ggplot(data=dat[ind.sel,], mapping=aes(x=unifvalue, linetype=case, col=case)) + geom_density(size=1.5, adjust = 1/3)+ theme_bw(base_size=24) + theme(axis.text=element_text(size=21)) + labs(x=expression(F[D[Q]](Q[obs])), linetype="Model", col="Model") + scale_x_continuous(breaks=seq(0,1,by=0.1), expand=c(0.001,0)) + scale_y_continuous(expand=c(0.001,0)) + geom_abline(slope=0, intercept=1) + expand_limits(x=c(0,1.01))
        plot(g.obj1)
        plot(g.obj2)
    }
}
plot.powspec.quantiles <- function(dat, sudriv, xlim=NA, ind.sel=NA, distr="norm"){
    ind.sel <- select.ind(sudriv=sudriv, xlim=xlim, ind.sel=ind.sel)
    n.case <- length(unique(dat$case))
    if(n.case > 1){
        cat("powspec function is not implemented for model comparison")
        return(NA)
    }
    dd <- dat[ind.sel,"quant"]
    len <- floor(length(dd)/20)
    ar <- numeric(length(dd))
    ar[1] <- rnorm(1)
    for(i in 2:length(dd)){
        ar[i] <- rnorm(1, 0.95*ar[i-1], sqrt(1-0.95^2))
    }
    sp=NULL;sp2=NULL
    for(i in 1:20){
        tmp <- spectrum(x=c(dd[((i-1)*len+1):(i*len)]),plot=FALSE)
        sp <- cbind(sp,tmp$spec)
        tmp2 <- spectrum(x=c(ar[((i-1)*len+1):(i*len)]),plot=FALSE)
        sp2 <- cbind(sp2,tmp2$spec)
        fr <- tmp$freq
    }
    sp <- rowMeans(sp)
    sp2 <- rowMeans(sp2)
    dat <- data.frame(spec.dens=c(sp,sp2), freq=c(fr,fr), type=rep(c("obs", "ar1"), each=length(fr)))
    g.obj1 <- ggplot(data=dat, mapping=aes(x=freq,y=spec.dens,col=type,shape=type)) + geom_point() + geom_line() + labs(x="Frequency [1/h]", y="Spectral density") + theme(axis.text=element_text(size=24)) + theme_bw(base_size=24) + scale_x_log10() + scale_y_log10()
    plot(g.obj1)
}
plot.dens.innovation <- function(dat, sudriv, xlim=NA, ind.sel=NA, distr="unif"){
    ind.sel <- select.ind(sudriv=sudriv, xlim=xlim, ind.sel=ind.sel)
    if(distr=="norm"){
        g.obj1 <- ggplot(data=dat[ind.sel,], mapping=aes(sample=innovation)) + stat_qq() + geom_abline(slope=1, intercept=0) + expand_limits(y=-4) + theme(axis.text=element_text(size=24)) + theme_bw(base_size=24)
        plot(g.obj1)
    }else{
        g.obj1 <- ggplot(data=dat[ind.sel,], mapping=aes(sample=innovation)) + stat_qq(distribution=stats::qunif) + geom_abline(slope=1, intercept=0) + scale_x_continuous(breaks=seq(0,1,by=0.1)) + scale_y_continuous(breaks=seq(0,1,by=0.1))+ expand_limits(y=0) + theme(axis.text=element_text(size=24)) + theme_bw(base_size=24)
        g.obj2 <- ggplot(data=dat[ind.sel,], mapping=aes(x=innovation)) + geom_density()+ theme(axis.text=element_text(size=24)) + theme_bw(base_size=24) + labs(x=expression(F[f(eta[i],eta[i-1])](eta[i]))) + scale_x_continuous(breaks=seq(0,1,by=0.1)) + expand_limits(x=0)
        plot(g.obj1)
        plot(g.obj2)
    }
}
plot.Qdet.quantiles <- function(dat, sudriv, xlim=NA, ind.sel=NA){
    dat[,"det"] <- dat[,"det"]/sudriv$layout$timestep.fac
    g.obj <- ggplot(data=dat[ind.sel,], mapping=aes(x=det, y=quant)) + geom_point()+ geom_abline(slope=0, intercept=0) + theme_bw(base_size=24) + theme(axis.text=element_text(size=24), axis.title=element_text(size=24)) + scale_x_log10() + labs(x=expression(Q["det"]*" (mm/h)"), y=expression(eta))
    plot(g.obj)
}
plot.pacf.quantiles <- function(dat, sudriv, xlim=NA, ind.sel=NA, lag.max=NULL, confidence=0.95){
    ind.sel <- select.ind(sudriv=sudriv, xlim=xlim, ind.sel=ind.sel)
    n.case <- length(unique(dat$case))
    if(n.case > 1){
        pac <- numeric()
        lag <- numeric()
        pac.obj <- pacf(dat$white.noise[ind.sel], plot=FALSE, lag.max=lag.max, na.action=na.pass)
        pac <- c(pac, c(1,pac.obj$acf))
        lag <- c(lag, c(0,pac.obj$lag))
        ind.new <- ind.sel
        for(i in 2:(n.case)){
            ind.new <- ind.new + sum(dat$case==dat$case[1])
            ind.sel <- c(ind.sel, ind.new)
            print(summary(ind.new))
            pac.obj <- pacf(dat$white.noise[ind.new], plot=FALSE, lag.max=lag.max, na.action=na.pass)
            pac <- c(pac, c(1,pac.obj$acf))
            lag <- c(lag, c(0,pac.obj$lag))
        }
    }else{
        pac.obj <- pacf(dat$white.noise[ind.sel], plot=FALSE, lag.max=lag.max, na.action=na.pass)
        pac <- c(1,pac.obj$acf)
        lag <- c(0,pac.obj$lag)
    }
    conf <- 2/sqrt(length(ind.sel))
    print(conf)
    dat.pac <- data.frame(pac=pac, lag=lag*sudriv$layout$timestep.fac, case=rep(unique(dat$case), each=length(pac.obj$lag)+1))
    g.obj <- ggplot(data=dat.pac, mapping=aes(x=lag, y=pac, shape=case, col=case)) + geom_point(size=2) + geom_line() + geom_hline(yintercept=conf, linetype="dotted") + geom_hline(yintercept=-1*conf, linetype="dotted") + geom_hline(yintercept=0) + scale_x_continuous(expand=c(0.001,0)) + theme_bw(base_size=24) + theme(axis.text=element_text(size=24), axis.title=element_text(size=24), legend.text=element_text(size=21)) + labs(x=ifelse(sudriv$layout$time.units=="days", "Lag [d]", "Lag [h]"), y=expression("PACF of "~eta[i]-E*"["*eta[i]*"|"*eta[i-1] * "]"), shape="Model", col="Model")
    plot(g.obj)
}
calc.metrics <- function(sudriv, xlim=NA, file=NA, ...){
    ind.sel <- select.ind(sudriv,xlim=xlim,ind.sel=NA)
    flash <- c(apply(sudriv$predicted$sample[,ind.sel],1,calc.flashiness))
    obs <- sudriv$observations[c(sudriv$layout$calib,sudriv$layout$pred)][ind.sel]
    flash.obs <- calc.flashiness(obs)
    det <- c(sampling_wrapper(sudriv, sample.par=FALSE, n.sample=1, sample.likeli=FALSE))[ind.sel]
    nse.det <- 1-sum((det-obs)^2)/sum((obs - mean(obs))^2)
    nse <- c(apply(sudriv$predicted$sample[,ind.sel],1,calc.nse,obs=obs))
    crps <- c(calc.crps(sudriv$predicted$sample[,ind.sel],obs=obs))
    crps <- crps/sudriv$layout$timestep.fac
    strmflw.tot <- c(apply(sudriv$predicted$sample[,ind.sel],1,sum))
    strmflw.err <- (sum(obs)-strmflw.tot)/sum(obs)*100
    df <- round(data.frame(nse.det=nse.det, nse.mu=mean(nse), nse.med=median(nse), nse.sd=sd(nse),
                     sferr.mu=mean(strmflw.err), sferr.med=median(strmflw.err), sferr.sd=sd(strmflw.err),
                     flash.mu=mean(flash), flash.med=median(flash), flash.sd=sd(flash), flash.obs=flash.obs
                     ), digits=3)
    if(is.na(file)){
        print(df)
    }else{
        write.table(df, file=file, sep="\t", ...)
    }
}
plot.flashiness <- function(dat, list.su, xlim=NA, ind.sel=NA){
    n.case <- length(unique(dat$case))
    ind.sel <- select.ind(list.su[[1]],xlim=xlim,ind.sel=ind.sel)
    flash <- c(apply(list.su[[1]]$predicted$sample[,ind.sel],1,calc.flashiness))
    n.samp <- nrow(list.su[[1]]$predicted$sample)
    if(n.case > 1){
        for(i in 2:(n.case)){
            flash <- c(flash,c(apply(list.su[[i]]$predicted$sample[,ind.sel],1,calc.flashiness)))
            n.samp <- c(n.samp,nrow(list.su[[i]]$predicted$sample))
        }
    }
    dat.flash <- data.frame(flash=flash, case=rep(unique(dat$case), times=n.samp))
    cat("FI: mean: ", with(dat.flash,tapply(flash,case,mean)), " median: ", with(dat.flash,tapply(flash,case,median))," sd: ", with(dat.flash,tapply(flash,case,sd)), "\n")
    flash.obs <- calc.flashiness(subset(dat, case==case[1])$obs[ind.sel])
    cat("Flashiness observed: ", flash.obs, "\n")
    g.obj <- ggplot(dat.flash, aes(x=case, y=flash)) + geom_boxplot() + labs(x="Model", y="Flashiness index [-]") + geom_hline(yintercept=flash.obs, col="red")+ theme_bw(base_size=24) ##+ theme(axis.text=element_text(size=24), axis.title=element_text(size=24))
    plot(g.obj)
}
plot.nse <- function(dat, list.su, xlim=NA, ind.sel=NA){
    n.case <- length(unique(dat$case))
    ind.sel <- select.ind(list.su[[1]],xlim=xlim,ind.sel=ind.sel)
    nse <- c(apply(list.su[[1]]$predicted$sample[,ind.sel],1,calc.nse,obs=list.su[[1]]$observations[c(list.su[[1]]$layout$calib,list.su[[1]]$layout$pred)][ind.sel]))
    n.samp <- nrow(list.su[[1]]$predicted$sample)
    if(n.case > 1){
        for(i in 2:(n.case)){
            nse <- c(nse,c(apply(list.su[[i]]$predicted$sample[,ind.sel],1,calc.nse,obs=list.su[[i]]$observations[c(list.su[[i]]$layout$calib,list.su[[i]]$layout$pred)][ind.sel])))
            n.samp <- c(n.samp,nrow(list.su[[i]]$predicted$sample))
        }
    }
    dat.nse <- data.frame(nse=nse, case=rep(unique(dat$case), times=n.samp))
    cat("NSE: mean: ", with(dat.nse,tapply(nse,case,mean)), " median: ", with(dat.nse,tapply(nse,case,median))," sd: ", with(dat.nse,tapply(nse,case,sd)), "\n")
    g.obj <- ggplot(dat.nse, aes(x=case, y=nse)) + geom_boxplot() + labs(x="Model", y="Nash-Sutcliffe Efficiency[-]")+ theme_bw(base_size=24) ##+ theme(axis.text=element_text(size=24), axis.title=element_text(size=24))
    plot(g.obj)
}
plot.crps <- function(dat, list.su, xlim=NA, ind.sel=NA){
    n.case <- length(unique(dat$case))
    ind.sel <- select.ind(list.su[[1]],xlim=xlim,ind.sel=ind.sel)
    crps <- c(calc.crps(list.su[[1]]$predicted$sample[,ind.sel],obs=list.su[[1]]$observations[c(list.su[[1]]$layout$calib,list.su[[1]]$layout$pred)][ind.sel]))
    n.samp <- ncol(list.su[[1]]$predicted$sample[,ind.sel])
    if(n.case > 1){
        for(i in 2:(n.case)){
            crps <- c(crps,c(calc.crps(list.su[[i]]$predicted$sample[,ind.sel],obs=list.su[[i]]$observations[c(list.su[[i]]$layout$calib,list.su[[i]]$layout$pred)][ind.sel])))
            n.samp <- c(n.samp,ncol(list.su[[i]]$predicted$sample[,ind.sel]))
        }
    }
    crps <- crps/list.su[[1]]$layout$timestep.fac
    dat.crps <- data.frame(crps=crps, obs=rep(list.su[[1]]$observations[c(list.su[[1]]$layout$calib,list.su[[1]]$layout$pred)][ind.sel], times=n.case), case=rep(unique(dat$case), times=n.samp))
    cat("CRPS: mean: ", with(dat.crps,tapply(crps,case,mean)), " median: ", with(dat.crps,tapply(crps,case,median))," sd: ", with(dat.crps,tapply(crps,case,sd)), "\n")
    g.obj <- ggplot(dat.crps, aes(x=case, y=crps)) + geom_boxplot() + labs(x="Model", y="CRPS [-]")+ theme_bw(base_size=24) ##+ theme(axis.text=element_text(size=24), axis.title=element_text(size=24))
    plot(g.obj)
    g.obj <- ggplot(dat.crps, aes(x=obs, y=crps)) + geom_point() + labs(x=expression(Q[obs]), y="CRPS [-]")+ theme_bw(base_size=24) ##+ theme(axis.text=element_text(size=24), axis.title=element_text(size=24))
    plot(g.obj)
}
plot.strmflw.err <- function(dat, list.su, xlim=NA, ind.sel=NA){
    n.case <- length(unique(dat$case))
    ind.sel <- select.ind(list.su[[1]],xlim=xlim,ind.sel=ind.sel)
    strmflw.tot <- c(apply(list.su[[1]]$predicted$sample[,ind.sel],1,sum))
    n.samp <- nrow(list.su[[1]]$predicted$sample)
    if(n.case > 1){
        for(i in 2:(n.case)){
            strmflw.tot <- c(strmflw.tot,c(apply(list.su[[i]]$predicted$sample[,ind.sel],1,sum)))
            n.samp <- c(n.samp,nrow(list.su[[i]]$predicted$sample))
        }
    }
    strmflw.obs <- sum(subset(dat, case==case[1])$obs[ind.sel])
    print(strmflw.obs)
    dat.strmflw.err <- data.frame(strmflw.err=(strmflw.obs-strmflw.tot)/strmflw.obs*100, case=rep(unique(dat$case), times=n.samp))
    cat("Strmflw Err: mean: ", with(dat.strmflw.err,tapply(strmflw.err,case,mean)), " median: ", with(dat.strmflw.err,tapply(strmflw.err,case,median))," sd: ", with(dat.strmflw.err,tapply(strmflw.err,case,sd)), "\n")
    g.obj <- ggplot(dat.strmflw.err, aes(x=case, y=strmflw.err)) + geom_boxplot() + labs(x="Model", y="Total streamflow error [%]") + geom_hline(yintercept=0, col="red")+ theme_bw(base_size=24) + theme(axis.text=element_text(size=24), axis.title=element_text(size=24))
    plot(g.obj)
}
plot.sd <- function(sudriv, fsd=1){
    Q=exp(seq(-6.9, 2.3, 0.1))
    par.likeli<- ifelse(as.logical(sudriv$likelihood$tran), exp(sudriv$likelihood$parameters), sudriv$likelihood$parameters)
    names(par.likeli) <- names(sudriv$likelihood$parameters)
    a <- par.likeli["C1Wv_Qstream_a_lik"]
    b <- par.likeli["C1Wv_Qstream_b_lik"]
    c <- par.likeli["C1Wv_Qstream_c_lik"]
    sd0 <- par.likeli["C1Wv_Qstream_sd0_lik"]
    Q0 <- par.likeli["C1Wv_Qstream_Q0_lik"]
    if(fsd==1) sdd <- a*sd0*((Q/Q0+b)/(1+b))^c
    if(fsd==2) sdd <- a*sd0*(Q/Q0)^c + b*Q0
    if(fsd==4){
        sdd <- numeric(length(Q))
        sdd[Q<Q0] <-  a*sd0*Q[Q<Q0]/Q0 + b*Q0
        sdd[Q>=Q0] <- b*Q0 + a*sd0/c*(c-1+(Q[Q>=Q0]/Q0)^c)
    }
    dat.sd <- data.frame(Q=Q, sd=sdd)
    g.obj1 <- ggplot(data=dat.sd, mapping=aes(x=Q, y=sd)) + geom_line() + scale_x_log10() + scale_y_log10() + theme_bw(base_size=24) + theme(axis.text=element_text(size=24), axis.title=element_text(size=24)) + labs(x="Q (mm/d)", y="Standard deviation (mm/d)")
    g.obj2 <- ggplot(data=dat.sd, mapping=aes(x=Q, y=sd)) + geom_line() + theme_bw(base_size=24) + theme(axis.text=element_text(size=24), axis.title=element_text(size=24)) + labs(x="Q (mm/h)", y="Standard deviation (mm/d)")
    plot(g.obj1)
    plot(g.obj2)
}
tau.time <- function(dat, sudriv, obspred="obspred", calc.S=FALSE){
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
plot.tau.time <- function(dat, sudriv, xlim=NA, ind.sel=NA, tme.orig="1985-01-01"){
    ind.sel <- select.ind(sudriv=sudriv, xlim=xlim, ind.sel=ind.sel)
    n.case <- length(unique(dat$case))
    if(n.case > 1){
        for(i in 2:(n.case)){
            ind.sel <- c(ind.sel, ind.sel + (i-1)*sum(dat$case==dat$case[1]))
        }
    }
    dat.obsdet <- rbind(dat[ind.sel,c("x","det")], setNames(dat[ind.sel,c("x", "obs")], nm=c("x", "det")))
    dat.obsdet <- rbind(dat.obsdet, setNames(dat[ind.sel,c("x", "tau.calc")], nm=c("x", "det")))
    dat.obsdet <- rbind(dat.obsdet, setNames(dat[ind.sel,c("x", "tau.obs")], nm=c("x", "det")))
    dat.obsdet <- rbind(dat.obsdet, setNames(dat[ind.sel,c("x", "tau.obs.rel")], nm=c("x", "det")))
    dat.obsdet <- rbind(dat.obsdet, setNames(dat[ind.sel,c("x", "S")], nm=c("x", "det")))
    dat.obsdet$var <- rep(c("det","obs","tau.det","tau.obs","tau.obs.rel","S"), each=length(ind.sel))
    p.all <- sudriv$input$inputobs[sudriv$layout$layout$time[c(sudriv$layout$calib, sudriv$layout$pred)], c(1,2)] ## ATTENTION: we assume that the input and the layou time are the same...
    ##p.all[,2] <- rollmean(p.all[,2], k=5, fill=0)
    p <- as.data.frame(p.all[ind.sel,])
    p1 <- as.POSIXct(x=as.numeric(p[,1])*60*60, origin=tme.orig)
    p[,1] <- p1
    p <- as.data.frame(cbind(p, rep("precip", nrow(p))), row.names=NULL)
    dat.obsdet <- rbind(dat.obsdet, setNames(p, nm=names(dat.obsdet)), make.row.names=FALSE)
    names(dat.obsdet) <- c("x", "y", "var")
    gobj1 <- ggplot(data=subset(dat.obsdet, var %in% c("obs","det")), mapping=aes(x=x, y=y,shape=var,colour=var)) + geom_point(size=0.5) + geom_line() + labs(x="", y="Discharge (mm/h)") + theme(axis.text=element_text(size=24)) + theme_bw(base_size=24) + scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
    gobj2 <- ggplot(data=subset(dat.obsdet, grepl("tau", var)), mapping=aes(x=x, y=y, col=var, shape=var, linetype=var)) + geom_point(size=0.5) + geom_line() + labs(x="", y=expression(tau)) + theme_bw(base_size=24) + theme(axis.text=element_text(size=24)) + scale_x_datetime(date_breaks="1 week", date_labels="%d %b") + scale_y_continuous(limits=c(0,50))
    gobj3 <- ggplot(data=subset(dat.obsdet, var %in% c("precip","S")), mapping=aes(x=x, y=y, col=as.numeric(y>0), shape=var)) + geom_point(size=0.5) + geom_line() + labs(x="", y="P and S (mm/h)") + theme_bw(base_size=24)+ theme(axis.text=element_text(size=24)) + scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
    weight.obs <- 1/(dat.obsdet$y[dat.obsdet$var=="tau.obs"]+1)
    ##weight.obs <- weight.obs/max(weight.obs)
    weight.det <- 1/(dat.obsdet$y[dat.obsdet$var=="tau.det"]+1)
    ##weight.det <- weight.det/max(weight.det)
    dat.obsdet$weight.obs <- NA
    dat.obsdet$weight.det <- NA
    dat.obsdet[dat.obsdet$var=="obs","weight.obs"] <- 5*weight.obs
    dat.obsdet[dat.obsdet$var=="obs","weight.det"] <- 5*weight.det
    gobj4 <- ggplot(data=subset(dat.obsdet, var=="obs"),mapping=aes(x=x, y=y)) + geom_point(alpha=weight.obs)+ geom_line(alpha = weight.obs) + labs(x="", y="Discharge (mm/h)") + theme(axis.text=element_text(size=24)) + theme_bw(base_size=24) + scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
    gobj5 <- ggplot(data=subset(dat.obsdet, var=="obs"), mapping=aes(x=x, y=y)) + geom_point(alpha=weight.det) +geom_line(alpha=weight.det) + labs(x="", y="Discharge (mm/h)") + theme(axis.text=element_text(size=24)) + theme_bw(base_size=24) + scale_x_datetime(date_breaks="1 week", date_labels="%d %b")

    ##deriv <- c(diff(dat$det),NA)/c(diff(p.all[,1]), NA)
    deriv <- c(diff(dat$S),NA)/c(diff(p.all[,1]), NA)
    dat.ptq <- data.frame(x=p.all[,1], p=p.all[,2], p1=c(NA,p.all[1:(nrow(p.all)-1),2]), tau.obs=dat$tau.obs, discharge=dat$S, deriv=deriv, deriv.rel=deriv/dat$S, deriv1=p.all[,2]*deriv/dat$S, deriv2=c(NA,p.all[1:(nrow(p.all)-1),2])*deriv/dat$S, P2=(p.all[,2]-0.3)/dat$S, P3=c(NA,NA,(p.all[,2]-0.3)/dat$S)[1:nrow(p.all)], tau.obs.rel=dat$tau.obs.rel, tau.calc=dat$tau.calc)

    gobj6 <- ggplot(dat.ptq, aes(x=p,y=tau.obs)) + geom_point(size=0.3) + scale_y_log10(limits=c(1e-1,NA)) + labs(x="Precipitation (mm/h)", y=expression(tau))

    gobj7 <- ggplot(subset(dat.ptq, deriv>0), aes(x=deriv1, y=tau.obs)) + geom_point(size=0.3) + scale_x_log10(limits=c(1e-5,1e1)) + scale_y_log10(limits=c(1e-1,NA)) + labs(x="(positive P*dS/dt/S)", y=expression(tau))

    gobj8 <- ggplot(dat.ptq[p.all[,2]==0 & dat.ptq$deriv<0,], aes(x=deriv*-1, y=tau.obs)) + geom_point(size=0.3) + scale_x_log10() + scale_y_log10(limits=c(1e-1,NA)) + labs(x="negative dS/dt when P=0", y=expression(tau))

    ## gobj9 <- ggplot(dat.ptq, aes(x=p1,y=tau.obs)) + geom_point(size=0.3) + scale_y_log10(limits=c(1e-1,NA)) + labs(x="Precipitation(i-1) (mm/h)", y="tau * S/sd")
    gobj9 <- ggplot(dat.ptq, aes(x=tau.calc, y=tau.obs)) + geom_point(size=0.3) + scale_x_log10(limits=c(1e0,NA)) + scale_y_log10() + labs(x="tau calc", y=expression(tau))

    ## gobj10 <- ggplot(dat.ptq, aes(x=P3, y=tau.obs)) + geom_point(size=0.3) + scale_x_log10(limits=c(1e-2,1e1)) + scale_y_log10(limits=c(1e-1,NA)) + labs(x="(P-0.3)(i-2)/S(i-2)", y=expression(tau))
    gobj10 <- ggplot(dat.ptq, aes(x=tau.calc, y=tau.obs.rel)) + geom_point(size=0.3) + scale_x_log10(limits=c(1e-1,NA)) + scale_y_log10(limits=c(1e-1,NA)) + labs(x="tau calc", y="tau*Qdet/sd")

    gobj11 <- ggplot(dat.ptq[p.all[,2]==0,], aes(x=discharge, y=tau.obs)) + geom_point(size=0.3) + scale_x_log10(limits=c(1e-5,NA)) + scale_y_log10() + labs(x="S when P=0", y=expression(tau))

    dat.tmp <- dat.ptq[ind.sel,]
    col      <-  pmin(pmax(log(dat.tmp$tau.calc/(dat.tmp$tau.obs+0.01)),-5),5)
    ##ind.out3 <-  dat.tmp$tau.calc/(dat.tmp$tau.obs+0.01) > 10

    grid.newpage()
    grid.draw(gtable_rbind(ggplotGrob(gobj1),ggplotGrob(gobj2),ggplotGrob(gobj3)))
    ## grid.newpage()
    ## grid.draw(gtable_rbind(ggplotGrob(gobj4),ggplotGrob(gobj5)))
    ## grid.newpage()
    ## grid.draw(gtable_rbind(gtable_cbind(ggplotGrob(gobj6),ggplotGrob(gobj7),ggplotGrob(gobj8)),
    ##                        gtable_cbind(ggplotGrob(gobj9),ggplotGrob(gobj10),ggplotGrob(gobj11))))

    ## col <- numeric(length(ind.sel))
    ## col[ind.out2] <- -1
    ## col[ind.out3] <- 1
    dat.tmp <- data.frame(x=rep((1:nrow(dat))[ind.sel],2), y=c(dat$obs[ind.sel], dat$det[ind.sel]), var=rep(c("obs","det"), each=length(ind.sel)), deviance=rep(col, 2))
    gobj1 <- ggplot(dat.tmp, aes(x=x,y=y,color=deviance)) + geom_point(data=subset(dat.tmp, var=="obs"), size=0.5) + geom_line(data=subset(dat.tmp, var=="det")) + labs(x="", y="Discharge (mm/h)") + theme(axis.text=element_text(size=24)) + theme_bw(base_size=24) +scale_colour_gradient2(low="blue", mid="grey", high="red")##+ scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
    gobj3 <- ggplot(data=data.frame(x=p.all[ind.sel,1], P=p.all[ind.sel,2]), mapping=aes(x=x, y=P, colour=as.numeric(P>0))) + geom_point(size=0.5) + geom_line() + labs(x="", y="P (mm/h)") + theme_bw(base_size=24)+ theme(axis.text=element_text(size=24), legend.position="none")## + scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
    plot(gobj1)
    plot(gobj3)
}

constrain_parameters_wrapper <- function(sudriv, mcmc.sample){
    sc <- mcmc.sample
    ## account for the constrained parameters
    lb <- c(sudriv$model$args$parLo[as.logical(sudriv$model$par.fit)], sudriv$likelihood$lb[as.logical(sudriv$likelihood$par.fit)])
    ub <- c(sudriv$model$args$parHi[as.logical(sudriv$model$par.fit)], sudriv$likelihood$ub[as.logical(sudriv$likelihood$par.fit)])
    for(i in 1:(dim(mcmc.sample)[3])){
        sc[,,i] <- constrain_parameters(mcmc.sample[,,i], lb[i], ub[i])
    }
    return(sc)
}
