prepare.timedepargs <- function(su,tag,which.timedep,remove.taumax,fix.taumax,fix.a,f_mean,sclshifts,fit.sd.ou=FALSE,hybrid=FALSE,mod=NULL,var=NULL,scaleshift=NULL){
  ## prepare the arguments needed for the infer.timedep function of Peter Reicherts framework
  if(hybrid & (is.null(mod) | is.null(var) | is.null(scaleshift))) stop("'mod', 'var', and 'scaleshift' have to be suppled when hybrid=TRUE")
  if(remove.taumax){
    su <- remove.taumax.Q(su)
  }else if(fix.taumax){
    su <- fix.taumax.Q(su)
  }
  if(fix.a) su <- fix.a.func(su)
  su$predicted <- NULL ## drop predictions to reduce object size, because su object is saved as argument of timedeppar result.
  
  if("Glo%Cmlt_P" %in% which.timedep){ ## make sure that it is among the fitted parameters
    su$model$par.fit[names(su$model$parameters) == "Glo%Cmlt_P"] <- 1
    su$model$prior$distdef[["Glo%Cmlt_P"]] <- c("normaltrunc","1","0.05","0.5","2")
  }else if("Glo%tStart_VL" %in% which.timedep){ ## make sure that it is among the fitted parameters
    su$model$par.fit[names(su$model$parameters)=="Glo%tStart_VL"] <- 1
    su$model$prior$distdef[["Glo%tStart_VL"]] <- c("normaltrunc","1","0.5","0","2")
  }
  param.ini <- as.list(c(su$model$parameters[as.logical(su$model$par.fit)], su$likelihood$parameters[as.logical(su$likelihood$par.fit)]) - ifelse(grepl("Cl[0-9]",tag), runif(1,0.1,0.2), 0))
  if(f_mean){
    mnprm <- paste0(which.timedep, "_fmean")
  }else{
    mnprm <- NA
  }
  param.log <- logical(0)
  par.tran <- as.logical(su$model$args$parTran[names(su$model$parameters) %in% which.timedep])
  names(par.tran) <- names(su$model$parameters)[names(su$model$parameters) %in% which.timedep]
  sds <- c("Glo%Cmlt_Dspl_SD" = 0.1,
           "Glo%CmltSmax_UR" = 0.3,
           "Glo%CmltSmax_IR" = 0.15,
           "Glo%Cmlt_P" = 0.1,
           "Glo%Cmlt_BeQq_UR" = 0.1,
           "GloTr%CmltSlOne_IR" = 0.3,
           "GloTr%CmltSlTwo_IR" = 0.3,
           "Glo%Cmlt_AlQq_FR" = 0.1,
           "Glo%Cmlt_AlQq_SR" = 0.1,
           "Glo%Cmlt_K_Qq_FR" = 0.2,
           "Glo%Cmlt_K_Qq_RR" = 0.2,
           "Glo%Cmlt_K_Qq_SR" = 0.2,
           "U1W%KpQq_FR"      = 0.05,
           "Glo%Cmlt_Pmax_ED" = 0.2,
           "Glo%tStart_VL" = 0.7,
           "GloTr%CmltKd_WR" = 0.1,
           "GloTr%CmltRs_WR" = 0.2,
           "Glo%Cmlt_E" = 0.1)
  ## special cases:
  if("Glo%CmltSmax_IR" %in% which.timedep){
    su$model$prior$distdef[["Glo%CmltSmax_IR"]][4] <- log(1) # set lower bound to 1 mm
  }else if("Glo%Cmlt_P" %in% which.timedep){
    ## redefine bounds of precipitation multiplier
    su$model$args$parLo[names(su$model$parameters)=="Glo%Cmlt_P"] <- 0.5
    su$model$args$parHi[names(su$model$parameters)=="Glo%Cmlt_P"] <- 2
  }else if("Glo%Cmlt_Pmax_ED" %in% which.timedep){
    ## relax the upper bound of this one
    su$model$args$parHi[names(su$model$parameters)=="Glo%Cmlt_Pmax_ED"] <- exp(6)
  }else if("Glo%Cmlt_K_Qq_RR" %in% which.timedep){
    ## relax the upper bound of this one
    su$model$args$parHi[names(su$model$parameters)=="Glo%Cmlt_K_Qq_RR"] <- exp(3)
  }
  param.ou.ini <- numeric()
  param.ou.fixed <- numeric()
  if(!all(which.timedep=="none")){
    for(td.curr in which.timedep){ ## construct param.ini, param.ou.fixed, param.ou.ini
      if(f_mean){ # re-parameterize mean as constant parameter
        if(td.curr %in% names(sclshifts)){
          param.ini[[td.curr]] <- cbind(1:nrow(su$input$inputobs), 0)
          param.ini[[paste0(td.curr,"_fmean")]] <- sigm.trans.inv(ifelse(td.curr=="Glo%tStart_VL", 1, as.numeric(su$model$parameters[td.curr])), sclshifts[[td.curr]][1],sclshifts[[td.curr]][2])
          mn.fix <- 0
          sd.fix <- sds[td.curr] ## fixed standard deviation of the rescaled OU-process
        }else{
          param.ini[[td.curr]] <- cbind(1:nrow(su$input$inputobs), ifelse(par.tran[td.curr], 0, 1))
          param.ini[[paste0(td.curr,"_fmean")]] <- as.numeric(su$model$parameters[td.curr])
          mn.fix <- ifelse(par.tran[td.curr], 0, 1)
          sd.fix <- sds[td.curr] ## fixed standard deviation of the rescaled OU-process
        }
        fix.new <- c(mn.fix, sd.fix, 1/(10*24*4)) # in 15 min units (gamma = 1/tau)
        names(fix.new) <- c(paste0(td.curr,"_mean"), paste0(td.curr,"_sd"), paste0(td.curr,"_gamma"))
        if(fit.sd.ou){
          param.ou.fixed <- c(param.ou.fixed, fix.new[c(1,3)])
          param.ou.ini   <- c(param.ou.ini, fix.new[2])
        }else{
          param.ou.fixed <- c(param.ou.fixed, fix.new)
        }
        ## make prior for fmean param, just to construct the range further down. This prior is removed again.
        su$model$prior$distdef[[paste0(td.curr,"_fmean")]] <- su$model$prior$distdef[[td.curr]]
      }else{
        warning("not reparametrizing mean (!f_mean) has not been maintained")
        if(td.curr %in% names(sclshifts)){
          if(td.curr=="Glo%tStart_VL"){
            param.ini[[td.curr]] <- cbind(1:nrow(su$input$inputobs), sigm.trans.inv(1,sclshifts[[td.curr]][1],sclshifts[[td.curr]][2]))
            ini.new <- sigm.trans.inv(1,sclshifts[[td.curr]][1],sclshifts[[td.curr]][2])
            fix.new <- sds[td.curr]
          }else{
            param.ini[[td.curr]] <- cbind(1:nrow(su$input$inputobs), sigm.trans.inv(as.numeric(su$model$parameters[td.curr]),sclshifts[[td.curr]][1],sclshifts[[td.curr]][2]))
            ini.new <- sigm.trans.inv(as.numeric(su$model$parameters[td.curr]),sclshifts[[td.curr]][1],sclshifts[[td.curr]][2])
            fix.new <- sds[td.curr]
          }
        }else{
          param.ini[[td.curr]] <- cbind(1:nrow(su$input$inputobs), as.numeric(su$model$parameters[td.curr]))
          ini.new <- as.numeric(su$model$parameters[td.curr])
          fix.new <- sds[td.curr]
        }
        names(ini.new) <- paste0(td.curr,"_mean")
        param.ou.ini <- c(param.ou.ini, ini.new)
        fix.new <- c(fix.new, 1/(10*24*4)) # in 15 min units (gamma = 1/tau)
        names(fix.new) <- c( paste0(td.curr,"_sd"), paste0(td.curr,"_gamma"))
        param.ou.fixed <- c(param.ou.fixed, fix.new)
      }
    }
    su$model$timedep$pTimedep <- names(su$model$parameters) %in% which.timedep # indicate which parameter is time-dependent
  }
  
  ## prepare param.range based on prior
  param.range <- param.ini
  for(i in names(param.range)){
    if(!grepl("lik", i)){
      dist <- su$model$prior$distdef[[i]]
    }else{
      dist <- su$likelihood$prior$distdef[[i]]
    }
    if(dist[1]=="uniform"){
      param.range[[i]] <- as.numeric(dist[c(2,3)])
      if(!any(sapply(which.timedep, grepl, x=i))) param.ini[[i]] <- max(min(param.ini[[i]], param.range[[i]][2]), param.range[[i]][1]) ## force the initial constant parameter within bounds
    }else if(grepl("trunc", dist[1])){
      param.range[[i]] <- as.numeric(dist[c(4,5)])
      if(!any(sapply(which.timedep, grepl, x=i))) param.ini[[i]] <- max(min(param.ini[[i]], param.range[[i]][2]), param.range[[i]][1])
    }else{
      param.range[[i]] <- NULL
    }
  }
  su$model$prior$distdef[paste0(which.timedep,"_fmean")] <- NULL
  if(f_mean & !all(which.timedep=="none")){
    for(i in which.timedep){
      if(par.tran[i]){
        param.range[[i]] <- c(log(1/5),log(5))
      }else{
        param.range[[i]] <- c(1/5,5)
      }
    }
  }
  for (td.curr in which.timedep){
    if(td.curr %in% names(sclshifts)){
      param.range[[td.curr]] <- c(-20,20)
      param.range[[paste0(td.curr,"_mean")]] <- c(-5,5)
      ## adapt boundaries of mean scaling factor
      if(f_mean) param.range[[paste0(td.curr,"_fmean")]] <- c(-20,20)
    }
  }
  if(hybrid){
    hybrid.args <- prepare.hybrid.args(sudriv=su, tag=tag, mod=mod, var=var, scaleshift=scaleshift)
  }else{
    hybrid.args <- NULL
  }
  print("range:")
  print(param.range)
  print("param.ini.timedep:")
  print(head(param.ini[[which.timedep]]))
  timdep.args <- list(param.ini = param.ini,
                      param.range = param.range,
                      param.log = param.log,
                      param.ou.ini = param.ou.ini,
                      param.ou.fixed = param.ou.fixed,
                      su = su,
                      mnprm=mnprm,
                      hybrid=hybrid,
                      hybrid.args=hybrid.args)
  return(timdep.args)
}
select.maxlikpars.timedep <- function(sudriv, res.timedep, scaleshift=NA, lik.not.post=FALSE){ # update sudriv object with maximum posterior timedependent parameters
  if(!all(is.na(scaleshift)) & is.null(rownames(scaleshift))) stop("no scaleshift rownames supplied")
  ## lik.not.post: select maximum likelihood parameter instead of maximum posterior.
  ind.timedep <- unlist(lapply(res.timedep$param.maxpost, length))>1
  ## get index of maximum posterior
  pm <- which.max(res.timedep$sample.logpdf[,ifelse(lik.not.post,"loglikeliobs","logposterior")])
  par <- res.timedep$sample.param.const[pm,]
  fmn <- grepl("_fmean", x = names(par))
  if(any(fmn)){
    print("tran:")
    print((names(par)[fmn] %in% rownames(scaleshift)))
    tran <- sudriv$model$args$parTran[names(sudriv$model$parameters)==gsub("_fmean","",names(par)[fmn])] == 1 | (names(par)[fmn] %in% rownames(scaleshift))
    fmn.val <- par[fmn]
  }
  
  ## update the maximum posterior time-dependent parameters
  if(sum(ind.timedep)>0){
    partd <- res.timedep$sample.param.timedep[[1]][-1,][pm,] ## this only works for one time-dependent parameter so far
    if(any(fmn)){
      if(tran){
        partd <- partd + fmn.val
      }else{
        partd <- partd * fmn.val
      }
    }
    parmat <- matrix(partd, nrow=length(partd))
    if(sum(sudriv$model$timedep$pTimedep)!=sum(ind.timedep)) stop("sudriv and res.timedeppar do not have the same number of timdep parameters")
    if(any(names(sudriv$model$parameters)[sudriv$model$timedep$pTimedep] != names(res.timedep$param.maxpost[ind.timedep]))) stop("pTimedep of sudriv and res.timedep do not have the same timedependent parameters")
    parmat <- as.matrix(parmat)
    colnames(parmat) <- NULL
    if(!all(is.na(scaleshift))){
      if(nrow(scaleshift)!=sum(ind.timedep) | ncol(scaleshift)!=2) stop("dimension of scaleshift is not right")
      for(i in 1:ncol(parmat)){
        parmat[,i] <- sigm.trans(parmat[,i], scale=scaleshift[i,1], shift=scaleshift[i,2])
        if(any(fmn)) par[which(fmn)[i]] <- sigm.trans(par[which(fmn)[i]], scale=scaleshift[i,1], shift=scaleshift[i,2]) #constant parameters
      }
    }
    ## force time course within bounds after addition or multiplication with fmean parameter
    lo <- sudriv$model$args$parLo[sudriv$model$timedep$pTimedep]
    hi <- sudriv$model$args$parHi[sudriv$model$timedep$pTimedep]
    for(i in 1:ncol(parmat)){
      parmat[,i] <- pmin(pmax(parmat[,i], lo[i]), hi[i])
      if(any(fmn)) par[which(fmn)[i]] <- pmin(pmax(par[which(fmn)[i]], lo[i]), hi[i]) #constant parameters
    }
    sudriv$model$timedep$par <- parmat
  }
  
  ## update maximum posterior of the constant parameters
  names(par) <- gsub("_fmean","",names(par))
  print("updated constant parameters to:")
  print(par)
  match.m <- match(names(par), names(sudriv$model$parameters))
  match.l <- match(names(par), names(sudriv$likelihood$parameters))
  sudriv$model$parameters[match.m[!is.na(match.m)]] <- par[!is.na(match.m)]
  sudriv$likelihood$parameters[match.l[!is.na(match.l)]] <- par[!is.na(match.l)]
  
  return(sudriv)
}
accept.frequ.get <- function(res, n.burnin=0){
  thin    <- 1; if ( !is.null(res$control) ) { if(!is.null(res$control$thin))    thin    <- res$control$thin } 
  n.adapt <- 0; if ( !is.null(res$control) ) { if(!is.null(res$control$n.adapt)) n.adapt <- min(res$control$n.adapt,res$control$n.iter-1) }
  freq <- list()
  for ( i in 1:length(res$sample.param.timedep) )
  {
    t <- res$sample.param.timedep[[i]][1,]
    sample <- signif(res$sample.param.timedep[[i]][-1,],digits=6)
    ind.range <- 1:nrow(sample)
    if ( floor(max(n.adapt,n.burnin)/thin)+1 < nrow(sample) ) ind.range <- (floor(max(n.adapt,n.burnin)/thin)+1):nrow(sample)
    freq.tmp <- rep(NA,length(t))
    for ( j in 1:length(t) ) freq.tmp[j] <- (length(unique(sample[ind.range,j]))-1)/length(ind.range)*100
    freq[[i]] <- freq.tmp
  }
  names(freq) <- names(res$sample.param.timedep)
  return(freq)
}

prepare.hybrid.args <- function(sudriv, tag, mod, var, scaleshift){
  layout.model.td <- sudriv$layout
  layout.model.td$layout <- data.frame(var=rep(var, each=nrow(sudriv$input$inputobs)), time=rep(sudriv$input$inputobs[,1], times=length(var)))
  layout.model.td$lump <- NA
  layout.model.td$calib <- 1:nrow(layout.model.td$layout)
  layout.model.td$pred.layout <- layout.model.td$layout
  ## get feature data (model states of this will be changed while iterating)
  data <- get.loess.input(sudriv=sudriv, tag=tag, vars=var, t.lim=c(-Inf,Inf), remove.na=FALSE, with.td=FALSE)
  data <- data[,"U5F1Wv_Ss1",drop=FALSE]
  ## run hybrid model
  model.td <- function(dat, mod){
    colnames(dat) <- gsub("U5F1Wv_Ss1", "x", colnames(dat))
    out <- dat[,"x"] < min(mod$x) | dat[,"x"] > max(mod$x)
    if(any(out)) warning(paste0("Input for empirical model outside range at ", sum(out), " positions. Truncating ..."))
    dat[,"x"] <- pmin(pmax(dat[,"x"], min(mod$x)), max(mod$x))
    y <- predict(mod, newdata=dat)
    if(any(!is.finite(y))){
      print("y is not finite at:")
      print(summary(dat[!is.finite(y),]))
    }
    return(matrix(y, ncol=1))
  }
  args <- list(data=data,
               data.time=data$time,
               model.td=model.td,
               args.model.td=list(mod=mod),
               lstm=FALSE,
               scaleshift=scaleshift,
               layout.model.td=layout.model.td,
               verbose=1)
  return(args)
}

plot.timedeppar.dynamics <- function(res.timedeppar, burn.in=0, plot=TRUE, file=NA, conf=c(0.6,0.8,0.9), time.info=list(tme.orig=NA,t0=NA,tme.units=NA,timestep.fac=NA), xlim=c(-Inf,Inf), xintercept=NULL, tag.red=NULL, applic=FALSE){
  ## this function plots the temporal dynamics of the time course of a parameter estimated with the infer.timedeppar function
  ## get minimum and maximum of intervals to shade (assuming xintercept is a list)
  if(!is.null(xintercept)){
    tmp <- sapply(xintercept, function(x) x[1], simplify=FALSE)
    xint.min <- unlist(tmp)
    attributes(xint.min) <- attributes(tmp[[1]])
    tmp <- sapply(xintercept, function(x) x[2], simplify=FALSE)
    xint.max <- unlist(tmp)
    attributes(xint.max) <- attributes(tmp[[1]])
  }
  
  su.tmp <- res.timedeppar$dot.args$sudriv
  td <- transform.timedep.par.sample(res.timedeppar$sample.param.timedep, res.timedeppar$sample.param.const, su.tmp,
                                     res.timedeppar$dot.args$mnprm, res.timedeppar$dot.args$scaleshift)
  for(tdcurr in names(td)){ # don't forget log transformation
    if(su.tmp$model$args$parTran[names(su.tmp$model$parameters)==tdcurr]==1){
      td[[tdcurr]][2:nrow(td[[tdcurr]]),] <- exp(td[[tdcurr]][2:nrow(td[[tdcurr]]),])
    }
  }
  ## transform all to data frames
  burn.in <- burn.in/res.timedeppar$control$thin
  td <- lapply(td, function(x) as.data.frame(t(x[c(1,(burn.in+2):nrow(x)),]))) # always keep the first row, since it is the time, not a sample
  ## make one table for all the parameters combined
  td <- td %>% enframe() %>% unnest(cols=c(value)) %>% rename(time=V1)
  ## transform time into time units
  if(!all(is.na(unlist(time.info)))){
    if(any(is.na(unlist(time.info)))) warning("NA in 'time.info', cannot deal with that...")
    td <- td %>% mutate(time=as.POSIXct(time.info$tme.orig)+
                                            ((time-1)*time.info$timestep.fac+time.info$t0)*60*60*ifelse(time.info$tme.units=="days",24,1))
  }
  ## and limit time axis
  td <- td %>% filter(time >= xlim[1] & time <=xlim[2])

  ## make the table longer
  td <- td %>% pivot_longer(c(-name,-time), names_to="k", values_to="value")
  n.samp <- length(unique(td$k))
  ## make quantiles based on confs
  qnts <- c((1-conf)/2, (1-(1-conf)/2)[length(conf):1])
  ## get quantiles for each timestep
  bounds <- td %>%
    nest(data=c(k,value)) %>%
    mutate(quant.tmp=map(data,~quantile(.$value, probs = qnts)), quant=map(quant.tmp, ~bind_rows(.) %>%
                                                                             pivot_longer(cols=everything(), names_to="nm.quant", values_to="val.quant"))) %>% unnest(quant)
  bounds <- bounds %>% select(-data,-quant.tmp)
  ## prepare the mapping of quantiles to confidence levels
  tmp <- names(quantile(rnorm(100), probs = qnts))
  mp <- c(conf, conf[length(conf):1])
  names(mp) <- tmp
  ## map quantiles to confidence levels
  bounds <- bounds %>% mutate(conf=mp[nm.quant], lw.up=case_when(nm.quant%in%names(mp)[1:(length(mp)/2)] ~ "lower",
                                                                 nm.quant%in%names(mp)[(length(mp)/2+1):length(mp)] ~ "upper"))
  ## make wider
  bounds.wide <- bounds %>% select(-nm.quant) %>% pivot_wider(names_from=lw.up, values_from=val.quant)
  ## manual alpha scale mapping
  alp <- pmin(1-bounds.wide$conf+0.1,1)
  names(alp) <- bounds.wide$conf
  tmp <- make.breaks(xlim)
  if(all(!is.finite(xlim))) xlim <- NULL
  gg <- ggplot(bounds.wide%>%mutate(conf=as.character(conf)), aes(x=time)) +
    geom_ribbon(mapping = aes(ymin=lower, ymax=upper, alpha=conf)) + scale_alpha_manual(values=alp) + 
    scale_x_datetime(date_breaks=tmp$brks, date_labels=tmp$frmt, limits=xlim) + 
    labs(alpha="Confidence", x="Time", y=ifelse(is.null(tag.red),"parameter",mylabeller.param.units(tag.red)), 
         caption=paste0("Based on ",n.samp," samples")) + theme_light() + 
    theme(plot.margin=unit(c(0,0.01,0.1,0.1), "cm"), text=element_text(size=12),
          legend.title=element_text(size=14), legend.text=element_text(size=14))
  if(applic) gg <- gg + geom_vline(xintercept=as.POSIXct("2009-05-19 12:00"), linetype="dashed", size=0.5, color="red")
  if(!is.null(xintercept)) gg <- gg + annotate("rect", xmin=xint.min, xmax=xint.max, ymin=0, ymax=Inf, fill="blue", alpha=0.2)
  if(plot){
    if(is.na(file)){
      plot(gg)
    }else{
      ggsave(file, gg)
    }
  }else{
    return(gg)
  }
}

transform.timedep.par.sample <- function(sample.in, s.cnst, sudriv, mnprm, scaleshift){
  if(!all(names(sample.in) == names(sudriv$model$parameters)[as.logical(sudriv$model$timedep$pTimedep)])) stop("order of timdependent parameters is wrong.")
  if(!all(is.na(mnprm))){ ## shift mean of some timedep-parameters for which the mean was re-parameterized as a constant parameter
    td <- names(sample.in)
    ind.mnprm <- which(td %in% gsub("_fmean","",mnprm))
    if(!all(gsub("_fmean","",mnprm) %in% td)) stop("some parameters of ", mnprm, " not found")
    #if(sum(ind.mnprm)!=1) stop(paste0("too many or too few parameters of ", mnprm, " found"))
    for(i in 1:length(ind.mnprm)){
      mn.curr <- td[ind.mnprm[i]]
      mn <- s.cnst[,paste0(mn.curr,"_fmean")]
      tran <- sudriv$model$args$parTran[names(sudriv$model$parameters)==mn.curr] == 1 | (mn.curr %in% gsub("_fmean","",rownames(scaleshift)))
      if(tran){
        sample.in[[ind.mnprm[i]]][-1,] <- sample.in[[ind.mnprm[i]]][-1,] + mn
      }else{
        sample.in[[ind.mnprm[i]]][-1,] <- sample.in[[ind.mnprm[i]]][-1,] * mn
      }
    }
  }
  if(!all(is.na(scaleshift))){ ## back-transform parameter with sigmoid transformation
    if(nrow(scaleshift)!=length(sample.in) | ncol(scaleshift)!=2) stop("dimension of scaleshift is not right")
    for(i in 1:length(sample.in)){
      sample.in[[i]][-1,] <- sigm.trans(sample.in[[i]][-1,], scale=scaleshift[i,1], shift=scaleshift[i,2])
    }
  }
  ## force time course within bounds after addition or multiplication with fmean parameter
  lo <- sudriv$model$args$parLo[sudriv$model$timedep$pTimedep]
  hi <- sudriv$model$args$parHi[sudriv$model$timedep$pTimedep]
  for(i in 1:length(sample.in)){
    sample.in[[i]][-1,] <- pmin(pmax(sample.in[[i]][-1,], lo[i]), hi[i])
  }
  return(sample.in)
}

analyze.clones <- function(timedep1, timedep2, brn.in1=1, brn.in2=1){
  nm.tdpar <- names(timedep1$sample.param.timedep)
  tdsamp1 <- timedep1$sample.param.timedep[[nm.tdpar]]
  tdsamp1 <- tdsamp1[c(1,(brn.in1+1):nrow(tdsamp1)),]
  tdsamp2 <- timedep2$sample.param.timedep[[nm.tdpar]]
  tdsamp2 <- tdsamp2[c(1,(brn.in2+1):nrow(tdsamp2)),]
  logpdf1 <- timedep1$sample.logpdf[(brn.in1+1):nrow(timedep1$sample.logpdf),]
  logpdf2 <- timedep2$sample.logpdf[(brn.in2+1):nrow(timedep2$sample.logpdf),]
  #prepare constpar samples
  cnst1 <- timedep1$sample.param.const[(brn.in1+1):nrow(timedep1$sample.param.const),]
  cnst2 <- timedep2$sample.param.const[(brn.in2+1):nrow(timedep2$sample.param.const),]
  cnst1 <- cnst1 %>% as_tibble %>% pivot_longer(everything(), names_to="par", values_to="value") %>% mutate(clone=1)
  cnst2 <- cnst2 %>% as_tibble %>% pivot_longer(everything(), names_to="par", values_to="value") %>% mutate(clone=2)
  cnst.data <- rbind(cnst1,cnst2) %>% mutate(clone=as.factor(clone))
  
  #transfrom samples
  fmean1 <- timedep1$sample.param.const[,grepl("fmean", colnames(timedep1$sample.param.const))]
  fmean1 <- fmean1[(brn.in1):length(fmean1)]
  fmean2 <- timedep2$sample.param.const[,grepl("fmean", colnames(timedep2$sample.param.const))]
  fmean2 <- fmean2[(brn.in2):length(fmean2)]
  tran1 <- exp(fmean1 + tdsamp1[2:nrow(tdsamp1),])
  tran2 <- exp(fmean2 + tdsamp2[2:nrow(tdsamp2),])
  rng1 <- apply(tran1, 2, range)
  rng2 <- apply(tran2, 2, range)
  #make data frame for time course
  bounds <- rbind(cbind(time=tdsamp1[1,], t(rng1), clone=1),cbind(time=tdsamp2[1,], t(rng2), clone=2))
  colnames(bounds) <- c("time", "lower", "upper", "clone")
  bounds <- bounds %>% as.data.frame %>% mutate(clone=as.factor(clone))##%>% pivot_longer(cols=c(lower,upper), names_to="type", values_to="value")
  #make data frame for logpost
  dat.logpost <- as.data.frame(rbind(cbind(logpdf1,clone=1), cbind(logpdf2, clone=2))) %>% mutate(clone=as.factor(clone))
  # plot time courses
  plt1 <- ggplot(data=bounds) + geom_ribbon(mapping=aes(x=time, ymin=lower, ymax=upper, fill=clone, group=clone), alpha=0.6)
  #ggsave(filename = "compare_clones.pdf", plt, width=8, height=5)
  # plot logposteriors
  plt2 <- ggplot(data=dat.logpost) + stat_density(mapping=aes(x=logposterior, fill=clone, color=clone, group=clone), alpha=0.6, trim=TRUE)
  plt3 <- ggplot(data=dat.logpost) + stat_density(mapping=aes(x=loglikeliobs, fill=clone, color=clone, group=clone), alpha=0.6, trim=TRUE)
  plt <- egg::ggarrange(plt1, plt2, plt3, ncol=1)
  ggsave(filename = "compare_clones.pdf", plt, width=8, height=5)
  
  ##plot posterior marginal for parameters
  gglist <- list()
  for(par.curr in unique(cnst.data$par)){
    dat.curr <- cnst.data %>% filter(par==par.curr)
    gg <- ggplot(data=dat.curr) + stat_density(mapping=aes(x=value, fill=clone, color=clone, group=clone), alpha=0.6, trim=TRUE)+
      labs(title =par.curr) + theme(legend.position = "none")
    gglist <- c(gglist, list(gg))
  }
  gg <- egg::ggarrange(plots=gglist, ncol=4, draw = FALSE)
  ggsave(filename = "compare_parameters.pdf", gg, width=10, height=10)
}

## Translate variable names for plotting
mylabeller.param <- function(labs){
  x <- c(dsplsd=expression(D),
         smaxur=expression(S[u*",max"]),
         beqqur=expression(beta[u]),
         kqqsr2=expression(k[g]),
         kpqqfr=expression(k[i]),
         pmaxed=expression(P[ex]),
         smaxir=expression(S[t*",max"]),
         kqqrr=expression(k[c]),
         cmlte=expression(phi[e]),
         alqqsr=expression(alpha[g]),
         cmltp=expression(phi[p]),
         kqqfr=expression(k[d]),
         kdwr=expression(lambda),
         rswr=expression(r[s]),
         sloneir=expression(S[t*",z1"]),
         sltwoir=expression(S[t*",z2"]),
         alqqfr=expression(alpha[d]))
  return(x[labs])
}

mylabeller.param.units <- function(labs, log=FALSE){
  x <- c(dsplsd=expression(D~"(-)"),
         smaxur=expression(S[u*",max"]~"(mm)"),
         beqqur=expression(beta[u]~"(-)"),
         kqqsr2=ifelse(log, expression(ln~k[g]^"*"~"(-)"),
                       expression(k[g]~"("*mm^{1-alpha[g]}~h^{-1}*")")),
         kpqqfr=expression(k[i]~"("*h^{-1}*")"),
         pmaxed=expression(P[ex]~"(mm"*h^{-1}*")"),
         smaxir=expression(S[t*",max"]~"(mm)"),
         kqqrr=expression(k[c]~"("*h^{-1}*")"),
         cmlte=expression(phi[e]~"(-)"),
         alqqsr=ifelse(log, expression(ln~alpha[g]^"*"~"(-)"), expression(alpha[g]~"(-)")),
         cmltp=expression(phi[p]~"(-)"),
         kqqfr=ifelse(log, expression(ln~k[d]^"*"~"(-)"),
                      expression(k[d]~"("*mm^{1-alpha[d]}~h^{-1}*")")),
         kdwr=expression(lambda~"("*h^{-1}*")"),
         rswr=expression(r[s]~"("*h^{-1}*")"),
         sloneir=expression(S[t*",z1"]~"(mm)"),
         sltwoir=ifelse(log, expression(ln~S[t*",z2"]^"*"~"(-)"), expression(S[t*",z2"]~"(mm)")),
         alqqfr=expression(alpha[d]~"(-)"))
  return(x[labs])
}

mylabeller.feat <- function(labs){
  x <- c(U5F1Wv_Ss1=expression(S[g]),
         U3F1Wv_Su1=expression(S[u]),
         C1Wv_Qstream="Streamflow",
         prec="Precip.",
         temp="Temp.",
         U3F1Wv_Si1=expression(S[t]),
         epot=expression(E[pot]),
         U1F1Wv_Sf1=expression(S[d]),
         Luft.Feuchte="Humidity",
         Wind.v=expression(v[Wind]))
  return(x[labs])
}

mylabeller.feat.units <- function(labs){
  x <- c(U5F1Wv_Ss1=expression(S[g]~"(mm)"),
         U3F1Wv_Su1=expression(S[u]~"(mm)"),
         C1Wv_Qstream="Streamflow",
         prec="Precip.",
         temp="Temp.",
         U3F1Wv_Si1=expression(S[t]),
         epot=expression(E[pot]),
         U1F1Wv_Sf1=expression(S[d]),
         Luft.Feuchte="Humidity",
         Wind.v=expression(v[Wind]))
  return(x[labs])
}

mylabeller.feat.many <- function(labs){
  mtch <- c(U5F1Wv_Ss1="S[g]",
            U3F1Wv_Su1="S[u]",
            C1Wv_Qstream="Streamflow",
            prec="Precip.",
            temp="Temp.",
            U3F1Wv_Si1="S[t]",
            epot="E[pot]",
            U1F1Wv_Sf1="S[d]",
            Luft.Feuchte="Humidity",
            Wind.v="v[Wind]")
  fn <- function(x,y){
    sapply(y, grepl, x=x)
  }
  mat <- mapply("fn", x=labs, MoreArgs=list(y=names(mtch)))
  mat <- apply(mat, 2, function(x) mtch[x])
  ex <- vector(mode="expression")
  if(is.list(mat)){
    for(i in 1:length(mat)){
      ex[i] <- parse(text=paste(mat[[i]], collapse="+"))
    }
  }else{
    for(i in 1:ncol(mat)){
      ex[i] <- parse(text=paste(mat[,i], collapse="+"))
    }
  }
  return(ex)
}
