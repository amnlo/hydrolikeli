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