remove.taumax.Q <- function(sudriv){
  ## removes the correlation in Q of the error model and removes taumax from the fitted parameters. This is manly needed because the time-dependent parameters replace the correlation in the error model.
  ind <- which(names(sudriv$likelihood$parameters)=="GLOB_Mult_Q_taumax_lik")
  sudriv$likelihood$tran[ind] <- 0
  sudriv$likelihood$lb[ind] <- 0
  sudriv$likelihood$ub[ind] <- 10000
  sudriv$likelihood$parameters["GLOB_Mult_Q_taumax_lik"] <- 0
  sudriv$likelihood$par.fit[ind] <- 0
  tmp <- dimnames(sudriv$parameter.sample)[[2]]=="GLOB_Mult_Q_taumax_lik"
  sudriv$parameter.sample <- sudriv$parameter.sample[,!tmp,]
  return(sudriv)
}
fix.taumax.Q <- function(sudriv){
  ## fixes taumax at the value that was inferred with the constant parameters
  ind <- which(names(sudriv$likelihood$parameters)=="GLOB_Mult_Q_taumax_lik")
  print("taumax is fixed at:")
  print(sudriv$likelihood$parameters[ind])
  sudriv$likelihood$par.fit[ind] <- 0
  tmp <- dimnames(sudriv$parameter.sample)[[2]]=="GLOB_Mult_Q_taumax_lik"
  sudriv$parameter.sample <- sudriv$parameter.sample[,!tmp,]
  return(sudriv)
}
fix.a.func <- function(sudriv){
  ## fixes a at the value that was inferred with the constant parameters
  ind <- which(names(sudriv$likelihood$parameters)=="C1Wv_Qstream_a_lik")
  print("a is fixed at:")
  print(sudriv$likelihood$parameters[ind])
  sudriv$likelihood$par.fit[ind] <- 0
  tmp <- dimnames(sudriv$parameter.sample)[[2]]=="C1Wv_Qstream_a_lik"
  sudriv$parameter.sample <- sudriv$parameter.sample[,!tmp,]
  
  ind <- which(names(sudriv$likelihood$parameters)=="C1Tc1_Qstream_a_lik")
  print("a is fixed at:")
  print(sudriv$likelihood$parameters[ind])
  sudriv$likelihood$par.fit[ind] <- 0
  tmp <- dimnames(sudriv$parameter.sample)[[2]]=="C1Tc1_Qstream_a_lik"
  sudriv$parameter.sample <- sudriv$parameter.sample[,!tmp,]
  
  ind <- which(names(sudriv$likelihood$parameters)=="C1Tc2_Qstream_a_lik")
  print("a is fixed at:")
  print(sudriv$likelihood$parameters[ind])
  sudriv$likelihood$par.fit[ind] <- 0
  tmp <- dimnames(sudriv$parameter.sample)[[2]]=="C1Tc2_Qstream_a_lik"
  sudriv$parameter.sample <- sudriv$parameter.sample[,!tmp,]
  
  return(sudriv)
}
release.b <- function(sudriv){
  ## makes the "b" of the residual error model of streamflow, atra and terb a fitted parameter
  sudriv$likelihood$par.fit[which(names(sudriv$likelihood$parameters)=="C1Wv_Qstream_b_lik")] <- 1
  sudriv$likelihood$par.fit[which(names(sudriv$likelihood$parameters)=="C1Tc1_Qstream_b_lik")] <- 1
  sudriv$likelihood$par.fit[which(names(sudriv$likelihood$parameters)=="C1Tc2_Qstream_b_lik")] <- 1
  ## invent a sample for "b", based on which the initial covariance of the jump distribution will be calculated by infer.timedeppar
  dims <- dim(sudriv$parameter.sample)
  b.samp.q <- rnorm(prod(dims[c(1,3)]), mean=-2, sd=0.3)
  b.samp.c <- rnorm(prod(dims[c(1,3)]), mean=-4, sd=0.3)
  tmp1 <- array(b.samp.q, dim=c(dims[1],1,dims[3]))
  tmp2 <- array(b.samp.q, dim=c(dims[1],1,dims[3]))
  # combine the invented sample with the existing one
  smp.new <- abind(sudriv$parameter.sample, tmp1, along=2)
  colnames(smp.new)[ncol(smp.new)] <- "C1Wv_Qstream_b_lik" #streamflow
  smp.new <- abind(smp.new, tmp2, along=2)
  colnames(smp.new)[ncol(smp.new)] <- "C1Tc1_Qstream_b_lik" # atrazine
  smp.new <- abind(smp.new, tmp2, along=2)
  colnames(smp.new)[ncol(smp.new)] <- "C1Tc2_Qstream_b_lik" # terbuthylazine
  
  sudriv$parameter.sample <- smp.new # replace
  return(sudriv)
}
redef.init.range <- function(sudriv, drop=0.9, jitter=0, init.range.orig=matrix(0,ncol=2)){
  if(jitter != 0 & is.na(init.range.orig[1])) warning("No init.range.orig applied in case of jitter")
  logpost <- quantile(sudriv$posterior.sample[nrow(sudriv$posterior.sample),], drop)
  s <- remove.chains(sudriv, brn.in=nrow(sudriv$posterior.sample)-2, logpost=logpost)$sample
  init.range <- t(apply(s[nrow(s),,], 1, range))
  lower.new <- init.range[,1] - jitter*(init.range[,1]-init.range.orig[,1])
  init.range[,2] <- init.range[,2] + jitter*(init.range.orig[,2]-init.range[,2])
  init.range[,1] <- lower.new
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
select.maxlikpars.optimized <- function(sudriv, optimized){
  ## select a parameter set of the multiple optimizations and add the optimizations to the sudriv object
  objective <- Inf
  sudriv$optimized <- list()
  j <- 1
  for(optim.curr in optimized){
    sudriv$optimized[[j]] <- list(status = optim.curr$status, iterations = optim.curr$iterations, objective = optim.curr$objective, solution = optim.curr$solution, message = optim.curr$message)
    if(optim.curr$objective < objective){
      objective <- optim.curr$objective
      solution  <- optim.curr$solution
    }
    j <- j + 1
  }
  sudriv$model$parameters[as.logical(sudriv$model$par.fit)] <- solution[1:sum(sudriv$model$par.fit)]
  sudriv$likelihood$parameters[as.logical(sudriv$likelihood$par.fit)] <- solution[(sum(sudriv$model$par.fit)+1):length(solution)]
  return(sudriv)
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
select.ind <- function(sudriv, xlim, ind.sel, calibpred="calib"){
  ## select index of layout that satisfies condition
  ## create data frame for ggplot-object
  L <- ifelse(calibpred=="calib", "layout", "pred.layout")
  if(is.na(xlim[1])) xlim <- c(-Inf, Inf)
  if("POSIXct" %in% class(xlim[1])){
    tme <- as.POSIXct(sudriv$layout$tme.orig) + sudriv$layout[[L]]$time*60*60*ifelse(sudriv$layout$time.units=="days",24,1)
    ind.sel <- which(tme >= xlim[1] & tme <= xlim[2])
    return(ind.sel)
  }else{
    if(xlim[1]=="pred") xlim <- range(sudriv$layout$pred.layout$time)
    if(xlim[1]=="calib") xlim <- range(sudriv$layout$layout$time[sudriv$layout$calib])
    if(is.na(ind.sel[1])){
      ind.sel <- which(sudriv$layout[[L]]$time >= xlim[1] & sudriv$layout[[L]]$time <= xlim[2])
    }else{
      ind.sel <- ind.sel[sudriv$layout[[L]]$time[ind.sel] >= xlim[1] & sudriv$layout[[L]]$time[ind.sel] <= xlim[2]]
    }
    return(ind.sel)
  }
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
get.loess.input <- function(sudriv, tag, vars, res=NULL, add.data=NULL, t.lim=NULL, remove.na=TRUE, with.td=TRUE){
  ## This function extracts the data needed to fit some linear and nonlinear models to the time-course of the time dependent parameter.
  ## if with.td=TRUE and 'res' is supplied, the time series of the parameter is taken from res. If 'res' is not supplied, the time series is taken form sudriv
  if(with.td){
    if(is.null(res)){
      if(is.null(sudriv$model$timedep)) stop("function 'get.loess.input' requires non-null sudriv$model$timedep")
      if(dim(sudriv$model$timedep$par)[2]>1) warning("'get.loess.input' is not (yet) implemented for multiple timedependent parameters")
      y.timedep <- c(sudriv$model$timedep$par)
    }else{
      if(length(res$sample.param.timedep)!=1) stop("none or multiple timedependent parameters found. This function cannot cope with that.")
      pm <- which.max(res$sample.logpdf[,"loglikeliobs"])
      y.timedep <- res$sample.param.timedep[[1]][-1,][pm,]
    }
  }
  if(!is.null(vars)){
    layout.states <- list(layout = data.frame(var=rep(vars, each=nrow(sudriv$input$inputobs)), time=rep(sudriv$input$inputobs[,1], length(vars)), stringsAsFactors=FALSE),
                          lump   = rep(NA, nrow(sudriv$input$inputobs)*length(vars)))
    y.all <- run.model(layout=layout.states, sudriv=sudriv, lump=FALSE)$original
    y.all <- cbind(layout.states$layout, y.all)
    y.all <- y.all %>% spread(var, y.all)
  }else{
    y.all <- data.frame(nothing99=rep(NA,nrow(sudriv$input$inputobs))) ## initialize y.all without model output
  }
  ## get the states to compare it to
  y.all <- y.all %>% mutate(prec = pmax(sudriv$input$inputobs[,"P"],0), prec=rollmean(prec, k=5, fill=0)) # moving average of precipitation
  y.all <- y.all %>% mutate(epot = pmax(sudriv$input$inputobs[,"Epot"],0))
  ##y.all <- y.all %>% mutate(temp = pmax(sudriv$input$inputobs[,"T"],0))
  if(with.td) y.all <- y.all %>% mutate(y.td = y.timedep)
  if("nothing99" %in% colnames(y.all)) y.all <- y.all %>% select(-nothing99)
  ## add the additional data in function argument
  if(!is.null(add.data)){
    ## assuming the time column of add.data is named time, transform it to time of sudriv object
    add.data <- add.data %>% mutate(time= as.numeric((time - as.POSIXct(sudriv$layout$tme.orig)))*ifelse(sudriv$layout$time.units=="days",1,24))
    ## interpolate it to existing data
    add.data <- apply(X=add.data%>%select(-time),2,FUN=function(y) approx(x=add.data$time, y=y, xout=y.all$time)$y)
    y.all <- cbind(y.all, add.data)
  }
  
  ## consistency check
  if(with.td) if(length(y.timedep) != nrow(y.all)) stop("dimension mismatch")
  ## lm1 <- lm(y.td ~ ., data=y.all%>%select(-time))
  y.all2 <- y.all
  ## limit the analysis to the period where we actually have data...
  tag.red <- gsub("_.*","",tag)
  if(is.null(t.lim)){
    if(tag.red %in% c("kdwr","rswr","sloneir","sltwoir")){# if it is a chemistry related parameter
      strt <- sudriv$layout$layout %>% slice(sudriv$layout$calib) %>% filter(var %in% c("C1Tc1_Qstream","C1Tc2_Qstream")) %>% select(time) %>% min
      end <- sudriv$layout$layout %>% slice(sudriv$layout$calib) %>% filter(var %in% c("C1Tc1_Qstream","C1Tc2_Qstream")) %>% select(time) %>% max
    }else{ #if it is a more water related parameter
      strt <- min(sudriv$layout$layout$time[sudriv$layout$calib])
      end <- max(sudriv$layout$layout$time[sudriv$layout$calib])
    }
  }else{
    strt <- t.lim[1]
    end <- t.lim[2]
  }
  cat("strt: ",strt,"\n")
  cat("end: ",end,"\n")
  y.all2 <- y.all2 %>% filter(time >= strt & time <= end)
  if(remove.na) y.all2 <- y.all2 %>% na.omit
  cat("dim data:\t",dim(y.all2),"\n")
  if(with.td){
    if(is.null(res) & sudriv$model$args$parTran[which(sudriv$model$timedep$pTimedep)[1]] == 1){
      ## in this case we want the time course in the original units, not the transformed one (e.g. for plotting)
      y.all2 <- y.all2 %>% mutate(y.td=exp(y.td))
    }
  }
  return(y.all2)
}
remove.constpar.su <- function(sudriv, param){
  if(length(param)!=1) stop("'param' must be of length 1")
  ## removes an arbitrary parameter from the fitted model parameters of a sudriv object
  ind <- which(names(sudriv$model$parameters)==param)
  sudriv$model$par.fit[ind] <- 0
  tmp <- dimnames(sudriv$parameter.sample)[[2]]==param
  sudriv$parameter.sample <- sudriv$parameter.sample[,!tmp,]
  return(sudriv)
}