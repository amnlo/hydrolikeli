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
get.loess.input <- function(sudriv, tag, vars, add.data, t.lim=NULL, remove.na=TRUE){
  ## This function extracts the data needed to fit some linear and nonlinear models to the time-course of the time dependent parameter.
  if(is.null(sudriv$model$timedep)) stop("function 'find.pattern.timedep' requires non-null sudriv$model$timedep")
  if(dim(sudriv$model$timedep$par)[2]>1) warning("'find.pattern.timedep' is not (yet) implemented for multiple timedependent parameters")
  nm.td <- names(sudriv$model$parameters)[sudriv$model$timedep$pTimedep]
  y.timedep <- c(sudriv$model$timedep$par)
  if(!is.null(vars)){
    layout.states <- list(layout = data.frame(var=rep(vars, each=nrow(sudriv$input$inputobs)), time=rep(sudriv$input$inputobs[,1], length(vars)), stringsAsFactors=FALSE),
                          lump   = rep(NA, nrow(sudriv$input$inputobs)*length(vars)))
    y.mod <- run.model(layout=layout.states, sudriv=sudriv, lump=FALSE)$original
    y.mod <- cbind(layout.states$layout, y.mod)
    y.mod <- y.mod %>% spread(var, y.mod)
  }else{
    y.mod <- data.frame(nothing99=rep(NA,nrow(sudriv$input$inputobs))) ## initialize y.mod without model output
  }
  ## get the states to compare it to
  ## if(is.null(layout.states)) layout.states <- list(layout=sudriv$layout$layout, lump=rep(NA, nrow(sudriv$layout$layout)))
  y.mod <- y.mod %>% mutate(prec = pmax(rollmean(sudriv$input$inputobs[,"P"], k=10*24, na.pad=TRUE),0))
  y.mod <- y.mod %>% mutate(epot = pmax(rollmean(sudriv$input$inputobs[,"Epot"], k=10*24, na.pad=TRUE),0))
  y.mod <- y.mod %>% mutate(temp = pmax(rollmean(sudriv$input$inputobs[,"T"], k=10*24, na.pad=TRUE),0))
  y.all <- y.mod %>% mutate(y.td = y.timedep)
  if("nothing99" %in% colnames(y.mod)) y.mod <- y.mod %>% select(-nothing99)
  ## add the additional data in function argument
  if(!is.null(add.data)){
    ## assuming the time column of add.data is named time, transform it to time of sudriv object
    add.data <- add.data %>% mutate(time= as.numeric((time - as.POSIXct(sudriv$layout$tme.orig)))*ifelse(sudriv$layout$time.units=="days",1,24))
    ## interpolate it to existing data
    add.data <- apply(X=add.data%>%select(-time),2,FUN=function(y) approx(x=add.data$time, y=y, xout=y.all$time)$y)
    y.all <- cbind(y.all, add.data)
  }
  
  ## consistency check
  if(length(y.timedep) != nrow(y.all)) stop("dimension mismatch")
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
  if(sudriv$model$args$parTran[which(sudriv$model$timedep$pTimedep)[1]] == 1){
    y.all2 <- y.all2 %>% mutate(y.td=exp(y.td))
  }
  return(y.all2)
}