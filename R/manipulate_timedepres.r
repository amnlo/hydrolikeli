remove.constpar <- function(res, param){
  ## this function takes a result of 'infer.timedeppar' of the timedeppar package, removes the parameter 'param' from the constant parameters that are inferred and returns the new result object.
  ##if(res$package != "timedeppar 0.5 2019-10-18") warning("not the version of result I expected...")
  res.new <- res
  for(p.curr in param){
    if(!(p.curr %in% names(res$param.ini))){warning(paste0("parameter ",p.curr," not found"));next}
    res.new$param.ini[[p.curr]] <- NULL
    res.new$param.range[[p.curr]] <- NULL
    res.new$sample.param.const <- res.new$sample.param.const[,colnames(res.new$sample.param.const)!=p.curr]
    res.new$param.maxpost[[p.curr]] <- NULL
    rmcnst1 <- rownames(res.new$cov.prop.const)==p.curr
    rmcnst2 <- colnames(res.new$cov.prop.const)==p.curr
    res.new$cov.prop.const <- res.new$cov.prop.const[!rmcnst1,!rmcnst2]
  }
  return(res.new)
}
fix.a.q.restart <- function(res){
  sudriv <- res$dot.args$sudriv
  ## fixes a of Q at the value that was inferred with the previous run with time-dependent parameters
  ind <- which(names(sudriv$likelihood$parameters)=="C1Wv_Qstream_a_lik")
  sudriv$likelihood$parameters["C1Wv_Qstream_a_lik"] <- res$sample.param.const[nrow(res$sample.param.const),"C1Wv_Qstream_a_lik"]
  print("a is fixed at:")
  print(sudriv$likelihood$parameters[ind])
  sudriv$likelihood$par.fit[ind] <- 0
  tmp <- dimnames(sudriv$parameter.sample)[[2]]=="C1Wv_Qstream_a_lik"
  sudriv$parameter.sample <- sudriv$parameter.sample[,!tmp,]
  res$dot.args$sudriv <- sudriv
  
  res$param.ini[["C1Wv_Qstream_a_lik"]] <- NULL
  res$sample.param.const <- res$sample.param.const[,colnames(res$sample.param.const) != "C1Wv_Qstream_a_lik"]
  res$param.maxpost[["C1Wv_Qstream_a_lik"]] <- NULL
  res$cov.prop.const <- res$cov.prop.const[rownames(res$cov.prop.const) != "C1Wv_Qstream_a_lik", colnames(res$cov.prop.const) != "C1Wv_Qstream_a_lik"]
  return(res)
}
release.a.q.restart <- function(res, ini.a=NULL, sd.a=0.1){
  ## This function relaxes the fixation of "a" of streamflow of a result of a timedeppar inference
  ## "sd.a" is the standard deviation of the jump distribution in the dimension of "a"
  sudriv <- res$dot.args$sudriv
  ind <- which(names(sudriv$likelihood$parameters)=="C1Wv_Qstream_a_lik")
  ## set the initial "a" if not provided in arguments
  if(is.null(ini.a)){
    ini.a <- sudriv$likelihood$parameters["C1Wv_Qstream_a_lik"]
  }
  ## Add "a" to fitted parameters of sudriv object
  sudriv$likelihood$par.fit[ind] <- 1
  res$dot.args$su <- sudriv
  
  print("a is freed at initial value:")
  print(ini.a)

  ## Insert initial value for "a"
  insrt.before <- which(names(res$param.ini)=="C1Tc1_Qstream_a_lik")
  if("C1Wv_Qstream_a_lik" %in% names(res$param.ini)) stop("trying to release a, but it is already fitted")
  res$param.ini <- c(res$param.ini[1:(insrt.before-1)], list("C1Wv_Qstream_a_lik"=ini.a), res$param.ini[insrt.before:length(res$param.ini)])

  ## Insert (fake) parameter sample for "a"
  insrt.before <- which(colnames(res$sample.param.const)=="C1Tc1_Qstream_a_lik")
  res$sample.param.const <- cbind(res$sample.param.const[,1:(insrt.before-1)], "C1Wv_Qstream_a_lik"=rnorm(nrow(res$sample.param.const),ini.a,sd.a), res$sample.param.const[,insrt.before:ncol(res$sample.param.const)])

  ## Insert (fake) parameter maxpost for "a"
  insrt.before <- which(names(res$param.maxpost)=="C1Tc1_Qstream_a_lik")
  res$param.maxpost <- c(res$param.maxpost[1:(insrt.before-1)], list("C1Wv_Qstream_a_lik"=ini.a), res$param.maxpost[insrt.before:length(res$param.maxpost)])

  ## Insert row/col of proposal covariance matrix for "a"
  insrt.before <- which(colnames(res$cov.prop.const)=="C1Tc1_Qstream_a_lik")
  res$cov.prop.const <- cbind(res$cov.prop.const[,1:(insrt.before-1)], "C1Wv_Qstream_a_lik"=0, res$cov.prop.const[,insrt.before:ncol(res$cov.prop.const)])
  res$cov.prop.const <- rbind(res$cov.prop.const[1:(insrt.before-1),], "C1Wv_Qstream_a_lik"=0, res$cov.prop.const[insrt.before:nrow(res$cov.prop.const),])
  res$cov.prop.const[insrt.before,insrt.before] <- sd.a^2
  
  return(res)
}
release.cmltp.restart <- function(res, ini.cmltp=NULL, sd.cmltp=0.1){
  ## This function relaxes the fixation of the multiplier parameter for precipitation of a result of a timedeppar inference
  ## "sd.cmltp" is the standard deviation of the jump distribution in the dimension of "cmltp"
  sudriv <- res$dot.args$sudriv
  ind <- which(names(sudriv$model$parameters)=="Glo%Cmlt_P")
  ## set the initial "cmltp" if not provided in arguments
  if(is.null(ini.cmltp)){
    ini.cmltp <- sudriv$model$parameters["Glo%Cmlt_P"]
  }
  ## Add "cmltp" to fitted parameters of sudriv object
  sudriv$model$par.fit[ind] <- 1
  ## Adapt the boundary to match the boundary in the timedep case
  sudriv$model$args$parLo[ind] <- 0.9
  sudriv$model$args$parHi[ind] <- 1.4
  res$dot.args$su <- sudriv

  print("cmltp is freed at initial value:")
  print(ini.cmltp)

  ## Insert initial value for "cmltp" (this is neglected in the current version of timedeppar)
  if("Glo%Cmlt_P" %in% names(res$sample.param.const)) stop("trying to release cmltp, but it is already fitted")
  res$param.ini <- c(list("Glo%Cmlt_P"=ini.cmltp), res$param.ini)
  print(res$param.ini[sapply(res$param.ini, length)==1])
  
  ## Insert (fake) parameter sample for "cmltp"
  res$sample.param.const <- cbind("Glo%Cmlt_P"=rnorm(nrow(res$sample.param.const),ini.cmltp,sd.cmltp), res$sample.param.const)
  
  ## Insert (fake) parameter maxpost for "cmltp"
  res$param.maxpost <- c(list("Glo%Cmlt_P"=ini.cmltp), res$param.maxpost)
  
  ## Insert row/col of proposal covariance matrix for "cmltp"
  res$cov.prop.const <- cbind("Glo%Cmlt_P"=0, res$cov.prop.const)
  res$cov.prop.const <- rbind("Glo%Cmlt_P"=0, res$cov.prop.const)
  res$cov.prop.const[1,1] <- sd.cmltp^2
  
  return(res)
}
remove.timedeppar <- function(res){
  ## This function removes all timedependent parameters from the inference result. It was
  ## created so that we can continue the inference from a result, not inferring the time series, but calculating it based on
  ## an empirical internal relationship (hybrid model)
  nm.td <- names(res$param.maxpost)[which(sapply(res$param.maxpost, length)>1)]
  for(i in nm.td){
    res$param.ini[[i]] <- NULL#mean(res$param.ini[[i]][,2])
  }
  for(i in nm.td){
    res$param.maxpost[[i]] <- NULL#mean(res$param.maxpost[[i]][,2])
  }
  res$param.ou.ini <- numeric()
  res$param.ou.fixed <- numeric()
  res$sample.param.timedep <- list()
  res$param.ou.maxpost <- numeric()
  res$cov.prop.ou <- NA

  return(res)
}
release.sd.ou.restart <- function(res, scale.prop.ou.ini=0.3){
  ## this function makes the standard deviation of the OU process a fitted parameter (if it was fixed before)
  nme <- grep("_sd", names(res$param.ou.fixed), value=TRUE)
  sd.ou <- res$param.ou.fixed[nme]
  if(length(sd.ou)>1) stop("cannot deal with multiple timedependent parameters")
  res$param.ou.fixed <- res$param.ou.fixed[names(res$param.ou.fixed)!=nme] # remove it from the fixed parameters
  res$param.ou.ini <- sd.ou
  res$sample.param.ou <- matrix(as.numeric(sd.ou), nrow=nrow(res$sample.param.const), ncol=1)
  colnames(res$sample.param.ou) <- nme
  res$param.ou.maxpost <- sd.ou
  mat <- matrix(as.numeric(sd.ou^2), nrow=1, ncol=1)
  colnames(mat) <- nme
  rownames(mat) <- nme
  res$cov.prop.ou <- list(mat)
  res$scale.prop.ou <- scale.prop.ou.ini
  return(res)
}
