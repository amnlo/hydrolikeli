remove.constpar <- function(res, param){
  ## this function takes a result of 'infer.timedeppar' of the timedeppar package, removes the parameter 'param' from the constant parameters that are inferred and returns the new result object.
  ##if(res$package != "timedeppar 0.5 2019-10-18") warning("not the version of result I expected...")
  res.new <- res
  for(p.curr in param){
    if(!(p.curr %in% names(res$param.ini))){warning(paste0("parameter ",pucrr," not found"));next}
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
  sudriv$likelihood$par.fit[ind] <- 1
  print("a is freed at initial value:")
  print(sudriv$likelihood$parameters[ind])

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
  
  ## Add "a" to fitted parameters of sudriv object
  su <- res$dot.args$su
  su$likelihood$par.fit[names(su$likelihood$parameters)=="C1Wv_Qstream_a_lik"] <- 1
  res$dot.args$su <- su
  
  return(res)
}
