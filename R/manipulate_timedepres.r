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
