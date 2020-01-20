s.m.mcmc.wrapper <- function(log.posterior, max.iter, sudriv, savepath, tag, drop=0, init.range=NULL, init.state=NULL, iter.step=30000, jitter=0, ...){
  done.iter <- 0
  init.range.orig <- init.range
  warm.up <- 0
  while(done.iter < max.iter){
    if((done.iter != 0 | warm.up != 0) & drop > 0) init.range <- redef.init.range(sudriv, drop=drop, jitter=ifelse(done.iter<(max.iter-2*iter.step),jitter,0), init.range.orig=init.range.orig)
    if(warm.up < 3){
      result.s.m = my.s.m.mcmc(log.posterior, max.iter=600, init.range=init.range, init.state=init.state, sudriv=sudriv, ...)
      warm.up <- warm.up + 1
      if(all(!is.finite(c(result.s.m$log.p)))) stop("all logposteriors are non-numeric")
    }else{
      result.s.m = my.s.m.mcmc(log.posterior, max.iter=iter.step, init.range=init.range, init.state=init.state, sudriv=sudriv, ...)
      done.iter <- done.iter + iter.step
    }
    s <- result.s.m$samples
    init.state <- s[,dim(s)[2],]
    if(drop > 0) init.state <- NULL
    sudriv$parameter.sample <- aperm(s, perm=c(2,3,1))
    colnames(sudriv$parameter.sample) <- c(names(sudriv$model$parameters)[as.logical(sudriv$model$par.fit)], names(sudriv$likelihood$parameters)[as.logical(sudriv$likelihood$par.fit)])##, paste("mu", 1:(1*n.mu), sep=""))
    sudriv$posterior.sample <- t(result.s.m$log.p)
    su <- sudriv
    dir.create(savepath)
    save(su, file=paste(savepath, "/su_", tag, ".RData", sep=""), version=2)
  }
  return(sudriv)
}
my.s.m.mcmc <- function (f, max.iter, n.walkers, n.dim, init.range, init.state=NULL, ...){
  chain.length = max.iter%/%n.walkers
  log.p = matrix(NA, nrow = n.walkers, ncol = chain.length)
  log.p.old = rep(NA, n.walkers)
  ensemble.old = matrix(NA, nrow = n.walkers, ncol = n.dim)
  ensemble.new = matrix(NA, nrow = n.walkers, ncol = n.dim)
  samples = array(NA, dim = c(n.walkers, chain.length, n.dim))
  mcmc.object = array(NA, dim = c(n.walkers, chain.length,n.dim + 1))
  if(is.null(init.state)){
    for (k in 1:n.walkers) {
      for (g in 1:n.dim) {
        ensemble.old[k, g] = runif(1, init.range[g, 1], init.range[g,
                                                                   2])
      }
      log.p.old[k] = f(ensemble.old[k, ], ...)
    }
  }else{
    ensemble.old <- init.state
    for(k in 1:n.walkers){
      log.p.old[k] <- f(init.state[k, ], ...)
    }
  }
  log.p[, 1] = log.p.old
  samples[, 1, ] = ensemble.old
  for (l in 2:chain.length) {
    for (n in 1:n.walkers) {
      z = ((runif(1) + 1)^2)/2
      a = sample((1:n.walkers)[-n], 1)
      par.active = ensemble.old[a, ]
      ensemble.new[n, ] = par.active + z * (ensemble.old[n,
                                                         ] - par.active)
      log.p.new = f(ensemble.new[n, ], ...)
      if (!is.finite(log.p.new)) {
        acc = 0
      }
      else {
        acc = z^(n.dim - 1) * exp(log.p.new - log.p.old[n])
      }
      test = runif(1)
      if (acc > test) {
        samples[n, l, ] = ensemble.new[n, ]
        ensemble.old[n, ] = ensemble.new[n, ]
        log.p[n, l] = log.p.new
        log.p.old[n] = log.p.new
      }
      else {
        samples[n, l, ] = ensemble.old[n, ]
        log.p[n, l] = log.p.old[n]
      }
    }
  }
  mcmc.list = list(samples = samples, log.p = log.p)
  return(mcmc.list)
}