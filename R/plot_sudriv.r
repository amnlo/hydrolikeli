gg_color_hue <- function(n,start) {
  hues = seq(start, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
plot.markov.hist <- function(sudriv, brn.in = 0, n=1e4, pridef = NULL, v.line=NULL, lower.logpost=NA, prior.only=FALSE, plot=TRUE, kl.div=TRUE, file.hist=NA, width=9, height=7, file.kl=NA, tag=NULL, scl="posterior", lab=""){
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
  ## account for the fact that df = df - 2
  a.re[grepl("_df_lik", a.re$param),"value"] <- a.re[grepl("_df_lik", a.re$param),"value"] + 2
  ## scale parameters to the desired time format
  if(!is.null(sudriv$model$par.time)){
    ## make sure the global multipliers have the correct info regarding time proportionality (hard-coded)
    nm <- names(su$model$parameters)
    sudriv$model$par.time[grepl("Pmax_ED", nm)] <- -1
    sudriv$model$par.time[grepl("K[p|_Q|d]", nm)] <- -1
    sudriv$model$par.time[grepl("Rs", nm)] <- -1
    sudriv$model$par.time[grepl("Start", nm)] <- -1
    sudriv$model$par.time[grepl("tau", nm)] <- 1
    n.prop <- c(names(sudriv$model$parameters)[sudriv$model$par.time==1], names(sudriv$likelihood$parameters)[sudriv$likelihood$par.time==1])
    n.invprop <- c(names(sudriv$model$parameters)[sudriv$model$par.time==-1], names(sudriv$likelihood$parameters)[sudriv$likelihood$par.time==-1])
    for(i in par.names){
      if(i %in% n.prop) a.re[a.re$param==i,"value"] <- a.re[a.re$param==i,"value"]*sudriv$layout$timestep.fac
      if(i %in% n.invprop) a.re[a.re$param==i,"value"] <- a.re[a.re$param==i,"value"]/sudriv$layout$timestep.fac
    }
  }
  ## insert prior into data frame to be plotted:
  if(!is.null(pridef)){
    pri.samp <- matrix(NA, nrow=n, ncol=length(par.names))
    colnames(pri.samp) <- par.names
    l.pri <- 1000
    a.pri <- data.frame(value=rep(NA, l.pri*length(par.names)), param=NA, walker=NA, y=NA, pri=1)
    j <- 1
    for(par.curr in par.names){
      uni <- FALSE
      cut.up <- FALSE
      if(grepl("_a_lik", par.curr) & grepl("lognorm",pridef[[par.curr]][1])){
        a.m <- as.numeric(pridef[[par.curr]][2])
        a.sd <- as.numeric(pridef[[par.curr]][3])
        a.sd2 <- sqrt(log(1+a.sd^2/a.m^2))
        a.m2 <- log(a.m) - 0.5*a.sd2^2
        a.m2 <- a.m2 - log(5)
        pridef[[par.curr]][2] <- as.character(exp(a.m2+a.sd2^2/2))
        pridef[[par.curr]][3] <- as.character(sqrt(exp(2*a.m2+a.sd2^2)*(exp(a.sd2^2)-1)))
      }
      if(which(par.curr==par.names) %in% which(as.logical(par.trans)) & pridef[[par.curr]][1] %in% c("normal", "Normal", "norm", "Norm", "normaltrunc")){
        cut.up <- TRUE
        if(grepl("trunc", pridef[[par.curr]][1])){
          pridef[[par.curr]][1] <- "lognormaltrunc"
          pridef[[par.curr]][4:5] <- as.character(exp(as.numeric(pridef[[par.curr]][4:5])))
          if(par.curr %in% n.prop) pridef[[par.curr]][4:5] <- as.character(as.numeric(pridef[[par.curr]][4:5])*sudriv$layout$timestep.fac)
          if(par.curr %in% n.invprop) pridef[[par.curr]][4:5] <- as.character(as.numeric(pridef[[par.curr]][4:5])/sudriv$layout$timestep.fac)
        }else{
          pridef[[par.curr]][1] <- "lognormal"
        }
        m <- as.numeric(pridef[[par.curr]][2])
        if(par.curr %in% n.prop) m <- m + log(sudriv$layout$timestep.fac) ## ATTENTION: this should also be done in case pardef is lognormal, but as of now we don't have any case of a non-transformed, lognormally distributed parameter that is dependent on the time units.
        if(par.curr %in% n.invprop) m <- m - log(sudriv$layout$timestep.fac)
        s <- as.numeric(pridef[[par.curr]][3])
        pridef[[par.curr]][2] <- exp(m + s^2/2)
        pridef[[par.curr]][3] <- as.numeric(pridef[[par.curr]][2])*sqrt(exp(s^2)-1)
      }
      if(which(par.curr==par.names) %in% which(as.logical(par.trans)) & pridef[[par.curr]][1] %in% c("uniform", "Uniform", "unif", "Unif")){
        uni <- TRUE
      }
      ##g.obj <- g.obj + stat_function(data=subset(a.re, param==par.curr), fun = calcpdf, args=list(distpar=pridef[[par.curr]], log=FALSE))
      mu <- as.numeric(pridef[[par.curr]][2])
      sd <- as.numeric(pridef[[par.curr]][3])
      if(is.na(sd)) sd <- mu
      pri.x.max <- ifelse(cut.up,qlnorm(min(1,sqrt(2*mu/sd)),m,s),mu+2*sd)
      if(cut.up & pridef[[par.curr]][1]=="lognormaltrunc") pri.x.max <- pmin(pri.x.max, as.numeric(pridef[[par.curr]][5]))
      rang  <- range(subset(a.re, param==par.curr)$value)
      pri.x.min <- pmax(mu - 2*sd,0)
      if(!prior.only){
        pri.x.max <- pmin(max(pri.x.max, rang[2]), rang[1]+10*(rang[2]-rang[1]))
        pri.x.min <- min(pri.x.min, rang[1])
      }
      pri.x <- seq(pri.x.min, pri.x.max, length.out=l.pri)
      if(!prior.only & scl=="posterior"){
        if(rang[1]>pri.x[1] & rang[2]<pri.x[l.pri]){
          pri.x <- seq(rang[1], rang[2], length.out=l.pri)
        }
      }
      if(uni){
        pri.x <- seq(exp(as.numeric(pridef[[par.curr]][2])), exp(as.numeric(pridef[[par.curr]][3])), length.out=l.pri)
        if(!prior.only & scl=="posterior"){
          if(rang[1]>pri.x[1] & rang[2]<pri.x[l.pri]){
            pri.x <- seq(rang[1], rang[2], length.out=l.pri)
          }
        }
      }
      pri.dens <- calcpdf(pri.x, distpar=pridef[[par.curr]], log=FALSE)
      if(grepl("_df_lik", par.curr)) pri.x <- pri.x + 2
      set.seed(9)
      pri.samp[,par.curr] <- rpdf(n=n,   distpar=pridef[[par.curr]])
      a.pri[(l.pri*(j-1)+1):(l.pri*j),] <- data.frame(value=pri.x, param=par.curr, walker=NA, y=pri.dens, pri=1)
      j <- j + 1
    }
    post.samp <- matrix(NA, nrow=sum(a.re$param==a.re$param[1]), ncol=length(par.names))
    colnames(post.samp) <- par.names
    for(pr in par.names){
      post.samp[,pr] <- a.re$value[a.re$param==pr]
    }
    ## prepare and calculate KL divergence
    if(kl.div){
      ##TotalKL <- KL.divergence(X=post.samp*matrix(rnorm(prod(dim(post.samp)),1,1e-4),nrow=nrow(post.samp)), Y=pri.samp)
      catch <- strsplit(su$settings$subcatchment, split="[0-9]")[[1]][1]
      ind <- gregexpr("[0-9]", su$settings$subcatchment)[[1]][1]
      splt <- strsplit(su$settings$subcatchment,split="")[[1]]
      reso <- ifelse(ind<0,"1h",paste0(splt[ind:length(splt)], collapse=""))
      cat(paste0("catchment ", "reso ", "errmod ", paste(par.names, collapse=" "), "\n", catch, " ", reso, " ", tag), file=file.kl, append=FALSE)
      for(pr in par.names){
        kl <- KL.divergence(X=as.numeric(post.samp[,pr])*rnorm(nrow(post.samp),1,1e-4), Y=as.numeric(pri.samp[,pr]))
        cat(paste0(" ",mean(kl)), file=file.kl, append=TRUE)
      }
    }
    ## add density of prior for plotting (nothing to do with KL divergence)
    a.re <- rbind(a.re, a.pri)
  }
  if(!is.null(v.line)){
    vl <- v.line[levels(as.factor(a.re$param))]
    vline.dat <- data.frame(param=levels(as.factor(a.re$param)), vl=vl)
  }
  ## create list with plots to be arranged with cowplot
  g.list <- list()
  ## prepare labels
  labs  <- par.names
  labs <- gsub("%", "", labs)
  labs <- gsub("_lik", "", labs)
  labs <- gsub("C1Wv_Qstream_", "", labs)
  labs <- gsub("GloCmlt_", "", labs)
  labs <- gsub("GloTrCmlt", "", labs)
  labs <- gsub("Qstream_", "", labs)
  labs <- gsub("U1W", "", labs)
  labs[labs=="a"] <- "a[q]~~\"[-]\""
  labs[labs=="b"] <- "b~~\"[-]\""
  labs[labs=="E"] <- "E~~\"[-]\""
  labs <- gsub("Dspl_SD", "D~~\"[-]\"", labs)
  labs <- gsub("Pmax_ED", "P[ex]~~\"[mm\"~h^{-1}*\"]\"", labs)
  labs <- gsub("GloCmltSmax_IR", "S[list(t,max)]~~\"[mm]\"", labs)
  labs <- gsub("GloCmltSmax_UR", "S[list(u,max)]~~\"[mm]\"", labs)
  labs <- gsub("BeQq_UR", "beta[u]~~\"[-]\"", labs)
  labs <- gsub("K_Qb_UR", "k[list(u,b)]~~\"[\"*h^{-1}*\"]\"", labs)
  labs <- gsub("K_Qq_FR", "k[d]~~\"[\"*h^{-1}*\"]\"", labs)
  labs <- gsub("K_Qq_RR", "k[c]~~\"[\"*h^{-1}*\"]\"", labs)
  labs <- gsub("K_Qq_SR", "k[g]~~\"[\"*h^{-1}*\"]\"", labs)
  labs <- gsub("AlQq_FR", "alpha[d]~~\"[-]\"", labs)
  labs <- gsub("AlQq_SR", "alpha[g]~~\"[-]\"", labs)
  labs <- gsub("K_Qb_SR", "k[list(g,b)]~~\"[\"*h^{-1}*\"]\"", labs)
  labs <- gsub("KpQq_FR", "k[i]~~\"[\"*h^{-1}*\"]\"", labs)
  labs <- gsub("Rs_WR", "r[s]~~\"[\"*h^{-1}*\"]\"", labs)
  labs <- gsub("Kd_WR", "lambda~~\"[\"*h^{-1}*\"]\"", labs)
  labs <- gsub("SlOne_IR", "S[z1]~~\"[mm]\"", labs)
  labs <- gsub("SlTwo_IR", "S[z2]~~\"[mm]\"", labs)
  labs <- gsub("C1Tc1_a", "a[atra]~~\"[-]\"", labs)
  labs <- gsub("C1Tc2_a", "a[terb]~~\"[-]\"", labs)
  labs <- gsub("GLOB_Mult_Q_taumax", "tau[q]~~\"[h]\"", labs)
  labs <- gsub("GLOB_Mult_T_taumax", "tau[c]~~\"[h]\"", labs)
  catch <- strsplit(su$settings$subcatchment, split="[0-9]")[[1]][1]
  ind <- gregexpr("[0-9]", su$settings$subcatchment)[[1]][1]
  splt <- strsplit(su$settings$subcatchment,split="")[[1]]
  reso <- ifelse(ind<0,"1h",paste0(splt[ind:length(splt)], collapse=""))
  j <- 1
  for(par.curr in par.names){
    den <- density(c(post.samp[,par.curr]))
    cls <- "Posterior"
    if(prior.only){den <- list();cls<-NULL}
    pri <- subset(a.re, pri==1 & param==par.curr)
    dat <- rbind(data.frame(value=den$x, y=den$y, class=cls), data.frame(value=pri[,"value"], y=pri[,"y"], class="Prior"))
    dat$class <- factor(dat$class, levels = c("Prior", "Posterior"))
    g.obj <- ggplot(data=dat, mapping=aes(x=value, y=y, fill=class, alpha=class))
    g.obj <- g.obj + geom_area() + theme_bw() + theme(legend.margin=margin(l=1,unit="in"), legend.text=element_text(size=rel(1.5)), legend.title=element_blank(), axis.title.y=element_blank(), axis.title.x=element_text(size=rel(1.5)), axis.text=element_text(size=rel(0.8)), plot.margin=unit(c(ifelse(j<=3,0.4,0.1),0.25,0.1,0.25),"in")) + xlab(label=parse(text=labs[j])) + scale_alpha_discrete(range=c(0.3,0.7)) + scale_fill_brewer(palette="Dark2") + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
    if(scl!="posterior") g.obj <- g.obj + coord_cartesian(ylim=c(0,max(subset(a.re, pri==1 & param==par.curr)$y,na.rm=TRUE)))
    leg <- get_legend(g.obj)
    dev.off()
    g.obj <- g.obj + theme(legend.position="none")
    g.list[[par.curr]] <- g.obj
    j <- j + 1
  }
  pg <- plot_grid(plotlist=g.list,rel_heights=c(1.12,1,1))#leg[[1]]))
  tag.sub <- gsub("P","*",tag)
  tag.sub <- gsub("mean", "", tag.sub)
  pg <- pg + draw_label(label=lab, x=0,y=1,hjust=0,vjust=1,size=18)
  if(plot){
    save_plot(file.hist, pg, ncol=3, base_width=5, base_height=ifelse(length(g.list)>6,9,6))
    dev.off()
  }else{
    return(pg)
  }
}
plot.markov.chain <- function(sudriv, brn.in = 0, thin=1, lower.logpost=NA){
  ## Visualizes marginal parameter distributions of Markov Chains
  ndim <- length(dim(sudriv$parameter.sample))
  if(ndim < 2){stop("parameter.sample might be missing..., plotting nothing.")}
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
  ## make nicer parameter names
  a.re$param <- gsub("%", "", a.re$param)
  a.re$param <- gsub("_lik", "", a.re$param)
  a.re$param <- gsub("C1Wv_Qstream_", "", a.re$param)
  a.re$param <- gsub("GloCmlt_", "", a.re$param)
  a.re$param <- gsub("GloTrCmlt", "", a.re$param)
  a.re$param <- gsub("Qstream_", "", a.re$param)
  ## actual plotting
  g.obj <- ggplot(data=a.re, mapping=aes(x=x,y=value, color=walker)) + geom_line() + facet_wrap("param", nrow=floor(sqrt(dim(a)[2])), scales="free") + theme(legend.position="none")
  plot(g.obj)
}

plot.cor <- function(sudriv, brn.in=0, thin=1, lower.logpost=NA, plot=TRUE){
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
  ## ATTENTION: divide "a" by 5 to represent new notation where sigma_0 was cancelled
  print(head(df))
  df[,grepl("a_lik", par.names)] <- df[,grepl("a_lik", par.names)]/5
  colnames(df) <- gsub("%", "", colnames(s))
  colnames(df) <- gsub("_lik", "", colnames(df))
  colnames(df) <- gsub("C1Wv_Qstream_", "", colnames(df))
  colnames(df) <- gsub("GloCmlt_", "", colnames(df))
  colnames(df) <- gsub("GloTrCmlt", "", colnames(df))
  colnames(df) <- gsub("Qstream_", "", colnames(df))
  colnames(df) <- gsub("U1W", "", colnames(df))
  labels <- colnames(df)
  labels <- gsub("Cmlt_E", "C[E]", labels)
  labels <- gsub("Smax_UR", "S[max]", labels)
  labels <- gsub("K_Qb_UR", "k[u]", labels)
  labels <- gsub("K_Qq_FR", "k[f]", labels)
  labels <- gsub("taumax", "tau[max]", labels)
  labels <- gsub("taumin", "tau[min]", labels)
  print(labels)
  ##labels <- substr(labels, start=nchar(labels)-10, stop=nchar(labels))
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
  pm <- ggpairs(df, lower=list(continuous=myfun), diag=list(continuous=myfun.diag), columnLabels=labels, labeller="label_parsed") + theme_bw(base_size=11)
  if(plot){
    print(pm)
  }else{
    return(pm)
  }
}

plot.predictions <- function(list.su, probs=NA, n.samp=0, sub.set="all", rand=TRUE, xlim=NA, ylim=NA, tme.orig="1000-01-01", lp.num.pred=NA, plt=TRUE, metrics=FALSE, capt.nsamp=FALSE, arrange=NA, plot.var=NA, scl=1, alp=1, loads.det=list(), app.hru.areas=list(), file=NA, type.band=c(par="sample.parunc",par.obs="sample"), type.realiz="par", xintercept=NULL, applic=FALSE, x.ax.tex=TRUE){
  ## ' xlim is a list with an element for each event, which is a vector of length 2: the starting and the end time for that event. The events listed in xlim are plotted side by side.
  translate.var <- c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream","U5F1Wv_Ss1","U5F1Wv_Su1","U3F1Tc1Lv1_Si1","U2F1Wv_Sr1","U3F1Wv_Sf1")
  translate.to <- c(paste0("Streamflow ", ifelse(list.su[[1]]$layout$time.units=="hours", "(mm/h)", "(mm/d)")), expression("Atrazine "*(mu*g/l)), expression("Terbuthylazine "*(mu*g/l)), expression(S[g]~"(mm)"), expression(S[u]~"(mm)"), expression("Atraz. conc. in"~S[t]~(mu*g/l)), expression(S[c]~"(mm)"), expression(S[d]~"(mm)"))
  if(!is.null(xintercept)){ #get minimum and maximum of intervals to shade (assuming xintercept is a list)
    tmp <- sapply(xintercept, function(x) x[1], simplify=FALSE)
    xint.min <- unlist(tmp)
    attributes(xint.min) <- attributes(tmp[[1]])
    tmp <- sapply(xintercept, function(x) x[2], simplify=FALSE)
    xint.max <- unlist(tmp)
    attributes(xint.max) <- attributes(tmp[[1]])
  }
  ## consistency checks
  if(length(type.band)==2 & !all(names(type.band)==c("par","par.obs"))) stop("if length of 'type.band' is 2, must have names 'par' and 'par.obs'")
  ## create data frame for ggplot-object
  if(!is.na(arrange[1]) & length(arrange)!=length(list.su)){warning("length of 'arrange' not equal to length of 'list.su'");return(NA)}
  if(is.na(arrange[1])){arrange <- rep(1,length(list.su));names(arrange) <- names(list.su)}
  if(is.null(names(list.su))) stop("list.su must be named list")
  if(!is.na(xlim[1]) & is.null(names(xlim))) stop("xlim must have named elements")
  area.catch <- 1182895 #m^2 (area of total catchment)
  n.case <- length(list.su)
  ## limit the samples to 'sub.set'
  if(sub.set!="all"){
    for(case.curr in 1:n.case){
      for(bnd in type.band){
        list.su[[case.curr]]$predicted[[bnd]] <- list.su[[case.curr]]$predicted[[bnd]][sub.set,]
      }
    }
  }
  if("C1Wv_Qstream" %in% plot.var){## Adapt streamflow units to timestep factor
    strmflw      <- grepl("Wv_Qstream", list.su[[1]]$layout$layout$var)
    strmflw.pred <- grepl("Wv_Qstream", list.su[[1]]$layout$pred.layout$var)
    for(case.curr in 1:n.case){ # adapt units of streamflow
      list.su[[case.curr]]$predicted$det[1,strmflw.pred] <- list.su[[case.curr]]$predicted$det[1,strmflw.pred]/list.su[[case.curr]]$layout$timestep.fac
      list.su[[case.curr]]$observations[strmflw] <- list.su[[case.curr]]$observations[strmflw]/list.su[[case.curr]]$layout$timestep.fac
      for(bnd in type.band){
        list.su[[case.curr]]$predicted[[bnd]][,strmflw.pred] <- list.su[[case.curr]]$predicted[[bnd]][,strmflw.pred]/list.su[[case.curr]]$layout$timestep.fac
      }
    }
  }
  sudriv <- list.su[[1]]
  ind.sel     <- which(sudriv$layout$pred.layout$var %in% plot.var)
  if(sum(ind.sel)==0){warning("no time period selected"); return(NA)}
  time     <- sudriv$layout$pred.layout$time[ind.sel]
  time.obs <- sudriv$layout$layout$time
  time     <- as.POSIXlt(x=tme.orig)+time*60*60*ifelse(sudriv$layout$time.units=="days",24,1)
  time.obs <- as.POSIXlt(x=tme.orig)+time.obs*60*60*ifelse(sudriv$layout$time.units=="days",24,1)
  obsval <- sudriv$observations
  dt <- sudriv$predicted$det[1,ind.sel]
  if(metrics){
    outside <- obsval > c(apply(sudriv$predicted$sample[,ind.sel], 2, quantile, probs=probs)[2,]) | obsval < c(apply(sudriv$predicted$sample[,ind.sel], 2, quantile, probs=probs)[1,])
    frc <- round(1 - sum(outside)/length(outside), 2)
    mbe <- (sum(obsval)-sum(dt))/sum(obsval)*100
    nse <- 1-sum((dt-obsval)^2)/sum((obsval - mean(obsval))^2)
    capt <- paste("MBE: ", round(mbe), "%, NSE: ", round(nse,2), ", Logpost calib: ", round(lp.num.pred[1]), ", Frac. in bounds: ", frc, sep="")
  }else if(capt.nsamp){
    n.captsamp <- nrow(sudriv$predicted[[type.band[1]]])
    capt <- paste0("Based on ",n.captsamp," samples")
  }else{
    capt <- NULL
  }
  atra <- FALSE
  atra.u3 <- FALSE
  terb <- FALSE
  terb.u3 <- FALSE
  ss <- matrix(NA, nrow=nrow(sudriv$predicted[[type.band[1]]]), ncol=length(ind.sel)*n.case)
  if(length(type.band)==2) ss2 <- ss
  n.water <- sum(sudriv$layout$pred.layout$var[ind.sel]=="C1Wv_Qstream")
  n.atr.u3 <- sum(sudriv$layout$pred.layout$var[ind.sel]=="U3F1Tm1_Qstrm")
  load.atra <- matrix(NA, nrow=nrow(sudriv$predicted[[type.band[1]]]), ncol=n.water*n.case)
  load.atra.u3 <- matrix(NA, nrow=nrow(sudriv$predicted[[type.band[1]]]), ncol=n.atr.u3*n.case)
  load.terb <- matrix(NA, nrow=nrow(sudriv$predicted[[type.band[1]]]), ncol=n.water*n.case)
  if(!is.na(probs[1])){# calculate uncertainty bands
    for(i in 1:n.case){
      sudriv <- list.su[[i]]
      ss.curr <- sudriv$predicted[[type.band[1]]][,ind.sel]
      if(length(type.band)==2) ss.curr2 <- sudriv$predicted[[type.band[2]]][,ind.sel]
      if(("C1Wv_Qstream" %in% plot.var) & ("C1Tc1_Qstream" %in% plot.var)){ #calculate total load of substance exported
        atra <- TRUE
        load.atra[,((i-1)*n.water+1):(i*n.water)] <- ss.curr[,sudriv$layout$pred.layout$var[ind.sel]=="C1Wv_Qstream"]*sudriv$layout$timestep.fac*area.catch*ss.curr[,sudriv$layout$pred.layout$var[ind.sel]=="C1Tc1_Qstream"] # timesetp.fac because streamflow was adapted above
      }
      if("U3F1Tm1_Qstrm" %in% plot.var){ #calculate total load of substance exported
        atra.u3 <- TRUE
        load.atra.u3[,((i-1)*n.atr.u3+1):(i*n.atr.u3)] <- ss.curr[,sudriv$layout$pred.layout$var[ind.sel]=="U3F1Tm1_Qstrm"]
      }
      if(("C1Wv_Qstream" %in% plot.var) & ("C1Tc2_Qstream" %in% plot.var)){ #calculate total load of substance exported
        terb <- TRUE
        load.terb[,((i-1)*n.water+1):(i*n.water)] <- ss.curr[,sudriv$layout$pred.layout$var[ind.sel]=="C1Wv_Qstream"]*sudriv$layout$timestep.fac*area.catch*ss.curr[,sudriv$layout$pred.layout$var[ind.sel]=="C1Tc2_Qstream"] # timestep.fac because streamflow was adapted above
      }
      ss[,((i-1)*length(ind.sel)+1):(i*length(ind.sel))] <- ss.curr
      if(length(type.band)==2) ss2[,((i-1)*length(ind.sel)+1):(i*length(ind.sel))] <- ss.curr2
      if((atra | atra.u3 | terb) & length(type.band)==2) warning("I cannot properly deal with load uncertainties when distinguishing parameter and residual uncertainty")
    }
    quants <- apply(ss, 2, quantile, probs=probs)
    if(length(type.band)==2){
      quants <- rbind(quants, apply(ss2, 2, quantile, probs=probs))
    }else{
      quants <- rbind(quants, quants)
    }
  }else{
    quants <- data.frame(rbind(NA,NA,NA,NA))
  }
  
      
  if(n.samp > 0){## plot actual realisations
    preds <- numeric()
    if(type.realiz=="par") smp <- type.band[type.realiz]
    for(i in 1:n.case){
      if(rand){
        ss <- list.su[[i]]$predicted[[smp]][sample(1:nrow(sudriv$predicted[[smp]]),n.samp),ind.sel,drop=FALSE]
      }else{
        ss <- list.su[[i]]$predicted[[smp]][1:min(n.samp, nrow(sudriv$predicted[[smp]])),ind.sel,drop=FALSE]
      }
      dms <- dim(ss)
      preds <- c(preds,array(t(ss), dim=c(prod(dms), 1)))
    }
    if(length(list.su)==1){
      stoch.realiz.name <- "Stoch. realiz."
    }else{
      stoch.realiz.name <- paste0(names(list.su), " stoch")
    }
    stoch <- data.frame(x=rep(time,n.case*n.samp), value=c(preds), var=rep(sudriv$layout$pred.layout[ind.sel,"var"], n.case*n.samp), simu=rep(paste0(stoch.realiz.name, ifelse(dms[1]>1,1:(dms[1]),"")), each = dms[2]), lower=c(preds), lower2=c(preds), upper=c(preds), upper2=c(preds))
  }else{
    stoch <- data.frame()
  }
  obs   <- data.frame(x=time.obs, value=obsval, var=sudriv$layout$layout[,"var"], simu="Observed", lower=obsval, lower2=obsval, upper=obsval, upper2=obsval)
  # expand dt if there are multiple models
  if(n.case>1){for(i in 2:n.case){dt <- c(dt,list.su[[i]]$predicted$det[1,ind.sel])}}
  det <-   data.frame(x=rep(time,n.case), value = c(dt), var=rep(sudriv$layout$pred.layout[ind.sel,"var"], n.case), simu=paste(rep(names(list.su),each=length(time)),sep=""), lower=c(quants[1,]), lower2=c(quants[3,]), upper=c(quants[2,]), upper2=c(quants[4,]))
  data.plot <- rbind(det, stoch, obs)
  ## combine the upper and the lower ribbons into one variable each, so that it can be passed to an aesthetics for plotting
  data.plot <- data.plot %>% pivot_longer(cols=c(lower, lower2), names_to="lw.nm", values_to="lw")
  data.plot <- data.plot %>% mutate(lw.nm=replace(lw.nm,lw.nm=="lower","par")) %>% mutate(lw.nm=replace(lw.nm,lw.nm=="lower2","par.obs"))
  data.plot <- data.plot %>% pivot_longer(cols=c(upper, upper2), names_to="up.nm", values_to="up")
  data.plot <- data.plot %>% mutate(up.nm=replace(up.nm,up.nm=="upper","par")) %>% mutate(up.nm=replace(up.nm,up.nm=="upper2","par.obs"))
  data.plot <- data.plot %>% filter(lw.nm==up.nm) %>% rename(typ=lw.nm) %>% select(-up.nm) %>% mutate(alp=replace(typ,typ=="par",0.2)) %>% mutate(alp=as.numeric(replace(alp,alp=="par.obs",0.6)))
  ## remove observational uncertainty for states for which we have none
  noerrmod <- plot.var[!sapply(plot.var, function(x) any(grepl(paste0(x,".*_lik"), names(list.su[[1]]$likelihood$parameters))))]
  data.plot <- data.plot %>% mutate(lw=replace(lw, var%in%noerrmod & typ=="par.obs", NA)) %>% mutate(up=replace(up, var%in%noerrmod & typ=="par.obs", NA))
  ## add NA observations for those states to keep plot colors consistent
  dummy.dat <- data.plot %>% filter(var==unique(var)[1])
  dummy.dat2 <- list()
  for(i in 1:length(noerrmod)){
    dummy.dat$var <- noerrmod[i]
    dummy.dat$simu <- "Observed"
    dummy.dat$value <- NA
    dummy.dat$lw <- NA
    dummy.dat$up <- NA
    dummy.dat$alp <- NA
    dummy.dat2 <- c(dummy.dat2, list(dummy.dat))
  }
  dummy.dat <- do.call(rbind, dummy.dat2)
  data.plot <- rbind(data.plot, dummy.dat)
    
  ## good so far...
  ## actual plotting
  n <- n.samp+1
  if(is.na(plot.var[1])){ # plot all states
    plot.var <- unique(data.plot$var)
  }
  i <- 1
  ij <- 1
  g.objs <- list()
  loads.atra <- list()
  loads.atra.u3 <- list()
  loads.terb <- list()
  loads.atra.all <- data.frame(x=numeric(), simu=character())
  loads.atra.u3.all <- data.frame(x=numeric(), simu=character())
  loads.terb.all <- data.frame(x=numeric(), simu=character())
  xlim.q <- xlim[!grepl("Total", names(xlim))]
  for(event.curr in xlim){
    tmp <- make.breaks(event.curr)
    brks <- tmp$brks
    frmt <- tmp$frmt
    print(i)
    j <- 1
    if(!grepl("Total", names(xlim)[i])){
      for(var.curr in plot.var){
        print(var.curr)
        for(panel.curr in unique(arrange)){# create the ggplot object for each panel
          print(panel.curr)
          ## get index of rows of su objects of current panel
          last <- j==length(unique(arrange))*length(plot.var)
          cases <- names(arrange[arrange==panel.curr])
          data.curr <- data.plot %>% filter(var==var.curr & x>=event.curr[1] & x<=event.curr[2])
          data.curr <- data.curr %>% mutate(typ=gsub("par.obs","residual",typ), typ=gsub("par","intrinsic",typ))
          g.obj <- ggplot(data=data.curr, mapping=aes(x=x,y=value,color=simu,linetype=simu))
          if(!is.na(probs[1])){
            g.obj <- g.obj + geom_ribbon(aes(ymin=lw,ymax=up,alpha=typ), data=data.curr, linetype=ifelse(length(cases)>1, "solid", 0)) + scale_alpha_manual(values=c("intrinsic"=0.5,"residual"=0.3))
          }
          g.obj <- g.obj + geom_line(data=data.curr%>%distinct(x,value,simu))
          if(applic) g.obj <- g.obj + geom_vline(xintercept=as.POSIXct("2009-05-19 12:00"), linetype="dashed", size=0.5, color="red")
          if(!is.null(xintercept)) g.obj <- g.obj + annotate("rect", xmin=xint.min, xmax=xint.max, ymin=0, ymax=Inf, fill="blue", alpha=0.2)
          g.obj <- g.obj + theme_light() + theme(text=element_text(size=12), plot.margin=unit(c(ifelse(j==1,0.1,0),0.01,ifelse(last&x.ax.tex,0.1,-0.3),ifelse(i==1,0.2,0.1)), "cm"), 
                                              legend.position=ifelse(i==length(xlim.q) & j==2,"right","none"), legend.title=element_text(size=14), legend.text=element_text(size=14)) + 
            labs(linetype="", color="", x="", y=translate.to[translate.var==var.curr], alpha="Stochast.") + 
            scale_x_datetime(date_breaks=brks, date_labels=frmt, limits=c(event.curr)) + 
            scale_y_continuous(expand=c(0.01,0)) + scale_color_viridis(discrete=TRUE)
          if(last & x.ax.tex){g.obj <- g.obj + theme(axis.text.x=element_text(size=8))}else{g.obj <- g.obj + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())}
          if(i!=1) g.obj <- g.obj + theme(axis.title.y=element_blank())
          if(!is.na(ylim[1])){
            if(!all(names(ylim) %in% plot.var)) stop("names for list ylim not found")
            g.obj <- g.obj + coord_cartesian(ylim=ylim[[var.curr]])
          }
          g.objs[[ij]] <- g.obj
          j <- j + 1
          ij <- ij + 1
        }
      }
    }
    if(atra){
      t1=time[sudriv$layout$pred.layout$var[ind.sel]=="C1Wv_Qstream"]
      evnt <- rep(t1,n.case)>event.curr[1] & rep(t1,n.case)<=event.curr[2]
      ll <- load.atra[,evnt]
      cat("here: ", which.max(rowSums(ll)), "\n")
      if(names(xlim)[1]=="E0") save(ll, file="loads_atra_E0.RData", version=2)
      loads.atra[[i]] <- apply(load.atra[,evnt], 1, function(x,list.su,t1,evnt) c(tapply(x, rep(names(list.su), each=length(t1))[evnt], sum)), list.su=list.su, t1=t1, evnt=evnt)
      if(n.case==1){
        loads.atra[[i]] <- array(loads.atra[[i]], dim=c(1,length(loads.atra[[i]])))## make sure that the dimension of loads.atra is stable if n.case > 1
        dimnames(loads.atra[[i]]) <- list(c(names(list.su)), NULL)
      }
      cat("here: ", which.max(as.numeric(loads.atra[[i]])), "\n")
      loads.atra.all <- rbind(loads.atra.all, data.frame(x=c(t(loads.atra[[i]])[,names(list.su)]), simu=rep(names(list.su), each=nrow(load.atra)), event=names(xlim)[i]))
    }
    if(atra.u3){
      t1=time[sudriv$layout$pred.layout$var[ind.sel]=="U3F1Tm1_Qstrm"]
      evnt <- rep(t1,n.case)>event.curr[1] & rep(t1,n.case)<=event.curr[2]
      loads.atra.u3[[i]] <- apply(load.atra.u3[,evnt], 1, function(x,list.su,t1,evnt) c(tapply(x, rep(names(list.su), each=length(t1))[evnt], sum)), list.su=list.su, t1=t1, evnt=evnt)
      if(n.case==1) loads.atra.u3[[i]] <- array(loads.atra.u3[[i]], dim=c(1,length(loads.atra.u3[[i]])))## make sure t
      loads.atra.u3.all <- rbind(loads.atra.u3.all, data.frame(x=c(t(loads.atra.u3[[i]])[,names(list.su)]), simu=rep(names(list.su), each=nrow(load.atra.u3)), event=names(xlim)[i]))
    }
    if(terb){
      t1=time[sudriv$layout$pred.layout$var[ind.sel]=="C1Wv_Qstream"]
      evnt <- rep(t1,n.case)>event.curr[1] & rep(t1,n.case)<=event.curr[2]
      loads.terb[[i]] <- apply(load.terb[,evnt], 1, function(x,list.su,t1,evnt) c(tapply(x, rep(names(list.su), each=length(t1))[evnt], sum)), list.su=list.su, t1=t1, evnt=evnt)
      if(n.case==1){
        loads.terb[[i]] <- array(loads.terb[[i]], dim=c(1,length(loads.terb[[i]])))## make sure that the dimension of loads.terb is stable if n.case > 1
        dimnames(loads.terb[[i]]) <- list(c(names(list.su)), NULL)
      }
      loads.terb.all <- rbind(loads.terb.all, data.frame(x=c(t(loads.terb[[i]])[,names(list.su)]), simu=rep(names(list.su), each=nrow(load.terb)), event=names(xlim)[i]))
    }
    i <- i + 1
  }
  if(atra & !is.null(app.hru.areas$atra)) loads.atra.all$x <- loads.atra.all$x/1000/1000 # convert from micro g to g
  if(terb & !is.null(app.hru.areas$atra)) loads.terb.all$x <- loads.terb.all$x/1000/1000 # convert from micro g to g
  print(app.hru.areas$atra)
  if(names(xlim)[1]=="E0") save(loads.atra.all, file="loads_atra_all.RData", version=2)
  if(atra & terb & !is.null(app.hru.areas$atra)){
    brks.rel <- c(0.001,seq(0.001,0.01,0.001), seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1))
    brks.abs <- c(0.01,seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10))
    brks.ms.rel <- c(0.01,seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1))
    loads.atra.rel <- loads.atra.all
    loads.atra.rel$x <- loads.atra.rel$x/8269.5*100
    loads.terb.rel <- loads.terb.all
    first.applic <- loads.terb.rel[,"event"] %in% c("E00", "E0", "E1","E2")
    loads.terb.rel[first.applic,"x"] <- loads.terb.rel[first.applic,"x"]/5594.5*100
    loads.terb.rel[!first.applic,"x"] <- loads.terb.rel[!first.applic,"x"]/(5594.5+4252)*100
    if(n.case==1){
      gg.atra.abs <- ggplot(data=loads.atra.all, aes(x=x, fill=event, color=event)) + geom_density(alpha=alp) + scale_x_log10(breaks=brks.abs, labels=every_nth(brks.abs, 5, inverse=TRUE)) + labs(x=expression("Exported atrazine (g)"), y="Probability (-)", fill="Event", color="Event")
      gg.atra <- ggplot(data=loads.atra.rel, aes(x=x, fill=event, color=event)) + geom_density(alpha=alp) + scale_x_log10(breaks=brks.rel, labels=every_nth(brks.rel, 5, inverse=TRUE)) + labs(x=expression("Exported atrazine (% of applied)"), y="Probability (-)", fill="Event", color="Event")
      gg.terb.abs <- ggplot(data=loads.terb.all, aes(x=x, fill=event, color=event)) + geom_density(alpha=alp) + scale_x_log10(breaks=brks.abs, labels=every_nth(brks.abs, 5, inverse=TRUE)) + labs(x=expression("Exported terbuthylazine (g)"), y="Probability (-)", fill="Event", color="Event")
      gg.terb <- ggplot(data=loads.terb.rel, aes(x=x, fill=event, color=event)) + geom_density(alpha=alp) + scale_x_log10(breaks=brks.rel, labels=every_nth(brks.rel, 5, inverse=TRUE)) + labs(x=expression("Exported terbuthylazine (% of applied)"), y="Probability (-)", fill="Event", color="Event")
    }else{
      gg.atra.abs <- ggplot(data=loads.atra.all, aes(x=x, y=simu, fill=event)) + geom_density_ridges(scale=scl, alpha=alp) + scale_x_log10(breaks=brks.abs, labels=every_nth(brks.abs, 5, inverse=TRUE)) + labs(x=expression("Exported atrazine (g)"), y="", fill="Event") + theme_ridges() + scale_y_discrete(expand = c(0.01, 0))
      gg.atra <- ggplot(data=loads.atra.rel, aes(x=x, y=simu, fill=event)) + geom_density_ridges(scale=scl, alpha=alp) + scale_x_log10(breaks=brks.rel, labels=every_nth(brks.rel, 5, inverse=TRUE)) + labs(x=expression("Exported atrazine (% of applied)"), y="", fill="Event") + theme_ridges() + scale_y_discrete(expand = c(0.01, 0))
      gg.terb.abs <- ggplot(data=loads.terb.all, aes(x=x, y=simu, fill=event)) + geom_density_ridges(scale=scl, alpha=alp) + scale_x_log10(breaks=brks.abs, labels=every_nth(brks.abs, 5, inverse=TRUE)) + labs(x=expression("Exported terbuthylazine (g)"), y="", fill="Event") + theme_ridges() + scale_y_discrete(expand = c(0.01, 0))
      gg.terb <- ggplot(data=loads.terb.rel, aes(x=x, y=simu, fill=event)) + geom_density_ridges(scale=scl, alpha=alp) + scale_x_log10(breaks=brks.rel, labels=every_nth(brks.rel, 5, inverse=TRUE)) + labs(x=expression("Exported terbuthylazine (% of applied)"), y="", fill="Event") + theme_ridges() + scale_y_discrete(expand = c(0.01, 0))
    }
    ## these are the masses exported from the hrus in the maximum posterior parameter set. They are used to calculate the fraction of the HRU contributions in the stochastic case
    flux.atra.tot <- loads.det$flux.atra.tot
    flux.terb.tot <- loads.det$flux.terb.tot
    ## areas in each HRU to which atrazine and terbuthylazine was applied
    atra.app.hru.areas <- app.hru.areas$atra
    terb.app.hru.areas <- app.hru.areas$terb
    loads.atra.tot <- subset(loads.atra.all, event==ifelse(length(xlim)>1,"Total",names(xlim)[1])) # convert total exports to absolute amount (g)
    loads.terb.tot <- subset(loads.terb.all, event==ifelse(length(xlim)>1,"Total",names(xlim)[1])) # convert total exports to absolute amount (g)
    x <- rep(NA, sum(table(loads.atra.tot$simu) * unlist(lapply(atra.app.hru.areas[names(list.su)], length))))
    loads.atra.hru <- data.frame(x=x, simu=NA, hru=NA)
    loads.atra.hru.rel <- data.frame(x=x, simu=NA, hru=NA)
    x <- rep(NA, sum(table(loads.terb.tot$simu) * unlist(lapply(terb.app.hru.areas[names(list.su)], length))))
    loads.terb.hru <- data.frame(x=x, simu=NA, hru=NA)
    loads.terb.hru.rel <- data.frame(x=x, simu=NA, hru=NA)
    rwind <- 1:(table(loads.atra.tot$simu)[1])
    for(case.curr in names(list.su)){
      print(case.curr)
      mdl <- loads.atra.tot$simu==case.curr
      print("mdl")
      print(sum(mdl))
      for(hru in 1:length(flux.atra.tot[[case.curr]])){
        print(hru)
        if(!(case.curr == (names(list.su)[1]) & hru==1)){
          rwind <- (rwind[length(rwind)]+1):(rwind[length(rwind)]+table(loads.atra.tot$simu)[names(list.su)==case.curr])
          print("rwind:")
          print(range(rwind))
        }
        if(hru==1){
          loads.atra.hru[rwind,"x"] <- NA#loads.atra.tot$x*8269.5/100*flux.atra.tot[hru]/sum(flux.atra.tot) # absolute mass exported per HRU (g)
          loads.atra.hru.rel[rwind,"x"] <- NA#loads.atra.hru[rwind,"x"] / (0.1*115.00287/1000/1000*17410) * 100 # divide by the absolute mass of atrazine sprayed on impervious areas (taken from the input)
          loads.terb.hru[rwind,"x"] <- NA#loads.terb.tot$x*(5594.5+4252)/100*flux.terb.tot[hru]/sum(flux.terb.tot) # absolute mass exported per HRU (g)
          loads.terb.hru.rel[rwind,"x"] <- NA#loads.terb.hru[rwind,"x"] / (0.1*115.00287/1000/1000*17410) * 100 # divide by the absolute mass of terbuthylazine sprayed on impervious areas (taken from the input)
        }else{
          loads.atra.hru[rwind,"x"] <- loads.atra.tot$x[mdl]*flux.atra.tot[[case.curr]][hru]/sum(flux.atra.tot[[case.curr]]) # absolute mass exported per HRU (g)
          loads.atra.hru.rel[rwind,"x"] <- loads.atra.tot$x[mdl]/8269.5*100*flux.atra.tot[[case.curr]][hru]/sum(flux.atra.tot[[case.curr]])/(atra.app.hru.areas[[case.curr]][hru]/sum(atra.app.hru.areas[[case.curr]],na.rm=TRUE)) # exported mass per HRU relative to applied mass in that HRU
          loads.terb.hru[rwind,"x"] <- loads.terb.tot$x[mdl]*flux.terb.tot[[case.curr]][hru]/sum(flux.terb.tot[[case.curr]]) # absolute mass exported per HRU (g)
          M.terb.hru <- terb.app.hru.areas[[case.curr]][["first"]][hru]/sum(terb.app.hru.areas[[case.curr]][["first"]],na.rm=TRUE)*5594.5 + terb.app.hru.areas[[case.curr]][["second"]][hru]/sum(terb.app.hru.areas[[case.curr]][["second"]],na.rm=TRUE)*4252
          loads.terb.hru.rel[rwind,"x"] <- loads.terb.tot$x[mdl]*100*flux.terb.tot[[case.curr]][hru]/sum(flux.terb.tot[[case.curr]]) / M.terb.hru # exported mass per HRU relative to applied mass in that HRU
          ## loads.terb.hru.rel[rwind,"x"] <- loads.terb.tot$x[mdl]/4252*100*flux.terb.tot[[case.curr]][hru]/sum(flux.terb.tot[[case.curr]])/(terb.app.hru.areas[[case.curr]][["second"]][hru]/sum(terb.app.hru.areas[[case.curr]][["second"]],na.rm=TRUE)) # exported mass per HRU relative to applied mass in that HRU
        }
        print(table(loads.atra.tot$simu[mdl]))
        loads.atra.hru[rwind,"simu"] <- loads.atra.tot$simu[mdl]
        loads.atra.hru.rel[rwind,"simu"] <- loads.atra.tot$simu[mdl]
        print(names(flux.atra.tot[[case.curr]])[hru])
        loads.atra.hru[rwind,"hru"] <- names(flux.atra.tot[[case.curr]])[hru]
        loads.atra.hru.rel[rwind,"hru"] <- names(flux.atra.tot[[case.curr]])[hru]
        loads.terb.hru[rwind,"simu"] <- loads.terb.tot$simu[mdl]
        loads.terb.hru.rel[rwind,"simu"] <- loads.terb.tot$simu[mdl]
        loads.terb.hru[rwind,"hru"] <- names(flux.terb.tot[[case.curr]])[hru]
        loads.terb.hru.rel[rwind,"hru"] <- names(flux.terb.tot[[case.curr]])[hru]
      }
    }
    ## remove the 'connected and drained' hru
    loads.atra.hru     %<>% filter(hru != "Connected and Drained")
    loads.atra.hru.rel %<>% filter(hru != "Connected and Drained")
    loads.terb.hru     %<>% filter(hru != "Connected and Drained")
    loads.terb.hru.rel %<>% filter(hru != "Connected and Drained")
    sm <- loads.atra.hru.rel %>% group_by(simu, hru) %>% summarise(mean_exprt=mean(x), q05=quantile(x,0.05,na.rm=TRUE), q95=quantile(x,0.95,na.rm=TRUE))
    write.table(sm, file="atra_export_summary.txt", row.names=FALSE)
    sm <- loads.terb.hru.rel %>% group_by(simu, hru) %>% summarise(mean_exprt=mean(x), q05=quantile(x,0.05, na.rm=TRUE), q95=quantile(x,0.95,na.rm=TRUE))
    write.table(sm, file="terb_export_summary.txt", row.names=FALSE)
    if(n.case==1){
      gg.atra.hru <- ggplot(data=loads.atra.hru, aes(x=x, fill=hru, color=str_wrap(hru,12))) + geom_density(alpha=alp) + scale_x_log10(breaks=brks.abs, labels=every_nth(brks.abs, 5, inverse=TRUE)) + labs(x=expression("Exported atrazine (g)"), y="Probability (-)", fill="HRU", color="HRU")
      gg.atra.hru.rel <- ggplot(data=loads.atra.hru.rel, aes(x=x, fill=hru, color=str_wrap(hru,12))) + geom_density(alpha=alp) + scale_x_log10(breaks=brks.ms.rel, labels=every_nth(brks.ms.rel, 5, inverse=TRUE)) + labs(x=expression("Exported atrazine (% of applied)"), y="Probability (-)", fill="HRU", color="HRU")
      gg.terb.hru <- ggplot(data=loads.terb.hru, aes(x=x, fill=hru, color=str_wrap(hru,12))) + geom_density(alpha=alp) + scale_x_log10(breaks=brks.abs, labels=every_nth(brks.abs, 5, inverse=TRUE)) + labs(x=expression("Exported terbuthylazine (g)"), y="Probability (-)", fill="HRU", color="HRU")
      gg.terb.hru.rel <- ggplot(data=loads.terb.hru.rel, aes(x=x, fill=hru, color=str_wrap(hru,12))) + geom_density(alpha=alp) + scale_x_log10(breaks=brks.ms.rel, labels=every_nth(brks.ms.rel, 5, inverse=TRUE)) + labs(x=expression("Exported terbuthylazine (% of applied)"), y="Probability (-)", fill="HRU", color="HRU")
    }else{
      gg.atra.hru <- ggplot(data=loads.atra.hru, aes(x=x, y=simu, fill=str_wrap(hru,12))) + geom_density_ridges(scale=scl, alpha=alp) + scale_x_log10(breaks=brks.abs, labels=every_nth(brks.abs, 5, inverse=TRUE)) + labs(x=expression("Exported atrazine (g)"), y="", fill="HRU") + theme_ridges() + scale_y_discrete(expand = c(0.01, 0))
      gg.atra.hru.rel <- ggplot(data=loads.atra.hru.rel, aes(x=x, y=simu, fill=str_wrap(hru,12))) + geom_density_ridges(scale=scl, alpha=alp) + scale_x_log10(breaks=brks.ms.rel, labels=every_nth(brks.ms.rel, 5, inverse=TRUE)) + labs(x=expression("Exported atrazine (% of applied)"), y="", fill="HRU") + theme_ridges() + scale_y_discrete(expand = c(0.01, 0))
      gg.terb.hru <- ggplot(data=loads.terb.hru, aes(x=x, y=simu, fill=str_wrap(hru,12))) + geom_density_ridges(scale=scl, alpha=alp) + scale_x_log10(breaks=brks.abs, labels=every_nth(brks.abs, 5, inverse=TRUE)) + labs(x=expression("Exported terbuthylazine (g)"), y="", fill="HRU") + theme_ridges() + scale_y_discrete(expand = c(0.01, 0))
      gg.terb.hru.rel <- ggplot(data=loads.terb.hru.rel, aes(x=x, y=simu, fill=str_wrap(hru,12))) + geom_density_ridges(scale=scl, alpha=alp) + scale_x_log10(breaks=brks.ms.rel, labels=every_nth(brks.ms.rel, 5, inverse=TRUE)) + labs(x=expression("Exported terbuthylazine (% of applied)"), y="", fill="HRU") + theme_ridges() + scale_y_discrete(expand = c(0.01, 0))
    }
  }
  if(plt){
    cat("plotting ...\n")
    if(!is.na(file)) pdf(file=file, width=10, height=7)
    if(!is.null(capt)){
      capt <- textGrob(capt,gp = gpar(fontface = 3, fontsize = 9),
      hjust = 1,
      x = 1)
    }
    egg::ggarrange(plots=g.objs, nrow=length(plot.var), byrow=FALSE, newpage=FALSE, bottom=capt)
    if(!is.na(file)) dev.off()
    if(atra & terb & !is.null(app.hru.areas$atra)){
      pub1 <- ggpubr::ggarrange(gg.atra,
                                gg.terb, legend="right", common.legend=TRUE, nrow=2, labels=c("a","b"))
      pub2 <- ggpubr::ggarrange(gg.atra.hru.rel,
                                gg.terb.hru.rel, legend="right", common.legend=TRUE, nrow=2, labels=c("c","d"))
      if(!is.na(file)){ggexport(pub1, pub2, filename=gsub(".pdf","_export.pdf",file), nrow=2)}
      ##grid.arrange(gg.atra.abs, gg.terb.abs, gg.atra.hru, gg.terb.hru, ncol=1, newpage=TRUE)
      pub1 <- ggpubr::ggarrange(gg.atra.abs,
                                gg.terb.abs, legend="right", common.legend=TRUE, nrow=2, labels=c("a","b"))
      pub2 <- ggpubr::ggarrange(gg.atra.hru,
                                gg.terb.hru, legend="right", common.legend=TRUE, nrow=2, labels=c("c","d"))
      if(!is.na(file)){ggexport(pub1, pub2, filename=gsub(".pdf","_export2.pdf",file), nrow=2)}
    }
  }else{
    warning("returning ggplot not implemented properly. Returning only part of the ggplot objects ...")
    return(g.objs)
  }
}
make.breaks <- function(limits){
  period <- (limits[2] - limits[1]) #duration of current event in days, used to calculate the breaks
  if(period <= 1){brks <- "12 hours"; frmt <- "%d.%m. %H:%M"}
  if(period > 1){brks <- "1 day"; frmt <- "%d.%m. %H:%M"}
  if(period > 3){brks <- "1 day"; frmt <- "%d.%m"}
  if(period > 5){brks <- "4 days"; frmt <- "%d.%m"}
  if(period > 20){brks <- "1 week"; frmt <- "%d.%m"}
  if(period > 60){brks <- "1 month"; frmt <- "%d.%m"}
  return(list(brks=brks, frmt=frmt))
}

plot.ts.quantiles <- function(dat, sudriv, xlim=NA, ind.sel=NA, precip=FALSE, plim=0, plot=TRUE){
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
    dat.rect <- dat.rect[-1,]
    for(i in 2:(table(dat.ind.sel$case)[1])){
      if(dat.ind.sel[i,n.p] & (!dat.ind.sel[i-1,n.p] | i==2)) dat.rect <- rbind(dat.rect,data.frame(from=dat.ind.sel[i,"x"], to=dat.ind.sel[i,"x"]))
      if((!dat.ind.sel[i,n.p] | i==(table(dat.ind.sel$case)[1])) &  dat.ind.sel[i-1,n.p]) dat.rect[nrow(dat.rect),2] <- dat.ind.sel[i-1,"x"]
    }
  }
  dat.rect$pr <- "Precipitation"
  dat.ind.sel[,"case"] <- gsub("P", "\u002A", dat.ind.sel[,"case"])
  g.obj1 <- ggplot() + geom_point(mapping=aes(x=x,y=quant), data=dat.ind.sel, size=0.8) + geom_line(mapping=aes(x=x,y=quant), data = dat.ind.sel)
  if(precip) g.obj1 <- g.obj1 + geom_rect(data=dat.rect, aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=pr), alpha=0.3) + scale_fill_manual(values=c("grey"))
  g.obj1 <- g.obj1 + labs(x="", y=expression(eta), fill="")+ theme_bw() + theme(axis.text=element_text(size=12)) + scale_x_datetime(labels = date_format("%b %y"))##+ scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
  if(plot){
    plot(g.obj1)
  }else{
    return(g.obj1)
  }
}
plot.ts.white.noise <- function(dat, sudriv, xlim=NA, ind.sel=NA, precip=FALSE, plim=0, plot=TRUE){
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
    dat.rect <- dat.rect[-1,]
    for(i in 2:(table(dat.ind.sel$case)[1])){
      if(dat.ind.sel[i,n.p] & ! dat.ind.sel[i-1,n.p]) dat.rect <- rbind(dat.rect,data.frame(from=dat.ind.sel[i,"x"], to=dat.ind.sel[i,"x"]))
      if((!dat.ind.sel[i,n.p] | i==(table(dat.ind.sel$case)[1])) &  dat.ind.sel[i-1,n.p]) dat.rect[nrow(dat.rect),2] <- dat.ind.sel[i-1,"x"]
    }
  }
  dat.rect$pr <- "Precipitation"
  dat.ind.sel[,"case"] <- gsub("P", "\u002A", dat.ind.sel[,"case"])
  g.obj1 <- ggplot() + geom_point(mapping=aes(x=x,y=white.noise, shape=case, color=case), data=dat.ind.sel, size=1.2) + geom_line(mapping=aes(x=x,y=white.noise, color=case), data = dat.ind.sel)
  if(precip) g.obj1 <- g.obj1 + geom_rect(data=dat.rect, aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=pr), alpha=0.3) + scale_fill_manual(values=c("grey"))
  g.obj1 <- g.obj1 + labs(x="", y=expression("Standardized innovations,"~~chi), shape="Error Model", color="Error Model", fill="")+ theme_bw() + theme(text=element_text(size=20)) ##+ scale_x_datetime(date_breaks="1 week", date_labels="%d %b")
  if(plot){
    plot(g.obj1)
  }else{
    return(g.obj1)
  }
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
plot.Qdet.quantiles <- function(dat, sudriv, xlim=NA, ind.sel=NA, plot=TRUE){
  dat[,"det"] <- dat[,"det"]/sudriv$layout$timestep.fac
  g.obj <- ggplot(data=dat[ind.sel,], mapping=aes(x=det, y=quant)) + geom_point()+ geom_abline(slope=0, intercept=0) + theme_bw(base_size=18) + theme(axis.text=element_text(size=24), axis.title=element_text(size=18)) + scale_x_log10() + labs(x=expression(Q["det"]*" [mm "*h^{-1}*"]"), y=expression(eta))
  if(plot){
    plot(g.obj)
  }else{
    return(g.obj)
  }
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
plot.results.summary <- function(files=NA,outpath="sudriv_output/"){
  murg <- read.table(files["murg"], sep=" ", header=TRUE)
  maimai <- read.table(files["maimai"], sep=" ", header=TRUE)
  murg$catchment <- "Murg"
  maimai$catchment <- "Maimai"
  dat <- rbind(maimai,murg)
  ## decide where to use the smoothed version and where the original one
  ## dat <- dat[-which(dat$catchment=="Maimai" & dat$reso=="1h" & dat$errmod %in% c("E3", "E4")),]
  ## dat <- dat[-which(dat$catchment=="Murg" & dat$reso=="1h" & dat$errmod =="E3aP"),]
  dat$reso <- gsub("h", "", dat$reso)
  xx <- dat$reso
  xx[xx=="6"] <- 2
  xx[xx=="24"] <- 3
  xx <- as.numeric(xx)
  dat$reso <- reorder(as.factor(dat$reso),X=xx)
  ##dat$reso <- as.factor(dat$reso)
  dat$errmod <- gsub("P", "", dat$errmod)
  dat$errmod <- gsub("mean", "", dat$errmod)
  dat$errmod[dat$errmod=="E3"] <- "'E3(\\u002A)'"
  dat$errmod[dat$errmod=="E3a"] <- "'E3a(\\u002A)'"
  dat$errmod[dat$errmod=="E4"] <- "'E4(\\u002A)'"
  dat$errmod[dat$errmod=="E4a"] <- "'E4a(\\u002A)'"
  xx <- dat$errmod
  xx[grepl("E1",xx)] <- 1
  xx[grepl("E2",xx)] <- 2
  xx[grepl("E3",xx)] <- 3
  xx[grepl("E3a",xx)] <- 4
  xx[grepl("E4",xx)] <- 5
  xx[grepl("E4a",xx)] <- 6
  xx <- as.numeric(xx)
  dat$errmod <- reorder(as.factor(dat$errmod),X=xx)
  dat$case <- paste(as.character(dat$catchment), as.character(dat$errmod))
  ## dat$meas <- 1-dat$reli^(1-dat$prec)
  ## dat$meas.valid <- 1-dat$reli.valid^(1-dat$prec.valid)
  notext <- element_blank()
  ## colours
  my_palette = c("#000000", brewer.pal(3, "Set1")[1], brewer.pal(9, "Greens")[c(7,5)], brewer.pal(9, "Purples")[c(7,4)])
  scale_colour_discrete = function(...) scale_colour_manual(..., values = palette())
  palette(my_palette)
  dat.simp <- subset(dat,!grepl("E4",errmod))
  #reliability
  g.reli <- ggplot(data=dat, aes(x=reso, y=reli)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=reli, group=case, colour=errmod, linetype=catchment))+geom_blank(aes(y=reli.valid))+geom_hline(yintercept=0)+labs(shape="Catchment", colour="Error Model")+labs(y="Reliability [-]")+scale_x_discrete(expand=c(0.1,0.1))+scale_y_continuous(expand=c(0,0))+theme_bw()
  g.reli.simp <- ggplot(data=dat.simp, aes(x=reso, y=reli, shape=catchment, colour=errmod)) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(aes(y=reli.valid))+geom_hline(yintercept=0)+labs(shape="Catchment", colour="Error Model")+labs(y="Reliability [-]")+theme_bw()
  g.reli.valid <- ggplot(data=dat, aes(x=reso, y=reli.valid)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=reli.valid, group=case, colour=errmod, linetype=catchment))+geom_blank(aes(y=reli))+geom_hline(yintercept=0)+scale_x_discrete(expand=c(0.1,0.1))+scale_y_continuous(expand=c(0,0))+theme_bw()
  g.reli.valid.simp <- ggplot(data=dat.simp, aes(x=reso, y=reli.valid, shape=catchment, colour=errmod)) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(aes(y=reli))+geom_hline(yintercept=0)+theme_bw()
  #spread
  g.spread <- ggplot(data=dat, aes(x=reso, y=spread)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=spread, group=case, colour=errmod, linetype=catchment))+geom_blank(aes(y=spread.valid))+labs(y="Spread [-]",x="Resolution [h]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  g.spread.simp <- ggplot(data=dat.simp, aes(x=reso, y=spread, shape=catchment, colour=errmod)) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(aes(y=spread.valid))+labs(y="Spread [-]",x="Resolution [h]")+theme_bw()
  g.spread.valid <- ggplot(data=dat, aes(x=reso, y=spread.valid)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=spread.valid, group=case, colour=errmod, linetype=catchment))+geom_blank(aes(y=spread))+labs(x="Resolution [h]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  g.spread.valid.simp <- ggplot(data=dat.simp, aes(x=reso, y=spread.valid, shape=catchment, colour=errmod)) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(aes(y=spread))+labs(x="Resolution [h]")+theme_bw()
  #vol.bias
  g.sferr <- ggplot(data=dat, aes(x=reso, y=sferr.med)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=sferr.med, group=case, colour=errmod, linetype=catchment))+geom_blank(aes(y=sferr.med.valid))+geom_hline(yintercept=0)+labs(x="Resolution [h]", y="Streamflow Error [%]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  g.sferr.simp <- ggplot(data=dat.simp, aes(x=reso, y=sferr.med, shape=catchment, colour=errmod)) + geom_hline(yintercept=0) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(aes(y=sferr.med.valid))+labs(x="Resolution [h]", y="Streamflow Error [%]")+theme_bw()
  g.sferr.valid <- ggplot(data=dat, aes(x=reso, y=sferr.med.valid)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=sferr.med.valid, group=case, colour=errmod, linetype=catchment))+geom_blank(aes(y=sferr.med))+geom_hline(yintercept=0)+labs(x="Resolution [h]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  g.sferr.valid.simp <- ggplot(data=dat.simp, aes(x=reso, y=sferr.med.valid, shape=catchment, colour=errmod)) +geom_hline(yintercept=0) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(aes(y=sferr.med))+labs(x="Resolution [h]")+theme_bw()
  ## Nash-Sutcliffe Efficiency
  g.nse <- ggplot(data=dat, aes(x=reso, y=nse.det)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=nse.det, group=case, colour=errmod, linetype=catchment))+geom_blank(mapping=aes(x=reso,y=nse.det.valid))+labs(shape="Catchment", colour="Error Model",linetype="Catchment",x="Resolution [h]",y=expression(widehat(E)[N*","*det]~" [-]"))+scale_colour_discrete(label=parse_format())+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  g.nse.simp <- ggplot(data=dat.simp, aes(x=reso, y=nse.det, shape=catchment, colour=errmod)) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(mapping=aes(x=reso,y=nse.det.valid))+labs(shape="Catchment", colour="Error Model",x="Resolution [h]",y=expression(widehat(E)[N*","*det]~" [-]"))+scale_colour_discrete(label=parse_format())+theme_bw()
  g.nse.valid <- ggplot(data=dat, aes(x=reso, y=nse.det.valid)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=nse.det.valid, group=case, colour=errmod, linetype=catchment))+geom_blank(mapping=aes(x=reso,y=nse.det))+labs(x="Resolution [h]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  g.nse.valid.simp <- ggplot(data=dat.simp, aes(x=reso, y=nse.det.valid, shape=catchment, colour=errmod)) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(mapping=aes(x=reso,y=nse.det))+labs(x="Resolution [h]")+theme_bw()
  ## flashiness index
  g.fi <- ggplot(data=dat, aes(x=reso, y=flash.obs-flash.med)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=flash.obs-flash.med,group=case,colour=errmod, linetype=catchment))+geom_blank(aes(x=reso,y=flash.obs.valid-flash.med.valid))+labs(shape="Catchment", linetype="Catchment", colour="Error Model")+geom_hline(yintercept=0)+scale_colour_discrete(label=parse_format())+labs(x="Resolution [h]", y=expression(I[F*","*obs]-~"median("~I[F]~")"))+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  g.fi.simp <- ggplot(data=dat.simp, aes(x=reso, y=flash.obs-flash.med, shape=catchment, colour=errmod)) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(aes(x=reso,y=flash.obs.valid-flash.med.valid))+labs(shape="Catchment", colour="Error Model")+geom_hline(yintercept=0)+scale_colour_discrete(label=parse_format())+labs(x="Resolution [h]", y=expression(I[F*","*obs]-~"median("~I[F]~")"))+theme_bw()
  g.fi.valid <- ggplot(data=dat, aes(x=reso, y=flash.obs.valid-flash.med.valid)) + geom_point(aes(shape=catchment,colour=errmod),size=1.3,alpha=1)+geom_path(aes(x=reso,y=flash.obs.valid-flash.med.valid,group=case,colour=errmod, linetype=catchment))+geom_blank(aes(x=reso,y=flash.obs-flash.med))+geom_hline(yintercept=0)+labs(x="Resolution [h]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  g.fi.valid.simp <- ggplot(data=dat.simp, aes(x=reso, y=flash.obs.valid-flash.med.valid, shape=catchment, colour=errmod)) + geom_jitter(width=0.2,height=0,size=1.3,alpha=1)+geom_blank(aes(x=reso,y=flash.obs-flash.med))+geom_hline(yintercept=0)+labs(x="Resolution [h]")+theme_bw()
  # put together flashiness index, reliability and spread of reduced plot
  leg <- get_legend(g.fi.simp+theme(legend.position="right",legend.margin=margin(t=8,unit="pt"),text=element_text(size=14)))
  leg <- plot_grid(leg,NULL,NULL, nrow=3, labels=c("","",""))
  pg <- plot_grid(g.fi.simp+theme(legend.position="none",plot.margin=unit(c(30,0,5.5,5.5),"pt"),text=element_text(size=14),axis.title.x=notext,axis.text.x=notext,axis.ticks.x=notext),g.fi.valid.simp+theme(legend.position="none",plot.margin=unit(c(30,5.5,5.5,0),"pt"),text=element_text(size=14),axis.title=notext,axis.text.x=notext,axis.ticks.x=notext),
                  g.reli.simp+theme(legend.position="none",text=element_text(size=14),axis.title.x=notext,axis.text.x=notext,axis.ticks.x=notext,plot.margin=unit(c(5.5,0,5.5,0),"pt")),g.reli.valid.simp+theme(legend.position="none",text=element_text(size=14),plot.margin=unit(c(5.5,5.5,5.5,0),"pt"),axis.title=notext,axis.text.x=notext,axis.ticks.x=notext),
                  g.spread.simp+theme(legend.position="none",text=element_text(size=14),plot.margin=unit(c(5.5,0,5.5,0),"pt")),g.spread.valid.simp+theme(legend.position="none",text=element_text(size=14),plot.margin=unit(c(5.5,5.5,5.5,0),"pt"),axis.title.y=notext), ncol=2, labels=c("Calibration","Validation","","","",""), rel_heights=c(1.2,1,1), label_x=c(0.25,0.28), align="v")
  save_plot(paste0(outpath,"plot_results1_simp.pdf"), plot_grid(pg, leg, ncol=2, rel_widths=c(1,0.25)), base_height=7.5,base_aspect_ratio=1.0)
  dev.off()
  # put together flashiness index, reliability and spread of full plot
  leg <- get_legend(g.fi+theme(legend.position="right",legend.margin=margin(t=8,unit="pt"),text=element_text(size=14)))
  better1 <- ggdraw() + draw_image("sudriv_output/bettersvg.svg",y=0.2,width=0.5,height=0.5)
  better2 <- ggdraw() + draw_image("sudriv_output/bettersvg.svg",y=0.3,width=0.5,height=0.5)
  leg <- plot_grid(leg,better1,better2, nrow=3, labels=c("","",""))
  pg <- plot_grid(g.fi+theme(legend.position="none",plot.margin=unit(c(30,0,5.5,5.5),"pt"),text=element_text(size=14),axis.title.x=notext,axis.text.x=notext,axis.ticks.x=notext),g.fi.valid+theme(legend.position="none",plot.margin=unit(c(30,5.5,5.5,0),"pt"),text=element_text(size=14),axis.title=notext,axis.text.x=notext,axis.ticks.x=notext),
                  g.reli+theme(legend.position="none",text=element_text(size=14),axis.title.x=notext,axis.text.x=notext,axis.ticks.x=notext,plot.margin=unit(c(5.5,0,5.5,0),"pt")),g.reli.valid+theme(legend.position="none",plot.margin=unit(c(5.5,5.5,5.5,0),"pt"),text=element_text(size=14),axis.title=notext,axis.text.x=notext,axis.ticks.x=notext),
                  g.spread+theme(legend.position="none",text=element_text(size=14),plot.margin=unit(c(5.5,0,5.5,0),"pt")),g.spread.valid+theme(legend.position="none",plot.margin=unit(c(5.5,5.5,5.5,0),"pt"),text=element_text(size=14),axis.title.y=notext), ncol=2, labels=c("Calibration","Validation","","","",""), rel_heights=c(1.2,1,1), label_x=c(0.25,0.28), align="v")
  save_plot(paste0(outpath,"plot_results1.pdf"), plot_grid(pg, leg, ncol=2, rel_widths=c(1,0.25)), base_height=7.5,base_aspect_ratio=1.0)
  dev.off()
  ## put together vol. bias and NSE of full plot
  leg <- get_legend(g.nse+theme(legend.position="right",legend.margin=margin(t=8,unit="pt"),text=element_text(size=14)))
  leg <- plot_grid(leg,NULL,NULL, nrow=3, labels=c("","",""))
  pg <- plot_grid(g.sferr+theme(legend.position="none",plot.margin=unit(c(30,0,5.5,5.5),"pt"),text=element_text(size=14),axis.title.x=notext,axis.text.x=notext,axis.ticks.x=notext),g.sferr.valid+theme(legend.position="none",plot.margin=unit(c(30,5.5,5.5,5.5),"pt"),text=element_text(size=14),axis.title=notext,axis.text.x=notext,axis.ticks.x=notext),
                  g.nse+theme(legend.position="none",plot.margin=unit(c(5.5,0,5.5,0),"pt"),text=element_text(size=14)),g.nse.valid+theme(legend.position="none",text=element_text(size=14),axis.title.y=notext), ncol=2, labels=c("Calibration","Validation","",""), rel_heights=c(1.2,1), label_x=c(0.25,0.28), align="v")
  save_plot(paste0(outpath,"plot_results2.pdf"), plot_grid(pg, leg, ncol=2, rel_widths=c(1,0.25)), base_height=7.5,base_aspect_ratio=1.0)
  dev.off()
  ## put together vol. bias and NSE of reduced plot
  leg <- get_legend(g.nse.simp+theme(legend.position="right",legend.margin=margin(t=8,unit="pt"),text=element_text(size=14)))
  leg <- plot_grid(leg,NULL,NULL, nrow=3, labels=c("","",""))
  pg <- plot_grid(g.sferr.simp+theme(legend.position="none",plot.margin=unit(c(30,0,5.5,5.5),"pt"),text=element_text(size=14),axis.title.x=notext,axis.text.x=notext,axis.ticks.x=notext),g.sferr.valid.simp+theme(legend.position="none",plot.margin=unit(c(30,5.5,5.5,5.5),"pt"),text=element_text(size=14),axis.title=notext,axis.text.x=notext,axis.ticks.x=notext),
                  g.nse.simp+theme(legend.position="none",plot.margin=unit(c(5.5,0,5.5,0),"pt"),text=element_text(size=14)),g.nse.valid.simp+theme(legend.position="none",text=element_text(size=14),axis.title.y=notext), ncol=2, labels=c("Calibration","Validation","",""), rel_heights=c(1.2,1), label_x=c(0.25,0.28), align="v")
  save_plot(paste0(outpath,"plot_results2_simp.pdf"), plot_grid(pg, leg, ncol=2, rel_widths=c(1,0.25)), base_height=7.5,base_aspect_ratio=1.0)
  dev.off()
}
plot.KL.summary <- function(file=NA,outpath="sudriv_output/"){
  dat <- read.table(file, sep=" ", header=TRUE)
  dat[dat[,"catchment"]=="waengi","catchment"] <- "Murg"
  dat[dat[,"catchment"]=="maimai","catchment"] <- "Maimai"
  ## decide where to use the smoothed version and where the original one
  ## dat <- dat[-which(dat$catchment=="Maimai" & dat$reso=="1h" & dat$errmod %in% c("E3", "E4")),]
  ## dat <- dat[-which(dat$catchment=="Murg" & dat$reso=="1h" & dat$errmod =="E3aP"),]
  dat$reso <- gsub("h", "", dat$reso)
  xx <- dat$reso
  xx[xx=="6"] <- 2
  xx[xx=="24"] <- 3
  xx <- as.numeric(xx)
  ##dat$reso <- factor(xx, levels=c(1,2,3))
  dat$reso <- reorder(as.factor(dat$reso),X=xx)
  dat$errmod <- gsub("P", "", dat$errmod)
  dat$errmod <- gsub("mean", "", dat$errmod)
  dat$errmod[dat$errmod=="E3"] <- "E3(\u002A)"
  dat$errmod[dat$errmod=="E3a"] <- "E3a(\u002A)"
  dat$errmod[dat$errmod=="E4"] <- "E4(\u002A)"
  dat$errmod[dat$errmod=="E4a"] <- "E4a(\u002A)"
  xx <- dat$errmod
  xx[grepl("E1",xx)] <- 1
  xx[grepl("E2",xx)] <- 2
  xx[grepl("E3",xx)] <- 3
  xx[grepl("E3a",xx)] <- 4
  xx[grepl("E4",xx)] <- 5
  xx[grepl("E4a",xx)] <- 6
  xx <- as.numeric(xx)
  dat$errmod <- reorder(as.factor(dat$errmod),X=xx)
  dat$case <- paste(as.character(dat$catchment), as.character(dat$errmod))
  ## dat$meas <- 1-dat$reli^(1-dat$prec)
  ## dat$meas.valid <- 1-dat$reli.valid^(1-dat$prec.valid)
  notext <- element_blank()
  ## colours
  my_palette = c("#000000", brewer.pal(3, "Set1")[1], brewer.pal(9, "Greens")[c(7,5)], brewer.pal(9, "Purples")[c(7,4)])
  scale_colour_discrete = function(...) scale_colour_manual(..., values = palette())
  palette(my_palette)
  ## Cmlt_E
  g.Cmlt_E <- ggplot(data=dat, aes(x=reso, y=U1W_Cmlt_E)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=U1W_Cmlt_E,colour=errmod,linetype=catchment,group=case))+labs(shape="Catchment",linetype="Catchment",colour="Error Model",title=expression(C[E]),y="KL-Divergence [-]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  # Smax_UR
  g.Smax_UR <- ggplot(data=dat, aes(x=reso, y=U1W_Smax_UR)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=U1W_Smax_UR,colour=errmod,linetype=catchment,group=case))+labs(title=expression(S[max]))+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  # U1W_K_Qb_UR
  g.K_Qb_UR <- ggplot(data=dat, aes(x=reso, y=U1W_K_Qb_UR)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=U1W_K_Qb_UR,colour=errmod,linetype=catchment,group=case))+labs(title=expression(k[u]),y="KL-Divergence [-]",x="Resolution [h]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  # U1W_K_Qq_FR
  g.K_Qq_FR <- ggplot(data=dat, aes(x=reso, y=U1W_K_Qq_FR)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=U1W_K_Qq_FR,colour=errmod,linetype=catchment,group=case))+labs(title=expression(k[f]),x="Resolution [h]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  # C1Wv_Qstream_a_lik
  g.a_lik<- ggplot(data=dat, aes(x=reso, y=C1Wv_Qstream_a_lik)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=C1Wv_Qstream_a_lik,colour=errmod,linetype=catchment,group=case))+labs(shape="Catchment", linetype="Catchment", colour="Error Model",title="a",y="KL-Divergence [-]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  # C1Wv_Qstream_b_lik
  g.b_lik<- ggplot(data=dat, aes(x=reso, y=C1Wv_Qstream_b_lik)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=C1Wv_Qstream_b_lik,colour=errmod,linetype=catchment,group=case))+labs(title="b")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  # C1Wv_Qstream_taumax_lik
  g.taumin_lik<- ggplot(data=dat, aes(x=reso, y=C1Wv_Qstream_taumin_lik)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=C1Wv_Qstream_taumin_lik,colour=errmod,linetype=catchment,group=case))+labs(title=expression(tau[min]),y="KL-Divergence [-]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  # C1Wv_Qstream_df_lik
  g.taumax_lik<- ggplot(data=dat, aes(x=reso, y=C1Wv_Qstream_taumax_lik)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=C1Wv_Qstream_taumax_lik,colour=errmod,linetype=catchment,group=case))+labs(title=expression(tau[max]))+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  # C1Wv_Qstream_gamma_lik
  g.gamma_lik<- ggplot(data=dat, aes(x=reso, y=C1Wv_Qstream_gamma_lik)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=C1Wv_Qstream_gamma_lik,colour=errmod,linetype=catchment,group=case))+labs(title=expression(gamma),y="KL-Divergence [-]",x="Resolution [h]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  # C1Wv_Qstream_taumin_lik
  g.df_lik<- ggplot(data=dat, aes(x=reso, y=C1Wv_Qstream_df_lik)) + geom_point(aes(shape=catchment,colour=errmod),size=1,alpha=1)+geom_path(aes(x=reso,y=C1Wv_Qstream_df_lik,colour=errmod,linetype=catchment,group=case))+labs(title="df",x="Resolution [h]")+scale_x_discrete(expand=c(0.1,0.1))+theme_bw()
  ## put together hydrological model parameters
  leg <- get_legend(g.Cmlt_E+theme(legend.position="right"))
  pg <- plot_grid(g.Cmlt_E+theme(legend.position="none",plot.margin=unit(c(5.5,0,5.5,5.5),"pt"),axis.title.x=notext),g.Smax_UR+theme(legend.position="none",plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt"),axis.title.x=notext,axis.title.y=notext),
                  g.K_Qb_UR+theme(legend.position="none",plot.margin=unit(c(5.5,0,5.5,0),"pt")),g.K_Qq_FR+theme(legend.position="none",plot.margin=unit(c(5.5,0,5.5,0),"pt"),axis.title.y=notext), ncol=2, align="v")
  save_plot(paste0(outpath,"plot_KL1.pdf"), plot_grid(pg, leg, ncol=2, rel_widths=c(1,0.25)), base_height=6,base_aspect_ratio=0.7)
  dev.off()
  ## put together the error model parameters
  leg <- get_legend(g.a_lik+theme(legend.position="right"))
  pg <- plot_grid(g.a_lik+theme(legend.position="none",plot.margin=unit(c(5.5,0,5.5,5.5),"pt"),axis.title.x=notext),g.b_lik+theme(legend.position="none",plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt"),axis.title.x=notext,axis.title.y=notext),
                  g.taumin_lik+theme(legend.position="none",axis.title.x=notext,plot.margin=unit(c(5.5,0,5.5,0),"pt")),g.taumax_lik+theme(legend.position="none",axis.title.x=notext,axis.title.y=notext),
                  g.gamma_lik+theme(legend.position="none",plot.margin=unit(c(5.5,0,5.5,0),"pt")),g.df_lik+theme(legend.position="none",axis.title.y=notext), ncol=2, align="v")
  save_plot(paste0(outpath,"plot_KL2.pdf"), plot_grid(pg, leg, ncol=2, rel_widths=c(1,0.25)), base_height=6,base_aspect_ratio=0.7)
  dev.off()
}
plot.dens.par <- function(list.su, brn.ins, covariates=NA, pars="C1Wv_Qstream_taumin_lik"){
  data <- data.frame()
  if(!is.na(covariates[1])){if(length(covariates) != length(list.su)){warning("length of covariates and list.su differ"); return(NA)}}
  j <- 1
  for(su.curr in list.su){
    for(par.curr in pars){
      value <- c(su.curr$parameter.sample[(brn.ins[j]):(dim(su.curr$parameter.sample)[1]),par.curr,])
      ## back-transform parameters to original scale
      trans <- par.curr %in% c(names(su.curr$model$parameters)[as.logical(su.curr$model$args$parTran)], names(su.curr$likelihood$parameters)[as.logical(su.curr$likelihood$tran)])
      if(trans) value <- exp(value)
      ## scale model parameters to the desired time format
      if(!is.null(su.curr$model$par.time)){
        if(par.curr %in% names(su.curr$model$parameters)[su.curr$model$par.time==1]) value <- value*su.curr$layout$timestep.fac
        if(par.curr %in% names(su.curr$model$parameters)[su.curr$model$par.time==-1]) value <- value/su.curr$layout$timestep.fac
      }
      ## scale likelihood parameters to the desired time format
      if(!is.null(su.curr$likelihood$par.time)){
        if(par.curr %in% names(su.curr$likelihood$parameters)[su.curr$likelihood$par.time==1]) value <- value*su.curr$layout$timestep.fac
        if(par.curr %in% names(su.curr$likelihood$parameters)[su.curr$likelihood$par.time==-1]) value <- value/su.curr$layout$timestep.fac
      }
      data <- rbind(data, data.frame(value=value, model=names(list.su)[[j]],par=par.curr,covariates=ifelse(is.na(covariates[1]),NULL,covariates[j])))
    }
    j <- j + 1
  }
  print(unique(data$par))
  data$par[which(data$par == "C1Wv_Qstream_taumin_lik")] <- "taumin"
  data$par[which(data$par == "C1Wv_Qstream_taumax_lik")] <- "taumax"
  data$model <- gsub("maimai.h1", "Maimai 1 h", data$model)
  data$model <- gsub("maimai.h6", "Maimai 6 h", data$model)
  data$model <- gsub("maimai.h24", "Maimai 24 h", data$model)
  data$model <- gsub("murg.h1", "Murg 1 h", data$model)
  data$model <- gsub("murg.h6", "Murg 6 h", data$model)
  data$model <- gsub("murg.h24", "Murg 24 h", data$model)
  xx <- data$model
  xx[grep("Murg 1 h",xx)] <- 1
  xx[grep("Murg 6 h",xx)] <- 2
  xx[grep("Murg 24 h",xx)] <- 3
  xx[grep("Maimai 1 h",xx)] <- 4
  xx[grep("Maimai 6 h",xx)] <- 5
  xx[grep("Maimai 24 h",xx)] <- 6
  xx <- as.numeric(xx)
  data$model <- reorder(as.factor(data$model),X=xx)
  ## g.obj <- ggplot(data, aes(x=value, fill=model)) + geom_density(alpha=0.5) + labs(x=expression(paste(tau[min]," [h]")), y="Density", fill="Temp. resolution") + theme(legend.position=c(0.7,0.7))
  g.obj <- ggplot(data,aes(x=value,y=model)) + geom_density_ridges(aes(fill=paste(model,par),alpha=covariates), scale=2)  + scale_y_discrete(expand=c(0.01,0)) + scale_fill_cyclical(breaks=c("Murg 1 h taumin","Murg 1 h taumax"), labels=c(`Murg 1 h taumin`=expression(tau[min]), `Murg 1 h taumax`=expression(tau[max])), values=c("#ff0000", "#0000ff"), name="Parameter", guide="legend") + theme_ridges(font_size=16) + scale_x_log10(breaks=c(20,50,100,200,500,1000), labels=c(20,50,100,200,500,1000), expand=c(0.01,0)) + labs(x="Value [h]", y="Case", alpha=expression(Xi[reli]))
  pdf("sudriv_output/taumin_timeres2.pdf", height=6)
  plot(g.obj)
  dev.off()
}
plot.eta.Qdet.cor <- function(list.su, tme.orig, brn.in=0, ylim=NA, outpath=""){
  cor1 <- plot.cor(sudriv=list.su[[1]], brn.in=brn.in, plot=FALSE)
  cor2 <- plot.cor(sudriv=list.su[[2]], brn.in=brn.in, plot=FALSE)
  dat1 <- pred.stats(list.su[1],tme.orig=tme.orig)
  dat2 <- pred.stats(list.su[2],tme.orig=tme.orig)
  eta1 <- plot.Qdet.quantiles(dat1, list.su[[1]], ind.sel=list.su[[1]]$layout$pred, plot=FALSE)
  eta2 <- plot.Qdet.quantiles(dat2, list.su[[2]], ind.sel=list.su[[2]]$layout$pred, plot=FALSE)
  if(!is.na(ylim[1])) eta1 <- eta1 + coord_cartesian(ylim=ylim)
  if(!is.na(ylim[1])) eta2 <- eta2 + coord_cartesian(ylim=ylim)
  notext <- element_blank()
  pg <- plot_grid(eta1 + theme(plot.margin=unit(c(40,10,5.5,5.5),"pt")),eta2 + theme(axis.title.y=notext, plot.margin=unit(c(40,0,5.5,10),"pt")),ggmatrix_gtable(cor1),ggmatrix_gtable(cor2),ncol=2,rel_heights=c(0.8,1),labels=c("E2", "E3", "", ""), label_size=20, label_x=c(0.5,0.5,0,0))
  save_plot(paste0(outpath,"figure6.png"), pg, base_height=12, base_aspect_ratio=1.3)
  dev.off()
}
plot.ts.eta <- function(su.daily, su.hourly, dat.daily, dat.hourly, tme.orig, ylim=NA, outpath=""){
  ts1 <- plot.ts.quantiles(dat.daily, su.daily, xlim="calib", precip=TRUE, plot=FALSE)
  ts2 <- plot.ts.quantiles(dat.hourly, su.hourly, xlim=c(4000,6300), precip=TRUE, plot=FALSE)
  if(!is.na(ylim[1])) ts1 <- ts1 + coord_cartesian(ylim=ylim)
  if(!is.na(ylim[1])) ts2 <- ts2 + coord_cartesian(ylim=ylim)
  notext <- element_blank()
  pg <- plot_grid(ts1+ theme(plot.margin=unit(c(25,20,5.5,5.5),"pt")), ts2 + theme(plot.margin=unit(c(25,20,5.5,5.5),"pt")),ncol=1,labels=c("Daily data", "Hourly data"), label_size=14, label_x=c(0.3,0.3))
  save_plot(paste0(outpath,"figure92.png"), pg, base_height=5, base_aspect_ratio=1.8)
  dev.off()
}
plot.white.noise.paper <- function(list.su, tme.orig, brn.in=0, ylim=NA, outpath=""){
  ## we assume that list.su contains "maimai1hE2", "maimai1hE3P", "waengiE2" and "waengiE3"
  ## and tme.orig is vector of length 2 for maimai1h2 and waengi
  dat1 <- pred.stats(list.su[1],tme.orig=tme.orig[1])
  dat2 <- pred.stats(list.su[2],tme.orig=tme.orig[1])
  dat1$precip <- list.su[[2]]$input$P.roll
  dat2$precip <- list.su[[2]]$input$P.roll
  dat.maimai <- rbind(dat1,dat2)
  dat1 <- pred.stats(list.su[3],tme.orig=tme.orig[2])
  dat2 <- pred.stats(list.su[4],tme.orig=tme.orig[2])
  dat1$precip <- list.su[[4]]$input$P.roll
  dat2$precip <- list.su[[4]]$input$P.roll
  dat.waengi <- rbind(dat1,dat2)
  noi.maimai <- plot.ts.white.noise(dat.maimai, list.su[[1]], xlim=c(11490,11778), precip=TRUE, plot=FALSE)
  noi.waengi <- plot.ts.white.noise(dat.waengi, list.su[[3]], xlim=c(18170,18350), precip=TRUE, plot=FALSE)
  ##qq.maimai <- ggplot(dat.maimai, aes(sample=white.noise, colour=case, shape=case)) + stat_qq() + stat_qq_line()
  ##qq.waengi <- ggplot(dat.waengi, aes(sample=white.noise, colour=case, shape=case)) + stat_qq() + stat_qq_line()
  pg <- plot_grid(noi.maimai + theme(plot.margin=unit(c(0.5,0.5,0.2,0.5),"in")),noi.waengi + theme(plot.margin=unit(c(0.5,0.5,0.2,0.5),"in")),ncol=1,labels=c("Maimai", "Murg"), label_size=18)
  save_plot(paste0(outpath,"figure_white_noise.pdf"), pg, base_height=12, base_aspect_ratio=1.3)
  dev.off()
}
compare.logpdfs <- function(lgs, file="plot_logpdfs.png"){
  ## This function compares and plots the logpdfs reached with different time dependent parameters.
  ## Input is a table of logpdfs reached for different parameters.
  lgs[,"var"] <- sapply(lgs[,"var"], function(x) strsplit(x, "_")[[1]][1])
  gg <- ggplot(data = lgs) + theme_bw() 
  gg.lik   <- gg + geom_point(aes(y=var,x=loglikeliobs)) + geom_vline(xintercept=as.numeric(lgs%>%filter(grepl("none",var))%>%select(loglikeliobs)%>%summarise(mn=mean(loglikeliobs)))) + labs(y="Parameter", x="Observational likelihood")
  gg.post  <- gg + geom_point(aes(y=var,x=logposterior)) + geom_vline(xintercept=as.numeric(lgs%>%filter(grepl("none",var))%>%select(logposterior)%>%summarise(mn=mean(logposterior)))) + labs(y="Parameter", x="Posterior")
  gg.ou    <- gg + geom_point(aes(y=var,x=logpdfou_timedeppar)) + labs(y="Parameter", x="logpdf of time-course")
  #pg <- plot_grid(gg.lik, gg.post, gg.ou, nrow = 3)
  save_plot(file, plot = gg.lik, base_height = 7, base_asp = 0.7, dpi=500)
}
plot.loess.scatter <- function(list.su, plot.which, mfrow=c(1,2), args.ggsave, ...){
  #pdf(...)
  #par(mfrow=mfrow, mar=c(5,5,1,2))
  j <- 1
  gg.list <- list()
  for(feat in names(plot.which)){
    tag.red <- grep(plot.which[j], names(list.su), value=TRUE)
    xlab <- mylabeller.feat.units(feat)
    ylab <- mylabeller.param.units(tag.red, log=TRUE)
    ## prepare data
    if(grepl("\\+",feat)){
      nms <- names(list.su[[tag.red]]$model$timedep$model.td.Nd)
      if(feat %in% nms){
        sm <- list.su[[tag.red]]$model$timedep$model.td.Nd[[feat]]
      }else{ # try to flip features around, maybe they were specified in the wrong order...
        tmp <- strsplit(feat, "\\+")[[1]]
        feat <- paste0(tmp[2:1], collapse="+")
        sm <- list.su[[tag.red]]$model$timedep$model.td.Nd[[feat]]
      }
      xlab <- mylabeller.feat.units(strsplit(feat, "\\+")[[1]])
    }else{
      sm <- list.su[[tag.red]]$model$timedep$model.td[[feat]]
    }
    if(is.null(sm)){
      warning(paste0("no model found for ", tag.red, " and ", feat))
    }else{
      r2 <- 1-sum((sm$residuals)^2)/sum((sm$y-mean(sm$y))^2)
      if(grepl("\\+",feat)){ # contour plot
        print(feat)
        x1 <- seq(min(sm$x[,1]), max(sm$x[,1]), length.out = 500)
        x2 <- seq(min(sm$x[,2]), max(sm$x[,2]), length.out = 500)
        dat <- expand.grid(x1, x2, KEEP.OUT.ATTRS = FALSE)
        colnames(dat) <- colnames(sm$x)
        y.pred <- predict(sm, newdata=dat)
        dat <- cbind(dat, y=y.pred)
        colnames(dat)[3] <- tag.red
        nms <- colnames(dat)
        first = sym(nms[1])
        second = sym(nms[2])
        third = sym(nms[3])
        gg.list[[j]] <- ggplot(data=dat, aes(x=!!first, y=!!second)) + 
          stat_contour(geom="polygon", aes(z=!!third, fill=..level..)) + geom_tile(aes(fill=!!third)) + stat_contour(aes(z=!!third), bins=15) + 
          scale_fill_continuous(..., type="viridis", guide=guide_colorbar(title=ylab)) + geom_point(data = as.data.frame(sm$x), shape=20, color="red", alpha=0.1)+
          labs(x=xlab[1], y=xlab[2], fill=ylab, subtitle=bquote(R^2 == .(round(r2,2))))
      }else{ # scatterplot
        #https://www.inwt-statistics.com/read-blog/smoothscatter-with-ggplot2-513.html
        dat <- as.data.frame(cbind(sm$x, sm$y, fitted=sm$fitted))
        nms <- colnames(dat)
        first <-  sym(nms[1])
        second <- sym(nms[2])
        gg.list[[j]] <- ggplot(data = dat, aes(x=!!first, y=!!second)) +
          stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200, show.legend=FALSE) +
          #scale_fill_continuous(low = "white", high = "dodgerblue4")+
          scale_fill_viridis(direction=-1, option="magma", begin=0.5)+
          #geom_point(color="red", alpha=0.1, shape=20)+
          geom_line(aes(x=!!first, y=fitted))+
          annotate("text", x=max(dat[,nms[1]]), y=max(dat[,nms[2]]), label=bquote(R^2 == .(round(r2,2))), hjust=1, vjust=1)+
          labs(x=xlab, y=ylab)
        #smoothScatter(x=sm$x, y=sm$y, nrpoints=1000, xlab=xlab, ylab=ylab)
        #lines(x=sm$x, y=sm$fitted, col="red")
        #title(sub=bquote(R^2 == .(round(r2, 2))))
      }
    }
    j <- j + 1
  }
  last <- egg::ggarrange(plots=gg.list, nrow=2, labels=paste0("(",letters[1:length(gg.list)], ")"))
  args.ggsave <- c(args.ggsave, list(plot=last))
  do.call(ggsave, args = args.ggsave)
  dev.off()
}