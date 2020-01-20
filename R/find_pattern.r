find.pattern.timedep <- function(sudriv, vars=NULL, validation_split=0.2, add.data=NULL, tag=""){
  ## This function compares the time course of the time dependent parameters to the model states, output (and potentially other variables) and identifies matching patterns.
  ## consistency checks:
  print(tag)
  tag.red <- gsub("_.*","",tag)
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
  if(tag.red %in% c("kdwr","rswr","sloneir","sltwoir","alqqfr")){# if it is a chemistry related parameter
    strt <- sudriv$layout$layout %>% slice(sudriv$layout$calib) %>% filter(var %in% c("C1Tc1_Qstream","C1Tc2_Qstream")) %>% select(time) %>% min
    end <- sudriv$layout$layout %>% slice(sudriv$layout$calib) %>% filter(var %in% c("C1Tc1_Qstream","C1Tc2_Qstream")) %>% select(time) %>% max
  }else{ #if it is a more water related parameter
    strt <- min(sudriv$layout$layout$time[sudriv$layout$calib])
    end <- max(sudriv$layout$layout$time[sudriv$layout$calib])
  }
  cat("strt: ",strt,"\n")
  cat("end: ",end,"\n")
  y.all2 <- y.all2 %>% filter(time >= strt & time <= end) %>% na.omit
  cat("dim data:\t",dim(y.all2),"\n")
  ## add some features
  y.all2[paste0(colnames(y.all2),"2")] <- y.all2^2
  #y.all2[paste0(colnames(y.all),"3")] <- sqrt(y.all)
  y.all2 <- y.all2 %>% select(-time2, -y.td2)
  ## lm2 <- lm(y.td ~ ., data=y.all2)
  if(sudriv$model$args$parTran[which(sudriv$model$timedep$pTimedep)[1]] == 1){
    y.all2 <- y.all2 %>% mutate(y.td=exp(y.td))
  }
  train <- (1:nrow(y.all2))[1:round(nrow(y.all2)*(1-validation_split))] ## train in beginning, test at end
  test <- (1:nrow(y.all2))[-train]
  test2 <- (1:nrow(y.all2))[1:round(nrow(y.all2)*validation_split)] ## train at end, test in beginning
  train2 <- (1:nrow(y.all2))[-test2]
  train3 <- (1:nrow(y.all2))[c(1:round(nrow(y.all2)*(0.5-0.5*validation_split)),round(nrow(y.all2)*(0.5+0.5*validation_split)):nrow(y.all2))] ## train at beginning and end, test in the middle
  test3 <- (1:nrow(y.all2))[-train3]
  
  ## =================================================================================================
  ## scatterplot between timedep par and explanatory variables
  pdf(paste0("../output/timedeppar/A1Str07h2x/",tag,"/plot_scatter.pdf"))
  mapply(function(x,y,nm,tag){
    dat <- data.frame(x=x,y=y) %>% arrange(x)
    smoothScatter(x=dat$x,y=dat$y,main=paste(tag,"&",nm),xlab=nm,ylab=tag.red,nrpoints=1000)
    sm <- loess(y~x,data=dat,span=1/4)
    pred <- predict(sm, newdata=dat)
    lines(x=dat$x, y=pred, col="red")
    title(sub=bquote(R^2 == .(round(1-var(pred-dat$y)/var(dat$y), 2))))
  }, y.all2, nm=colnames(y.all2), MoreArgs=list(y=y.all2[,"y.td"], tag=tag.red))
  dev.off()
  
  ## =================================================================================================
  ## scale the data for fitting glms later on
  print("y.all2:")
  print(summary(y.all2))
  boxes <- apply(y.all2, 2, function(x){
    print(summary(x))
    if(all(x>=0)){
      res <- tryCatch({
        unlist(boxcoxfit(x,lambda2=TRUE)[c("lambda","beta.normal","sigmasq.normal")])
      },error=function(cond){
        cond
      },warning={}
      )
      if(inherits(res,"error")){
        return(c(NA,NA,NA,NA))
      }else{
        return(res)
      }
    }else{
      return(c(NA,NA,NA,NA))
    }
  }
  )
  y.all2scaled <- y.all2
  for(i in 1:ncol(y.all2)){
    if(!(colnames(y.all2)[i]=="y.td")){
      if(any(is.na(boxes[,i]))) next
      if(boxes[1,i]!=0){
        y.all2scaled[,i] <- (((y.all2scaled[,i]+boxes[2,i])^boxes[1,i] - 1)/boxes[1,i] - boxes[3,i])/sqrt(boxes[4,i])
      }else{
        y.all2scaled[,i] <- (log(y.all2scaled[,i]+boxes[2,i]) - boxes[3,i])/sqrt(boxes[4,i])
      }
    }
  }
  
  ## =================================================================================================
  ## clustering analysis
  # clust.dat <- y.all2scaled
  # clust.dat <- clust.dat[,which(!grepl("2",colnames(clust.dat)))]
  # clust.dat <- as.data.frame(scale(clust.dat))
  # save(clust.dat,file = paste0("../output/timedeppar/A1Str07h2x/",tag,"/tddata.RData"))
  # print(summary(clust.dat))
  # print("starting clustering analysis ...")
  # clust <- dbscan(y.all2scaled %>% select(-time), eps=1, minPts=15)
  # print(clust)
  # plt <- fviz_cluster(clust, data = y.all2scaled %>% select(-time), stand = FALSE, ellipse = FALSE, show.clust.cent = FALSE, geom = "point",palette = "jco", ggtheme = theme_classic())
  # ggsave(filename = paste0("../output/timedeppar/A1Str07h2x/",tag,"/plot_cluster.png"))
  # print("done.")
  
  ## =================================================================================================
  ## look at some additional things
  cors <- apply(y.all2scaled, 2, cor, y=y.all2scaled[,"y.td"])
  cat("correlations:\n")
  print(cors)
  ## test for stationarity
  stati <- apply(y.all2scaled, 2, function(x) adf.test(x)$p.value)
  names(stati) <- colnames(y.all2scaled)
  cat("p-values under null hypothesis of non-stationarity:\n")
  print(stati)
  pdf(paste0("../output/timedeppar/A1Str07h2x/",tag,"/plot_crosscorr.pdf"))
  par(mfrow=c(3,3))
  mapply(function(y,x,lag.max,nm,plot,tag) ccf(x=y,y=x,lag.max=lag.max,plot=plot,main=paste(tag.red,"&",nm),ylim=c(-0.6,0.6)), y.all2scaled, nm=colnames(y.all2scaled), MoreArgs=list(y=y.all2scaled[,"y.td"], lag.max=2*7*24*4, plot=TRUE, tag=tag)) ## 2 weeks max lag
  dev.off()
  
  # ================================================================================
  ## fit linear models
  sudriv$model$timedep$empir.model$glm <- NULL
  if(is.null(sudriv$model$timedep$empir.model$glm)){
    hlf <- ncol(as.data.frame(y.all2scaled)%>%select(-time,-y.td))%/%3
    tmp <- colnames(as.data.frame(y.all2scaled)%>%select(-time,-y.td))[1:hlf]
    frm <- as.formula(paste0("y.td ~ ",paste(tmp,collapse="+")))
    print(paste0("y.td ~ ",paste(tmp,collapse="+")))
    cf <- tryCatch(
      {coef(glm(frm, family=Gamma, data=as.data.frame(y.all2scaled)[train,]%>%select(-time), na.action="na.exclude")) ## get coefficients of first half of columns as starting coefficients
      },error=function(cond){message("fitting glm failed:");message(cond);return(NULL)})
    if(!is.null(cf)){
      glm2 <- glm(y.td ~ ., family=Gamma, data=as.data.frame(y.all2scaled)[train,]%>%select(-time), na.action="na.exclude", start=c(cf,rep(0,ncol(y.all2scaled)-2-hlf)))
    }
    tmp <- colnames(as.data.frame(y.all2scaled)%>%select(-time,-y.td))
    tmp <- tmp[!grepl("2", tmp)]
    frm <- as.formula(paste0("y.td~",paste(tmp,collapse="*")))
    cat("only linear terms formula: ", paste0("y.td~",paste(tmp,collapse="*")),"\n")
    linmod  <- lm(formula=frm, data=as.data.frame(y.all2scaled)[train,]%>%select(-time), na.action="na.exclude")
    linmod2 <- lm(formula=frm, data=as.data.frame(y.all2scaled)[train2,]%>%select(-time), na.action="na.exclude")
    linmod3 <- lm(formula=frm, data=as.data.frame(y.all2scaled)[train3,]%>%select(-time), na.action="na.exclude")
  }else{
    glm2 <- sudriv$model$timedep$empir.model$glm
  }
  ## if(is.null(sudriv$model$timedep$empir.model$nn)){
  ##     nn   <- neuralnet(y.td ~ ., data=as.data.frame(y.all2scaled)[train,]%>%select(-time)%>%na.omit, hidden=c(ncol(y.all2)-1,ncol(y.all2)%/%2), threshold=0.36, lifesign="full", act.fct="tanh")
  ##     sudriv$model$timedep$empir.model$nn <- nn
  ## }else{
  ##     nn <- sudriv$model$timedep$empir.model$nn
  ## }
  if(!is.null(cf)) pred.glm    <- predict.glm(glm2, newdata=as.data.frame(y.all2scaled)%>%select(-time), type="response")
  pred.linmod  <- predict.lm(linmod,  newdata=as.data.frame(y.all2scaled)%>%select(-time), type="response")
  pred.linmod2 <- predict.lm(linmod2, newdata=as.data.frame(y.all2scaled)%>%select(-time), type="response")
  pred.linmod3 <- predict.lm(linmod3, newdata=as.data.frame(y.all2scaled)%>%select(-time), type="response")
  
  ## convert time to date
  y.all2 <- y.all2 %>% mutate(time = as.POSIXct(sudriv$layout$tme.orig) + time * ifelse(sudriv$layout$time.units=="hours", 1, 24) * 60 *60)
  ## obs and pred data for first split
  compr <- data.frame(time=y.all2[,"time"], value=y.all2[,"y.td"], predobs="obs")
  #if(!is.null(cf)) compr <- rbind(compr, data.frame(time=y.all2[,"time"], value=pred.glm, predobs="pred.glm")) ## data frame for comparison of predictions
  compr <- rbind(compr, data.frame(time=y.all2[,"time"], value=pred.linmod, predobs="pred.linmod"))
  
  ## obs and pred data for second split
  compr2 <- data.frame(time=y.all2[,"time"], value=y.all2[,"y.td"], predobs="obs")
  #if(!is.null(cf)) compr2 <- rbind(compr2, data.frame(time=y.all2[,"time"], value=pred.glm, predobs="pred.glm")) ## data frame for comparison of predictions
  compr2 <- rbind(compr2, data.frame(time=y.all2[,"time"], value=pred.linmod2, predobs="pred.linmod"))
  
  ## obs and pred data for third split
  compr3 <- data.frame(time=y.all2[,"time"], value=y.all2[,"y.td"], predobs="obs")
  #if(!is.null(cf)) compr2 <- rbind(compr2, data.frame(time=y.all2[,"time"], value=pred.glm, predobs="pred.glm")) ## data frame for comparison of predictions
  compr3 <- rbind(compr3, data.frame(time=y.all2[,"time"], value=pred.linmod3, predobs="pred.linmod"))
  
  
  gg0 <- ggplot(y.all2, aes(x=time, y=y.td)) + geom_point() + labs(y=paste(tag,"observed"))
  gg1 <- ggplot(compr, aes(x=time, y=value, colour=predobs)) + geom_point(size=0.5) + labs(x="", y=tag, colour="") + theme(legend.position="top") + geom_vline(xintercept=y.all2[round(nrow(y.all2)*(1-validation_split)),"time"])
  gg2 <- ggplot(compr2, aes(x=time, y=value, colour=predobs)) + geom_point(size=0.5)  + labs(x="", y=tag, colour="") + theme(legend.position="none") + geom_vline(xintercept=y.all2[round(nrow(y.all2)*validation_split),"time"])
  gg3 <- ggplot(compr3, aes(x=time, y=value, colour=predobs)) + geom_point(size=0.5)  + labs(x="", y=tag, colour="") + theme(legend.position="none") + geom_vline(xintercept=y.all2[c(round(nrow(y.all2)*(0.5-0.5*validation_split)),round(nrow(y.all2)*(0.5+0.5*validation_split))),"time"])
  cw <- plot_grid(gg1,gg2,gg3,nrow=3)
  save_plot(filename=paste0("../output/timedeppar/A1Str07h2x/",tag,"/plot_predobs.pdf"), plot=cw, base_height=8)
  # ==============================================================================
  return(sudriv)
}
