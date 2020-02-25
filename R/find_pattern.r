find.pattern.timedep <- function(sudriv, vars=NULL, res=NULL, validation_split=0.2, add.data=NULL, tag="", Nd=1, keep.data=FALSE){
  ## This function compares the time course of the time dependent parameters to the model states, output (and potentially other variables) and identifies matching patterns.
  if(grepl("iniQE1",tag)){
    tag.ext <- "_iniQE1"
    iniQE.curr <- "QE1"
  }else{
    tag.ext <- ""
    iniQE.curr <- "QE3"
  }
  tag.red <- gsub("_.*","",tag)
  
  ## Prepare data
  y.all2 <- get.loess.input(sudriv=sudriv, tag=tag, vars=vars, res=res, add.data=add.data, remove.na=FALSE)

  ## Write data
  save(y.all2,file = paste0("../output/timedeppar/A1Str07h2x/",paste0(tag.red,tag.ext),".RData"))

  ## =================================================================================================
  ## scatterplot between timedep par and explanatory variables
  dir <- paste0("../output/timedeppar/A1Str07h2x/",tag,"/") 
  file.remove(paste0(dir,"r2_loess.txt"))
  conn <- file(paste0(dir,"r2_loess.txt"), "w")
  pdf(paste0(dir,"plot_scatter.pdf"))
  par(mfrow=c(2,2))
  loess.models <- mapply(function(x,y,nm,tag){
    dat <- data.frame(x=x,y=y) %>% arrange(x)
    smoothScatter(x=dat$x,y=dat$y,main=paste(tag,"&",nm),xlab=nm,ylab=tag.red,nrpoints=1000)
    sm <- loess(y~x,data=dat,span=1)
    pred <- predict(sm, newdata=dat)
    lines(x=dat$x, y=pred, col="red")
    r2 <- 1-sum((pred-dat$y)^2)/sum((dat$y-mean(dat$y))^2)
    write(c(tag,iniQE.curr,nm,r2), file=conn, ncolumns=4, append=TRUE)
    title(sub=bquote(R^2 == .(round(r2, 2))))
    return(sm)
  }, y.all2, nm=colnames(y.all2), MoreArgs=list(y=y.all2[,"y.td"], tag=tag.red), SIMPLIFY = FALSE)
  dev.off()
  sudriv$model$timedep$model.td <- loess.models
  sudriv$model$timedep$model.td.data <- y.all2
  if(Nd>1){ ## make N-dimensional loess analysis
    if(Nd>2) warning("Fitting all combinations of 3 or more variables might take a lot of time ...")
    data <- y.all2 %>% select(-time)
    data <- data[,c("y.td", grep("y.td", colnames(data), value=TRUE, invert=TRUE))]
    loess.Nd <- function(x, data){
      fm <- as.formula(paste0("y.td~",paste(colnames(data)[x], collapse="+")))
      ls <- tryCatch( loess(fm, data=data, span=1, control=loess.control(trace.hat="approximate")), error=function(x){message(x);return(NA)}, 
                      warning=function(x){message(x);return(NA)}) # trace.hat=approximate to be used for 1000 or more data points
      if(!is.na(ls[1])){
        pred <- predict(ls, newdata=data)
        r2 <- 1-sum((pred-data$y.td)^2)/sum((data$y.td-mean(data$y.td))^2)
      }else{
        r2 <- NA
      }
      nm <- paste(colnames(data)[x], collapse="+")
      print(nm)
      write(c(tag.red,iniQE.curr,nm,r2), file=conn, ncolumns=4, append=TRUE)
      if(!keep.data){ls$x <- NULL; ls$y <- NULL; ls$residual <- NULL} # save some space
      ls <- list(ls)
      names(ls) <- nm # make sure to keep the names of the feature combination (+)
      return(ls)
    }
    sm <- combn(2:ncol(data), Nd, loess.Nd, data=data, simplify=FALSE)
    sm <- unlist(sm, recursive=FALSE) # removes the uppermost layer of lists
    sudriv$model$timedep$model.td.Nd <- sm
  }
  close(conn)
  ## adapt the global table with loess results
  par.curr.dat <- read.table(paste0(dir,"r2_loess.txt"))
  names(par.curr.dat) <- c("parameter","iniQE","feature","r2")
  loess.dat <- read.table(paste0(dir,"../loess_overview_data.txt"), header=TRUE)
  loess.dat <- loess.dat %>% filter(!(parameter %in% par.curr.dat[,1] & iniQE==iniQE.curr)) %>% rbind(., par.curr.dat)
  write.table(loess.dat, paste0(dir,"../loess_overview_data.txt"), quote=FALSE)
  
  ## =================================================================================================
  ## scale the data for fitting glms later on
  if(FALSE){
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
  # cors <- apply(y.all2scaled, 2, cor, y=y.all2scaled[,"y.td"])
  # cat("correlations:\n")
  # print(cors)
  # ## test for stationarity
  # stati <- apply(y.all2scaled, 2, function(x) adf.test(x)$p.value)
  # names(stati) <- colnames(y.all2scaled)
  # cat("p-values under null hypothesis of non-stationarity:\n")
  # print(stati)
  # pdf(paste0("../output/timedeppar/A1Str07h2x/",tag,"/plot_crosscorr.pdf"))
  # par(mfrow=c(3,3))
  # mapply(function(y,x,lag.max,nm,plot,tag) ccf(x=y,y=x,lag.max=lag.max,plot=plot,main=paste(tag.red,"&",nm),ylim=c(-0.6,0.6)), y.all2scaled, nm=colnames(y.all2scaled), MoreArgs=list(y=y.all2scaled[,"y.td"], lag.max=2*7*24*4, plot=TRUE, tag=tag)) ## 2 weeks max lag
  # dev.off()
  
  # ================================================================================
  ## fit linear models
  if(FALSE){
    train <- (1:nrow(y.all2))[1:round(nrow(y.all2)*(1-validation_split))] ## train in beginning, test at end
    train2 <- (1:nrow(y.all2))[-test2]
    train3 <- (1:nrow(y.all2))[c(1:round(nrow(y.all2)*(0.5-0.5*validation_split)),round(nrow(y.all2)*(0.5+0.5*validation_split)):nrow(y.all2))] ## train at beginning and end, test in the middle

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
  }
  # ==============================================================================
  return(sudriv)
}
