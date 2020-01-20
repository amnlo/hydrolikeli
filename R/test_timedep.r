test.timedeppar <- function(sudriv){
  ## run with constant parameters
  cat("running no timedep pars...\n")
  res.dummy <- run.engine(sudriv)
  const <- res.dummy$y
  cat("done\n")
  
  layout.ur <- list(layout = data.frame(var=rep("U3F1Wv_Su1", nrow(sudriv$input$inputobs)), time=sudriv$input$inputobs[,1], stringsAsFactors=FALSE),
                    lump   = rep(NA, nrow(sudriv$input$inputobs)))
  layout.ur.sr <- list(layout = data.frame(var=rep(c("U3F1Wv_Su1","U3F1Wv_Ss1"), each=nrow(sudriv$input$inputobs)), time=rep(sudriv$input$inputobs[,1],2), stringsAsFactors=FALSE),
                       lump   = rep(NA, 2*nrow(sudriv$input$inputobs)))
  layout.q <- list(layout = data.frame(var=rep("C1Wv_Qstream", nrow(sudriv$input$inputobs)), time=sudriv$input$inputobs[,1], stringsAsFactors=FALSE),
                   lump   = rep(NA, nrow(sudriv$input$inputobs)))
  ## run with time-dependent parameter, which is constant
  ## td <- pmin(pmax(randou(mean=0.95,sd=0.01,tau=10000,t=1:nrow(sudriv$input$inputobs))$y, 0.51), 0.999)
  maxUR.orig <- as.numeric(sudriv$model$parameters["Glo%CmltSmax_UR"])
  td <- rep(maxUR.orig, nrow(sudriv$input$inputobs))
  sudriv$model$timedep$par <- matrix(td,ncol=1)
  sudriv$model$timedep$pTimedep <- names(sudriv$model$parameters) %in% c("Glo%CmltSmax_UR") # indicate which parameter is time-dependent
  cat("running with constant timedep pars...\n")
  timedep.const <- run.engine(sudriv)$y
  res.dummy1 <- run.model(layout=layout.ur, sudriv=sudriv, lump=FALSE)$original
  td <- rep(3, nrow(sudriv$input$inputobs))
  sudriv$model$timedep$par <- matrix(td,ncol=1)
  res.dummy2 <- run.model(layout=layout.ur, sudriv=sudriv, lump=FALSE)$original
  td <- rep(maxUR.orig, nrow(sudriv$input$inputobs))
  sudriv$model$timedep$par <- matrix(td,ncol=1)
  res.dummy3 <- run.model(layout=layout.ur, sudriv=sudriv, lump=FALSE)$original
  cat("done\n")
  
  sudriv$input$inputobs[1,"P"] <- 10 ## to test the first value of the timedep parameter
  ## running with different values for Pmax_ED
  td <- rep(log(15), nrow(sudriv$input$inputobs)) # this should be above max precip (12)
  sudriv$model$timedep$par <- matrix(td,ncol=1)
  sudriv$model$timedep$pTimedep <- names(sudriv$model$parameters) %in% c("Glo%Cmlt_Pmax_ED") # indicate which parameter is time-dependent
  print(names(sudriv$model$parameters)[sudriv$model$timedep$pTimedep])
  cat("running with constant timedep pars...\n")
  res.maxED15 <- run.model(layout=layout.q, sudriv=sudriv, lump=FALSE)$original
  td <- rep(log(19), nrow(sudriv$input$inputobs)) # this should be above max precip (12)
  sudriv$model$timedep$par <- matrix(td,ncol=1)
  res.maxED19 <- run.model(layout=layout.q, sudriv=sudriv, lump=FALSE)$original
  td <- rep(log(11), nrow(sudriv$input$inputobs))
  sudriv$model$timedep$par <- matrix(td,ncol=1)
  res.maxED11 <- run.model(layout=layout.q, sudriv=sudriv, lump=FALSE)$original
  lrgthan11 <- which(sudriv$input$inputobs[,"P"]>=11)[1]
  td <- pmin(pmax(rnorm(nrow(sudriv$input$inputobs), log(15), 0.05),log(13)),log(19.9)) # this should be above max precip (12)
  sudriv$model$timedep$par <- matrix(td,ncol=1)
  res.maxEDvar <- run.model(layout=layout.q, sudriv=sudriv, lump=FALSE)$original
  ##td.keep <- ifelse(rollmean(sudriv$input$inputobs[,"P"],k=2,fill=c(0,NULL,0),align="right")>1e-9,log(15),log(1))
  ##td.keep <- ifelse(c(10,sudriv$input$inputobs[1:(nrow(sudriv$input$inputobs)-1),"P"])>0,log(15),log(1)) # this should always be above precipitation
  cat("running with first value different ...\n")
  td.keep <- c(log(5), rep(log(15), nrow(sudriv$input$inputobs)-1))
  sudriv$model$timedep$par <- matrix(td.keep,ncol=1)
  res.maxEDevade <- run.model(layout=layout.q, sudriv=sudriv, lump=FALSE)$original
  cat("done\n")
  
  ## run with simple step function as time-dependent parameter
  td <- rep(maxUR.orig, nrow(sudriv$input$inputobs))
  td[2000:length(td)] <- td[2000:length(td)] - log(2) # adapt second half of time series of parameter
  sudriv$model$timedep$par <- matrix(td,ncol=1)
  sudriv$model$timedep$pTimedep <- names(sudriv$model$parameters) %in% c("Glo%CmltSmax_UR") # indicate which parameter is time-dependent
  cat("running with step function timedep par...\n")
  step <- run.engine(sudriv)$y
  res.dummy.step <- run.model(layout=layout.ur.sr, sudriv=sudriv, lump=FALSE)$original
  pdf("stepfunction_test.pdf")
  plot(res.dummy.step[layout.ur.sr$layout$var=="U3F1Wv_Su1"])
  points(res.dummy.step[layout.ur.sr$layout$var=="U3F1Wv_Ss1"], col="red")
  dev.off()
  cat("done\n")
  
  ## run with chaotic timedep par
  set.seed(11)
  td <- rnorm(nrow(sudriv$input$inputobs), mean=log(100), sd=1)
  td <- pmin(pmax(td, log(30)), log(300))
  cat("chaos Smax_UR:\n")
  print(summary(td))
  sudriv$model$timedep$par <- matrix(td,ncol=1)
  sudriv$model$timedep$pTimedep <- names(sudriv$model$parameters) %in% c("Glo%CmltSmax_UR") # indicate which parameter is time-dependent
  cat("running with chaotic parameters ...\n")
  chaos <- run.engine(sudriv)$y
  print(length(chaos))
  cat("done\n")
  
  
  cat("test sucessful if all following are true...\n")
  print(length(const)==length(timedep.const))
  dff <- const - timedep.const
  print(all(dff==0))
  print(all(res.dummy1==res.dummy3))
  print(!all(res.dummy1==res.dummy2))
  print(all(res.maxED15==res.maxED19))
  print(!all(res.maxED15==res.maxED11))
  print(all(res.maxED15[1:(lrgthan11-1)]==res.maxED11[1:(lrgthan11-1)]))
  print(all(res.maxED15==res.maxEDvar))
  print("evade:")
  print(!all(res.maxED15==res.maxEDevade))
  ## cat(sum(res.maxED15!=res.maxEDevade)," / ",length(res.maxEDevade),"\n")
  ## dif <- res.maxED15-res.maxEDevade
  ## print(summary(dif))
  ## print(summary(which(dif!=0)))
  ## cat("max dif at: ",which.max(abs(dif)),"\n")
  ## print(sort(dif)[1:50])
  ## cat("max precip at: ",which.max(sudriv$input$inputobs[,"P"]),"\n")
  ## cat("max precip: ",max(sudriv$input$inputobs[,"P"]),"\n")
  ## cat("precip larger0: ",sum(sudriv$input$inputobs[,"P"]>0),"\n")
  ## print(cbind(7850:7880,sudriv$input$inputobs[7850:7880,"P"]))
  ## print(cbind(7850:7880,td.keep[7850:7880]))
  ## print(cbind(7850:7880,dif[7850:7880]))
  
  ## test is first half of timeseries is equal and second is not (step function)
  print(head(sudriv$model$outnames))
  qs <- result_index_var(res.dummy, file.o=NA, variables=c("C1Wv_Qstream"), outnames=sudriv$model$outnames)
  print(qs)
  q.const <- const[(qs[[1]][1]):(qs[[1]][2])]
  q.step  <- step[(qs[[1]][1]):(qs[[1]][2])]
  q.chaos <- chaos[(qs[[1]][1]):(qs[[1]][2])]
  cat("total streamflow for chaos: ", round(sum(q.chaos),5),"\n")
  print(summary((q.const - q.step)[1:1999]))
  print(summary((q.const - q.step)[2000:length(td)]))
  print(summary(q.const - q.chaos))
}
