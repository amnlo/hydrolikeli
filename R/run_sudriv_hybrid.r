run.sudriv.hybrid <- function(sudriv, ..., layout.model.td=NULL, lump=FALSE, td.ini=NULL, iter.max=50, tol=1e-3, verbose=0){
  ## This function runs the sudriv object with a timedependent parameter that is estimated based on a certain
  ## model that is required as input. The procedure is iterative: for a certain time-course of the parameter, the
  ## sudriv model is run and the states of the model are obtained over time. These states are then used as the input
  ## for the model that calculates a new time-course of the parameter, and so on until convergence is reached.
  
  if(is.null(td.ini)) td.ini <- sudriv$model$timedep$par
  if(is.null(layout.model.td)) layout.model.td <- sudriv$layout
  sudriv$model$timedep$par <- td.ini
  converged <- FALSE
  iter <- 0
  states.old <- 0
  while(!converged & iter < iter.max){
    
    states.new <- run.model(layout=layout.model.td, sudriv=sudriv, lump=lump)$original
    tmp <- model.td.wrapper(states.new, layout=layout.model.td$layout, ...)
    ## insert the new time series of the parameter into the sudriv object
    sudriv$model$timedep$par[(tmp$cut.beg+1):nrow(sudriv$model$timedep$par),] <- tmp$pred.td
    err <- sum(abs(states.new-states.old))
    states.old <- states.new
    
    converged <- err < tol
    iter <- iter + 1
    if(verbose>0){
      cat("Iteration:\t", iter, "\n")
      cat("Absolute error:\t", err, "\n")
    }
  }
  if(iter>=iter.max) warning("maximal number of iterations reached")
  # run final simulation with converged states and parameters
  y.mod <- run.model(layout=sudriv$layout, sudriv=sudriv, lump=lump) # note that the final output is for the layout of the sudriv object, not the layout we need for the internal states (layout.model.td)
  return(y.mod)
}

model.td.wrapper <- function(runmodel.out, layout, data, data.time, lags, col.names, model.td, args.model.td){
  ## this function calculates the time-course of the timedependent parameters based on the output of a run of the sudriv model
  inp.combined <- combine.states.inp(runmodel.out=runmodel.out, data=data, data.time=data.time, layout=layout)
  if(!is.null(col.names)) inp.combined <- inp.combined[,col.names] # sort columns according to col.names
  ## arrange data in 3D arrays by considering lags
  X <- 1:(nrow(inp.combined)-lags+1)
  tmp <- lapply(X = X, FUN = function(x,data) data[x:(x+lags-1),], data=inp.combined)
  tmp <- array(unlist(tmp), dim=c(lags,ncol(inp.combined),length(tmp)))
  dimnames(tmp)[[2]] <- colnames(inp.combined)
  inp.combined <- aperm(tmp, c(3,1,2))
  ## make sure that data is multiple of batch size, discard data at beginning
  cut.batch <- nrow(inp.combined)%%args.model.td$batch_size + 1
  inp.combined <- inp.combined[cut.batch:nrow(inp.combined),,]
  cat("predicting timdep par with lstm ...\n")
  pred.td <- do.call(model.td, c(list(X=inp.combined), args.model.td))
  cat("done.\n")
  ## keep track of the number of timesteps that were cut at the beginning of the series
  cut.beg <- lags-1+cut.batch-1
  return(list(pred.td=pred.td, cut.beg=cut.beg))
}

combine.states.inp <- function(runmodel.out, data, data.time, layout){
  ## this function combines the output of the sudriv model and the environmental input to produce the input for the model of the timedependent parameter
  ## 'data' is the scaled feature matrix used for prediction: it should contain the scaling info that was used in training
  runmodel.out <- data.frame(val=runmodel.out, time=layout$time, var=layout$var)
  runmodel.out <- as.data.frame(spread(runmodel.out, key=-val, value=val))
  if(nrow(runmodel.out)!=nrow(data)) stop(paste0("n.rows of model output: ",nrow(runmodel.out), " and feature data: ", nrow(data), " not the same"))
  if(!all(runmodel.out[,"time"] == data.time)) stop("time steps do not agree completely")
  runmodel.out <- runmodel.out %>% select(-time)
  if(!all(colnames(runmodel.out) %in% colnames(data))) stop("model output found for states that are not in features")
  ## scale data in same way as was done for training the LSTM
  model.out.scld <- myscale.boxcox(runmodel.out, box.param=box.param[,colnames(runmodel.out)])$dat.scaled
  model.out.scld2 <- scale(model.out.scld, center=attr(data, "scaled:center")[colnames(runmodel.out)], scale=attr(data, "scaled:scale")[colnames(runmodel.out)])
  ## insert the modelled states into the X data set
  srt <- match(colnames(model.out.scld2), colnames(data))
  data[,srt] <- model.out.scld2
  return(data)
}