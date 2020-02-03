run.sudriv.hybrid <- function(sudriv, model.td.wrapper, ..., layout.model.td=NULL, lump=FALSE, td.ini=NULL, iter.max=50, tol=1e-3, verbose=0){
  ## This function runs the sudriv object with a timedependent parameter that is estimated based on a certain
  ## model that is required as input. The procedure is iterative: for a certain time-course of the parameter, the
  ## sudriv model is run and the states of the model are obtained over time. These states are then used as the input
  ## for the model that calculates a new time-course of the parameter, and so on until convergence is reached.
  ## model.td.wrapper receives as a first argument the ouput of the run.model function
  
  if(is.null(td.ini)) td.ini <- sudriv$model$timedep$par
  if(is.null(layout.model.td)) layout.model.td <- sudriv$layout
  sudriv$model$timedep$par <- td.ini
  converged <- FALSE
  iter <- 0
  states.old <- 0
  while(!converged & iter < iter.max){
    
    states.new <- run.model(layout=layout.model.td, sudriv=sudriv, lump=lump)$original
    td.new <- model.td.wrapper(states.new, layout=layout.model.td$layout, ...)
    sudriv$model$timedep$par <- td.new
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

model.td.wrapper <- function(runmodel.out, inp.env, model.td, param, layout, col.names){
  ## this function calculates the time-course of the timedependent parameters based on the output of a run of the sudriv model
  inp.combined <- combine.states.inp(runmodel.out=runmodel.out, inp.env=inp.env, layout=layout)
  if(!is.null(col.names)) inp.combined <- inp.combined[,col.names] # sort columns according to col.names
  pred.td <- model.td(inp.combined, param)
  return(pred.td)
}

combine.states.inp <- function(runmodel.out, inp.env, layout){
  ## this function combines the output of the sudriv model and the environmental input to produce the input for the model of the timedependent parameter
  model.out.data <- runmodel.out
  model.out.data <- data.frame(val=model.out.data, time=layout$time, var=layout$var)
  model.out.data <- as.data.frame(spread(model.out.data, key=-val, value=val))
  data <- cbind(inp.env, model.out.data)
  return(data)
}

model.td <- function(x, param){
  return(matrix(-6, nrow=21601))
}
