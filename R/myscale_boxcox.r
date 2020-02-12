myscale.boxcox <- function(data, box.param=NULL){
  if(is.null(box.param)){  
    box.param <- apply(data, 2, function(x){
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
  }
  datascaled <- data
  for(i in 1:ncol(data)){
    if(any(is.na(box.param[,i]))) next
    if(box.param[1,i]!=0){
      datascaled[,i] <- (((datascaled[,i]+box.param[2,i])^box.param[1,i] - 1)/box.param[1,i] - box.param[3,i])/sqrt(box.param[4,i])
    }else{
      datascaled[,i] <- (log(datascaled[,i]+box.param[2,i]) - box.param[3,i])/sqrt(box.param[4,i])
    }
  }
  return(list(dat.scaled=datascaled, box.param=box.param))
}