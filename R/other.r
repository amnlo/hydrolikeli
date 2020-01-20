interp.precip <- function(time.p, precip, layout, avr.prev=TRUE){## ATTENTION: behaviour is unclear when precipitation contains missings.
    if(avr.prev){
        p.intp <- approx(x=time.p, y=precip, xout=layout$time)$y
    }else{
        t.tmp1 <- time.p - 0.01*abs(diff(time.p[1:2]))
        t.tmp2 <- time.p + 0.01*abs(diff(time.p[1:2]))
        times.p.double <- c(time.p[1]+diff(time.p[1:2]), c(t(cbind(t.tmp1,t.tmp2))))
        p.p.double     <- c(rep(precip, each=2), precip[length(precip)])
        p.intp <- approx(x=times.p.double, y=p.p.double, xout=layout$time)$y
    }
    return(p.intp)
}


every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) # credit: user "adamdsmith" on stackoverflow: https://stackoverflow.com/questions/34533472/insert-blanks-into-a-vector-for-e-g-minor-tick-labels-in-r/34533473#34533473
  {
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
      } else {
        x[1:nth != 1]
        }
    } else {
      if(empty) {
        x[1:nth != 1] <- ""
        x
        } else {
          x[1:nth == 1]
        }
    }
}
