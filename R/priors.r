param.logprior <- function(const.par){
  ## this is the hard-coded joint prior for the constant parameters in the 'timedeppar' package
  distdef <- list(`Glo%Cmlt_P` = c("normaltrunc","1","0.05","0.5","2"),
                  `Glo%Cmlt_E` = c("normaltrunc", "0", "0.2","-1","1"),
                  `Glo%Cmlt_Dspl_SD`=c("lognormaltrunc","0.7","0.3","0.5","1"),
                  `Glo%Cmlt_Pmax_ED`= c("normaltrunc","2.3","0.2","0.7","2.99"),
                  `Glo%CmltSmax_IR`= c("normaltrunc","2.7","0.1","2.302585","3.912023"),
                  `Glo%CmltSmax_UR`= c("normaltrunc","5.5","0.5","3.4","6.9"),
                  `Glo%Cmlt_BeQq_UR`= c("normaltrunc","0.9","1","-2.3","1.8"),
                  `Glo%Cmlt_K_Qq_RR`= c("normaltrunc","-0.5","0.5","-2.5","1"),
                  `Glo%Cmlt_K_Qq_FR`= c("normaltrunc","-2.2","0.3","-3.5","-1"),
                  `Glo%Cmlt_AlQq_FR`= c("normaltrunc","0.2","0.2","0","1.2"),
                  `Glo%Cmlt_K_Qq_SR`= c("normaltrunc","-6","1","-12","-5"),
                  `Glo%Cmlt_AlQq_SR`= c("normaltrunc","0.7","0.3","0","1.8"),
                  `GloTr%CmltSlOne_IR`= c("normaltrunc","3","0.5","0","7"),
                  `GloTr%CmltKd_WR`= c("normaltrunc","-6.8","0.4","-9.21034","-1"),
                  `GloTr%CmltSlTwo_IR`= c("normaltrunc","3.8","0.5","0","7"),
                  `GloTr%CmltRs_WR`= c("normaltrunc","-5","0.5","-9.21034","-1"),
                  `U1W%KpQq_FR`  = c("normaltrunc","-2","0.5","-4","1"),
                  `GLOB_Mult_Q_taumax_lik`=c("normaltrunc","4","1","0","6"),
                  `C1Wv_Qstream_a_lik`=c("exponential","1"),
                  `C1Wv_Qstream_b_lik`=c("normaltrunc","-0.5","0.5","-4","0.5","NA"),
                  `C1Tc1_Qstream_a_lik`=c("exponential","1"),
                  `C1Tc2_Qstream_a_lik`=c("exponential","1"))
  fmean <- names(const.par)[grep("_fmean", names(const.par))]
  if(length(fmean)>0){
    fmean <- gsub("_fmean","",fmean)
    ind.fmean <- names(distdef) %in% fmean
    names(distdef)[ind.fmean] <- paste0(names(distdef)[ind.fmean],"_fmean")
    if("Glo%Cmlt_Dspl_SD" %in% fmean){ ## consider that this parameter is sigmoid transformed
      distdef[["Glo%Cmlt_Dspl_SD_fmean"]] <- c("normal", "-0.4987", "1.6741")
    }else if("Glo%tStart_VL" %in% fmean){ ## consider that this parameter is sigmoid transformed
      distdef[["Glo%tStart_VL_fmean"]] <- c("normal", "0.3252", "0.8542")
    }else if("Glo%Cmlt_P" %in% fmean){ ## consider that this parameter is sigmoid transformed
      distdef[["Glo%Cmlt_P_fmean"]] <- c("normal", "-1.386294", "0.1")
    }else if("Glo%CmltSmax_IR" %in% fmean){ ## consider that the lower bound of the prior was changed here
      distdef[["Glo%CmltSmax_IR_fmean"]] <- c("normaltrunc","2.7","0.1","0","3.912023")
    }
  }
  distdef  <- distdef[names(distdef) %in% names(const.par)]
  distdef  <- distdef[match(names(const.par), names(distdef))]
  ## =======================================================
  srp <- grepl("CmltSl.{3}_IR",names(const.par))
  cnst.test <- names(const.par)[!srp]
  if(!all(cnst.test %in% names(distdef))) stop("prior for some fitted parameters is not implemented")
  ## calculate logprior
  mvprior <- TRUE
  pri <- list(dist="indep", mean=0, sd=1, cor=0, cor.inv=NA, log=TRUE, distdef=distdef)
  if(mvprior){
    pri$distdef <- pri$distdef[!srp]
  }else{
    pri$distdef <- pri$distdef
  }
  if(mvprior){
    args.pdf       <- c(list(z=const.par[!srp]), pri)
    sorp.pars            <- const.par[srp]
    if((any(sorp.pars <= 0)) | any(sorp.pars >= 7)){
      logpri.sorp <- -Inf
    }else{
      logpri.sorp          <- log(dlnorm.rplus(exp(sorp.pars), meanlog=c(3,3.8), varlog=matrix(c(0.52,0.48,0.48,0.52), ncol=2)))
    }
  }else{
    args.pdf       <- c(list(z=const.par), pri)
    logpri.sorp          <- 0
  }
  logprior      <- do.call(calcpdf_mv, args.pdf) + logpri.sorp
  return(logprior)
}
param.ou.logprior <- function(oupar){
  distdef.mn <- list(`Glo%Cmlt_P` = c("normaltrunc","1","0.05","0.5","2"),
                     `Glo%Cmlt_E` = c("normaltrunc", "0", "0.2","-1","1"),
                     `Glo%Cmlt_Dspl_SD`=c("lognormaltrunc","0.7","0.3","0.5","1"),
                     `Glo%tStart_VL`=c("lognormaltrunc","1.2","0.4","0","2"),
                     `Glo%Cmlt_Pmax_ED`= c("normaltrunc","2.3","0.2","0.7","2.99"),
                     `Glo%CmltSmax_IR`= c("normaltrunc","2.7","0.1","0","3.912023"), # note that the lower bound was relaxed here
                     `Glo%CmltSmax_UR`= c("normaltrunc","5.5","0.5","3.4","6.9"),
                     `Glo%Cmlt_BeQq_UR`= c("normaltrunc","0.9","1","-2.3","1.8"),
                     `Glo%Cmlt_K_Qq_RR`= c("normaltrunc","-0.5","0.5","-2.5","1"),
                     `Glo%Cmlt_K_Qq_FR`= c("normaltrunc","-2.2","0.3","-3.5","-1"),
                     `Glo%Cmlt_AlQq_FR`= c("normaltrunc","0.2","0.2","0","1.2"),
                     `Glo%Cmlt_K_Qq_SR`= c("normaltrunc","-6","1","-12","-5"),
                     `Glo%Cmlt_AlQq_SR`= c("normaltrunc","0.7","0.3","0","1.8"),
                     `GloTr%CmltSlOne_IR`= c("normaltrunc","3","0.5","0","7"),
                     `GloTr%CmltKd_WR`= c("normaltrunc","-6.8","0.4","-9.21034","-1"),
                     `GloTr%CmltSlTwo_IR`= c("normaltrunc","3.8","0.5","0","7"),
                     `GloTr%CmltRs_WR`= c("normaltrunc","-5","0.5","-9.21034","-1"),
                     `U1W%KpQq_FR`  = c("normaltrunc","-2","0.5","-4","1"),
                     `GLOB_Mult_Q_taumax_lik`=c("normaltrunc","4","1","0","6"),
                     `C1Wv_Qstream_a_lik`=c("exponential","1"),
                     `C1Wv_Qstream_b_lik`=c("normaltrunc","-0.5","0.5","-4","0.5","NA"),
                     `C1Tc1_Qstream_a_lik`=c("exponential","1"),
                     `C1Tc2_Qstream_a_lik`=c("exponential","1"))
  ## is the parameter transformed?
  tran <- c(`Glo%Cmlt_P` = FALSE,
            `Glo%Cmlt_E` = TRUE,
            `Glo%Cmlt_Dspl_SD`= FALSE,
            `Glo%tStart_VL`= FALSE,
            `Glo%Cmlt_Pmax_ED`= TRUE,
            `Glo%CmltSmax_IR`= TRUE, # note that the lower bound was relaxed here
            `Glo%CmltSmax_UR`= TRUE,
            `Glo%Cmlt_BeQq_UR`= TRUE,
            `Glo%Cmlt_K_Qq_RR`= TRUE,
            `Glo%Cmlt_K_Qq_FR`= TRUE,
            `Glo%Cmlt_AlQq_FR`= TRUE,
            `Glo%Cmlt_K_Qq_SR`= TRUE,
            `Glo%Cmlt_AlQq_SR`= TRUE,
            `GloTr%CmltSlOne_IR`= TRUE,
            `GloTr%CmltKd_WR`= TRUE,
            `GloTr%CmltSlTwo_IR`= TRUE,
            `GloTr%CmltRs_WR`= TRUE,
            `U1W%KpQq_FR`  = TRUE,
            `GLOB_Mult_Q_taumax_lik` = TRUE,
            `C1Wv_Qstream_a_lik`= TRUE,
            `C1Tc1_Qstream_a_lik`= TRUE,
            `C1Tc2_Qstream_a_lik`= TRUE)
  pp.all <- gsub("_mean","",names(oupar))
  pp.all <- gsub("_sd","",pp.all)
  pp.all <- gsub("_gamma","",pp.all)
  distdef <- list()
  for(pp in pp.all){
    sd.fac <- NA
    if(pp=="Glo%Cmlt_BeQq_UR"){
      sd.fac <- 0.05 # for this parameter, we do not expect large changes in time
    }else if(pp=="Glo%Cmlt_P"){
      sd.fac <- 1 # for this parameter, we expect that a large part of the possible change varies also in time
    }else{
      sd.fac <- 0.2
    }
    distdef.sd.ou <- list(c("normaltrunc", "0", as.character(as.numeric(distdef.mn[[pp]][3])*sd.fac), "0", "10"))
    names(distdef.sd.ou) <- paste0(pp,"_sd")
    names(distdef.mn) <- paste0(names(distdef.mn), "_mean") ## use the same prior for the mean of the OU process as for the constant parameters, if reparameterization is happening, the mean is not fitted and this not used
    distdef.mn.ou <- distdef.mn[names(distdef.mn) %in% names(oupar)]
    if(pp=="Glo%Cmlt_Dspl_SD"){ ## consider that this parameter is sigmoid transformed
      distdef.mn.ou <- list(c("normal", "-0.4987", "1.6741"))
      names(distdef.mn.ou) <- "Glo%Cmlt_Dspl_SD_mean"
      distdef.sd.ou <- list(c("normaltrunc", "0", as.character(as.numeric(distdef.mn.ou[[1]][3])*sd.fac), "0", "10"))
      names(distdef.sd.ou) <- "Glo%Cmlt_Dspl_SD_sd"
    }else if(pp=="Glo%tStart_VL"){ ## consider that this parameter is sigmoid transformed
      distdef.mn.ou <- list(c("normal", "0.3252", "0.8542"))
      names(distdef.mn.ou) <- "Glo%tStart_VL_mean"
      distdef.sd.ou <- list(c("normaltrunc", "0", as.character(as.numeric(distdef.mn.ou[[1]][3])*sd.fac), "0", "10"))
      names(distdef.sd.ou) <- "Glo%tStart_VL_sd"
    }else if(pp=="Glo%Cmlt_P"){ ## consider that this parameter is sigmoid transformed
      distdef.mn.ou <- list(c("normal", "-1.386294", "0.1"))
      names(distdef.mn.ou) <- "Glo%Cmlt_P_mean"
      distdef.sd.ou <- list(c("normaltrunc", "0", as.character(as.numeric(distdef.mn.ou[[1]][3])*sd.fac), "0", "10"))
      names(distdef.sd.ou) <- "Glo%Cmlt_P_sd"
    }
    distdef  <- c(distdef, distdef.mn.ou, distdef.sd.ou)
  }
  distdef <- distdef[names(oupar)] # use only the ou parameters that are fitted
  distdef  <- distdef[match(names(oupar), names(distdef))]
  pri <- list(dist="indep", mean=0, sd=1, cor=0, cor.inv=NA, log=TRUE, distdef=distdef)
  args.pdf       <- c(list(z=oupar), pri)
  logprior <- do.call(calcpdf_mv, args.pdf)
  return(logprior)
}
