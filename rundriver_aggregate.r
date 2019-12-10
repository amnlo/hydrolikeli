op.sys <- Sys.info()["sysname"]
if(op.sys == "Linux"){.libPaths("/local/ammannlo/personalR/library")}else{Sys.setlocale("LC_TIME", "English")}
options(repos="https://stat.ethz.ch/CRAN/")
if ( ! require(ecosim) )     { install.packages("ecosim");     library(ecosim) }
if ( ! require(pryr) )     { install.packages("pryr");     library(pryr) }
if ( ! require(compare) )     { install.packages("compare");     library(compare) }
if ( ! require(zoo) )     { install.packages("zoo");     library(zoo) }
if ( ! require(reshape) )     { install.packages("reshape");     library(reshape) }
if ( ! require(tidyr) )     { install.packages("tidyr");     library(tidyr) }
if ( ! require(dplyr) )     { install.packages("dplyr");     library(dplyr) }
if ( ! require(magrittr) )     { install.packages("magrittr");     library(mgrittr) }
if ( ! require(stringr) )     { install.packages("stringr");     library(stringr) }
if ( ! require(nloptr) )     { install.packages("nloptr");     library(nloptr) }
if ( ! require(truncnorm) )     { install.packages("truncnorm");     library(truncnorm) }
if ( ! require(compositions) )     { install.packages("compositions");     library(compositions) }
if ( ! require(rjson) )     { install.packages("rjson");     library(rjson) }
if ( ! require(ggplot2) )     { install.packages("ggplot2");     library(ggplot2) }
if ( ! require(ggridges) )     { install.packages("ggridges");     library(ggridges) }
if ( ! require(GGally) )     { install.packages("GGally");     library(GGally) }
if ( ! require(grid) )     { install.packages("grid");     library(grid) }
if ( ! require(gridExtra) )     { install.packages("gridExtra");     library(gridExtra) }
if ( ! require(cowplot) )     { install.packages("cowplot");     library(cowplot) }
if ( ! require(egg) )     { install.packages("egg");     library(egg) }
if ( ! require(ggpubr) )     { install.packages("ggpubr");     library(ggpubr) }
if ( ! require(skewt) )     { install.packages("skewt");     library(skewt) }
if ( ! require(scoringRules) )     { install.packages("scoringRules");     library(scoringRules) }
if ( ! require(parallel) )     { install.packages("parallel");     library(parallel) }
if ( ! require(doParallel) )     { install.packages("doParallel");     library(doParallel) }
if ( ! require(dream) )     { install.packages("dream", repos="http://R-Forge.R-project.org");     library(dream) }
#if ( ! require(devtools) )     { install.packages("devtools");     library(devtools) }
# --------------------------------
# Change this:
reload 				<- TRUE
reload.timedep		<- TRUE
sample.posterior 	<- FALSE
optimize.posterior 	<- FALSE
nparallel 			<- 1
plot.chains 		<- FALSE
pred 				<- FALSE
calc.met            <- FALSE
summarize.results 	<- FALSE
plot.states 		<- TRUE
set.par             <- NA#"../../A1Str07h2x/v05_connred3aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExc0/opt_par_v05_connred3aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExc0.txt"
plot.uncert		 	<- FALSE
sensitivity.analysis <- FALSE
place.terb.early    <- TRUE ## ATTENTION: is this what you want?
analyze.likelihood <- FALSE
pritag             <- "Exc0"
errmod.conc		   <- "E1"
struct             <- "Str07h2xTimedep"
modNam			   <- "MexpH4" # only used for plot.markov.hist
tag 			   <- paste0("v05_connred3aggr_QE3P7c07i_TT",errmod.conc,"b02c08i_SE1c08_prior",pritag)
if(reload.timedep) tag	   <- "kqqsr2_QE1_adptintrv"
sett  			   <- paste0("settings_",struct,"_A1_QE3c07i_T",errmod.conc,"b02c08i_SE1c08.json")
start.par          <- NA#paste0("../../A1","Str07h2x","/v05_connred3aggr_QE1P7c07i_TT",errmod.conc,"b02c08i_SE1c08_prior","Exc4","/opt_par_v05_connred3aggr_QE1P7c07i_TT",errmod.conc,"b02c08i_SE1c08_prior","Exc4",".txt")
vary			   <- 1.1 # vary the initial parameters by this factor to create the starting interval
time.units <- "hours"
maxtime <- 8*60*60 # h * min * sec for optimization
timestep.fac <- 0.25
brn.in <- 250
thin <- 10
max.iter  <- 100000
iter.step <- 100000
n.walkers <- 50
jitter <- 0
n.sample <- 500
cross.valid  <- FALSE
vars.valid <- c("C3Wv_Qstream", "C3Tc1_Qstream", "C3Tc2_Qstream")
vars.calib <- c("C4Wv_Qstream", "C4Tc1_Qstream", "C4Tc2_Qstream", "HYPERSTATE_sorption")
compare.models <- NA#c(#MexpH1="settings_Str070ts_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#MexpH2="settings_Str074ts_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#MexpH3a="settings_Str0725ts_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#MexpH3b="settings_Str0726ts_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#MprxH4="settings_Str076ts_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#MexpH4_ns="settings_Str07h2x_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#MexpH4Corr="settings_Str07h2x_A1_QE3c07i_TE2b02c08i_SE1c08.json",
					#MexpH4="settings_Str07h2x_A1_QE3c07i_TE1b02c08i_SE1c08.json")
					#Mkqqfr="settings_Str07h2xTimedep_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#Mkqqrr="settings_Str07h2xTimedep_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#Mkpqqfr="settings_Str07h2xTimedep_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#Mkdwr="settings_Str07h2xTimedep_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#Mrswr="settings_Str07h2xTimedep_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#Mdsplsd="settings_Str07h2xTimedep_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#Malqqfr="settings_Str07h2xTimedep_A1_QE3c07i_TE1b02c08i_SE1c08.json",
					#Mnone="settings_Str07h2xTimedep_A1_QE3c07i_TE1b02c08i_SE1c08.json")
					#Msmaxur="settings_Str07h2xTimedep_A1_QE3c07i_TE1b02c08i_SE1c08.json")
tag.mult <- c(#MexpH1="v05_connred3aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExb0",
			  #MexpH2="v05_connred3aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExc0",
			  #MexpH3a="v05_connred2aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExb",
			  #MexpH3b="v05_connred2aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExb",
			  #MprxH4="v05_proxim2aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExc0",
			  #MexpH4_ns="v05_nosprconnred2aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExb",
			  #MexpH4Corr="v05_connred3aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExb0",
			  #MexpH4="v05_connred3aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExc0")
			  #Mkqqfr="kqqfr_QE3",
			  Mkqqrr="kqqrr_QE3",
			  #Mkpqqfr="kpqqfr_QE3",
			  Mkdwr="kdwr_QE3",
			  Mrswr="rswr_QE3",
			  #Mdsplsd="dsplsd_QE3",
			  #Malqqfr="alqqfr_QE3",
			  Mnone="none_QE3")
			  #Msmaxur="smaxur_QE3")
export <- FALSE
arrange=rep(1,length(compare.models))
names(arrange) <- names(compare.models)
plot.var=c("C1Wv_Qstream", "C1Tc1_Qstream", "C1Tc2_Qstream")
if(grepl("dsplsd",tag)){
	scaleshift = matrix(c(0.5,0.5), nrow=1)
	rownames(scaleshift) <- "Glo%Cmlt_Dspl_SD_fmean"
}else if(grepl("tstartvl",tag)){
	scaleshift = matrix(c(2,0), nrow=1)
	rownames(scaleshift) <- "Glo%tStart_VL_fmean"
}else if(grepl("cmltp",tag)){
	scaleshift = matrix(c(1.4-0.9,0.9), nrow=1)
	rownames(scaleshift) <- "Glo%Cmlt_P_fmean"
}else{
	scaleshift=NA
}
# ------------------------------
tme.orig <- "2008-01-01 00:00:00"
reproduce.synthetic <- FALSE
dir.driver <- "../../../../SuperFlex/Driver/sudriv_package/R"
source("source.functions.r")
source.functions(dir.driver)
source("../../../../Likelihood/src/hydrolikeli_herb.r")
source("../../../../Likelihood/src/plot.markov.r")
options(stringsAsFactors=FALSE)
su <- sudriv.setup(settings=sett)
su$settings$file.par <- paste("p",su$settings$structure,su$settings$par.tag,".par",sep="")

if(reload.timedep){
	savepath <- paste("../output/timedeppar/", su$settings$subcatchment, su$settings$structure, "/", tag, sep="")
}else{
	savepath <- paste("../output/", su$settings$subcatchment, su$settings$structure, "/", tag, sep="")
}
dir.create(savepath, recursive=TRUE)
if(op.sys == "Linux"){
    pth.base <- "\"/local/ammannlo/PhD/Herbicides/Workprog/Modelling/"
    f.path.config <- paste(pth.base, "config/flexConfigEschMain_", su$settings$subcatchment, "_", su$settings$structure, "_", su$settings$tracer, ".dat", "\"", sep = "")
    f.path.hru 	    <- paste("/local/ammannlo/PhD/Herbicides/Workprog/Modelling/config/", su$settings$subcatchment, "_", su$settings$input.tag, ".dat", sep="")
    f.path.hyperfun <- paste("/local/ammannlo/PhD/Herbicides/Workprog/Modelling/config/hyperfun_", su$settings$subcatchment, "_", su$settings$structure,".RData", sep="")
}else{
    pth.base <- "Q:\\Abteilungsprojekte\\siam\\ammannlo\\PhD\\Herbicides\\Workprog\\Modelling\\"
    f.path.config <- paste(pth.base, "config\\flexConfigEschMain_", su$settings$subcatchment, "_", su$settings$structure, "_", su$settings$tracer, ".dat", sep = "")
    f.path.hru <- paste(pth.base, "config\\", su$settings$subcatchment, "_", su$settings$input.tag, ".dat", sep="")
    f.path.hyperfun <- paste(pth.base, "config\\hyperfun_", su$settings$subcatchment, "_", su$settings$structure,".RData", sep="")
}
file.config <- file("path_config.txt")
write(f.path.config, file.config)
close(file.config)
su <- model.setup(su, settings=sett, writeO=TRUE, f.path.hru=f.path.hru, f.path.transK=NA, f.path.hyperfun=ifelse(grepl("SE[0-9]",tag) | reload.timedep, f.path.hyperfun, NA)) ## ATTENTION: HRU file is read here! It has to agree with the object that is loaded. (f.path.config).
if(!reload){
	## Add hyperstates layout and observations for sorption
	if(grepl("SE[0-9]", tag)){
		hyp.lay <- data.frame(var="HYPERSTATE_sorption", time=c(1,2,3,4,5,6))
		su$layout$layout <- rbind(su$layout$layout, hyp.lay)
		su$observations <- c(su$observations, rep(2.1,6))##c(1.5,1.5,1.5,2.1,4.0,13)) #sorption coefficient calculated by Camenzuli (l/kg)
	}
	su <- parameters.setup(su, with.names=TRUE, settings=su$settings, replace.param=TRUE)
	su <- optimizer.setup(su, settings=sett)
	su <- likelihood.setup(su, settings=sett, replace.param=TRUE)
	su <- construct.prior(su, file=paste0("../sudriv_input/priordef_S07",pritag,".txt"))
	su$layout$calib <- 1:nrow(su$layout$layout)
	su$input$inputobs[,1] <- su$input$inputobs[,1]*24
	su$layout$layout[,"time"] <- su$layout$layout[,"time"]*24
	if(grepl("SE[0-9]",tag)) su$layout$layout[su$layout$layout$var=="HYPERSTATE_sorption","time"] <- 12108+c(0,3,7,15,30,60)*24
	su <- regularize.layout(su, threshold=40000, set.dt=0.25, remove=NA, aggregate=TRUE)
    rem <- numeric()
    if(grepl("A3",sett)) rem <- which((su$layout$layout$var=="C1Wv_Qstream" & ((su$layout$layout$time >= 11095 & su$layout$layout$time <= 11101.01) | (su$layout$layout$time >= 14940 & su$layout$layout$time <= 14945))) | (su$layout$layout$var=="C1Tc1_Qstream" & (su$layout$layout$time >= 12607.36 & su$layout$layout$time <= 12607.54))) # remove some strange data points in A3
	rem <- c(rem,which(grepl("Tc2_Qstream", su$layout$layout$var) & (su$layout$layout$time/24)>515)) # exclude observations of terbuthylazine after e2
	if(length(rem)==0) rem <- NA
	print("rem:")
	print(length(rem))
	print(summary(rem))
	su <- regularize.layout(su, threshold=Inf, remove=rem)
	print("finished regularizing layout")
	if(cross.valid){
		su$layout$pred  <- su$layout$calib[su$layout$layout$var[su$layout$calib] %in% vars.valid]
		su$layout$calib <- su$layout$calib[su$layout$layout$var[su$layout$calib] %in% vars.calib]
	}
	su$layout$timestep.fac <- timestep.fac
	su$layout$time.units <- time.units
	if(grepl("mode", tag) & !grepl("mode", su$settings$f.likeli)){warning("tag contains \"mode\", but name of likelihood doesn't...")}
	su$layout$tme.orig <- tme.orig
	print("checking layout...")
	su$layout <- check.layout(su$layout)
	print("done")
	##remove time points from calibration with non-positive timestep
	##rmve <- su$layout$calib[which(diff(su$layout$layout$time[su$layout$calib])<=0 & diff(as.numeric(as.factor(su$layout$layout$var[su$layout$calib])))==0)]
	##while(length(rmve)>0){
	##	su$layout$calib <- su$layout$calib[!(su$layout$calib %in% rmve)]
	##	rmve <- su$layout$calib[which(diff(su$layout$layout$time[su$layout$calib])<=0 & diff(as.numeric(as.factor(su$layout$layout$var[su$layout$calib])))==0)]
	##}
	
	## set some fixed parameters and bounds
	## set some fixed parameters and bounds
	print("hard-coded change to following parameters:")
	print(su$model$parameters[grepl("U[2-9]W%Msmt(QE|E)_S0_[^SR]", names(su$model$parameters))])
	su$model$parameters[grepl("U[2-9]W%Msmt(QE|E)_S0_[^SR]", names(su$model$parameters))] <- -4.60517
	su$model$parameters[grepl("U[2-9]W%BetaE_UR", names(su$model$parameters))] <- 0.01 # this is the same as line above but not transformed
	if(!grepl("7h2y", sett)) su$model$parameters[grepl("Glo%Cmlt_Sb_FR", names(su$model$parameters))] <- -11.51293
	su$model$args$parHi[grepl("Glo%CmltSmax_IR", names(su$model$parameters))] <- log(51)
	if(grepl("Str07h2w",sett)) su$model$args$parHi[grepl("Glo%Cmlt_K_Qb_SR", names(su$model$parameters))] <- log(0.11)

	
    su$layout$calib <- su$layout$calib[!is.na(su$observations[su$layout$calib])]
	if(length(list.files(path=savepath, pattern="sol_.*RData"))==1){ #start sampling from result of optimizer
		load(paste(savepath,"/sol_",tag,".RData",sep=""))
		tran <- c(as.logical(su$model$args$parTran[as.logical(su$model$par.fit)]), as.logical(su$likelihood$tran[as.logical(su$likelihood$par.fit)]))
		prs <- ifelse(tran, log(sol[2:nrow(sol),1]), sol[2:nrow(sol),1])
		init.range <- cbind(prs-abs(prs)*0.05, prs+abs(prs)*0.05)
	}else{ # start sampling from inital parameter set
		if(!is.na(start.par[1])){# read from file
			ini <- read.table(paste(savepath,"/",start.par,sep=""), header=TRUE)
			np <- nrow(ini) - sum(grepl("_lik", ini[,1]))
			prs.mod <- subset(ini[1:np,], parameter %in% names(su$model$parameters)[as.logical(su$model$par.fit)])
			prs.lik <- subset(ini[(np+1):nrow(ini),], parameter %in% names(su$likelihood$parameters)[as.logical(su$likelihood$par.fit)])
			prs.mod[,2] <- ifelse(as.logical(su$model$args$parTran[which(names(su$model$parameters) %in% prs.mod[,"parameter"])]), log(prs.mod[,2]), prs.mod[,2])
			prs.lik[,2] <- ifelse(as.logical(su$likelihood$tran[which(names(su$likelihood$parameters) %in% prs.lik[,"parameter"])]), log(prs.lik[,2]), prs.lik[,2])
			prs.m <- su$model$parameters[as.logical(su$model$par.fit)]
			prs.m[which(names(prs.m) %in% prs.mod[,"parameter"])] <- prs.mod[,2]
			prs.l <- su$likelihood$parameters[as.logical(su$likelihood$par.fit)]
			prs.l[which(names(prs.l) %in% prs.lik[,"parameter"])] <- prs.lik[,2]
			prs <- c(prs.m,prs.l)
		}else{ # use default ones
			prs <- c(su$model$parameters[as.logical(su$model$par.fit)], su$likelihood$parameters[as.logical(su$likelihood$par.fit)])
		}
		init.range <- cbind(prs/vary, prs*vary)
	}
	init.range[prs==0,1] <- -1
	init.range[prs==0,2] <- 1
    init.range <- cbind(pmin(init.range[,1], init.range[,2]), pmax(init.range[,1], init.range[,2]))
    init.range[,1] <- pmax(init.range[,1], c(su$model$args$parLo[as.logical(su$model$par.fit)], su$likelihood$lb[as.logical(su$likelihood$par.fit)]))
    init.range[,2] <- pmin(init.range[,2], c(su$model$args$parHi[as.logical(su$model$par.fit)], su$likelihood$ub[as.logical(su$likelihood$par.fit)]))
	su$layout$pred.layout <- data.frame(var=rep(c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"),each=nrow(su$input$inputobs)), time=rep(su$input$inputobs[,1],3))
    P.roll <- rollmean(su$input$inputobs[,2], k=7, fill=0)
	P.roll[P.roll<1e-8] <- 0
    ##P.roll <- su$input$inputobs[,2]
    P.roll.calib <- interp.precip(su$input$inputobs[,1], P.roll, su$layout$layout)
    P.roll.pred  <- interp.precip(su$input$inputobs[,1], P.roll, su$layout$pred.layout)
    su$input$P.roll <- P.roll.calib
	su$input$P.roll.pred <- P.roll.pred	
	print("done preparing sudriv")
	
	##lyt <- su$layout
	##save(lyt, file=paste0(savepath,"/constr_layout.RData"))
	##obsr <- su$observations
	##save(obsr, file=paste0(savepath,"/constr_observations.RData"))
}else{
	##hyperfun <- su$model$hyperfun
	load(paste(savepath,"/su_",tag,".RData",sep=""))
	##su$model$hyperfun <- hyperfun #ATTENTION: overwrite the hyperfunction of the loaded object
	##pritag <- "Exa2"
	##su <- construct.prior(su, file=paste0("../sudriv_input/priordef_S07",pritag,".txt")) # ATTENTION overvrite prior of the loaded object
	##tag <- "v05_sprbckeros0wtns0_QE3P7c07i_TTE2b02c08i_SE1c08_priorDx4"
	##ATTENTION: update (overwrite) the likelihood and sampling functions of the loaded object
	##su$likelihood$f.likeli <- LogLikelihoodHydrology_la9esimp_fast_skewt
	##su$likelihood$f.sample <- LogLikelihoodHydrology_la9esimp_skewt_sample
	##su$model$args$nout <- 216
    ##outnames <- as.character(read.table("outnames.txt")[,1])
    ##outnames <- gsub("%", "_", outnames)
    ##outnames <- gsub("\\[", "", outnames)
    ##outnames <- gsub("\\]", "", outnames)
	##su$model$outnames <- outnames
	##su$model$args$parLo[grepl("U[2-9]T[1-2]%Sl.{3}_RR", names(su$model$parameters))] <- -15
	##su$model$parameters[grepl("U[2-9]T[1-2]%Sl.{3}_RR", names(su$model$parameters))] <- -12
	##su$model$args$parLo[grepl("U[2-9]T[1-2]%Sl.{3}_RR", names(su$model$parameters))] <- -15
	##su$model$parameters["GloTr%CmltKd_SR"] <- 11.7 # this is the value estimated from Atrazine conc in Rhine
	##su <- likelihood.setup(su, settings=sett, replace.param=TRUE)
	if(is.null(su$optimized)){
		if(!reload.timedep) init.range <- redef.init.range(su)
		##init.range["Glo%Cmlt_AlQq_FR",] <- c(0,0.1)
	}else{
		load(paste(savepath,"/sol_",tag,".RData",sep=""))
	}
	if(!is.null(su$parameter.sample)) init.state <- t(su$parameter.sample[dim(su$parameter.sample)[1],,])
	##init.state[,"Glo%Cmlt_AlQq_FR"] <- runif(nrow(init.state),0.9,1.1)


	## add pred.layout since this is new
	##su$likelihood$f.sample <- LogLikelihoodHydrology_la9esimp_skewt_sample
	#pred.layout <- subset(su$layout$layout, var=="C1Wv_Qstream")
	#a <- pred.layout
	#a$var <- "C1Tc1_Qstream"
	#pred.layout <- rbind(pred.layout, a)
	#a$var <- "C1Tc2_Qstream"
	#pred.layout <- rbind(pred.layout, a)
	#pred.layout <- su$layout$layout[su$layout$calib,]
	#su$layout$pred.layout <- pred.layout
	su$layout$pred.layout <- data.frame(var=rep(c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"),each=nrow(su$input$inputobs)), time=rep(su$input$inputobs[,1],3))
    P.roll <- rollmean(su$input$inputobs[,2], k=7, fill=0)
	P.roll[P.roll<1e-8] <- 0
    P.roll.pred  <- interp.precip(su$input$inputobs[,1], P.roll, su$layout$pred.layout)
	su$input$P.roll.pred <- P.roll.pred
	
	## set some fixed parameters and bounds
	print("hard-coded change to following parameters:")
	print(su$model$parameters[grepl("U[2-9]W%Msmt(QE|E)_S0_[^SR]", names(su$model$parameters))])
	su$model$parameters[grepl("U[2-9]W%Msmt(QE|E)_S0_[^SR]", names(su$model$parameters))] <- -4.60517 # this did not work because of Msmt"h"
	su$model$parameters[grepl("U[2-9]W%BetaE_UR", names(su$model$parameters))] <- 0.01 # this is the same as line above but not transformed
	if(!grepl("7h2y", sett)) su$model$parameters[grepl("Glo%Cmlt_Sb_FR", names(su$model$parameters))] <- -11.51293
	su$model$args$parHi[grepl("Glo%CmltSmax_IR", names(su$model$parameters))] <- log(51)
	if(grepl("Str07h2w",sett)) su$model$args$parHi[grepl("Glo%Cmlt_K_Qb_SR", names(su$model$parameters))] <- log(0.11)

}
if(!(grepl("A3",sett))){ #in this case, there is only one application date, since terb was not applied in this subcatchment
	if(place.terb.early){# move the 2nd application of terbuthylazine earlier to 3 June, 12:00
		col.app <- which(grepl("App_A[0-9]_Terb.3",colnames(su$input$inputobs)))[1]
		ind.terb2 <- which(su$input$inputobs[,col.app] > 0)
		if(length(ind.terb2) != 2){
			print("ind.terb2:")
			print(ind.terb2)
			stop("strange terb input encountered")
		}
		ind.terb2 <- ind.terb2[2]
		may28 <- which.min(abs(as.POSIXlt(su$layout$tme.orig) + su$input$inputobs[,"julday"]*60*60 - as.POSIXlt("2009-06-03 12:00:00")))
		if(ind.terb2 != may28){
			su$input$inputobs[may28,"P"] <- su$input$inputobs[ind.terb2,"P"]
			su$input$inputobs[may28,grepl("App_A[0-9]_Terb",colnames(su$input$inputobs))] <- su$input$inputobs[ind.terb2,grepl("App_A[0-9]_Terb",colnames(su$input$inputobs))]
			su$input$inputobs[ind.terb2,"P"] <- 0
			su$input$inputobs[ind.terb2,grepl("App_A[0-9]_Terb",colnames(su$input$inputobs))] <- 0
		}
	}else{ # move it to 10th of june
		col.app <- which(grepl("App_A[0-9]_Terb.1",colnames(su$input$inputobs)))[1]
		print(col.app)
		ind.terb2 <- which(su$input$inputobs[,col.app] > 0)
		if(length(ind.terb2) != 2) stop("strange terb input encountered")
		ind.terb2 <- ind.terb2[2]
		june10 <- which.min(abs(as.POSIXlt(su$layout$tme.orig) + su$input$inputobs[,"julday"]*60*60 - as.POSIXlt("2009-06-10 12:00:00")))
		if(ind.terb2 != june10){
			su$input$inputobs[june10,"P"] <- su$input$inputobs[ind.terb2,"P"]
			su$input$inputobs[june10,grepl("App_A[0-9]_Terb",colnames(su$input$inputobs))] <- su$input$inputobs[ind.terb2,grepl("App_A[0-9]_Terb",colnames(su$input$inputobs))]
			su$input$inputobs[ind.terb2,"P"] <- 0
			su$input$inputobs[ind.terb2,grepl("App_A[0-9]_Terb",colnames(su$input$inputobs))] <- 0
		}	
	}
}
if(sample.posterior){
    ##library(MCMCEnsembleSampler)
	print(cbind(names(c(su$model$parameters[as.logical(su$model$par.fit)], su$likelihood$parameters[as.logical(su$likelihood$par.fit)])), init.range))
    n.dim <- sum(c(su$model$par.fit, su$likelihood$par.fit))
    su <- s.m.mcmc.wrapper(logposterior, max.iter=max.iter, sudriv=su, savepath=savepath, tag=tag, init.range=init.range, init.state=NULL, iter.step=iter.step, drop=0.8, jitter=jitter, n.walkers=n.walkers, n.dim=n.dim, prior=TRUE, apprx=TRUE, weight.equally=FALSE, mvprior=TRUE)
	##dir.create(savepath, showWarnings=FALSE, recursive=TRUE)
    ##save(su, file=paste(savepath, "/su_", tag, ".RData", sep=""))
}else if(optimize.posterior){
	lb <- c(su$model$args$parLo[as.logical(su$model$par.fit)], su$likelihood$lb[as.logical(su$likelihood$par.fit)]) + 1e-3
	ub <- c(su$model$args$parHi[as.logical(su$model$par.fit)], su$likelihood$ub[as.logical(su$likelihood$par.fit)]) - 1e-3
	if(nparallel > 1){
		cl = makeCluster(nparallel, methods = FALSE)
		cl = makeCluster(nparallel, methods = FALSE)
		clusterEvalQ(cl, .libPaths("/local/ammannlo/personalR/library"))
		registerDoParallel(cl)
		clusterSetRNGStream(cl = cl, iseed = NULL) # Random number generator on the cluster
		optimized <- foreach(i=1:nparallel, .packages=c("nloptr", "truncnorm", "rjson", "skewt")) %dopar% {
			source("../../../../Likelihood/src/hydrolikeli_herb.r")
			source("../../../../Likelihood/src/plot.markov.r")
			su <- model.setup(su, settings=sett, writeO=TRUE, f.path.transK=NA)
			neg.logposterior <- function(...){-1*logposterior(...)}
			args <- list(eval_f           = neg.logposterior,
						lb               = lb,
						ub               = ub,
						 opts             = list(algorithm =  "NLOPT_G_MLSL_LDS",#"NLOPT_LN_BOBYQA",#"NLOPT_GN_DIRECT_L",
												 local_opts= list(algorithm="NLOPT_LN_NELDERMEAD", xtol_rel=1e-5, maxtime=10*60),
												 maxtime   = maxtime,
												 maxeval   = -1),
						sudriv           = su,
						prior            = TRUE,
						apprx            = TRUE,
						weight.equally   = FALSE)
			if(reload){
				x0 <- pmax(pmin(su$optimized[[i]]$solution, ub-1e-8), lb+1e-8)
			}else{
				x0 <- runif(nrow(init.range), init.range[,1], init.range[,2])
			}
			args.new <- c(list(x0=x0), args)
			do.call(nloptr, args=args.new)
		}
		stopCluster(cl)
	}else{
		neg.logposterior <- function(...){-1*logposterior(...)}
		x0 <- c(su$model$parameters[as.logical(su$model$par.fit)], su$likelihood$parameters[as.logical(su$likelihood$par.fit)])
		low  <- lb
		high <- ub
		## read initial point in parameter space
		if(reload){
			x0 <- sol[2:nrow(sol),1]
		}else{
			if(length(list.files(path=savepath, pattern="opt_par_ini"))==1){
				ini <- read.table(paste(savepath,"/opt_par_ini.txt",sep=""), header=TRUE)
				for(i in 1:sum(su$model$par.fit)){##find initial point
					x0[i] <- ini[ini[,1]==names(x0[i]),2]
				}
				for(j in 1:sum(su$likelihood$par.fit)){##find initial point
					x0[i+j] <- ini[ini[,1]==names(x0[i+j]),2]
				}
			}
		}
		if(reload | length(list.files(path=savepath, pattern="opt_par_ini"))==1){
			np <- sum(su$model$par.fit)
			x0[1:np] <- ifelse(as.logical(su$model$args$parTran[as.logical(su$model$par.fit)]), log(x0[1:np]), x0[1:np])
			x0[(np+1):length(x0)] <- ifelse(as.logical(su$likelihood$tran[as.logical(su$likelihood$par.fit)]), log(x0[(np+1):length(x0)]), x0[(np+1):length(x0)])
		}
		##print(cbind(init.range[,1],x0,init.range[,2]))
		##ll <- numeric(1000)
		##for(i in 1:10){
		##	x0 <- runif(length(x0), min=pmax(init.range[,1],low), max=pmin(init.range[,2],high))
		##	ll[i] <- neg.logposterior(x0=x0, sudriv=su, prior=TRUE, apprx=TRUE, weight.equally=FALSE)
		##}
		print(x0)
		optimized <- nloptr(x0               = x0,
						 eval_f           = neg.logposterior,
						 lb               = lb,
						 ub               = ub,
						 opts             = list(algorithm =  "NLOPT_GN_CRS2_LM",#"NLOPT_LN_BOBYQA",#"NLOPT_GN_DIRECT_L",
												 #local_opts= list(algorithm="NLOPT_LN_NELDERMEAD", xtol_rel=1e-5, maxtime=10*60),
												 maxtime   = maxtime,
												 maxeval   = -1,
												 xtol_rel = -1),
						 sudriv           = su,
						 prior            = TRUE,
						 apprx            = TRUE,
						 weight.equally   = FALSE)
		optimized <- list(optimized)
	}
	su <- select.maxlikpars.optimized(su, optimized)
	nmes <- names(c(su$model$parameters[as.logical(su$model$par.fit)], su$likelihood$parameters[as.logical(su$likelihood$par.fit)]))
	j <- 1
	sol <- matrix(nrow=length(su$optimized[[1]]$solution)+1, ncol=length(su$optimized))
	for(optim.curr in su$optimized){
		su$optimized[[j]]$solution <- optim.curr$solution
		sol.trans <- ifelse(as.logical(c(su$model$args$parTran[as.logical(su$model$par.fit)],su$likelihood$tran[as.logical(su$likelihood$par.fit)])), exp(optim.curr$solution), optim.curr$solution)
		names(sol.trans) <- nmes
		su$optimized[[j]]$solution.trans <- sol.trans
		sol[,j] <- c(optim.curr$objective, sol.trans)
		j <- j + 1
	}
	save(su, file=paste(savepath, "/su_opt_", tag, ".RData", sep = ""), version=2)
	rownames(sol) <- c("objective", names(su$optimized[[1]]$solution.trans))
	save(sol, file=paste(savepath, "/sol_", tag, ".RData", sep=""), version=2)
}

if(plot.chains){
	su <- select.maxlikpars(su)
	ind.fit <- as.logical(c(su$model$par.fit, su$likelihood$par.fit))
	a <- c(su$model$parameters, su$likelihood$parameters)[ind.fit]
	nm <- c(names(su$model$parameters), names(su$likelihood$parameters))[ind.fit]
	b <- ifelse(as.logical(c(su$model$args$parTran, su$likelihood$tran))[ind.fit], exp(a), a)
	write.table(data.frame(parameter=nm, value=as.numeric(b)), file=paste(savepath, "/opt_par_", tag, ".txt", sep=""), row.names=FALSE)
    pridef <- c(su$model$prior$distdef, su$likelihood$prior$distdef)
    lower.logpost <- NA
    res <- 300
    height <- 7 #inches
    width <- 9 #inches
    plot.markov.hist(sudriv=su, brn.in=brn.in, pridef=pridef,lower.logpost=lower.logpost, plot=TRUE, file.hist=paste(savepath, "/histpri", tag, ".png", sep=""), width=width*res, height=height*res, prior.only=TRUE, kl.div=FALSE, lab=modNam)
    #plot.markov.hist(sudriv=su, brn.in=brn.in, pridef=pridef,lower.logpost=lower.logpost, plot=TRUE, file.hist=paste(savepath, "/hist", tag, ".png", sep=""), width=width*res, height=height*res, prior.only=FALSE, kl.div=FALSE, lab=modNam)
    plot.markov.hist(sudriv=su, brn.in=brn.in, pridef=pridef,lower.logpost=lower.logpost, plot=TRUE, file.hist=paste(savepath, "/hist", modNam, ".pdf", sep=""), width=width, height=height, prior.only=FALSE, kl.div=FALSE, lab=modNam)
    #png(paste(savepath, "/cor", tag, ".png", sep=""), width=width*res, height=height*res, res=res*0.8)
    #plot.cor(sudriv=su, brn.in=brn.in, thin=thin, lower.logpost=lower.logpost)
    #dev.off()
    png(paste(savepath, "/chains", tag, ".png", sep=""), width=width*res, height=height*res, res=res)
    plot.markov.chain(sudriv=su, brn.in=brn.in, thin=1,  lower.logpost=lower.logpost )
    dev.off()
}

if(pred){
	su <- select.maxlikpars(su)
	su$predicted$det    <- sampling_wrapper(su, sample.par=FALSE, n.sample=1, sample.likeli=FALSE)
	su$predicted$sample.parunc <- sampling_wrapper(su, brn.in=brn.in, sample.par=TRUE, n.sample=n.sample, sample.likeli=FALSE)
	##su$predicted$sample <- sampling_wrapper(su, brn.in=brn.in, sample.par=TRUE, n.sample=n.sample, sample.likeli=TRUE)
    save(su, file=paste(savepath, ifelse(place.terb.early,"","/terb_late"), "/su_", tag, ".RData", sep=""), version=2)
}
if(calc.met){
	calc.metrics(su, xlim="calib", file.out=paste0(savepath, ifelse(place.terb.early,"","/terb_late"), "/metrics_",tag,".txt"), vars=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
}

if(summarize.results){
    ## plot stats in prediction period
    ## extend prediction to calibration period to make compatible with new plotting routines
    ## su$predicted$sample <- cbind(matrix(NA, nrow(su$predicted$sample), length(su$layout$calib)), su$predicted$sample)
    ## su$predicted$det <- cbind(matrix(NA, 1, length(su$layout$calib)), su$predicted$det)
    list.su <- list(E3=su)
    dat  <- pred.stats(list.su,tme.orig=su$layout$tme.orig)
	dat1 <- pred.stats(list.su[1],tme.orig=su$layout$tme.orig)
    ## dat <- tau.time(dat=dat, sudriv=su, obspred="obspred")

    ## plot predictions
    ## xlim1 <- su$layout$layout$time[su$layout$calib][round(11500/su$layout$timestep.fac)]#valid was 3500
    ## xlim <- c(xlim1, xlim1+900/su$layout$timestep.fac)
    pdf(paste(savepath, "/predictions4", tag, ".pdf", sep=""), height=9, width=12)
    plot.predictions(list.su, probs=c(0.05,0.95), xlim="calib", lp.num.pred=dat$lp.num.pred, tme.orig=su$layout$tme.orig, metrics=TRUE)
    plot.predictions(list.su, probs=c(0.05,0.95), xlim="pred", lp.num.pred=dat$lp.num.pred, tme.orig=su$layout$tme.orig, metrics=TRUE)
    ## plot.predictions(list.su, probs=c(0.05,0.95), xlim=xlim, lp.num.pred=dat$lp.num.pred, tme.orig=tme.orig)
    ## plot.predictions(list.su, probs=c(0.05,0.95), xlim=xlim, ylim=c(0,1))
    ## plot.predictions(list.su, probs=c(0.05,0.95), xlim=xlim, ylim=c(0,1), n.samp=1, arrange=c(E1=1,E2=2,E3a=3,E3ab=4, E4a=5))
    dev.off()
    calc.metrics(su, dat=dat1, xlim="calib", newline=FALSE, file1=paste(savepath, "/metrics", tag, ".txt", sep=""),file2=paste("../sudriv_output/table_results2_",su$settings$subcatchment,".txt",sep=""),errmod=tag, catchment=su$settings$subcatchment, append=TRUE)
    calc.metrics(su, dat=dat1, xlim="pred", newline=TRUE, file1=paste(savepath, "/metrics", tag, ".txt", sep=""), file2=paste("../sudriv_output/table_results2_",su$settings$subcatchment,".txt",sep=""), append=TRUE)

    pdf(paste(savepath, "/stats_valid", tag, ".pdf", sep=""), width=12)
    ## plot.ts.quantiles(dat, su, xlim=c(24000,25000))
    plot.ts.quantiles(dat, su, xlim="pred")
    ## plot.ts.tau(dat, su, xlim=c(4000,5000), ind.sel=su$layout$pred)
    plot.ts.tau(dat, su, xlim="pred")
    ## plot.ts.quantiles(dat, su, xlim=c(35000,40000))
    plot.dens.quantiles(dat, su, ind.sel=su$layout$pred, distr="unif")
    plot.dens.quantiles(dat, su, ind.sel=su$layout$pred)
    plot.dens.innovation(dat, su, ind.sel=su$layout$pred)
    plot.Qdet.quantiles(dat, su, ind.sel=su$layout$pred)
    plot.pacf.quantiles(dat, su, ind.sel=su$layout$pred, lag.max=10, confidence=0.95)
    plot.flashiness(dat, list.su, xlim="pred")
    plot.strmflw.err(dat, list.su, xlim="pred")
    plot.nse(dat, list.su, xlim="pred")
    plot.crps(dat, list.su, xlim="pred")
    ##plot.tau.time(dat, su, ind.sel=su$layout$pred)
    dev.off()
    pdf(paste(savepath, "/stats_calib", tag, ".pdf", sep=""), width=12)
    ## plot.ts.quantiles(dat, su, xlim=c(10000,11000))
    plot.ts.quantiles(dat, su, xlim="calib")
    ## plot.ts.tau(dat, su, ind.sel=su$layout$calib)
    plot.ts.tau(dat, su, xlim="calib")
    plot.powspec.quantiles(dat, su, ind.sel=su$layout$calib)
    ## plot.ts.quantiles(dat, su, xlim=c(35000,40000))
    plot.dens.quantiles(dat, su, ind.sel=su$layout$calib, distr="unif")
    plot.dens.quantiles(dat, su, ind.sel=su$layout$calib)
    plot.dens.innovation(dat, su, ind.sel=su$layout$calib)
    plot.Qdet.quantiles(dat, su, ind.sel=su$layout$calib)
    plot.sd(su, fsd=2)
    plot.pacf.quantiles(dat, su, ind.sel=su$layout$calib, lag.max=10, confidence=0.95)
    plot.flashiness(dat, list.su, xlim="calib")
    plot.strmflw.err(dat, list.su, xlim="calib")
    plot.nse(dat, list.su, xlim="calib")
    plot.crps(dat, list.su, xlim="calib")
    ## if prediction was made for the calibration phase:
    ## plot.flashiness(dat, list.su)
    ## plot.strmflw.err(dat, list.su)
    dev.off()
}

if(plot.states){
    allsub <- ifelse(su$settings$subcatchment == "AAll", TRUE, FALSE)
    if(!is.na(set.par[1])){
		x0 <- c(su$model$parameters[as.logical(su$model$par.fit)], su$likelihood$parameters[as.logical(su$likelihood$par.fit)])
		ini <- read.table(paste(savepath,"/",set.par,sep=""), header=TRUE)
		for(i in 1:sum(su$model$par.fit)){##find initial point
			x0[i] <- ini[ini[,1]==names(x0[i]),2]
		}
		for(j in 1:sum(su$likelihood$par.fit)){##find initial point
			x0[i+j] <- ini[ini[,1]==names(x0[i+j]),2]
		}
		np <- sum(su$model$par.fit)
		x0[1:np] <- ifelse(as.logical(su$model$args$parTran[as.logical(su$model$par.fit)]), log(x0[1:np]), x0[1:np])
		x0[(np+1):length(x0)] <- ifelse(as.logical(su$likelihood$tran[as.logical(su$likelihood$par.fit)]), log(x0[(np+1):length(x0)]), x0[(np+1):length(x0)])
		su$model$parameters[as.logical(su$model$par.fit)] <- x0[1:np]
		su$likelihood$parameters[as.logical(su$likelihood$par.fit)] <- x0[(np+1):length(x0)]
    }else{
		if(reload.timedep){
			load(paste0("../output/timedeppar/result_",tag,".RData"))
			if(!("result" %in% ls())) result <- res
			##result$param.maxpost <- result$param.maxpost[unlist(lapply(result$param.maxpost, length))<=1]
			su <- select.maxlikpars.timedep(su, result, scaleshift=scaleshift, lik.not.post=TRUE)
		}else{
			if(!is.null(su$parameter.sample)) su <- select.maxlikpars(su)
		}
    }
	print("selected parameters:")
	print(su$model$parameters[su$model$par.fit==1])
	##su$model$parameters["U3T1%SlTwo_SR"] <- 6.5##5.5
	##su$model$parameters["U3T1%SlOne_SR"] <- 6.5##5.5
	##su$model$parameters["GloTr%CmltSlOne_UR"] <- 6.5
	##su$model$parameters["GloTr%CmltSlTwo_UR"] <- 3
	##su$model$parameters[grepl("U[2-9]T[1-2]%SlOne_SR", names(su$model$parameters))] <- 10
	##su$model$parameters[grepl("U[2-9]T[1-2]%SlTwo_SR", names(su$model$parameters))] <- -2
	##su$model$args$parHi[grepl("U[2-9]T[1-2]%Kd_SR", names(su$model$parameters))] <- 4
	##su$model$hyperfun <- hyperfun #ATTENTION: overwrite the hyperfunction of the loaded object
	##su$model$parameters["GloTr%CmltKd_SR"] <- -9.54
	##print(su$model$parameters[grepl("U[2-9]T2%Kd_IR",names(su$model$parameters))])
	##su$model$parameters[grepl("U[2-9]T2%Kd_IR",names(su$model$parameters))] <- log(0.5)
	##su$model$parameters["U3T2%SlOne_SR"] <- 6.5##5.5
	##su$model$parameters["GloTr%CmltRs_WR"] <- log(4e-3)
	##su$model$parameters["GloTr%CmltSlOne_IR"] <- log(20)
	##su$model$parameters["GloTr%CmltSlTwo_IR"] <- log(60)
    ##su$model$parameters["Glo%CmltSmax_IR"] <- log(10)
    ##su$model$parameters["Glo%tStart_VL"] <- 1
	##su$model$args$parHi[names(su$model$parameters)=="U3W%K_Qb_SR"] <- log(10)
	##su$model$parameters["U3W%K_Qb_SR"] <- log(3.5)
	##su$model$parameters["U3W%K_Qq_SR"] <- log(2)
	##su$model$parameters["U3W%Alfa_Qq_SR"] <- log(2)
	##su$model$parameters["Glo%Cmlt_AlQq_FR"] <- log(1.5*1.20137439397018)
	##su$model$parameters["Glo%Cmlt_K_Qq_SR"] <- log(2*6.2e-06)
	##su$model$parameters["Glo%Cmlt_K_Qq_FR"] <- log(1.5*0.046536493860088)
	##su$model$parameters["U3W%Smax_UR"] <- log(1)
	##su$model$parameters["U4W%Dspl_RE"] <- 0.8
	##su$model$parameters["Glo%Cmlt_Dspl_SD"] <- 0.9
	##su$model$parameters["Glo%Cmlt_K_Qq_FR"] <- log(0.3)
	##su$model$args$parHi[names(su$model$parameters)=="Glo%Cmlt_K_Qq_RR"] <- log(4.5)
	##su$model$parameters["Glo%Cmlt_K_Qq_RR"] <- log(4.5)
	##su$model$parameters["Glo%Cmlt_Pmax_ED"] <- log(11.5)
	##su$model$parameters[grepl("Alfa_Qq_FR",names(su$model$parameters))] <- log(3)
	##su$model$parameters[grepl("K_Qb_FR",names(su$model$parameters))] <- 0.01
	##su$model$parameters["Glo%Cmlt_K_Qb_UR"] <- log(9e-6)
	##su$model$parameters["Glo%Cmlt_K_Qb_SR"] <- log(0.5e-2)
	##su$model$parameters["GloTr%CmltKd_WR"] <- log(1e-3)
	##load("../output/A1Str07h2b/v05_spraydrift_QE1c07i_TTE1c07i_priorC2/model.parameters.RData")
	##su$model$parameters <- model.parameters
    ##su$model$timedep$par[,1] <- pmin(su$model$timedep$par[,1],log(2))
    ##tag <- paste0(tag, "_max2")
    var.obs.q <- c("C1Wv_Qstream","C2Wv_Qstream","C3Wv_Qstream","C4Wv_Qstream","C5Wv_Qstream","P") ## observations and the model output of those are plotted
    var.obs.atra <- c("C1Tc1_Qstream", "C2Tc1_Qstream", "C3Tc1_Qstream", "C4Tc1_Qstream", "C5Tc1_Qstream")
    var.obs.terb <- c("C1Tc2_Qstream", "C2Tc2_Qstream", "C3Tc2_Qstream", "C4Tc2_Qstream", "C5Tc2_Qstream")
	if(grepl("SE[0-9]",sett)){
		var.obs.sorp <- c("HYPERSTATE_sorption")
	}
    if(!allsub){
		var.obs.q <- var.obs.q[1]
		var.obs.atra <- var.obs.atra[1]
		var.obs.terb <- var.obs.terb[1]
		var.obs.mult <- c(var.obs.q[1], var.obs.atra[1], var.obs.terb[1], "P")
    }else{
		var.obs.mult <- c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream","C3Wv_Qstream","C3Tc1_Qstream","C3Tc2_Qstream","P")
	}
	if(reload.timedep & !grepl("none", tag)) var.obs.mult <- c(var.obs.mult, "timedep")
    if(grepl("Str07", sett)){
    	if(allsub){
			var.mod1 <- c("U1F5Wv_Sf1", "U5F5Wv_Si1", "U5F5Wv_Su1") ## only the model output of those are plotted
			var.mod2 <- c("U2F5Wv_Sr1", "U3F2Wv_Sf1", "U3F4Wv_Ss1", "U5F5Wv_Ss1")
			if(grepl("Str072", sett)){
				var.mod2 <- c("U2F5Wv_Sr1", "U2F4Wv_Ss1", "U5F5Wv_Ss1")
			}
			var.modAcCon<- c("U2F1Tc1Lv1_Si1", "U2F1Tc1Lv2_Si1", "U2F1Tc1Lv1_Sr1", "U2F1Tc1Lv2_Sr1")
			var.modAcDrn<- c("U3F1Tc1Lv1_Si1", "U3F1Tc1Lv2_Si1", "U3F1Tc1Lv1_Sf1", "U3F1Tc1Lv2_Sf1")
			var.modTcCon<- c("U2F1Tc2Lv1_Si1", "U2F1Tc2Lv2_Si1", "U2F1Tc2Lv1_Sr1", "U2F1Tc2Lv2_Sr1")
			var.modTcDrn<- c("U3F1Tc2Lv1_Si1", "U3F1Tc2Lv2_Si1", "U3F1Tc2Lv1_Sf1", "U3F1Tc2Lv2_Sf1")
			var.modTcGwt<- c("U3F1Tc1Lv1_Ss1", "U3F1Tc1Lv2_Ss1", "U5F1Tc1Lv1_Ss1", "U5F1Tc1Lv2_Ss1")
			if(grepl("Str071", sett)){
				var.modTcCon<- c("U2F1Tc1Lv1_Si1", "U2F1Tc2Lv1_Si1")
				var.modTcDrn<- c("U3F1Tc1Lv1_Si1", "U3F1Tc2Lv1_Si1")
				var.modTcGwt<- c("U5F1Tc1Lv1_Ss1", "U5F1Tc2Lv1_Ss1")
			}
			var.mod3 <- c("U1F2Tm1_Qstrm", "U2F2Tm1_Qstrm", "U3F2Tm1_Qstrm", "U4F2Tm1_Qstrm", "U5F2Tm1_Qstrm")
			var.mod4 <- c("U1F2Tm2_Qstrm", "U2F2Tm2_Qstrm", "U3F2Tm2_Qstrm", "U4F2Tm2_Qstrm", "U5F2Tm2_Qstrm")
			var.mod5 <- c("U1F2Wv_Qstrm", "U2F2Wv_Qstrm", "U3F2Wv_Qstrm", "U4F2Wv_Qstrm", "U5F2Wv_Qstrm")
			var.modM <- c("U3F2Wv_Si1", "U3F2Tc1Lv1_Si1", "U3F2Tc2Lv1_Si1", "U3F2Tc1Lv2_Si1", "U3F2Tc2Lv2_Si1")
			var.modM2 <- c("U2F2Wv_Si1", "U2F2Tc1Lv1_Si1", "U2F2Tc2Lv1_Si1", "U2F2Tc1Lv2_Si1", "U2F2Tc2Lv2_Si1")
			var.modM5 <- c("U5F2Wv_Si1", "U5F2Tc1Lv1_Si1", "U5F2Tc2Lv1_Si1", "U5F2Tc1Lv2_Si1", "U5F2Tc2Lv2_Si1")
			var.modMU<- c("U3F2Wv_Su1", "U3F2Tc1Lv1_Su1", "U3F2Tc2Lv1_Su1", "U3F2Tc1Lv2_Su1", "U3F2Tc2Lv2_Su1")
			var.modMS<- c("U3F2Wv_Ss1", "U3F2Tc1Lv1_Ss1", "U3F2Tc2Lv1_Ss1", "U3F2Tc1Lv2_Ss1", "U3F2Tc2Lv2_Ss1")
			var.modQind<- c("U1F2Wv_Qq_FR", "U2F2Wv_Q_RR", "U3F2Wv_Qq_FR", "U2F2Wv_Qq_SR", "U3F2Wv_Qq_SR", "U4F2Wv_Qq_SR", "U5F2Wv_Qq_SR")
			var.modQsur<- c("U1F2Wv_Qq_FR", "U2F2Wv_Qe_UR", "U3F2Wv_Qe_UR")
			var.modAind<- c("U1F2Tm1_Qq_FR", "U2F2Tm1_Q_RR", "U3F2Tm1_Qq_FR", "U2F2Tm1_Qq_SR", "U3F2Tm1_Qq_SR", "U4F2Tm1_Qq_SR", "U5F2Tm1_Qq_SR")
			var.modAsur<- c("U2F2Tm1_Qe_UR", "U3F2Tm1_Qe_UR")
			var.modTind<- c("U1F2Tm2_Qq_FR", "U2F2Tm2_Q_RR", "U3F2Tm2_Qq_FR", "U2F2Tm2_Qq_SR", "U3F2Tm2_Qq_SR", "U4F2Tm2_Qq_SR", "U5F2Tm2_Qq_SR")
			var.modTsur<- c("U2F2Tm2_Qe_UR", "U3F2Tm2_Qe_UR")
			pp.mod1  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod1)
			pp.mod2  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod2)
			pp.modAcCon <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modAcCon)
			pp.modAcDrn <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modAcDrn)
			pp.modTcCon <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTcCon)
			pp.modTcDrn <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTcDrn)
			pp.modTcGwt <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTcGwt)
			pp.mod3  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod3)
			pp.mod4  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod4)
			pp.mod5  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod5)
			pp.modM  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modM)
			pp.modM2  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modM2)
			pp.modM5  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modM5)
			pp.modMU <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modMU)
			pp.modMS <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modMS)
			pp.modQind <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modQind)
			pp.modQsur <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modQsur)
			pp.modAind <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modAind)
			pp.modTind <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTind)
			pp.modAsur <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modAsur)
			pp.modTsur <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTsur)
		}else{
			proxim <- grepl("Str076",sett)
			noRR <- grepl("Str070",sett) | grepl("Str074",sett) | grepl("Str0725",sett)
			noFR <- grepl("Str0726",sett)
			var.mod1 <- c("U1F1Wv_Sf1", "U3F1Wv_Si1", "U3F1Wv_Su1") ## only the model output of those are plotted
			var.mod2 <- c(switch(noRR+1, "U2F1Wv_Sr1", NULL), switch(noFR+1, "U3F1Wv_Sf1", NULL), "U3F1Wv_Ss1", "U3F1Wv_Ss1")
			##var.mod25 <- c("U2F1Wv_Su1", "U4F1Wv_Su1", "U2F1Wv_Sr1", "U4F1Wv_Sr1")
			if(grepl("Str072ts", sett)){
				var.mod2 <- c("U2F1Wv_Sf1", "U2F1Wv_Ss1", "U5F1Wv_Ss1")
			}
			var.modAcCon<- c("U2F1Tc1Lv1_Si1", "U2F1Tc1Lv2_Si1", "U2F1Tc1Lv1_Su1", "U2F1Tc1Lv2_Su1")
			var.modAcDrn<- c("U3F1Tc1Lv1_Si1", "U3F1Tc1Lv2_Si1", "U3F1Tc1Lv1_Su1", "U3F1Tc1Lv2_Su1")
			var.modAcConDrn<- c("U4F1Tc1Lv1_Si1", "U4F1Tc1Lv2_Si1", "U4F1Tc1Lv1_Su1", "U4F1Tc1Lv2_Su1")
			var.modTcCon<- c("U2F1Tc2Lv1_Si1", "U2F1Tc2Lv2_Si1", "U2F1Tc2Lv1_Su1", "U2F1Tc2Lv2_Su1")
			var.modTcDrn<- c("U3F1Tc2Lv1_Si1", "U3F1Tc2Lv2_Si1", "U3F1Tc2Lv1_Su1", "U3F1Tc2Lv2_Su1")
			var.modTcConDrn<- c("U4F1Tc2Lv1_Si1", "U4F1Tc2Lv2_Si1", "U4F1Tc2Lv1_Su1", "U4F1Tc2Lv2_Su1")
			var.modTcGwt<- c("U3F1Tc1Lv1_Ss1", "U3F1Tc1Lv2_Ss1", ifelse(proxim,"U4F1Tc1Lv1_Ss1","U5F1Tc1Lv1_Ss1"), ifelse(proxim,"U4F1Tc1Lv2_Ss1","U5F1Tc1Lv2_Ss1"))
			var.mod3 <- c("U1F1Tm1_Qstrm", "U2F1Tm1_Qstrm", "U3F1Tm1_Qstrm", "U4F1Tm1_Qstrm", switch(proxim+1, "U5F1Tm1_Qstrm", NULL))
			var.mod4 <- c("U1F1Tm2_Qstrm", "U2F1Tm2_Qstrm", "U3F1Tm2_Qstrm", "U4F1Tm2_Qstrm", switch(proxim+1, "U5F1Tm2_Qstrm", NULL))
			var.mod5 <- c("U1F1Wv_Qstrm", "U2F1Wv_Qstrm", "U3F1Wv_Qstrm", "U4F1Wv_Qstrm", switch(proxim+1,"U5F1Wv_Qstrm", NULL), "C1Wv_Qstream")
			var.modM <- c("U3F1Wv_Si1", "U3F1Tc1Lv1_Si1", "U3F1Tc2Lv1_Si1", "U3F1Tc1Lv2_Si1", "U3F1Tc2Lv2_Si1")
			var.modM2 <- c("U2F1Wv_Si1", "U2F1Tc1Lv1_Si1", "U2F1Tc2Lv1_Si1", "U2F1Tc1Lv2_Si1", "U2F1Tc2Lv2_Si1")
			var.modM5 <- c("U5F1Wv_Si1", "U5F1Tc1Lv1_Si1", "U5F1Tc2Lv1_Si1", "U5F1Tc1Lv2_Si1", "U5F1Tc2Lv2_Si1")
			var.modMU<- c("U3F1Wv_Su1", "U3F1Tc1Lv1_Su1", "U3F1Tc2Lv1_Su1", "U3F1Tc1Lv2_Su1", "U3F1Tc2Lv2_Su1")
			var.modMS<- c("U3F1Wv_Ss1", "U3F1Tc1Lv1_Ss1", "U3F1Tc2Lv1_Ss1", "U3F1Tc1Lv2_Ss1", "U3F1Tc2Lv2_Ss1")
			if(grepl("Str071", sett)){ # consider that Lv2 is not existing
				var.modTcCon<- c("U2F1Tc1Lv1_Si1", "U2F1Tc2Lv1_Si1")
				var.modTcDrn<- c("U3F1Tc1Lv1_Si1", "U3F1Tc2Lv1_Si1")
				var.modTcGwt<- c("U5F1Tc1Lv1_Ss1", "U5F1Tc2Lv1_Ss1")
				var.modM <- c("U3F1Wv_Si1", "U3F1Tc1Lv1_Si1", "U3F1Tc2Lv1_Si1")
				var.modM2 <- c("U2F1Wv_Si1", "U2F1Tc1Lv1_Si1", "U2F1Tc2Lv1_Si1")
				var.modM5 <- c("U5F1Wv_Si1", "U5F1Tc1Lv1_Si1", "U5F1Tc2Lv1_Si1")
				var.modMU<- c("U3F1Wv_Su1", "U3F1Tc1Lv1_Su1", "U3F1Tc2Lv1_Su1")
				var.modMS<- c("U3F1Wv_Ss1", "U3F1Tc1Lv1_Ss1", "U3F1Tc2Lv1_Ss1")
			}
			if(proxim) var.modM5 <- gsub("U5","U4",var.modM5)
			var.modQind<- c("U1F1Wv_Qq_FR", switch(noRR+1, "U2F1Wv_Q_RR", NULL), switch(noFR+1,"U3F1Wv_Qq_FR", NULL), "U2F1Wv_Qq_SR", "U3F1Wv_Qq_SR", "U4F1Wv_Qq_SR", switch(proxim+1,"U5F1Wv_Qq_SR", NULL))
			var.modQsur<- c("U1F1Wv_Qq_FR", "U2F1Wv_Qe_UR", "U3F1Wv_Qe_UR")
			var.modAind<- c("U1F1Tm1_Qq_FR", switch(noRR+1, "U2F1Tm1_Q_RR", NULL), switch(noFR+1,"U3F1Tm1_Qq_FR", NULL), "U2F1Tm1_Qq_SR", "U3F1Tm1_Qq_SR", "U4F1Tm1_Qq_SR", switch(proxim+1,"U5F1Tm1_Qq_SR", NULL))
			var.modAsur<- c("U2F1Tm1_Qe_UR", "U3F1Tm1_Qe_UR")
			var.modTind<- c("U1F1Tm2_Qq_FR", switch(noRR+1,"U2F1Tm2_Q_RR", NULL), switch(noFR+1,"U3F1Tm2_Qq_FR", NULL), "U2F1Tm2_Qq_SR", "U3F1Tm2_Qq_SR", "U4F1Tm2_Qq_SR", switch(proxim+1,"U5F1Tm2_Qq_SR", NULL))
			var.modTsur<- c("U2F1Tm2_Qe_UR", "U3F1Tm2_Qe_UR")
			pp.mod1  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod1)
			pp.mod2  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod2)
			##pp.mod25 <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod25)
			if(!grepl("Str071", sett)){
				pp.modAcCon <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modAcCon)
				pp.modAcDrn <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modAcDrn)
				pp.modAcConDrn <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modAcConDrn)
				pp.modTcConDrn <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTcConDrn)
			}
			pp.modTcCon <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTcCon)
			pp.modTcDrn <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTcDrn)
			pp.modTcGwt <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTcGwt)
			pp.mod3  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod3)
			pp.mod4  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod4)
			pp.mod5  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod5)
			pp.modM  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modM)
			pp.modM2  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modM2)
			pp.modM5  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modM5)
			pp.modMU <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modMU)
			pp.modMS <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modMS)
			pp.modQind <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modQind)
			pp.modQsur <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modQsur)
			pp.modAind <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modAind)
			pp.modTind <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTind)
			pp.modAsur <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modAsur)
			pp.modTsur <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modTsur)
		}
    }else if(grepl("Str02", sett)){
        var.mod <- c("U1F1Wv_Sf1", "U2F1Wv_Su1", "U2F1Wv_Ss1") ## only the model output of those are plotted
		var.modT<- c("U2F1Tc1Lv1_Su1", "U2F1Tc2Lv1_Su1")
        pp.mod  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.mod)
        pp.modT  <- prepare.plot.layout(sudriv=su, var.obs=c("C1Wv_Qstream"), var.mod=var.modT)
	}
    E00 <- as.POSIXct(c("2009-05-19 00:00", "2009-05-20 12:00"))
    E0 <- as.POSIXct(c("2009-05-21 12:00", "2009-05-23 00:00"))
	E1 <- as.POSIXct(c("2009-05-26 12:00", "2009-05-28 00:00"))
    E1long <- as.POSIXct(c("2009-05-26 12:00", "2009-06-05 00:00"))
    E1tail <- as.POSIXct(c("2009-05-28 00:00", "2009-06-05 00:00"))
    E2 <- as.POSIXct(c("2009-06-06 00:00", "2009-06-07 18:00"))
    E2long <- as.POSIXct(c("2009-06-04 00:00", "2009-06-09 18:00"))
    E25 <- as.POSIXct(c("2009-06-19 00:00", "2009-06-21 00:00"))
    E255 <- as.POSIXct(c("2009-06-23 00:00", "2009-06-26 00:00"))
    E3 <- as.POSIXct(c("2009-06-25 00:00", "2009-06-30 00:00"))
    E35<- as.POSIXct(c("2009-06-30 00:00", "2009-07-10 00:00"))
	E4 <- as.POSIXct(c("2009-07-10 00:00", "2009-07-25 18:00"))
    Ebeg <- as.POSIXct(c("2009-03-01 00:00", "2009-05-15 00:00"))
    Etmp <- as.POSIXct(c("2009-05-01 00:00", "2009-05-15 00:00"))
    Etmp2 <- as.POSIXct(c("2009-08-01 00:00", "2009-08-15 00:00"))
    Etmp3 <- as.POSIXct(c("2009-08-03 12:00", "2009-08-04 00:00"))
    Eapp <- as.POSIXct(c("2009-05-15 00:00", "2009-07-30 00:00"))
	Eapplong <- c(as.POSIXct("2009-05-01"),Eapp[2])
	Eapplonglong <- c(as.POSIXct("2009-05-01"),"2009-10-30 00:00")
    Eapp.tobi <- as.POSIXct(c("2009-05-20 00:00", "2009-07-23 00:00"))	
	Eapp <- Eapp.tobi
    Ewint <- as.POSIXct(c("2008-10-01 00:00", "2009-03-31 00:00"))
	## xlim <- as.POSIXct(c("2026-03-01", "2026-08-01"))
    pp.q <- prepare.plot.layout(sudriv=su, var.obs=var.obs.q)
    pp.atra <- prepare.plot.layout(sudriv=su, var.obs=var.obs.atra)
    pp.terb <- prepare.plot.layout(sudriv=su, var.obs=var.obs.terb)
	pp.mult <- prepare.plot.layout(sudriv=su, var.obs=var.obs.mult)
	if(grepl("SE[0-9]",sett)){
		pp.sorp <- prepare.plot.layout(sudriv=su, var.obs=var.obs.sorp)
	}
	if(is.null(su$model$hru.areas)){
		hru.areas <- NA
	}else{
		hru.areas <- su$model$hru.areas
	}
    # hru.areas = matrix(c(5, 12776, 24512, 3920, 130057.8,
						# 0, 11588, 40852, 3148, 69701.53,
						# 7660, 10632, 126820, 2136, 99174.42,
						# 2230, 19472, 54300, 4828, 86515.46,
						# 7515, 48752, 50160, 28556, 337583.4),byrow=TRUE,ncol=5)
	# hru.areas <- matrix(colSums(hru.areas), ncol=5)
	ind.fit <- as.logical(c(su$model$par.fit, su$likelihood$par.fit))
	a <- c(su$model$parameters, su$likelihood$parameters)[ind.fit]
	nm <- c(names(su$model$parameters), names(su$likelihood$parameters))[ind.fit]
	b <- ifelse(as.logical(c(su$model$args$parTran, su$likelihood$tran))[ind.fit], exp(a), a)
	write.table(data.frame(parameter=nm, value=as.numeric(b)), file=paste(savepath, "/opt_par_", tag, ".txt", sep=""), row.names=FALSE)
	
    translate.var <- c("C1Wv_Qstream", "P", "C2Wv_Qstream", "C3Wv_Qstream", "C4Wv_Qstream", "C5Wv_Qstream", "C1Tc1_Qstream", "C1Tc2_Qstream", "C2Tc1_Qstream", "C2Tc2_Qstream", "C3Tc1_Qstream", "C3Tc2_Qstream", "C4Tc1_Qstream", "C4Tc2_Qstream", "C5Tc1_Qstream", "C5Tc2_Qstream", "U1F1Wv_Qstrm", "U2F1Wv_Qstrm", "U3F1Wv_Qstrm", "U4F1Wv_Qstrm", "U5F1Wv_Qstrm", "U1F1Tm1_Qstrm", "U2F1Tm1_Qstrm", "U3F1Tm1_Qstrm", "U4F1Tm1_Qstrm", "U5F1Tm1_Qstrm", "U1F1Tm2_Qstrm", "U2F1Tm2_Qstrm", "U3F1Tm2_Qstrm", "U4F1Tm2_Qstrm", "U5F1Tm2_Qstrm", "U1F1Wv_Qq_FR", "U2F1Wv_Q_RR", "U3F1Wv_Qq_FR", "U1F1Tm1_Qq_FR", "U2F1Tm1_Q_RR", "U3F1Tm1_Qq_FR", "U1F1Tm2_Qq_FR", "U2F1Tm2_Q_RR", "U3F1Tm2_Qq_FR", "U3F1MaT1_Si1Lv0", "U3F1MaT1_Si1Lv1", "U3F1MaT1_Si1Lv2", "U3F1MaT2_Si1Lv0", "U3F1MaT2_Si1Lv1", "U3F1MaT2_Si1Lv2", "HYPERSTATE_sorption")
    translate.to <- c("C1Wv_Qstream", expression("Precip."~(mm/15*"min")), "C2Wv", "C3Wv", "C4Wv", "C5Wv", expression("Atrazine"~(mu*g/l)), expression("Terbuthylazine"~(mu*g/l)), "C2Tc1", "C2Tc2", "C3Tc1", "C3Tc2", "C4Tc1", "C4Tc2", "C5Tc1", "C5Tc2", "Impervious", "Shortcut", "Drained", ifelse(grepl("Str076", struct), "Rest", "SC and Drained"), "Rest", "Impervious", "Shortcut", "Drained", ifelse(grepl("Str076", struct), "Rest", "SC and Drained"), "Rest", "Impervious", "Shortcut", "Drained", ifelse(grepl("Str076", struct), "Rest", "SC and Drained"), "Rest", "Impervious reservoir", "Shortcut reservoir", "Drainage reservoir", "Impervious", "Shortcut", "Drainage",  "Impervious", "Shortcut", "Drainage", "Dissolved", "Fast sorbed", "Slow sorbed", "Dissolved", "Fast sorbed", "Slow sorbed", "Apparent~K[d]~(l/kg)")
	
    plot.results(layout.mod=pp.q$layout.mod$layout, y.mod=pp.q$y.mod, layout.obs=pp.q$layout.obs, y.obs=pp.q$y.obs, variables=var.obs.q, file=paste(savepath,"/states_q_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.q$layout.mod$layout, y.mod=pp.q$y.mod, layout.obs=pp.q$layout.obs, y.obs=pp.q$y.obs, variables=var.obs.q, file=paste(savepath,"/states_q_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.q$layout.mod$layout, y.mod=pp.q$y.mod, layout.obs=pp.q$layout.obs, y.obs=pp.q$y.obs, variables=var.obs.q, file=paste(savepath,"/states_q_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.q$layout.mod$layout, y.mod=pp.q$y.mod, layout.obs=pp.q$layout.obs, y.obs=pp.q$y.obs, variables=var.obs.q, file=paste(savepath,"/states_q_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.q$layout.mod$layout, y.mod=pp.q$y.mod, layout.obs=pp.q$layout.obs, y.obs=pp.q$y.obs, variables=var.obs.q, file=paste(savepath,"/states_q_E25_", tag, ".pdf", sep=""), xlim=E25, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.q$layout.mod$layout, y.mod=pp.q$y.mod, layout.obs=pp.q$layout.obs, y.obs=pp.q$y.obs, variables=var.obs.q, file=paste(savepath,"/states_q_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.q$layout.mod$layout, y.mod=pp.q$y.mod, layout.obs=pp.q$layout.obs, y.obs=pp.q$y.obs, variables=var.obs.q, file=paste(savepath,"/states_q_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.q$layout.mod$layout, y.mod=pp.q$y.mod, layout.obs=pp.q$layout.obs, y.obs=pp.q$y.obs, variables=var.obs.q, file=paste(savepath,"/states_q_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.atra$layout.mod$layout, y.mod=pp.atra$y.mod, layout.obs=pp.atra$layout.obs, y.obs=pp.atra$y.obs, variables=var.obs.atra, file=paste(savepath,"/states_atra_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.atra$layout.mod$layout, y.mod=pp.atra$y.mod, layout.obs=pp.atra$layout.obs, y.obs=pp.atra$y.obs, variables=var.obs.atra, file=paste(savepath,"/states_atra_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.atra$layout.mod$layout, y.mod=pp.atra$y.mod, layout.obs=pp.atra$layout.obs, y.obs=pp.atra$y.obs, variables=var.obs.atra, file=paste(savepath,"/states_atra_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.atra$layout.mod$layout, y.mod=pp.atra$y.mod, layout.obs=pp.atra$layout.obs, y.obs=pp.atra$y.obs, variables=var.obs.atra, file=paste(savepath,"/states_atra_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.atra$layout.mod$layout, y.mod=pp.atra$y.mod, layout.obs=pp.atra$layout.obs, y.obs=pp.atra$y.obs, variables=var.obs.atra, file=paste(savepath,"/states_atra_E25_", tag, ".pdf", sep=""), xlim=E25, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.atra$layout.mod$layout, y.mod=pp.atra$y.mod, layout.obs=pp.atra$layout.obs, y.obs=pp.atra$y.obs, variables=var.obs.atra, file=paste(savepath,"/states_atra_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.atra$layout.mod$layout, y.mod=pp.atra$y.mod, layout.obs=pp.atra$layout.obs, y.obs=pp.atra$y.obs, variables=var.obs.atra, file=paste(savepath,"/states_atra_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.atra$layout.mod$layout, y.mod=pp.atra$y.mod, layout.obs=pp.atra$layout.obs, y.obs=pp.atra$y.obs, variables=var.obs.atra, file=paste(savepath,"/states_atra_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.terb$layout.mod$layout, y.mod=pp.terb$y.mod, layout.obs=pp.terb$layout.obs, y.obs=pp.terb$y.obs, variables=var.obs.terb, file=paste(savepath,"/states_terb_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.terb$layout.mod$layout, y.mod=pp.terb$y.mod, layout.obs=pp.terb$layout.obs, y.obs=pp.terb$y.obs, variables=var.obs.terb, file=paste(savepath,"/states_terb_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.terb$layout.mod$layout, y.mod=pp.terb$y.mod, layout.obs=pp.terb$layout.obs, y.obs=pp.terb$y.obs, variables=var.obs.terb, file=paste(savepath,"/states_terb_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.terb$layout.mod$layout, y.mod=pp.terb$y.mod, layout.obs=pp.terb$layout.obs, y.obs=pp.terb$y.obs, variables=var.obs.terb, file=paste(savepath,"/states_terb_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_Etmp_", tag, ".pdf", sep=""), xlim=Etmp, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_Etmp2_", tag, ".pdf", sep=""), xlim=Etmp2, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_Etmp3_", tag, ".pdf", sep=""), xlim=Etmp3, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_Eapplong_", tag, ".pdf", sep=""), xlim=Eapplong, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_Eapplonglong_", tag, ".pdf", sep=""), xlim=Eapplonglong, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_Ebeg_", tag, ".pdf", sep=""), xlim=Ebeg, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E1long_", tag, ".pdf", sep=""), xlim=E1long, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E1tail_", tag, ".pdf", sep=""), xlim=E1tail, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E2long_", tag, ".pdf", sep=""), xlim=E2long, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E25_", tag, ".pdf", sep=""), xlim=E25, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E255_", tag, ".pdf", sep=""), xlim=E255, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E35_", tag, ".pdf", sep=""), xlim=E35, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    plot.results(layout.mod=pp.mult$layout.mod$layout, y.mod=pp.mult$y.mod, layout.obs=pp.mult$layout.obs, y.obs=pp.mult$y.obs, variables=var.obs.mult, file=paste(savepath,"/states_mult_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
	if(grepl("SE[0-9]",sett)){
		plot.results(layout.mod=pp.sorp$layout.mod$layout, y.mod=pp.sorp$y.mod, layout.obs=pp.sorp$layout.obs, y.obs=pp.sorp$y.obs, variables=var.obs.sorp, file=paste(savepath,"/states_sorp_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    }
	if(grepl("Str07", sett)){
		plot.results(layout.mod=pp.modM2$layout.mod$layout, y.mod=pp.modM2$y.mod, layout.obs=pp.modM2$layout.obs, y.obs=pp.modM2$y.obs, variables=var.modM2, file=paste(savepath,"/states_massIRU2_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modM$layout.mod$layout, y.mod=pp.modM$y.mod, layout.obs=pp.modM$layout.obs, y.obs=pp.modM$y.obs, variables=var.modM, file=paste(savepath,"/states_massIRU3_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modM5$layout.mod$layout, y.mod=pp.modM5$y.mod, layout.obs=pp.modM5$layout.obs, y.obs=pp.modM5$y.obs, variables=var.modM5, file=paste(savepath,"/states_massIRU5_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modM$layout.mod$layout, y.mod=pp.modM$y.mod, layout.obs=pp.modM$layout.obs, y.obs=pp.modM$y.obs, variables=var.modM, file=paste(savepath,"/states_massIR_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modM$layout.mod$layout, y.mod=pp.modM$y.mod, layout.obs=pp.modM$layout.obs, y.obs=pp.modM$y.obs, variables=var.modM, file=paste(savepath,"/states_massIR_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modM$layout.mod$layout, y.mod=pp.modM$y.mod, layout.obs=pp.modM$layout.obs, y.obs=pp.modM$y.obs, variables=var.modM, file=paste(savepath,"/states_massIR_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modM$layout.mod$layout, y.mod=pp.modM$y.mod, layout.obs=pp.modM$layout.obs, y.obs=pp.modM$y.obs, variables=var.modM, file=paste(savepath,"/states_massIR_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modM$layout.mod$layout, y.mod=pp.modM$y.mod, layout.obs=pp.modM$layout.obs, y.obs=pp.modM$y.obs, variables=var.modM, file=paste(savepath,"/states_massIR_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMU$layout.mod$layout, y.mod=pp.modMU$y.mod, layout.obs=pp.modMU$layout.obs, y.obs=pp.modMU$y.obs, variables=var.modMU, file=paste(savepath,"/states_massUR_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMU$layout.mod$layout, y.mod=pp.modMU$y.mod, layout.obs=pp.modMU$layout.obs, y.obs=pp.modMU$y.obs, variables=var.modMU, file=paste(savepath,"/states_massUR_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMU$layout.mod$layout, y.mod=pp.modMU$y.mod, layout.obs=pp.modMU$layout.obs, y.obs=pp.modMU$y.obs, variables=var.modMU, file=paste(savepath,"/states_massUR_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMU$layout.mod$layout, y.mod=pp.modMU$y.mod, layout.obs=pp.modMU$layout.obs, y.obs=pp.modMU$y.obs, variables=var.modMU, file=paste(savepath,"/states_massUR_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMU$layout.mod$layout, y.mod=pp.modMU$y.mod, layout.obs=pp.modMU$layout.obs, y.obs=pp.modMU$y.obs, variables=var.modMU, file=paste(savepath,"/states_massUR_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMU$layout.mod$layout, y.mod=pp.modMU$y.mod, layout.obs=pp.modMU$layout.obs, y.obs=pp.modMU$y.obs, variables=var.modMU, file=paste(savepath,"/states_massUR_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMS$layout.mod$layout, y.mod=pp.modMS$y.mod, layout.obs=pp.modMS$layout.obs, y.obs=pp.modMS$y.obs, variables=var.modMS, file=paste(savepath,"/states_massSR_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMS$layout.mod$layout, y.mod=pp.modMS$y.mod, layout.obs=pp.modMS$layout.obs, y.obs=pp.modMS$y.obs, variables=var.modMS, file=paste(savepath,"/states_massSR_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMS$layout.mod$layout, y.mod=pp.modMS$y.mod, layout.obs=pp.modMS$layout.obs, y.obs=pp.modMS$y.obs, variables=var.modMS, file=paste(savepath,"/states_massSR_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMS$layout.mod$layout, y.mod=pp.modMS$y.mod, layout.obs=pp.modMS$layout.obs, y.obs=pp.modMS$y.obs, variables=var.modMS, file=paste(savepath,"/states_massSR_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMS$layout.mod$layout, y.mod=pp.modMS$y.mod, layout.obs=pp.modMS$layout.obs, y.obs=pp.modMS$y.obs, variables=var.modMS, file=paste(savepath,"/states_massSR_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modMS$layout.mod$layout, y.mod=pp.modMS$y.mod, layout.obs=pp.modMS$layout.obs, y.obs=pp.modMS$y.obs, variables=var.modMS, file=paste(savepath,"/states_massSR_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, parameters=su$model$parameters, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modQind$layout.mod$layout, y.mod=pp.modQind$y.mod, layout.obs=pp.modQind$layout.obs, y.obs=pp.modQind$y.obs, variables=var.modQind, file=paste(savepath,"/states_fluxQind_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=FALSE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modQind$layout.mod$layout, y.mod=pp.modQind$y.mod, layout.obs=pp.modQind$layout.obs, y.obs=pp.modQind$y.obs, variables=var.modQind, file=paste(savepath,"/states_fluxQind_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=FALSE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modQind$layout.mod$layout, y.mod=pp.modQind$y.mod, layout.obs=pp.modQind$layout.obs, y.obs=pp.modQind$y.obs, variables=var.modQind, file=paste(savepath,"/states_fluxQind_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=FALSE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modQind$layout.mod$layout, y.mod=pp.modQind$y.mod, layout.obs=pp.modQind$layout.obs, y.obs=pp.modQind$y.obs, variables=var.modQind, file=paste(savepath,"/states_fluxQind_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=FALSE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modQind$layout.mod$layout, y.mod=pp.modQind$y.mod, layout.obs=pp.modQind$layout.obs, y.obs=pp.modQind$y.obs, variables=var.modQind, file=paste(savepath,"/states_fluxQind_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=FALSE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modQind$layout.mod$layout, y.mod=pp.modQind$y.mod, layout.obs=pp.modQind$layout.obs, y.obs=pp.modQind$y.obs, variables=var.modQind, file=paste(savepath,"/states_fluxQind_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=FALSE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modQsur$layout.mod$layout, y.mod=pp.modQsur$y.mod, layout.obs=pp.modQsur$layout.obs, y.obs=pp.modQsur$y.obs, variables=var.modQsur, file=paste(savepath,"/states_fluxQsur_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=FALSE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modAind$layout.mod$layout, y.mod=pp.modAind$y.mod, layout.obs=pp.modAind$layout.obs, y.obs=pp.modAind$y.obs, variables=var.modAind, file=paste(savepath,"/states_fluxAind_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modAind$layout.mod$layout, y.mod=pp.modAind$y.mod, layout.obs=pp.modAind$layout.obs, y.obs=pp.modAind$y.obs, variables=var.modAind, file=paste(savepath,"/states_fluxAind_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modAind$layout.mod$layout, y.mod=pp.modAind$y.mod, layout.obs=pp.modAind$layout.obs, y.obs=pp.modAind$y.obs, variables=var.modAind, file=paste(savepath,"/states_fluxAind_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modAind$layout.mod$layout, y.mod=pp.modAind$y.mod, layout.obs=pp.modAind$layout.obs, y.obs=pp.modAind$y.obs, variables=var.modAind, file=paste(savepath,"/states_fluxAind_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modAind$layout.mod$layout, y.mod=pp.modAind$y.mod, layout.obs=pp.modAind$layout.obs, y.obs=pp.modAind$y.obs, variables=var.modAind, file=paste(savepath,"/states_fluxAind_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAind$layout.mod$layout, y.mod=pp.modAind$y.mod, layout.obs=pp.modAind$layout.obs, y.obs=pp.modAind$y.obs, variables=var.modAind, file=paste(savepath,"/states_fluxAind_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modAsur$layout.mod$layout, y.mod=pp.modAsur$y.mod, layout.obs=pp.modAsur$layout.obs, y.obs=pp.modAsur$y.obs, variables=var.modAsur, file=paste(savepath,"/states_fluxAsur_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modTind$layout.mod$layout, y.mod=pp.modTind$y.mod, layout.obs=pp.modTind$layout.obs, y.obs=pp.modTind$y.obs, variables=var.modTind, file=paste(savepath,"/states_fluxTind_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modTsur$layout.mod$layout, y.mod=pp.modTsur$y.mod, layout.obs=pp.modTsur$layout.obs, y.obs=pp.modTsur$y.obs, variables=var.modTsur, file=paste(savepath,"/states_fluxTsur_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.mod3$layout.mod$layout, y.mod=pp.mod3$y.mod, layout.obs=pp.mod3$layout.obs, y.obs=pp.mod3$y.obs, variables=var.mod3, file=paste(savepath,"/states_fluxA_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, write.load=TRUE, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod3$layout.mod$layout, y.mod=pp.mod3$y.mod, layout.obs=pp.mod3$layout.obs, y.obs=pp.mod3$y.obs, variables=var.mod3, file=paste(savepath,"/states_fluxA_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod3$layout.mod$layout, y.mod=pp.mod3$y.mod, layout.obs=pp.mod3$layout.obs, y.obs=pp.mod3$y.obs, variables=var.mod3, file=paste(savepath,"/states_fluxA_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod3$layout.mod$layout, y.mod=pp.mod3$y.mod, layout.obs=pp.mod3$layout.obs, y.obs=pp.mod3$y.obs, variables=var.mod3, file=paste(savepath,"/states_fluxA_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod3$layout.mod$layout, y.mod=pp.mod3$y.mod, layout.obs=pp.mod3$layout.obs, y.obs=pp.mod3$y.obs, variables=var.mod3, file=paste(savepath,"/states_fluxA_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod3$layout.mod$layout, y.mod=pp.mod3$y.mod, layout.obs=pp.mod3$layout.obs, y.obs=pp.mod3$y.obs, variables=var.mod3, file=paste(savepath,"/states_fluxA_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod4$layout.mod$layout, y.mod=pp.mod4$y.mod, layout.obs=pp.mod4$layout.obs, y.obs=pp.mod4$y.obs, variables=var.mod4, file=paste(savepath,"/states_fluxT_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, write.load=TRUE, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod4$layout.mod$layout, y.mod=pp.mod4$y.mod, layout.obs=pp.mod4$layout.obs, y.obs=pp.mod4$y.obs, variables=var.mod4, file=paste(savepath,"/states_fluxT_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod4$layout.mod$layout, y.mod=pp.mod4$y.mod, layout.obs=pp.mod4$layout.obs, y.obs=pp.mod4$y.obs, variables=var.mod4, file=paste(savepath,"/states_fluxT_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod4$layout.mod$layout, y.mod=pp.mod4$y.mod, layout.obs=pp.mod4$layout.obs, y.obs=pp.mod4$y.obs, variables=var.mod4, file=paste(savepath,"/states_fluxT_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod4$layout.mod$layout, y.mod=pp.mod4$y.mod, layout.obs=pp.mod4$layout.obs, y.obs=pp.mod4$y.obs, variables=var.mod4, file=paste(savepath,"/states_fluxT_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, timestep.fac=su$layout$timestep.fac, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod5$layout.mod$layout, y.mod=pp.mod5$y.mod, layout.obs=pp.mod5$layout.obs, y.obs=pp.mod5$y.obs, variables=var.mod5, file=paste(savepath,"/states_fluxQ_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod5$layout.mod$layout, y.mod=pp.mod5$y.mod, layout.obs=pp.mod5$layout.obs, y.obs=pp.mod5$y.obs, variables=var.mod5, file=paste(savepath,"/states_fluxQ_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod5$layout.mod$layout, y.mod=pp.mod5$y.mod, layout.obs=pp.mod5$layout.obs, y.obs=pp.mod5$y.obs, variables=var.mod5, file=paste(savepath,"/states_fluxQ_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod5$layout.mod$layout, y.mod=pp.mod5$y.mod, layout.obs=pp.mod5$layout.obs, y.obs=pp.mod5$y.obs, variables=var.mod5, file=paste(savepath,"/states_fluxQ_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod5$layout.mod$layout, y.mod=pp.mod5$y.mod, layout.obs=pp.mod5$layout.obs, y.obs=pp.mod5$y.obs, variables=var.mod5, file=paste(savepath,"/states_fluxQ_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod5$layout.mod$layout, y.mod=pp.mod5$y.mod, layout.obs=pp.mod5$layout.obs, y.obs=pp.mod5$y.obs, variables=var.mod5, file=paste(savepath,"/states_fluxQ_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod5$layout.mod$layout, y.mod=pp.mod5$y.mod, layout.obs=pp.mod5$layout.obs, y.obs=pp.mod5$y.obs, variables=var.mod5, file=paste(savepath,"/states_fluxQ_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=FALSE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod1$layout.mod$layout, y.mod=pp.mod1$y.mod, layout.obs=pp.mod1$layout.obs, y.obs=pp.mod1$y.obs, variables=var.mod1, file=paste(savepath,"/states_mod1_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod1$layout.mod$layout, y.mod=pp.mod1$y.mod, layout.obs=pp.mod1$layout.obs, y.obs=pp.mod1$y.obs, variables=var.mod1, file=paste(savepath,"/states_mod1_Ebeg_", tag, ".pdf", sep=""), xlim=Ebeg, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod1$layout.mod$layout, y.mod=pp.mod1$y.mod, layout.obs=pp.mod1$layout.obs, y.obs=pp.mod1$y.obs, variables=var.mod1, file=paste(savepath,"/states_mod1_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod1$layout.mod$layout, y.mod=pp.mod1$y.mod, layout.obs=pp.mod1$layout.obs, y.obs=pp.mod1$y.obs, variables=var.mod1, file=paste(savepath,"/states_mod1_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod1$layout.mod$layout, y.mod=pp.mod1$y.mod, layout.obs=pp.mod1$layout.obs, y.obs=pp.mod1$y.obs, variables=var.mod1, file=paste(savepath,"/states_mod1_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod1$layout.mod$layout, y.mod=pp.mod1$y.mod, layout.obs=pp.mod1$layout.obs, y.obs=pp.mod1$y.obs, variables=var.mod1, file=paste(savepath,"/states_mod1_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod1$layout.mod$layout, y.mod=pp.mod1$y.mod, layout.obs=pp.mod1$layout.obs, y.obs=pp.mod1$y.obs, variables=var.mod1, file=paste(savepath,"/states_mod1_E25_", tag, ".pdf", sep=""), xlim=E25, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod1$layout.mod$layout, y.mod=pp.mod1$y.mod, layout.obs=pp.mod1$layout.obs, y.obs=pp.mod1$y.obs, variables=var.mod1, file=paste(savepath,"/states_mod1_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod1$layout.mod$layout, y.mod=pp.mod1$y.mod, layout.obs=pp.mod1$layout.obs, y.obs=pp.mod1$y.obs, variables=var.mod1, file=paste(savepath,"/states_mod1_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod2$layout.mod$layout, y.mod=pp.mod2$y.mod, layout.obs=pp.mod2$layout.obs, y.obs=pp.mod2$y.obs, variables=var.mod2, file=paste(savepath,"/states_mod2_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod2$layout.mod$layout, y.mod=pp.mod2$y.mod, layout.obs=pp.mod2$layout.obs, y.obs=pp.mod2$y.obs, variables=var.mod2, file=paste(savepath,"/states_mod2_Ebeg_", tag, ".pdf", sep=""), xlim=Ebeg, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod2$layout.mod$layout, y.mod=pp.mod2$y.mod, layout.obs=pp.mod2$layout.obs, y.obs=pp.mod2$y.obs, variables=var.mod2, file=paste(savepath,"/states_mod2_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod2$layout.mod$layout, y.mod=pp.mod2$y.mod, layout.obs=pp.mod2$layout.obs, y.obs=pp.mod2$y.obs, variables=var.mod2, file=paste(savepath,"/states_mod2_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod2$layout.mod$layout, y.mod=pp.mod2$y.mod, layout.obs=pp.mod2$layout.obs, y.obs=pp.mod2$y.obs, variables=var.mod2, file=paste(savepath,"/states_mod2_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        ##plot.results(layout.mod=pp.mod25$layout.mod$layout, y.mod=pp.mod25$y.mod, layout.obs=pp.mod25$layout.obs, y.obs=pp.mod25$y.obs, variables=var.mod25, file=paste(savepath,"/states_mod3_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod2$layout.mod$layout, y.mod=pp.mod2$y.mod, layout.obs=pp.mod2$layout.obs, y.obs=pp.mod2$y.obs, variables=var.mod2, file=paste(savepath,"/states_mod2_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod2$layout.mod$layout, y.mod=pp.mod2$y.mod, layout.obs=pp.mod2$layout.obs, y.obs=pp.mod2$y.obs, variables=var.mod2, file=paste(savepath,"/states_mod2_E25_", tag, ".pdf", sep=""), xlim=E25, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod2$layout.mod$layout, y.mod=pp.mod2$y.mod, layout.obs=pp.mod2$layout.obs, y.obs=pp.mod2$y.obs, variables=var.mod2, file=paste(savepath,"/states_mod2_E3_", tag, ".pdf", sep=""), xlim=E3, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.mod2$layout.mod$layout, y.mod=pp.mod2$y.mod, layout.obs=pp.mod2$layout.obs, y.obs=pp.mod2$y.obs, variables=var.mod2, file=paste(savepath,"/states_mod2_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcCon$layout.mod$layout, y.mod=pp.modTcCon$y.mod, layout.obs=pp.modTcCon$layout.obs, y.obs=pp.modTcCon$y.obs, variables=var.modTcCon, file=paste(savepath,"/states_modTcCon_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcCon$layout.mod$layout, y.mod=pp.modTcCon$y.mod, layout.obs=pp.modTcCon$layout.obs, y.obs=pp.modTcCon$y.obs, variables=var.modTcCon, file=paste(savepath,"/states_modTcCon_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcCon$layout.mod$layout, y.mod=pp.modTcCon$y.mod, layout.obs=pp.modTcCon$layout.obs, y.obs=pp.modTcCon$y.obs, variables=var.modTcCon, file=paste(savepath,"/states_modTcCon_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcCon$layout.mod$layout, y.mod=pp.modTcCon$y.mod, layout.obs=pp.modTcCon$layout.obs, y.obs=pp.modTcCon$y.obs, variables=var.modTcCon, file=paste(savepath,"/states_modTcCon_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcCon$layout.mod$layout, y.mod=pp.modTcCon$y.mod, layout.obs=pp.modTcCon$layout.obs, y.obs=pp.modTcCon$y.obs, variables=var.modTcCon, file=paste(savepath,"/states_modTcCon_E25_", tag, ".pdf", sep=""), xlim=E25, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAcCon$layout.mod$layout, y.mod=pp.modAcCon$y.mod, layout.obs=pp.modAcCon$layout.obs, y.obs=pp.modAcCon$y.obs, variables=var.modAcCon, file=paste(savepath,"/states_modAcCon_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAcConDrn$layout.mod$layout, y.mod=pp.modAcConDrn$y.mod, layout.obs=pp.modAcConDrn$layout.obs, y.obs=pp.modAcConDrn$y.obs, variables=var.modAcConDrn, file=paste(savepath,"/states_modAcConDrn_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAcDrn$layout.mod$layout, y.mod=pp.modAcDrn$y.mod, layout.obs=pp.modAcDrn$layout.obs, y.obs=pp.modAcDrn$y.obs, variables=var.modAcDrn, file=paste(savepath,"/states_modAcDrn_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAcDrn$layout.mod$layout, y.mod=pp.modAcDrn$y.mod, layout.obs=pp.modAcDrn$layout.obs, y.obs=pp.modAcDrn$y.obs, variables=var.modAcDrn, file=paste(savepath,"/states_modAcDrn_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAcDrn$layout.mod$layout, y.mod=pp.modAcDrn$y.mod, layout.obs=pp.modAcDrn$layout.obs, y.obs=pp.modAcDrn$y.obs, variables=var.modAcDrn, file=paste(savepath,"/states_modAcDrn_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAcDrn$layout.mod$layout, y.mod=pp.modAcDrn$y.mod, layout.obs=pp.modAcDrn$layout.obs, y.obs=pp.modAcDrn$y.obs, variables=var.modAcDrn, file=paste(savepath,"/states_modAcDrn_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAcDrn$layout.mod$layout, y.mod=pp.modAcDrn$y.mod, layout.obs=pp.modAcDrn$layout.obs, y.obs=pp.modAcDrn$y.obs, variables=var.modAcDrn, file=paste(savepath,"/states_modAcDrn_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAcDrn$layout.mod$layout, y.mod=pp.modAcDrn$y.mod, layout.obs=pp.modAcDrn$layout.obs, y.obs=pp.modAcDrn$y.obs, variables=var.modAcDrn, file=paste(savepath,"/states_modAcDrn_E25_", tag, ".pdf", sep=""), xlim=E25, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modAcDrn$layout.mod$layout, y.mod=pp.modAcDrn$y.mod, layout.obs=pp.modAcDrn$layout.obs, y.obs=pp.modAcDrn$y.obs, variables=var.modAcDrn, file=paste(savepath,"/states_modAcDrn_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcConDrn$layout.mod$layout, y.mod=pp.modTcConDrn$y.mod, layout.obs=pp.modTcConDrn$layout.obs, y.obs=pp.modTcConDrn$y.obs, variables=var.modTcConDrn, file=paste(savepath,"/states_modTcConDrn_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcDrn$layout.mod$layout, y.mod=pp.modTcDrn$y.mod, layout.obs=pp.modTcDrn$layout.obs, y.obs=pp.modTcDrn$y.obs, variables=var.modTcDrn, file=paste(savepath,"/states_modTcDrn_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcDrn$layout.mod$layout, y.mod=pp.modTcDrn$y.mod, layout.obs=pp.modTcDrn$layout.obs, y.obs=pp.modTcDrn$y.obs, variables=var.modTcDrn, file=paste(savepath,"/states_modTcDrn_E00_", tag, ".pdf", sep=""), xlim=E00, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcDrn$layout.mod$layout, y.mod=pp.modTcDrn$y.mod, layout.obs=pp.modTcDrn$layout.obs, y.obs=pp.modTcDrn$y.obs, variables=var.modTcDrn, file=paste(savepath,"/states_modTcDrn_E0_", tag, ".pdf", sep=""), xlim=E0, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcDrn$layout.mod$layout, y.mod=pp.modTcDrn$y.mod, layout.obs=pp.modTcDrn$layout.obs, y.obs=pp.modTcDrn$y.obs, variables=var.modTcDrn, file=paste(savepath,"/states_modTcDrn_E1_", tag, ".pdf", sep=""), xlim=E1, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcDrn$layout.mod$layout, y.mod=pp.modTcDrn$y.mod, layout.obs=pp.modTcDrn$layout.obs, y.obs=pp.modTcDrn$y.obs, variables=var.modTcDrn, file=paste(savepath,"/states_modTcDrn_E2_", tag, ".pdf", sep=""), xlim=E2, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcDrn$layout.mod$layout, y.mod=pp.modTcDrn$y.mod, layout.obs=pp.modTcDrn$layout.obs, y.obs=pp.modTcDrn$y.obs, variables=var.modTcDrn, file=paste(savepath,"/states_modTcDrn_E25_", tag, ".pdf", sep=""), xlim=E25, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcDrn$layout.mod$layout, y.mod=pp.modTcDrn$y.mod, layout.obs=pp.modTcDrn$layout.obs, y.obs=pp.modTcDrn$y.obs, variables=var.modTcDrn, file=paste(savepath,"/states_modTcDrn_E4_", tag, ".pdf", sep=""), xlim=E4, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcGwt$layout.mod$layout, y.mod=pp.modTcGwt$y.mod, layout.obs=pp.modTcGwt$layout.obs, y.obs=pp.modTcGwt$y.obs, variables=var.modTcGwt, file=paste(savepath,"/states_modTcGwt_Eapp_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
		plot.results(layout.mod=pp.modTcGwt$layout.mod$layout, y.mod=pp.modTcGwt$y.mod, layout.obs=pp.modTcGwt$layout.obs, y.obs=pp.modTcGwt$y.obs, variables=var.modTcGwt, file=paste(savepath,"/states_modTcGwt_Ebeg_", tag, ".pdf", sep=""), xlim=Ebeg, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    }else if(grepl("Str02", sett)){
        plot.results(layout.mod=pp.mod$layout.mod$layout, y.mod=pp.mod$y.mod, layout.obs=pp.mod$layout.obs, y.obs=pp.mod$y.obs, variables=var.mod, file=paste(savepath,"/states_mod_", tag, ".pdf", sep=""), xlim=Ebeg, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
        plot.results(layout.mod=pp.modT$layout.mod$layout, y.mod=pp.modT$y.mod, layout.obs=pp.modT$layout.obs, y.obs=pp.modT$y.obs, variables=var.modT, file=paste(savepath,"/states_modT_", tag, ".pdf", sep=""), xlim=Eapp, per.area=TRUE, hru.areas=hru.areas, translate.var=translate.var, translate.to=translate.to)
    }
}

if(plot.uncert){
    E00 <- as.POSIXct(c("2009-05-19 00:00", "2009-05-20 12:00"))
    E0 <- as.POSIXct(c("2009-05-21 12:00", "2009-05-23 00:00"))
	#E0 <- as.POSIXct(c("2009-05-21 12:00", "2009-05-22 12:00"))
    E0zoom <- as.POSIXct(c("2009-05-21 22:30", "2009-05-22 03:00"))
	E1 <- as.POSIXct(c("2009-05-26 12:00", "2009-05-27 12:00"))
    E2 <- as.POSIXct(c("2009-06-06 06:00", "2009-06-07 18:00"))
    E25 <- as.POSIXct(c("2009-06-19 00:00", "2009-06-21 00:00"))
    E3 <- as.POSIXct(c("2009-06-25 08:00", "2009-06-29 00:00"))
    E35<- as.POSIXct(c("2009-07-05 00:00", "2009-07-10 00:00"))
	E4 <- as.POSIXct(c("2009-07-10 00:00", "2009-07-25 18:00"))
    Ebeg <- as.POSIXct(c("2009-03-01 00:00", "2009-05-15 00:00"))
    Eapp <- as.POSIXct(c("2009-05-15 00:00", "2009-07-30 00:00"))
    Ewint <- as.POSIXct(c("2008-10-01 00:00", "2009-03-31 00:00"))
    Eapp.tobi <- as.POSIXct(c("2009-05-20 00:00", "2009-07-23 00:00"))
	
	#load(paste0(savepath,"/constr_layout.RData"))
	#load(paste0(savepath,"/constr_observations.RData"))
	#su$layout <- lyt ## ATTENTION: the layout is overwritten to include terbuthylazine after event e2
	#su$observations <- obsr	
    
	list.su <- list("MexpH4"=su)
	dir.create(paste(savepath,"/predictions",sep=""))
	tag.save <- tag
	tag <- paste0(tag,ifelse(place.terb.early,"","_terb_late"))
	pdf(paste(savepath, "/predictions/", tag, "_E2E3E4E5.pdf", sep=""), width=10, height=7)
	#fil <- paste(savepath, "/predictions/", tag, "_E1E2E3E4.pdf", sep="")
    plot.predictions(list.su, probs=c(0.05,0.95), xlim=list(E2=E1,E3=E2,E4=E3,E5=E4,Total=Eapp.tobi), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_E00.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), n.samp=1, xlim=list(E00=E00), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_E0.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), n.samp=1, xlim=list(E0=E0), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_E1.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), n.samp=1, xlim=list(E1=E1), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_E2.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), n.samp=1, xlim=list(E2=E2), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_E25.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), n.samp=1, xlim=list(E25=E25), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_E3.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), n.samp=1, xlim=list(E3=E3), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_E4.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), n.samp=1, xlim=list(E4=E4), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_Eapp.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), xlim=list(Eapp=Eapp), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_Eapp_tobi.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), xlim=list(Eapp=Eapp.tobi), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	pdf(paste(savepath, "/predictions/", tag, "_Ebeg.pdf", sep=""), width=10, height=7)
    plot.predictions(list.su, probs=c(0.05,0.95), xlim=list(Eapp=Ebeg), tme.orig=su$layout$tme.orig, metrics=FALSE, plot.var=c("C1Wv_Qstream","C1Tc1_Qstream","C1Tc2_Qstream"))
    dev.off()
	tag <- tag.save
	##plot.predictions(list.su, probs=c(0.05,0.95), xlim="pred", lp.num.pred=dat$lp.num.pred, tme.orig=su$layout$tme.orig, metrics=TRUE)
    ## plot.predictions(list.su, probs=c(0.05,0.95), xlim=xlim, lp.num.pred=dat$lp.num.pred, tme.orig=tme.orig)
    ## plot.predictions(list.su, probs=c(0.05,0.95), xlim=xlim, ylim=c(0,1))
    ## plot.predictions(list.su, probs=c(0.05,0.95), xlim=xlim, ylim=c(0,1), n.samp=1, arrange=c(E1=1,E2=2,E3a=3,E3ab=4, E4a=5))
    #dev.off()
}

if(sensitivity.analysis){
	##load("../output/A1Str07h2b/v05_spraydrift_QE1c07i_TTE1c07i_priorC2/model.parameters.RData")
	##load("../output/A1Str07h2b/v05_spraydrift_QE1c07i_TTE1c07i_priorC2/input.RData")
	##su$model$parameters <- model.parameters
	##su$parameter.sample <- NULL
	##su$input <- input
    ##tag <- paste0(tag, "_T1")
	##sensitivity_analysis(sudriv=su, vary.all=list("Glo%Cmlt_E"=log(c(1.0,1.5)),
													##"Glo%Cmlt_K_Qq_FR"=log(c(0.005,0.2))))
													##"BetaQq_UR"=log(c(1.5,4))))
													##"U[^1]T[1-2]%SlOne_FR"=log(c(3,15)),
													##"U[^1]T[1-2]%SlTwo_FR"=log(c(3,15)),
													##"GloTr%CmltSlOne_RR"=log(c(3,15)),
													##"GloTr%CmltSlTwo_RR"=log(c(3,15))))
													##"U[^1]W%nRes_FR"=log(c(1.9,2.9))))
													##"K_Qb_FR"=c(0.001,0.01)))
													##"Alfa_Qq_FR"=log(c(1,2.5)),
													##"Glo%Cmlt_Dspl_SD"=c(0.6,0.9),
													##"Glo%Cmlt_Pmax_ED"=log(c(9,13)),
													##"Glo%CmltSmax_IR"=log(c(4,14))))
													##"Glo%CmltSmax_UR"=log(c(150,250)),
													##"Glo%Cmlt_K_Qq_RR"=log(c(0.8,1.8)),
													##"Glo%Cmlt_K_Qq_SR"=log(c(0.00001, 0.001)),
													##"GloTr%CmltSlOne_IR"=log(c(3,15)),
													##"GloTr%CmltKd_WR"=log(c(0.0006,0.0015)),
													##"GloTr%CmltSlTwo_IR"=log(c(10,60)),
													##"GloTr%CmltRs_WR"=log(c(0.001,0.02)),
													##"U1W%KpQq_FR"=log(c(0.05,0.2))))
}

if(reproduce.synthetic){


    su2 <- su
    ##su2$layout$pred <- su$layout$calib
    ##su2$likelihood$f.sample <- LogLikelihoodHydrology_la9_skewt_sample_synth
    ##load("realiz_eta_correlated_OU.RData")
    a <- sampling_wrapper(su2, sample.par=FALSE, n.sample=1, biased=grepl("la9", tag) & !grepl("mu0", tag), sample.likeli=TRUE)
    su2$observations[su$layout$calib] <- c(a)[su$layout$calib]
    ## try to infer the parameters of the synthetic realization
    n.dim <- sum(c(su2$model$par.fit, su2$likelihood$par.fit))
    library(MCMCEnsembleSampler)
    log.posterior <- logposterior
    max.iter <- 10000
    n.walkers <- 100
    prs <- c(su2$model$parameters[as.logical(su2$model$par.fit)], su2$likelihood$parameters[as.logical(su2$likelihood$par.fit)])
    init.range <- cbind(prs/2, prs*2)
    init.range[prs<0,] <- init.range[prs<0,c(2,1)]

    s <- remove.chains(su2, brn.in=0, logpost=NA)$sample
    init.range <- t(apply(s[nrow(s),,], 1, quantile, probs=c(0.1,0.9)))
    result.s.m = s.m.mcmc(log.posterior, max.iter=max.iter, n.walkers=n.walkers, n.dim=n.dim, init.range=init.range, sudriv=su2, prior=TRUE, apprx=TRUE)
    s <- result.s.m$samples
    su2$parameter.sample <- aperm(s, perm=c(2,3,1))
    colnames(su2$parameter.sample) <- c(names(su2$model$parameters)[as.logical(su2$model$par.fit)], names(su2$likelihood$parameters)[as.logical(su2$likelihood$par.fit)])##, paste("mu", 1:(1*n.mu), sep=""))
    su2$posterior.sample <- t(result.s.m$log.p)
    save(su2, file=paste("sudriv_output/", su2$settings$subcatchment, "/eval_lik/", tag, "/repr_synth", tag, ".RData", sep=""))

    ## plot the histograms of the sample, together with the true value (used to generate the time series)
    prs <- ifelse(as.logical(c(su2$model$args$parTran, su2$likelihood$tran)), exp(c(su2$model$parameters, su2$likelihood$parameters)), c(su2$model$parameters, su2$likelihood$parameters))
    v.line <- prs[as.logical(c(su2$model$par.fit, su2$likelihood$par.fit))]
    names(v.line) <- colnames(su2$parameter.sample)
    thin <- 1
    lower.logpost <- NA
    res <- 300
    height <- 7 #inches
    width <- 9 #inches
    png(paste("sudriv_output/", su$settings$subcatchment, "/eval_lik/", tag, "/reprhist5", tag, ".png", sep=""), width=width*res, height=height*res, res=res)
    plot.markov.hist(sudriv=su2, brn.in=0, v.line=v.line)
    dev.off()
    png(paste("sudriv_output/", su$settings$subcatchment, "/eval_lik/", tag, "/reprchain5", tag, ".png", sep=""), width=width*res, height=height*res, res=res)
    plot.markov.chain(sudriv=su2, brn.in=0, thin=thin)
    dev.off()
    png(paste("sudriv_output/", su$settings$subcatchment, "/eval_lik/", tag, "/reprcor5", tag, ".png", sep=""), width=width*res, height=height*res, res=res)
    plot.cor(sudriv=su2, brn.in=0, thin=thin)
    dev.off()

    ## calculate confidence intervals and draw boxplot
    brn.in <- 100
    samp.brn <- su2$parameter.sample[brn.in:nrow(su2$parameter.sample),,]
    mm <- which(su2$posterior.sample[brn.in:nrow(su2$posterior.sample),]==max(su2$posterior.sample[brn.in:nrow(su2$posterior.sample),]), arr.ind =TRUE)[1,]
    v.line.hours <- v.line
    collapsed <- matrix(aperm(samp.brn, c(1,3,2)), nrow=prod(dim(samp.brn)[c(1,3)]), ncol=dim(samp.brn)[2])
    for(i in 1:ncol(collapsed)){
        if(as.logical(c(su2$model$args$parTran, su2$likelihood$tran)[as.logical(c(su2$model$par.fit, su2$likelihood$par.fit))])[i]){
            collapsed[,i] <- exp(collapsed[,i])
        }
        if(i %in% c(3,4)){collapsed[,i] <- collapsed[,i]/24; v.line.hours[i] <- v.line.hours[i]/24}#transformation to hours
        if(i == 7){collapsed[,i] <- collapsed[,i]*24; v.line.hours[i] <- v.line.hours[i]*24}
    }
    quants <- apply(collapsed, 2, quantile, c(0.025,0.975))
    means <- apply(collapsed, 2, mean)
    mm <- matrix(su2$posterior.sample[brn.in:nrow(su2$posterior.sample),], nrow=prod(dim(samp.brn)[c(1,3)]))
    ml <- collapsed[which(c(mm)==max(mm))[1],]
    colnames(quants) <- colnames(su2$parameter.sample)

    ## save the sythetic parameters (true parameters)
    synth.par <- c(su$model$parameters, su$likelihood$parameters)
    save(synth.par, file=paste("sudriv_output/", su2$settings$subcatchment, "/eval_lik/", tag, "/synth_par", tag, ".RData", sep=""))

    ## make prediction
    su2 <- select.maxlikpars(su2)
    m <- which(su2$posterior.sample == max(su2$posterior.sample))[1]
    su2$predicted$sample <- sampling_wrapper(su2, brn.in=75, sample.par=TRUE, n.sample=500, biased=grepl("la9", tag) & !grepl("mu0", tag), sample.likeli=TRUE, rep.mu.times=rep.mu.times)
    su2$predicted$det    <- sampling_wrapper(su2, sample.par=FALSE, n.sample=1, biased=grepl("la9", tag) & !grepl("mu0", tag), sample.likeli=FALSE, rep.mu.times=rep.mu.times)
    save(su2, file=paste("sudriv_output/", su2$settings$subcatchment, "/eval_lik/", tag, "/repr_synth", tag, ".RData", sep=""))

    ## plot stats in prediction period
    ##mus <- c(su2$parameter.sample[m,(sum(c(su2$model$par.fit, su2$likelihood$par.fit))+1):(dim(su2$parameter.sample)[2])])
    mus <- NA
    mu=mus
    su2$layout$pred <- su$layout$calib
    su2$layout$calib <- numeric(0)
    dat <- pred.stats(su2, mu=mus, biased=is.numeric(mus[1]), rep.mu.times=rep.mu.times)

    ## plot predictions
    pdf(paste("sudriv_output/", su2$settings$subcatchment, "/eval_lik/", tag, "/synth_predictions", tag, ".pdf", sep=""), width=12)
    plot.predictions(su2, xlim=c(24000,25000), probs=c(0.05,0.95), lp.num.pred=dat$lp.num.pred)
    plot.predictions(su2, xlim=c(20000,26000), probs=c(0.05,0.95), lp.num.pred=dat$lp.num.pred)
    ## plot.predictions(su2, xlim=c(35000,40000), probs=c(0.05, 0.95), lp.num.pred=dat$lp.num.pred)
    ## plot.predictions(su2, xlim=c(35000,40000), probs=c(0.05, 0.95), lp.num.pred=dat$lp.num.pred)
    plot.predictions(su2, xlim=c(24000,25000), n.samp=2)
    plot.predictions(su2, xlim=c(20000,26000), n.samp=2)
    ## plot.predictions(su2, xlim=c(35000,40000), n.samp=2)
    ## plot.predictions(su2, xlim=c(35000,40000), n.samp=2)
    dev.off()

    pdf(paste("sudriv_output/", su2$settings$subcatchment, "/eval_lik/", tag, "/stats", tag, ".pdf", sep=""), width=12)
    plot.ts.quantiles(dat, su2, xlim=c(24000,25000))
    ## plot.ts.quantiles(dat, su2, xlim=c(35000,40000))
    plot.dens.quantiles(dat, su2, ind.sel=su2$layout$pred, distr="unif")
    plot.dens.quantiles(dat, su2, ind.sel=su2$layout$pred)
    plot.dens.innovation(dat, su2, ind.sel=su2$layout$pred)
    plot.Qdet.quantiles(dat, su2, ind.sel=su2$layout$pred)
    plot.pacf.quantiles(dat, su2, ind.sel=su2$layout$pred, lag.max=10)
    ## plot.dens.quantiles(dat, su2, ind.sel=su2$layout$calib[!auto2.hyp])
    ## plot.dens.quantiles(dat, su2, ind.sel=su2$layout$calib[auto2.hyp])
    ## plot.Qdet.quantiles(dat, su2, ind.sel=su2$layout$calib[!auto2.hyp])
    ## plot.Qdet.quantiles(dat, su2, ind.sel=su2$layout$calib[auto2.hyp])
    ## plot.pacf.quantiles(dat, su2, ind.sel=su2$layout$calib[!auto2.hyp], lag.max=10)
    ## plot.pacf.quantiles(dat, su2, ind.sel=su2$layout$calib[auto2.hyp], lag.max=10)
    ##plot.dens.quantiles(dat, su2, ind.sel=su2$layout$calib[auto2.hyp])
    ##plot.Qdet.quantiles(dat, su2, ind.sel=su2$layout$calib[!auto2.hyp])
    dev.off()
    pdf(paste("sudriv_output/", su2$settings$subcatchment, "/eval_lik/", tag, "/stats_calib", tag, ".pdf", sep=""), width=12)
    plot.ts.quantiles(dat, su2, xlim=c(24000,25000))
    ## plot.ts.quantiles(dat, su2, xlim=c(35000,40000))
    plot.dens.quantiles(dat, su2, ind.sel=su2$layout$calib, distr="unif")
    plot.dens.quantiles(dat, su2, ind.sel=su2$layout$calib)
    plot.dens.innovation(dat, su2, ind.sel=su2$layout$calib)
    plot.Qdet.quantiles(dat, su2, ind.sel=su2$layout$calib)
    plot.pacf.quantiles(dat, su2, ind.sel=su2$layout$calib, lag.max=10)
    ## plot.dens.quantiles(dat, su2, ind.sel=su2$layout$calib[!auto2.hyp])
    ## plot.dens.quantiles(dat, su2, ind.sel=su2$layout$calib[auto2.hyp])
    ## plot.Qdet.quantiles(dat, su2, ind.sel=su2$layout$calib[!auto2.hyp])
    ## plot.Qdet.quantiles(dat, su2, ind.sel=su2$layout$calib[auto2.hyp])
    ## plot.pacf.quantiles(dat, su2, ind.sel=su2$layout$calib[!auto2.hyp], lag.max=10)
    ## plot.pacf.quantiles(dat, su2, ind.sel=su2$layout$calib[auto2.hyp], lag.max=10)
    ##plot.dens.quantiles(dat, su2, ind.sel=su2$layout$calib[auto2.hyp])
    ##plot.Qdet.quantiles(dat, su2, ind.sel=su2$layout$calib[!auto2.hyp])
    dev.off()
}
if(!is.na(compare.models[1])){
	# load sudriv objects
	list.su <- list()
	for(mod.curr in names(compare.models)){
		su <- sudriv.setup(settings=compare.models[mod.curr])
		savepath <- paste0("../output/", ifelse(reload.timedep,"timedeppar/",""), su$settings$subcatchment, su$settings$structure, "/", tag.mult[mod.curr])
		load(paste0(savepath,"/su_",tag.mult[mod.curr],".RData"))
		# if(mod.curr==(names(compare.models)[1])){
			# load(paste0(savepath,"/constr_layout.RData"))
			# load(paste0(savepath,"/constr_observations.RData"))
			# su$layout <- lyt ## ATTENTION: the layout is overwritten to include terbuthylazine after event e2
			# su$observations <- obsr
		# }
		list.su[[mod.curr]] <- su
	}
	names(list.su) <- names(compare.models)
    E00 <- as.POSIXct(c("2009-05-19 00:00", "2009-05-20 12:00"))
    E0 <- as.POSIXct(c("2009-05-21 12:00", "2009-05-23 00:00"))
	#E0 <- as.POSIXct(c("2009-05-21 12:00", "2009-05-22 12:00"))
    E0zoom <- as.POSIXct(c("2009-05-21 22:30", "2009-05-22 03:00"))
	E1 <- as.POSIXct(c("2009-05-26 12:00", "2009-05-27 12:00"))
    E2 <- as.POSIXct(c("2009-06-06 06:00", "2009-06-07 18:00"))
    E25 <- as.POSIXct(c("2009-06-19 00:00", "2009-06-21 00:00"))
    E3 <- as.POSIXct(c("2009-06-25 08:00", "2009-06-29 00:00"))
    E35<- as.POSIXct(c("2009-07-05 00:00", "2009-07-10 00:00"))
	E4 <- as.POSIXct(c("2009-07-10 00:00", "2009-07-25 18:00"))
    Ebeg <- as.POSIXct(c("2009-03-01 00:00", "2009-05-15 00:00"))
    Eapp <- as.POSIXct(c("2009-05-15 00:00", "2009-07-30 00:00"))
    Ewint <- as.POSIXct(c("2008-10-01 00:00", "2009-03-31 00:00"))
    Eapp.tobi <- as.POSIXct(c("2009-05-20 00:00", "2009-07-23 00:00"))
	if(export){
		tmp <- read.app.hru.areas(compare.models=compare.models, tag.mult=tag.mult)
		print("tmp: ")
		print(tmp)
		fil <- paste0("../output/",ifelse(reload.timedep,"timedeppar/",""),"modelcomparison/",paste0(names(compare.models),collapse=""),"stoch_E1E2E3E4.pdf")
		plot.predictions(list.su=list.su, probs=c(0.05,0.95), xlim=list(E1=E0,E2=E1,E3=E2,E4=E3,E5=E4,Total=Eapp.tobi), tme.orig=tme.orig,arrange=arrange,plot.var=plot.var, scl=1, alp=0.8, loads.det=tmp$loads.det, app.hru.areas=tmp$app.hru.areas, file=fil)
		#dev.off()
	}
	pdf(paste0("../output/",ifelse(reload.timedep,"timedeppar/",""),"modelcomparison/",paste0(names(compare.models),collapse=""),"_E0.pdf"))
	plot.predictions(list.su=list.su, xlim=list(E1=E0), tme.orig=tme.orig,arrange=arrange,plot.var=plot.var)
	dev.off()
	pdf(paste0("../output/",ifelse(reload.timedep,"timedeppar/",""),"modelcomparison/",paste0(names(compare.models),collapse=""),"_E0zoom.pdf"))
	plot.predictions(list.su=list.su, xlim=list(E1=E0zoom), tme.orig=tme.orig,arrange=arrange,plot.var=plot.var)
	dev.off()
	#png(paste0("../output/",ifelse(reload.timedep,"timedeppar/",""),"modelcomparison/",paste0(names(compare.models),collapse=""),"_E0.png"), height=6, width=6, units="in", res=300)
	#plot.predictions(list.su=list.su, xlim=list(E1=E0), tme.orig=tme.orig,arrange=arrange,plot.var=plot.var)
	#dev.off()
	pdf(paste0("../output/",ifelse(reload.timedep,"timedeppar/",""),"modelcomparison/",paste0(names(compare.models),collapse=""),"_E1E2E4E5.pdf"), width=10, height=7)
	plot.predictions(list.su=list.su, xlim=list(E1=E0,E2=E1,E4=E3,E5=E4), tme.orig=tme.orig,arrange=arrange,plot.var=plot.var)
	dev.off()
	pdf(paste0("../output/",ifelse(reload.timedep,"timedeppar/",""),"modelcomparison/",paste0(names(compare.models),collapse=""),"_E2E3E4E5.pdf"), width=10, height=7)
	plot.predictions(list.su=list.su, xlim=list(E2=E1,E3=E2,E4=E3,E5=E4), tme.orig=tme.orig,arrange=arrange,plot.var=plot.var)
	dev.off()
	pdf(paste0("../output/",ifelse(reload.timedep,"timedeppar/",""),"modelcomparison/",paste0(names(compare.models),collapse=""),"_Eapp.pdf"), width=10, height=7)
	plot.predictions(list.su=list.su, xlim=list(Eapp=Eapp), tme.orig=tme.orig,arrange=arrange,plot.var=plot.var)
	dev.off()
	pdf(paste0("../output/",ifelse(reload.timedep,"timedeppar/",""),"modelcomparison/",paste0(names(compare.models),collapse=""),"_Ebeg.pdf"), width=10, height=7)
	plot.predictions(list.su=list.su, xlim=list(Ebeg=Ebeg), tme.orig=tme.orig,arrange=arrange,plot.var=plot.var)
	dev.off()
	##png(paste0("../output/modelcomparison/",paste0(names(compare.models),collapse=""),"_E1E2E3E4.png"))
	##plot.predictions(list.su=list.su, xlim=list(E1=E0,E2=E1,E4=E3,E5=E4), tme.orig=tme.orig,arrange=arrange,plot.var=plot.var)
	##dev.off()
}
if(analyze.likelihood){
## here we assume that a loaded su object exists
    su <- select.maxlikpars(su)
    x0 <- c(su$model$parameters[as.logical(su$model$par.fit)], su$likelihood$parameters[as.logical(su$likelihood$par.fit)])
    logpost <- logposterior(x0=x0, sudriv=su, prior=FALSE, mode=TRUE, apprx=FALSE, verbose=FALSE, auto=NA, weight.equally=FALSE)
    ll <- c(logpost$loglik[su$layout$layout$var[su$layout$calib]=="C1Wv_Qstream"])
    summary(ll)
    sort(ll)[1:20]
}
