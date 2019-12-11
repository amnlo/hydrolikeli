## This script tests the time-dependent parameter inference of Peter with my workflow
op.sys <- Sys.info()["sysname"]
if(op.sys == "Linux"){.libPaths("/local/ammannlo/personalR/library")}else{Sys.setlocale("LC_TIME", "English")}
options(repos="https://stat.ethz.ch/CRAN/")
if ( ! require(ecosim) )     { install.packages("ecosim");     library(ecosim) }
if ( ! require(rjson) )     { install.packages("rjson");     library(rjson) }
if ( ! require(dplyr) )     { install.packages("dplyr");     library(dplyr) }
if ( ! require(tidyr) )     { install.packages("tidyr");     library(tidyr) }
if ( ! require(tseries) )     { install.packages("tseries");     library(tseries) }
if ( ! require(stats) )     { install.packages("stats");     library(stats) }
if ( ! require(zoo) )     { install.packages("zoo");     library(zoo) }
if ( ! require(ggplot2) )     { install.packages("ggplot2");     library(ggplot2) }
if ( ! require(car) )     { install.packages("car");     library(car) }
if ( ! require(cowplot) )     { install.packages("cowplot");     library(cowplot) }
if ( ! require(geoR) )     { install.packages("geoR");     library(geoR) }
if ( ! require(IDPmisc) )     { install.packages("IDPmisc");     library(IDPmisc) }
if ( ! require(neuralnet) )     { install.packages("neuralnet");     library(neuralnet) }
if ( ! require(compositions) )     { install.packages("compositions");     library(compositions) }
if ( ! require(graphics) )     { install.packages("graphics");     library(graphics) }
library(timedeppar)
library(truncnorm)

dir.driver <- "../../../../SuperFlex/Driver/sudriv_package/R"
source("source.functions.r")
source.functions(dir.driver)
source("../../../../Likelihood/src/hydrolikeli_herb.r")
source("../../../../Likelihood/src/plot.markov.r")
n.iter <- 20000

test 		 <- FALSE

remove.taumax<- TRUE
fix.taumax   <- FALSE
f_mean    	 <- TRUE
infer 		 <- TRUE
restart 	 <- FALSE
adapt.intrv  <- FALSE # this should be used only together with restart=TRUE
continue 	 <- FALSE
save.su 	 <- FALSE
plot		 <- FALSE
find.pattern <- FALSE
table.logpdf <- FALSE
verbose 	 <- 1

if(remove.taumax){
		vrs <- "1"
}else if(fix.taumax){
		vrs <- "3fix"
}else{
		vrs <- "3"
}

tag.vrs <- paste0("QE",vrs)

control <- list(n.interval = 50, # what is the characteristic time scale of the variability you want to keep?
				n.timedepini= 10,
				n.constpertimedep=10,
				n.save = 100,
				thin = 50)

run.these <- 1#[-c(1,5,8,9,10,11,16)]
for(cse in run.these){
	cases <- list(cse, `1`=list(c("Glo%Cmlt_Dspl_SD","Glo%CmltSmax_UR"), "dsplsd_smaxur"),
			`2`=list(c("Glo%CmltSmax_UR"), "smaxur"),
			`3`=list(c("Glo%Cmlt_BeQq_UR"), "beqqur"),
			`4`=list(c("Glo%Cmlt_K_Qq_SR"), "kqqsr2"),
			`5`=list(c("U1W%KpQq_FR"), "kpqqfr"),
			`6`=list(c("Glo%tStart_VL"), "tstartvl"),
			`7`=list(c("Glo%CmltSmax_IR"), "smaxir"),
			`8`=list(c("Glo%Cmlt_K_Qq_RR"), "kqqrr"),
			`9`=list(c("Glo%Cmlt_K_Qq_FR"), "kqqfr"),
			`10`=list(c("GloTr%CmltKd_WR"), "kdwr"),
			`11`=list(c("GloTr%CmltRs_WR"), "rswr"),
			`12`=list(c("Glo%Cmlt_Pmax_ED"), "pmaxed"),
			`13`=list(c("Glo%Cmlt_E"), "cmlte"),
			`14`=list(c("GloTr%CmltSlOne_IR"), "sloneir"),
			`15`=list(c("GloTr%CmltSlTwo_IR"), "sltwoir"),
			`16`=list(c("Glo%Cmlt_AlQq_FR"), "alqqfr"),
			`17`=list(c("Glo%Cmlt_AlQq_SR"), "alqqsr"),
			`18`=list(c("Glo%Cmlt_P"), "cmltp"),
			`19`=list(c("none","none")))
	sel <- do.call(switch, cases)
	which.timedep <- sel[[1]]
	tag <- sel[[2]]
	tag <- paste0(tag,"_",tag.vrs)
	if(cse %in% c(2,3,4,6,7,9,13,14,16,18)) tag <- paste0(tag,"_adptintrv")
	##if(cse %in% c()) tag <- paste0(tag,"_adptintrv")
	##tag <- paste0(tag,"_tautd")

	errmod.conc		   <- "E1"
	struct             <- "Str07h2xTimedep"
	sett  			   <- paste0("settings_",struct,"_A1_QE3c07i_T",errmod.conc,"b02c08i_SE1c08.json")
	su <- sudriv.setup(settings=sett)
	su$settings$file.par <- paste("p",su$settings$structure,su$settings$par.tag,".par",sep="")
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
	su <- model.setup(su, settings=sett, writeO=TRUE, f.path.hru=f.path.hru, f.path.transK=NA, f.path.hyperfun=f.path.hyperfun) ## ATTENTION: HRU file is read here! It has to agree with the object that is loaded. (f.path.config).
	load("../output/A1Str07h2x/v05_connred3aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExc0/su_v05_connred3aggr_QE3P7c07i_TTE1b02c08i_SE1c08_priorExc0.RData")
	sclshifts <- list("Glo%Cmlt_Dspl_SD" = c(1-0.5, 0.5),
	                  "Glo%tStart_VL" = c(2, 0),
	                  "Glo%Cmlt_P" = c(1.4-0.9,0.9))
	if(any(which.timedep %in% names(sclshifts))){
	  scaleshift = matrix(unlist(sclshifts[names(sclshifts) %in% which.timedep]), nrow=sum(names(sclshifts) %in% which.timedep), byrow=TRUE)
	  rownames(scaleshift) <- paste0(names(sclshifts)[names(sclshifts) %in% which.timedep],ifelse(f_mean,"_fmean",""))
	  print(scaleshift)
	}else{
	  scaleshift <- NA
	}
	if(infer){
	  ags <- prepare.timedepargs(su=su,
	                             which.timedep = which.timedep,
	                             remove.taumax = remove.taumax,
	                             fix.taumax = fix.taumax,
	                             f_mean = f_mean,
	                             sclshifts = sclshifts)
	}
	if(restart==TRUE){
	  load(paste0("../output/timedeppar/result_",tag,".RData"))
		if(!("result" %in% ls())) result <- res
		## retrieve sudriv object that was saved as argument of previous infer.timedeppar
		result$loglikeli <- wrap.loglik
		if(adapt.intrv){ ## adapt the density of the intervals to the apparent acceptence frequency of the previous result
			tag <- paste0(tag, "_adptintrv")
			control$splitmethod="weighted"
			acf <- accept.frequ.get(result)
			control$interval.weights <- list()
			for(i in names(acf)){
			  control$interval.weights[[i]] <- -1*acf[[i]]+max(acf[[i]])+10
			  if(length(result$control$interval.weights)>0) control$interval.weights[[i]] <- control$interval.weights[[i]]*result$control$interval.weights[[i]]
			  miin <- min(control$interval.weights[[i]])
			  control$interval.weights[[i]] <- rollmean(control$interval.weights[[i]], k=500, fill=c(miin,NA,miin))
			  print("weights:")
			  print(summary(control$interval.weights[[i]]))
			}
		}
		result <- infer.timedeppar(loglikeli = wrap.loglik,
							 task = "restart",
							 n.iter = n.iter,
							 control = control,
							 file.save=paste0("../output/timedeppar/result_",tag),
							 res.infer.timedeppar=result,
							 #param.logprior = param.logprior,
							 #param.ou.logprior = param.ou.logprior,
							 verbose = verbose,
							 logposterior = result$dot.args$logposterior,
							 sudriv = result$dot.args$su,
							 scaleshift = result$dot.args$scaleshift,
							 mnprm=result$dot.args$mnprm)
		save(result, file=paste0("../output/timedeppar/result_",tag,".RData"), version=2)
	}else if(continue==TRUE){
		load(paste0("../output/timedeppar/result_",tag,".RData"))
		if(!("result" %in% ls())) result <- res
		result$loglikeli <- wrap.loglik
		## retrieve sudriv object that was saved as argument of previous infer.timedeppar
		result <- infer.timedeppar(loglikeli = wrap.loglik,
							 task = "continue",
							 n.iter = n.iter,
							 control = control,
							 file.save=paste0("../output/timedeppar/result_",tag),
							res.infer.timedeppar=result,
							 verbose = verbose,
							 logposterior = result$dot.args$logposterior,
							 sudriv = result$dot.args$sudriv,
							 scaleshift = result$dot.args$scaleshift,
							 mnprm=result$dot.args$mnprm)
		save(result, file=paste0("../output/timedeppar/result_",tag,".RData"), version=2)
	}else if(infer){
	  su <- ags$su
		## get covariance matrix of the posterior sample of loaded sudriv object
		ch.len <- dim(su$parameter.sample)[1]
		b <- aperm(su$parameter.sample[(ch.len-500):ch.len,,], perm = c(1,3,2))
		samp <- matrix(as.numeric(b), ncol = dim(su$parameter.sample)[2])
		## drop parameters that are time-dependent
		dm <- dimnames(su$parameter.sample)[[2]]
		colnames(samp) <- dm
		for(td.curr in which.timedep){
		  if(f_mean){# add a fictional sample of the fmean parameter
		    if(!(td.curr %in% dm)){
		      if(td.curr=="Glo%tStart_VL"){
		        samp <- cbind(samp, rtruncnorm(nrow(samp),0.01,1.99,as.numeric(su$model$prior$distdef[[td.curr]][2]),as.numeric(su$model$prior$distdef[[td.curr]][3])))
		      }else if(td.curr=="Glo%Cmlt_P"){
		        samp <- cbind(samp, rtruncnorm(nrow(samp),0.9,1.4,as.numeric(su$model$prior$distdef[[td.curr]][2]),as.numeric(su$model$prior$distdef[[td.curr]][3])))
		      }else{
		        samp <- cbind(samp, rnorm(nrow(samp),as.numeric(su$model$prior$distdef[[td.curr]][2]),as.numeric(su$model$prior$distdef[[td.curr]][3])))
		      }
		      dm <- c(dm, paste0(td.curr,"_fmean"))
		      colnames(samp) <- dm
		    }else{
		      dm <- gsub(td.curr,paste0(td.curr,"_fmean"),dm)
		      colnames(samp) <- dm
		    }
		    if(td.curr %in% names(sclshifts)){ # adapt the sample to consider the sigmoid transformation
		      samp[,td.curr==gsub("_fmean","",colnames(samp))] <- sigm.trans.inv(samp[,td.curr==gsub("_fmean","",colnames(samp))], sclshifts[[td.curr]][1], sclshifts[[td.curr]][2])
		    }
		  }
		  print(summary(ags$param.ini[[td.curr]]))
		}
		samp <- samp[,match(names(ags$param.ini)[!(names(ags$param.ini) %in% which.timedep)], dm)]
		keep <-  (colnames(samp) %in% colnames(samp)[!(colnames(samp) %in% which.timedep)])
		cov.prop.const.ini <- cov(samp[keep,keep])
		print("param.ini.const:")
		print(ags$param.ini[sapply(ags$param.ini, length) <=1])
		print(ags$mnprm)
		result <- infer.timedeppar(loglikeli = wrap.loglik,
								 param.ini = ags$param.ini,
								 param.range = ags$param.range,
								 param.log = ags$param.log,
								 param.logprior = param.logprior,
								 param.ou.ini = ags$param.ou.ini,
								 param.ou.fixed = ags$param.ou.fixed,
								 param.ou.logprior = param.ou.logprior,
								 n.iter = n.iter,
								 cov.prop.const.ini = cov.prop.const.ini,
								scale.prop.ou.ini = rep(0.005,length(which.timedep)),
								control = control,
								 file.save=paste0("../output/timedeppar/result_",tag),
								 verbose = verbose,
								 logposterior = logposterior,
								 sudriv = su,
								 scaleshift = scaleshift,
								 mnprm=ags$mnprm)
		save(result, file=paste0("../output/timedeppar/result_",tag,".RData"), version=2)
	}
	#if(adapt.intrv & !grepl("adptintrv",tag)) tag <- paste0(tag, "_adptintrv")
	if(plot){
		a <- tryCatch({
			load(paste0("../output/timedeppar/result_",tag,".RData"))
			},error=function(cond){
				cond
			},warning=function(cond){
				cond
			}
		)
		print(a)
		print(inherits(a,"error"))
		print(inherits(a,"warning"))
		if(!(inherits(a, "error") | inherits(a,"warning"))){
			if(!("result" %in% ls())) result <- res
			result$loglikeli <- wrap.loglik
			dir.create(paste0("../output/timedeppar/A1Str07h2x/",tag), recursive=TRUE)
			pdf(paste0("../output/timedeppar/A1Str07h2x/",tag,"/plot_",tag,".pdf"))
			if(grepl("tautd",tag)){
			chns <- c(12000,13000,14000)
			}else{
			chns <- c(4000,6000,8000)
			}
			plot(result, nrow=2, chains.at=chns)
			dev.off()	
			## compare change in loglikeli with change in sd_ou to inform a good prior on sd_ou
			## loglik <- result$sample.logpdf[,"loglikeliobs"]
			## sd_ou <- result$sample.param.ou[,paste0(which.timedep,"_sd")]
			## x <- 1:length(loglik)
			## lm.loglik <- lm(loglik~x)
			## lm.sd_ou <- lm(sd_ou~x)
			## tabl <- data.frame(slp.loglik=coef(lm.loglik)[2], slp.sd_ou=coef(lm.sd_ou)[2], icpt.sd_ou =coef(lm.sd_ou)[1], mnou=mean(result$sample.param.ou[,paste0(which.timedep,"_mean")]), last.sd.ou=mean(tail(sd_ou, n=10)))
			## tabl <- signif(tabl, 3)
			## tabl$param <- which.timedep
			## write.table(tabl, file="../output/timedeppar/sd_ou_loglik.txt", append=ifelse(cse==run.these[1],FALSE,TRUE),row.names=FALSE,col.names=ifelse(cse==run.these[1],TRUE,FALSE))
		}
	}
	if(save.su){
		a <- tryCatch({
			load(paste0("../output/timedeppar/result_",tag,".RData"))
			},error=function(cond){
				cond
			},warning=function(cond){
				cond
			}
		)
		if(!(inherits(a, "error") | inherits(a,"warning"))){
			if(!("result" %in% ls())) result <- res
			result$loglikeli <- wrap.loglik
			su <- result$dot.args$sudriv ## retrieve sudriv object that was saved as argument of previous infer.timedeppar
			##result$param.maxpost <- result$param.maxpost[unlist(lapply(result$param.maxpost, length))<=1]
			su <- select.maxlikpars.timedep(sudriv=su, res.timedep=result, scaleshift=result$dot.args$scaleshift, lik.not.post=TRUE)
			print("selected parameters:")
			print(su$model$parameters[su$model$par.fit==1])
			su$parameter.sample.timdedep <- result$sample.param.timedep
			su$parameter.sample.const   <- result$sample.param.const
			su$predicted$det    <- sampling_wrapper_timedep(su, sample.par=FALSE, n.sample=1, sample.likeli=FALSE, scaleshift=result$dot.args$scaleshift, mnprm=NULL)
			dir.create(paste0("../output/timedeppar/A1Str07h2x/",tag), recursive=TRUE)
			save(su, file=paste0("../output/timedeppar/A1Str07h2x/",tag,"/su_",tag,".RData"), version=2)
		}
	}
	if(find.pattern){
		load(paste0("../output/timedeppar/A1Str07h2x/", tag, "/su_", tag, ".RData"))
		vars <- c("C1Wv_Qstream", "U5F1Wv_Ss1", "U3F1Wv_Su1")
		su <- find.pattern.timedep(sudriv=su, vars=vars, scaleshift=scaleshift, tag=tag)
		##save(su, file=paste0("../output/timedeppar/A1Str07h2x/",tag,"/su_",tag,".RData"), version=2)
	}
	if(table.logpdf){
		a <- tryCatch({
			load(paste0("../output/timedeppar/result_",tag,".RData"))
			},error=function(cond){
				message(cond)
				next
			},warning=function(cond){
				message(cond)
			},finally={}
		)
		if(!(inherits(a, "error") | inherits(a,"warning"))){
			if(!("result" %in% ls())) result <- res
			if(grepl("none",tag)){
				logpdf.curr <- cbind(as.data.frame(result$sample.logpdf),logprior_oupar=NA,logpdfou_timedeppar=NA,var=tag)
			}else{
				logpdf.curr <- cbind(as.data.frame(result$sample.logpdf),var=tag)
			}
			colnames(logpdf.curr) <- c("logposterior","loglikeliobs","logprior_constpar","logprior_oupar","logpdfou_timedeppar","var")
			write.table(x = logpdf.curr, file="sample_logpdf.txt", append=ifelse(cse==run.these[1],FALSE,TRUE),row.names=FALSE,col.names=ifelse(cse==run.these[1],TRUE,FALSE))
		}
	}
	if(test){# test the implementation of the time-dependent parameters
		dummy <- test.timedeppar(sudriv=su)
	}
	rm(result,res,su)
}
if(table.logpdf){
	lgs <- read.table("sample_logpdf.txt", header=TRUE)
	compare.logpdfs(lgs)
}
