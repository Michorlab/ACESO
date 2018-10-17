
#' @import data.table
NULL
#' @import stats
NULL
#' @import utils
NULL
#' @import mrgsolve
NULL


#' @export
read.PKdata=function(file){
  PKdata=read.csv(file,header = T)

  if(length(setdiff(c('CONC', 'Time'),colnames(PKdata)))>0){
    missing=setdiff(c('CONC', 'Time'),colnames(PKdata))
    stop(paste("\n Column",missing,"is missing"))
  }
  return(PKdata)
}


#' @export
pk.function=function(model=cmt1_iv,dosing_schedule=NULL,parameters=NULL,variability=NULL,tend,delta=0.00005,pk.function=NULL,missing_times=NULL,output.time=NULL,output.conc=NULL,scale=1,...){
  if(!is.null(output.time)){
   if(!is.numeric(output.time) | !is.numeric(output.conc)) stop('output.time and output.conc must have numeric values')
    pk <- approxfun(output.time, output.conc)
  }else if(!is.null(pk.function)){
    pk=pk.function
  }else{
    if(!is.null(parameters)){
      model <- param(model, parameters)
    }
    if(!is.null(variability)){
      model=omat(model,variability)
    }
    if(!is.null(missing_times) & !is.null(dosing_schedule)){ #TODO: meter for si tengo load+maint
      dosing_interval=dosing_schedule@data$ii
      dose=dosing_schedule@data$amt
      tamt=seq(dosing_schedule@data$time[1],tend,dosing_interval)
      tamt=tamt[!tamt %in% missing_times]
      dosing_schedule <- ev(amt = dose, time=tamt)
    }
    out <-(mrgsim(model,events = dosing_schedule,end = tend,delta=delta,...))
    out@data$CP=out$CP*scale
    pk <- approxfun(out$time, out$CP)
  }

  return(pk)
}


mrgsolve.function=function(model=cmt1_iv,dosing_schedule,tend,parameters=NULL,variability=NULL,...){
  if(!is.null(parameters)){
    model <- param(model, parameters)
  }
  if(!is.null(variability)){
    model=omat(model,variability)
  }
  out <-as.data.frame(mrgsim(model,events = dosing_schedule,start=tend,end = tend,...))
  #pk <- approxfun(out$time, out$CP)
  return(out$CP[2])
}

#' @export
easy.mrgsim=function(model=cmt1_iv,dosing_schedule=NULL,parameters=NULL,variability=NULL,precision=0.1,tend,scale=1,missing_times=NULL,...){
  if(!is.null(parameters)){
    model <- param(model, parameters)
  }
  if(!is.null(variability)){
    model=omat(model,variability)
  }
  if(!is.null(missing_times)){
    dosing_interval=dosing_schedule@data$ii
    dose=dosing_schedule@data$amt
    tamt=seq(dosing_schedule@data$time[1],tend,dosing_interval)
    tamt=tamt[!tamt %in% missing_times]
    dosing_schedule <- ev(amt = dose, time=tamt)
  }
  out <-(mrgsim(model,events = dosing_schedule,end = tend,delta=precision,...))
  out@data$CP=out$CP*scale
  return(out)
}




# PK=read.csv('D:\\Mis Documentos\\ITZIAR\\Boston\\PK\\Erlotinib_1D_sparse.csv',header = T)
# mod <- mread("pk1", modlib())
# Estimate.PK(PK_data=PK,model=mod,initial_estimates=c(CL=1, V=100, KA=0.1))
# Estimate.PK(PK_data=PK,model=mod,initial_estimates=c(CL=1, V=100, KA=0.1),method='NLS',weighted = T)
# Estimate.PK(PK_data=PK,model=mod,initial_estimates=c(CL=1, V=100, KA=0.1),method='MLE')
#
# library(dplyr)
# data(Indometh)
# colnames(Indometh)=c('Subject', 'time', 'CONC')
# obs <- as.data.frame(Indometh) %>%
#   mutate(evid = 0, cmt = 0, ID = as.numeric(Subject))
#
# dose <- distinct(obs, ID) %>%
#   mutate(amt = 25, time = 0, CONC = NA, evid = 1, cmt = 2)
#
# #Put it back together
# PK <- bind_rows(obs, dose) %>%
#   arrange(ID, time) %>%
#   mutate(Subject = NULL)
#
# mod <- mread_cache("pk2", modlib())
#
# Estimate.PK(PK_data=PK,model=mod,initial_estimates=c(CL = 2, V2 = 4, Q = 5, V3 = 10),method='MLE',log.yaxis = T)
# fit=Estimate.PK(PK_data=PK,model=mod,initial_estimates=c(CL = 2, V2 = 4, Q = 5, V3 = 10),method='NLS',log.yaxis = T, weighted = F)
# exp(coef(fit))
# fit=Estimate.PK(PK_data=PK,model=mod,initial_estimates=c(CL = 2, V2 = 4, Q = 5, V3 = 10),method='NLS',log.yaxis = T, weighted = T)
# exp(coef(fit))

#' @export
Estimate.PK=function(PK_data,DV="CONC",model,initial_estimates,method=c('NLS','MLE'),log.yaxis=F,weighted=F){

  # Maximum likelihood
  mle_fun <- function(p, .mod, .data,.yobs, pred = FALSE,log.transform){

    p <- lapply(p,exp)

    sigma=p$sigma
    p$sigma=NULL

    .mod <- update(.mod, param=p)

    out <- mrgsim(.mod, data = .data, obsonly=TRUE)

    if(pred) return(as.data.frame(out))

    mle=-sum(dnorm((.yobs),(out$CP),(sigma), log=TRUE))
    return(mle)
  }

  # nonlinear least squares
  nls_pred <- function(p, .mod, .data,.yobs,pred = FALSE,log.transform=F){
   # .data <- cbind(.data,exp(p))
     p <- lapply(p,exp)
    .mod <- update(.mod, param=p)
    out <- mrgsim(.mod, data = .data, obsonly = TRUE)
    if(pred) return(as.data.frame(out))
    return( out$CP)
  }

  method <- match.arg(method)

  obj=nls_pred
  if(method=='MLE') obj=mle_fun

  obs  <- PK_data[PK_data$evid==0,]
  dose <- PK_data[PK_data$evid==1,]
  yobs <- obs[,DV]
  if(!weighted){
    weighted=rep(1,length(yobs))
  }
  else{
    weighted=1/yobs^2
  }

  theta <- log(initial_estimates)
  if(method=="NLS"){
    fit <- minpack.lm::nlsLM(yobs~nls_pred(p=eval(parse(text=paste0("c(",paste(names(theta),"=",names(theta),collapse = ","),")")))
                             ,.mod=model,.data=PK_data),
               start = theta, weights =  weighted,
               control=nls.control(maxiter=100,minFactor=1/4096,warnOnly = T),
              na.action = na.pass)

    est=(coef(fit))
   result=fit
  }else{
    theta=c(theta,sigma=0)
    fit <- minqa::newuoa(theta, obj, .mod = model, .data=PK_data, .yobs = yobs,log.transform=log.transform)
    est <- setNames(fit$par, names(theta))
    OFV=obj(est,model,PK_data,yobs,log.transform=log.transform)
    AIC=2*OFV+2*length(theta)
    final_est=exp(est)
    final_est=final_est[-length(final_est)]
    names(final_est)=names(initial_estimates)
    result=list(parameter_estimates=final_est,OFV=OFV,AIC=AIC)
  }

  model=update(model,end=max(PK_data$time), delta=0.5)
  pred <-as.data.frame(obj(est,model,PK_data,yobs,pred=TRUE,log.transform=log.transform))


  p<-ggplot() +
    geom_line(data=pred, aes(time,CP),lwd=2, col="firebrick") +
    geom_point(data=PK_data, aes(time,CONC), size=2, col="darkslateblue")

  if(log.yaxis==T) p<-p+scale_y_continuous(trans="log")

  print(p)

  return(result)
}

