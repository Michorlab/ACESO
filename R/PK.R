
#' @import data.table
NULL
#' @import methods
NULL
#' @import stats
NULL
#' @import utils
NULL
#' @import drc
NULL
#' @importFrom mrgsolve mrgsim param omat ev mcode
NULL
#' @importFrom minpack.lm nlsLM
NULL

#' read.PKdata
#'
#' Read the pharmacokinetic data
#'
#' @param file name of the file where to drug concentrations over time are stored.
#' @return A data.frame with the pharmacokinetic data 
#' @export
read.PKdata=function(file){
  PKdata=read.csv(file,header = T)

  if(length(setdiff(c('CONC', 'Time'),colnames(PKdata)))>0){
    missing=setdiff(c('CONC', 'Time'),colnames(PKdata))
    stop(paste("\n Column",missing,"is missing"))
  }
  return(PKdata)
}

#' pk.function
#'
#' Makes an interpolation of Pharmacokinetic (PK) data or a simulated PK model over time
#'
#' @param model Structural PK model definition (mrgsolve model specification)
#' @param dosing_schedule Dosing schedule to simulate the PK model
#' @param parameters Estimates of the pharmacokinetic parameters of the model.
#' @param variability variability terms for the specified PK parameters.
#' @param tend Final time for the simulation of the PK model.
#' @param delta a number indicating the increment of a time sequence to simulate from 0 to tend. A small value is needed for a proper interpolation. Default: 0.00005.
#' @param pk.function user defined pk function to interpolate instead of interpolating the simulation results of a specified PK model
#' @param missing_times numeric vector to indicate the times where the patient missed some of the doses specified in the dosing_schedule argument.
#' @param output.time vector of times to interpolate 
#' @param output.conc vector of drug concentrations to interpolate
#' @param scale scaling parameter for y axis
#' @param scale.x scaling parameter for x axis
#' @param ... additional arguments for the mrgsim function from mrgsolve package
#' @return This function returns an interpolation of PK data or a simulated PK model 

#' @export
#' @examples
#' \dontrun{
#' # mrgsolve model specification for a one compartment model and intravenuos administration:
#'  code_iv <- '$PARAM @annotated
#'  TVCL   :  8 : Clearance (volume/time)
#'  TVV    : 80 : Central volume (volume)
#'  
#'  $CMT  @annotated
#'CENT : Central compartment
#'
#'
#' $MAIN
#' double CL = exp(log(TVCL) + ETA1);
#' double V = exp(log(TVV)  + ETA2);
#'
#' $OMEGA @labels ETA1 ETA2
#' 0 0
#' 
#' $GLOBAL
#' #define CP (CENT)
#' 
#' $PKMODEL ncmt = 1, depot = FALSE
#' 
#' $CAPTURE @annotated
#' CP : Plasma concentration (mass/volume)'
#' 
#' library(mrgsolve)
#' 
#' #Read the code:
#' cmt1_iv<- mcode("mymodel", code_iv)%>% Req(CP) 
#' 
#' #Specify dosing schedule: 150mg once a day for a month (30 days)
#' e1 <- ev(amt = 150, ii = 1, addl = 30, time=0)
#' 
#' #Define or change the parameters for the model. Use a list:
#' newpar <- list(TVCL=7.14*24,TVV=155)
#' # * TVCL: Typical value for the clearance (volume/time) 
#' # * TVV: Typical value for the volumen of distribution
#' 
#' #Pk function to be used during the simulation of the evolutionary process:
#' pk=pk.function(model=cmt1_iv,dosing_schedule=e1,tend=30,parameters = newpar)
#'
#'}
pk.function=function(model,dosing_schedule=NULL,parameters=NULL,variability=NULL,tend,delta=0.00005,pk.function=NULL,missing_times=NULL,output.time=NULL,output.conc=NULL,scale=1,scale.x=1,...){
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
    out@data$time=out$time*scale.x
    pk <- approxfun(out$time, out$CP)
  }

  return(pk)
}


mrgsolve.function=function(model,dosing_schedule,tend,parameters=NULL,variability=NULL,...){
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

#' Simulate from a model object
#'
#' Extended version of mrgsim function from mrgsolve package to compute predictions of PK models following mrgsolve model specifications.
#'
#' @param model Structural PK model definition (mrgsolve model specification)
#' @param dosing_schedule Dosing schedule to simulate the PK model
#' @param parameters Estimates of the pharmacokinetic parameters of the model.
#' @param variability variability terms for the specified PK parameters.
#' @param tend Final time for the simulation of the PK model.
#' @param precision number: increment of the sequence to simulate from 0 to tend.
#' @param missing_times numeric vector to indicate the times where the patient missed some of the doses specified in the dosing_schedule argument.
#' @param scale scaling parameter for y axis
#' @param scale.x scaling parameter for x axis
#' @param ... additional arguments for the mrgsim function from mrgsolve package
#' @return This function computes predictions of PK models following mrgsolve model specifications.
#' @export
#' @examples
#' \dontrun{
#' # mrgsolve model specification for a one compartment model and intravenuos administration:
#'  code_iv <- '$PARAM @annotated
#'  TVCL   :  8 : Clearance (volume/time)
#'  TVV    : 80 : Central volume (volume)
#'  
#'  $CMT  @annotated
#'CENT : Central compartment
#'
#'
#' $MAIN
#' double CL = exp(log(TVCL) + ETA1);
#' double V = exp(log(TVV)  + ETA2);
#'
#' $OMEGA @labels ETA1 ETA2
#' 0 0
#' 
#' $GLOBAL
#' #define CP (CENT)
#' 
#' $PKMODEL ncmt = 1, depot = FALSE
#' 
#' $CAPTURE @annotated
#' CP : Plasma concentration (mass/volume)'
#' 
#' library(mrgsolve)
#' 
#' #Read the code:
#' cmt1_iv<- mcode("mymodel", code_iv)%>% Req(CP) 
#' 
#' #Specify dosing schedule: 150mg once a day for a month (30 days). The dose is given at time 0.
#' d1 <- ev(amt = 150, ii = 1, addl = 30, time=0)
#' 
#' #Define or change the parameters for the model. Use a list:
#' newpar <- list(TVCL=7.14*24,TVV=155)
#' # * TVCL: Typical value for the clearance (volume/time) 
#' # * TVV: Typical value for the volumen of distribution
#' # To see the names of the parameters in the model:
#' param(cmt1_iv)
#' 
#' #Simulate the model with easy.mrgsim function 
#' easy.mrgsim(model=cmt1_iv,dosing_schedule=d1,tend=30,delta=0.1,parameters = newpar) %>% plot
#'
#' #Now simulate another dosign schedule: 1600mg once a week for 45 days. The dose is given at day 2.
#' d2 <- ev(amt = 1600, ii = 7, addl = 30, time=2)
#' 
#' easy.mrgsim(model=cmt1_iv,dosing_schedule=d2,tend=45,delta=0.1,parameters = newpar) %>% plot
#'
#' #Save the results as a data.frame
#' Simulation2<-as.data.frame(
#' easy.mrgsim(model=cmt1_iv,dosing_schedule=d2,tend=45,delta=0.1,parameters = newpar))
#'}
easy.mrgsim=function(model,dosing_schedule=NULL,parameters=NULL,variability=NULL,precision=0.1,tend,scale=1,scale.x=1,missing_times=NULL,...){
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
  out@data$time=out$time*scale.x
  return(out)
}




#' Estimate.PK
#'
#' Fit Pharmacokinetic (PK) data and estimate PK parameters. The use of other software like NONMEM or MONOLIX is encoraged for estimation of parameters from different individuals.
#'
#' @param PK_data Data to fit (drug concentrations over time)
#' @param x.name Column name of the independent variable (time) in the x axis. "time" by default.
#' @param y.name Column name of the dependent variable (drug concentrations) in the y axis. "CONC" by default.
#' @param model Structural PK model definition (mrgsolve model specification).
#' @param initial_estimates list with the initial estimates for the PK parameters.
#' @param log.yaxis logical to indicate if the values for the dependent variable should be log transformed. FALSE by default.
#' @param weighted logical. If TRUE a weighted NLS method is applied.
#' @return This function returns the final estimates for the PK paramters and the plot of the observation values and the prediction of the model
#' @export
#' @examples
#' \dontrun{
#' data(PK_example)
#' head(PK.example.data)
#' #Read code from model repository: 1 compartment model with extravascular administration
#' cmt1_oral<- mread("PK_models/1cmt_oral")%>% Req(CP) 
#' #See parameter names in the model:
#' param(cmt1_oral)
#' #Estimate parameters:
#  Estimate.PK(PK_data=PK.example.data,model=cmt1_oral,initial_estimates=c(TVCL=1, TVV=100, TVKA=0.1))
#' 
#' }

Estimate.PK=function(PK_data,x.name="time",y.name="CONC",model,initial_estimates,log.yaxis=F,weighted=F){
 CP<-NULL
 CONC<-NULL
 
 if(y.name!='CONC') colnames(PK_data)[colnames(PK_data)=='CONC']='CONC_2'
 colnames(PK_data)[colnames(PK_data)==y.name]='CONC'
 
 colnames(PK_data)[colnames(PK_data)==x.name]='time'

  # nonlinear least squares
  nls_pred <- function(p, .mod, .data,.yobs){
   # .data <- cbind(.data,exp(p))
     p <- lapply(p,exp)
    .mod <- mrgsolve::update(.mod, param=p)
    out <- mrgsim(.mod, data = .data, obsonly = TRUE)
    #if(pred) return(as.data.frame(out))
    return( out$CP)
  }

  obj=nls_pred
  obs  <- PK_data[PK_data$evid==0,]
  dose <- PK_data[PK_data$evid==1,]
  yobs <- obs[,"CONC"]
  if(!weighted){
    weighted=rep(1,length(yobs))
  }
  else{
    weighted=1/yobs^2
  }

  theta <- log(initial_estimates)
  
  fit <- minpack.lm::nlsLM(yobs~nls_pred(p=eval(parse(text=paste0("c(",paste(names(theta),"=",names(theta),collapse = ","),")")))
                                         ,.mod=model,.data=PK_data),
                           start = theta, weights =  weighted,
                           control=nls.control(maxiter=100,minFactor=1/4096,warnOnly = T),
                           na.action = na.pass)
  
  est=(coef(fit))
  result=fit
 
 
  model=mrgsolve::update(model,param=exp(est),end=max(PK_data$time), delta=0.5)
  pred <-as.data.frame(mrgsolve::mrgsim(model,data=PK_data[PK_data$amt!=0,],end=max(PK_data$time)))


  p<-ggplot() +
    geom_line(data=pred, aes(time,CP),lwd=2, col="firebrick") +
    geom_point(data=PK_data, aes(time,CONC), size=2, col="darkslateblue")

  if(log.yaxis==T) p<-p+scale_y_continuous(trans="log")

  print(p)

  return(result)
}


