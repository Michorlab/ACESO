
#' curve.fit
#'
#' A general model fitting function for single drug effects using drm from the drc package.
#'
#' @param data Concentration-effect dataframe.
#' @param resp Name of the column with the response values. Default is "Birth_rate".
#' @param conc Name of the column with the drug concentration values. Default to "CONC".
#' @param fct nonlinear function. Currently available functions can be found in drc package (use \code{\link[drc]{getMeanFunctions}} for a full list).
#'  By default the four-parameter log-logistic model LL.4 is used. More examples include the three- and five-parameter log-logistic models LL.3, LL.5,the Weibull model W1.4, etc.
#' @param ... Additional arguments for the model fitting function. See \code{\link[drc]{drm}}for more info.
#' @return This function returns the results of the model fitting function
#' @export
#' @examples
#' \dontrun{
#' data(AU565_dataset) 
#' #Fit the data with a 4 parameter log-logistic function:
#' curve.fit(AU565_dataset,resp="Viable.cells",conc="CONC",fct=drc::LL.4())
#' }
curve.fit=function(data,resp="Birth_rate",conc='CONC',fct=drc::LL.4(),...){

  colnames(data)[colnames(data)=='effect']='effect_2'
  if(conc!='CONC') colnames(data)[colnames(data)=='CONC']='CONC_2'
  colnames(data)[colnames(data)==resp]='effect'
  colnames(data)[colnames(data)==conc]='CONC'


  if(!all(is.na(data$effect))){
    #if(length(unique(data$effect))==1){ stop("The response must vary over the concentrations")}
    fit = drc::drm(
      effect~CONC, data=data, fct=fct, ...)
    colnames(fit$data)[2]=resp
  }else{
    stop(paste(resp,' column is not numeric'))
  }

  #plot(fit, type = "all")
  return(fit)
}

#' best.singlefit
#'
#' Function where linear and several non-linear models are used to fit the curve of single drug effects using drm from the drc package.
#'
#' @param data Concentration-effect dataframe.
#' @param resp Name of the column with the response values. Default is Birth_rate.
#' @param conc Name of the column with the drug concentration values.
#' @param IC string for supplying the information criterion to be used. "AIC" and "BIC" are the two options. "AIC" by default.
#' @param type a character string specifying the data type (parameter estimation will depend on the data type as different log likelihood function will be used).
#' @param ... Additional arguments for the selection of the best model fitting function. See ?model.select for more info.
#' @param compare logical indicating if the result of all the compared models should be returned (TRUE) or only the result of the model with the lowest AIC/BIC (FALSE, default) 
#' @return This function returns the best model fitted to the data.
#' @export
#' @examples
#' \dontrun{
#' data(AU565_dataset) 
#' #Fit the data with a 4 parameter log-logistic function:
#' best.singlefit(AU565_dataset,resp="Viable.cells",conc="CONC")
#' # The best model is the 3 parameter log-logistic function. See ?drc::LL.3
#' #To see the objective function and AIC of all the compared functions, 
#' # use the compare argument of the function:
#' best.singlefit(AU565_dataset,resp="Viable.cells",conc="CONC",compare=T)
#' }
best.singlefit=function(data,resp="Birth_rate",conc='CONC',type='continuous',IC='AIC',compare=FALSE, ...){

  if(conc!='CONC') colnames(data)[colnames(data)=='CONC']='CONC_2'
  colnames(data)[colnames(data)==conc]='CONC'

  formu=paste0(resp,'~CONC')
  icfct=AIC
  if(IC=='BIC') icfct=BIC

  if(!all(is.na(eval(parse(text=paste0('data$',resp)))))){

    #M1=suppressWarnings(model.select(eval(parse(text=formu)),data=data, list(LL.4(),L.3(),L.4(),L.5(),LL.2(),LL.3(), LL.5(), W1.2(), W1.3(), W1.4(),W2.4(), EXD.2(),EXD.3(),G.2(),G.3(),G.4(), LN.3(),LN.4()),type=type,icfct=icfct,...))
    M1=suppressWarnings(model.select(eval(parse(text=formu)),data=data, list('LL.4','L.3','L.4','L.5','LL.2','LL.3','LL.5', 'W1.2', 'W1.3', 'W1.4', 'W2.4', 'EXD.2','EXD.3','G.2','G.3','G.4', 'LN.3','LN.4'),type=type,icfct=icfct,...))
    
    #linear model
    L1=eval(parse(text=paste0('lm(',formu, ',data=data)')))
    ic=AIC(L1)
    if(IC=='BIC') ic=BIC(L1)

    if(M1[,'IC'][1]<ic){
      fit = eval(parse(text=paste0('drc::drm(',formu,', data=data,fct=',rownames(M1)[1],'(), type = "continuous")')))
      #drc::drm(effect~CONC, data=data, fct=eval(parse(text =paste0(rownames(M1)[1],'()'))), type = "continuous",...)
    }else{
      fit=L1
      fit$fct$fct=function(CONC, parm){
        parm[1]+parm[2]*CONC
      }
    }

  }else{
    stop(paste(resp,' column is not numeric'))
  }

  if(compare==T){ return(M1)}
  return(fit)
}

#' linear.fit
#'
#' Function to fit linear models using lm function.
#'
#' @param data Concentration-effect dataframe.
#' @param resp Name of the column with the response values. Default is Birth_rate.
#' @param conc Name of the column with the drug concentration values.
#' @param ... Additional arguments for the lm function.
#' @return This function returns the results of the model fitting function

#' @export
linear.fit=function(data,resp="Birth_rate",conc="CONC",...){

  colnames(data)[colnames(data)=='effect']='effect_2'
  if(conc!='CONC') colnames(data)[colnames(data)=='CONC']='CONC_2'
  colnames(data)[colnames(data)==resp]='effect'
  colnames(data)[colnames(data)==conc]='CONC'

  if(!all(is.na(data$effect))){
    fit = lm(effect~CONC, data=data,...)
    colnames(fit$model)[1]=resp
  }else{
    stop(paste(resp,' column is not numeric'))
  }
  return(fit)
}

#' modelfit.plot
#'
#' Function to plot the results of fitted models.
#'
#' @param model a model object for which prediction is desired.
#' @param linear.model logic argument to specify if the fitted model is linear (TRUE) or non-linear (FALSE).
#' @param linecol color for the prediction line.
#' @param log.x logical. Logaritmically transform x axis values? FALSE by default.
#' @return A graph with the real observations plus the prediction of a previously fitted model.

#' @export
#' @examples
#' \dontrun{
#' data(AU565_dataset) 
#' #Fit the data with a 4 parameter log-logistic function:
#' fit<-curve.fit(AU565_dataset,resp="Viable.cells",conc="CONC",fct=drc::LL.4())
#' #To see the name of other nonlinear function write \code{\link[drc]{getMeanFunctions}}
#' modelfit.plot(fit)
#' modelfit.plot(fit,log.x=F)
#' }
modelfit.plot=function(model, linear.model=FALSE,log.x=FALSE,linecol='firebrick'){
  CONC<-x<-y<-NULL
  fitdata=model$data
  effect=colnames(fitdata)[2]
  colnames(fitdata)[2]="effect"
  conc=colnames(fitdata)[1]
  colnames(fitdata)[1]="CONC"

  if(linear.model){
    fitdata=model$model
    effect=colnames(fitdata)[1]
    colnames(fitdata)[1]="effect"
  }

  mml <- data.frame(CONC = seq(min(fitdata[,1]), max(fitdata[,1]), length.out = 1000))
  fit.ggplot=data.frame(y=predict(model, newdata=mml),x=mml[,conc])
  plot_fit<-ggplot2::ggplot(data=fitdata,ggplot2::aes(x=CONC,y=effect))+ggplot2::geom_point(size=1.5)+
    ggplot2::geom_line(data=fit.ggplot,ggplot2::aes(x=(x),y=(y)),size=1.3,col=linecol)+ggplot2::xlab("Concentrations")+ggplot2::ylab(effect)+
    theme_minimal()
  if(log.x){
    plot_fit<-plot_fit+scale_x_continuous(trans="log")
  } 

  return(plot_fit)
}

#' Multiple.singlefit
#'
#' Multiple concentration-response fit
#'
#' @param data Concentration-response dataframe.
#' @param resp Response to be fitted (Y axis). Default is Birth_rate.
#' @param conc Name of the column with the drug concentration values.
#' @param fct nonlinear function. Currently available functions can be found in drc package (use \code{\link[drc]{getMeanFunctions}} for a full list).
#'  By default the four-parameter log-logistic model LL.4 is used. More examples include the three- and five-parameter log-logistic models LL.3, LL.5,the Weibull model W1.4, etc.
#' @param linear.model logical value indicating whether the curve should be fitted using a linear model. If TRUE, the fct argument is ignored. FALSE by default.
#' @param log.x logical. Logaritmically transform x axis values? FALSE by default.
#' @param ... Additional arguments for the model fitting function. See \code{\link[drc]{drm}} for more info.
#' @return This function returns the coefficient estimates of the fitted model and the corresponding plots for each cell line.
#' @export
#' @examples
#' \dontrun{
#' #Effect of Alvocidib on the cell counts of BT-20 and MCF7 cell lines.
#' 
#' filename=system.file("extdata","2cell_lines.txt",package="ACESO")
#' 
#' growth_data=read.cellcount.data(filename,sep="\t")
#' 
#' #Calculate net growth rate assuming an exponential growth of the cells:
#' growth_data<-net_growth_rate(growth_data)
#' 
#' #Fit the data with a 4 parameter log-logistic function:
#' Multiple.singlefit(growth_data,resp="Net_growth",conc="CONC",fct=drc::LL.4())
#' 
#' Multiple.singlefit(growth_data,resp="Net_growth",conc="CONC",fct=drc::LL.4(),log.x=T)
#' }
Multiple.singlefit=function(data,resp="Birth_rate",conc='CONC',fct=drc::LL.4(),linear.model=FALSE,log.x=FALSE,...){

  Type<-NULL
  Cell.line<-NULL
  CONC<-NULL
  x<-NULL
  y<-NULL
  
  fit.result<- vector("list", length = length(unique(data$Cell.line))*length(unique(data$Type)))
  plot_fit<- vector("list", length = length(unique(data$Cell.line))*length(unique(data$Type)))
  l=0

  if(!(any(colnames(data)=="Type"))){
    data$Type=0
  }
  #Is there a cell line column?
  if(!(any(colnames(data)=="Cell.line"))){
    data$Cell.line='Cell line 1'
  }

  for(i in unique(data$Cell.line)){
    datacl=subset(data,Cell.line==i)
    for(j in unique(datacl$Type)){
      datacltype=subset(datacl,Type==j)
      if(linear.model==FALSE){
        fit=curve.fit(data=datacltype,fct=fct,resp=resp,conc=conc,...)
        fit$Type=j
        fit$Cell.line=i
        #For the plot
        effect=colnames(fit$data)[2]
        colnames(fit$data)[2]="effect"
        fitdata=fit$data
        colnames(fit$data)[2]=effect
      }else{
        fit=linear.fit(data=datacltype,resp=resp,...)
        fit$Type=j
        fit$Cell.line=i
        fit$fct$fct=function(CONC, parm){
          parm[1]+parm[2]*CONC
        }
        #For the plot
        effect=colnames(fit$model)[1]
        colnames(fit$model)[1]="effect"
        fitdata=fit$model
        colnames(fit$model)[1]=effect
      }

      fit.result[[l+1]]=fit

      #Plot
      mml <- data.frame(CONC = seq(min(datacltype[,conc]), max(datacltype[,conc]), length.out = 1000))
      fit.ggplot=data.frame(y=predict(fit, newdata=mml),x=mml$CONC)
      p<-ggplot2::ggplot(data=fitdata[,c("CONC","effect")],ggplot2::aes(x=(CONC),y=effect))+ggplot2::geom_point(size=1.5)+
        ggplot2::geom_line(data=fit.ggplot,ggplot2::aes(x=(x),y=(y)),size=1.3,col='firebrick')+ggplot2::xlab("Concentrations")+ggplot2::ylab(effect)+
        ggplot2::ggtitle(paste0(i,': ', 'Type ',j))+theme_minimal()+ theme(text = element_text(size=14))
      
      if(log.x){
        p<-p+scale_x_continuous(trans="log")
      } 
      plot_fit[[l+1]]<-p
      
      l=l+1
    }
  }

  fit.result[sapply(fit.result, is.null)] <- NULL
  plot_fit[sapply(plot_fit, is.null)] <- NULL
  do.call("grid.arrange", c(plot_fit))
  return(fit.result)
}


#' Multiple.best.singlefit
#'
#'  Function where linear and several non-linear models are used to fit the curve of single drug effects on multiple cell lines at the same time.
#'
#' @param data Concentration-response dataframe.
#' @param resp Response to be fitted (Y axis). Default is Birth_rate.
#' @param conc Name of the column with the drug concentration values.
#' @param IC string for supplying the information criterion to be used. "AIC" and "BIC" are the two options. "AIC" by default.
#' @param type a character string specifying the data type (parameter estimation will depend on the data type as different log likelihood function will be used).
#' @param ... Additional arguments for the selection of the best model fitting function. See ?model.select for more info.
#' @return This function returns the best model fitted to each cell line data.
#' @export
#' @examples
#' \dontrun{
#' #Effect of Alvocidib on the cell counts of BT-20 and MCF7 cell lines.
#' 
#' filename=system.file("extdata","2cell_lines.txt",package="ACESO")
#' 
#' growth_data=read.cellcount.data(filename,sep="\t")
#' 
#' #Calculate net growth rate assuming an exponential growth of the cells:
#' growth_data<-net_growth_rate(growth_data)
#' 
#' #Fit the data with a 4 parameter log-logistic function:
#' Multiple.best.singlefit(growth_data,resp="Net_growth",conc="CONC")
#' }
Multiple.best.singlefit=function(data,resp="Birth_rate",conc='CONC',IC='AIC',type='continuous',...){
  #parameters=c()
 Type<-NULL
 Cell.line<-NULL
 CONC<-x<-y<-NULL
 
   if(!(any(colnames(data)=="Type"))){
    data$Type=0
  }
  #Is there a cell line column?
  if(!(any(colnames(data)=="Cell.line"))){
    data$Cell.line='Cell line 1'
  }

  fit.result<- vector("list", length = length(unique(data$Cell.line))*length(unique(data$Type)))
  plot_fit<- vector("list", length = length(unique(data$Cell.line))*length(unique(data$Type)))
  l=0

  for(i in unique(data$Cell.line)){
    datacl=subset(data,Cell.line==i)
    for(j in unique(datacl$Type)){
      datacltype=subset(datacl,Type==j)
      fit=best.singlefit(data=datacltype,resp=resp,conc=conc,IC=IC,type=type,...)
      fit$Type=j
      fit$Cell.line=i
      if(class(fit)!='lm'){
        #For the plot
        effect=colnames(fit$data)[2]
        colnames(fit$data)[2]="effect"
        fitdata=fit$data
        colnames(fit$data)[2]=effect
      }else{
        #For the plot
        effect=colnames(fit$model)[1]
        colnames(fit$model)[1]="effect"
        fitdata=fit$model
        colnames(fit$model)[1]=effect
      }
      fit.result[[l+1]]=fit
      #Plot
      mml <- data.frame(CONC = seq(min(datacltype[,conc]), max(datacltype[,conc]), length.out = 1000))
      fit.ggplot=data.frame(y=predict(fit, newdata=mml),x=mml$CONC)
      plot_fit[[l+1]]<-ggplot2::ggplot(data=fitdata[,c("CONC","effect")],ggplot2::aes(x=(CONC),y=effect))+ggplot2::geom_point(size=1.5)+
        ggplot2::geom_line(data=fit.ggplot,ggplot2::aes(x=(x),y=(y)),size=1.3,col='firebrick')+ggplot2::xlab("Concentrations")+ggplot2::ylab(effect)+
        ggplot2::ggtitle(paste0(i,': ', 'Type ',j))+theme_minimal()+ theme(text = element_text(size=14))
      l=l+1
    }
  }

  fit.result[sapply(fit.result, is.null)] <- NULL
  plot_fit[sapply(plot_fit, is.null)] <- NULL
  do.call("grid.arrange", c(plot_fit))
  return(fit.result)
}

#' model.select
#'
#'  Comparison of different models using the following criteria: the log likelihood value, Akaike's or Bayes information criterion (AIC) or the estimated residual standard error.
#'
#' @param formula a symbolic description of the model to be fit. Either of the form 'response ~ dose' or as a data frame with response values in first column and dose values in second column.
#' @param data Concentration-response dataframe.
#' @param fctList a list of dose-response functions to be compared.
#' @param type a character string specifying the data type (parameter estimation will depend on the data type as different log likelihood function will be used).
#' @param nested logical. TRUE results in F tests between adjacent models (in 'fctList'). Only sensible for nested models.
#' @param sorted character string determining according to which criterion the model fits are ranked.
#' @param icfct function supplying the information criterion to be used. "AIC" and "BIC" are the two options. "AIC" by default.
#' @param ... Additional arguments for the model fitting function. See \code{\link[drc]{drm}} for more info.
#' @return This function compares different models using the following criteria: the log likelihood value, Akaike's or Bayes information criterion (AIC) or the estimated residual standard error.
#' @export
model.select=function (formula, data,fctList = NULL, type='continuous',nested = FALSE, sorted = c("IC","Res var", "Lack of fit", "no"), icfct = AIC,...)
{
  #options(show.error.messages = FALSE) #The user can see in the final matrix if an error has ocurred or not when fitting the model
  sorted <- match.arg(sorted)
  if (!is.logical(nested)) {
    stop("'nested' argument takes only the values: FALSE, TRUE")
  }
  lenFL <- length(fctList)
  for(i in 1:lenFL){
    object=try(drc::drm(formula, data=data, fct=eval(parse(text=paste0("drc::",fctList[[i]],"()"))),type=type,...),silent=TRUE)
    if(!inherits(object, "try-error")){
      object$fct$name=fctList[i]
      fctList=fctList[-i] #Elimino el que ha funcionado de la lista
      break #Si no falla, seguimos
    }
  }
  lenFL <- length(fctList)
  contData <- identical(type, "continuous")
  nestedInd <- 3 + contData + nested
  mc <- match.call()
  
  retMat <- matrix(0, lenFL + 1, 3 + contData + nested)
  retMat[1, 1] <- logLik(object)
  retMat[1, 2] <- icfct(object)
  retMat[1, 3] <- modelFit(object)[2, 5]
  if (contData) {
    tryRV <- try(summary(object)$resVar, silent = TRUE)
    if (!inherits("tryRV", "try-error")) {
      retMat[1, 4] <- tryRV
    }
    else {
      retMat[1, 4] <- NA
    }
  }
  if (nested) {
    retMat[1, nestedInd] <- NA
  }
  fctList2 <- rep("", lenFL + 1)
  fctList2[1] <- object$fct$name
  if (!is.null(fctList)) {
    prevObj <- object
    for (i in 1:lenFL) {
      #tempObj <- try(update(object, fct = fctList[[i]]),  silent = TRUE) 
      #fctList2[i + 1] <- fctList[[i]]$name
      tempObj <- try(update(object, fct = eval(parse(text=paste0("drc::",fctList[[i]],"()")))),  silent = TRUE) 
      fctList2[i + 1] <- fctList[i]
      if (!inherits(tempObj, "try-error")) {
        retMat[i + 1, 1] <- logLik(tempObj)
        retMat[i + 1, 2] <- icfct(tempObj)
        retMat[i + 1, 3] <- modelFit(tempObj)[2, 5]
        if (contData) {
          tryRV2 <- try(summary(tempObj)$resVar, silent = TRUE)
          if (!inherits("tryRV2", "try-error")) {
            retMat[i + 1, 4] <- tryRV2
          }
          else {
            retMat[i + 1, 4] <- NA
          }
        }
        if (nested) {
          retMat[i + 1, nestedInd] <- anova(prevObj, 
                                            tempObj, details = FALSE)[2, 5]
        }
      }
      else {
        retMat[i + 1, ] <- NA
      }
      prevObj <- tempObj
    }
  }
  rownames(retMat) <- as.vector(unlist(fctList2))
  cnames <- c("logLik", "IC", "Lack of fit")
  if (contData) {
    cnames <- c(cnames, "Res var")
  }
  if (nested) {
    cnames <- c(cnames, "Nested F test")
  }
  colnames(retMat) <- cnames
  
  if (sorted != "no") {
    return(retMat[order(retMat[, sorted]), -3])
  }
  else {
    return(retMat[,-3])
  }
}
