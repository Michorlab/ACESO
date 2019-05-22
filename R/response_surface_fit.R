#' @import mgcv
NULL
#' @importFrom braidReports responseMap
NULL

#### RESPONSE SURFACE FIT ####
#' resp.surface.fit
#'
#' Fits a generalized additive model (gam) or  a polynomial surface (loess) to data.
#'
#' @param data Concentration-effect dataframe.
#' @param resp Name of the column with the response values. Default is Birth_rate.
#' @param conc1 Name of the column with the concentration values for drug 1.
#' @param conc2 Name of the column with the concentration values for drug 2.
#' @param Drug1.name string to specify the name of drug 1.
#' @param Drug2.name string to specify the name of drug 2.
#' @param method array to specify if a generalized additive model (gam) or  a polynomial surface (loess) model must be fitted to the data.
#' @return This function returns a fitted gam object (see \code{\link[mgcv]{gamObject}} for a detailed description) or an object of class "loess".
#' @export
#' @examples 
#' \dontrun{
#' data(Dactolisib_Trametinib_rates)
#' gam.model=resp.surface.fit(GD,resp='Birth_rate',conc1='CONC',conc2='CONC2')
#' print(gam.model)
#' }


resp.surface.fit=function(data,resp='Birth_rate',conc1='CONC',conc2='CONC2',Drug1.name='Drug 1',Drug2.name='Drug 2',method=c('gam','loess')){
  method <- match.arg(method)

  GD=data[,c('Cell.line',conc1,conc2,'Type',resp)]
  GD=unique(GD)
  CONC=GD[,conc1]
  CONC2=GD[,conc2]

  if(method=='gam'){  #generalized additive model (GAM)
    b1 <- mgcv::gam(GD[,resp] ~ s(CONC, CONC2))  #spline based smooths
    b2 <- mgcv::gam(GD[,resp] ~ te(CONC, CONC2)) #tensor product smooths

    if(AIC(b1)<AIC(b2)){
      return(b1)
    }else{
      return(b2)
    }
  }else{ #Local Polynomial Regression Fitting (loess)
    colnames(GD)[colnames(GD)==resp]='effect'
    loess.model = loess(effect~CONC+CONC2,data=GD, control=loess.control(surface="direct"),span=1) #Fit a polynomial surface
    return(loess.model)
  }

}

#' Multiple.resp.surface.fit
#'
#' Fits a generalized additive model (gam) or  a polynomial surface (loess) to data of multiple cell lines.
#'
#' @param data Concentration-effect dataframe.
#' @param resp Name of the column with the response values. Default is Birth_rate.
#' @param conc1 Name of the column with the concentration values for drug 1.
#' @param conc2 Name of the column with the concentration values for drug 2.
#' @param Drug1.name string to specify the name of drug 1.
#' @param Drug2.name string to specify the name of drug 2.
#' @param logscale determines whether the input variables will be plotted on a logarithmic scale, or a standard linear scale. Because this function is intended primarily for visualizing combined action dose-response, a logarithmic dose-pair space is the default.
#' @param method array to specify if a generalized additive model (gam) or  a polynomial surface (loess) model must be fitted to the data.
#' @param title title for the drug.
#' @return This function returns a fitted gam object (see \code{\link[mgcv]{gamObject}} for a detailed description) or an object of class "loess".
#' @seealso \code{\link[mgcv]{gamObject}}
#' @export
#' @examples 
#' \dontrun{
#' data(Dactolisib_Trametinib_rates)
#' gam.model=Multiple.resp.surface.fit(GD,resp='Birth_rate',conc1='CONC',conc2='CONC2', 
#' title=", GAM with raw data")
#' print(gam.model)
#' }
Multiple.resp.surface.fit=function(data,resp='Birth_rate',conc1='CONC',conc2='CONC2',Drug1.name='Drug 1',Drug2.name='Drug 2',logscale=T,method=c('gam','loess'),title=""){
 Cell.line<-NULL
 Type<-NULL
 CONC<-NULL
 rmap<-NULL
 
  method <- match.arg(method)

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
      fit=resp.surface.fit(data=datacltype,resp=resp,conc1=conc1,conc2=conc2,Drug1.name=Drug1.name,Drug2.name=Drug2.name,method=method)
      fit$Type=j
      fit$Cell.line=i
      fit.result[[l+1]]=fit
      #Plot
      GD=datacltype[,c('Cell.line',conc1,conc2,'Type',resp)]
      GD=unique(GD)
      CONC=GD[,conc1]
      CONC2=GD[,conc2]
      GD$predict=predict(fit, newdata =data.frame(CONC=CONC,CONC2=CONC2))
      rmap <<- braidReports::responseMap(predict~CONC+CONC2,GD,logscale=logscale)
      p=ResponseSurface.plot(rmap=rmap,xl=Drug1.name,yl=Drug2.name,zl=resp)

      plot_fit[[l+1]]<-p+ ggplot2::ggtitle(paste(i,':', 'Type',j,title))+theme_minimal()+theme(text = element_text(size=12))
      l=l+1
    }
  }

  fit.result[sapply(fit.result, is.null)] <- NULL
  plot_fit[sapply(plot_fit, is.null)] <- NULL
  do.call("grid.arrange", c(plot_fit))
  return(fit.result)
}


