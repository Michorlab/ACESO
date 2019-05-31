#################### HIGHEST SINGLE AGENT ###################
#' Fill the dataset with the predicted effect under Highest Single Agent (HSA) model
#'
#' @param data input dataset
#' @param resp Name of the column in the dataset with the effect/response to be fitted.
#' @param conc1 Name of the column in the dataset with the concentrations for drug 1. 'CONC' by default.
#' @param conc2 Name of the column in the dataset with the concentrations for drug 1. 'CONC2' by default.
#' @param Emax.effect argument to indicate whether the maximum effect of the drug in the dataset is associated with the minimum or the maximum values of the resp column.
#'
#' @return The input dataset with a new column called 'HSA_response' to store the expected effect of combination under HSA model
#' @export
#'
#' @examples 
#' \dontrun{
#' data(Dactolisib_Trametinib_rates)
#' head(GD)
#' GD=HSA(GD,resp='Birth_rate',conc1='CONC',conc2='CONC2')
#' head(GD)
#' rmap <- braidReports::responseMap(HSA_response~CONC+CONC2,GD,logscale=T)
#' plot.ResponseSurface(rmap,xl="[Dactolisib] (µM)",yl="[Trametinib] (µM)",
#' zl="Birth rate of \n sensitive cells (1/h)",title="HSA")
#' }
HSA<-function(data,resp='Birth_rate',conc1='CONC',conc2='CONC2',Emax.effect=c('min','max')){

  data$HSA_response=NULL
  GD=data[,c('Cell.line',conc1,conc2,'Type',resp)]
  GD=unique(GD)
  Emax.effect <- match.arg(Emax.effect)
  #Monotherapy curves
  drug1.resp<- GD[GD[,conc2]==0,resp]  #GD[GD$CONC2==0,resp]
  drug2.resp <-GD[GD[,conc1]==0,resp]

  ref.mat <- matrix(NA,ncol=length(drug2.resp),nrow=length(drug1.resp))
  for (i in 1:length(drug1.resp)) {
    for (j in 1:length(drug2.resp)) {
      if(Emax.effect=='min'){
        ref.mat[i, j] <- ifelse(drug1.resp[i] < drug2.resp[j], #lower response is associated with a higher effect
                                drug1.resp[i], drug2.resp[j])
      }else{
        ref.mat[i, j] <- ifelse(drug1.resp[i] > drug2.resp[j], #higher response is associated with a higher effect
                                drug1.resp[i], drug2.resp[j])
      }

    }
  }

  GD$HSA_response=as.vector(t(ref.mat))
  data=merge(data,GD)
  data=data[order(data[,conc1],data[,conc2]),]
  return(data)

}

