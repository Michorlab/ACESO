################# LOEWE ADDITIVITY ###########################

resLoewe <- function(a, b, EmaxA, EmaxB, E0A, E0B, nA, nB, EC50A, EC50B){
  inversehill <- function(E, Emax, E0, n, EC50){
    d <- ((Emax - E0) / (E - E0) - 1)^ (1 / n)
    x <- EC50 / d
    return(x)
  }
  loewefun <- function(x, a, b, EmaxA, EmaxB, E0A, E0B, nA, nB, EC50A, EC50B){
    A <- a / (inversehill(x, EmaxA, E0A, nA, EC50A))
    if (!is.finite(A)) A <- 0
    B <- b / (inversehill(x, EmaxB, E0B, nB, EC50B))
    if (!is.finite(B)) B <- 0
    return(A + B - 1)
  }
  res <- uniroot(loewefun, sort(c(E0A, E0B, EmaxA, EmaxB))[c(2, 3)] + c(1e-9, -1e-9),
                 a = a, b = b, EmaxA = EmaxA, EmaxB = EmaxB, E0A = E0A,
                 E0B = E0B, nA = nA, nB = nB, EC50A = EC50A, EC50B = EC50B)$root
  return(res)
}
#' Compute the predicted effect under Loewe model
#'
#' @param d1 dose of drug 1
#' @param d2 dose of drug 2
#' @param model hill model for drug 1 and drug 2
#' @param same.Emax logical value to indicate whether the hill model for drug 1 and drug 2 share the same Emax value (classical Loewe model)
#'
#' @return expected effect of combination under Loewe additive model
#' @export
#'
#' @examples NULL
try.loewe <- function(d1, d2, model,same.Emax=F){
  loewe.E0.Emax <- function(d1, d2, model){
    resLoewe(d1, d2, EmaxA = model[3], EmaxB = model[3],
             E0A = model[4], E0B = model[4],
             nA = model[1], nB = model[2],
             EC50A = model[5], EC50B = model[6])
  }
  loewe.E0 <- function(d1, d2, model){
    resLoewe(d1, d2, EmaxA = model[3], EmaxB = model[4],
             E0A = model[5], E0B = model[5],
             nA = model[1], nB = model[2],
             EC50A = model[6], EC50B = model[7])
  }
  
  if(same.Emax){
    t <- try(loewe.E0.Emax(d1, d2, model), silent = TRUE)
  }else{
    t <- try(loewe.E0(d1, d2, model), silent = TRUE)
  }
  if (class(t) == "try-error") res <- NA
  else res <- t
  return(res)
}

try.loewe.separate.models <- function(d1, d2, mod1, mod2){
  loewe <- function(d1, d2, mod1, mod2){
    resLoewe(d1, d2, EmaxA = mod1[2], EmaxB = mod2[2],
             E0A = mod1[3], E0B = mod2[3],
             nA = mod1[1], nB = mod2[1],
             EC50A = mod1[4], EC50B = mod2[4])
  }
  t <- try(loewe(d1, d2, mod1, mod2), silent = TRUE)
  if (class(t) == "try-error") res <- NA
  else res <- t
  return(res)
}

#' Fill the dataset with the predicted effect under Loewe model
#'
#' @param data input data
#' @param resp Name of the column in the dataset with the effect/response to be fitted.
#' @param conc1 Name of the column in the dataset with the concentrations for drug 1. 'CONC' by default.
#' @param conc2 Name of the column in the dataset with the concentrations for drug 1. 'CONC2' by default.
#' @param same.Emax logical value to indicate whether the hill model for drug 1 and drug 2 share the same Emax value (classical Loewe model). FALSE by default.
#' @param Emax.effect argument to indicate whether the maximum effect of the drug in the dataset is associated with the minimum or the maximum values of the resp column.
#'
#' @return The input dataset with a new column called 'loewe_additivity' to store the expected effect of combination under Loewe additive model
#' @export
#'
#' @examples 
#' \dontrun{
#' data(Dactolisib_Trametinib_rates)
#' head(GD)
#' GD=Loewe(data=GD,resp = 'Birth_rate')
#' head(GD)
#' rmap <-  braidReports::responseMap(loewe_additivity~CONC+CONC2,GD,logscale=T)#' 
#' plot.ResponseSurface(rmap,xl="[Dactolisib] (µM)",yl="[Trametinib] (µM)",
#' zl="Birth rate of \n sensitive cells (1/h)",title="Loewe")
#' }
Loewe<-function(data,resp='Birth_rate',conc1='CONC',conc2='CONC2',same.Emax=F,Emax.effect=c('min','max')){
  
  variable<-NULL
  data$loewe_additivity=NULL
  GD=data[,c('Cell.line',conc1,conc2,'Type',resp)]
  GD=unique(GD)
  Emax.effect <- match.arg(Emax.effect)
  
  drug1.model=curve.fit(data=GD[GD[,conc2]==0,],resp=resp,conc=conc1)
  drug2.model=curve.fit(data=GD[GD[,conc1]==0,],resp=resp,conc=conc2)
  #Marginal data: data with the effect of only one of the drugs
  Marginal.data=GD[GD[,conc1]==0 | GD[,conc2]==0,]
  Marginal.data=data.table::melt(Marginal.data, measure.vars=c(conc1,conc2))
  first.row=Marginal.data[1,]
  second.row=first.row
  second.row$variable=conc2
  Marginal.data=Marginal.data[Marginal.data$value!=0,]
  Marginal.data=rbind(first.row,second.row,Marginal.data)
  Marginal.data=Marginal.data[order(Marginal.data$variable,Marginal.data$value),]
  
  formu=paste0(resp,'~value')
  if(!same.Emax){
    #same E0
    combi <- drc::drm(formu,variable, data=Marginal.data, fct = drc::LL.4(),
                      pmodels=list(~variable-1, ~variable-1, ~1, ~variable-1))
    drug1.model$coefficients=coef(combi)[c(1,3,5,6)]
    drug2.model$coefficients=coef(combi)[c(2,4,5,7)]
  }else{
    #same E0 y Emax
    combi <- drc::drm(formu,variable, data=Marginal.data, fct = drc::LL.4(),
                      pmodels=list(~variable-1, ~1, ~1, ~variable-1))
    drug1.model$coefficients=coef(combi)[c(1,3,4,5)]
    drug2.model$coefficients=coef(combi)[c(2,3,4,6)]
  }
  
  row.conc <- unique(GD[,conc1])[-1]
  col.conc <- unique(GD[,conc2])[-1]
  # Initialize matrix
  loewe.mat <- matrix(NA,ncol=length(col.conc)+1,nrow=length(row.conc)+1)
  loewe.mat[,1]=fitted(drug1.model)
  loewe.mat[1,]=fitted(drug2.model)
  
  for (i in 1:length(col.conc)) {
    for (j in 1:length(row.conc)) {
      x1 <- row.conc[j]
      x2 <- col.conc[i]
      
      res<-try.loewe(x1, x2, model=coef(combi),same.Emax=same.Emax)
      #res<-try.loewe.separate.models(x1, x2, mod1=coef(drug1.model),mod2=coef(drug2.model))
      
      if (!is.na(res)) {
        loewe.mat[j + 1, i + 1] <- res
      } else {
        y.loewe1<-dr.function(drug1.model,CONC=x1+x2)
        y.loewe2<-dr.function(drug2.model,CONC=x1+x2)
        if(Emax.effect=='min'){
          loewe.mat[j + 1, i + 1] <- ifelse(y.loewe1 < y.loewe2, y.loewe1, y.loewe2)
        }else{
          loewe.mat[j + 1, i + 1] <- ifelse(y.loewe1 > y.loewe2, y.loewe1, y.loewe2)
        }
        
      }
      
    }
  }
  
  GD$loewe_additivity=as.vector(t(loewe.mat))
  data=merge(data,GD)
  data=data[order(data[,conc1],data[,conc2]),]
  return(data)
}


