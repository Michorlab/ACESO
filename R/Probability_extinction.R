#Probability of extintion:

#' Pext
#'
#' Probability of extinction at time Text of a single resistant cell produced at time t
#' @param t time of production of a resistant cell clone
#' @param Text time of extinction
#' @param type_i Type-i S4 object
#' @param i resistant clone identifier
#' @param approximation logical argument indicating if an approximation of the numerical integration method must be used or not. Default to TRUE for faster computation.
#' @param int.function integration function. Possible options are "integrate" or "pracma". The default option is integrate function in R. If "pracma" is selected, the numerical integration methods from pracma package are used (more robust but slower option).
#' @return Pext function returns the probability of extinction (between 0 and 1) at time T of a single resistant cell produced at time t
#' @export
#'

Pext <- function(t,Text,type_i,i,approximation=T,int.function=c("integrate","pracma")){
  int.function <- match.arg(int.function)
  if(approximation){
    if(int.function=="integrate"){
      P.extinction<-Pext.aprox(t,Text,type_i,i)
    }else if(int.function=="pracma"){
      P.extinction<-Pext.pracma.aprox(t,Text,type_i,i)
    }
  }else{
    if(int.function=="integrate"){
      P.extinction<-Pext.integrate(t,Text,type_i,i)
    }else if(int.function=="pracma"){
      P.extinction<-Pext.pracma(t,Text,type_i,i)
    }

  }
}
Pext.pracma <- function(t,Text,type_i,i){
  wi <- Vectorize(function(n,time,bi,di){di(n+time)-bi(n+time)},"n")

  InnerIntegralVec <- Vectorize(function(tau,t,bi,di){
    di(tau+t)*exp(pracma::integral(wi,xmin=0,xmax=(tau),time=t,bi=bi,di=di,reltol = 1e-5)) #integrate fails, better use integral from pracma package
  },"tau")

  Q<-integrate(InnerIntegralVec,lower=0,upper=Text-t,t=t,bi=type_i[[i]]@bi,di=type_i[[i]]@di)$value
  return(Q/(1+Q))
}

Pext.pracma.aprox <- function(t,Text,type_i,i){

  wi <- Vectorize(function(n,time,bi,di){di(n+time)-bi(n+time)},"n")
  int_wi<- approxfun(seq(0,Text,1), sapply(seq(0,Text,1),function(x) pracma::integral(wi,xmin=0,xmax=x,time=t,bi=type_i[[i]]@bi,di=type_i[[i]]@di)))

  InnerIntegralVec <- Vectorize(function(tau,t,bi,di){
    di(tau+t)*exp(int_wi(tau))
  },"tau")

  Q<-integrate(InnerIntegralVec,lower=0,upper=Text-t,t=t,bi=type_i[[i]]@bi,di=type_i[[i]]@di)$value
  return(Q/(1+Q))
}


Pext.integrate <- function(t,Text,type_i,i){
  wi <- Vectorize(function(n,time,bi,di){di(n+time)-bi(n+time)},"n")

  InnerIntegralVec <- Vectorize(function(tau,t,bi,di){
    di(tau+t)*exp(integrate(wi,lower=0,upper=(tau),time=t,bi=bi,di=di)$value) #integrate fails, better use integral from pracma package
  },"tau")

  Q<-integrate(InnerIntegralVec,lower=0,upper=Text-t,t=t,bi=type_i[[i]]@bi,di=type_i[[i]]@di)$value
  return(Q/(1+Q))
}

Pext.aprox <- function(t,Text,type_i,i){
  wi <- Vectorize(function(n,time,bi,di){di(n+time)-bi(n+time)},"n")
  int_wi<- approxfun(seq(0,Text,1), sapply(seq(0,Text,1),function(x) integrate(wi,lower=0,upper=x,time=t,bi=type_i[[i]]@bi,di=type_i[[i]]@di,subdivisions = 400L,stop.on.error = F)$value))

  InnerIntegralVec <- Vectorize(function(tau,t,bi,di){
    di(tau+t)*exp(int_wi(tau))
  },"tau")

  Q<-integrate(InnerIntegralVec,lower=0,upper=Text-t,t=t,bi=type_i[[i]]@bi,di=type_i[[i]]@di,stop.on.error = F)$value
  return(Q/(1+Q))
}

#' Pext_preexisting_resist
#'
#' Probability of extinction at time t of a resistant cell type when there is pre-existing resistance at treatment initiation
#' @param t time of production of a resistant cell clone
#' @param type_i Type-i S4 object
#' @param i resistant clone identifier
#' @return Pext_preexisting_resist function returns the probability of extinction (between 0 and 1) at time t of a single resistant cell produced at time t
#' @export

Pext_preexisting_resist <- function(t,type_i,i){
  wi <- Vectorize(function(n,bi,di){di(n)-bi(n)},"n")

  InnerIntegralVec <- Vectorize(function(tau,bi,di){
    di(tau)*exp(integrate(wi,lower=0,upper=(tau),bi=bi,di=di,stop.on.error = F)$value)
  },"tau")

  Q<-integrate(InnerIntegralVec,lower=0,upper=t,bi=type_i[[i]]@bi,di=type_i[[i]]@di)$value
  P0_res<-(Q/(1+Q))^type_i[[i]]@Ni
  return(P0_res)
}
