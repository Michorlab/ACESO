#' @importFrom pracma integral
NULL

#' EN_type0_cells
#'
#' Expected number of sensitive type 0 cells at time t after treatment initiation
#'
#' @param t Final time
#' @param T0 Initial time. Default 0.
#' @param type_0 Type-0 S4 object
#' @param ui vector containing the mutation probabilities of the resistant cell clones. Default is 0 (single-type birth-death process, no resistant clones).
#' @param int.function integration function. Possible options are "integrate" or "pracma". The default option is integrate function in R. If "pracma" is selected, the numerical integration methods from pracma package are used (more robust but slower option).
#' @return EN_type0_cells returns the number of sensitive cells at time t after treatment initiation
#' @export
#' @examples
#' \dontrun{
#' #Birth rate of type 0 cells as a function of time:
#' b0=function(time){0.05*sin(0.1*time)+0.1}
#' #Create Type-0 S4 object structure with the parameters of type 0 sensitive cells
#' Type0 <-define.Type0.cells(N0=1,birth_rate = b0,death_rate= 0.08)
#' #Call EN_type0_cells function to calculate the number of cells at time 250.
#' EN_type0_cells(t=250,type_0=Type0)
#' EN_type0_cells(T0=0,t=50,type_0=Type0)*EN_type0_cells(T0=50,t=250,type_0=Type0)
#'}
#'
EN_type0_cells<-function(T0=0,t,type_0,ui=0,int.function=c("integrate","pracma")){
  int.function <- match.arg(int.function)
  if(int.function=="integrate"){
    n0=EN_type0_cells(T0=T0,t=t,type_0=type_0,ui=ui)
  }else if(int.function=="pracma"){
    n0=EN_type0_cells.pracma(T0=T0,t=t,type_0=type_0,ui=ui)
  }
  return(n0)
}

f0<-Vectorize(function(t,b0,d0){b0(t)-d0(t)},"t")

EN_type0_cells_integrate<-function(T0=0,t,type_0,ui=0){
  integral_0<-integrate(f0,lower=T0,upper=t,b0=type_0@b0,d0=type_0@d0,subdivisions = 400L,stop.on.error = F)$value
  n0=type_0@N0*exp(integral_0)
  return(n0)
}
EN_type0_cells.pracma<-function(T0=0,t,type_0,ui=0){
  integral_0<-pracma::integral(f0,xmin=T0,xmax=t,b0=type_0@b0,d0=type_0@d0)
  n0=type_0@N0*exp(integral_0)
  return(n0)
}

