
#' En_resistant_cells_i
#'
#' Expected number of resistant cells in clone i at time t after treatment initiation
#'
#' @param t Final time
#' @param type_0 Type-0 S4 object
#' @param type_i list with Type-i S4 object
#' @param  i number of cell type
#' @param type_icr Type-icr S4 object for the i cell type
#' @param approximation logical argument indicating if an approximation of the numerical integration method must be used or not. Default to TRUE for faster computation.
#' @return En_resistant_cells returns the number of resistant cells in clone i at time t after treatment initiation
#' @export
#' @examples
#' \dontrun{
#' #Birth rate of type 0 cells as a function of time:
#'b0=function(time){0.05*sin(0.1*time)+0.1}
#'#Create Type-0 S4 object structure with the parameters of type 0 sensitive cells
#'Type0 <-define.Type0.cells(N0=100,birth_rate = b0,death_rate= 0.14)

#'#Birth rate of type 1 cells as a function of time:
#'b1=function(time){0.05*sin(0.1*time)+0.12}
#'#Create Type-i S4 object structure with the parameters of type i resistant cells
#'Type1 <- define.Typei.cells(Ni=0,birth_rate = b1,death_rate  = 0.09,mutation_rate=10^-3)
#'#Group the different type i resistant cells
#'type_i<-list(Type1)
#'#  #Call En_resistant_cells_i function
#'ni_t<-En_resistant_cells_i(t=100,type_0=Type0,type_i=type_i,i=1,approximation=F)
#' }
En_resistant_cells_i<-function(t,type_0,type_i,i,type_icr=NULL,approximation=T){
 # int.function <- match.arg(int.function)
  if(approximation){
      En_res<-En_resistant_cells_i.aprox(t,type_0,type_i,i,type_icr)
  }else{
      En_res<-En_resistant_cells_i.integrate(t,type_0,type_i,i,type_icr)
  }
  return(En_res)
}

fi<-Vectorize(function(s,bi,di){bi(s)-di(s)},"s")

En_resistant_cells_i.integrate<-function(t,type_0,type_i,i,type_icr=NULL){

  if(is.null(type_icr)) ui_cr=0
  else{ui_cr=type_icr@ui_cr}

  InnerIntegral = function(y,bi,di,ui,type_0,i,ui_cr){
    integral_0<-integrate(f0,lower=0,upper=y,b0=type_0@b0,d0=type_0@d0, subdivisions = 400L)
    n0=type_0@N0*exp(integral_0$value)
    exp(-integrate(fi,lower=0,upper=y,bi=bi,di=di,subdivisions = 400L)$value)*ui(y)*type_0@b0(y)*n0
  }

  #ui<-sapply(1:length(type_i),function(x)(type_i[[x]]@ui))
  ni<- (type_i[[i]]@Ni +
          integrate(Vectorize(InnerIntegral,"y") , subdivisions = 200,
                    lower=0, upper=t,bi=type_i[[i]]@bi,di=type_i[[i]]@di,ui=type_i[[i]]@ui,type_0=type_0,i=i,ui_cr=ui_cr)$value)/exp(-integrate(fi,lower=0,upper=t,bi=type_i[[i]]@bi,di=type_i[[i]]@di,subdivisions = 400)$value)
  return(ni)
}


En_resistant_cells_i.pracma<-function(t,type_0,type_i,i){

  InnerIntegral = function(y,bi,di,ui,type_0,i){
    integral_0<-pracma::integral(f0,xmin=0,xmax=y,b0=type_0@b0,d0=type_0@d0)
    n0=type_0@N0*exp(integral_0)
    exp(-pracma::integral(fi,0,y,bi=bi,di=di))*ui(y)*type_0@b0(y)*n0
  }

  #ui<-sapply(1:length(type_i),function(x)(type_i[[x]]@ui))
  ni<- (type_i[[i]]@Ni +
          pracma::quad(Vectorize(InnerIntegral,"y") ,
                       0, t,tol=10^-5,bi=type_i[[i]]@bi,di=type_i[[i]]@di,ui=type_i[[i]]@ui,type_0=type_0,i=i))/exp(-pracma::integral(fi,0,t,bi=type_i[[i]]@bi,di=type_i[[i]]@di))
  return(ni)
}

En_resistant_cells_i.aprox<-function(t,type_0,type_i,i,type_icr=NULL){
  if(t==0) return(type_i[[i]]@Ni)
  if(is.null(type_icr)) ui_cr=0
  else{ui_cr=type_icr@ui_cr}
  #ui<-sapply(1:length(type_i),function(x)(type_i[[x]]@ui))
  int_f0<- approxfun(seq(0,t+1,1), sapply(seq(0,t+1,1),function(x) integrate(f0,lower=0,upper=x,b0=type_0@b0,d0=type_0@d0, subdivisions = 400L,stop.on.error = F)$value))
  int_fi<- approxfun(seq(0,t+1,1), sapply(seq(0,t+1,1),function(x) integrate(fi,lower=0,upper=x,bi=type_i[[i]]@bi,di=type_i[[i]]@di,subdivisions = 400L,stop.on.error = F)$value))

  InnerIntegral = function(y,ui,type_0,i){
    n0=type_0@N0*exp(int_f0(y))
    exp(-int_fi(y))*ui(y)*type_0@b0(y)*n0
  }

  ni<- (type_i[[i]]@Ni +
          integrate(Vectorize(InnerIntegral,"y") , subdivisions = 400L,
                    lower=0, upper=t,ui=type_i[[i]]@ui,type_0=type_0,i=i,stop.on.error = F)$value)/exp(-int_fi(t))
  return(ni)
}
