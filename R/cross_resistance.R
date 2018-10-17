
fi_cr<-Vectorize(function(s,bi_cr,di_cr){bi_cr-di_cr},"s")
###### Expected number cross-resistant cells: asumo que los birth y death rates no cambian con el tiempo
#' En_crossresistant_cells_i
#'
#' Expected number of cross-resistant cells arising from resistant clone i at time t after treatment initiation
#'
#' @param t Final time
#' @param type_icr Type-icr S4 object
#' @param type_0 Type-0 S4 object
#' @param type_i list with Type-i S4 object
#' @param  i cell type number from which the cross-resistance is originating
#' @param approximation logical argument indicating if an approximation of the numerical integration method must be used or not. Default to TRUE for faster computation.
#' @return En_crossresistant_cells_i returns the number of cross-resistant cells that arise from resistant clone i at time t after treatment initiation
#' @export
En_crossresistant_cells_i<-function(t,type_icr,type_0,type_i,i,approximation=T){
  if(approximation){
    ni<-En_crossresistant_cells_i.approx(t,type_icr,type_0,type_i,i)
  }else{
    ni<-En_crossresistant_cells_i.integrate(t,type_icr,type_0,type_i,i)
  }
  return(ni)
}

En_crossresistant_cells_i.integrate<-function(t,type_icr,type_0,type_i,i){

  InnerIntegral = function(y,type_icr,type_0,type_i,i){
    ni=En_resistant_cells_i(y,type_0=type_0,type_i=type_i,i=i,type_icr=type_icr)
    exp(-integrate(fi_cr,lower=0,upper=y,bi_cr=type_icr@bi_cr,di_cr=type_icr@di_cr,subdivisions = 400L)$value)*type_icr@ui_cr*type_i[[i]]@bi(y)*ni
  }

  ni<- (type_icr@Ni_cr +
          integrate(Vectorize(InnerIntegral,"y") , subdivisions = 200,
                    lower=0, upper=t,type_icr=type_icr,type_0=type_0,type_i=type_i,i=i)$value)/exp(-integrate(fi_cr,lower=0,upper=t,bi_cr=type_icr@bi_cr,di_cr=type_icr@di_cr,subdivisions = 400)$value)
  return(ni)
}

En_crossresistant_cells_i.approx<-function(t,type_icr,type_0,type_i,i){
  if(t==0) return(type_icr@Ni_cr)
  int_fi_cr<- approxfun(seq(0,t+1,1), sapply(seq(0,t+1,1),function(x) integrate(fi_cr,lower=0,upper=x,bi_cr=type_icr@bi_cr,di_cr=type_icr@di_cr,subdivisions = 400L)$value))

  InnerIntegral = function(y,type_icr,type_0,type_i,i){
    ni=En_resistant_cells_i.aprox(t=y,type_0=type_0,type_i=type_i,i=i,type_icr=type_icr)
    exp(-int_fi_cr(y))*type_icr@ui_cr*type_i[[i]]@bi(y)*ni
  }

  ni<- (type_icr@Ni_cr +
          integrate(Vectorize(InnerIntegral,"y") , subdivisions = 400L,
                    lower=0, upper=t,type_icr=type_icr,type_0=type_0,type_i=type_i,i=i)$value)/exp(-int_fi_cr(t))
  return(ni)
}
