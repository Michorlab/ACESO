#' En_resistant_cells
#'
#' Expected number of resistant cells in the different clones over time after treatment initiation
#'
#' @param N Number of different resistant cell types
#' @param t Final time
#' @param type_0 Type-0 S4 object
#' @param type_i Type-i S4 object
#' @param type_icr Type-icr S4 object for the i cell type
#' @param approximation logical argument indicating if an approximation of the numerical integration method must be used or not. Default to TRUE for faster computation.
#' @return En_resistant_cells returns a matrix with the number of resistant cells over time in the different clones
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
#'Type1 <- define.Typei.cells(Ni=0,birth_rate = b1,death_rate  = 0.09,mutation_rate=10^-4)
#'#Birth rate of type 2 cells as a function of time:
#'b2=function(time){0.05*sin(0.1*time)+0.12}
#'#Create Type-i S4 object structure with the parameters of type i resistant cells
#'Type2 <- define.Typei.cells(Ni=0,birth_rate = b2,death_rate  = 0.06,mutation_rate=10^-4)
#'Call En_resistant_cells function
#'En_resistant_cells(N=2,t=100,type_0=Type0,type_i=list(Type1,Type2))
#' }
#'
En_resistant_cells<-function(N,t,type_0,type_i,type_icr=NULL,approximation=T){
  if(approximation){
    ni<-En_resistant_cells.aprox(N,t,type_0,type_i,type_icr)
  }else{
    ni<-En_resistant_cells.integrate(N,t,type_0,type_i,type_icr)
  }
  return(ni)
}

En_resistant_cells.integrate<-function(N,t,type_0,type_i,type_icr=NULL){


  InnerIntegral = function(y,bi,di,ui,type_0,i,ui_cr){
    exp(-pracma::integral(fi,lower=0,upper=y,bi=bi,di=di))*ui(y)*type_0@b0(y)*EN_type0_cells(y,type_0 =type_0)
  }

  ni<-matrix(NA,ncol=t+1,nrow=N)
  #ui<-sapply(1:length(type_i),function(x)(type_i[[x]]@ui))

  #Iterate over time:
  for(i in 1:N){
    if(is.null(type_icr)) ui_cr=0
    else{ui_cr=type_icr[[i]]@ui_cr}

    ni[i,]<- sapply(0:t,function(j){
      (type_i[[i]]@Ni +
         integrate(Vectorize(InnerIntegral,"y") ,
                   lower=0, upper=j,bi=type_i[[i]]@bi,di=type_i[[i]]@di,ui=type_i[[i]]@ui,N0=type_0@N0,b0=type_0@b0,d0=type_0@d0,i=i,ui_cr=ui_cr)$value)/exp(-integrate(fi,0,j,bi=type_i[[i]]@bi,di=type_i[[i]]@di)$value)
    })
  }
  return(ni)
}

En_resistant_cells.aprox<-function(t,N,type_0,type_i,type_icr=NULL){

  if(!is.list(type_i)) stop("type_i must be a list of resistant cell types")
  type.icr=NULL

  ni<-matrix(NA,ncol=1,nrow=N)

  #Iterate over different cell types:
  for(i in 1:N){
    #Select which type_icr corresponds to i
    if(!is.null(type_icr) & is.list(type_icr)) type.icr=type_icr[unlist(lapply(type_icr, function(x) (x@mother_cell==i)))][[1]]
    ni[i,1]<- En_resistant_cells_i.aprox(t=t,type_0 = type_0,type_i=type_i,type_icr=type.icr,i=i)

  }

  colnames(ni)=paste0("t=",t)
  return(ni)
}
