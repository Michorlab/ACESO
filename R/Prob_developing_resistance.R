
#Probability of resistance:
#' Prob_resistance
#'
#' Probability of resistance at time t
#' @param t time of production of a resistant cell clone
#' @param type_0 Type-0 S4 object
#' @param type_i Type-i S4 object
#' @param N number of resistant cell clones
#' @return returns the probability that there exists at least one resistant cell of any type at time T after treatment iniciation
#' @export
#' @examples
#' \dontrun{
#' #Birth rate of sensitive cells as a function of time:
#'  b0 <- function(time){0.05*sin(0.1*time)+0.14}
#'  #Birth rate of resistant cells as a function of time:
#'  b1 <- function(time){0.05*sin(0.1*time)+0.12}
#'  #Create Type-0 S4 object structure with the parameters of sensitive cells
#'  Type0 <-define.Type0.cells(N0=100,birth_rate = b0,death_rate= 0.14)
#'  #Create Type-i S4 object structure with the parameters of resistant cells
#'  Type1 <- define.Typei.cells(Ni=0,birth_rate = b1,death_rate= 0.1,mutation_rate=10^-3)
#'  Group the resistant cells types in a list
#'  Type_i<-list(Type1)
#'  Call the function:
#'  Prob_resistance(t=100,type_0=Type0,type_i=Type_i,N=1)
#'  }
#'
Prob_resistance.pracma_aprox<-function(t,type_0,type_i,N){

  int_f0<- approxfun(seq(0,t+1,0.5), sapply(seq(0,t+1,0.5),function(x) pracma::integral(f0,xmin=0,xmax=x,b0=type_0@b0,d0=type_0@d0)))


  PrIntegral <- Vectorize(function(t,Text,type_0,type_i,N){
    sum(unlist(lapply(1:N,function(i){
      # integral_0<-pracma::integral(f0,xmin=0,xmax=t,b0=type_0@b0,u=sum(ui),d0=type_0@d0)
      # n0=type_0@N0*exp(integral_0)
      n0=type_0@N0*exp(int_f0(t))
      #EN_type0_cells(t=t,type_0=type_0,ui=ui)*ui[i]*type_0@b0(t)*(1-Pext(t,Text=Text,type_i=type_i,i=i))
      n0*type_i[[i]]@ui(t)*type_0@b0(t)*(1-Pext.pracma.aprox(t,Text=Text,type_i=type_i,i=i))
    }),use.names = FALSE))
  },"t")

  P0sens<-exp(-integrate(PrIntegral,lower=0,upper=t,Text=t,type_0=type_0,type_i=type_i,N=N)$value)

  P0_res<-sapply(1:N,function(x){
    if(type_i[[x]]@Ni!=0){ #Only calculate if Ni!=0
      Pext_preexisting_resist(t=t,type_i=type_i,i=x)
    }else{
      P0_res=1
    }
  })

  1-P0sens*P0_res
}

Prob_resistance.pracma<-function(t,type_0,type_i,N){
  PrIntegral <- Vectorize(function(t,Text,type_0,type_i,N) {
    #    ui<-sapply(1:length(type_i),function(x)(type_i[[x]]@ui))
    sum(unlist(lapply(1:N,function(i){

      integral_0<-pracma::integral(f0,xmin=0,xmax=t,b0=type_0@b0,d0=type_0@d0)
      #integral_0<-quadgk_f0(f0,xmin=0,xmax=t,b0=type_0@b0,u=sum(ui),d0=type_0@d0)
      n0=type_0@N0*exp(integral_0)
      n0*type_i[[i]]@ui(t)*type_0@b0(t)*(1-Pext.pracma(t,Text=Text,type_i=type_i,i=i))
    }),use.names = FALSE))
  },"t")

  P0sens<-exp(-integrate(PrIntegral,lower=0,upper=t,Text=t,type_0=type_0,type_i=type_i,N=N)$value)

  P0_res<-sapply(1:N,function(x){
    if(type_i[[x]]@Ni!=0){ #Only calculate if Ni!=0
      Pext_preexisting_resist(t=t,type_i=type_i,i=x)
    }else{
      P0_res=1
    }
  })

  1-P0sens*P0_res
}


Prob_resistance<-function(t,type_0,type_i,N){

  PrIntegral <- Vectorize(function(t,Text,type_0,type_i,N) {
    #ui<-sapply(1:length(type_i),function(x)(type_i[[x]]@ui))
    sum(unlist(lapply(1:N,function(i){
      integral_0<-integrate(f0,lower=0,upper=t,b0=type_0@b0,d0=type_0@d0)$value
      n0=type_0@N0*exp(integral_0)
      n0*type_i[[i]]@ui(t)*type_0@b0(t)*(1-Pext(t,Text=Text,type_i=type_i,i=i))
    }),use.names = FALSE))
  },"t")

  P0sens<-exp(-integrate(PrIntegral,lower=0,upper=t,Text=t,type_0=type_0,type_i=type_i,N=N)$value)

  P0_res<-sapply(1:N,function(x){
    if(type_i[[x]]@Ni!=0){ #Only calculate if Ni!=0
      Pext_preexisting_resist(t=t,type_i=type_i,i=x)
    }else{
      P0_res=1
    }
  })

  1-P0sens*P0_res
}



Prob_resistance.aprox<-function(t,type_0,type_i,N){
  #ui<-sapply(1:length(type_i),function(x)(type_i[[x]]@ui))
  int_f0<- approxfun(seq(0,t+1,1), sapply(seq(0,t+1,1),function(x) integrate(f0,lower=0,upper=x,b0=type_0@b0,d0=type_0@d0,subdivisions = 400L,stop.on.error = F)$value))


  PrIntegral <- Vectorize(function(t,Text,type_0,type_i,N) {
    # ui<-sapply(1:length(type_i),function(x)(type_i[[x]]@ui))
    sum(unlist(lapply(1:N,function(i){
      #integral_0<-integrate(f0,lower=0,upper=t,b0=type_0@b0,u=sum(ui),d0=type_0@d0)$value
      n0=type_0@N0*exp(int_f0(t))
      #EN_type0_cells(t=t,type_0=type_0,ui=ui)*ui[i]*type_0@b0(t)*(1-Pext(t,Text=Text,type_i=type_i,i=i))
      n0*type_i[[i]]@ui(t)*type_0@b0(t)*(1-Pext.aprox(t,Text=Text,type_i=type_i,i=i))
    }),use.names = FALSE))
  },"t")

  P0sens<-exp(-integrate(PrIntegral,lower=0,upper=t,Text=t,type_0=type_0,type_i=type_i,N=N)$value)

  P0_res<-sapply(1:N,function(x){
    if(type_i[[x]]@Ni!=0){ #Only calculate if Ni!=0
      Pext_preexisting_resist(t=t,type_i=type_i,i=x)
    }else{
      P0_res=1
    }
  })

  P0_res=prod(P0_res)

  1-P0sens*P0_res
}
