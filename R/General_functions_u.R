

setClass("Type-0", representation(N0 = "numeric", b0 = "function", d0="function"))
setClass("Type-i", representation(Ni = "numeric", bi = "function", di="function",ui="function"))
setClass("Type-icr", representation(Ni_cr = "numeric", bi_cr = "numeric", di_cr="numeric",ui_cr="numeric",mother_cell="numeric"))

f0<-Vectorize(function(t,b0,d0){b0(t)-d0(t)},"t")
fi<-Vectorize(function(s,bi,di){bi(s)-di(s)},"s")
fi_cr<-Vectorize(function(s,bi_cr,di_cr){bi_cr-di_cr},"s")




####################################################################################
#Example: #' @examples
# b0=function(time){0.05*sin(0.1*time)+0.1}
# Type0 <-define.Type0.cells(N0=100,birth_rate = b0,death_rate= 0.14)
#
# b1=function(time){0.05*sin(0.1*time)+0.12}
# Type1 <- define.Typei.cells(Ni=0,birth_rate = b1,death_rate= 0.09,mutation_rate=10^-3)
# type_i<-list(Type1)
#
# ni_t<-En_resistant_cells_i(t=100,type_0=Type0,type_i=type_i,i=1)
# ni_t<-En_resistant_cells_i.aprox(t=100,type_0=Type0,type_i=type_i,i=1)
#
# EN_type0_cells(t=100,type_0=Type0,ui=10^-3)
#
# type_icr=define.Typeicr.cells(Ni_0=0,birth_rate=0.02,death_rate= 0.01,mutation_rate=10^-6)
# En_crossresistant_cells_i(t=100,type_icr=type_icr,type_0=Type0,type_i=type_i,i=1)
# #3.570655e-05
# En_crossresistant_cells_i.aprox(t=100,type_icr=type_icr,type_0=Type0,type_i=type_i,i=1)
# #3.570169e-05
# En_resistant_cells_i(t=100,type_0=Type0,type_i=type_i,i=1,type_icr=type_icr)
# En_resistant_cells_i.aprox(t=100,type_0=Type0,type_i=type_i,i=1,type_icr=type_icr)
# #
# ni<-sapply(0:100,function(t)En_resistant_cells_i(t=t,type_0=Type0,type_i=type_i,i=1))
#  plot(0:100,ni,type="l")
# ni<-sapply(0:100,function(t)En_resistant_cells_i.aprox(t=t,type_0=Type0,type_i=type_i,i=1))
# plot(0:100,ni,type="l")
# #All resistant cells together:
# ni<-En_resistant_cells.aprox(N=1,t=100,type_0=Type0,type_i=type_i)
#
# # #Probability of resistance
# b0 <- function(time){0.05*sin(0.1*time)+0.14}
# b1=function(time){0.05*sin(0.1*time)+0.12}
# Type0 <-define.Type0.cells(N0=100,birth_rate = b0,death_rate= 0.14)
# Type1 <- define.Typei.cells(Ni =0 , birth_rate = b1, death_rate= 0.1,mutation_rate=10^-3)
# Prob_resistance(t=100,type_0=Type0,type_i=list(Type1),N=1)
# Prob_resistance.aprox(t=100,type_0=Type0,type_i=list(Type1),N=1)
# Prob_resistance.pracma(t=100,type_0=Type0,type_i=list(Type1),N=1) #slower, but works for flat integrands
# Prob_resistance.pracma_aprox(t=100,type_0=Type0,type_i=list(Type1),N=1)
# PR_T<-matrix(NA,ncol=100)
# for(j in 1:100){
#
#   PR_T[j]=Prob_resistance(t=j,type_0=Type0,type_i=list(Type1),N=1)
#
# }
#
# PR=sapply(0:100,function(j){Prob_resistance(t=j,type_0=Type0,type_i=list(Type1),N=1)})
# plot(1:100,PR_T)

##############################################################################################
#' define.Type0.cells
#'
#' Define Type 0 cells: sensitive cell lines
#' @param N0 cell population at time 0.
#' @param birth_rate birth rate function. It can be a numeric value, a user created function or the result of a model fitting function.
#' @param death_rate death rate function. It can be a numeric value, a user created function or the result of a model fitting function.
#' @param scale scaling parameter. Numeric value.
#' @param pk.function the name of the function(s) which computes the pharmacokinetic profile of the drug(s) that affects type 0 cells.
#' @return a S4 object with all the above information gathered together.
#' @export
#' @examples
#' define.Type0.cells(N0=1,birth_rate = 0.001,death_rate  = 0.005)
define.Type0.cells=function(N0,birth_rate,death_rate,scale=1,pk.function='pk'){
  #Define birth rate:
  #in case of a constant
  if(is.numeric(birth_rate)){
    if(birth_rate<0){stop('Birth rate parameter must be positive')}
    b0.function=paste0('
    b0=function(t,scale=',scale,'){
    BR_predict=',birth_rate,'*scale
    return(BR_predict)
    }')
    eval(parse(text = b0.function))
  }else if(is.function(birth_rate)){
    b0=birth_rate
  }
  else{
    b0.function=
      paste0('b0=function(t,model=',birth_rate,',scale=',scale,'){\n',
             paste0(paste0('CONC',seq(1,length(pk.function))),'=',pk.function,'(t)\n',collapse=''),
             paste0('BR_predict=dr.function(model,',
                    sub('1','',paste0(paste0('CONC',seq(1,length(pk.function))),'=',paste0('CONC',seq(1,length(pk.function))),collapse=',')),')
    BR_predict=BR_predict*scale
    return(BR_predict)
    }'))
    eval(parse(text = b0.function))

  }
  #Define death rate:
  if(is.numeric(death_rate)){
    if(death_rate<0){stop('Death rate parameter must be positive')}
    d0.function=paste0('
    d0=function(t,scale=',scale,'){
    DR_predict=',death_rate,'*scale
    return(DR_predict)
    }')
    eval(parse(text = d0.function))
  }else if(is.function(death_rate)){
    d0=death_rate
  }else{
    d0.function=
      paste0('d0=function(t,model=',death_rate,',scale=',scale,'){\n',
             paste0(paste0('CONC',seq(1,length(pk.function))),'=',pk.function,'(t)\n',collapse=''),
             paste0('DR_predict=dr.function(model,',
                    sub('1','',paste0(paste0('CONC',seq(1,length(pk.function))),'=',paste0('CONC',seq(1,length(pk.function))),collapse=',')),')
    DR_predict=DR_predict*scale
    return(DR_predict)
    }'))

    eval(parse(text = d0.function))
  }

  #Define Type-0 cell type:

  methods::new("Type-0", N0 = N0, b0=b0, d0 = d0)
}

#' define.Typei.cells
#'
#' Define Type i cells: resistant cell line
#' @param Ni cell population at time 0.
#' @param birth_rate birth rate function. It can be a numeric value, a user created function or the result of a model fitting function.
#' @param death_rate death rate function. It can be a numeric value, a user created function  or the result of a model fitting function.
#' @param mutation_rate numeric value or function for the mutation rate.
#' @param scale scaling parameter. Numeric value.
#' @param pk.function the name of the function which computes the pharmacokinetic profile of the drug that affects type i cells.
#' @return a S4 object with all the above information gathered together.
#' @export
#' @examples
#' define.Typei.cells(Ni=0,birth_rate = 0.001,death_rate  = 0.005,mutation_rate=0.0001)

define.Typei.cells=function(Ni=0,birth_rate,death_rate,mutation_rate,scale=1,pk.function='pk'){ #TODO:Ni_0
  #Define birth rate:
  #in case of a constant
  if(is.numeric(birth_rate)){
    bi.function=paste0('
    bi=function(t,scale=',scale,'){
    BR_predict=',birth_rate,'*scale
    return(BR_predict)
    }')
    eval(parse(text = bi.function))
  }else if(is.function(birth_rate)){
    bi=birth_rate
  }else{ #in case of a model
    bi.function=paste0('
    bi=function(t,model=',birth_rate,',scale=',scale,'){
    CONC=',pk.function,'(t)
    BR_predict=dr.function(model,CONC=CONC)
    BR_predict=BR_predict*scale
    return(BR_predict)
    }')
    eval(parse(text = bi.function))  #TODO: en caso de combinacion???
  }
  #Define death rate:
  if(is.numeric(death_rate)){
    di.function=paste0('
    di=function(t,scale=',scale,'){
      DR_predict=',death_rate,'*scale
      return(DR_predict)
    }')
    eval(parse(text = di.function))
  }else if(is.function(death_rate)){
    di=death_rate
  }else{ #in case of a model
    di.function=paste0('
    di=function(t,model=',death_rate,',scale=',scale,'){
      CONC=',pk.function,'(t)
      DR_predict=dr.function(model,CONC=CONC)
      DR_predict=DR_predict*scale
      return(DR_predict)
    }')
    eval(parse(text = di.function))
  }
  #Define mutation rate:
  if(is.numeric(mutation_rate)){
    ui.function=paste0('
                       ui=function(t){
                       ui_predict=',mutation_rate,'
                       return(ui_predict)
                       }')
    eval(parse(text = ui.function))
  }else if(is.function(mutation_rate)){
    ui=mutation_rate
  }else{ #in case of a model
    ui.function=paste0('
    ui=function(t,model=',mutation_rate,'){
                       CONC=',pk.function,'(t)
                       ui_predict=dr.function(model,CONC=CONC)
                       return(ui_predict)
  }')
    eval(parse(text = ui.function))
}

  #Define Type-i cell type:
  methods::new("Type-i", Ni = Ni, bi=bi, di = di,ui=ui)
}

#' define.Typeicr.cells
#'
#' Define Type icr cells: cross-resistant cell line
#' @param Ni_0 cell population at time 0.
#' @param birth_rate numeric value for the birth rate of cross-resistant cells.
#' @param death_rate numeric value for the death rate of cross-resistant cells.
#' @param mutation_rate numeric value for mutation rate.
#' @param scale scaling parameter. Numeric value.
#' @param mother_cell number corresponding to the mother cell type from which type icr cell has emerged.
#' @return a S4 object with all the above information gathered together.
#' @export
#' @examples
#' define.Typeicr.cells(Ni_0=0,birth_rate = 0.001,death_rate  = 0.005,mutation_rate=0.0001)
#'
define.Typeicr.cells=function(Ni_0=0,birth_rate,death_rate,mutation_rate,mother_cell=1, scale=1){
  if(!is.numeric(mutation_rate)){
    stop('Mutation rate should have a numeric value')
  }
  if(!is.numeric(birth_rate)){
    stop('birth_rate should have a numeric value')
  }
  if(!is.numeric(death_rate)){
    stop('death_rate should have a numeric value')
  }
  
  #Define Type-icr cell type:
  methods::new("Type-icr", Ni_cr = Ni_0, bi_cr=birth_rate, di_cr = death_rate,ui_cr=mutation_rate,mother_cell=mother_cell)
}
