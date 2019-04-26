

#' read.celldeath.file
#'
#' Read the csv file where the apoptosis assay data is saved.
#'
#' @param filename name of the csv file. string.
#' @param sep the field separator character. Comma by default.
#' @param combination boolean argument to specify if the file has drug combination data, that is, more that one drug concentration column. Default to FALSE.
#' @param column.name string value to specify the column name where the total number or the fraction of death cells numbers are saved.
#' @return A dataframe with the relevant information about the cell viability assay.
#' @export
read.celldeath.file=function(filename,sep=",",combination=F,column.name="Apoptotic.fraction"){
  Apoptosis_time0<-Apoptotic.fraction<-Time<-Cell.line<-Type<-Replicate<-NULL
  cell_death_data=read.csv(filename,header = T, na='.',sep=sep)
  "."<-"callate"
  #Apoptotic.fraction: Annexin.Vplus: sum of early and late apoptotic cells (V+/Pi- + V+/Pi+)
  if(length(setdiff(c('CONC', column.name,'Time'),colnames(cell_death_data)))>0){
    missing=setdiff(c('CONC', column.name,'Time'),colnames(cell_death_data))
    stop(paste("\n Column",missing,"is missing"))
  }else if(combination){
    if(length(setdiff(c('CONC2'),colnames(cell_death_data)))>0){
      stop(paste("\n Column CONC2 is missing"))
    }
    if(!any(cell_death_data$CONC==0 | cell_death_data$CONC2==0)){
      stop('The dataset should provide information about each drug effect alone AND in combination')
    }
  }

  # Is there a replicate column?
  if(!any(colnames(cell_death_data)=="Replicate")) cell_death_data$Replicate=1
  # Is there a cell Type column?
  if(!(any(colnames(cell_death_data)=="Type"))){
    cell_death_data$Type=0
  }
  #Is there a cell line column?
  if(!(any(colnames(cell_death_data)=="Cell.line"))){
    cell_death_data$Cell.line='Cell line 1'
  }

  # #New column with time 0 data
  dtdata=data.table::as.data.table(cell_death_data)
  if(any(cell_death_data$Time==0)){
    dtdata[,Apoptosis_time0:=.(Apoptotic.fraction[Time==0]),by = .(Cell.line,Type,Replicate)]
    Time0=dtdata[dtdata$Time==unique(dtdata$Time)[is.finite(unique(dtdata$Time))][2],]
    Time0$Time=0
    Time0[,column.name]=Time0$Apoptosis_time0
    dtdata=dtdata[dtdata$Time!=0,]
    dtdata=rbind(Time0,dtdata)
    #   cell_death_data=as.data.frame(dtdata)
  }
  return(dtdata)

}

#' calculate_death_rate
#'
#' Function to calculate the death rate from the apoptosis assay data.
#'
#' @param cell_death_data name of the dataframe where the apoptosis assay data is saved.
#' @param net_growth_data name of the dataframe where the growth kinetics data is saved.
#' @param d_start initial parameter value for the death rate to be used in the least squared estimation process.
#' @param column.name string value to specify the column name where the total number or the fraction of death cells numbers are saved.
#' @param Apoptotic.fraction if TRUE, the fraction of death cells is used for the calculation of the death rate. If FALSE (default), the total number of death cells is used.
#' @param ... Additional arguments for nls function
#' @return A dataframe where the death rate for each cell type is included as a column to the information used as input.
#' @export
calculate_death_rate=function(cell_death_data,net_growth_data=NULL,d_start=0.005,Apoptotic.fraction=F,column.name="Apoptotic.cells",...){
  . <- "callate"
  Birth_rate<- Death_rate<-Cell.line<-Type <-NULL
  if(!is.null(net_growth_data)) cell_death_data=merge(cell_death_data,net_growth_data)
  data=data.table::as.data.table(cell_death_data)
  data$Death_rate=NULL
  #data[, Death_rate:= as.numeric(Death_rate)]
  if(Apoptotic.fraction==T){
    eval(parse(text=paste0('data[,Death_rate:=coef(nls(',column.name,'~(100*d*(-1+exp(Net_growth*Time)))/(-d+(Net_growth+d)*exp(Net_growth*Time)),start=list(d=d_start),...))[1],by=.(Cell.line,Type,',noquote(paste(colnames(data)[grep("CONC",colnames(data))],collapse=",")),')]')))
  }else{ #TODO: use column name?
    eval(parse(text=paste0('data[,Death_rate:=coef(nls(Apoptotic.cells~((Cell_Count_0*d/Net_growth)*(exp(Net_growth*Time)-1)),start=list(d=d_start),...))[1],by=.(Cell.line,Type,',noquote(paste(colnames(data)[grep("CONC",colnames(data))],collapse=",")),')]')))
  }

  data$Birth_rate=data$Net_growth+data$Death_rate
  # Cell clearance
  if(any(data$Birth_rate<0)){
    #data$Death_rate=data$Death_rate+(-min(data$Birth_rate))
    data[,Death_rate:=Death_rate+(-min(Birth_rate)[min(Birth_rate)<0]),by=.(Cell.line,Type)]
    data[,Birth_rate:=Birth_rate+(-min(Birth_rate)[min(Birth_rate)<0]),by=.(Cell.line,Type)]
  }

  data=as.data.frame(data)
  data=data[,colnames(net_growth_data)]
  #data$Replicate=NULL
  #data=unique(data[,c('Cell.line','Type','CONC','Net_growth','Birth_rate','Death_rate')])
  data=data[order(data$Cell.line,data$Type,data$Time,data$CONC),]
  return(data)
}
