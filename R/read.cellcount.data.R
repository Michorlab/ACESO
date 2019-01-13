
#' @import gridExtra
NULL

check_controls=function(dtdata,combination=F){
  Control<-Viable.cells<-CONC<-Cell.line<-Type<-Time<-Replicate<-CONC2<-NULL
  . <- "callate"
  # Is there a replicate column?
  if(!any(colnames(dtdata)=="Replicate")) dtdata$Replicate=1
  # Is there a cell Type column?
  if(!(any(colnames(dtdata)=="Type"))){
    dtdata$Type=0
  }
  #Is there a cell line column?
  if(!(any(colnames(dtdata)=="Cell.line"))){
    dtdata$Cell.line='Cell line 1'
  }

  #dtdata=data.table::as.data.table(data)
  if(!any(colnames(dtdata)=="Control") & any(dtdata$CONC==0)){
    if(!combination){
      dtdata[,Control:=.(Viable.cells[CONC==0]),by = .(Cell.line,Type,Time,Replicate)]
    }else{
      dtdata[,Control:=.(Viable.cells[CONC==0 & CONC2==0]),by = .(Cell.line,Type,Time,Replicate)]
    }
  }else if(any(colnames(dtdata)=="Control")){
    # Rewrite controls as dose=0
    CONC0=dtdata[,.SD[1],by = .(Cell.line,Type,Time,Replicate)]
    CONC0$CONC<-0
    CONC0$Viable.cells<-CONC0$Control
    if(combination){
      CONC0$CONC2<-0
    }
    #Add CONC0 to data
    dtdata=rbind(dtdata,CONC0)
    dtdata=unique(dtdata)
    if(!combination){
      dtdata=dtdata[order(dtdata$Cell.line,dtdata$Type,dtdata$Time,dtdata$CONC),]
    }else{
      dtdata=dtdata[order(dtdata$Cell.line,dtdata$Type,dtdata$Time,dtdata$CONC,dtdata$CONC2),]
    }

  }else if(!any(colnames(dtdata)=="Control") & !any(dtdata$CONC==0)){
    stop("No cell count associated with concentration 0")
  }
  data=as.data.frame(dtdata)
  return(data)
}


#' read.cellcount.data
#'
#' Read the csv file where the cell viability assay data is gathered.
#'
#' @param filename name of the csv file. string.
#' @param sep the field separator character. Comma by default.
#' @param combination boolean argument to specify if the file has drug combination data, that is, more that one drug concentration column. Default to FALSE.
#' @return A dataframe with the relevant information about the cell viability assay.
#' @export
read.cellcount.data=function(filename,sep=",",combination=F){
  Cell_Count_0<-Viable.cells<-Time<-Cell.line<-Type<-Replicate<-NULL
  . <- "callate"
  dataset=read.csv(filename,header=T,na='.',sep=sep)
  if(length(setdiff(c('CONC', 'Viable.cells','Time'),colnames(dataset)))>0){
    missing=setdiff(c('CONC', 'Viable.cells','Time'),colnames(dataset))
    stop(paste("\n Column",missing,"is missing"))
  }else if(combination){
    if(length(setdiff(c('CONC2'),colnames(dataset)))>0){
     stop(paste("\n Column CONC2 is missing"))
    }
    if(!any(dataset$CONC==0 | dataset$CONC2==0)){
      stop('The dataset should provide information about each drug effect alone AND in combination')
    }
  }
  dtdata=data.table::as.data.table(dataset)
  if(any(dataset$Time==0)){
    dtdata[,Cell_Count_0:=.(Viable.cells[Time==0]),by = .(Cell.line,Type,Replicate)]
    #dtdata[,Cell_Count_0:=mean(Viable.cells[Time==0]),by = .(Cell.line,Type)]
    dtdata=dtdata[dtdata$Time!=0,]
  }

  dataset=check_controls(dtdata=dtdata,combination=combination)

  return(dataset)

}



