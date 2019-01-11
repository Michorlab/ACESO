
#' net_growth_rate
#'
#' Function to calculate the net growth rate from the cell viability assay data.
#'
#' @param Inputdata name of the dataframe where the cell viability data is saved.
#' @param death_rate numeric value for the death rate. Default to NULL.
#' @param birth_rate numeric value for the birth rate. Default to NULL.
#' @param time0_data logical argument to specify if the cell viability data for time 0 values should be used or not. Default to TRUE.
#' @return The Inputdata dataframe with the net growth rate values for each cell type added.
#' @export
#' @examples 
#' \dontrun{
#' data(Dactolisib_Trametinib_combination)
#' growth_data<-net_growth_rate(Dactolisib_Trametinib_combination)
#' library(ggplot2)
#' ggplot(data=growth_data,aes(x=(CONC),y=Net_growth,col=factor(CONC2)))+geom_point()+
#' geom_line(linetype="dashed")+
#' theme(text = element_text(size=14))+xlab(expression(paste("[Dactolisib] (",mu,"M)", sep="")))+
#' scale_colour_discrete(name=expression(paste("[Trametinib] (",mu,"M)", sep="")))+
#' theme_classic()+ylab("Net growth (1/h)")
#'  
#' #Calculate the birth rate from the net growth rate assuming a death of 0.021 1/h
#' ###Needed for dataset with more than 1 cell line
#' d0=vector("list", length = length(unique(growth_data$Cell.line))) 
#' ###To identify each death rate with the name of the cell line
#' names(d0)=unique(as.character(growth_data$Cell.line)) 
#' ###Introduce the value of 0.021 1/h
#' d0[[1]]=0.021
#' growth_data<-net_growth_rate(growth_data,death_rate = d0)
#'
#' ggplot(data=growth_data,aes(x=(CONC),y=Birth_rate,col=factor(CONC2)))+geom_point()+
#' geom_line(linetype="dashed")+
#' theme(text = element_text(size=14))+xlab(expression(paste("[Dactolisib] (",mu,"M)", sep="")))+
#' scale_colour_discrete(name=expression(paste("[Trametinib] (",mu,"M)", sep="")))+
#' theme_classic()+ylab("Birth rate (1/h)")
#' }

net_growth_rate=function(Inputdata,death_rate=NULL,birth_rate=NULL,time0_data=T){
  Viable.cells<-NULL
  Cell_Count_0<-NULL
  Time<-NULL
  CONC<-NULL
  Type<-NULL
  "."<-"callate"

  #For multiple cell lines:
  for(l in unique(Inputdata$Cell.line)){
    data=Inputdata[Inputdata$Cell.line==l,]

    #If Cell count is measured for a unique time
    if(length(unique(data$Time))==1){
      Net_growth=with(data,log(Viable.cells/Cell_Count_0)/Time)
      data$Net_growth=Net_growth
      gd=data.table::as.data.table(data)
      #gd[,mean(Net_growth),by=.(CONC,CONC2,Type)]
      eval(parse(text=paste0('gd[,Net_growth:=mean(Net_growth),by=.(Type,',noquote(paste(colnames(data)[grep("CONC",colnames(data))],collapse=",")),')]')))
    }else{
      # Do you want to use time 0 data for model fitting?
      if(time0_data==T){
        Time0=data[data$Time==unique(data$Time)[is.finite(unique(data$Time))][1],]
        Time0$Time=0
        Time0$Viable.cells=Time0$Cell_Count_0
        gd=rbind(Time0,data)
      }else if(time0_data==F){
        gd=data[data$Time!=0,]
      }

      #For every dose and cell type: fit the linear curve
      #data$Net_growth<-NA
      gd=data.table::as.data.table(gd)
      gd[,Net_growth:=coef(lm(log(Viable.cells)~Time))[2],by=.(CONC,Type)]
      #for(i in unique(gd$CONC)[is.finite(unique(gd$CONC))]){
      #fit<-lm(log(Viable.cells)~Time, data=gd[gd$CONC==i,])
      #gd$Net_growth[data$CONC==i]<-coef(fit)[2]
      #}
    }
    data=as.data.frame(gd)
    #From net growth to birth and death rates
    data$Death_rate=NA
    data$Birth_rate=NA
    if(!is.null(death_rate)){
      if(length(death_rate)!=length(unique(Inputdata$Cell.line))) stop("Number of supplied death rates and cell lines are not equal")
      for(i in 0:(length(death_rate[[l]])-1)){
        data$Death_rate[data$Type==i]<-death_rate[[l]][i+1] #TODO:d[[1]]['Type0']
      }
      data$Birth_rate<-data$Net_growth+data$Death_rate
      #Mean of birth rates in case of multiple times:
      #dtdata=data.table::as.data.table(edited_data)
      #dtdata[,BR:= mean(Birth_rate),by = .(CONC)]
      #edited_data=as.data.frame(dtdata)
    }else if(!is.null(birth_rate)){
      if(length(birth_rate)!=length(unique(Inputdata$Cell.line))) stop("Number of supplied birth rates and cell lines are not equal")
      for(i in 0:(length(birth_rate[[l]])-1)){
        data$Birth_rate[data$Type==i]<-birth_rate[[l]][i+1]
      }
      data$Death_rate<-data$Birth_rate-data$Net_growth

    }
    if(which(unique(Inputdata$Cell.line) %in% l)==1){
      edited_data=data
    }else{
      edited_data=data.table::rbindlist(list(edited_data,data),use.names=T,fill=T)
    }

  }

  return(edited_data)

}


growth_rate_inhibition=function(data){
  Viable.cells<-NULL
  Cell_Count_0<-NULL
  Control<-NULL

  GR=with(data,
          2^(log2(Viable.cells/Cell_Count_0)/log2(Control/Cell_Count_0))-1)

  return(GR)
}

