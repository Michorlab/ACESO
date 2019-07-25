## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warnings= FALSE
)

## ------------------------------------------------------------------------
library(ACESO)
library(mrgsolve)
library(ggplot2)

## ------------------------------------------------------------------------
growth_data=read.cellcount.data(system.file("extdata", "cell_viability_assay.txt", package = "ACESO"), sep=";")
head(growth_data)

## ------------------------------------------------------------------------
#Exploraty plot:  viable cells over time
library(ggplot2)
ggplot(data=growth_data,aes(x=Time,y=Viable.cells,col=factor(CONC)))+geom_point()+
  geom_smooth(se=F,method="lm")+facet_wrap(~Type)+scale_y_continuous(trans="log")+ylab("Viable cells (log scale)")+scale_colour_discrete(name="Dose (µM)")

## ------------------------------------------------------------------------
cell_death=read.celldeath.file(system.file("extdata", "apoptosis_assay.txt", package = "ACESO"), column.name="Apoptotic.fraction",sep=";")
head(cell_death)

## ------------------------------------------------------------------------
library(dplyr)
my_sum <- cell_death %>%
  group_by(CONC,Time,Type2) %>%
  summarise(
    n=n(),
    mean=mean(Apoptotic.fraction),
    sd=sd(Apoptotic.fraction)
  )

ggplot(cell_death[cell_death$Time!=0,])+
  geom_bar(data=my_sum[my_sum$Time!=0,],aes(x=as.factor(Time), y=mean,fill=as.factor(CONC),color=as.factor(CONC)), stat="identity", alpha=0.7, position=position_dodge(0.9)) +
  geom_jitter(aes(x=as.factor(Time), y=Apoptotic.fraction,color=as.factor(CONC)), size=0.9, position=position_dodge(0.9))+
  facet_wrap(~factor(Type2,levels=c("sensitive","resistant")))+
  xlab("Time (hours)")+ylab("Fraction of death cells")+
  geom_errorbar(data=my_sum[my_sum$Time!=0,], aes(x=as.factor(Time), ymin=mean-sd, ymax=mean+sd,color=as.factor(CONC)),width=.2, position=position_dodge(0.9))+
  theme_bw()+scale_color_brewer(name="Erlotinib concentration (µM)",palette="Spectral")+scale_fill_brewer(name="Erlotinib concentration (µM)",palette="Spectral")+
  theme(text = element_text(size=14))+theme(legend.position="top") 

## ------------------------------------------------------------------------
growth_data<-net_growth_rate(growth_data,time0_data = F)
head(growth_data) #See the values in Net_growth column

## ------------------------------------------------------------------------
all.data=calculate_death_rate(net_growth_data=growth_data,cell_death_data =cell_death, column.name = "Apoptotic.fraction",Apoptotic.fraction = T)

## ------------------------------------------------------------------------
ggplot(all.data,aes(x=CONC,y=Death_rate))+geom_point()+scale_y_continuous(limits = c(0, 0.05))+facet_wrap(~Type)+ylab("Death rates (1/h)")+xlab("Drug concentration (µM)")

ggplot(all.data,aes(x=CONC,y=Birth_rate))+geom_point()+facet_wrap(~Type)+
  ylab("Birth rates (1/h)")+xlab("Drug concentration (µM)")

## ------------------------------------------------------------------------
#Fit birth rates
BR.fit=Multiple.best.singlefit(all.data,resp="Birth_rate")
BR.fit
#Fit death rates
DR.fit=Multiple.best.singlefit(all.data,resp="Death_rate")
DR.fit

## ------------------------------------------------------------------------
best.singlefit(all.data[all.data$Type==1,],resp="Birth_rate",compare=T)
#Some error messages might appear because not all the models tested are able to fit the data. Ignore them.

## ------------------------------------------------------------------------
B1.fit=Multiple.singlefit(all.data[all.data$Type==1,],fct=drc::LL.4(),resp="Birth_rate")

## ------------------------------------------------------------------------
B1.fit=Multiple.singlefit(all.data[all.data$Type==1,],resp="Birth_rate",linear.model = T)

## ------------------------------------------------------------------------
model_library(list=T)

cmt1_ev <- mread("1cmt_ev", model_library()) %>% Req(CP)
see(cmt1_ev)

## ------------------------------------------------------------------------
param(cmt1_ev)

## ------------------------------------------------------------------------
newpar <-  list('TVCL' = 3.95, #L/h
          'TVV'  = 233,  #L
          'TVKA' = 0.95) #h-1

## ------------------------------------------------------------------------
e1 <-  ev(amt = 150, ii = 0, addl = 0, time=0)
# amt: amount
# ii: dosing interval
# addl: additional doses
# time: time when the dose is given.

## ------------------------------------------------------------------------
easy.mrgsim(model=cmt1_ev,dosing_schedule=e1,delta=0.1,tend=48,parameters = newpar) %>% plot

## ------------------------------------------------------------------------
e2 <-  ev(amt = 150, ii = 1, addl = 30, time=0) #150 mg every day for 30 days
easy.mrgsim(model=cmt1_ev,dosing_schedule=e2,delta=0.1,tend=30,parameters = list(TVCL=3.95*24,TVV=233,TVKA=0.95*24)) %>% plot(scales=list(cex=1.5))


## ------------------------------------------------------------------------
easy.mrgsim(model=cmt1_ev,dosing_schedule=e2,delta=0.1,tend=30,parameters = list(TVCL=3.95*24,TVV=233,TVKA=0.95*24),scale=1000/429.9) %>% plot(scales=list(cex=1.5))


## ------------------------------------------------------------------------
pk=pk.function(model=cmt1_ev,dosing_schedule=e2,tend=30,parameters = list(TVCL=3.95*24,TVV=233,TVKA=0.95*24),scale=1000/429.9)

## ----eval=FALSE----------------------------------------------------------
#  pk2=function(t, Cp0,k, time_interval,time_first_dose=0){
#    n = floor((t-time_first_dose)/time_interval) + 1
#    Cp<-Cp0*(1-exp(-k*n*time_interval))*
#      exp(-k*((t-time_first_dose)-(n-1)*time_interval))/(1-exp(-k*time_interval))
#    return(Cp)
#  }

## ------------------------------------------------------------------------
Type0 <-define.Type0.cells(N0=10^6,birth_rate = 'BR.fit[[1]]',death_rate= 0.0085,scale=24,pk.function = 'pk')
#The function returns a S4 object with all the information gathered together.
Type0
#For resistant cells and additional argument is needed: mutation_rate
Type1 <-define.Typei.cells(Ni=0,birth_rate = 'B1.fit[[1]]',death_rate  = 0.0035,mutation_rate=10^-8,scale=24, pk.function = 'pk')

## ------------------------------------------------------------------------
plot(seq(0,30,0.1),Type0@b0(seq(0,30,0.1)),type='l',xlab='Time',ylab='b0')
plot(seq(0,30,0.1),Type1@bi(seq(0,30,0.1)),type='l',xlab='Time',ylab='b1')

## ------------------------------------------------------------------------
 EN_type0_cells(t=10,type_0=Type0,ui=c(Type1@ui))
 EN_type0_cells(t=10,type_0=Type0,ui=c(Type1@ui),int.function="pracma")

## ------------------------------------------------------------------------
sapply(seq(0,30,1),function(i){ EN_type0_cells(t=i,type_0=Type0,ui=c(Type1@ui))})

## ----eval=FALSE----------------------------------------------------------
#  En_resistant_cells(N=1,t=10,type_0=Type0,type_i=list(Type1))

## ----eval=FALSE----------------------------------------------------------
#  sapply(c(0,2,4,8),function(i){En_resistant_cells(N=1,t=i,type_0=Type0,type_i=list(Type1))})
#  

## ---- eval=FALSE---------------------------------------------------------
#  #Prob_resistance(t=10,type_0=Type0,type_i=list(Type1),N=1)
#  #Not executed because it is a very time consuming function:

