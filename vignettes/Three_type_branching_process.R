## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warnings= FALSE
)

## ------------------------------------------------------------------------
library(ACESO)
library(ggplot2)

## ------------------------------------------------------------------------
data(Alpelisib_Trametinib_combination)
head(Alpelisib_Trametinib_combination)

DrugA=as.character(Alpelisib_Trametinib_combination$Drug.Name[1])
DrugB=as.character(Alpelisib_Trametinib_combination$Drug2.Name[1])
print(c(DrugA,DrugB))

## ------------------------------------------------------------------------
ggplot(data=Alpelisib_Trametinib_combination,aes(x=(CONC),y=Viable.cells,col=factor(CONC2)))+
  geom_point(size=1.5)+
  xlab(paste0("[",DrugA,"] (µM)"))+scale_colour_discrete(name=paste0("[",DrugB,"] (µM)"))+
  theme_classic()+theme(text = element_text(size=14))

## ------------------------------------------------------------------------
growth_data<-net_growth_rate(Alpelisib_Trametinib_combination)
head(growth_data)

## ------------------------------------------------------------------------
ggplot(data=growth_data,aes(x=(CONC),y=Net_growth,col=factor(CONC2)))+geom_point()+geom_line(linetype="dashed")+
  theme(text = element_text(size=14))+xlab(paste0("[",DrugA,"] (µM)"))+
  scale_colour_discrete(name=paste0("[",DrugB,"] (µM)"))+
  theme_classic()+ylab("Net growth (1/h)")+theme(text = element_text(size=14))

## ------------------------------------------------------------------------
d0=vector("list", length = length(unique(growth_data$Cell.line)))
names(d0)=unique(as.character(growth_data$Cell.line))
d0[[1]]=-min(growth_data$Net_growth)
d0

## ------------------------------------------------------------------------
growth_data<-net_growth_rate(growth_data,death_rate = d0)
ggplot(data=growth_data,aes(x=(CONC),y=Birth_rate,col=factor(CONC2)))+geom_point()+geom_line(linetype="dashed")+
  theme(text = element_text(size=14))+xlab(paste0("[",DrugA,"] (µM)"))+
  scale_colour_discrete(paste0("[",DrugB,"] (µM)"))+
  theme_classic()+ylab("Birth rate of \n sensitive cells, b0 (1/h)")+theme(text = element_text(size=14))


## ------------------------------------------------------------------------
#We are going to remove the columns that are not necessary for the analysis of the data in this example:
GD=growth_data[,c('Cell.line','CONC','CONC2','Net_growth','Type','Birth_rate','Death_rate')]
GD=unique(GD)

## ------------------------------------------------------------------------
rmap <- responseMap(Birth_rate~CONC+CONC2,GD,logscale=T,interpolate=FALSE)
DifferenceSurface.plot(rmap,zcenter=max(growth_data$Birth_rate)/2,
                       xl=paste0("[",DrugA,"] (µM)"),
                       yl=paste0("[",DrugB,"] (µM)"),
                       zl="Birth rate of \n sensitive cells, b0 (1/h)",
                       mid="yellow",low="hotpink1",high="darkturquoise")


## ------------------------------------------------------------------------
rmap <- responseMap(Birth_rate~CONC+CONC2,GD,logscale=T)

ResponseSurface.plot(rmap,xl=paste0("[",DrugA,"] (µM)"),
                     yl=paste0("[",DrugB,"] (µM)"),
                     zl="Birth rate of \n sensitive cells, b0 (1/h)",
                     palette=c("hotpink1","yellow","darkturquoise"))


## ------------------------------------------------------------------------
gam.model=Multiple.resp.surface.fit(GD, title=", GAM",Drug1.name=paste0("[",DrugA,"] (µM)"),Drug2.name=paste0("[",DrugB,"] (µM)"))

## ------------------------------------------------------------------------
 data(Alpelisib_sensitive)
 head(Alpelisib_sensitive)
 ggplot(data=Alpelisib_sensitive,aes(x=CONC,y=Birth_rate))+geom_point()+theme_bw()


## ------------------------------------------------------------------------
#Some error messages might appear because not all the models tested are able to fit the data. Ignore them.
B1=Multiple.best.singlefit(data=Alpelisib_sensitive,resp='Birth_rate')
B1

## ------------------------------------------------------------------------
 data(Trametinib_sensitive)
 head(Trametinib_sensitive)
 ggplot(data=Trametinib_sensitive,aes(x=CONC2,y=Birth_rate))+geom_point()+theme_bw()


## ------------------------------------------------------------------------
B2=Multiple.best.singlefit(data=Trametinib_sensitive,resp='Birth_rate',conc = "CONC2")
B2

## ------------------------------------------------------------------------
best.singlefit(Trametinib_sensitive,resp="Birth_rate",conc = "CONC2",compare=T)

## ------------------------------------------------------------------------
library(mrgsolve)
#Alpelisib
cmt1_oral<- mread("1cmt_ev", model_library())
#Trametinib: 
cmt2_oral<- mread("2cmt_ev", model_library())
see(cmt2_oral)


## ------------------------------------------------------------------------
e2 <- ev(amt = 2, time=0, ii=1, addl=30)
easy.mrgsim(model=cmt2_oral,dosing_schedule=e2,delta=0.1,tend=30,parameters = list(TVCL=4.91*24,TVV=214,TVKA=2.05*24,Vp=568,Q=60*24),scale=1000/615.4,Req="CP") %>% plot(scales=list(cex=1.5))

## ------------------------------------------------------------------------
Type0 <-define.Type0.cells(N0=10^6,birth_rate = 'gam.model[[1]]',death_rate= 0.03,scale=24,pk.function = c('pk1','pk2'))

Type1 <-define.Typei.cells(Ni=1000,birth_rate = 'B1[[1]]',death_rate  = 0.015,mutation_rate=10^-7,scale=24,pk.function='pk1')

Type2 <-define.Typei.cells(Ni=1000,birth_rate = 'B2[[1]]',death_rate  = 0.01,mutation_rate=10^-7,scale=24,pk.function='pk2')

## ------------------------------------------------------------------------
#Alpelisib: 300mg/day
e1 <- ev(amt = 300, ii = 1, addl = 60, time=0)
pk1=pk.function(model=cmt1_oral,dosing_schedule=e1,tend=50,parameters = list(TVCL=11.5*24,TVV=118,TVKA=0.784*24,ALAG=0.489/24),scale=1000/441.47) 

#Trametinib: 2mg/day
e2 <- ev(amt = 2, time=0, ii=1, addl=60)
pk2=pk.function(model=cmt2_oral,dosing_schedule=e2,tend=50,parameters = list(TVCL=4.91*24,TVV=214,TVKA=2.05*24,Vp=568,Q=60*24),scale=1000/615.4) 

## ------------------------------------------------------------------------
EN_type0_cells(t=7,type_0=Type0,ui=c(Type1@ui,Type2@ui))

## ------------------------------------------------------------------------
En_resistant_cells(N=2,t=7,type_0=Type0,type_i=list(Type1,Type2))

