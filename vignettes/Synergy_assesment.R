## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warnings= FALSE
)

## ------------------------------------------------------------------------
library(ACESO)

## ------------------------------------------------------------------------
data(Dactolisib_Trametinib_rates)
head(GD)

DrugA=as.character(GD$Drug.Name[1])
DrugB=as.character(GD$Drug2.Name[1])
print(c(DrugA,DrugB))

## ------------------------------------------------------------------------
rmap <- responseMap(Birth_rate~CONC+CONC2,GD,logscale=T,interpolate=FALSE)
DifferenceSurface.plot(rmap,zcenter=max(GD$Birth_rate)/2,
                       xl=" [Dactoloisib] (µM)",
                       yl="[Trametinib] (µM)",
                       zl="Birth rate of \n sensitive cells, b0 (1/h)",
                       mid="yellow",low="hotpink1",high="darkturquoise")


## ------------------------------------------------------------------------
rmap <- responseMap(Birth_rate~CONC+CONC2,GD,logscale=T)

ResponseSurface.plot(rmap,xl=" [Dactoloisib] (µM)",
                     yl="[Trametinib] (µM)",
                     zl="Birth rate of \n sensitive cells, b0 (1/h)",
                     palette=c("hotpink1","yellow","darkturquoise"))


## ------------------------------------------------------------------------
GD=Loewe(data=GD,resp = 'Birth_rate')

## ------------------------------------------------------------------------
rmap_loewe <- responseMap(loewe_additivity~CONC+CONC2,GD,logscale=T)

ResponseSurface.plot(rmap_loewe,xl=" [Dactoloisib] (µM)",
                     yl="[Trametinib] (µM)",
                     zl="Birth rate of \n sensitive cells, b0 (1/h)",
                     palette=c("hotpink1","yellow","darkturquoise"))

## ------------------------------------------------------------------------
GD$diffLoewe=(GD$loewe_additivity-GD$Birth_rate)

p=SynergyMatrix.plot(GD,resp="diffLoewe")
p+ggplot2::labs(x="Dactolisib concentration (µM)", y="Trametinib concentration (µM)")


## ------------------------------------------------------------------------
GD=HSA(GD,resp = 'Birth_rate')

GD$diffHSA=(GD$HSA_response-GD$Birth_rate)

p2=SynergyMatrix.plot(GD,resp="diffHSA")
p2+ggplot2::labs(x="Dactolisib concentration (µM)", y="Trametinib concentration (µM)")


