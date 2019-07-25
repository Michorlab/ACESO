## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warnings= FALSE
)

## ------------------------------------------------------------------------
#' Load libraries
library(ACESO)
library(mrgsolve)
#install mrgsolve from github: devtools::install_github("metrumresearchgroup/mrgsolve")

## ------------------------------------------------------------------------
model_library(list=T)

cmt1_iv <- mread("1cmt_iv", model_library()) %>% Req(CP)
#See the code of the model
see(cmt1_iv)

## ------------------------------------------------------------------------
param(cmt1_iv)

## ------------------------------------------------------------------------
newpar <-  list('TVCL' = 171.4, #L/h
          'TVV'  = 153.5) 
cmt1_iv<-param(cmt1_iv,newpar)

param(cmt1_iv)

## ------------------------------------------------------------------------
e1 <- ev(amt = 150, ii = 1, addl = 6, time=0)
e1

## ------------------------------------------------------------------------
easy.mrgsim(model=cmt1_iv,dosing_schedule=e1,delta=0.1,tend=7) %>% plot

## ------------------------------------------------------------------------
e2 <- seq(e1, wait = 7, e1, wait = 7)
e2
easy.mrgsim(model=cmt1_iv,dosing_schedule=e2,delta=0.1,tend=7*4) %>% plot

## ------------------------------------------------------------------------
e3 <- ev(amt = 50, ii = 1, addl = 6)
e4 <- seq(e1, e3, e1)
e4
easy.mrgsim(model=cmt1_iv,dosing_schedule=e4,delta=0.1,tend=7*3) %>% plot


## ------------------------------------------------------------------------
e5 <- ev(amt = 150, ii = 1, addl = (3*7)-1)
e6 <- seq(e5, wait = 7, e5)
easy.mrgsim(model=cmt1_iv,dosing_schedule=e6,delta=0.1,tend=7*4*2) %>% plot

## ------------------------------------------------------------------------
#High dose: 1600mg once every week
load<- ev(amt = 1600, ii = 7, addl = 0, time=0)
#Low dose: 50mg once every day, after the high dose
maintenance<- ev(amt = 50, ii = 1, addl = 5,time=1)
easy.mrgsim(model=cmt1_iv,dosing_schedule=load+maintenance,delta=0.1,tend=7) %>% plot

## ------------------------------------------------------------------------
load_maint <- ev_rep(load+maintenance, n=4)
head(load_maint)
easy.mrgsim(model=cmt1_iv,data=load_maint,delta=0.1,tend=7*4) %>% plot

## ------------------------------------------------------------------------
Nindividuals=5
easy.mrgsim(model=cmt1_iv,dosing_schedule=load+maintenance,delta=0.1,tend=7,
            variability = dmat(0.1, 0.2),nid=Nindividuals) %>% plot


## ------------------------------------------------------------------------
pk<-as.data.frame(easy.mrgsim(model=cmt1_iv,dosing_schedule=load+maintenance,delta=0.1,tend=7,
                              variability = dmat(0.1, 0.2),nid=Nindividuals) )
head(pk)

library(ggplot2)
ggplot(data=pk,aes(x=time,y=CP,col=as.factor(ID)))+geom_line(size=1.3)

## ------------------------------------------------------------------------
inf <- ev(amt = 150, time=0, rate=48)
easy.mrgsim(model=cmt1_iv,dosing_schedule =inf,delta=0.1,tend=7*2) %>% plot

## ------------------------------------------------------------------------
model_library(list=T)

cmt1_oral <- mread("1cmt_ev", model_library()) %>% Req(CP) 

oral1 <- ev(amt = 150, time=0)

easy.mrgsim(model=cmt1_oral,dosing_schedule=oral1,delta=0.1,tend=7) %>% plot

oral2 <- ev(amt = 150, ii=1, addl=6, time=0)

easy.mrgsim(model=cmt1_oral,dosing_schedule=oral2,delta=0.1,tend=7*3) %>% plot


## ------------------------------------------------------------------------
#See the parameters of the model
param(cmt1_oral)
#The parameter associated with the lag time is ALAG
easy.mrgsim(model=cmt1_oral,dosing_schedule=oral1,delta=0.1,tend=7*3,parameters=list(ALAG=2)) %>% plot


