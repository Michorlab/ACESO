
#Heaviside step function: to avoid negative results when drug concentrations are high.
HeavisideStep.f = function(x){
  x[x<0]=0
  return(x)
}

dr.function=function(model,CONC,CONC2=NULL){ #TODO: ampliar para combinaciones con mas de dos farmacos
  if(!is.null(model$fct$fct)){
    #pm <- t(model$parmMat)
    pm <-t(as.matrix(model$coefficients))
    dr_predict=model$fct$fct(CONC,pm)
    #if(model$dataList$names$orName!='Net_growth'){ #Para eliminar valores negativos, pero en caso de Net_growth puedo tener valores negativos
    dr_predict=HeavisideStep.f(dr_predict)
    #}
  }else if(inherits(model, "gam")){ #GAM model
    dr_predict=predict.gam.simplified(object=model,newdata=data.frame(CONC=CONC,CONC2=CONC2))
    #dr_predict=HeavisideStep.f(dr_predict)
  }else{
    if(is.null(CONC2)){
      dr_predict=predict(model,newdata=data.frame(CONC)) #Using this option is very slow in the integration part
    }else{
      dr_predict=predict(model,newdata=data.frame(CONC=CONC,CONC2=CONC2))
    }

    dr_predict=HeavisideStep.f(dr_predict) #TODO: no eliminar valores negativos en net_growth
  }
  return(dr_predict)
}

predict.gam.simplified=function(object,newdata){
  nb <- length(object$coefficients)
  Terms <- list(delete.response(object$pterms))
  b.size=length(newdata[[1]])
  n.smooth <- length(object$smooth)
  drop.ind <- attr(object$nsdf, "drop.ind")
  fit=array(0, b.size)
  start <- 1
  stop <- start + b.size - 1
  X <- matrix(0, b.size, nb + length(drop.ind))
  for (i in 1:length(Terms)) {
    mf <- model.frame(Terms[[i]], newdata, xlev = object$xlevels)

    Xp <- model.matrix(Terms[[i]], mf, contrasts = object$contrasts)

    if (object$nsdf[i] > 0)
      X[, 1:object$nsdf[i]] <- Xp

  }
  if (!is.null(drop.ind))
    X <- X[, -drop.ind]

  if (n.smooth)
    for (k in 1:n.smooth) {
      Xfrag <- PredictMat(object$smooth[[k]], newdata)
      X[, object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
    }

  fit[start:stop]<- X %*% object$coefficients
  return(as.numeric(fit))
}
