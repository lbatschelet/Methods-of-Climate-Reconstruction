## Function converted from the former sourced template fragment.
## Uses and updates objects stored in run.env.


compute_diagnostics <- function(run.env){
  local.env <- new.env(parent = run.env)
  evalq({
    result<-list()
######################################################
############## Residuals
########################################################################

# calculate residuals and AR1 -----------------------------------------------------
if(do.residuals==T){
  if(do.field==T){
    rx<-target-pred.calib.scaled
    residualsx<-apply(rx,2,sd)#rbind(residuals,rx)
    result<-c(result,list(residualsx=residualsx))
    if(do.residual.ar1==T){
      residual.ar1sx<-apply(rx,2,function(x) cor(x,c(x[-1],NA),use="complete.obs"))
      result<-c(result,list(residual.ar1sx=residual.ar1sx))
    }
  }
  if(do.index==T){
    rxc<-y.vec-y.neu.pred.calib.scaled
    residuals.contx<-sd(rxc)#c(residuals.cont,rxc)
    result<-c(result,list(residuals.contx=residuals.contx))
    if(write.residuals.all==T){
      add<-ifelse(ens==1,F,T)
      write.table(t(round(rxc,3)),residuals.all.file,quote=F,col.names=F,row.names=F,sep=";",append=add)
    }
    if(do.residual.ar1){
      residual.ar1s.contx<-acf(rxc,lag.max=1,plot=F)[[1]][2]
      result<-c(result,list(residual.ar1s.contx=residual.ar1s.contx))
    }
  }
}


######################################################
############## Verification
########################################################################


# Do verification ------------------------------------------------------
if (do.verif==T){
  
  if(do.field==T){
    estimated<-pred.scaled[double.years,][verif.years,]
    climatology <- target
    measured <- target.orig[verif.years,]
    
    if (do.re==T){
      modelIDe <- error.valid.res(measured,estimated,climatology)     
      res <- matrix(modelIDe$RE, ncol=dim(target)[2], nrow=dim(rows)[1], byrow=T)
    }
    if (do.ce==T){
      modelIDe <- error.valid.ce(measured,estimated)     
      ces <- matrix(modelIDe$RE, ncol=dim(target)[2], nrow=dim(rows)[1], byrow=T)
    }
    
    if (do.r2==T){
      modelIDe <- error.valid.r2(measured,estimated)     
      r2 <- matrix(modelIDe$RE, ncol=dim(target)[2], nrow=dim(rows)[1], byrow=T)
    }
    
  }
  
  if(do.index==T){
    estimated.cont<-y.neu.pred.scaled[double.years][verif.years]
    climatology.cont<-target.mean.calib
    measured.cont<-target.mean[verif.years]
    
    if(do.re.index==T){
      re.cont.calc<-error.valid.res.cont(measured.cont,estimated.cont,climatology.cont)
      re.cont<-rep(re.cont.calc$RE, length(dataRows)) 
    }
    if(do.ce.index==T){
      ce.cont.calc<-error.valid.ce.cont(measured.cont,estimated.cont)
      ce.cont<-rep(ce.cont.calc$RE, length(dataRows))  
    }
    
    if(do.r2.index==T){
      r2.cont.calc<-error.valid.r2.cont(measured.cont,estimated.cont)
      r2.cont<-rep(r2.cont.calc$RE, length(dataRows)) 
    }
    
    #rmse
    if(do.rmse.all==T){
      rmse.allx<-error.valid.rmse.cont(estimated.cont,measured.cont)$rmse
      result<-c(result,list(rmse.allx=rmse.allx))
    }
    
  }
  
}

  }, envir = local.env)
  if(exists('estimated', envir = local.env, inherits = FALSE)){
    assign('estimated', get('estimated', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('estimated.cont', envir = local.env, inherits = FALSE)){
    assign('estimated.cont', get('estimated.cont', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('res', envir = local.env, inherits = FALSE)){
    assign('res', get('res', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('ces', envir = local.env, inherits = FALSE)){
    assign('ces', get('ces', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('r2', envir = local.env, inherits = FALSE)){
    assign('r2', get('r2', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('re.cont', envir = local.env, inherits = FALSE)){
    assign('re.cont', get('re.cont', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('ce.cont', envir = local.env, inherits = FALSE)){
    assign('ce.cont', get('ce.cont', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('r2.cont', envir = local.env, inherits = FALSE)){
    assign('r2.cont', get('r2.cont', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('residualsx', envir = local.env, inherits = FALSE)){
    assign('residualsx', get('residualsx', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('residual.ar1sx', envir = local.env, inherits = FALSE)){
    assign('residual.ar1sx', get('residual.ar1sx', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('residuals.contx', envir = local.env, inherits = FALSE)){
    assign('residuals.contx', get('residuals.contx', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('residual.ar1s.contx', envir = local.env, inherits = FALSE)){
    assign('residual.ar1s.contx', get('residual.ar1s.contx', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('rmse.allx', envir = local.env, inherits = FALSE)){
    assign('rmse.allx', get('rmse.allx', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists("result", envir = local.env, inherits = FALSE)){
    merge.result.payloads(run.env, local.env$result)
  }
  invisible(NULL)
}





