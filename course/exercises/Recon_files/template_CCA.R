## Function converted from the former sourced template fragment.
## Uses and updates objects stored in run.env.


run_cca_reconstruction <- function(run.env){
  local.env <- new.env(parent = run.env)
  evalq({
    result<-list()
## CCA recon code translate from Matlab from Jianghao's scripts
rows<-as.matrix(dataRows)
if(is.null(dim(proxy.matrix.full))==T) proxy.matrix.full<-as.matrix(proxy.matrix.full)

overlap.years.his<-which(proxy.years %in% proxy.calib.years)
double.years<-which(alldataRows %in% overlap.years.his)
x.all<-proxy.matrix.full[alldataRows,dataCol]

s.calib<-x.all[double.years,][calib.years,]

cca_options$dp_max <- min(dim(s.calib)[2],cca_options$dp_max)
cca_options$model<-model
if (cca_options$loadparas==F) cca_options$posfile<-posfile

early        = 1:ceiling(length(calib.years)/2)
late         = (ceiling(length(calib.years)/2)+1):length(calib.years)

dp_max<-cca_options$dp_max

indices<-proxy.calib.window.years
indices[early]<-1
indices[late]<-2
cca_options$indices<-indices

if(cca_options$loadparas){
  params<-readMat(paste0(cca_options$parafileroot,model,".mat"))
  dt<-params$dt
  dp<-params$dp
  dcca<-params$dcca
  if(isTRUE(all.equal(as.vector(params$pattern),match(models$years[[model]],(startyear:endyear))))==F){
    print(paste("!!! Years from loaded CCA parameters and model years do not Match for model",model,"!!!"))
  }
}else{
  cca_cv_out<- cca_cv(target, s.calib, cca_options)
  dt<-cca_cv_out$dt_opt
  dp<-cca_cv_out$dp_opt
  dcca<-cca_cv_out$dcca_opt
}

if(cca_options$saveparas){
  writeMat(paste0(out.suffix,"/params_",model,".mat"),dp=dp,dt=dt,dcca=dcca,pattern=match(models$years[[model]],(startyear:endyear)))
}
 
if(cca_options$sample_params){
  dccas<-dcca+seq((cca_options$dcca_sample*-1),cca_options$dcca_sample)
  dccas<-dccas[which(dccas>0)]
  dts<-dt+seq((cca_options$dt_sample*-1),cca_options$dt_sample)
  dts<-dts[which(dts>0)]
  dps<-dp+seq((cca_options$dp_sample*-1),cca_options$dp_sample)
  dps<-dps[which(dps>0)]
}

if(cca_options$sample_params){
  dt<-sample(dts,1)
  dp<-sample(dps,1)
  dccas<-dccas[which(dccas<=min(dp,dt))]
  dcca<-sample(dccas,1)
}

if(dp>dim(s.calib)[2]) dp<-dim(s.calib)[2]

# Compute regression coefficients per pattern 
cca_bp_out<-cca_bp(s.calib,target,dp,dt,dcca);
B<-cca_bp_out$B
Tm<-cca_bp_out$Tm
Pm<-cca_bp_out$Pm
Ts<-cca_bp_out$Ts
Ps<-cca_bp_out$Ps


# Perform CCA reconstruction per pattern
field_rc = (x.all - repmat(Pm,dim(x.all)[1],1))%*%diag(1.0/Ps)%*%t(B)%*%diag(Ts)
pred.scaled  = field_rc  + repmat(Tm,dim(x.all)[1],1)
recon<- pred.scaled[reconRows,]
pred.calib.scaled<-pred.scaled[double.years,][calib.years,]

#### add AR(1) noise based on the residuals to the recon -----------------
if (add.arnoise==T){
  if(arnoise.version %in% c("simple","gene")==F) stop ("no correct version for arnoise addition provided")
  if(arnoise.version=="simple"){
    ##version 1 (simple addition)
    crusd<-apply(target,2,sd)
    crumean<-apply(target,2,mean)
   if(add.arnoise.spat){
     pred.scaled.noise<-add.residualnoise.cov(pred.scaled,target,double.years,calib.years,MCiterations,crusd,crumean,noiseRows,do.pc.arnoise.spat,latweights)
     if(MCiterations==1) pred.scaled.noise<-pred.scaled.noise[1,,]
   }else{
     pred.scaled.noise<-add.residualnoise(pred.scaled,target,double.years,calib.years,MCiterations,crusd,crumean,noiseRows)
   }
  }else{
    #version 2 (wahl & smerdon 2012)
    pred.scaled.noise<-add.residualnoise.gene.cca(pred.scaled[double.years,][calib.years,],target,MCiterations,s.calib,x.all,calib.years,double.years,noiseRows,dp,dt,dcca)
    if(MCiterations==1) pred.scaled.noise<-pred.scaled.noise[1,,]
  }
  if(MCiterations==1){
    if(length(reconRows)>1){
      recon<-array(pred.scaled.noise[1:length(reconRows),],dim=c(length(reconRows),dim(pred.scaled.noise)[2]))
    }else{
      recon<- array(pred.scaled.noise,dim=c(1,length(pred.scaled.noise)))
    }
  }else{
    recon<-pred.scaled.noise[1,1:length(reconRows),]
  }
  
}



# prepare early verif years -----------------------------------------------
if(do.verif.early==T | do.ens.scores==T){
  early.rows<-which(alldataRows %in% early.years)
  field_rc = (x.all[early.rows,] - repmat(Pm,length(early.rows),1))%*%diag(1.0/Ps)%*%t(B)%*%diag(Ts)
  pred.verif.scaled  = field_rc  + repmat(Tm,length(early.rows),1)
  if(do.ens.scores==T){
    if(add.arnoise==T){
      if(MCiterations>1){
        #fill directly the entire enseble
        verif.ensemblex<-aperm(pred.scaled.noise[,(length(reconRows)+1):length(noiseRows),],c(2,3,1))
        result<-c(result,list(verif.ensemblex=verif.ensemblex))
      }else{
        verif.ensx<-array(pred.scaled.noise[(length(reconRows)+1):length(noiseRows),],dim=c(length(early.rows),dim(pred.scaled.noise)[2]))
        result<-c(result,list(verif.ensx=verif.ensx))
      }
    }else{
      verif.ensx<-pred.verif.scaled
      result<-c(result,list(verif.ensx=verif.ensx))
    }
  }
}

  }, envir = local.env)
  if(exists('rows', envir = local.env, inherits = FALSE)){
    assign('rows', get('rows', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('double.years', envir = local.env, inherits = FALSE)){
    assign('double.years', get('double.years', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('recon', envir = local.env, inherits = FALSE)){
    assign('recon', get('recon', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('pred.scaled', envir = local.env, inherits = FALSE)){
    assign('pred.scaled', get('pred.scaled', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('pred.calib.scaled', envir = local.env, inherits = FALSE)){
    assign('pred.calib.scaled', get('pred.calib.scaled', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('pred.scaled.noise', envir = local.env, inherits = FALSE)){
    assign('pred.scaled.noise', get('pred.scaled.noise', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('pred.verif.scaled', envir = local.env, inherits = FALSE)){
    assign('pred.verif.scaled', get('pred.verif.scaled', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists("result", envir = local.env, inherits = FALSE)){
    merge.result.payloads(run.env, local.env$result)
  }
  invisible(NULL)
}




