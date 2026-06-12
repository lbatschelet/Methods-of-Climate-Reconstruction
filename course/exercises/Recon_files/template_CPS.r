## Function converted from the former sourced template fragment.
## Uses and updates objects stored in run.env.


run_cps_index_reconstruction <- function(run.env){
  local.env <- new.env(parent = run.env)
  evalq({
    result<-list()
rows<-as.matrix(c(dataRows))

if(is.null(dim(proxy.matrix.full))==T) proxy.matrix.full<-as.matrix(proxy.matrix.full)
if(is.null(dim(proxy.matrix.calib))==T) proxy.matrix.calib<-as.matrix(proxy.matrix.calib)

overlap.years.his<-which(proxy.years %in% proxy.calib.years)

double.years<-which(alldataRows %in% overlap.years.his)

composite<-recon.cps(proxy.matrix.calib,proxy.matrix.full,dataRows,dataCol,target.mean.calib,alldataRows,sample.cps.weight,cps.weighting,cps.weight.exponent,minsc.cps,maxsc.cps,calib.years,double.years,cps.weigh.distances,do.cps.cor.weight)

y.neu.pred.scaled<-composite

y.neu.pred.calib.scaled<-composite[double.years][calib.years]

recon.cont<- composite[reconRows]

#### add AR(1) noise based on the residuals to the recon -----------------
##version 1 (simple addition)
if (add.arnoise==T){
  if(arnoise.version %in% c("simple","gene")==F) stop ("no correct version for arnoise addition provided")
  if(arnoise.version=="simple"){
    crusd<-sd(target.mean.calib)
    crumean<-mean(target.mean.calib)
    y.neu.pred.scaled.noise<-add.residualnoise.index(composite,target.mean.calib,double.years,calib.years,MCiterations,crusd,crumean,noiseRows)
  }else{
    #version 2 (wahl & smerdon 2012)
    crusd<-sd(target.mean.calib)
    crumean<-mean(target.mean.calib)
    y.neu.pred.scaled.noise<-add.residualnoise.gene.cps(y.neu.pred.calib.scaled,target.mean.calib,MCiterations,crusd,crumean,proxy.matrix.calib,proxy.matrix.full,dataRows,dataCol,alldataRows,sample.cps.weight,cps.weight.exponent,calib.years,double.years,noiseRows)
  }
  recon.cont<-as.vector(y.neu.pred.scaled.noise[1,1:length(reconRows)])
}

# early verif years -------------------------------------------------------

if(do.verif.early==T | do.ens.scores==T){
  early.rows<-which(alldataRows %in% early.years)
  y.neu.pred.verif.scaled<-y.neu.pred.scaled[early.rows]
  if(do.ens.scores==T){
    if(add.arnoise==T){
      if(MCiterations>1){
        #fill directly the entire enseble
        verif.ensemble.rosmx<-t(y.neu.pred.scaled.noise[,(length(reconRows)+1):length(noiseRows)])
        result<-c(result,list(verif.ensemble.rosmx=verif.ensemble.rosmx))
      }else{
        verif.ensx.rosm<-y.neu.pred.scaled.noise[,(length(reconRows)+1):length(noiseRows)]
        result<-c(result,list(verif.ensx.rosm=verif.ensx.rosm))
      }
    }else{
      verif.ensx.rosm<-y.neu.pred.verif.scaled
      result<-c(result,list(verif.ensx.rosm=verif.ensx.rosm))
      
    }
  }
}

y.vec <- target.mean.calib


  }, envir = local.env)
  if(exists('rows', envir = local.env, inherits = FALSE)){
    assign('rows', get('rows', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('double.years', envir = local.env, inherits = FALSE)){
    assign('double.years', get('double.years', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('recon.cont', envir = local.env, inherits = FALSE)){
    assign('recon.cont', get('recon.cont', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('y.vec', envir = local.env, inherits = FALSE)){
    assign('y.vec', get('y.vec', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('y.neu.pred.scaled', envir = local.env, inherits = FALSE)){
    assign('y.neu.pred.scaled', get('y.neu.pred.scaled', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('y.neu.pred.calib.scaled', envir = local.env, inherits = FALSE)){
    assign('y.neu.pred.calib.scaled', get('y.neu.pred.calib.scaled', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('y.neu.pred.scaled.noise', envir = local.env, inherits = FALSE)){
    assign('y.neu.pred.scaled.noise', get('y.neu.pred.scaled.noise', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists('y.neu.pred.verif.scaled', envir = local.env, inherits = FALSE)){
    assign('y.neu.pred.verif.scaled', get('y.neu.pred.verif.scaled', envir = local.env, inherits = FALSE), envir = run.env)
  }
  if(exists("result", envir = local.env, inherits = FALSE)){
    merge.result.payloads(run.env, local.env$result)
  }
  invisible(NULL)
}



