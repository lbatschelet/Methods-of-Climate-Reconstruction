## Function converted from the former sourced template fragment.
## Uses and updates objects stored in run.env.


run_cps_field_reconstruction <- function(run.env){
  local.env <- new.env(parent = run.env)
  evalq({
    result<-list()
recon<-array(dim=c(length(reconRows),dim(target.orig)[2]))
pred.scaled<-array(dim=c(length(alldataRows),dim(target.orig)[2]))
pred.calib.scaled<-array(dim=c(length(calib.years),dim(target.orig)[2]))
if(do.verif.early==T) pred.verif.scaled<-array(dim=c(length(early.years),dim(target.orig)[2]))

overlap.years.his<-which(proxy.years %in% proxy.calib.years)
double.years<-which(alldataRows %in% overlap.years.his)

rows<-as.matrix(dataRows)
if(is.null(dim(proxy.matrix.full))==T) proxy.matrix.full<-as.matrix(proxy.matrix.full)
if(is.null(dim(proxy.matrix.calib))==T) proxy.matrix.calib<-as.matrix(proxy.matrix.calib)

if (add.arnoise==T) pred.scaled.noise<-array(dim=c(MCiterations,length(reconRows),dim(target)[2]))

if(do.ens.scores==T & add.arnoise==T & MCiterations>1){
    verif.ensemblex<-verif.ensemble
}
if(do.ens.scores==T & MCiterations==1){
  verif.ensx<-verif.ensemble[,,1]
  verif.ensx[]<-NA
}

exportvariables.cps<-exportvariables
exportvariables.cps<-c(exportvariables.cps,"rows","overlap.years.his","double.years")

##loop over grid points -----
#for (gp in 1:dim(target.orig)[2]){

results.cpsfield<-foreach (gp = 1:dim(target.orig)[2],.export=exportvariables.cps,.packages = export.packs) %op.dopar.cps.field% {

    result.cpsfield<-list()
    
  norecon<-F
  dataColuse<-dataCol
  if(few.na==T & do.radius==T & minproxies.radius>length(dataCol)){
    norecon<-T
    recon[,gp]<-NA
    if(do.verif.early==T) pred.verif.scaled[,gp]<-NA 
  }
  # select proxies to be used for cps based on radius (only if the minimum acceptable n proxies is smnaller than the # of proxies in the nest)
  if(do.radius==T & minproxies.radius<=length(dataCol)){
    selcps<-which(distances[dataCol,gp]<=search.radius)
    if(length(selcps)<minproxies.radius){
      if(few.na==T){
      #set all outputs to NA because of too few proxies available
        norecon<-T
        recon[,gp]<-NA
        if(do.verif.early==T)  pred.verif.scaled[,gp]<-NA        
      }else{
        selcps<-order(distances[dataCol,gp])[1:minproxies.radius]
      }
    }
    dataColuse<-dataCol[selcps]
  }
  
  if(norecon==F){
    composite<-recon.cps(proxy.matrix.calib,proxy.matrix.full,dataRows,dataColuse,target[,gp],alldataRows,sample.cps.weight,cps.weighting,cps.weight.exponent,minsc.cps,maxsc.cps,calib.years,double.years,cps.weigh.distances,do.cps.cor.weight)
    
    y.neu.pred.scaled<-composite
    
    reconx<-composite[reconRows]
    result.cpsfield<-c(result.cpsfield,list(pred.scaledx=composite))
    result.cpsfield<-c(result.cpsfield,list(pred.calib.scaledx=composite[double.years][calib.years]))
    
#### add AR(1) noise based on the residuals to the recon -----------------
    ##version 1 (simple addition)
    if (add.arnoise==T & add.arnoise.spat==F){
      if(arnoise.version %in% c("simple","gene")==F) stop ("no correct version for arnoise addition provided")
      if(arnoise.version=="simple"){
        crusd<-sd(target[,gp])
        crumean<-mean(target[,gp])
        pred.scaled.noise.gp<-add.residualnoise.index(composite,target[,gp],double.years,calib.years,MCiterations,crusd,crumean,noiseRows)
      }
      if(arnoise.version=="gene"){
        #version 2 (wahl & smerdon 2012)
        crusd<-sd(target[,gp])
        crumean<-mean(target[,gp])
        pred.scaled.noise.gp<-add.residualnoise.gene.cps(composite[double.years][calib.years],target[,gp],MCiterations,crusd,crumean,proxy.matrix.calib,proxy.matrix.full,dataRows,dataColuse,alldataRows,sample.cps.weight,cps.weight.exponent,calib.years,double.years,noiseRows)
      }
      result.cpsfield<-c(result.cpsfield,list(pred.scaled.noisex=pred.scaled.noise.gp))#[,1:length(reconRows)]))
      reconx<-as.vector(pred.scaled.noise.gp[1,1:length(reconRows)])
    }
    
    # early verif years -------------------------------------------------------
    if(do.verif.early==T| do.ens.scores==T){
      early.rows<-which(alldataRows %in% early.years)
      result.cpsfield<-c(result.cpsfield,list(pred.verif.scaledx=y.neu.pred.scaled[early.rows]))
      if(do.ens.scores==T){
        if(add.arnoise==T & add.arnoise.spat==F){
          if(MCiterations>1){
            verif.ensemblexx<-t(pred.scaled.noise.gp[,(length(reconRows)+1):length(noiseRows)])
            result.cpsfield<-c(result.cpsfield,list(verif.ensemblexx=verif.ensemblexx))
          }else{
            verif.ensxx<-pred.scaled.noise.gp[,(length(reconRows)+1):length(noiseRows)]
            result.cpsfield<-c(result.cpsfield,list(verif.ensxx=verif.ensxx))
          }
        }
        if(add.arnoise==F){
          verif.ensxx<-pred.verif.scaled[,gp]
          result.cpsfield<-c(result.cpsfield,list(verif.ensxx=verif.ensxx))
        }
      }
    }
  }

  result.cpsfield<-c(result.cpsfield,list(reconx=reconx))
                     
  result.cpsfield
}
print(paste('combining parallel results for member',ens,'...',date(),Sys.getpid()))

for(gp in 1:dim(target.orig)[2]){
 
  recon[,gp]<- results.cpsfield[[gp]]$reconx
  pred.scaled[,gp]<- results.cpsfield[[gp]]$pred.scaledx
  pred.calib.scaled[,gp]<-results.cpsfield[[gp]]$pred.calib.scaledx
  
  if (add.arnoise==T  & add.arnoise.spat==F){
    pred.scaled.noise[,,gp]<-results.cpsfield[[gp]]$pred.scaled.noisex
  }
  if(do.verif.early){
    pred.verif.scaled[,gp]<-results.cpsfield[[gp]]$pred.verif.scaledx
  }
    if(do.ens.scores==T){
     if(add.arnoise==T & add.arnoise.spat==F & MCiterations>1){
       verif.ensemblex[,gp,]<-results.cpsfield[[gp]]$verif.ensemblexx
     }
     if((add.arnoise==T & add.arnoise.spat==F & MCiterations==1) | add.arnoise==F){
       verif.ensx[,gp]<-results.cpsfield[[gp]]$verif.ensxx
     }
    }
}

if("verif.ensemblexx" %in% names(results.cpsfield[[1]]))  result<-c(result,list(verif.ensemblex=verif.ensemblex))
if("verif.ensxx" %in% names(results.cpsfield[[1]]))  result<-c(result,list(verif.ensx=verif.ensx))

if(add.arnoise==T & add.arnoise.spat==T){
  crusd<-apply(target,2,sd)
  crumean<-apply(target,2,mean)
  pred.scaled.noise<-add.residualnoise.cov(pred.scaled,target,double.years,calib.years,MCiterations,crusd,crumean,noiseRows,do.pc.arnoise.spat,latweights)
  if(MCiterations==1){
    recon<-array(pred.scaled.noise[,1:length(reconRows),],dim=c(length(reconRows),dim(pred.scaled.noise)[3])) 
  }else{
    recon<-pred.scaled.noise[1,1:length(reconRows),]
  }
}

if(do.ens.scores==T & add.arnoise==T & add.arnoise.spat==T){
  if(MCiterations>1){
    #fill directly the entire enseble
    verif.ensemblex<-aperm(pred.scaled.noise[,(length(reconRows)+1):length(noiseRows),],c(2,3,1))
    result<-c(result,list(verif.ensemblex=verif.ensemblex))
  }else{
    verif.ensx<-array(pred.scaled.noise[,(length(reconRows)+1):length(noiseRows),],dim=c(length(early.rows),dim(pred.scaled.noise)[3]))
    result<-c(result,list(verif.ensx=verif.ensx))
  }
}


rm(results.cpsfield)

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



