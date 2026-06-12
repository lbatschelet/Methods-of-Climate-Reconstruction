## Function converted from the former sourced template fragment.
## Uses and updates objects stored in run.env.


run_pcr_reconstruction <- function(run.env){
  local.env <- new.env(parent = run.env)
  evalq({
    result<-list()
rows<-as.matrix(dataRows)
if(is.null(dim(proxy.matrix.full))==T) proxy.matrix.full<-as.matrix(proxy.matrix.full)

# scaling and pcs of proxy data -------------------------------------------
#### scale over all overlap years of proxies but take pcs only over calib period

 overlap.years.his<-which(proxy.years %in% proxy.calib.years)
 double.years<-which(alldataRows %in% overlap.years.his)
 x.all<-proxy.matrix.full[alldataRows,dataCol]
 if(is.null(dim(x.all))==T) x.all<-as.matrix(x.all)
  if(sample.weights==F){
  x.scaled<-scale(x.all)
 }else{
  x.scaled<-scale(x.all,scale=scalefactors*apply(x.all,2,sd))
 }
if(length(dataCol)>1){
if(do.pc.opt==T){
  pc.x.all<-f.pc.opt(x.scaled[double.years,][calib.years,],pc.opt,pc.sample)
  north.npc.all.sx<-pc.x.all$npc
}else{
  if(npc.fix==T){
    pc.x.all<-f.pc.nfixed(x.scaled[double.years,][calib.years,],npc.s)  
  }else{
    pc.x.all<-f.pc(x.scaled[double.years,][calib.years,], proz.pc.s)
    if(sample.pcs){
      if(is.null(dim(pc.x.all$h))){
        north.npc.all.sx<-1
      }else{
        north.npc.all.sx<-dim(pc.x.all$h)[2]
      }
    }
  }
}

u.mat<-as.matrix(pc.x.all$h)
u.neu<-x.scaled %*% pc.x.all$a

}else{ #in case there is only 1 proxy used
  u.mat<-x.scaled[double.years,][calib.years]
  u.neu<-x.scaled
  north.npc.all.sx<-NA
}

u1.mat <- cbind(1, u.mat)
u1.neu <- cbind(1, u.neu)

###field recon -----------
if (do.field==T){

#  Principal component regression: ----------------------------------------
b1.mat <- solve(crossprod(u1.mat)) %*% t(u1.mat) %*% pc.target.h[calib.years,]

# Praediktion an der Stelle x.neu
y.neu.pred1 <- u1.neu %*% b1.mat
##back-transformation of PCs
y.neu.pred <-y.neu.pred1 %*% t(pc.target.a)


# Scale results -----------------------------------------------------------
if(do.var.adj==T){
  sdf.c<-apply(y.neu.pred[double.years,][calib.years,],2,sd)/apply(target,2,sd)
  pred.scaled<-scale(y.neu.pred,scale=sdf.c,center=F)
  mf.c<-apply(pred.scaled[double.years,][calib.years,],2,mean)-apply(target,2,mean)
  
  pred.scaled<-scale(pred.scaled,center=mf.c,scale=F)
}else{
  pred.scaled <-y.neu.pred/rep(latweights,each=dim(y.neu.pred)[1])
  #pred.scaled<-y.neu.pred
  pred.scaled <-pred.scaled*rep(instr.scale.factor,each=dim(y.neu.pred)[1])
  pred.scaled<-pred.scaled+rep(instr.center.factor,each=dim(y.neu.pred)[1])
}
  
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
      crusd<-apply(target,2,sd)
      crumean<-apply(target,2,mean)
      latwf<-rep(latweights,each=dim(y.neu.pred)[1])
      pred.scaled.noise<-add.residualnoise.gene(y.neu.pred1[double.years,][calib.years,],pc.target.h[calib.years,],MCiterations,u1.mat,u1.neu,pc.target.a,crusd,crumean,calib.years,double.years,noiseRows,instr.scale.factor,instr.center.factor,latwf,do.var.adj)
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
  pred.verif.scaled<-pred.scaled[early.rows,]
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

###end field recon ------------
}

# spatial mean recon ------------------------------------

if(do.index==T){

  y.vec <- target.mean.calib
b1.vec <- solve(crossprod(u1.mat)) %*% t(u1.mat) %*% y.vec
y.neu.pred <- as.vector(u1.neu %*% b1.vec)

#scaling
if(do.var.adj==T){
  sdf.c<-sd(y.neu.pred[double.years][calib.years])/sd(y.vec)
  y.neu.pred.scaled<-scale(y.neu.pred,scale=sdf.c,center=F)
  mf.c<-mean(y.neu.pred.scaled[double.years][calib.years])-mean(y.vec)
  
  y.neu.pred.scaled<-scale(y.neu.pred.scaled,center=mf.c,scale=F)
}else{
    y.neu.pred.scaled<-y.neu.pred
}

recon.cont<-y.neu.pred.scaled[reconRows]

y.neu.pred.calib.scaled<-y.neu.pred.scaled[double.years][calib.years]

#### add AR(1) noise based on the residuals to the recon -----------------
##version 1 (simple addition)
if (add.arnoise==T){
  if(arnoise.version %in% c("simple","gene")==F) stop ("no correct version for arnoise addition provided")
  if(arnoise.version=="simple"){
    crusd<-sd(target.mean.calib)
    crumean<-mean(target.mean.calib)
    y.neu.pred.scaled.noise<-add.residualnoise.index(y.neu.pred.scaled,target.mean.calib,double.years,calib.years,MCiterations,crusd,crumean,noiseRows)
  }else{
    #version 2 (wahl & smerdon 2012)
    crusd<-sd(target.mean.calib)
    crumean<-mean(target.mean.calib)
    y.neu.pred.scaled.noise<-add.residualnoise.gene.index(y.neu.pred[double.years][calib.years],target.mean.calib,MCiterations,u1.mat,u1.neu,crusd,crumean,noiseRows)
  }
  recon.cont<-as.vector(y.neu.pred.scaled.noise[1,1:length(reconRows)])
}

}

#early verif years -------------
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

if(do.pc.opt==T | sample.pcs==T){
result<-c(result,list(north.npc.all.sx=north.npc.all.sx))
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




