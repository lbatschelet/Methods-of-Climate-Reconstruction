### create ensemble of reconstructions by adding noise 
### noise calculation and addition adapted from Wahl & Smerdon 2012

library(waveslim)


add.residualnoise.gene<-function(TmatU.ret.anom.cal.fitted,crupccalib,MCiterations,u1.mat,u1.neu,pc.target.a,crusd,crumean,calib.years,double.years,reconRows,instr.scale.factor,instr.center.factor,latwf,do.var.adj){
  
  #####CALCULATION OF RPC RESIDUALS and MONTE CARLO CREATION OF NEW RESIDUALS and CREATION OF NEW INSTRUMENTAL PC SERIES BASED ON MC RESIDUALS
  
  TmatUresiduals<-crupccalib-TmatU.ret.anom.cal.fitted
  
  N<-nrow(TmatUresiduals)
  n<-ncol(TmatUresiduals)
  
  #DETERMINATION OF ACFs FOR RPC SERIES
  Data1<-array(rep(NA,N*n), dim=c(N,n))
  for (l in 1:ncol(TmatUresiduals)){
    Data1[,l]<-acf(TmatUresiduals[,l][!is.na(TmatUresiduals[,l])],N,plot=F)[[1]][1:N]}
  
  #CALCULATION OF SCALING FACTORS TO BE USED IN "NM" LOOP
  std_network<-apply(TmatUresiduals,2,sd)
  std_networkmatrix<-array(dim=c(N,n))
  for (o in 1:N){
    std_networkmatrix[o,]<-std_network}
  
  #RED NOISE SERIES SIMULATION LOOP
  TmatUresiduals_sim<-array (rep(NA,MCiterations*N*n), dim=c(N,MCiterations,n) )
  for (q in 1:MCiterations) {
    for (g in 1:n) {
      TmatUresiduals_sim[,q,g]<-hosking.sim(N,Data1[,g]) } # end "g" loop
    
    TmatUresiduals_sim[,q,]<-std_networkmatrix*TmatUresiduals_sim[,q,] } # end "q" loop
  pred.scaled.noisens<-array (dim=c(MCiterations,length(reconRows),dim(pc.target.a)[1]) )
  
  pred.scaled.noisens<-foreach (w = 1:MCiterations, .combine=function(...) abind(..., along=0),.multicombine = T) %op.dopar.addnoise% {
  #for (w in 1:MCiterations) {

    b1.mat <- solve(crossprod(u1.mat)) %*% t(u1.mat) %*% (TmatU.ret.anom.cal.fitted+TmatUresiduals_sim[,w,])
    y.neu.pred1<- u1.neu %*% b1.mat
    y.neu.pred <-y.neu.pred1 %*% t(pc.target.a)
    
    # Scale results -----------------------------------------------------------
   if(do.var.adj==T){
      sdf.c<-apply(y.neu.pred[double.years,][calib.years,],2,sd)/crusd
      pred.scaled<-scale(y.neu.pred,scale=sdf.c,center=F)
      mf.c<-apply(pred.scaled[double.years,][calib.years,],2,mean)-crumean
      pred.scaled<-scale(pred.scaled,center=mf.c,scale=F)
    }else{
       #pred.scaled<-y.neu.pred
       pred.scaled <-y.neu.pred/latwf
       pred.scaled <-pred.scaled*rep(instr.scale.factor,each=dim(y.neu.pred)[1])
       pred.scaled<-pred.scaled+rep(instr.center.factor,each=dim(y.neu.pred)[1])
     }
    if(w%%10==0) print (paste("AR noise iteration",w))
    out<-pred.scaled[reconRows,]
    if(length(reconRows)==1) out<-t(out)
    out
  }
  if(length(dim(pred.scaled.noisens))==2) pred.scaled.noisens<-array(pred.scaled.noisens,dim=c(1,dim(pred.scaled.noisens)[1],dim(target)[2]))
  
  pred.scaled.noisens
}


################
########## now for index recon --------------


add.residualnoise.gene.index<-function(TmatU.ret.anom.cal.fitted,crupccalib,MCiterations,u1.mat,u1.neu,crusd,crumean,reconRows){
  
  sds.TmatU.ret.anom.cal.fitted.vector<-sqrt(var(TmatU.ret.anom.cal.fitted))
  sdfp<-sds.TmatU.ret.anom.cal.fitted.vector/sqrt(var(crupccalib))
  TmatU.ret.anom.cal.fitted<-scale(TmatU.ret.anom.cal.fitted,scale=sdfp,center=F)
  
  TmatUresiduals<-crupccalib-TmatU.ret.anom.cal.fitted
  
  N<-nrow(TmatUresiduals)
  
  #DETERMINATION OF ACFs FOR RPC SERIES
    Data1<-acf(TmatUresiduals[!is.na(TmatUresiduals)],N,plot=F)[[1]][1:N]
  
  #CALCULATION OF SCALING FACTORS TO BE USED IN "NM" LOOP
  std_network<-sd(TmatUresiduals)
  
  #RED NOISE SERIES SIMULATION LOOP
  TmatUresiduals_sim<-array (rep(NA,MCiterations*N), dim=c(N,MCiterations) )
  for (q in 1:MCiterations) {
      TmatUresiduals_sim[,q]<-hosking.sim(N,Data1)
    
    TmatUresiduals_sim[,q]<-std_network*TmatUresiduals_sim[,q] } # end "q" loop
  pred.scaled.noisens<-array (dim=c(MCiterations,length(reconRows)))

  for (w in 1:MCiterations) {
    
    
    b1.vec <- solve(crossprod(u1.mat)) %*% t(u1.mat) %*% (TmatU.ret.anom.cal.fitted+TmatUresiduals_sim[,w])
    y.neu.pred<- u1.neu %*% b1.vec
    pred.scaled<-y.neu.pred
    pred.scaled.noisens[w,]<-pred.scaled[reconRows]
    
    
  }
  pred.scaled.noisens
}



################
########## now for CPS --------------

add.residualnoise.gene.cps<-function(TmatU.ret.anom.cal.fitted,crupccalib,MCiterations,crusd,crumean,proxy.matrix.calib,proxy.matrix.full,dataRows,dataCol,alldataRows,sample.cps.weight,cps.weight.exponent,calib.years,double.years,reconRows){
    
  sds.TmatU.ret.anom.cal.fitted.vector<-sqrt(var(TmatU.ret.anom.cal.fitted))
  sdfp<-sds.TmatU.ret.anom.cal.fitted.vector/sqrt(var(crupccalib))
  TmatU.ret.anom.cal.fitted<-scale(TmatU.ret.anom.cal.fitted,scale=sdfp,center=F)
  
  TmatUresiduals<-crupccalib-TmatU.ret.anom.cal.fitted
  
  N<-nrow(TmatUresiduals)
  
  #DETERMINATION OF ACFs FOR RPC SERIES
  Data1<-acf(TmatUresiduals[!is.na(TmatUresiduals)],N,plot=F)[[1]][1:N]
  
  #CALCULATION OF SCALING FACTORS TO BE USED IN "NM" LOOP
  std_network<-sd(TmatUresiduals)
  
  #RED NOISE SERIES SIMULATION LOOP
  TmatUresiduals_sim<-array (rep(NA,MCiterations*N), dim=c(N,MCiterations) )
  for (q in 1:MCiterations) {
    TmatUresiduals_sim[,q]<-hosking.sim(N,Data1)
    TmatUresiduals_sim[,q]<-std_network*TmatUresiduals_sim[,q] } # end "q" loop
    pred.scaled.noisens<-array (dim=c(MCiterations,length(reconRows)) )
    for (w in 1:MCiterations) {
      
      y.neu.pred<-recon.cps(proxy.matrix.calib,proxy.matrix.full,dataRows,dataCol,(TmatU.ret.anom.cal.fitted+TmatUresiduals_sim[,w]),alldataRows,sample.cps.weight,cps.weighting,cps.weight.exponent,minsc.cps,maxsc.cps,calib.years,double.years)
          
       # Scale results -----------------------------------------------------------
      
      sdf.c<-sd(y.neu.pred[double.years][calib.years])/crusd
      pred.scaled<-scale(y.neu.pred,scale=sdf.c,center=F)
      mf.c<-mean(pred.scaled[double.years][calib.years])-crumean
      pred.scaled<-scale(pred.scaled,center=mf.c,scale=F)
      
      pred.scaled.noisens[w,]<-pred.scaled[reconRows]
      
    
    }
  pred.scaled.noisens
}










################
########## now for CCA --------------
add.residualnoise.gene.cca<-function(TmatU.ret.anom.cal.fitted,crupccalib,MCiterations,s.calib,x.all,calib.years,double.years,reconRows,dp,dt,dcca){
  TmatUresiduals<-crupccalib-TmatU.ret.anom.cal.fitted
  
  N<-nrow(TmatUresiduals)
  n<-ncol(TmatUresiduals)
  
  #DETERMINATION OF ACFs FOR RPC SERIES
  Data1<-array(rep(NA,N*n), dim=c(N,n))
  for (l in 1:ncol(TmatUresiduals)){
    Data1[,l]<-acf(TmatUresiduals[,l][!is.na(TmatUresiduals[,l])],N,plot=F)[[1]][1:N]} # end "l" loop
  
  #CALCULATION OF SCALING FACTORS TO BE USED IN "NM" LOOP
  std_network<-apply(TmatUresiduals,2,sd)
  std_networkmatrix<-array(dim=c(N,n))
  for (o in 1:N){
    std_networkmatrix[o,]<-std_network}
  
  #RED NOISE SERIES SIMULATION LOOP
  TmatUresiduals_sim<-array (rep(NA,MCiterations*N*n), dim=c(N,MCiterations,n) )
  for (q in 1:MCiterations) {
    for (g in 1:n) {
      TmatUresiduals_sim[,q,g]<-hosking.sim(N,Data1[,g]) } # end "g" loop
      TmatUresiduals_sim[,q,]<-std_networkmatrix*TmatUresiduals_sim[,q,] } # end "q" loop
      pred.scaled.noisens<-array (dim=c(MCiterations,length(reconRows),dim(pc.target.a)[1]) )
      pred.scaled.noisens<-foreach (w = 1:MCiterations, .export = c("cca_bp","svd.matlab","repmat"), .combine=function(...) abind(..., along=0),.multicombine = T) %op.dopar.addnoise% {
      
      TT.resid<-(TmatU.ret.anom.cal.fitted+TmatUresiduals_sim[,w,])
      cca_bp_out.resid<-cca_bp(s.calib,TT.resid,dp,dt,dcca);
      B.resid<-cca_bp_out.resid$B
      Tm.resid<-cca_bp_out.resid$Tm
      Pm.resid<-cca_bp_out.resid$Pm
      Ts.resid<-cca_bp_out.resid$Ts
      Ps.resid<-cca_bp_out.resid$Ps
      
      field_rc.resid = (x.all - repmat(Pm.resid,dim(x.all)[1],1))%*%diag(1.0/Ps.resid)%*%t(B.resid)%*%diag(Ts.resid)
      pred.scaled.resid  = field_rc.resid  + repmat(Tm.resid,dim(x.all)[1],1)
      
      out<-pred.scaled.resid[reconRows,]
      if(length(reconRows)==1) out<-t(out)
      
      if(w%%10==0) print (paste("AR noise iteration",w))
      out
    }
  if(length(dim(pred.scaled.noisens))==2) pred.scaled.noisens<-array(pred.scaled.noisens,dim=c(1,dim(pred.scaled.noisens)[1],dim(target)[2]))
  
  pred.scaled.noisens
}

