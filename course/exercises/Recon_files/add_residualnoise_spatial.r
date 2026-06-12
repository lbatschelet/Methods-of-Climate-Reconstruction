add.residualnoise.spat<-function(pred.scaled,target,double.years,calib.years,MCiterations,crusd,crumean,reconRows){
  N<-dim(pred.scaled)[1]
  residuals<-target-pred.scaled[double.years,][calib.years,]
  
  lat.weight.factor<-(cos(lats.2d*pi/180)^latweight.cos.power)
  
  resid.latw<-residuals*rep(lat.weight.factor,each=dim(residuals)[1])#lat.weight.factor
  resid.pc<-f.pc(resid.latw,0.9)
  
  rhos<-apply(resid.pc$h,2,fastacf)
  vars<-apply(resid.pc$h,2,var)
  
  pred.scaled.noise<-foreach (s = 1:MCiterations, .export = "create.ar1.var.ts", .combine=function(...) abind(..., along=0),.multicombine = T) %op.dopar.addnoise% {
  
    pcs.noisy<-array(dim=c(N,dim(resid.pc$h)[2]))
    for(i in 1:dim(resid.pc$a)[2]){
      pcs.noisy[,i]<-create.ar1.var.ts(N,rho = rhos[i],varn = vars[i]) 
    }
    resid.back.noisy<-pcs.noisy %*% t(resid.pc$a)
    resid.back.noisy.nw<-resid.back.noisy/rep(lat.weight.factor,each=dim(resid.back.noisy)[1])#lat.weight.factor
    resid.back.noisy.scaled<-resid.back.noisy.nw
    
    for(gp in 1:dim(target)[2]){
      sdf<-sd(resid.back.noisy.nw[double.years[calib.years],gp])/sd(residuals[,gp])
      resid.back.noisy.scaled[,gp]<-resid.back.noisy.nw[,gp]/sdf
      mf<-mean(resid.back.noisy.scaled[double.years[calib.years],gp])-mean(residuals[,gp])
      resid.back.noisy.scaled[,gp]<-resid.back.noisy.scaled[,gp]-mf
    }
    
  pred.scaled.noise<-pred.scaled[reconRows,]+resid.back.noisy.scaled[reconRows,]
  out<-pred.scaled.noise
  if(length(reconRows)==1) out<-t(out)
  out
    
  }
  if(length(dim(pred.scaled.noise))==2) pred.scaled.noise<-array(pred.scaled.noise,dim=c(1,length(reconRows),dim(target)[2]))
    
  pred.scaled.noise
}

