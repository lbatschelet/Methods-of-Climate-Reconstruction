
library(MASS)
#library(mvnfast)
library(mgcv)

##add residualnoise based on residual covariance. translated from Jianghaos script, 2016/07/18, RN
add.residualnoise.cov<-function(pred.scaled,target,double.years,calib.years,MCiterations,crusd,crumean,reconRows,do.pc.arnoise.spat,latweights){

  n<-length(reconRows)
  residuals<-target-pred.scaled[double.years,][calib.years,]
  
  res.ac<-apply(residuals,2,fastacf)
  
  if(do.pc.arnoise.spat){
    resid.latw<-residuals*rep(latweights,each=dim(residuals)[1])
    resid.pc<-f.pc(resid.latw,0.9)
    covRes.pc<-cov2(resid.pc$h)
  }else{
      covRes<-cov2(residuals)
  }


  pred.scaled.noise<-foreach (s = 1:MCiterations, .combine=function(...) abind(..., along=0),.multicombine = T) %op.dopar.addnoise% {
  if(do.pc.arnoise.spat){
    yerr.pc<-mgcv::rmvn(n=n, mu=rep(0,dim(resid.pc$h)[2]),V = covRes.pc)
    yerr<-yerr.pc %*% t(resid.pc$a)
  } else{
    yerr<-mgcv::rmvn(n=n, mu=rep(0,dim(residuals)[2]),V = covRes)
  } 
  if(n>5){  #only make the ar stuff if the ts is longer than 5 years
    for(i in 1:dim(residuals)[2]){
      noise<-yerr[1,i]
      for(t in 2:n){
        noise[t]<-res.ac[i]*noise[t-1]+yerr[t,i]
      }
      sdf<-sqrt(var(residuals[,i]))/sqrt(var(noise))
      noise<-noise*sdf
      mf<-mean(residuals[,i])-mean(noise)
      yerr[,i]<-noise+mf
    }
  }    
    pred.scaled.noise<-pred.scaled[reconRows,]+yerr
    out<-pred.scaled.noise
    if(length(reconRows)==1) out<-t(out)
    out
    
  }
  if(length(dim(pred.scaled.noise))==2) pred.scaled.noise<-array(pred.scaled.noise,dim=c(1,length(reconRows),dim(target)[2]))
  
  pred.scaled.noise
}
