### create ensemble of reconstructions by adding noise to the reconstructed time series.
## noise has the mean, sd and AR1 coefficient of the calibration residuals

add.residualnoise<-function(y.neu.pred,target,double.years,calib.years,MCiterations,crusd,crumean,reconRows){
  N<-dim(y.neu.pred)[1]
  residuals<-target-y.neu.pred[double.years,][calib.years,]
  sd.res<-sqrt(apply(residuals,2,var))
  ars.int<-apply(residuals,2,function(x) cor(x,c(x[-1],NA),use="complete.obs"))

  pred.scaled.noise<-foreach (s = 1:MCiterations, .combine=function(...) abind(..., along=0),.multicombine = T)%op.dopar.addnoise% {
  #for (s in 1:MCiterations){
    y.neu.pred.noise<-array (dim=c(dim(y.neu.pred)[1],dim(y.neu.pred)[2]) )
    for(i in 1:dim(y.neu.pred)[2]){
      noise<-arima.sim(n=N,list(ar=ars.int[i]))
      noise<-noise/sqrt(var(noise))*sd.res[i]
      noise<-noise-mean(noise)
      y.neu.pred.noise[,i]<-y.neu.pred[,i]+noise
    }
    out<-y.neu.pred.noise[reconRows,]
    if(length(reconRows)==1) out<-t(out)
    out
  }
  pred.scaled.noise
}



## now for an index recon -----------------------

add.residualnoise.index<-function(y.neu.pred,target,double.years,calib.years,MCiterations,crusd,crumean,reconRows){
  N<-length(y.neu.pred)
  residuals<-target-y.neu.pred[double.years][calib.years]
  sd.res<-sqrt(var(residuals))
  ars.int<-cor(residuals,c(residuals[-1],NA),use="complete.obs")
  pred.scaled.noise<-array (dim=c(MCiterations,length(reconRows)) )
  for (s in 1:MCiterations){
    noise<-arima.sim(n=N,list(ar=ars.int))
    noise<-noise/sqrt(var(noise))*sd.res
    noise<-noise-mean(noise) #+mean(residuals)
    y.neu.pred.noise<-y.neu.pred+noise 
    pred.scaled.noise[s,]<-y.neu.pred.noise[reconRows]
  }
  pred.scaled.noise
}
