##cps recon function

recon.cps<-function(proxy.matrix.calib,proxy.matrix.full,dataRows,dataCol,target,alldataRows,sample.cps.weight,cps.weighting,cps.weight.exponent,minsc.cps,maxsc.cps,calib.years,double.years,cps.weigh.distances=F,do.cps.cor.weight=T){

  #correlate proxies with target
  if(do.cps.cor.weight==T){
      cors<- cor(proxy.matrix.calib[,dataCol],target,use="pairwise.complete.obs")
  }else{
    cors<-rep(1,length(dataCol))
  }

  # adjust sign -------------------------------------------------------------
  proxy.matrix.full.pos<-proxy.matrix.full[,dataCol]
  if(length(dataCol)>1){
    for( i in 1:dim(proxy.matrix.full.pos)[2]){
      proxy.matrix.full.pos[,i]<-proxy.matrix.full.pos[,i]*sign(cors[i])
    }
  }else{
    proxy.matrix.full.pos<-proxy.matrix.full.pos*sign(cors) 
  }
    
  # make the composite -----------------------------------------------------
  if(length(dataCol)>1){
    x.all<-proxy.matrix.full.pos[alldataRows,]
    x.scaled<-scale(x.all)
    #define weighting factor
    if(sample.cps.weight==F){  #case with nod additional random sampling factor: weight the proxies by their target correlation
      weights.cps<-abs(cors)
      x.scaled.w<-x.scaled
      for(i in 1:dim(x.scaled)[2]){
        x.scaled.w[,i]<-x.scaled[,i]*weights.cps[i]
      }
      if(cps.weigh.distances==T){  #weight by distance
        x.scaled.w<-x.scaled.w %*% (dpn[dataCol,dataCol]/max(dpn[dataCol,dataCol]))
      }
      composite<-apply(x.scaled.w,1,mean,na.rm=T)
    }else{
      x.scaled.w<-x.scaled
      if(cps.weighting=="c"){  #weigh by exponent as in cook 2010
        weights.cps<-abs(cors)^cps.weight.exponent
        for(i in 1:dim(x.scaled)[2]){
          x.scaled.w[,i]<-abs(x.scaled[,i])^weights.cps[i]*sign(x.scaled[,i])
        }
      }else{  #weigh by factor between x and 1/x 
        weights.cps<-abs(cors)*runif(dim(x.scaled)[2],min=minsc.cps,max=maxsc.cps)
        for(i in 1:dim(x.scaled)[2]){
          x.scaled.w[,i]<-x.scaled[,i]*weights.cps[i]
        }
      }
      if(cps.weigh.distances==T){
        x.scaled.w<-x.scaled.w %*% dpn[dataCol,dataCol]
      }
      composite<-apply(x.scaled.w,1,mean,na.rm=T)
    }
  }else{   #case with only one record
    composite<-proxy.matrix.full.pos[alldataRows]
  }
   
  # #scale the composite ----------------------------------------------------
  s<-sd(composite[double.years][calib.years])/sd(target)
  composite.scaled<-composite/s
  c<-mean(composite.scaled[double.years][calib.years])-mean(target)
  composite.scaled<-composite.scaled-c
  return(composite.scaled)
}
  