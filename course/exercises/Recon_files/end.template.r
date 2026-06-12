## Function converted from the former sourced template fragment.
## Returns the updated result list built from objects in run.env.


pack_results <- function(run.env){
  local.env <- new.env(parent = run.env)
  result <- evalq({
# first the field recon stuff--------------

if(do.field==T){

#recdata[which(is.na(recdata))]<-0
result<-c(result,list(pred.scaled=pred.scaled))#ensemble.sum<-ensemble.sum+pred.scaled

result<-c(result,list(recon=recon))#ensemble[ens,,]<-recon#[reconRows,]

}

###The spatial mean recons ---------------------
if(do.index==T){
  result<-c(result,list(y.neu.pred.scaled=y.neu.pred.scaled))#ensemble.sum.rosm<-ensemble.sum.rosm+y.neu.pred.scaled
  result<-c(result,list(recon.cont=recon.cont))#all.sigma[models$years[[model]],ens]<-recon.cont#[reconRows]
}

#  now the verif stuff ---------------------
############################################
if(do.field==T & do.verif==T){
if(do.proz.pos==T){
proz.pos<-apply(res,1,function(x) length(which(x>0))/dim(res)[2]*100)
result<-c(result,list(proz.pos=proz.pos))#all.proz.pos[models$years[[model]],ens]<-proz.pos
}
 
  if(do.re==T){
#resdata[which(is.na(resdata))]<-0
  result<-c(result,list(res=res))#res.sum<-res.sum+res

#if(write.quantiles.verif==T) res.ensemble[ens,,]<-res
}



# CE
if(do.ce==T){
  
  #cesdata[which(is.na(cesdata))]<-0
  result<-c(result,list(ces=ces))#ces.sum<-ces.sum+ces
  
  #if(write.quantiles.verif==T) ces.ensemble[ens,,]<-ces
}

# r2
if(do.r2==T){

#r2data[which(is.na(r2data))]<-0
  result<-c(result,list(r2=r2))#r2.sum<-r2.sum+r2

#if(write.quantiles.verif==T) r2.ensemble[ens,,]<-r2
}


}

if(do.index==T & do.verif==T){
  if(do.ce.index==T){
    result<-c(result,list(ce.cont=ce.cont))#all.ce.cont[models$years[[model]],ens]<-ce.cont
  }
  
  
  if(do.re.index==T){
    result<-c(result,list(re.cont=re.cont))#all.re.cont[models$years[[model]],ens]<-re.cont
  }
  
  
  if(do.r2.index==T){
    result<-c(result,list(r2.cont=r2.cont))#all.r2.cont[models$years[[model]],ens]<-r2.cont
  }

}
    result
  }, envir = local.env)
  result
}



