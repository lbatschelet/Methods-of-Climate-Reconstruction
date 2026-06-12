cca_cv<-function(temp,prox,cca_options){
##cca_CV translated from matlab function 2017/05/19 RN [dp_opt,dt_opt,dcca_opt,misfit] = cca_cv(temp,prox,cca_options)
# Choose an optimal set of CCA truncation parameters.
#
# Three truncation parameters are required for CCA reconstruction:
  #     dp:   Truncation parameter of the proxy matrix
#     dt:   Truncation parameter of the temperature matrix
#     dcca: Truncation parameter of the covariance matrix of proxy and temperature
#
#
#     dcca by design is smaller than min(dp,dt). The optimal selection of
#     dt, dp and dcca is based on the k-fold cross-validation RMSE.
#
#     The CCA reconstruction method is inspired by Smerdon et al. [2010a]
#
#
# Input:
  #     temp: Temperature matrix
#     prox: Proxy (pseudoproxy) matrix
#           Notice that temperature and proxy matrix should have the same
#           time dimension
#        k: Number of folds of K-fold cross-validation
#     indices: Indices of K-fold cross-validation
#     noise:  Index of noise level, used to save RMSE output
#
# Output:
  #     dp_opt:   Optimal choice of dp
#     dt_opt:   Optimal choice of dt
#     dcca_opt: Optimal chocie of dcca
#     misfit:   A struct containing the useful error information as follows:
  #     rmse:  RMSE of different folds
#     error: Truncated version of rmse
#
# See also TRUNCATION_CRITERIA, CCA_BP

  svd.matlab<-function(x){
    svdx<-svd(x,nu=max(nrow(x),ncol(x)),nv=max(nrow(x),ncol(x)))
    U<-svdx$u
    S<-x
    S[]<-0
    if(length(svdx$d)>1){
      s1<-diag(svdx$d)
      S[1:dim(s1)[1],1:dim(s1)[2]]<-s1
    }else{
      S[1]<-svdx$d
    }
    V<-svdx$v
    list(U=U,S=S,V=V)
  }

  
fopts    = names(cca_options)
K        = cca_options$K
indices  = cca_options$indices
weights  = cca_options$weights
# noise    = cca_options$noise
dt_max   = cca_options$dt_max
dp_max   = cca_options$dp_max
method   = cca_options$method
model   = cca_options$model
posfile= cca_options$posfile


dp_max<-min(cca_options$dp_max,length(which(indices==1)))
dp_range   = 1:dp_max

dt_max<-min(cca_options$dt_max,length(which(indices==1)))
dt_range   = 1:dt_max


#if (method=='smerdon10'){
  rmse = array(dim=c(min(dp_max,dt_max),dp_max,dt_max,K))
  for (k in 1:K){
    print(paste0('Fold-',k))
   # Define training and test data
    test <- (indices == k); train = (indices != k)
    ntest = length(test)
    
   # Temperature matrix
    temp_train    = temp[train,]
    temp_test     = temp[test,]
    Tds<- scale(temp_train)
    Tm<-apply(temp_train,2,mean)
    Ts<-apply(temp_train,2,sd)
   # Proxy matrix
    prox_train    = prox[train,]
    Pds<- scale(prox_train)
    Pm<-apply(prox_train,2,mean)
    Ps<-apply(prox_train,2,sd)

    svdp<-svd.matlab(t(Pds))
    UP<-svdp$U
    SP<-svdp$S #diag(svdp$d)
    VP<-svdp$V
    
    svdt<-svd.matlab(t(Tds))
    UT<-svdt$U
    ST<-svdt$S
    VT<-svdt$V
    
        # Cross Covariance matrix
    #for (dt in dt_range){
    results.ccapar<-foreach (dt = dt_range) %dopar% {
      n<-1
      rmses<-c()
      poss<-c()
      for (dp in dp_range){
        write.table(paste0('----nest ',model,' dp = ',dp,' dt = ',dt,' k = ',k, ', dp_max = ', (max(dp_range)),date(),'----'),posfile,append = T,row.names = FALSE,col.names = FALSE)
        dcca_max = min(dt,dp,cca_options$dcca_max)
        for (dcca in 1:dcca_max){
          FX          = t(VT[,1:dt])%*%VP[,1:dp]
          svdf<-svd.matlab(FX)
          UF<-svdf$U
          SF<-svdf$S
          VF<-svdf$V
        
         # Estimate CV regression matrix
          if(dt==1){
            utx<-t(t(UT[,1:dt]))
          }else{
            utx<-UT[,1:dt]
          }
          if(dp==1){
            B_cv       = utx%*%ST[1:dt,1:dt]%*%UF[1:dt,1:dcca]%*%SF[1:dcca,1:dcca]%*%t(VF[1:dp,1:dcca])%*%diag(1.0/diag(SP[1:dp,1:dp],nrow=length(SP[1:dp,1:dp])))%*%t(UP[,1:dp])
          }else{
            B_cv       = utx%*%ST[1:dt,1:dt]%*%UF[1:dt,1:dcca]%*%SF[1:dcca,1:dcca]%*%t(VF[1:dp,1:dcca])%*%diag(1.0/diag(SP[1:dp,1:dp]))%*%t(UP[,1:dp])
          }
          
         # Predict temperature over verification sample
          temp_pred  = (prox-t(t(rep(1,ntest)))%*%Pm)%*%diag(1.0/Ps)%*%t(B_cv)%*%diag(Ts)+ t(t(rep(1,ntest)))%*%Tm
          
         # Calculate the misfit (RMSE)
          mse                = apply((temp_pred[test,] - temp_test)^2,2,mean)
          rmses[n] = sqrt(mse%*%t(t(weights)))
          poss[n]=(k-1)*dim(rmse)[1]*dim(rmse)[2]*dim(rmse)[3]+(dt-1)*dim(rmse)[1]*dim(rmse)[2]+(dp-1)*dim(rmse)[1]+dcca
          n<-n+1
        }
      }
      list(rmses=rmses,poss=poss)
    }
  #  stopCluster(cl)
    
    for(i in seq_along(results.ccapar)){
       rmse[results.ccapar[[i]]$poss]<-results.ccapar[[i]]$rmses
    } 
  }
  
  rm(results.ccapar)
  
  mean_misfit = sqrt(apply(rmse^2,c(1,2,3),mean,na.rm=T))
  min_RMSE = min(mean_misfit,na.rm=T)
  bestopts = which(mean_misfit==min_RMSE,arr.ind=T)
  dcca_opt<-bestopts[1]
  dp_opt<-bestopts[2]
  dt_opt<-bestopts[3]

print('Computation is done!')

list(dp_opt=dp_opt,dt_opt=dt_opt,dcca_opt=dcca_opt,rmse=rmse,min_RMSE=min_RMSE)

}