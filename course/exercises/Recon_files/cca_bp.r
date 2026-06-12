cca_bp<-function(P,TT,dp,dt,dcca){
###cca_bp translated from matlab function 2017/05/19 RN function [Tp,B,Tm,Pm,Ts,Ps] = cca_bp(P,T,dp,dt,dcca)

# function to make predictions using CCA-BP
#
# Input parameters:
  # P - matrix of predictor vectors (columns)
# T - matrix of predictand vectors (columns)
# Note: P and T should have the same row dimension
# dp and dt are reduction dimensions for P and T respectively, 
# per Barnett and Preisendorfer's CCA version
# dcca is the number of CCA modes to be used for prediction
# Note: dcca cannot exceed min(dp,dt)
# 
# Output parameters:
# Tp - in-sample prediction for T
# B  - regression matrix
# Tm - average of T in row dimension
# Pm - average of P in row dimension
# Ts - std of T in row dimension
# Ps - std of P in row dimension
#
# For explanations see the paper 
# Smerdon, J.E., A. Kaplan, D. Chang, and M.N. Evans, 2010: 
# A pseudoproxy evaluation of the CCA and RegEM methods for reconstructing 
# climate fields of the last millennium, J. Climate, in press.
# http://www.ldeo.columbia.edu/~jsmerdon/papers/2010a_jclim_smerdonetal.pdf
#
#  A.Kaplan, J.E.Smerdon, August 2010 
# Modified by J. Wang, USC, April 2012

  nc = dim(P)[1]
  # Time-standardized matrices
  
  Pds<- scale(P)
  Pm<-apply(P,2,mean)
  Ps<-apply(P,2,sd)
  
  Tds<- scale(TT)
  Tm<-apply(TT,2,mean)
  Ts<-apply(TT,2,sd)
  
  svdp<-svd.matlab(t(Pds))
  UP<-svdp$U
  SP<-svdp$S
  VP<-svdp$V
  
  svdt<-svd.matlab(t(Tds))
  UT<-svdt$U
  ST<-svdt$S
  VT<-svdt$V
  
  # Cross Covariance matrix
  FX <- t(VT[,1:dt])%*%VP[,1:dp]
  svdf<-svd.matlab(FX)
  UF<-svdf$U
  SF<-svdf$S
  VF<-svdf$V
  
  if(dt==1){
    utx<-t(t(UT[,1:dt]))
  }else{
    utx<-UT[,1:dt]
  }
  if(dp==1){
    B       = utx%*%ST[1:dt,1:dt]%*%UF[1:dt,1:dcca]%*%SF[1:dcca,1:dcca]%*%t(VF[1:dp,1:dcca])%*%diag(1.0/diag(SP[1:dp,1:dp],nrow=length(SP[1:dp,1:dp])))%*%t(UP[,1:dp])
  }else{
    B       = utx%*%ST[1:dt,1:dt]%*%UF[1:dt,1:dcca]%*%SF[1:dcca,1:dcca]%*%t(VF[1:dp,1:dcca])%*%diag(1.0/diag(SP[1:dp,1:dp]))%*%t(UP[,1:dp])
  }
  
  
  Tdsp = B%*%t(Pds)
  
  Tp = t(Tdsp)%*%diag(Ts) + t(t(rep(1,dim(TT)[1])))%*%Tm
  
  list(Tp=Tp,B=B,Tm=Tm,Pm=Pm,Ts=Ts,Ps=Ps)
}