#####File from SOM of Werner & Tingley 2015(CP)
##cov width: temporal mean of of the ensemble 5%-95% spread at each location
##covRate: Fraction of years where the target is within hte 90% ensemble range at each location
##--> if below 0.9, the error range is too narrow, otherwise too wide.
##REli: similar to difference between CI and empirical CI i.e. abs(0.9-CovRate)
##potCRPS: Similar to mean ensemble spread
###CRPS: Reli + potCRPS


MakeConfIntervals <- function( x, probs=c(.05, .95)){
  # Construct confidence intervals, can be symmetric or asymmetric. If passing a
  # single number to probs it is symmetric. When passing two numbers it can be
  # forced to be asymmetric, though I do not know why one would like that...
  if ( length(probs) == 1) {probs <- (1 + c(-1,1)*probs)/2}
  if( dim( x)[2] > 1){
    CI <- matrix(NA, dim(x)[1],2)
    CI[,1] <- apply( x, 1, quantile, probs=probs[1], na.rm=TRUE)
    CI[,2] <- apply( x, 1, quantile, probs=probs[2], na.rm=TRUE)
    return( CI)
  } else {
    warning( "Use with a *matrix* of input data!")
    return( x)
  }
}

CoverageRate <- function( x, Target, probs = c(.05, .95)){
  CI <- MakeConfIntervals( x, probs = probs)
  return( mean(rowSums(Target < CI) * rowSums(Target > CI) ) )
}

CRPS <- function( XEns, x.a){
  # Continous Ranked Probability Score as in Gneiting and Raftery "Proper
  # Scoring Rules", JAmStat (2007), using the form of eq. (20) and hoping for
  # Xens to be 1. independent and 2. a large enough ensemble that the expected
  # value of the evaluated differences is close enough to the mean.
  # See H. Hersbach, JAmetSoc (2002) Sec 4b for the implementation (also used at
  # ECMWF?), especially look at (24--27) for the notation used here.
  x.i <- sort( XEns)
  p.i <- seq(0,length(x.i))/length(x.i)
  idx.eval <- ( x.i < x.a)
  alpha.idx <- which(idx.eval)
  beta.idx <- which(!idx.eval)
  alpha <- c( x.i[alpha.idx], x.a)
  beta <- c( x.a, x.i[beta.idx])
  return( sum( diff(alpha)*p.i[alpha.idx+1]**2) + sum( (1-p.i[beta.idx])**2*diff(beta)) )
}

avgCRPS <- function( x, y){
  if( length( y) != dim( x)[1]){
    warning("incompatible dimensions for Target 'y' and Forecast 'x'!")
    return(FALSE)
  } else {
    XEns <- t( apply( x, 1, sort)) #sorts the ensemble recon at each time step
    p.i <- (0:dim(x)[2])/dim(x)[2]
    idx.eval <- XEns < y # which recons are smaller than the target
    alpha.idx <- which( idx.eval, arr.ind=TRUE) # which realizations are smaller than target
    beta.idx <- which( !idx.eval, arr.ind=TRUE) #which realizations are not smaller than target (larger)
    alpha <- matrix( NA,dim(XEns)[1],dim(XEns)[2]+1) #
    alpha[alpha.idx] <- XEns[alpha.idx]
    alpha[cbind(1:dim(XEns)[1],rowSums( idx.eval) + 1)] <- y #alpha: in each year c(all recons that are small than target, the target, NAs)
    a.i.mat <- t( apply(alpha,1,diff)) #difference between these values for each year
    a.i.mat[ is.na(a.i.mat)] <- 0
    beta <- matrix( NA,dim(XEns)[1],dim(XEns)[2]+1)
    beta[beta.idx+matrix(c(0, 1), dim(beta.idx)[1], 2, byrow = TRUE)] <- XEns[beta.idx]
    beta[cbind(1:dim(XEns)[1],rowSums( idx.eval) + 1)] <- y #beta: in each year c(NA,the target, all recons that are larger than the target)
    b.i.mat <- t( apply(beta,1,diff)) #difference between these values
    b.i.mat[ is.na(b.i.mat)] <- 0
    # the index i goes from 1..N for a.i and 0..N-1 for b.i !
    g.i <- c( 0, colMeans(a.i.mat)) + c( colMeans( b.i.mat), 0) #colMeans(a.i.mat)= average difference between lowest and second lowest ense member etc.. (for all that are below y). same for b.i. but for those that are avove
    o.i <- c( colMeans(b.i.mat), 0) / g.i
    Reli <- sum( g.i*( o.i - p.i)**2, na.rm=TRUE)
    CRPSpot <- sum( g.i * o.i * (1 - o.i), na.rm=TRUE)
    return( list( CRPSmean = Reli+CRPSpot, CRPSpot=CRPSpot, Reli=Reli))
  }
}


SGenEntropy <- function( x, y){
  # This only works if the number of ensemble members is bigger than the number
  # of times (or locations) we are evaluating over, dim(x)[1] < dim(x)[2]
  if( length( y) != dim( x)[1]){
    warning("incompatible dimensions for Target 'y' and Forecast 'x'!")
    return(FALSE)
  } else {
    Sigma.P <- cov( t(x) )
    mu.P <- rowMeans( x)
    return(-log(prod(diag(chol( Sigma.P)))) - (y - mu.P) %*% solve(Sigma.P) %*% (y - mu.P))
  }
}

IntervalScore <- function( x, y, alpha){
  if( !is.null(dim(x)) ){
    CI <- MakeConfIntervals( x, probs=alpha)
    return(CI[,2] - CI[,1] + 2/alpha*( (CI[,1] - y)*(CI[,1] > y) + (CI[,2] - y)*(CI[,2] < y) ))
  } else {
    CI <- quantile( x, probs=(1 + c(-1,1)*alpha)/2)
    return(CI[2] - CI[1] + 2/alpha*( (CI[1] - y)*(CI[1] > y) + (CI[2] - y)*(CI[2] < y) ))
  }
}

