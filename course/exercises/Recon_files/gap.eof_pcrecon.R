gap.eof.neu.pcrecon <- function(dat, method="ols",eofmax=min(dim(dat))-1,err = c(1,1e-3), plot = TRUE, min.num = 40, max.iter = 500, keep.obs = TRUE) {
##############################################################################
## DESCRIPTION: computes EOF's, also for missing data matrixes, it "optimally
##              interpolates" missing values using the
##              EOF (SVD) method by Beckers et al. (2003)
##              this method needs no a priori information and is parameter free
## -------------------------------------------------------------------------
## ARGUMENTS
## Required: dat         (Data matrix, NA's must be included :-))
## Optional: 
            # method    method used for the regression.must be "Ols" or "lms"
##              err         (convergence criterion within a EOF n search,  
##                        percentage change from one to next "loop",abs. value [units])
##           cross.fract (fraction of data used to do cross-validation) 0.05=5%
##           min.num     (minimum number of points used for cross validation)
##           plot        (if TRUE, then graphical result output is given)
##           max.iter    (max. number of iterations within EOF reconstruction)
##           keep.obs    (keep the measured values, do not suppress noise in measurements)  
## -------------------------------------------------------------------------
## VALUE: interp.field, noEOF, cross.MSE, expl.var, singv, singv.zero,
##        EOF, PC, PCnorm, U, D
## -------------------------------------------------------------------------
## other functions needed: eof
## -------------------------------------------------------------------------
## DETAILS: 1) the routine makes the assumption, that the data is noise free
##             hence, the routine is thought for a priori homogenized data
##          2) if there are no missing values, the routine gives you the optimal
##             EOF truncation and the EOF's of the data  
##  
## CAUTION: the data matrix MUST be given with timesteps in ROWS and stations/xy
##          grid in COLUMNS
##                    .. .. -> xy
##                     | .. .. ..
##                     v .. .. ..
##                     t .. .. ..
##
##          Note: the data in the input data matrix should be normalised and
##                approximately Gaussian, this is a problem especially for
##                daily data as e.g. precip, snow etc..   
##          
##          if the number of timesteps is larger than stations/xypoints, the
##          routine automatically transposes the input matrix, the computation
##          is done without difference, also the output is ok  
## -------------------------------------------------------------------------
## Author: Simon Scherrer (NCCR@MCH), Date: 28 April 03
## Last modified: 8 December 2005, 6 August 2008 --> keep.obs option included
##############################################################################
  
# check for largest dimension
# if timesteps > xy dimension then swap

orig.dim1 <- dim(dat)[1]
orig.dim2 <- dim(dat)[2]

#identfiy complete records used for infilling
nnas<-apply(dat,2,function(x) length(which(is.na(x))))
nonas<-which(nnas==0)
missvals<-apply(dat,2,function(x) which(is.na(x)))

if (orig.dim1 <= orig.dim2) dat <- t(dat)
  
nrow <- dim(dat)[1]
ncol <- dim(dat)[2]

# determine eff. num. of cross-validation points
n <- round(min(0.01*nrow*ncol+min.num,0.03*nrow*ncol))
                                        
orig <- dat
num.na <- sum(is.na(orig))

# determine missing.coordinate
miss.coord <- numeric(0)

for (i in 1:nrow){ # loop over rows
  for (j in 1:ncol){ # loop over columns
    if (is.na(orig[i,j])) {
      # store missing coordinate pairs
      miss.coord <- rbind(miss.coord,cbind(i,j))
    }
  }
}

#for pc rec don't set NAs to 0
orig0<-dat
orig0[is.na(orig0)] <- 0.

cat("\nYour data set has ",num.na, "NA's.\nThat is a rate of ",
    num.na/(nrow*ncol)*100,"%\n\n")
cat("NA points are:\n")
for (i in 1:length(miss.coord[,1])){
  cat("[",miss.coord[i,1],";",miss.coord[i,2],"]")
}

dum <- array(NA,c(nrow,ncol))
diff <- array(NA,c(max.iter,max.iter))
# x can't be a 4 dim vector, because of memory problems (needs sev. GByte RAM)
# EOF "n" array x and result array res
x <- xold <- array(NA,c(nrow,ncol))
res <-  array(NA,c(eofmax+1,nrow,ncol))
x <- orig


######################################################################
# SVD EOF m reconstruction for missing values, result is saved in res array

cat("\n\nBegin EOF reconstruction for MISSING VALUES\n\n")

for (m in 1:eofmax){
   EOFnum <- m
   i <- 1
   diff[,1] <- 10*err  # set a first diff value
   diff.old <- 1e6
   while (((diff[m,i]-diff.old)/diff.old) < -err[1]/100 &&  diff[m,i] > err[2] && i <= max.iter) {
#        c <- La.svd(x)
#        dum <- cbind(c$u[,1:EOFnum],c$u[,(EOFnum+1):dim(orig)[2]]*0.) %*%
#               diag(c$d) %*% cbind(c$vt) #cbind(c$vt[,1:EOFnum],c$vt[,(EOFnum+1):dim(orig)[2]]*0.)
#      
     npc<-EOFnum
     do.t<-ifelse (orig.dim1 <= orig.dim2,T,F)
     if(method=="ols"){
       dum<-regrecon(x,do.t,npc,nonas,missvals)
     }else{
       dum<-lmsregrecon(x,do.t,npc,nonas,missvals)
     }
       
       
         xold <- x
         xold[is.na(xold)] <- 0.
       # change only values of missing and cross validation points
       x[miss.coord] <- dum[miss.coord]       
              
     # convergence criterion determination, a simple MAE
     i <- i+1
     diff[m,i] <- abs(mean(x[miss.coord] - xold[miss.coord]))
     diff.old <- diff[m,i-1]  
   }
   res[m,,] <- x
   
   if (i < 500) {
   cat("EOF",m,"reconstruction converged with abs(mean(diff))=",diff[m,i]," after ",i-1,"loop(s)\n")
   }
   else {
     cat("!!!!!   CAUTION:   !!!!!\n")
     cat("EOF",m," did NOT CONVERGE after ",i-1,"loops.\n")
     cat("Mean difference was:",diff[m,i],". Error limit is",err,
         " currently.\n")
     cat("Consider a lower error (err) limit!\n\n")
   }
}

cv <- TRUE
if (cv==TRUE){
  
  # determine a subset of "orig" cross-validation points
  # to determine, how much EOF's should be retained
  # put them to zero too, as the NA values have been above

  cross.coord <- numeric(0)
  k <- 0

  repeat {
    nrow.cand <- round(runif(1, min=1, max=nrow))
    ncol.cand <- round(runif(1, min=1, max=ncol))
    if (!is.na(dat[nrow.cand,ncol.cand]) && dat[nrow.cand,ncol.cand]!=0) { # dat in orig vector with NA's
       res[1,nrow.cand,ncol.cand] = 0.
       # store successful coordinate pairs
       cross.coord <- rbind(cross.coord,cbind(nrow.cand,ncol.cand))
       # increase counter
       k <- k+1
      if (k == n+1) break 
    }
  }
}

cat("\n\n",n,"cross-validation points determined:\n")

for (i in 1:length(cross.coord[,1])){
cat("[",cross.coord[i,1],";",cross.coord[i,2],"]")
}

dum <- res2d <- best.recon <- array(NA,c(nrow,ncol))
diff <- array(NA,c(max.iter,max.iter))
mse <- array(NA,ncol)

######################################################################
# repeat reconstruction process in a cross-validation manner with
# the res array as start values, determine opt. number of EOFs and error 

cat("\n\nBegin EOF reconstruction for CROSS-VALIDATION and error determination\n")

for (m in 1:eofmax){
   EOFnum <- m
   i <- 1
   diff[,1] <- 10*err  # set a first diff value
   diff.old <- 1e6
   res2d <- res[m,,] 
   while (((diff[m,i]-diff.old)/diff.old) < -err[1]/100 &&  diff[m,i] > err[2] && i <= max.iter) {
#        c <- La.svd(res2d)
#        dum <- cbind(c$u[,1:EOFnum],c$u[,(EOFnum+1):dim(orig)[2]]*0.) %*%
#               diag(c$d) %*% cbind(c$vt) #cbind(c$vt[,1:EOFnum],c$vt[,(EOFnum+1):dim(orig)[2]]*0.)
       npc<-EOFnum
       do.t<-ifelse (orig.dim1 <= orig.dim2,T,F)
       if(method=="ols"){
         dum<-regrecon.all(res2d,do.t,npc,nonas)
       }else{
         dum<-lmsregrecon.all(res2d,do.t,npc,nonas,missvals)
       }
       
       
       
       
              
       old.guess <- res2d
       # change only cross-validation points
       res2d[cross.coord] <- dum[cross.coord]       
              
     # convergence criterion determination, a simple MAE
     i <- i+1
     diff[m,i] <- abs(mean(res2d[cross.coord] - old.guess[cross.coord]))
     diff.old <- diff[m,i-1]  
   }
   
   if (i < 500) {
   cat("EOF",m,"reconstruction converged with abs(mean(diff))=",diff[m,i]," after ",i-1,"loop(s)\n")
   }
   else {
     cat("!!!!!   CAUTION:   !!!!!\n")
     cat("EOF",m," did NOT CONVERGE after ",i-1,"loops.\n")
     cat("Mean difference was:",diff[m,i],". Error limit is",err,
         " currently.\n")
     cat("Consider a lower error (err) limit!\n\n")
   }

   # determine cross validation RMS error for EOF m
   mse[EOFnum] <- mean((res2d[cross.coord] - orig[cross.coord])^2)
   cat("MSE of EOF",m,"reconstruction is: ",mse[m],"\n")
   
   # define m+1 starting values
   for (i in 1:length(cross.coord[,1])){
     res[m+1,cross.coord[i,1],cross.coord[i,2]] <- dum[cross.coord[i,1],cross.coord[i,2]]
   }
}

   # determine the best reconstruction
   min <- min(mse,na.rm=T)
   best.num <- 1
   repeat {
     if (mse[best.num] != min) best.num <- best.num + 1
     else break
   }

   best.recon <- res[best.num,,]
   # re-input the coordinates that were used for cross-validation
   best.recon[cross.coord] <- orig[cross.coord]
 
   cat("\nThe best \"filler\" retains ",best.num," EOF's!\n\n")

   if (keep.obs==FALSE) {
   # compute best reconstruction, suppressing noise 
#    tmp <- La.svd(x)
#    best.recon <- cbind(tmp$u[,1:best.num],tmp$u[,(best.num+1):dim(orig)[2]]*0.) %*%
#                        diag(tmp$d) %*% cbind(tmp$vt)
     stop("currently not implemented: keep.obs==F")
   }

   best.svd <- eof(best.recon,plot=F)

   cat("The best filled data set first ",best.num," EOF explain (%):\n")
   cat(best.svd$expl.var[1:best.num],"\n\n")
   cat("The zero filled data set first ",best.num," EOF explain (%):\n")
   eof.zero <- eof(orig0,plot=F)
   cat(eof.zero$expl.var[1:best.num],"\n\n")

   # output values
   r <- list(interp.field = best.recon, noEOF = best.num, convergenceEOF = diff, cross.MSE = mse,
             expl.var = best.svd$expl.var, singv = best.svd$singv,
             singv.zero = eof.zero$singv, EOF = best.svd$EOF, PC = best.svd$PC,
             PCnorm = best.svd$PCnorm, U = best.svd$U, D = best.svd$D)

   cat("Consider the following output variables:\n")
   cat(names(r),"\n")

   if (plot == T) {

   # some plots
   tick.lab <- c(7,7,7)
   par(mfrow=c(2,2),las=0,mar=c(3,3.5,2,2),mgp=c(2,0.7,0),pty="s",cex=0.8,lab=tick.lab,tck=0.02)
   # plot MSE as a function of EOF
   plot(mse,xlab="EOF number",ylab="cross-validation: mean squared error (MSE)",
   main="cross-validated MSE error [units^2]",panel.first=grid(),type="b")
   points(mse,pch=19)
   abline(v=best.num,col=2,lty=2)
   #text(best.num+1,min(mse)+range(mse)/3,best.num)
   
   # plot singular values for EOD interpolation and zero filled fields
   plot(best.svd$singv[1:best.num],type="b", xlab="EOF number", ylab="singular value",
        main="sing. values (line = EOF interp., dots = 0 filled)",
        panel.first=grid(),pch=19)#,log="y")
   points(eof.zero$singv[1:best.num],type="b")

   #plot sing. value differences
   plot(best.svd$singv[1:best.num]-eof.zero$singv[1:best.num],type="b", xlab="EOF number", ylab="sing. value difference",
        main="sing. value diff. (EOF interp. - 0 filled)",
        panel.first=grid(),pch=19)
   abline(h=0,lty=2)
 
   
   # plot of expl. variance
   # North test plot, first n EOFs
   num <- best.num
   dof <- min(orig.dim1,orig.dim2)
   plot(best.svd$singv[1:num]^2/sum(best.svd$singv^2)*100.,
   main = "Explained variance and error (North crit.)",
   type = "n", ylab = "explained variance (%)",
   xlab = "EOF order",xlim=c(1,num),
   panel.first=abline(v=seq(1,best.num,1), h=c(seq(0.1,1,0.1),seq(2,9,1),seq(10,100,10)),lty=3,col="grey85"),log="y",
   ylim=c(1e-1,(best.svd$singv[1]^2+sqrt(2/dof)*
         (best.svd$singv[1]^2))/sum(best.svd$singv^2)*100.+5),axes=F)
   axis(1,1:best.num)
   axis(2,0.2)
   axis(2,c(0.1,0.5))
   axis(2,c(1,2,5,10,20,50,100))
   axis(2,100)

   lines(best.svd$singv[1:num]^2/sum(best.svd$singv^2)*100.,
         lty = 3)
     
   for (i in 1:num) {
        lines(rep(i, 2), c((best.svd$singv[i]^2+sqrt(2/dof)*
                         best.svd$singv[i]^2)/sum(best.svd$singv^2)*100.,
        (best.svd$singv[i]^2 - sqrt(2/dof)*best.svd$singv[i]^2)/
                           sum(best.svd$singv^2)*100.),
        lty = 2, col = "darkgrey")
        lines(c(i - 0.25, i + 0.25), rep((best.svd$singv[i]^2+
              sqrt(2/dof)*best.svd$singv[i]^2)/sum(best.svd$singv^2)*100.,
              2), lwd = 2, col = "darkgrey")
        lines(c(i - 0.25, i + 0.25), rep((best.svd$singv[i]^2-
              sqrt(2/dof)*best.svd$singv[i]^2)/sum(best.svd$singv^2)*100.,
              2), lwd = 2, col = "darkgrey")
    }

   points(best.svd$singv[1:num]^2/sum(best.svd$singv^2)*100.)
   points(best.svd$singv[1:num]^2/sum(best.svd$singv^2)*100.,
          pch = 20, cex = 0.8, col = "red")
   box()

   } # end if (plot == T)
   
   class(r) <- "geof"
   r

} # end function 
##############################################################################







eof <- function(x, center = FALSE, scale. = FALSE, plot = T) {
##############################################################################
## DESCRIPTION: computes "normal" EOF's (no missing values)
## -------------------------------------------------------------------------
## ARGUMENTS
## Required: x (Data matrix)
## Optional: center (demean), scale. = (normalize), tol = tolerance display only
## -------------------------------------------------------------------------
## VALUE:
## -------------------------------------------------------------------------
## other functions needed: 
## -------------------------------------------------------------------------
## DETAILS:
##
## CAUTION: the data matrix MUST be given with timesteps in ROWS and stations/xy
##          grid in COLUMNS
##                    .. .. -> xy
##                     | .. .. ..
##                     v .. .. ..
##                     t .. .. ..
##
##          the function computes the SVD for the optimal matrix(M,N) (M > N)   
## -------------------------------------------------------------------------
## Author: Simon Scherrer (NCCR@MCH), Date: 23 April 03
## Last modified: 30 April 2003
##############################################################################
  
   # check for bigger dimension and make the row dimension > the colum dimension
   dim.orig1 <- dim(x)[1]
   dim.orig2 <- dim(x)[2]
   
   x <- as.matrix(x)
   x <- scale(x, center = center, scale = scale.) # normalise or center the data
    
   if (dim.orig1 < dim.orig2) x <- t(x)
   
   s <- La.svd(x)		                	   # SVD decomposition
 
   u <- s$u
   d <- diag(s$d)
   vT<- s$vt


   if (dim.orig1 >= dim.orig2) {
     
    s$d <- s$d/sqrt(max(1, nrow(x) - 1))           # standard deviation
    ev <- s$d^2                                    # singular values -> eigenvalues
    expl.var <- ev/sum(ev)*100.			   # expl. variance vector (%)
    PC <- x %*% t(vT)				   # PC's renorm on data
    # u %*% diag(s$d)
    PCsvd <- scale(u, scale=T)    		   # PC's directly from SVD decomp.
						   # but standardized
    r <- list(expl.var = expl.var, singv = s$d, eigenv = ev, PC=PC, EOF = vT, PCnorm = PCsvd,
              sdev = s$d, U = u, D = d, vT = vT)
   } 

   if (dim.orig1 < dim.orig2) {
     
    s$d <- s$d/sqrt(max(1, nrow(x) - 1))           # standard deviation
    ev <- s$d^2                                    # singular values -> eigenvalues
    expl.var <- ev/sum(ev)*100.			   # expl. variance vector (%)
    PC <- t(x)  %*% s$u          	      	   # PC's renorm on data
    PCsvd <- scale(t(vT), scale=T)	           # PC's directly from SVD decomp.
						   # but standardized
    r <- list(expl.var = expl.var, singv = s$d, eigenv = ev, EOF = s$u, PC = PC, PCnorm = PCsvd,
              sdev = s$d, U = u, D = d, vT = vT)
   } 


   if (plot == T) {
   tick.lab <- c(10,10,10)
   par(mfrow=c(1,1),lab = tick.lab)

   # plot of expl. variance
   # North test plot, first n EOFs
   num <- 12
   dof <- min(dim.orig1,dim.orig2)
   plot(ev[1:num]/sum(ev)*100.,
   main = "Explained variance and error (North (1982) criterion)",
   type = "n", ylab = "explained variance (%)",
   xlab = "EOF order",xlim=c(1,num),
   ylim=c(0,(ev[1]+sqrt(2/dof)*ev[1])/sum(ev)*100.+5))

   lines(ev[1:num]/sum(ev)*100.,lty = 3,panel.first=grid())

   for (i in 1:num) {
        lines(rep(i, 2), c((ev[i]+sqrt(2/dof)*ev[i])/sum(ev)*100.,
        (ev[i] - sqrt(2/dof)*ev[i])/sum(ev)*100.),
        lty = 2, col = "darkgrey")
        lines(c(i - 0.25, i + 0.25), rep((ev[i]+sqrt(2/dof)*
              ev[i])/sum(ev)*100., 2), lwd = 2, col = "darkgrey")
        lines(c(i - 0.25, i + 0.25), rep((ev[i]-sqrt(2/dof)*
              ev[i])/sum(ev)*100., 2), lwd = 2, col = "darkgrey")
    }

   points(ev[1:num]/sum(ev)*100.)
   points(ev[1:num]/sum(ev)*100., pch = 20, cex = 0.8, col = "red")
   box()

   } # end if (plot == T)

   cat("Consider the following output variables:\n")
   cat(names(r),"\n")
   class(r) <- "eof"
   r
}
## end function eof






##regression recon--------------------
regrecon<-function(data,do.t=F,npc=1,nonas=NULL,missvals){
  if(do.t==T) data<-t(data)
  data.fill<-data
  if(is.null(nonas)){
    nnas<-apply(data,2,function(x) length(which(is.na(x))))
    nonas<-which(nnas==0)
  }
  
  if(length(nonas)<npc) npc<-length(nonas)
  
  for(pr in 1:dim(data)[2]){
    if((pr %in% nonas)==F){     ##(length(miss)>0){
      miss<-missvals[[pr]]  #which(is.na(data[,pr]))
      cn<-La.svd(data[,nonas])
      y<-t(t(data[-miss,pr]))
      
      xs<-cn$u  %*% diag(cn$d)
      xs<-cbind(1,xs[,1:npc])
      xsc<-xs[-miss,]
      xsr<-xs[miss,]
      
      recb<-solve(crossprod(xsc)) %*% t(xsc) %*% y
    #  rec<-xsc %*% recb
      recna<-xsr %*% recb
      data.fill[miss,pr]<-recna
      # pl.mts(cbind(data.fill[,pr],data[,pr]))
#       test<-data.fill[,pr]
#       test[-miss]<-rec
#       pl.mts(cbind(data.fill[,pr],data[,pr],test))
#       
    }
  }
  if(do.t==T) data.fill<-t(data.fill)
  data.fill
}


##regression recon no missing values--------------------
regrecon.all<-function(data,do.t=F,npc=1,nonas=NULL){
  if(do.t==T) data<-t(data)
  data.fill<-data
  
  if(is.null(nonas)){
    nnas<-apply(data,2,function(x) length(which(is.na(x))))
    nonas<-which(nnas==0)
  }
  
  if(length(nonas)<npc) npc<-length(nonas)
  
  for(pr in 1:dim(data)[2]){
    if((pr %in% nonas)==F){ 
      cn<-La.svd(data[,nonas])
      y<-t(t(data[,pr]))
      
      xs<-cn$u  %*% diag(cn$d)
      xs<-cbind(1,xs[,1:npc])
      recb<-solve(crossprod(xs)) %*% t(xs) %*% y
      recna<-xs %*% recb
      data.fill[,pr]<-recna
      # pl.mts(cbind(data.fill[,pr],data[,pr]))
      #       test<-data.fill[,pr]
      #       test[-miss]<-rec
      #       pl.mts(cbind(data.fill[,pr],data[,pr],test))
      #       
    }
  }
  if(do.t==T) data.fill<-t(data.fill)
  data.fill
}



##lms regression recon--------------------
lmsregrecon<-function(data,do.t=F,npc=1,nonas,missvals){
  if(do.t==T) data<-t(data)
  data.fill<-data
  if(is.null(nonas)){
    nnas<-apply(data,2,function(x) length(which(is.na(x))))
    nonas<-which(nnas==0)
  }
  
  if(length(nonas)<npc) npc<-length(nonas)
  
  for(pr in 1:dim(data)[2]){
    if((pr %in% nonas)==F){     ##(length(miss)>0){
      miss<-missvals[[pr]]  #which(is.na(data[,pr]))
      cn<-La.svd(data[,nonas])
      y<-t(t(data[-miss,pr]))
      
      xs<-cn$u  %*% diag(cn$d)
      xs<- t(t(xs[,1:npc]))# xs<-cbind(1,xs[,npc])
      xsc<-xs[-miss,]
      xsr<-xs[miss,]
      if(length(miss)==1) xsr<-t(xsr)
      
      recb<-lmsreg(x=xsc,y=y,method=lms)
    #  rec<-cbind(1,xsc) %*% recb$coefficients
      recna<-cbind(1,xsr) %*% recb$coefficients
      data.fill[miss,pr]<-recna
      # pl.mts(cbind(data.fill[,pr],data[,pr]))
      #       test<-data.fill[,pr]
      #       test[-miss]<-rec
      #       pl.mts(cbind(data.fill[,pr],data[,pr],test))
      #       
    }
  }
  if(do.t==T) data.fill<-t(data.fill)
  data.fill
}


##lms regression recon--------------------
lmsregrecon.all<-function(data,do.t=F,npc=1,nonas=NULL,missvals){
  if(do.t==T) data<-t(data)
  data.fill<-data
  if(is.null(nonas)){
    nnas<-apply(data,2,function(x) length(which(is.na(x))))
    nonas<-which(nnas==0)
  }
  
  if(length(nonas)<npc) npc<-length(nonas)
  
  for(pr in 1:dim(data)[2]){
    if((pr %in% nonas)==F){     ##(length(miss)>0){
      miss<-missvals[[pr]]  #which(is.na(data[,pr]))
      cn<-La.svd(data[,nonas])
      y<-t(t(data[,pr]))
      
      xs<-cn$u  %*% diag(cn$d)
      xs<- t(t(xs[,1:npc]))# xs<-cbind(1,xs[,npc])

      recb<-lmsreg(x=xs,y=y,method=lms)
      recna<-cbind(1,xs) %*% recb$coefficients
      data.fill[,pr]<-recna
      # pl.mts(cbind(data.fill[,pr],data[,pr]))
      #       test<-data.fill[,pr]
      #       test[-miss]<-rec
      #       pl.mts(cbind(data.fill[,pr],data[,pr],test))
      #       
    }
  }
  if(do.t==T) data.fill<-t(data.fill)
  data.fill
}