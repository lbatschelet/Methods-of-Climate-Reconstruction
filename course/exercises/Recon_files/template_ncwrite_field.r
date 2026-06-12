## Function converted from the former sourced template fragment.
## Uses the objects stored in run.env and writes outputs in place.


write_field_output <- function(run.env){
  local.env <- new.env(parent = run.env)
  evalq({
## write to nc files (formerly in end template)

if(members.out==T & (write.quantiles.recon==F | bigjob==T)){
  
  if(length(dataRows)>1){
    recdata.3d<-vecttomat(recon,grid.dims,nona)
    if(onecol==T) recdata.3d<-array(recdata.3d[1,1,],dim=c(1,1,dim(recdata.3d)[3]))
  }else{
    recdata.3d<-vecttomat.1d(recon,grid.dims,nona)
    if(onecol==T) recdata.3d<-array(recdata.3d[1,1],dim=c(1,1,dim(recdata.3d)[2]))
  }
  
  
  if(bigjob==F){
    if(nblocks[[model]]==1){
      ncvar_put(out.nc,out.nc$var[[1]]$name,recdata.3d,c(1,1,ens,blocks[[model]][1]),c(length(lons),length(lats),1,length(dataRows)))
    }else{
      for (b in 1:nblocks[model]) {
        selyears<-which(models$years[[model]] %in% blocks[[model]][[b]])
        ncvar_put(out.nc,out.nc$var[[1]]$name,recdata.3d[,,selyears],c(1,1,ens,blocks[[model]][[b]][1]),c(length(lons),length(lats),1,length(blocks[[model]][[b]])))
      }
    }
  }else{
    for(y in 1:length(dataRows)){
      # out.nc<-open.nc(nc.outfiles.recon[dataRows[y]],write=T)
      if(length(dataRows)>1){ 
        ncvar_put(get(paste("out.nc",y,sep=".")),get(paste("out.nc",y,sep="."))$var[[1]]$name,recdata.3d[,,y],c(1,1,ens,1),c(length(lons),length(lats),1,1))
      }else{
        ncvar_put(get(paste("out.nc",y,sep=".")),get(paste("out.nc",y,sep="."))$var[[1]]$name,recdata.3d,c(1,1,ens,1),c(length(lons),length(lats),1,1))
      }
      # close.nc(out.nc)
    }
  }
}


if(do.field==T & do.verif==T){

  # RE
  if(do.re==T){
    if(members.out==T & write.quantiles.verif==F){
      if(length(dataRows)>1){
        resdata.3d<-vecttomat(res,grid.dims,nona)
        if(onecol==T) resdata.3d<-array(resdata.3d[1,1,],dim=c(1,1,dim(resdata.3d)[3]))
      }else{
        resdata.3d<-vecttomat.1d(res,grid.dims,nona)
        if(onecol==T) resdata.3d<-array(resdata.3d[1,1],dim=c(1,1,dim(resdata.3d)[2]))
      }
      
      # if(bigjob==F){
      if(nblocks[[model]]==1){
        ncvar_put(res.nc,res.nc$var[[1]]$name,resdata.3d,c(1,1,ens,blocks[[model]][1]),c(length(lons),length(lats),1,length(dataRows)))
        
      }else{
        for (b in 1:nblocks[model]) {
          selyears<-which(models$years[[model]] %in% blocks[[model]][[b]])
          ncvar_put(res.nc,res.nc$var[[1]]$name,resdata.3d[,,selyears],c(1,1,ens,blocks[[model]][[b]][1]),c(length(lons),length(lats),1,length(blocks[[model]][[b]])))
        }
      }

    }
  }
    
   if(do.ce==T){
    
      if(members.out==T & write.quantiles.verif==F){
        if(length(dataRows)>1){
          cesdata.3d<-vecttomat(ces,grid.dims,nona)
          if(onecol==T) cesdata.3d<-array(cesdata.3d[1,1,],dim=c(1,1,dim(cesdata.3d)[3]))
        }else{
          cesdata.3d<-vecttomat.1d(ces,grid.dims,nona)
          if(onecol==T) cesdata.3d<-array(cesdata.3d[1,1],dim=c(1,1,dim(cesdata.3d)[2]))
        }
      # if(bigjob==F){
        cesdata.3d[which(cesdata.3d==-Inf)]<-NA
        if(nblocks[[model]]==1){
          ncvar_put(ces.nc,ces.nc$var[[1]]$name,cesdata.3d,c(1,1,ens,blocks[[model]][1]),c(length(lons),length(lats),1,length(dataRows)))
          
        }else{
          for (b in 1:nblocks[model]) {
            selyears<-which(models$years[[model]] %in% blocks[[model]][[b]])
            ncvar_put(ces.nc,ces.nc$var[[1]]$name,cesdata.3d[,,selyears],c(1,1,ens,blocks[[model]][[b]][1]),c(length(lons),length(lats),1,length(blocks[[model]][[b]])))
          }
        }

      }
   }
  
  if(do.r2==T){
    if(members.out==T & write.quantiles.verif==F){
      if(length(dataRows)>1){
        r2data.3d<-vecttomat(r2,grid.dims,nona)
        if(onecol==T) r2data.3d<-array(r2data.3d[1,1,],dim=c(1,1,dim(r2data.3d)[3]))
      }else{
        r2data.3d<-vecttomat.1d(r2,grid.dims,nona)
        if(onecol==T) r2data.3d<-array(r2data.3d[1,1],dim=c(1,1,dim(r2data.3d)[2]))
      }
      
      # if(bigjob==F){
      if(nblocks[[model]]==1){
        ncvar_put(r2.nc,r2.nc$var[[1]]$name,r2data.3d,c(1,1,ens,blocks[[model]][1]),c(length(lons),length(lats),1,length(dataRows)))
      }else{
        for (b in 1:nblocks[model]) {
          selyears<-which(models$years[[model]] %in% blocks[[model]][[b]])
          ncvar_put(r2.nc,r2.nc$var[[1]]$name,r2data.3d[,,selyears],c(1,1,ens,blocks[[model]][[b]][1]),c(length(lons),length(lats),1,length(blocks[[model]][[b]])))
        }
      }
    }
    
  }
  
}
    
    
    
    
    

  }, envir = local.env)
  invisible(NULL)
}

