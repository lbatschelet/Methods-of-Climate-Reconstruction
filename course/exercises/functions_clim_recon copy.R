# read in CRUTEM NetCDF gridded temperature 
# set syr and eyr to cut period 
read_cru <- function(syr,eyr,crupath,crufilename){
  nc=nc_open(paste0(crupath,crufilename), write=F)
  print(nc)
  t2m_tmp=ncvar_get(nc, "temperature_anomaly") # for CRU temp
  t <- ncvar_get(nc, "time")
  tunits <- ncatt_get(nc, "time", "units")
  tustr <- strsplit(tunits$value, " ")
  origin <- paste(unlist(tustr)[3],unlist(tustr)[4])
  datesd <- as.POSIXct(t*24*60*60, origin=origin) # adress time values * times seconds and insert "origin" 
  dates <- format(datesd, "%Y-%m-%d %H:%M")  # gives dates in the format y-m-d-hh 
  dates
  yr <- unlist(strsplit(dates, "-"))[seq(1,length(dates)*3,3)]
  mon <- unlist(strsplit(dates, "-"))[seq(2,length(dates)*3,3)]
  timepos <- which((yr>(syr-1)) & (yr<(eyr+1)))
  t2m <- t2m_tmp[,,timepos] # from year 1901-2015
  lonlist=ncvar_get(nc, "longitude") 
  #lonlist[lonlist > 180] <- lonlist[lonlist > 180] - 360 # in case of lon 0 to 360
  latlist=ncvar_get(nc, "latitude") 
  nc_close(nc)
  cru <- list(data=t2m,lon=lonlist,lat=latlist,
                 time=dates[timepos])
  return(cru)
}



# read in GHCN version 4 (temperature)
# set syr and eyr to cut period 
# set refsyr and refeyr to require that series have any data in 
# this period. can be used to extract only long time series
read_ghcn_refyr <- function(syr,eyr,refsyr,refeyr,ghcntemppath,ghcnfilename,ghcnmetafilename){
  years <- syr:eyr
  ryrs <- refsyr:refeyr
  stationfile <- paste(ghcntemppath,ghcnfilename,sep='')
  ww <- c(11, -1 , 8, -1, 9, -1, 6, -1, 30, -1, 4, 1, -1, 4, 2, 2, 2, 2, 1, 2, 16, 1)
  stations <- read.fwf(stationfile, ww, header=F, fill=T, na.string='-999')
  colnames(stations) <- c('ID', 'LATITUDE', 'LONGITUDE', 'STNELEV', 'NAME', 'GRELEV', 'POPCLS', 'POPSIZ', 'TOPO', 'STVEG', 'STLOC', 'OCNDIS', 'AIRSTN', 'TOWNDIS', 'GRVEG', 'POPCSS')
  stations[stations == -9999 ] <- NA
  ww2 <- c(11, 4, 4, rep(c(5,1,1,1),12)) 
  tmp <- as.matrix(read.fwf(paste(ghcntemppath,ghcnmetafilename,sep=''), ww2, header=F, fill=T, na.string='-9999'))
  colnames(tmp) <- c('ID', 'YEAR', 'ELEMENT', 'VALUE1', 'DMFLAG1', 'QCFLAG1', 'DSFLAG1', 'VALUE2', 'DMFLAG2', 'QCFLAG2', 'DSFLAG2', 'VALUE3', 'DMFLAG3', 'QCFLAG3', 'DSFLAG3', 'VALUE4', 'DMFLAG4', 'QCFLAG4', 'DSFLAG4', 'VALUE5', 'DMFLAG5', 'QCFLAG5', 'DSFLAG5', 'VALUE6', 'DMFLAG6', 'QCFLAG6', 'DSFLAG6', 'VALUE7', 'DMFLAG7', 'QCFLAG7', 'DSFLAG7', 'VALUE8', 'DMFLAG8', 'QCFLAG8', 'DSFLAG8', 'VALUE9', 'DMFLAG9', 'QCFLAG9', 'DSFLAG9', 'VALUE10', 'DMFLAG10', 'QCFLAG10', 'DSFLAG10', 'VALUE11', 'DMFLAG11', 'QCFLAG11', 'DSFLAG11', 'VALUE12', 'DMFLAG12', 'QCFLAG12', 'DSFLAG12')
  # ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v4/README
  tmp2 <- sort(unique(tmp[,1]))
  tmp3 <- sort(unique(stations[,1]))
  tmp4 <- match(tmp2,tmp3)
  tmp5 <- which(tmp4>0)
  stat.i <- tmp2[tmp5]
  tmp.arr <- array(NA,c(length(years), 12, length(stat.i)))
  for (i in 1:length(seq(stat.i))){
    if (i %% 100 == 0) { print(i) }
    s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
    y.i <- which(years %in% tmp[s.i,2])
    r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,2] == years[x])))
    ry.i <- which(ryrs %in% tmp[s.i,2])
    if (length(ry.i) > 0) {
      #print(paste(i,'of',length(seq(stat.i))))
      tmp.arr[y.i,,i] <- tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]
    } else {
      if (length(y.i)==1) {
        tmp.arr[y.i,,i] <- array(NA,length(tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]))
      } else {
        tmp.arr[y.i,,i] <- array(NA,dim(tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]))
      }
    }  
  }
  stations.new <- stations[which(stations[,'ID'] %in% stat.i),]
  # check if at least 80% of data is available and NOT missing (NA)
  #mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])
  ## check if at least 1% of data is available and NOT missing (NA)
  mask <- apply(!is.na(tmp.arr), 3, sum) > 0.01 * prod(dim(tmp.arr)[1:2])
  ghcn.data.tmp <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
  ghcn.data=array(as.numeric(ghcn.data.tmp),c(dim(ghcn.data.tmp)[1],dim(ghcn.data.tmp)[2]))
  ghcn <- list(data=ghcn.data/100, 
               lon=stations.new[mask,'LONGITUDE'], 
               lat=stations.new[mask,'LATITUDE'], 
               names=gsub(' *$', '', as.character(stations.new[mask,'NAME'])), 
               height=stations.new[mask,'STNELEV'], 
               time=seq(min(years) + 1/24, by=1/12, length=nrow(ghcn.data)))
  return(ghcn)
}



# find grid boxes matching station coordinates
# station and grid must be lists with $lon and $lat
getgridboxnum <- function(station, grid) {
  m <- k <- l <- rep(NA,length(station$lon))
  gridatprox.arr <- array(NA,c(length(station$lon), length(station$time)))
  for(i in 1:(length(station$lon))){
    #if (i %% 100 == 0) { 
    #  print(paste("station number:",i))
    #}
    slon <- station$lon[i]
    slat <- station$lat[i]
    glon <- grid$lon[!is.na(grid$lon)]
    glat <- grid$lat[!is.na(grid$lat)]
    k[i]=which(round(abs(glon-slon),digits=3)==min(round(abs(glon-slon),digits=3)))[1] 
    l[i]=which(round(abs(glat-slat),digits=3)==min(round(abs(glat-slat),digits=3)))[1]
  }
  m <- list(xgridnum=k, ygridnum=l)
  invisible(m)
}  






###filter time series with custom filter:
tsfilt<-function(x,width=31,method="loess",cut.end=T){
  #possibilities
  #"loess" --> loessfilt function
  #"spline2" --> splinesmoother2 from Dave
  # "gauss" --> gaussfilter
  #"rm" --> running mean
  #"hamming" --> hamming
  #"bw" --> butterworth
  if(method=="loess"){
    filtered<-loessfilt(x,width,cut.end=cut.end)
  }
  if(method=="spline2"){
    if(length(which(is.na(x)))>0){
      if(is.null(dim(x))==T){
        filtered<-splinesmoother2.nas(x,width,cut.end=cut.end)
      }else{ 
        filtered<-ts(apply(x,2,function(y) splinesmoother2.nas(y,width,cut.end)),start=start(x)[1])
      }
    }else{
      filtered<-splinesmoother2(x,width)
      if(cut.end==T){
        sx<-fy(x)+floor(width/2)-1
        ex<-ly(x)-floor((width-0.5)/2)+1
        if(is.null(dim(x))==T){
          window(filtered,end=sx)<-NA
          window(filtered,start=ex)<-NA 
        }else{
          for (i in 1:dim(x)[2]){
            window(filtered[,i],end=sx[i])<-NA
            window(filtered[,i],start=ex[i])<-NA
          }
        }
      }
    }
  }
  if(method=="gauss"){
    if(is.null(dim(x))==T){
      filtered<-gauss.na(x,width)
    }else{ 
      filtered<-ts(apply(x,2,function(y) gauss.na(y,width)),start=start(x)[1])
    }
  }
  if(method=="rm"){
    if(is.null(dim(x))==T){
      filtered<-rollmean.na(x,width)
    }else{ 
      if(length(which(is.na(x)))>0){
        filtered<-ts(apply(x,2,function(y) rollmean.na(y,width)),start=start(x)[1])
      }else{
        filtered<-rollmean(x,width)
      }
    } 
  }
  if(method=="hamming"){
    require(oce)
    filtered<-stats::filter(x,makeFilter("hamming", 50, asKernel=FALSE))
  }
  
  
  if(method=="bw"){
    if(is.null(dim(x))==T){
      filtered<-butterfilt.na(x,width)
    }else{ 
      filtered<-tsapply(x,2,function(y) butterfilt.na(y,width))
    }
    if(cut.end==T){
      sx<-start(x)[1]+floor(width/2)-1
      ex<-end(x)[1]-floor((width-0.5)/2)+1
      window(filtered,end=sx)<-NA
      window(filtered,start=ex)<-NA
    }
  }
  filtered
}


##apply a function to a timeseries-matrix. the outcome is a time series with the same tsp
##fun should look something like this: function(x) mean(x,na.rm=T)
tsapply<-function(x,dim,fun){
  out<-ts(apply(x,dim,fun),start=start(x)[1])
}

#butterworth filter removing missing values
butterfilt.na<-function(y,tsc,type="low",order=2){
  #require(signal)
  # dt is expected to be equal to 1 sampling unit!!
  #first get rid of NAs at begining and nend
  sy<-min(which(!is.na(y)))
  ey<-max(which(!is.na(y)))
  x<-y[sy:ey]
  mx<-mean(x)
  #anomalies to full peroid mean may help to reduce end effects in some cases
  x<-x-mx
  nx<-tsc*order
  x2<-c(rep(mean(x[1:tsc]),nx),x,rep(mean(x[(length(x)-tsc+1):length(x)]),nx))
  bf <- signal:::butter(order, 1/tsc, type=type)
  b1 <- signal:::filtfilt(bf, x2)
  b1 <- b1[-c(1:nx,(length(b1)-nx+1):length(b1))]
  #pl.mts(cbind(x2,b1))
  b1<-b1+mx
  z<-y
  z[sy:ey]<-b1
  return(z)
  #detach(package:signal,unload=TRUE,force=TRUE)
  #unloadNamespace("signal")
  #  freqz(b1)
  #  zplane(bf)
}

