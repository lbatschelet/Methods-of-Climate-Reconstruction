source(paste0(recon.files.path,"gap.eof.R"))
source(paste0(recon.files.path,"gap.eof_pcrecon.R"))

##some required functions ----------
scaletoperiod<-function(data,start,end){
  if(length(dim(data))>0){
    sds<-apply(window(data,start=start,end=end),2,sd,na.rm=T)
  }else{
    sds<-sd(window(data,start=start,end=end),na.rm=T)
  }
  scaled<-ts(scale(data,scale=sds,center=F),start=start(data)[1])
  if(length(dim(scaled))>0){
    ms<-apply(window(scaled,start=start,end=end),2,mean,na.rm=T)
  }else{
    ms<-mean(window(scaled,start=start,end=end),na.rm=T)
  }
  scaled<-scale(scaled,scale=F,center=ms)
  return(scaled)
}

recon.cps.infill<-function(s.dat,target.index){
  #correlate proxies with target
  cors<- cor(s.dat[,-target.index],s.dat[,target.index],use="pairwise.complete.obs")
  
  # adjust sign -------------------------------------------------------------
  s.seasons.his.pos<-s.dat[,-target.index]
  if(dim(s.dat[,-target.index])[2]>1){
    for( i in 1:dim(s.dat[,-target.index])[2]){
      s.seasons.his.pos[,i]<-s.seasons.his.pos[,i]*sign(cors[i])
    }
  }else{
    s.seasons.his.pos<-s.seasons.his.pos*sign(cors) 
  }
  
  # make the composite -----------------------------------------------------
  if(dim(s.dat[,-target.index])[2]>1){
    x.scaled<-scale(s.seasons.his.pos)
    weights.cps<-abs(cors)
    composite<-apply(x.scaled,1,weighted.mean,w=weights.cps,na.rm=T)
  }else{
    composite<-s.seasons.his.pos
  }
  
  #scale the composite ----------------------------------------------------
  s<-sd(composite)/sd(s.dat[,target.index],na.rm=T)
  composite.scaled<-composite/s
  c<-mean(composite.scaled)-mean(s.dat[,target.index],na.rm=T)
  composite.scaled<-composite.scaled-c
  
  composite.scaled
}


# prepare data -----------------
timesforinfill<-proxy.table.calib[,1]
fill.start<-timesforinfill[1]
data.in<-ts(proxy.table.calib[,-1],start=fill.start)

m.orig<-apply(window(data.in,calib.start,calib.end),2,mean,na.rm=T)
sd.orig<-apply(window(data.in,calib.start,calib.end),2,sd,na.rm=T)
data.long<-scaletoperiod(data.in,calib.start,calib.end)

data<-window(data.long,start=fill.start,end=calib.end)


# actual infilling ------------

#percentage of missing values
print(c("percentage of missing values in fill-in period:",round(length(which(is.na(data)))/length(data)*100,2)))

## EOF infilling from Scherer and Appenzeller 2006
if(infill.method=="gap.eof"){
  data.filled<-gap.eof.neu(data)
  if(dim(data)[2]>dim(data)[1]) data.filled[[1]]<-t(data.filled[[1]])
  data.filled<-data.filled[[1]]
}

## CPS
if(infill.method=="cps"){
  data.filled<-data
  for(pr in 1:dim(data)[2]){
    if(length(which(is.na(data[,pr])))>0){
      infilled<-recon.cps.infill(data,pr)
      data.filled[which(is.na(data.filled[,pr])),pr]<-infilled[which(is.na(data.filled[,pr]))]
    }
  #  print(pr)
  }
}

## PC recon infilling (not just using SVD transformations)
if(infill.method=="pc.recon"){
  data.filled<-gap.eof.neu.pcrecon(data,method="ols",eofmax=30)
  if(dim(data)[2]>dim(data)[1]) data.filled.pcr[[1]]<-t(data.filled.pcr[[1]])
  data.filled<-data.filled[[1]]
}

## robust (lms) regression infilling (s. Van Ommen Antarctic2k recon)
if(infill.method=="lms.recon"){
  data.filled<-gap.eof.neu.pcrecon(data,method="lms",eofmax=20)
  if(dim(data)[2]>dim(data)[1]) data.filled.lms[[1]]<-t(data.filled.lms[[1]])
  data.filled<-data.filled[[1]]
}

### DINEOF of Taylor et al. 2013

if(infill.method=="dineof"){
  ###library(devtools)
  ###install_github("marchtaylor/sinkr")
  library(sinkr)
  din <- dineof(data)
  data.filled <- din$Xa 
}

layout(1)

# #Plots and checks
# pal <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
# 
# image(data.filled, col=pal(100))
# 
# image(data, col=pal(100))
# image(data[-(1:90),], col=pal(100))
# image(data.filled[-(1:90),], col=pal(100))
# 
bads<-0
problems<-numeric()


#test for reasonable infilling
for(pr in 1:dim(data)[2]){
  if(length(which(is.na(data[,pr])))>0){
    mn<-mean(data.filled[,pr][which(is.na(data[,pr]))])
    mo<-mean(data.filled[,pr][-which(is.na(data[,pr]))])
    sdo<-sd(data.filled[,pr][-which(is.na(data[,pr]))])
    if(mn>(mo+2*sdo) | mn<(mo-2*sdo)){
      print(paste0("potential infilling problem with proxy",pr,": (infilled values outside 2 std.dev. range)"))
      problems<-c(problems,pr)
      bads<-bads+1
    }
  }
}
# 
# pr<-0
# pr<-pr+1
# plot(cbind(data.filled[,problems[pr]],data[,problems[pr]]),plot.type="s", col=1:2)
# #pl.mts(cbind(data.filled.dineof[,problems[pr]],data.filled.lms[[1]][,problems[pr]],data.filled.pcr[[1]][,problems[pr]],data.fill.cps[,problems[pr]],data.filled[,problems[pr]],data[,problems[pr]]))
# colnames(data)[pr]


proxies.fill<-ts(data.filled,start=start(data)[1])

#scale back to original mean and sd
data.out1<-t(t(proxies.fill)*sd.orig)
data.out1<-ts(t(t(data.out1)+m.orig),start=start(data.in)[1])

data.out<-cbind(timesforinfill,data.out1)
colnames(data.out)<-colnames(proxy.table.calib)

## Check plots
# layout(1)
# i<-0
# 
# i<-i+1
# plot(cbind(data.out1[,i],data.in[,i]),plot.type="s", col=1:2)
# print(i)


