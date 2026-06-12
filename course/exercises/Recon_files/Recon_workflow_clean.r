ptm<-proc.time()

################################################################################
## Part I: Setup and input preparation
################################################################################

## 1. Validate option combinations ----------------------------------------
validate.workflow.options<-function(env=parent.frame()){
  get.opt<-function(name) get(name,envir=env,inherits=FALSE)
  set.opt<-function(name,value) assign(name,value,envir=env)
  is.true<-function(name) isTRUE(get.opt(name))
  disable.opt<-function(name,reason){
    if(isTRUE(get.opt(name))){
      set.opt(name,FALSE)
      message(sprintf("Disabling %s: %s",name,reason))
    }
  }
  
  if(is.true("do.cps") && is.true("do.cca")){
    stop("Only one reconstruction method can be active. Set either do.cps or do.cca to TRUE, not both.")
  }
  if(!is.true("do.field") && !is.true("do.index")){
    stop("Nothing to reconstruct. Set do.field and/or do.index to TRUE.")
  }
  if(get.opt("calib.start") > get.opt("calib.end")){
    stop("Invalid calibration period: calib.start must be smaller than or equal to calib.end.")
  }
  if(get.opt("startyear") > get.opt("endyear")){
    stop("Invalid reconstruction period: startyear must be smaller than or equal to endyear.")
  }
  if(is.true("do.verif.early") || is.true("do.ens.scores")){
    if(get.opt("early.start") > get.opt("early.end")){
      stop("Invalid early verification period: early.start must be smaller than or equal to early.end.")
    }
    if(get.opt("early.end") >= get.opt("calib.start")){
      stop("The early verification period must end before the calibration period starts.")
    }
  }
  if(is.true("do.ens.scores") && get.opt("nens") < 2){
    stop("do.ens.scores requires nens >= 2.")
  }
  if(is.true("do.cca")){
    if(!is.true("do.field") || is.true("do.index")){
      stop("CCA currently supports field reconstructions only. Use do.field=TRUE and do.index=FALSE.")
    }
    if(isTRUE(get.opt("cca_options")$loadparas) && isTRUE(get.opt("cca_options")$saveparas)){
      stop("cca_options$loadparas and cca_options$saveparas cannot both be TRUE in the same run.")
    }
    if(isTRUE(get.opt("cca_options")$saveparas) && get.opt("nens") != 1){
      stop("cca_options$saveparas should be used with nens == 1 so each model writes one parameter file.")
    }
    if(isTRUE(get.opt("cca_options")$loadparas)){
      parafileroot<-get.opt("cca_options")$parafileroot
      if(is.null(parafileroot) || !nzchar(parafileroot)){
        stop("cca_options$loadparas requires a non-empty cca_options$parafileroot.")
      }
      if(length(Sys.glob(paste0(parafileroot,"*.mat"))) == 0){
        stop(sprintf("No CCA parameter files found for cca_options$parafileroot='%s'.",parafileroot))
      }
    }
  }
  if(is.true("npc.fix") && is.true("sample.pcs")){
    stop("npc.fix and sample.pcs cannot both be TRUE. Disable PC sampling when using a fixed number of PCs.")
  }
  if(is.true("do.residual.ar1") && !is.true("do.residuals")){
    stop("do.residual.ar1 requires do.residuals == TRUE.")
  }
  if(is.true("write.residuals.all") && !is.true("do.residuals")){
    stop("write.residuals.all requires do.residuals == TRUE.")
  }
  if(is.true("add.arnoise") && get.opt("MCiterations") > 1 && is.true("dopar")){
    stop("MCiterations > 1 with add.arnoise == TRUE currently requires dopar == FALSE. The noise-only ensemble reuse relies on serial ensemble execution.")
  }
  if(is.true("do.calibvar")){
    if(get.opt("mincl") < 1 || get.opt("mincl") > get.opt("calib.length")){
      stop("Invalid calibration sampling setup: mincl must be between 1 and calib.length.")
    }
    if(get.opt("maxcl") < get.opt("mincl") || get.opt("maxcl") > get.opt("calib.length")){
      stop("Invalid calibration sampling setup: maxcl must be between mincl and calib.length.")
    }
    if(get.opt("verif.block.length") < 1){
      stop("verif.block.length must be at least 1.")
    }
  } else {
    calib.years.fix.values<-get.opt("calib.years.fix")
    if(length(calib.years.fix.values) == 0){
      stop("calib.years.fix must contain at least one calibration year index when do.calibvar == FALSE.")
    }
    if(anyDuplicated(calib.years.fix.values)){
      stop("calib.years.fix must not contain duplicate indices.")
    }
    if(any(calib.years.fix.values < 1 | calib.years.fix.values > get.opt("calib.length"))){
      stop("calib.years.fix must only contain indices between 1 and calib.length.")
    }
  }
  
  if(!is.true("do.field")) disable.opt("members.out","members are only written for field reconstructions.")
  if(!is.true("do.field")) disable.opt("do.ens.mean.nc","ensemble-mean NetCDF output is only written for field reconstructions.")
  if(is.true("write.quantiles.recon") && is.true("bigjob")){
    set.opt("write.quantiles.recon",FALSE)
    message("Disabling write.quantiles.recon: quantile writing is not supported when bigjob == TRUE.")
  }
  if(is.true("do.proz.pos") && (!is.true("do.field") || !is.true("do.re"))){
    set.opt("do.proz.pos",FALSE)
    message("Disabling do.proz.pos: positive-RE percentages require do.field == TRUE and do.re == TRUE.")
  }
  if(is.true("sample.radius") && !is.true("do.radius")){
    set.opt("sample.radius",FALSE)
    message("Disabling sample.radius: radius sampling requires do.radius == TRUE.")
  }
  if(is.true("few.na") && !is.true("do.radius")){
    set.opt("few.na",FALSE)
    message("Disabling few.na: the radius fallback is only used when do.radius == TRUE.")
  }
  if(is.true("sample.cps.weight") && !is.true("do.cps.cor.weight")){
    set.opt("sample.cps.weight",FALSE)
    message("Disabling sample.cps.weight: CPS weight sampling requires do.cps.cor.weight == TRUE.")
  }
  
  invisible(NULL)
}

## 2. Runtime defaults and package setup ----------------------------------
if(!exists("bigjob")) bigjob<-F
if(!exists("write.nests")) write.nests<-F
if(!exists("do.filter")) do.filter<-F
if(!exists("filterlength")) filterlength<-30
if(!exists("recon.files.path")) recon.files.path<-""
if(nzchar(recon.files.path) && !grepl("[/\\\\]$",recon.files.path)){
  recon.files.path<-paste0(recon.files.path,"/")
}
validate.workflow.options(environment())

library(foreach)
library(doParallel)
library(parallel)
if(add.arnoise) library(abind)

## 3. Workflow templates and sourced code ---------------------------------

template.file<-paste0(recon.files.path,"template.r")
template.file.cps<-paste0(recon.files.path,"template_CPS.r")
template.file.cps.field<-paste0(recon.files.path,"template_CPS_field.r")
template.file.cca<-paste0(recon.files.path,"template_CCA.R")
template.file.verif<-paste0(recon.files.path,"template_verif.r")
end.template.file<-paste0(recon.files.path,"end.template.r")
ncwrite.template.file<-paste0(recon.files.path,"template_ncwrite_field.r")

source(template.file)
source(template.file.cps)
source(template.file.cps.field)
source(template.file.cca)
source(template.file.verif)
source(end.template.file)
source(ncwrite.template.file)

if(do.cps==T){
  cps.function.file<-paste0(recon.files.path,"CPS_recon_function.r")
}
if(do.cca){
  source(paste0(recon.files.path,"cca_bp.r"))
  source(paste0(recon.files.path,"cca_cv.r"))
}
## 4. Load shared helper routines ----------------------------------------
source(paste0(recon.files.path,"R-Functions_recons.r"))
if(add.arnoise==T){
  if(arnoise.version=="simple"){
    source(paste0(recon.files.path,"add_residualnoise.r"))
  }
  if(arnoise.version=="gene"){
    source(paste0(recon.files.path,"add_residualnoise_gene.r"))
  }
  if(add.arnoise.spat){
  source(paste0(recon.files.path,"add_residualnoise_spatial_cov.r"))
  }
}

## 5. Align ensemble size for pure AR-noise experiments -------------------
if(add.arnoise==T & MCiterations>1){
  do.rmse<-F
  do.rmse.all<-F
  do.proz.pos<-F
  do.res.calib.verif.mean<-F
  do.residuals<-F
  do.residual.ar1<-F
  write.residuals.all<-F
  nens<-MCiterations
}


## 6. Output file paths and run log ----------------------------------------
if(file.exists(out.suffix)==F) dir.create(out.suffix,"0775")
output.paths<-build.output.paths(
  out.suffix=out.suffix,
  do.verif.early=do.verif.early,
  do.ens.scores=do.ens.scores,
  do.field=do.field,
  do.index=do.index,
  do.residuals=do.residuals,
  do.residual.ar1=do.residual.ar1
)
list2env(output.paths,envir=environment())
rm(output.paths)

## 6a. Write runtime options -----------------------------------------------
write.options(optionsfile)


## 7. Common reconstruction time axis -------------------------------------
firstdate<-paste(startyear,"-01-01",sep="")
tim<-as.Date(firstdate, "%Y-%m-%d")
for (t in (startyear+1):endyear){
tim<-c(tim,as.Date(paste(t,"01","01",sep="-"), "%Y-%m-%d"))
}
times<-as.numeric(tim)

reconyears<-startyear:endyear

library(ncdf4)

## 8. Load and preprocess target data -----------------------------------

if(read.instr==TRUE){
  onecol<-F
  if(ncfile==T){
    nc <- nc_open(targetfile)
    grid <- ncvar_get(nc,nc.var)
    lons<-ncvar_get(nc,nc.lon)
    lats<-ncvar_get(nc,nc.lat)
    nc_close(nc)
    if(length(dim(grid))==3){
      grid.dims<-dim(grid)[1:2]
      grid.2d<-array(grid,dim=c(grid.dims[1]*grid.dims[2],dim(grid)[3]))
    }else{
      if(length(dim(grid))==2){
        grid.2d<-t(grid)
      }else{
        if (do.field==TRUE) stop ("Target file does not contain a grid. Field recon not possible")
        onecol<-T
        target<-as.vector(grid)
      }
    }
    if(onecol==F){
      na.cells<-apply(grid.2d,1,function(x) (length(which(!is.na(x[-1])))==dim(grid.2d)[2]-1))
      nona<-which(na.cells==T & apply(grid.2d,1,sd,na.rm=T)>0.1)
      target<-t(grid.2d[nona,])
   } 
    target<-ts(target,start=nc.startyear)
    if(do.verif.early==T | do.ens.scores==T) target.early<-window(target,start=early.start,end=early.end)
    target<-window(target,start=calib.start,end=calib.end)
  }else{
    target<-as.matrix(read.table(targetfile,sep=targetfile.sep,header=targetfile.head))
    if(do.verif.early==T  | do.ens.scores==T) target.early<-target[which(target[,1]<=early.end & target[,1]>=early.start ),-1]
    target<-target[which(target[,1]<=calib.end),]
    target<-target[which(target[,1]>=calib.start),]
    if(file.exists(nonafile)==T) nona<-as.matrix(read.table(nonafile))
    grid.dims<-c(length(lons),length(lats),dim(target)[1])
    nc.varname<-"var1"
    if(dim(target)[2]==2){
      if (do.field==TRUE) stop ("Target file does not contain a grid. Field recon not possible")
      onecol<-T
      lons<-1
      lats<-1
    }
    target<-ts(target[,-1],start=target[1,1])
  }
}else{
  
  onecol<-F
  
  targetwindow<-which(targetdata$time>=calib.start & targetdata$time<=calib.end)
  targetdata$data <- targetdata$data[,,targetwindow]
  grid.dims<-dim(targetdata$data)[1:2]
  grid.2d<-array(targetdata$data,dim=c(grid.dims[1]*grid.dims[2],dim(targetdata$data)[3]))
  
  #na.cells<-apply(grid.2d,1,function(x) (length(which(!is.na(x[-1])))==dim(grid.2d)[2]-1))
  na.cells<-apply(grid.2d,1,function(x) (length(which(!is.na(x)))>2)) #more than 1 value needed for infilling...
  nona<-which(na.cells==T) # & apply(grid.2d,1,sd,na.rm=T)>0.1)
  cru<-t(grid.2d[nona,])
  
  
  #fill in missing data
  if (length(which(is.na(cru)))>0){
    print("Instrumental data need to be infilled!")
    print(paste0("They have ",length(which(is.na(cru)))," mising values"))
    print(paste0("This is ",round(length(which(is.na(cru)))/length(cru)*100,0),"% of the data points"))
    Sys.sleep(3)
    
    proxy.table.calib <- cbind(calib.start:calib.end,cru)
    
    source(paste0(recon.files.path,"infill_proxies_EOF_CPS_PCR_LMS_DINEOF_for_recon.r"))
    
    cru<-data.out[,-1]
    
    print("...instrumental data infilling finished")
  }  
  
  target <- cru
  lons<-targetdata$lon
  lats<-targetdata$lat
  
  if(do.verif.early==T  | do.ens.scores==T) target.early<-window(target,early.start,early.end)
  target<-window(target,calib.start,calib.end)
}

if(onecol==T & do.index==T){
  target.mean<-target
  if(do.verif.early==T | do.ens.scores==T)  target.early.mean<-target.early
}

if(length(dim(lons))==2){
  lons.2d<-as.vector(lons)
  lats.2d<-as.vector(lats)
  lons<-unique(lons.2d)
  lats<-unique(lats.2d)
}else{
  lons.2d<-rep(lons,length(lats))
  lats.2d<-rep(lats,each=length(lons))
}


## 7a. Create the latitude-weighted target mean ---------------------------
if(onecol==F & do.index==T){
  latweights<-cos(lats.2d*pi/180)^latweight.cos.power
  if(length(nona)>1) latweights<-latweights[nona]
  target.w<-target*rep(latweights,each=dim(target)[1])
  target.mean<-ts(apply(target.w,1,sum)/sum(latweights),start=start(target.w)[1])
  wmeanout<-target.mean
  if(do.verif.early==T| do.ens.scores==T){
    target.early.w<-target.early*rep(latweights,each=dim(target.early)[1])
    target.early.mean<-ts(apply(target.early.w,1,sum)/sum(latweights),start=start(target.early.w)[1])
    wmeanout<-tsapply(cbind(target.mean,target.early.mean),1,function(x) mean(x,na.rm=T))
  }
  write.table(cbind(time(wmeanout),wmeanout),paste(out.suffix,"/Target_field_mean_latw_",out.suffix,".txt",sep=""),quote=F,sep=";",row.names=F,col.names=c("Year","lat weighted target field mean"))
}


target.orig<-target

## 8. Create requested output NetCDF containers ---------------------------
print("creating output netcdf files...")
if(bigjob==T & members.out==T & do.field==T){
  members.dir<-paste(out.suffix,"/members",sep="")
  if(file.exists(members.dir)==F) dir.create(members.dir,"0775")
  nc.outfiles.recon<-paste(out.suffix,"/members/output_",reconyears,"_",out.suffix,".nc",sep="")
  for(y in seq_along(reconyears)){
    if(file.exists(nc.outfiles.recon[y])==F) make.nc.nc4(nc.outfiles.recon[y],lons,lats,(1:nens),times[y],nc.varname,nc.varunit)
  }
}
if(members.out==T & do.field==T){
  if(bigjob==F){
    nlevels<-ifelse(write.quantiles.recon==T,length(quantiles),nens)
    if(file.exists(nc.outfile.recon)==F) make.nc.nc4(nc.outfile.recon,lons,lats,(1:nlevels),times,nc.varname,nc.varunit)
    out.nc<-nc_open(nc.outfile.recon,write=T)
  }
  if(do.verif==T){
    nlevels.verif<-ifelse(write.quantiles.verif==T,length(quantiles),nens)
    if(do.re==T){
      if(file.exists(nc.outfile.res)==F) make.nc.nc4(nc.outfile.res,lons,lats,(1:nlevels.verif),times,paste(nc.varname,"_REs",sep=""),"unitless")
      res.nc<-nc_open(nc.outfile.res,write=T)
    }
    if(do.ce==T){
      if(file.exists(nc.outfile.ces)==F) make.nc.nc4(nc.outfile.ces,lons,lats,(1:nlevels.verif),times,paste(nc.varname,"_CEs",sep=""),"unitless")
      ces.nc<-nc_open(nc.outfile.ces,write=T)
    }
    if(do.r2==T){
      if(file.exists(nc.outfile.r2)==F) make.nc.nc4(nc.outfile.r2,lons,lats,(1:nlevels.verif),times,paste(nc.varname,"_r2",sep=""),"unitless")
      r2.nc<-nc_open(nc.outfile.r2,write=T)
    }
  }
}
if(do.ens.mean.nc==T & do.field==T){
  if(file.exists(nc.outfile.ens.mean)==F) make.nc.nc4(nc.outfile.ens.mean,lons,lats,1,times,paste(nc.varname,"_Ens_mean",sep=""),nc.varunit)
  ens.mean.nc<-nc_open(nc.outfile.ens.mean,write=T)
  if(do.verif==T){
    if(do.re==T){
      if(file.exists(nc.outfile.ens.mean.res)==F) make.nc.nc4(nc.outfile.ens.mean.res,lons,lats,1,times,paste(nc.varname,"_Ens_mean_REs",sep=""),"unitless")
      ens.mean.res.nc<-nc_open(nc.outfile.ens.mean.res,write=T)
    }
    if(do.ce==T){
      if(file.exists(nc.outfile.ens.mean.ces)==F) make.nc.nc4(nc.outfile.ens.mean.ces,lons,lats,1,times,paste(nc.varname,"_Ens_mean_CEs",sep=""),"unitless")
      ens.mean.ces.nc<-nc_open(nc.outfile.ens.mean.ces,write=T)
    }
    if(do.r2==T){
      if(file.exists(nc.outfile.ens.mean.r2)==F) make.nc.nc4(nc.outfile.ens.mean.r2,lons,lats,1,times,paste(nc.varname,"_Ens_mean_r2",sep=""),"unitless")
      ens.mean.r2.nc<-nc_open(nc.outfile.ens.mean.r2,write=T)
    }
  }
}
if(do.rmse==T & do.field==T){
  if(file.exists(nc.outfile.rmse)==F) make.nc.nc4(nc.outfile.rmse,lons,lats,1,times,paste(nc.varname,"_RMSE",sep=""),nc.varunit)
  rmse.nc<-nc_open(nc.outfile.rmse,write=T)
  if(file.exists(nc.outfile.rmse.scaled)==F) make.nc.nc4(nc.outfile.rmse.scaled,lons,lats,1,times,paste(nc.varname,"_RMSE",sep=""),nc.varunit)
  rmse.scaled.nc<-nc_open(nc.outfile.rmse.scaled,write=T)
}
if(do.verif.early==T & do.field==T){
  if(file.exists(nc.outfile.res.early)==F) make.nc.nc4(nc.outfile.res.early,lons,lats,1,times,paste(nc.varname,"_Early_verif_REs",sep=""),"unitless")
  early.verif.res.nc<-nc_open(nc.outfile.res.early,write=T)
  if(file.exists(nc.outfile.ces.early)==F)  make.nc.nc4(nc.outfile.ces.early,lons,lats,1,times,paste(nc.varname,"_Early_verif_CEs",sep=""),"unitless")
  early.verif.ces.nc<-nc_open(nc.outfile.ces.early,write=T)
  if(file.exists(nc.outfile.r2.early)==F)  make.nc.nc4(nc.outfile.r2.early,lons,lats,1,times,paste(nc.varname,"_Early_verif_r2",sep=""),"unitless")
  early.verif.r2.nc<-nc_open(nc.outfile.r2.early,write=T)
  if(file.exists(nc.outfile.res.early.scaled)==F)  make.nc.nc4(nc.outfile.res.early.scaled,lons,lats,1,times,paste(nc.varname,"_Early_verif_REs_scaled",sep=""),"unitless")
  early.verif.res.nc.scaled<-nc_open(nc.outfile.res.early.scaled,write=T)
  if(file.exists(nc.outfile.ces.early.scaled)==F)  make.nc.nc4(nc.outfile.ces.early.scaled,lons,lats,1,times,paste(nc.varname,"_Early_verif_CEs_scaled",sep=""),"unitless")
  early.verif.ces.nc.scaled<-nc_open(nc.outfile.ces.early.scaled,write=T)
  if(file.exists(nc.outfile.r2.early.scaled)==F)  make.nc.nc4(nc.outfile.r2.early.scaled,lons,lats,1,times,paste(nc.varname,"_Early_verif_r2_scaled",sep=""),"unitless")
  early.verif.r2.nc.scaled<-nc_open(nc.outfile.r2.early.scaled,write=T)
}
if(do.ens.scores==T){
  ens.scores.names<-c("xc","rmse","CovR90","avgCRPS","potCRPS","Reli","Cov90w")
  source(paste0(recon.files.path,"ScoringRules.R"))
  if(do.field==T){
    if(file.exists(nc.outfile.ens.scores)==F)make.nc.nc4(nc.outfile.ens.scores,lons,lats,ens.scores.names,times,paste(nc.varname,"_ensemble_scores",sep=""),"various")
    ens.scores.nc<-nc_open(nc.outfile.ens.scores,write=T)
  }
}
print("...done")

## 9. Load, infill, and validate proxy data -------------------------------


if(read.proxies==T){
  proxy.table.input <- as.matrix(read.table(proxyfile, header=proxyfile.head,sep=proxyfile.sep)) 
}else{
  proxy.table.input <- cbind(proxydata$time,proxydata$data)
}

infill.start<-ifelse(do.verif.early==T | do.ens.scores==T,early.start,calib.start)
proxy.table.infill<-proxy.table.input[which(proxy.table.input[,1]>=infill.start & proxy.table.input[,1]<=calib.end),]


# in case there are missing values in the proxy data during the calibration or early verification period, fill them in

if(length(which(is.na(proxy.table.infill)))>0){
  if(fill.in.proxies){
    print("data infilling started...")
    proxy.table.calib<-proxy.table.infill
    
    source(paste0(recon.files.path,"infill_proxies_EOF_CPS_PCR_LMS_DINEOF_for_recon.r"))
    proxy.table.input[which(proxy.table.input[,1]>=infill.start & proxy.table.input[,1]<=calib.end),] <- data.out
    proxy.table.calib<-proxy.table.input[which(proxy.table.input[,1]>=calib.start & proxy.table.input[,1]<=calib.end),]
    
    print("...data infilling finished")
    
  }else{
    stop("Missing values in proxy matrix during calibration or early verification period. change period, proxies, or set fill.in.proxies to TRUE")
  }
}

proxy.table.calib<-proxy.table.input[which(proxy.table.input[,1]>=calib.start & proxy.table.input[,1]<=calib.end),]


proxy.table.full<-proxy.table.input[which(proxy.table.input[,1]<=max(endyear,calib.end)),]
proxy.startyear<-ifelse(do.verif.early==T | do.ens.scores==T,min(startyear,early.start),startyear)
proxy.table.full<-proxy.table.full[which(proxy.table.full[,1]>=proxy.startyear),]

proxy.table.recon<-proxy.table.input[which(proxy.table.input[,1]<=endyear),]
if(is.null(dim(proxy.table.recon))) proxy.table.recon<-t(as.matrix(proxy.table.recon))
proxy.table.recon<-proxy.table.recon[which(proxy.table.recon[,1]>=startyear),]
if(is.null(dim(proxy.table.recon))) proxy.table.recon<-t(as.matrix(proxy.table.recon))

recon.years <- startyear:endyear
missing.recon.years <- recon.years[!recon.years %in% proxy.table.recon[,1]]
if(length(missing.recon.years) > 0){
  stop(paste0(
    "Proxy file must contain every year in the reconstruction period (startyear:endyear). Missing year(s): ",
    paste(missing.recon.years, collapse=", ")
  ))
}

proxy.coverage.recon <- apply(proxy.table.recon[,-1,drop=FALSE],1,function(x) any(!is.na(x)))
if(any(!proxy.coverage.recon)){
  empty.proxy.years <- proxy.table.recon[which(proxy.coverage.recon==FALSE),1]
  stop(paste0(
    "Proxy file must provide at least one non-missing proxy record per year in the reconstruction period (startyear:endyear). Year(s) with zero available proxies: ",
    paste(empty.proxy.years, collapse=", ")
  ))
}

proxy.years <- proxy.table.full[,1]

proxies<-colnames(proxy.table.full)[-1]

if(do.verif.early==T | do.ens.scores==T){
 early.years<-which(proxy.table.full[,1]>=early.start & proxy.table.full[,1]<=early.end)
}


## 10. Initialize run-level storage ---------------------------------------

anzpcs.models<-numeric()

scalefactors.all<-array(dim=c(nens,dim(proxy.table.full)[2]-1),data=NA)
eof.truncations.all<-array(dim=c(nens,2),data=NA)
calib.years.all<-list()

if(do.cps==T){
  if(sample.cps.weight==T) allcpsweights<-numeric()
  if(do.field==T & sample.radius==T) allsearchradii<-numeric()
}


## 11. Define models / nests ----------------------------------------------
dir.create(paste(out.suffix,"/nests",sep=""),"0775")

print("defining nests...")
if(nest.blocks==F){
  md<-proxy.table.recon[,-1]
  if(is.null(dim(md))) md<-t(as.matrix(md))
  models<-define.models(md)
}else{
  endy.blocks<-min((endyear+1),ifelse(do.verif.early==T,early.start,calib.start))
  nocalibblock<-ifelse(endyear<ifelse(do.verif.early==T,early.start,calib.start),T,F) 
  
  models<-define.models.blocks(proxy.table.recon,nest.blocks.l,endy.blocks,calib.end,nocalibblock)
  empty.blocks<-which(vapply(models$proxies,length,integer(1))==0)
  if(length(empty.blocks)>0){
    empty.block.msg<-vapply(empty.blocks,function(i){
      block.years<-proxy.table.recon[models$years[[i]],1]
      if(length(block.years)==1){
        paste0("block ",i," (year ",block.years,")")
      }else{
        paste0("block ",i," (years ",block.years[1],"-",block.years[length(block.years)],")")
      }
    },character(1))
    stop(paste(
      "nest.blocks cannot be used with the current settings because no proxy record covers the full length of",
      paste(empty.block.msg,collapse=", "),
      "Reduce nest.blocks.l, disable nest.blocks, or inspect proxy coverage for these years."
    ))
  }
}
models$allyears<-list()
for(i in seq_along(models$proxies)){
  models$allyears[[i]]<-which(!is.na(apply(proxy.table.full[,-1],1,function(y) sum(y[models$proxies[[i]]]))))
}


if(sample.proxies==F){
  proxy.table.full.member<-proxy.table.full
  proxy.table.calib.member<-proxy.table.calib
  nona.years<-1:dim(proxy.table.full.member)[1]
}
##identify blocks of consecutive years for each model
nblocks<-numeric()
blocks<-list()
for (model in 1:models$n){
  diffs<-diff(models$years[[model]])
  jumps<-which(diffs>1)
  if(length(jumps)>0){
    jumps<-c(jumps,length(models$years[[model]]))
    nblocks[model]<-length(jumps)
    s<-1
    blocks[[model]]<-list() 
    for (b in 1:nblocks[model]) {
      blocks[[model]][[b]]<-models$years[[model]][s:jumps[b]]
      s<-jumps[b]+1
    }
  }else{
    nblocks[model]<-1
    blocks[[model]]<-models$years[[model]]
  }
}

modelsname<-paste(out.suffix,"/nests/Nest_years_",out.suffix,".txt",sep="")
proxiesname<-paste(out.suffix,"/nests/Nest_proxies_",out.suffix,".txt",sep="")
for(i in 1:models$n){
  write.table(t(models$years[[i]]+startyear-1),modelsname,quote=F,sep=";",col.name=F,row.names=F,append=ifelse(i==1,F,T))
  write.table(t(models$proxies[[i]]),proxiesname,quote=F,sep=";",col.name=F,row.names=F,append=ifelse(i==1,F,T))
}
print("...done")

## 12. Precompute shared method inputs ------------------------------------
if(do.pc.opt==T | (sample.pcs==T & npc.fix==F)){
  north.npc.all.s<-array(dim=c(models$n,nens))
  north.npc.all.q<-array(dim=c(models$n,nens))
}


## 12a. Instrumental target PCs for PCR/CCA -------------------------------
if(sample.pcs==F) pc.sample<-0 ##allow using the f.pc.north function also for non-sampled experiments also needed for proxy pcs
if(do.field==T){
    target.orig.scaled<-scale(target.orig)
    instr.center.factor <- attr(target.orig.scaled,"scaled:center")
    instr.scale.factor <- attr(target.orig.scaled,"scaled:scale")
    latweights<-cos(lats.2d*pi/180)^latweight.cos.power
	 if(length(nona)>1) latweights<-latweights[nona]
    target.orig.latw<-target.orig.scaled*rep(latweights,each=dim(target.orig)[1])
    pc.target <- prcomp(target.orig.latw)
    if(do.pc.opt==T) npc.target<-get.npc.opt(target.orig.latw,pc.opt,100)

  eigwert <- pc.target$sdev^2  
  eigwert.relat <- cumsum(eigwert)/sum(eigwert)

if(sample.pcs==F){
  if(do.pc.opt==F){
    if(npc.fix==F){
     auswahl <- eigwert.relat < proz.pc.q 
     auswahl <- 1:(length(auswahl[auswahl])+1)
   }else{
     auswahl <- 1:npc.q
   }
  }else{
    auswahl <- 1: npc.target
    north.npc.all.q[]<-npc.target
  }
 pc.target.h <- as.matrix(pc.target$x[,auswahl])
 pc.target.a <- as.matrix(pc.target$rotation[,auswahl])
  }
}

## 13. CPS-specific setup -------------------------------------------------
if(do.cps==T){
  
  source(cps.function.file)

  if(do.field==T & cps.weigh.distances==T){
    ##read the proxy coordinates
    proxy.coords<-as.matrix(read.table(coord.file,sep=";",header=T))
    lats.proxies<-proxy.coords[1,]
    lons.proxies<-proxy.coords[2,]
    
    #change lon values of proxies
    if(changelon.proxies=="a"){
      negs<-which(lons.proxies<0)
      lons.proxies[negs]<-lons.proxies[negs]+360
    }
    if(changelon.proxies=="b"){
      negs<-which(lons.proxies>180)
      lons.proxies[negs]<-lons.proxies[negs]-360
    }
    
    ##get the distances from each proxy to each grid cell
    if(do.radius==T){
      distances<-array(dim=c(dim(proxy.table.full)[2]-1,dim(target.orig)[2]))
      for(pr in 1:(dim(proxy.table.full)[2]-1)){
        for(gp in 1:dim(target.orig)[2]){
          distances[pr,gp]<-distanz(lats.proxies[pr],lons.proxies[pr],lats.2d[gp],lons.2d[gp])
        }
      }
    }
  }
  if(cps.weigh.distances==T){
    distances.proxies<-array(dim=(c(dim(proxy.table.full)[2]-1,dim(proxy.table.full)[2]-1)))
    for(i in seq_along(lons.proxies)){
      for(j in seq_along(lons.proxies)){
        if(i==j){
          distances.proxies[i,j]<-0
        }else{
          distances.proxies[i,j]<-distanz(lat.a = lats.proxies[i],lon.a = lons.proxies[i],lat.b = lats.proxies[j],lon.b = lons.proxies[j])
        }
      }
    }
    
    dpn<-distances.proxies/sum(distances.proxies)#*(dim(proxy.table.full)[2]-1)
  }
}

## 14. Allocate output containers -----------------------------------------
if(do.field==T){
  if(do.rmse==T){
    rmse<-array(dim=c(length(startyear:endyear),dim(target.orig)[2]))
    rmse.scaled<-array(dim=c(length(startyear:endyear),dim(target.orig)[2]))
  }
  if(do.residuals==T) res.sds<-array(dim=c(length(startyear:endyear),dim(target.orig)[2]))
  if(do.residual.ar1==T) res.ar1s<-array(dim=c(length(startyear:endyear),dim(target.orig)[2]))
}
if(do.index==T){
  if(do.rmse==T){
    rmse.mean<-numeric()
    rmse.mean.scaled<-numeric()
  }
  if(do.residuals==T) res.sds.cont<-numeric()
  if(do.residual.ar1==T) res.ar1s.cont<-numeric()
  all.sigma<-array(dim=c(length(startyear:endyear),nens))
}
if(do.verif==T){
  if(do.re.index==T) all.re.cont<-array(dim=c(length(startyear:endyear),nens))
  if(do.ce.index==T) all.ce.cont<-array(dim=c(length(startyear:endyear),nens))
  if(do.r2.index==T) all.r2.cont<-array(dim=c(length(startyear:endyear),nens))
}
if(do.proz.pos==T) all.proz.pos<-array(dim=c(length(startyear:endyear),nens))
if(do.res.calib.verif.mean==T){
  calib.means.re<-numeric()
  verif.means.re<-numeric()
  calib.means.ce<-numeric()
  verif.means.ce<-numeric()
  calib.means.r2<-numeric()
  verif.means.r2<-numeric()
}
if(do.verif.early==T){
  if(do.index==T){
    skill.early.cont<-numeric()
    skill.early.scaled.cont<-numeric()
    skill.early.ce.cont<-numeric()
    skill.early.ce.scaled.cont<-numeric()
    skill.early.r2.cont<-numeric()
    skill.early.r2.scaled.cont<-numeric()
  }
}
if(do.ens.scores==T){
  ens.scores.cont<-array(dim=c(length(startyear:endyear),length(ens.scores.names)))
}
if(do.rmse.all==T){
  rmse.all<-array(dim=c(length(startyear:endyear),nens))
}

if(sample.proxies==T) all.proxy.selections<-array(dim=c(nens,dim(proxy.table.full)[2]-1),data=NA)

## 15. Pre-sample ensemble-level parameters -------------------------------

print("preparing ensemble members...")

for (ens in 1:nens){

  if(sample.proxies==T & sample.proxies.nest==F){
      if(nproxies.ret<0){
        nsamp<-round((dim(proxy.table.full)[2]-1)*nproxies.ret/(-100),0)
      }else{
        nsamp<-nproxies.ret
      }
      sel<-sample(1:(dim(proxy.table.full)[2]-1),nsamp)
      all.proxy.selections[ens,(-sel)]<-1
  }
  
  if(sample.weights==T){
     if(sample.proxies==T & sample.proxies.nest==F){
     scalefactors.orig<-runif(dim(proxy.table.full)[2]-length(sel)-1,min=1/maxsc,max=1/minsc)
     scalefactors.all[ens,(-sel)]<-scalefactors.orig
      }else{
      scalefactors.orig<-runif(dim(proxy.table.full)[2]-1,min=1/maxsc,max=1/minsc)
      scalefactors.all[ens,]<-scalefactors.orig
     }
  }
  
  if(sample.pcs==T){
      if(do.pc.opt==F){
        proz.pc.s<-runif(1,min=proz.pc.min.s,max=proz.pc.max.s)
        proz.pc.q<-runif(1,min=proz.pc.min.q,max=proz.pc.max.q)
        eof.truncations.all[ens,]<-c(proz.pc.s,proz.pc.q)
      }else{
        if(do.field==T){
          ns<-npc.target+sample((-1*pc.sample):(pc.sample),1)
          if(ns<1) ns<-1
  	      if(ns>dim(pc.target$x)[2]) ns<-dim(pc.target$x)[2]
          eof.truncations.all[ens,2]<-ns
        }
      }
  }
  
  ### for variable calib.verif period lengths
  if(do.calibvar==T){
      if(calib.variable==T){
		  if(mincl!=maxcl){
			cl<-sample(mincl:maxcl,1)
		  }else{
			cl<-mincl
		  }
      }else{
      cl<-calib.sub.length
      }
      nverifblocks<-ceiling((calib.length-cl)/verif.block.length)
      if((calib.length-cl)%%verif.block.length==0){
        vblocklengths<-rep(verif.block.length,nverifblocks)
      }else{
        vblocklengths<-c(rep(verif.block.length,nverifblocks-1),(calib.length-cl)%%verif.block.length)
      }
      avyears<-1:calib.length;vy<-numeric()
      for(i in 1:nverifblocks){
        blockstart<-sample(1:length(avyears),1)
        ivy<-avyears[blockstart:(blockstart+vblocklengths[i]-1)]
        if(length(which(is.na(ivy)))>0) ivy[which(is.na(ivy))]<-avyears[1:length(which(is.na(ivy)))]
        vy<-c(vy,ivy)
        avyears<-avyears[-which(avyears%in%ivy)]
      }
      calib.years.all[[ens]]<-avyears #sort(cy)
  
    if (ens==1) write.table("Calib_period_years",calibyearsname, sep=";", quote=F, col.names=F, row.names=F)
    write.table(t(calib.years.all[[ens]]),calibyearsname, sep=";", quote=F, col.names=F, row.names=F,append=T) 
  }
  
  ##CPS parameters
  if(do.cps==T){
    if(sample.cps.weight==T){
        cps.weight.exponent<-sample(cps.weightfactors,1)
        allcpsweights[[ens]]<-cps.weight.exponent
    }
    if(do.field==T & sample.radius==T){
        search.radius<-sample(min.radius:max.radius,1)
        allsearchradii[ens]<-search.radius
    }
  }
}

print("...done")



## 16. CCA-specific setup --------------------------------------------------
if(do.cca){
  library(R.matlab)

  weights0=cos(lats.2d*pi/180)
  weights=weights0/sum(weights0)

  cca_options$weights      = weights

  # if(cca_options$loadparas==F){
  #   write.table(paste0("start foreach loop ",date()),posfile,row.names = FALSE,col.names = FALSE)
  # }

}

## 17. Parallel backend ----------------------------------------------------
done.cl<-F
if(dopar==T | dopar.cps.field==T | dopar.ens.scores==T | (do.cca==T & cca_options$loadparas==F)){
  numCores <- detectCores()
  numCores<-min(maxcores,numCores)
  if(dopar) print(c("cores used:",numCores))
  cl <- makeCluster(numCores)#,outfile=posfile)
  registerDoParallel(cl)
  done.cl<-T
}

`%op.dopar.ens%` <- if (dopar) `%dopar%` else `%do%`
`%op.dopar.cps.field%` <- if (dopar.cps.field) `%dopar%` else `%do%`
`%op.dopar.ens.scores%` <- if (dopar.ens.scores) `%dopar%` else `%do%`
`%op.dopar.addnoise%` <- if (dopar.addnoise) `%dopar%` else `%do%`

export.packs<-"ncdf4"
if(add.arnoise==T & arnoise.version=="gene") export.packs<-c(export.packs,"waveslim")
if(add.arnoise) export.packs<-c(export.packs,"abind")
if(add.arnoise.spat) export.packs<-c(export.packs,"MASS","mvnfast")
if(do.cca) export.packs<-c(export.packs,"R.matlab")
if((do.cps==T & do.field==T) | do.ens.scores==T | add.arnoise==T) export.packs<-c(export.packs,"foreach","doParallel","parallel")

################################################################################
## Part II: Loop over models and ensemble members
################################################################################

## 18. Model loop ----------------------------------------------------------
for (model in 1:models$n){
#model<-1

# Initialize model-specific containers and bookkeeping --------------------
if(do.field==T){
  ensemble.sum<-array(dim=c(length(models$allyears[[model]]),dim(target)[2]),data=0)
  if(do.verif==T){
    if(do.re==T) res.sum<-array(dim=c(length(models$years[[model]]),dim(target)[2]),data=0)
    if(do.ce==T) ces.sum<-array(dim=c(length(models$years[[model]]),dim(target)[2]),data=0)
    if(do.r2==T) r2.sum<-array(dim=c(length(models$years[[model]]),dim(target)[2]),data=0)
  }
  if(write.quantiles.recon==T){
    ensemble<-array(dim=c(nens,length(models$years[[model]]),dim(target)[2]))
  }
  if(write.quantiles.verif==T & do.verif==T){
    res.ensemble<-array(dim=c(nens,length(models$years[[model]]),dim(target)[2]))
    ces.ensemble<-array(dim=c(nens,length(models$years[[model]]),dim(target)[2]))
    r2.ensemble<-array(dim=c(nens,length(models$years[[model]]),dim(target)[2]))
  }
}

if(do.ens.scores==T){
  if (do.field==T) verif.ensemble<-array(NA,c(length(early.years),dim(target)[2],nens))
  if (do.index==T) verif.ensemble.rosm<-array(NA,c(length(early.years),nens))
}

if(do.index==T){
  ensemble.sum.rosm<-vector(length=length(models$allyears[[model]]))
  ensemble.sum.rosm[]<-0
}

if(do.field==T) predicted<-array(dim=dim(target.orig),data=0)
if(do.index==T) predicted.mean<-rep(0,length(target.mean))
if(do.field==T & do.verif.early==T) predicted.early<-array(dim=dim(target.early),data=0)
if(do.index==T & do.verif.early==T) predicted.early.cont<-rep(0,length(target.early.mean))
if(do.field==T & do.residuals==T) residuals<-matrix(NA,nens,dim(target.orig)[2])
if(do.index==T & do.residuals==T) residuals.cont<-vector(length=nens)
if(do.field==T & do.residual.ar1==T) residual.ar1s<-matrix(NA,nens,dim(target.orig)[2])
if(do.index==T & do.residual.ar1==T) residual.ar1s.cont<-vector(length=nens)
if(do.index==T & do.res.calib.verif.mean==T){
  calib.means<-array(NA,c(calib.length,nens))
  verif.means<-array(NA,c(calib.length,nens))
}

count<-0

dataRows <- models$years[[model]]
alldataRows <- models$allyears[[model]]
reconRows <- match(dataRows,alldataRows)

if (add.arnoise==T){
  #define the years that are needed
  if(do.ens.scores==T){
    early.rows<-which(alldataRows %in% early.years)
    noiseRows<-c(reconRows,early.rows)
  }else{
    noiseRows<-reconRows
  }
}

if(bigjob==T & members.out==T & do.field==T){
  for(y in seq_along(dataRows)){
    assign(paste("out.nc",y,sep="."),nc_open(nc.outfiles.recon[dataRows[y]],write=T))
  }
}

if(write.residuals.all==T){
  residuals.all.file<-paste(out.suffix,"/nests/Residuals_all_",model,"_",out.suffix,".txt",sep="")
}

if(sample.proxies.nest==T) allsels<-array(dim=c(nens,dim(proxy.table.full)[2]-1),data=NA)

# Prepare and run the ensemble-member loop --------------------------------

# Prepare variables to export to foreach workers
sizes<-c();notexp<-c();needexp<-c();rem<-c();result<-c()
sizes<-( sapply(ls(),function(x){object.size(get(x))})) 
notexp<-which(sizes>10^5)
needexp<-c("target","target.orig","pc.target.a",
           "instr.scale.factor","instr.center.factor","residuals","residual.ar1s")
if(do.cps) needexp<-c(needexp,"recon.cps")
if(do.ens.scores==T & do.field==T) needexp<-c(needexp,"verif.ensemble")
if(do.ens.scores==T & do.index==T) needexp<-c(needexp,"verif.ensemble.rosm")

rem<-na.omit(match(needexp,ls()[notexp]))
if(length(rem)>0) notexp<-notexp[-na.omit(match(needexp,ls()[notexp]))]
exportvariables<-ls()[-notexp]
exportvariables<-c(exportvariables,"exportvariables")
if(do.pc.opt==T & pc.opt=="rulen" & exists("cis")) exportvariables<-c(exportvariables,"cis")
calib.years<-list()

# Run the ensemble-member loop
  results<-foreach (ens = 1:nens,.export=exportvariables,.noexport=exportvariables,.packages = export.packs) %op.dopar.ens% {

  result<-list();count<-0
  
#for(ens in 1:nens){
#ens<-1

#set parameters -------------------------------
  #### sample proxy data
if(sample.proxies==T){
 if(sample.proxies.nest==T){
  if(allow.fullset==T & nproxies.ret==1){
    nsamp<-nproxies.ret
    sel<-sample(c(0,models$proxies[[model]]),nsamp)  
    if(sel==0){
      allsels[ens,models$proxies[[model]]]<-1
      }else{
        allsels[ens,-sel]<-1
      }
  }else{
    if(nproxies.ret<0){
      nsamp<-round(length(models$proxies[[model]])*nproxies.ret/(-100),0)
    }else{
      nsamp<-nproxies.ret
    }
    sel<-sample(models$proxies[[model]],nsamp)  
    allsels[ens,-sel]<-1
  }
 }else{
   sel<-which(is.na(all.proxy.selections[ens,]))
 }
 if(sel[1]==0 | length(sel)==0){
   proxy.table.full.member<-proxy.table.full
   proxy.table.calib.member<-proxy.table.calib
   all.proxy.selections[ens,]<-1   
   sel<-numeric()
 }else{
   proxy.table.full.member<-proxy.table.full[,-(sel+1)]
   proxy.table.calib.member<-proxy.table.calib[,-(sel+1)]
   all.proxy.selections[ens,(-sel)]<-1   
 }
}else{
sel<-numeric()
}

  

##sample weights
if(sample.weights==T){
 if(sample.proxies==T){
   if(length(sel)>0){
     scalefactors.orig<-scalefactors.all[ens,(-sel)]
   }else{
     scalefactors.orig<-scalefactors.all[ens,]
   }
 }else{
  scalefactors.orig<-scalefactors.all[ens,]
 }
}

## sample PCs
if(sample.pcs==T){
 proz.pc.s<-eof.truncations.all[ens,1]
 if (do.field==T){
   if (do.pc.opt==F){
      proz.pc.q<-eof.truncations.all[ens,2]
      auswahl <- eigwert.relat < proz.pc.q 
      auswahl <- 1:(length(auswahl[auswahl])+1)
      north.npc.all.qx<-(length(auswahl[auswahl])+1)
   }else{
     ns<-eof.truncations.all[ens,2]
     auswahl<-1:(ns)
     north.npc.all.qx<-ns
   }
      pc.target.h <- as.matrix(pc.target$x[,auswahl])
      pc.target.a <- as.matrix(pc.target$rotation[,auswahl])
      result<-c(result,list(north.npc.all.qx=north.npc.all.qx))
 }

 
}

if(do.cps==T & do.field==T & sample.radius==T) search.radius<-allsearchradii[ens]
if(do.cps==T & sample.cps.weight==T) cps.weight.exponent<-allcpsweights[[ens]]


### for variable calib.verif period lengths
if(do.calibvar==T){
  calib.years<-calib.years.all[[ens]]
}else{
  calib.years<-calib.years.fix  
}

verif.years<-(1:calib.length)[-calib.years]
result<-c(result,list(calib.years=calib.years))
result<-c(result,list(verif.years=verif.years))

# Build calibration and reconstruction matrices for this member
proxy.matrix.calib<-proxy.table.calib.member[calib.years,]
proxy.matrix.calib <- proxy.matrix.calib[,-1]

proxy.calib.window.years <- proxy.table.calib.member[calib.years,1] 
proxy.calib.years<-proxy.table.calib.member[,1]

proxy.matrix.full <-proxy.table.full.member[,-1]

# Select calibration-period target data
if(do.field==T) target<-target.orig[calib.years,]
if(do.index==T) target.mean.calib<-target.mean[calib.years]


# Select the proxies available for this model/member

if(length(sel)>0){
	dataCol <- which((1:(dim(proxy.table.full)[2]-1))[-sel] %in% models$proxies[[model]])
}else{
	dataCol <- which((1:(dim(proxy.table.full)[2]-1)) %in% models$proxies[[model]])
}

## If no proxy survives selection, skip this member -----------------------
if(length(dataCol)>0){
 
  if(sample.weights==T){
    scalefactors<-scalefactors.orig[dataCol]
  }


# Run the chosen reconstruction method ------------------------------------
  ##if only AR noise-based uncertainties, all will be done for ens==1, afterwards only fill in the members
  if(add.arnoise==T & MCiterations>1 & ens>1){
  if(do.field==T) recon<-pred.scaled.noise[ens,1:length(reconRows),]
    if(do.index==T) recon.cont<-as.vector(y.neu.pred.scaled.noise[ens,1:length(reconRows)])
  }else{
    if(do.cps==F){
      if(do.cca==F){
        run_pcr_reconstruction(environment())
      }else{
        run_cca_reconstruction(environment())
      }
    }else{
      if(do.field==T) run_cps_field_reconstruction(environment())
      if(do.index==T) run_cps_index_reconstruction(environment())
    }
  }
  compute_diagnostics(environment())

  if(do.field==T & do.verif==T){
    result<-c(result,list(pred.calib.scaled=pred.calib.scaled))
    result<-c(result,list(estimated=estimated))
  }
  if(do.index==T & do.verif==T){
    result<-c(result,list(y.neu.pred.calib.scaled=y.neu.pred.calib.scaled))
    result<-c(result,list(estimated.cont=estimated.cont))
  }
  
  if(do.verif.early==T){
    if(do.field==T) result<-c(result,list(pred.verif.scaled=pred.verif.scaled))
    if (do.index==T) result<-c(result,list(y.neu.pred.verif.scaled=y.neu.pred.verif.scaled))
  }

  result<-c(result,list(rows=rows))
  
# Run shared post-processing for this member ------------------------------
 
  result<-pack_results(environment())

# Finalize this ensemble member -------------------------------------------

 count<-1#count<-count+1
 result<-c(result,list(count=count))
 
}else{
  print("no proxies selected in this member and model!!")
}
 #cat('member',ens,'model',model,'/',models$n,date(),Sys.getpid(),'\n')
 

#}

 
 result

}


# Combine member results for this model -----------------------------------
#cat('combining parallel results...\n')

for(ens in 1:nens){
  if("pred.calib.scaled" %in% names(results[[1]])) predicted[results[[ens]]$calib.years,]<-predicted[results[[ens]]$calib.years,]+results[[ens]]$pred.calib.scaled
  if("estimated" %in% names(results[[1]])) predicted[results[[ens]]$verif.years,]<-predicted[results[[ens]]$verif.years,]+results[[ens]]$estimated
  if("y.neu.pred.calib.scaled" %in% names(results[[1]])){
    predicted.mean[results[[ens]]$calib.years]<-predicted.mean[results[[ens]]$calib.years]+results[[ens]]$y.neu.pred.calib.scaled
    if(do.res.calib.verif.mean==T) calib.means[results[[ens]]$calib.years,ens]<-results[[ens]]$y.neu.pred.calib.scaled
  }
  if("estimated.cont" %in% names(results[[1]])){
    predicted.mean[results[[ens]]$verif.years]<-predicted.mean[results[[ens]]$verif.years]+results[[ens]]$estimated.cont
    if(do.res.calib.verif.mean==T) verif.means[results[[ens]]$verif.years,ens]<-results[[ens]]$estimated.cont
  }
  if("pred.verif.scaled" %in% names(results[[1]])) predicted.early<-predicted.early+results[[ens]]$pred.verif.scaled
  if("y.neu.pred.verif.scaled" %in% names(results[[1]])) predicted.early.cont<-predicted.early.cont+results[[ens]]$y.neu.pred.verif.scaled
  if("verif.ensx" %in% names(results[[1]]) & MCiterations==1) verif.ensemble[,,ens]<-results[[ens]]$verif.ensx
  if("verif.ensx.rosm" %in% names(results[[1]]) & MCiterations==1) verif.ensemble.rosm[,ens]<-results[[ens]]$verif.ensx.rosm
  if("pred.scaled" %in% names(results[[1]])) ensemble.sum<-ensemble.sum+results[[ens]]$pred.scaled
  if("recon" %in% names(results[[1]])){
    recon<-results[[ens]]$recon
    if(write.quantiles.recon==T & bigjob==F) ensemble[ens,,]<-recon
  }
  if("y.neu.pred.scaled" %in% names(results[[1]])) ensemble.sum.rosm<-ensemble.sum.rosm+results[[ens]]$y.neu.pred.scaled
  if("recon.cont" %in% names(results[[1]])) all.sigma[models$years[[model]],ens]<-results[[ens]]$recon.cont
  if("proz.pos" %in% names(results[[1]])) all.proz.pos[models$years[[model]],ens]<-results[[ens]]$proz.pos
  if("res" %in% names(results[[1]])){
    res<-results[[ens]]$res
    res.sum<-res.sum+res
    if(write.quantiles.verif==T) res.ensemble[ens,,]<-res
  }
  if("ces" %in% names(results[[1]])){
    ces<-results[[ens]]$ces
    ces.sum<-ces.sum+ces
    if(write.quantiles.verif==T) ces.ensemble[ens,,]<-ces
  }
  if("r2" %in% names(results[[1]])){
    r2<-results[[ens]]$r2
    r2.sum<-r2.sum+r2
    if(write.quantiles.verif==T) r2.ensemble[ens,,]<-r2
  }
  if("ce.cont" %in% names(results[[1]])) all.ce.cont[models$years[[model]],ens]<-results[[ens]]$ce.cont
  if("r2.cont" %in% names(results[[1]])) all.r2.cont[models$years[[model]],ens]<-results[[ens]]$r2.cont
  if("re.cont" %in% names(results[[1]])) all.re.cont[models$years[[model]],ens]<-results[[ens]]$re.cont
  if("residualsx" %in% names(results[[1]])) residuals[ens,]<-results[[ens]]$residualsx
  if("residuals.contx" %in% names(results[[1]])) residuals.cont[ens]<-results[[ens]]$residuals.contx
  if("residual.ar1sx" %in% names(results[[1]])) residual.ar1s[ens,]<-results[[ens]]$residual.ar1sx
  if("residual.ar1s.contx" %in% names(results[[1]])) residual.ar1s.cont[ens]<-results[[ens]]$residual.ar1s.contx
  if("rmse.allx" %in% names(results[[1]])) rmse.all[models$years[[model]],ens]<-results[[ens]]$rmse.allx
  if("north.npc.all.sx" %in% names(results[[ens]])) north.npc.all.s[model,ens]<-results[[ens]]$north.npc.all.sx
  if("north.npc.all.qx" %in% names(results[[1]])) north.npc.all.q[model,ens]<-results[[ens]]$north.npc.all.qx
  count<-count+results[[ens]]$count
  if(do.field) write_field_output(environment())
}
rows<-results[[1]]$rows
if("verif.ensemblex" %in% names(results[[1]])) verif.ensemble<-results[[1]]$verif.ensemblex
if("verif.ensemble.rosmx" %in% names(results[[1]])) verif.ensemble.rosm<-results[[1]]$verif.ensemble.rosmx

#cat('...done \n')

rm(results)
gc()

# Write model-level outputs and diagnostics -------------------------------

if(do.ens.mean.nc==T & do.field==T){
  ensemble.mean<-ensemble.sum/count
  if(do.re==T) res.mean<-res.sum/count
  if(do.ce==T) ces.mean<-ces.sum/count
  if(do.r2==T) r2.mean<-r2.sum/count
  ens.mean.3d<-vecttomat(ensemble.mean,grid.dims,nona)
  if(do.re==T) res.mean.3d<-vecttomat(res.mean,grid.dims,nona)
  if(do.ce==T) ces.mean.3d<-vecttomat(ces.mean,grid.dims,nona)
  if(do.r2==T) r2.mean.3d<-vecttomat(r2.mean,grid.dims,nona)
  if(nblocks[model]==1){
    ncvar_put(ens.mean.nc,ens.mean.nc$var[[1]]$name,ens.mean.3d[,,reconRows],c(1,1,1,blocks[[model]][1]),c(length(lons),length(lats),1,length(dataRows)))
    if(do.verif==T){
      if(do.re==T) ncvar_put(ens.mean.res.nc,ens.mean.res.nc$var[[1]]$name,res.mean.3d,c(1,1,1,blocks[[model]][1]),c(length(lons),length(lats),1,length(dataRows)))
      if(do.ce==T) ncvar_put(ens.mean.ces.nc,ens.mean.ces.nc$var[[1]]$name,ces.mean.3d,c(1,1,1,blocks[[model]][1]),c(length(lons),length(lats),1,length(dataRows)))
      if(do.r2==T) ncvar_put(ens.mean.r2.nc,ens.mean.r2.nc$var[[1]]$name,r2.mean.3d,c(1,1,1,blocks[[model]][1]),c(length(lons),length(lats),1,length(dataRows)))
    }
  }else{
    for (b in 1:nblocks[model]) {
      selyears<-which(models$years[[model]] %in% blocks[[model]][[b]])
      ncvar_put(ens.mean.nc,ens.mean.nc$var[[1]]$name,ens.mean.3d[,,reconRows[selyears]],c(1,1,1,blocks[[model]][[b]][1]),c(length(lons),length(lats),1,length(blocks[[model]][[b]])))
      if(do.verif==T){
        if(do.re==T) ncvar_put(ens.mean.res.nc,ens.mean.res.nc$var[[1]]$name,res.mean.3d[,,selyears],c(1,1,1,blocks[[model]][[b]][1]),c(length(lons),length(lats),1,length(blocks[[model]][[b]])))
        if(do.ce==T) ncvar_put(ens.mean.ces.nc,ens.mean.ces.nc$var[[1]]$name,ces.mean.3d[,,selyears],c(1,1,1,blocks[[model]][[b]][1]),c(length(lons),length(lats),1,length(blocks[[model]][[b]])))
        if(do.r2==T) ncvar_put(ens.mean.r2.nc,ens.mean.r2.nc$var[[1]]$name,r2.mean.3d[,,selyears],c(1,1,1,blocks[[model]][[b]][1]),c(length(lons),length(lats),1,length(blocks[[model]][[b]])))
      }
    }
  }
}

if(do.index==T){
  ens.mean.rosm<-ensemble.sum.rosm/count
}

if(do.field==T & members.out==T){
  if(write.quantiles.recon==T & bigjob==F){
    print("getting recon quantiles...")
    ensemble.quantiles<-apply(ensemble,c(2,3),function(x) quantile(x,quantiles,na.rm=T))
    ensemble.quantiles.3d<-vecttomat.quantiles(ensemble.quantiles,nlevels,grid.dims,dim(ensemble)[2],nona=nona)
    if(nblocks[[model]]==1){
      ncvar_put(out.nc,out.nc$var[[1]]$name,ensemble.quantiles.3d,c(1,1,1,blocks[[model]][1]),c(length(lons),length(lats),nlevels,length(dataRows)))
    }else{
      for (b in 1:nblocks[model]) {
        selyears<-which(models$years[[model]] %in% blocks[[model]][[b]])
        ncvar_put(out.nc,out.nc$var[[1]]$name,ensemble.quantiles.3d[,,,selyears],c(1,1,1,blocks[[model]][[b]][1]),c(length(lons),length(lats),nlevels,length(blocks[[model]][[b]])))
      }
    }
    print("...done")
  }
  if(do.verif==T & write.quantiles.verif==T){
    print("getting verif quantiles...")
    if(do.re==T){
      res.ensemble.quantiles<-apply(res.ensemble,c(2,3),function(x) quantile(x,quantiles,na.rm=T))
      res.ensemble.quantiles.3d<-vecttomat.quantiles(res.ensemble.quantiles,nlevels.verif,grid.dims,dim(res.ensemble)[2],nona=nona)
      if(nblocks[[model]]==1){
        ncvar_put(res.nc,res.nc$var[[1]]$name,res.ensemble.quantiles.3d,c(1,1,1,blocks[[model]][1]),c(length(lons),length(lats),nlevels.verif,length(dataRows)))
      }else{
        for (b in 1:nblocks[model]) {
          selyears<-which(models$years[[model]] %in% blocks[[model]][[b]])
          ncvar_put(res.nc,res.nc$var[[1]]$name,res.ensemble.quantiles.3d[,,,selyears],c(1,1,1,blocks[[model]][[b]][1]),c(length(lons),length(lats),nlevels.verif,length(blocks[[model]][[b]])))
        }
      }
    }
    if(do.ce==T){
      ces.ensemble.quantiles<-apply(ces.ensemble,c(2,3),function(x) quantile(x,quantiles,na.rm=T))
      ces.ensemble.quantiles.3d<-vecttomat.quantiles(ces.ensemble.quantiles,nlevels.verif,grid.dims,dim(ces.ensemble)[2],nona=nona)
      if(nblocks[[model]]==1){
        ncvar_put(ces.nc,ces.nc$var[[1]]$name,ces.ensemble.quantiles.3d,c(1,1,1,blocks[[model]][1]),c(length(lons),length(lats),nlevels.verif,length(dataRows)))
      }else{
        for (b in 1:nblocks[model]) {
          selyears<-which(models$years[[model]] %in% blocks[[model]][[b]])
          ncvar_put(ces.nc,ces.nc$var[[1]]$name,ces.ensemble.quantiles.3d[,,,selyears],c(1,1,1,blocks[[model]][[b]][1]),c(length(lons),length(lats),nlevels.verif,length(blocks[[model]][[b]])))
        }
      }
    }
    if(do.r2==T){
      r2.ensemble.quantiles<-apply(r2.ensemble,c(2,3),function(x) quantile(x,quantiles,na.rm=T))
      r2.ensemble.quantiles.3d<-vecttomat.quantiles(r2.ensemble.quantiles,nlevels.verif,grid.dims,dim(r2.ensemble)[2],nona=nona)
      if(nblocks[[model]]==1){
        ncvar_put(r2.nc,r2.nc$var[[1]]$name,r2.ensemble.quantiles.3d,c(1,1,1,blocks[[model]][1]),c(length(lons),length(lats),nlevels.verif,length(dataRows)))
      }else{
        for (b in 1:nblocks[model]) {
          selyears<-which(models$years[[model]] %in% blocks[[model]][[b]])
          ncvar_put(r2.nc,r2.nc$var[[1]]$name,r2.ensemble.quantiles.3d[,,,selyears],c(1,1,1,blocks[[model]][[b]][1]),c(length(lons),length(lats),nlevels.verif,length(blocks[[model]][[b]])))
        }
      }
    }
    print("...done")
  }
}

if(do.index==T & do.res.calib.verif.mean==T){
  calib.means.mean<-apply(calib.means,1,mean,na.rm=T)
  verif.means.mean<-apply(verif.means,1,mean,na.rm=T)
  calib.means.re[models$years[[model]]]<-error.valid.res.cont(target.mean,calib.means.mean,mean(target.mean))$RE
  verif.means.re[models$years[[model]]]<-error.valid.res.cont(target.mean,verif.means.mean,mean(target.mean))$RE
  calib.means.ce[models$years[[model]]]<-error.valid.ce.cont(target.mean,calib.means.mean)$RE
  verif.means.ce[models$years[[model]]]<-error.valid.ce.cont(target.mean,verif.means.mean)$RE
  calib.means.r2[models$years[[model]]]<-error.valid.r2.cont(target.mean,calib.means.mean)$RE
  verif.means.r2[models$years[[model]]]<-error.valid.r2.cont(target.mean,verif.means.mean)$RE
}

if(do.rmse==T){
  if(do.field==T){
    predicted<-predicted/count
    modelIDe<-error.valid.rmse(target.orig,predicted)$rmse
    rmse[models$years[[model]],]<-matrix(modelIDe, ncol=dim(target)[2], nrow=dim(rows)[1], byrow=T)
    predicted.scaled<-scaletots.period(ts(predicted,start=calib.start),target.orig,calib.start,calib.end)
    modelIDe<-error.valid.rmse(target.orig,predicted.scaled)$rmse
    rmse.scaled[models$years[[model]],]<-matrix(modelIDe, ncol=dim(target)[2], nrow=dim(rows)[1], byrow=T)
  }
  if(do.index==T){
    predicted.mean<-predicted.mean/count
    rmse.mean[models$years[[model]]]<-matrix(error.valid.rmse.cont(target.mean,predicted.mean)$rmse, ncol=1, nrow=dim(rows)[1], byrow=T)
    predicted.mean.scaled<-scaletots.period(ts(predicted.mean,start=calib.start),ts(target.mean,start=calib.start),calib.start,calib.end)
    rmse.mean.scaled[models$years[[model]]]<-matrix(error.valid.rmse.cont(target.mean,predicted.mean.scaled)$rmse, ncol=1, nrow=dim(rows)[1], byrow=T)
  }
}


if(do.verif.early==T){
 if(do.field==T){
   predicted.early<-predicted.early/count
   if(calib.start-early.end==1){
    pred.tot.early<-rbind(predicted.early,predicted) 
   }else{
    pred.tot.early<-rbind(predicted.early,array(dim=c((calib.start-early.end-1),dim(target.orig)[2]),data=1),predicted) 
   }
   predicted.early.scaled<-scaletots.period(ts(pred.tot.early,start=early.start),target.orig,calib.start,calib.end)
   predicted.early.scaled<- predicted.early.scaled[1:(early.end-early.start+1),]
     modelIDe <- error.valid.res(target.early,predicted.early,target.orig)     
     modelIDe.scaled <- error.valid.res(target.early,predicted.early.scaled,target.orig)     
     re.3d<-vecttomat.1d(modelIDe$RE,grid.dims,nona)
     re.3d.scaled<-vecttomat.1d(modelIDe.scaled$RE,grid.dims,nona)
     if(onecol==T){
       re.3d<-re.3d[1]
       re.3d.scaled<-re.3d.scaled[1]
     }
     for(y in seq_along(models$years[[model]])){
       ncvar_put(early.verif.res.nc,early.verif.res.nc$var[[1]]$name,re.3d,c(1,1,1, models$years[[model]][y]),c(length(lons),length(lats),1,1))
       ncvar_put(early.verif.res.nc.scaled,early.verif.res.nc.scaled$var[[1]]$name,re.3d.scaled,c(1,1,1, models$years[[model]][y]),c(length(lons),length(lats),1,1))
     }
     modelIDe <- error.valid.ce(target.early,predicted.early)     
     modelIDe.scaled <- error.valid.ce(target.early,predicted.early.scaled)     
     ce.3d<-vecttomat.1d(modelIDe$RE,grid.dims,nona)
     ce.3d.scaled<-vecttomat.1d(modelIDe.scaled$RE,grid.dims,nona)
     if(onecol==T){
       ce.3d<-ce.3d[1]
       ce.3d.scaled<-ce.3d.scaled[1]
     }
     for(y in seq_along(models$years[[model]])){
       ncvar_put(early.verif.ces.nc,early.verif.ces.nc$var[[1]]$name,ce.3d,c(1,1,1, models$years[[model]][y]),c(length(lons),length(lats),1,1))
       ncvar_put(early.verif.ces.nc.scaled,early.verif.ces.nc.scaled$var[[1]]$name,ce.3d.scaled,c(1,1,1, models$years[[model]][y]),c(length(lons),length(lats),1,1))
     }
     modelIDe <- error.valid.r2(target.early,predicted.early)     
     modelIDe.scaled <- error.valid.r2(target.early,predicted.early.scaled)    
     r2.3d<-vecttomat.1d(modelIDe$RE,grid.dims,nona)
     r2.3d.scaled<-vecttomat.1d(modelIDe.scaled$RE,grid.dims,nona)
     if(onecol==T){
       r2.3d<-r2.3d[1]
       r2.3d.scaled<-r2.3d.scaled[1]
     }
     for(y in seq_along(models$years[[model]])){
       ncvar_put(early.verif.r2.nc,early.verif.r2.nc$var[[1]]$name,r2.3d,c(1,1,1, models$years[[model]][y]),c(length(lons),length(lats),1,1))
       ncvar_put(early.verif.r2.nc.scaled,early.verif.r2.nc.scaled$var[[1]]$name,r2.3d.scaled,c(1,1,1, models$years[[model]][y]),c(length(lons),length(lats),1,1))
     }
 }

 if(do.index){
   predicted.early.cont<-predicted.early.cont/count
   pred.tot.early.cont<-c(predicted.early.cont,predicted.mean)
   if(calib.start-early.end==1){
    pred.tot.early.cont<-c(predicted.early.cont,predicted.mean) 
   }else{
    pred.tot.early.cont<-c(predicted.early.cont,rep(1,(calib.start-early.end-1)),predicted.mean) 
   }
   predicted.early.cont.scaled<-scaletots.period(ts(pred.tot.early.cont,start=early.start),target.mean,calib.start,calib.end)
   predicted.early.cont.scaled<- predicted.early.cont.scaled[1:(early.end-early.start+1)]
     modelIDe.cont <- error.valid.res.cont(target.early.mean,predicted.early.cont,target.mean) 
     modelIDe.cont.scaled <- error.valid.res.cont(target.early.mean,predicted.early.cont.scaled,target.mean) 
     skill.early.cont[models$years[[model]]] <- matrix(modelIDe.cont$RE, ncol=1, nrow=length(models$years[[model]]), byrow=T)
     skill.early.scaled.cont[models$years[[model]]] <- matrix(modelIDe.cont.scaled$RE, ncol=1, nrow=length(models$years[[model]]), byrow=T)
     modelIDe.cont <- error.valid.ce.cont(target.early.mean,predicted.early.cont) 
     modelIDe.cont.scaled <- error.valid.ce.cont(target.early.mean,predicted.early.cont.scaled) 
     skill.early.ce.cont[models$years[[model]]] <- matrix(modelIDe.cont$RE, ncol=1, nrow=length(models$years[[model]]), byrow=T)
     skill.early.ce.scaled.cont[models$years[[model]]] <- matrix(modelIDe.cont.scaled$RE, ncol=1, nrow=length(models$years[[model]]), byrow=T)
     modelIDe.cont <- error.valid.r2.cont(target.early.mean,predicted.early.cont) 
     modelIDe.cont.scaled <- error.valid.r2.cont(target.early.mean,predicted.early.cont.scaled) 
     skill.early.r2.cont[models$years[[model]]] <- matrix(modelIDe.cont$RE, ncol=1, nrow=length(models$years[[model]]), byrow=T)
     skill.early.r2.scaled.cont[models$years[[model]]] <- matrix(modelIDe.cont.scaled$RE, ncol=1, nrow=length(models$years[[model]]), byrow=T)
 } 
}

if(do.ens.scores==T){
  print("calculate ensemble scores...")
  if(do.field==T){
    TIdx <- seq(1, length(early.years))
    XC.Mat <- matrix(NA,dim(verif.ensemble)[2],dim(verif.ensemble)[3])
    ME.Mat <- matrix(NA,dim(verif.ensemble)[2],dim(verif.ensemble)[3])
    CovRate  <- rep(NA, dim(verif.ensemble)[2])
    CovWidth <- rep(NA, dim(verif.ensemble)[2])
    CRPS.avg <- matrix(NA, dim(target)[2], 3)
    exportvariables.ens.scores<-c("verif.ensemble","target.early","TIdx")
    results.ens.scores<-foreach (locIdx = 1:dim(target)[2],.export=exportvariables.ens.scores,.noexport=exportvariables.ens.scores) %op.dopar.ens.scores% {
      result.ens.scores<-list()
      ME.loc <- (colMeans(( verif.ensemble[TIdx, locIdx,] - target.early[TIdx, locIdx] )**2) )**.5
      result.ens.scores<-c(result.ens.scores,list(ME.Mat=ME.loc))
      XC.loc <- apply( verif.ensemble[TIdx , locIdx,],2, cor, y=target.early[TIdx, locIdx], use="pair")
      result.ens.scores<-c(result.ens.scores,list(XC.Mat=XC.loc))
      CovRate.loc <- CoverageRate( verif.ensemble[ TIdx,locIdx,], target.early[ TIdx,locIdx], probs=.9)
      result.ens.scores<-c(result.ens.scores,list(CovRate=CovRate.loc))
      CovWidth.loc <- mean(diff( apply(verif.ensemble[ TIdx, locIdx,],1,quantile,probs=c(.05,.95),na.rm=T)))
      result.ens.scores<-c(result.ens.scores,list(CovWidth=CovWidth.loc))
      CRPS.loc <- unlist(avgCRPS( verif.ensemble[ TIdx,locIdx,], target.early[ TIdx,locIdx]) )
      result.ens.scores<-c(result.ens.scores,list(CRPS.avg=CRPS.loc))
    }
    for ( locIdx in 1:dim(target)[2] ){
      ME.Mat[ locIdx,]<-results.ens.scores[[locIdx]]$ME.Mat
      XC.Mat[ locIdx,]<-results.ens.scores[[locIdx]]$XC.Mat
      CovRate[ locIdx]<-results.ens.scores[[locIdx]]$CovRate
      CovWidth[ locIdx]<-results.ens.scores[[locIdx]]$CovWidth
      CRPS.avg[ locIdx,]<-results.ens.scores[[locIdx]]$CRPS.avg
    }
    ErrMeas.df <- data.frame(xc = rowMeans( XC.Mat), rmse = rowMeans( ME.Mat), CovR90=CovRate, avgCRPS=CRPS.avg[,1], potCRPS=CRPS.avg[,2], Reli=CRPS.avg[,3], Cov90w=CovWidth)
    ens.scores.3d<-vecttomat(t(ErrMeas.df),grid.dims,nona)
    if(onecol==T){
      ens.scores.3d<-array(ens.scores.3d[1,1,],dim=c(1,1,dim(ens.scores.3d)[3]))
    }
    for(y in seq_along(models$years[[model]])){
      ncvar_put(ens.scores.nc,ens.scores.nc$var[[1]]$name,ens.scores.3d,c(1,1,1, models$years[[model]][y]),c(length(lons),length(lats),length(ens.scores.names),1))
    }
  }
  if(do.index==T){
    TIdx <- seq(1, length(early.years))
    ME.Mat <- (colMeans(( verif.ensemble.rosm[TIdx,] - target.early.mean[TIdx] )**2) )**.5
    XC.Mat <- apply( verif.ensemble.rosm[TIdx ,],2, cor, y=target.early.mean[TIdx], use="pair")
    CovRate <- CoverageRate( verif.ensemble.rosm[ TIdx,], target.early.mean[ TIdx], probs=.9)
    CovWidth <- mean(diff( apply(verif.ensemble.rosm[ TIdx, ],1,quantile,probs=c(.05,.95),na.rm=T)))
    CRPS.avg <- unlist(avgCRPS( verif.ensemble.rosm[ TIdx,], target.early.mean[ TIdx]) )
    for(y in seq_along(models$years[[model]])){
      ens.scores.cont[models$years[[model]][y],] <- c(mean( XC.Mat),mean( ME.Mat),CovRate,CRPS.avg,CovWidth)
    }
  }
  print("... done")
}

if(do.field==T & do.residuals==T){
  res.sds[models$years[[model]],]<-matrix(apply(residuals,2,median), ncol=dim(target)[2], nrow=length(models$years[[model]]), byrow=T)
}
if(do.index==T & do.residuals==T){
  res.sds.cont[models$years[[model]]]<-matrix(median(residuals.cont), ncol=1, nrow=length(models$years[[model]]), byrow=T)
}

if(do.field==T & do.residual.ar1==T){
  if(nens==1) residual.ar1s<-t(as.matrix(residual.ar1s))
  res.ar1s[models$years[[model]],]<-matrix(apply(residual.ar1s,2,median), ncol=dim(target)[2], nrow=length(models$years[[model]]), byrow=T)
}
if(do.index==T & do.residual.ar1==T){
  res.ar1s.cont[models$years[[model]]]<-matrix(median(residual.ar1s.cont), ncol=1, nrow=length(models$years[[model]]), byrow=T)
}

# Write sampled proxy selections for this model ---------------------------
if(sample.proxies==T & sample.proxies.nest==T & write.proxy.selection==T){
  allselsname<-paste(out.suffix,"/nests/Selected_proxies_nest-",model,"_",out.suffix,".txt",sep="")
  write.table(cbind(proxies,t(allsels)),allselsname, sep=";", quote=F, col.names=F, row.names=F)
}

if(bigjob==T & members.out==T & do.field==T){
  for(y in seq_along(dataRows)){
    nc_close(get(paste("out.nc",y,sep=".")))
  }
}


# Finalize this model -----------------------------------------------------

time<-(proc.time()-ptm)[3]/60

#write.table(paste(round(time,2),' minuten'),posfile,append = T,row.names = FALSE,col.names = FALSE)

cat("model",model,"/",models$n,' ',time,'minuten \n')

}

if(done.cl) stopCluster(cl)


################################################################################
## Part III: Finalize and write global outputs
################################################################################

## 19. Final output writing and cleanup -----------------------------------
if(do.rmse==T){
  if(do.index==T){
    write.table(cbind(startyear:endyear,round(rmse.mean,3)),rmsename.cont,quote=F,sep=";",col.names=F,row.names=F)
    write.table(cbind(startyear:endyear,round(rmse.mean.scaled,3)),rmsename.cont.scaled,quote=F,sep=";",col.names=F,row.names=F)
  }
  if(do.field==T){
    rmse.3d<-vecttomat(rmse,grid.dims,nona)
    ncvar_put(rmse.nc,rmse.nc$var[[1]]$name,round(rmse.3d,3),c(1,1,1,1),c(length(lons),length(lats),1,dim(rmse)[1]))
    nc_close(rmse.nc)
    rmse.scaled.3d<-vecttomat(rmse.scaled,grid.dims,nona)
    ncvar_put(rmse.scaled.nc,rmse.scaled.nc$var[[1]]$name,round(rmse.scaled.3d,3),c(1,1,1,1),c(length(lons),length(lats),1,dim(rmse.scaled)[1]))
    nc_close(rmse.scaled.nc)
  }
}


if(do.verif.early==T){
  if(do.field==T){
    nc_close(early.verif.res.nc)
    nc_close(early.verif.ces.nc)
    nc_close(early.verif.r2.nc)
    nc_close(early.verif.res.nc.scaled)
    nc_close(early.verif.ces.nc.scaled)
    nc_close(early.verif.r2.nc.scaled)
  }
  if(do.index==T){
    write.table(cbind(startyear:endyear,round(skill.early.cont,3)),skill.early.cont.name.re,quote=F,sep=";",col.names=F,row.names=F) 
    write.table(cbind(startyear:endyear,round(skill.early.scaled.cont,3)),skill.early.scaled.cont.name.re,quote=F,sep=";",col.names=F,row.names=F)
    write.table(cbind(startyear:endyear,round(skill.early.ce.cont,3)),skill.early.cont.name.ce,quote=F,sep=";",col.names=F,row.names=F) 
    write.table(cbind(startyear:endyear,round(skill.early.ce.scaled.cont,3)),skill.early.scaled.cont.name.ce,quote=F,sep=";",col.names=F,row.names=F)
    write.table(cbind(startyear:endyear,round(skill.early.r2.cont,3)),skill.early.cont.name.r2,quote=F,sep=";",col.names=F,row.names=F) 
    write.table(cbind(startyear:endyear,round(skill.early.r2.scaled.cont,3)),skill.early.scaled.cont.name.r2,quote=F,sep=";",col.names=F,row.names=F)
  }
}
if(do.ens.scores==T & do.index==T){
  write.table(cbind(startyear:endyear,round(ens.scores.cont,3)),ens.scores.cont.name,quote=F,sep=";",col.names=c("Year",ens.scores.names),row.names=F)
}
if(do.ens.scores==T & do.field==T){
  nc_close(ens.scores.nc)
}

if(do.field==T & do.residuals==T){
  write.table(cbind(startyear:endyear,round(res.sds,3)),res.sdname,quote=F,sep=";",col.names=F,row.names=F)
}
if(do.field==T & do.residual.ar1==T){
  write.table(cbind(startyear:endyear,round(res.ar1s,3)),res.ar1name,quote=F,sep=";",col.names=F,row.names=F)
}

if(do.index==T){
  write.table(cbind(startyear:endyear,round(all.sigma,3)),sigmaname, sep=";", quote=F, col.names=F, row.names=F)
  if(do.residuals==T) write.table(cbind(startyear:endyear,round(res.sds.cont,3)),res.sdname.cont,quote=F,sep=";",col.names=F,row.names=F)
  if(do.residual.ar1==T) write.table(cbind(startyear:endyear,round(res.ar1s.cont,3)),res.ar1name.cont,quote=F,sep=";",col.names=F,row.names=F)
}

if(do.verif==T){
  if(do.proz.pos==T) write.table(cbind(startyear:endyear,round(all.proz.pos,2)),prozposname, sep=";", quote=F, col.names=F, row.names=F)
  if(do.re.index==T) write.table(cbind(startyear:endyear,round(all.re.cont,2)),allrecontname, sep=";", quote=F, col.names=F, row.names=F)
  if(do.ce.index==T) write.table(cbind(startyear:endyear,round(all.ce.cont,2)),allcecontname, sep=";", quote=F, col.names=F, row.names=F)
  if(do.r2.index==T) write.table(cbind(startyear:endyear,round(all.r2.cont,2)),allr2contname, sep=";", quote=F, col.names=F, row.names=F)
}

if(sample.weights==T & do.cps==F) write.table(cbind(proxies,t(round(scalefactors.all,3))),scalefactorsname, sep=";", quote=F, col.names=F, row.names=F)

if(do.cps==F){
  # if(do.pc.opt==T | (sample.pcs==T & npc.fix==F)){
  #   write.table(t(north.npc.all.s),north.name.proxies,quote=F,col.names=F,row.names=F,sep=";")
  #   write.table(t(north.npc.all.q),north.name.instr,quote=F,col.names=F,row.names=F,sep=";")
  # }
  if(sample.pcs==T) write.table(cbind(c("PC_truncation_proxies","PC_truncation_Instr"),t(round(eof.truncations.all,3))),eoftruncationsname, sep=";", quote=F, col.names=F, row.names=F)
}

if(sample.proxies==T & sample.proxies.nest==F) write.table(cbind(proxies,t(all.proxy.selections)),proxyselectionname, sep=";", quote=F, col.names=F, row.names=F)

if(do.cps==T & sample.cps.weight==T) write.table(t(allcpsweights),cpsweights.name,quote=F,col.names=F,row.names=F)
if(do.cps==T & sample.radius==T & do.field==T) write.table(t(allsearchradii),radii.name,quote=F,col.names=F,row.names=F)

if(members.out==T & bigjob==F & do.field==T) nc_close(out.nc)
if(do.ens.mean.nc==T & do.field==T) nc_close(ens.mean.nc)
if(do.verif==T){
  if(do.re==T){
    if(members.out==T & do.field==T) nc_close(res.nc)
    if(do.ens.mean.nc==T & do.field==T) nc_close(ens.mean.res.nc)
  }
  if(do.ce==T){
    if(members.out==T & do.field==T) nc_close(ces.nc)
    if(do.ens.mean.nc==T & do.field==T) nc_close(ens.mean.ces.nc)
  }
  if(do.r2==T){
    if(members.out==T & do.field==T) nc_close(r2.nc)
    if(do.ens.mean.nc==T & do.field==T) nc_close(ens.mean.r2.nc)
  }
  if(do.index==T & do.res.calib.verif.mean==T){
    write.table(cbind(startyear:endyear,round(calib.means.re,3)),calib.mean.res.name,quote=F,sep=";",row.names=F,col.names=c("Year","Calibration_years_ensemble_mean_RE"))
    write.table(cbind(startyear:endyear,round(verif.means.re,3)),verif.mean.res.name,quote=F,sep=";",row.names=F,col.names=c("Year","Verification_years_ensemble_mean_RE"))
    write.table(cbind(startyear:endyear,round(calib.means.ce,3)),calib.mean.ces.name,quote=F,sep=";",row.names=F,col.names=c("Year","Calibration_years_ensemble_mean_RE"))
    write.table(cbind(startyear:endyear,round(verif.means.ce,3)),verif.mean.ces.name,quote=F,sep=";",row.names=F,col.names=c("Year","Verification_years_ensemble_mean_RE"))
    write.table(cbind(startyear:endyear,round(calib.means.r2,3)),calib.mean.r2.name,quote=F,sep=";",row.names=F,col.names=c("Year","Calibration_years_ensemble_mean_R2"))
    write.table(cbind(startyear:endyear,round(verif.means.r2,3)),verif.mean.r2.name,quote=F,sep=";",row.names=F,col.names=c("Year","Verification_years_ensemble_mean_R2"))
  }
}
if(do.rmse.all==T) write.table(cbind(startyear:endyear,round(rmse.all,3)),rmse.all.name,quote=F,sep=";",row.names=F,col.names=F)


time<-(proc.time()-ptm)[3]/60
cat(time,' minuten','\n')


## End of workflow --------------------------------------------------------
time<-(proc.time()-ptm)[3]/60
cat(time,' minuten','\n')


















