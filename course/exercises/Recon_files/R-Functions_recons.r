## Helper functions used by the clean reconstruction workflow.
## Sections are ordered by topic to keep the active helper file navigable.

## I/O and output helpers -------------------------------------------------
## Functions for NetCDF creation, option logging, and reshaping outputs.

make.nc.nc4<-function(filename,lons,lats,levs,times,varname,varunit){
  londim <- ncdim_def("lon","degrees_east",as.double(lons)) 
  latdim <- ncdim_def("lat","degrees_north",as.double(lats))
  lev.values <- if(is.character(levs)) seq_along(levs) else as.double(levs)
  levdim <- ncdim_def("lev","ensemble_member",lev.values)
  tunits3 <- "days since 1970-01-01 00:00:00.0 -0:00"
  timedim <- ncdim_def("time",tunits3,as.double(times))
  
  tmp_def <- ncvar_def(varname,varunit,list(londim,latdim,levdim,timedim),NA,prec="single")
  
  ncout <- nc_create(filename,tmp_def,force_v4=TRUE)
  nc_close(ncout)
}  

write.options<-function(optionsfile){
  if(file.exists(optionsfile)==T) file.remove(optionsfile)
  write.table(t(c("out.suffix",out.suffix)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("read.instr",read.instr)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  if(exists("targetfile")) write.table(t(c("targetfile",targetfile)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  if(exists("ncfile")) write.table(t(c("ncfile",ncfile)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("nc.var",nc.var)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("nc.lon",nc.lon)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("nc.lat",nc.lat)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("nc.startyear",nc.startyear)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("lons",lons)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("lats",lats)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("nonafile",nonafile)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("nc.varname",nc.varname)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("grid.dims",grid.dims)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("read.proxies",read.proxies)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  if(exists("proxyfile"))write.table(t(c("proxyfile",proxyfile)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  if(proxyfile.sep=="\t"){
    write.table(t(c("proxyfile.sep","tab")),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  }else{
    write.table(t(c("proxyfile.sep",proxyfile.sep)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  }
  write.table(t(c("do.field",do.field)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("calib.start",calib.start)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("calib.end",calib.end)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("calib.length",calib.length)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.verif.early",do.verif.early)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("early.start",early.start)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("early.end",early.end)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("startyear",startyear)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("endyear",endyear)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("nens",nens)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("dopar",dopar)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("maxcores",maxcores)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("members.out",members.out)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("write.quantiles.recon",write.quantiles.recon)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("write.quantiles.verif",write.quantiles.verif)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("quantiles",quantiles)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.ens.mean.nc",do.ens.mean.nc)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("bigjob",bigjob)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("write.nests",write.nests)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.filter",do.filter)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("filterlength",filterlength)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("nest.blocks",nest.blocks)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("sample.proxies",sample.proxies)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("nproxies.ret",nproxies.ret)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("sample.proxies.nest",sample.proxies.nest)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("write.proxy.selection",write.proxy.selection)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("allow.fullset",allow.fullset)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("sample.pcs",sample.pcs)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.pc.opt",do.pc.opt)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("pc.opt",pc.opt)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("pc.sample",pc.sample)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("proz.pc.min.s",proz.pc.min.s)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("proz.pc.max.s",proz.pc.max.s)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("proz.pc.min.q",proz.pc.min.q)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("proz.pc.max.q",proz.pc.max.q)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("proz.pc.s",proz.pc.s)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("proz.pc.q",proz.pc.q)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("npc.fix",npc.fix)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("npc.s",npc.s)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("npc.q",npc.q)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("latweight.cos.power",latweight.cos.power)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("sample.weights",sample.weights)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("minsc",minsc)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("maxsc",maxsc)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.var.adj",do.var.adj)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.calibvar",do.calibvar)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("calib.variable",calib.variable)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("mincl",mincl)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("maxcl",maxcl)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("calib.sub.length",calib.sub.length)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("verif.block.length",verif.block.length)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("calib.years.fix",calib.years.fix)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("add.arnoise",add.arnoise)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("arnoise.version",arnoise.version)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("add.arnoise.spat",add.arnoise.spat)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.pc.arnoise.spat",do.pc.arnoise.spat)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("MCiterations",MCiterations)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("dopar.addnoise",dopar.addnoise)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.index",do.index)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.verif",do.verif)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.re",do.re)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.ce",do.ce)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.r2",do.r2)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.proz.pos",do.proz.pos)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.re.index",do.re.index)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.ce.index",do.ce.index)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.r2.index",do.r2.index)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.rmse",do.rmse)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.rmse.all",do.rmse.all)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.res.calib.verif.mean",do.res.calib.verif.mean)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.residuals",do.residuals)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.residual.ar1",do.residual.ar1)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("write.residuals.all",write.residuals.all)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.ens.scores",do.ens.scores)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("dopar.ens.scores",dopar.ens.scores)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.cps",do.cps)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.cps.cor.weight",do.cps.cor.weight)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("sample.cps.weight",sample.cps.weight)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("cps.weighting",cps.weighting)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("cps.weightfactors",cps.weightfactors)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("minsc.cps",minsc.cps)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("maxsc.cps",maxsc.cps)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("cps.weigh.distances",cps.weigh.distances)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("coord.file",coord.file)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("changelon.proxies",changelon.proxies)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.radius",do.radius)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("sample.radius",sample.radius)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("min.radius",min.radius)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("max.radius",max.radius)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("search.radius",search.radius)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("minproxies.radius",minproxies.radius)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("few.na",few.na)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  write.table(t(c("do.cca",do.cca)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  if(do.cca){
    write.table(t(c("cca_options$dp_max",cca_options$dp_max)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$dt_max",cca_options$dt_max)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$dcca_max",cca_options$dcca_max)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$K",cca_options$K)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$loadparas",cca_options$loadparas)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$saveparas",cca_options$saveparas)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$sample_params",cca_options$sample_params)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$dcca_sample",cca_options$dcca_sample)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$dt_sample ",cca_options$dt_sample )),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$dp_sample",cca_options$dp_sample)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
    write.table(t(c("cca_options$parafileroot",cca_options$parafileroot)),optionsfile, sep="\t",quote=F, col.names=F, row.names=F,append=T)
  }
  
}
build.output.paths<-function(out.suffix,do.verif.early,do.ens.scores,do.field,do.index,do.residuals,do.residual.ar1){
  output.paths<-list(
    nc.outfile.recon=paste(out.suffix,"/output_",out.suffix,".nc",sep=""),
    nc.outfile.res=paste(out.suffix,"/REs_",out.suffix,".nc",sep=""),
    nc.outfile.ces=paste(out.suffix,"/CEs_",out.suffix,".nc",sep=""),
    nc.outfile.r2=paste(out.suffix,"/r2_",out.suffix,".nc",sep=""),
    nc.outfile.ens.mean=paste(out.suffix,"/Ens-mean_output_",out.suffix,".nc",sep=""),
    nc.outfile.ens.mean.res=paste(out.suffix,"/Ens-mean_REs_",out.suffix,".nc",sep=""),
    nc.outfile.ens.mean.ces=paste(out.suffix,"/Ens-mean_CEs_",out.suffix,".nc",sep=""),
    nc.outfile.ens.mean.r2=paste(out.suffix,"/Ens-mean_r2_",out.suffix,".nc",sep=""),
    nc.outfile.rmse=paste(out.suffix,"/rmse_",out.suffix,".nc",sep=""),
    nc.outfile.rmse.scaled=paste(out.suffix,"/rmse_scaled_",out.suffix,".nc",sep=""),
    sigmaname=paste(out.suffix,"/Recon_spatmean_ROSM_",out.suffix,".txt",sep=""),
    allrecontname=paste(out.suffix,"/REs_spatmean_ROSM_",out.suffix,".txt",sep=""),
    allcecontname=paste(out.suffix,"/CEs_spatmean_ROSM_",out.suffix,".txt",sep=""),
    allr2contname=paste(out.suffix,"/r2_spatmean_ROSM_",out.suffix,".txt",sep=""),
    prozposname=paste(out.suffix,"/proz_pos_",out.suffix,".txt",sep=""),
    scalefactorsname=paste(out.suffix,"/Scale-factors_",out.suffix,".txt",sep=""),
    eoftruncationsname=paste(out.suffix,"/EOF-truncations_",out.suffix,".txt",sep=""),
    caliblengthsname=paste(out.suffix,"/Calib-period_lengths_",out.suffix,".txt",sep=""),
    proxyselectionname=paste(out.suffix,"/Selected_proxies_",out.suffix,".txt",sep=""),
    optionsfile=paste(out.suffix,"/Log_",out.suffix,".txt",sep=""),
    calibyearsname=paste(out.suffix,"/Calib-period_years_",out.suffix,".txt",sep=""),
    rmsename.cont=paste(out.suffix,"/rmse_spatmean_ROSM_",out.suffix,".txt",sep=""),
    rmsename.cont.scaled=paste(out.suffix,"/rmse_scaled_spatmean_ROSM_",out.suffix,".txt",sep=""),
    calib.mean.res.name=paste(out.suffix,"/Calibration_mean_REs_",out.suffix,".txt",sep=""),
    verif.mean.res.name=paste(out.suffix,"/Verification_mean_REs_",out.suffix,".txt",sep=""),
    calib.mean.ces.name=paste(out.suffix,"/Calibration_mean_CEs_",out.suffix,".txt",sep=""),
    verif.mean.ces.name=paste(out.suffix,"/Verification_mean_CEs_",out.suffix,".txt",sep=""),
    calib.mean.r2.name=paste(out.suffix,"/Calibration_mean_R2_",out.suffix,".txt",sep=""),
    verif.mean.r2.name=paste(out.suffix,"/Verification_mean_R2_",out.suffix,".txt",sep=""),
    rmse.all.name=paste(out.suffix,"/RMSEs_all_",out.suffix,".txt",sep=""),
    north.name.proxies=paste(out.suffix,"/Opt_PC_truncation_proxies_",out.suffix,".txt",sep=""),
    north.name.instr=paste(out.suffix,"/Opt_PC_truncation_instr_",out.suffix,".txt",sep=""),
    radii.name=paste(out.suffix,"/Search_radii_",out.suffix,".txt",sep=""),
    cpsweights.name=paste(out.suffix,"/CPS-weights_",out.suffix,".txt",sep=""),
    posfile=paste(out.suffix,"/foreach_out.txt",sep="")
  )
  if(do.verif.early==T){
    output.paths$skill.early.cont.name.re<-paste(out.suffix,"/RE_spatmean_early_verif_period_",out.suffix,".txt",sep="")
    output.paths$skill.early.scaled.cont.name.re<-paste(out.suffix,"/RE_spatmean_early_verif_period_scaled_",out.suffix,".txt",sep="")
    output.paths$skill.early.cont.name.ce<-paste(out.suffix,"/CE_spatmean_early_verif_period_",out.suffix,".txt",sep="")
    output.paths$skill.early.scaled.cont.name.ce<-paste(out.suffix,"/CE_spatmean_early_verif_period_scaled_",out.suffix,".txt",sep="")
    output.paths$skill.early.cont.name.r2<-paste(out.suffix,"/r2_spatmean_early_verif_period_",out.suffix,".txt",sep="")
    output.paths$skill.early.scaled.cont.name.r2<-paste(out.suffix,"/r2_spatmean_early_verif_period_scaled_",out.suffix,".txt",sep="")
    output.paths$nc.outfile.res.early<-paste(out.suffix,"/REs_early_verif_period_",out.suffix,".nc",sep="")
    output.paths$nc.outfile.ces.early<-paste(out.suffix,"/CEs_early_verif_period_",out.suffix,".nc",sep="")
    output.paths$nc.outfile.r2.early<-paste(out.suffix,"/r2_early_verif_period_",out.suffix,".nc",sep="")
    output.paths$nc.outfile.res.early.scaled<-paste(out.suffix,"/REs_early_verif_period_scaled_",out.suffix,".nc",sep="")
    output.paths$nc.outfile.ces.early.scaled<-paste(out.suffix,"/CEs_early_verif_period_scaled_",out.suffix,".nc",sep="")
    output.paths$nc.outfile.r2.early.scaled<-paste(out.suffix,"/r2_early_verif_period_scaled_",out.suffix,".nc",sep="")
  }
  if(do.ens.scores==T){
    if(do.field==T) output.paths$nc.outfile.ens.scores<-paste(out.suffix,"/ensemble_scores_",out.suffix,".nc",sep="")
    if(do.index==T) output.paths$ens.scores.cont.name<-paste(out.suffix,"/ensemble_scores_spatmean_",out.suffix,".txt",sep="")
  }
  if(do.residuals==T){
    output.paths$res.sdname<-paste(out.suffix,"/Residual_SDs_",out.suffix,".txt",sep="")
    output.paths$res.sdname.cont<-paste(out.suffix,"/Residual_SDs_spatmean",out.suffix,".txt",sep="")
  }
  if(do.residual.ar1==T){
    output.paths$res.ar1name<-paste(out.suffix,"/Residual_AR1s_",out.suffix,".txt",sep="")
    output.paths$res.ar1name.cont<-paste(out.suffix,"/Residual_AR1s_spatmean",out.suffix,".txt",sep="")
  }
  if(length(grep("&",output.paths$posfile))>0) output.paths$posfile<-gsub("&","",output.paths$posfile)
  output.paths
}
merge.result.payloads<-function(run.env,payloads){
  if(is.null(payloads)) return(invisible(NULL))
  payload.names<-names(payloads)
  if(is.null(payload.names)) return(invisible(NULL))
  keep.named<-nzchar(payload.names)
  if(!any(keep.named)) return(invisible(NULL))
  payloads<-payloads[keep.named]
  if(exists("result",envir=run.env,inherits=FALSE)){
    result.current<-get("result",envir=run.env,inherits=FALSE)
  }else{
    result.current<-list()
  }
  current.names<-names(result.current)
  if(is.null(current.names)){
    result.current<-payloads
  }else{
    keep.current<-!(current.names %in% names(payloads))
    result.current<-c(result.current[keep.current],payloads)
  }
  assign("result",result.current,envir=run.env)
  invisible(NULL)
}

## Model and nest definition helpers --------------------------------------
## Functions for defining proxy combinations and nest-block structures.

define.models<-function(s.his){
  all.combinations<-apply(s.his,1,function(x) which(!is.na(x)))
  if(identical(dim(all.combinations),rev(dim(s.his)))==T){
    anz.models<-1
    models.years<-list(1:dim(s.his)[1])
    models.proxies<-list(1:dim(s.his)[2])
  }else{
    if(is.null(dim(all.combinations))){
      combinations<-unique(all.combinations)
      anz.models<-length(combinations)
      models.proxies<-list()
      for(i in 1:anz.models){
        models.proxies[[i]]<-as.vector(combinations[[i]])
      }
      models.years<-list(length=length(models.proxies))
      for(i in 1:anz.models){
        models.years[[i]]<-which(sapply(all.combinations,function(x) identical(as.vector(x),models.proxies[[i]]))==TRUE)
      }
    }else{
      anz.models<-1
      models.years<-list(1:dim(s.his)[1])
      models.proxies<-list(all.combinations[,1])
    }
  }
  models.order<-order(as.vector(do.call("c",models.years)))
  list(proxies=models.proxies,years=models.years,order=models.order,n=anz.models)
}

#same but with pre-defined nest lengths (blocks), allowing only proxies with no NAs in each block
define.models.blocks<-function(proxy.table.recon,nest.blocks.l,endy.blocks,calib.end,nocalibblock){
  x<-proxy.table.recon[,-1,drop=FALSE]
  tim<-proxy.table.recon[,1]
  bsy<-tim[1]
  bey<-bsy+nest.blocks.l-1
  anz.models<-0
  models.proxies<-list()
  models.years<-list()
  n<-1
  while(bey<endy.blocks){
    models.years[[n]]<-(which(tim==bsy):which(tim==bey))
    block<-x[models.years[[n]],,drop=FALSE]
    block.na<-colSums(is.na(block))
    models.proxies[[n]]<-which(block.na==0)
    anz.models<-anz.models+1
    bsy<-bey+1
    bey<-bsy+nest.blocks.l-1
    n<-n+1
  }  
  #last block before calib period start may be shorter
  if(bsy<endy.blocks){
    bey<-endy.blocks-1
    models.years[[n]]<-(which(tim==bsy):which(tim==bey))
    block<-x[models.years[[n]],,drop=FALSE]
    block.na<-colSums(is.na(block))
    models.proxies[[n]]<-which(block.na==0)
    anz.models<-anz.models+1
    n<-n+1
  }
  #nest for calib preiod
  if(nocalibblock==F){
    bsy<-endy.blocks
    bey<-calib.end
    models.years[[n]]<-(which(tim==bsy):which(tim==bey))
    block<-x[models.years[[n]],,drop=FALSE]
    block.na<-colSums(is.na(block))
    models.proxies[[n]]<-which(block.na==0)
    anz.models<-anz.models+1
    n<-n+1
    
    ##now the blocks for the post-calibration reconstructions
    bsy<-calib.end+1
    bey<-bsy+nest.blocks.l-1
    while(bsy<=max(tim)){
      if(bey>max(tim)) bey<-max(tim)
      models.years[[n]]<-(which(tim==bsy):which(tim==bey))
      block<-x[models.years[[n]],,drop=FALSE]
      block.na<-colSums(is.na(block))
      models.proxies[[n]]<-which(block.na==0)
      anz.models<-anz.models+1
      n<-n+1
      bsy<-bey+1
      bey<-bsy+nest.blocks.l-1
    }
  }
  models.order<-order(as.vector(do.call("c",models.years)))
  list(proxies=models.proxies,years=models.years,order=models.order,n=anz.models)
}

## Scaling and verification helpers ---------------------------------------
## Functions for scaling series and computing verification statistics.

##scale to another timeseries over a given period
scaletots.period<-function(data,targetdata,start,end){
  sy<-start(data)[1];ey<-end(data)[1]
  if(length(dim(data))>0){
    if(length(dim(targetdata))>0){
      sdf.c<-apply(window(data,start=start,end=end),2,sd)/apply(window(targetdata,start=start,end=end),2,sd)
    }else{
      sdf.c<-apply(window(data,start=start,end=end),2,sd)/sd(window(targetdata,start=start,end=end))
    }
  }else{
    sdf.c<-sd(window(data,start=start,end=end))/sd(window(targetdata,start=start,end=end))
  }
  data.scaled<-ts(scale(data,scale=sdf.c,center=F),start=sy,end=ey)
  if(length(dim(data))>0){
    if(length(dim(targetdata))>0){
      mf.c<-apply(window(data.scaled,start=start,end=end),2,mean)-apply(window(targetdata,start=start,end=end),2,mean)
    }else{
      mf.c<-apply(window(data.scaled,start=start,end=end),2,mean)-mean(window(targetdata,start=start,end=end))
    }
  }else{
    mf.c<-mean(window(data.scaled,start=start,end=end))-mean(window(targetdata,start=start,end=end))
  }
  data.scaled<-scale(data.scaled,center=mf.c,scale=F)
  return(data.scaled)
}

sum2<-function(x){
  return (sum(x^2))}

##apply a function to a timeseries-matrix. the outcome is a time series with the same tsp
##fun should look something like this: function(x) mean(x,na.rm=T)
tsapply<-function(x,dim,fun){
out<-ts(apply(x,dim,fun),start=start(x)[1])
}

error.valid.res<-function(measured, estimated,climatology)
{
  diffvals <-  measured - estimated
       diffvals2 <- t(t(measured) - apply(climatology,2,mean))
       Red.err <- (1-(apply(diffvals,2,sum2)/apply(diffvals2,2,sum2)))
       result<-list(RE=Red.err)
return(result)
}

error.valid.ce<-function(measured, estimated)
{
  diffvals <-  measured - estimated
       diffvals2 <- t(t(measured) - apply(measured,2,mean))
       ce <- (1-(apply(diffvals,2,sum2)/apply(diffvals2,2,sum2)))
       result<-list(RE=ce)
return(result)
}

error.valid.r2<-function(measured, estimated){
  r2<-c()
  for(i in 1:dim(measured)[2]){
    r2[i]<-cor(measured[,i],estimated[,i])^2
  }
  result<-list(RE=r2)
  
}

error.valid.rmse<-function(measured, estimated)
{
  diffvals <-  measured - estimated
       absdiffvals <- abs(diffvals)
	 mse <- apply(absdiffvals^2,2,mean,na.rm=T)
       rmse <- sqrt(mse)
       result<-list(rmse=rmse)
return(result)
}

error.valid.res.cont<-function(measured, estimated, climatology)
{
  diffvals <-  measured - estimated
       diffvals2 <- t(t(measured) - mean(climatology))
       Red.err <- (1-(sum2(diffvals)/sum2(diffvals2)))
       result<-list(RE=Red.err)
return(result)
}

error.valid.ce.cont<-function(measured, estimated)
{
  diffvals <-  measured - estimated
       diffvals2 <- t(t(measured) - mean(measured))
       ce<- (1-(sum2(diffvals)/sum2(diffvals2)))
       result<-list(RE=ce)
return(result)
}

error.valid.r2.cont<-function(measured, estimated)
{
       r2<-cor(measured,estimated)^2
       result<-list(RE=r2)
return(result)
}

error.valid.rmse.cont<-function(measured, estimated)
{
  diffvals <-  measured - estimated
       absdiffvals <- abs(diffvals)
	 mse <- mean(absdiffvals^2,na.rm=T)
       rmse <- sqrt(mse)
       result<-list(rmse=rmse)
return(result)
}

## Principal component and truncation helpers -----------------------------
## Functions and lookup data for EOF/PC selection and truncation rules.

if(file.exists(paste0(recon.files.path,"CIs_for_ruleN.RData"))){
  load(paste0(recon.files.path,"CIs_for_ruleN.RData"))
  rulen.cits.ns<-20:150
  rulen.cits.ps<-2:500
}

f.pc <- function(z.mat, var.prozent)
{ pc.z <- prcomp(z.mat,center = F)
eigwert <- pc.z$sdev^2
eigwert.relat <- cumsum(eigwert)/sum(eigwert)
auswahl <- eigwert.relat < var.prozent 
auswahl <- 1:(length(auswahl[auswahl])+1)
h <- pc.z$x[,auswahl]
list(a=pc.z$rotation[,auswahl], h=h) 
}

f.pc.nfixed <- function(z.mat, npc){
  pc.z <- prcomp(z.mat,center = F)
  if(npc>dim(pc.z$x)[2]) npc<-dim(pc.z$x)[2]-1 #allow all but one PCs as a max.
  auswahl <- 1:npc
  h <- pc.z$x[,auswahl]
  list(a=pc.z$rotation[,auswahl], h=h) 
}

find.npc.northrot<-function(x){
  z<-prcomp(x,center = F)
  p <- ncol(x)
  n<-nrow(x)
  eigspc<-z$sdev^2
  eigsLo.pc <- eigspc * (1 - sqrt(2/n))
  eigsHi.pc <- eigspc * (1 + sqrt(2/n))
  eigsd<-eigsLo.pc[-length(eigsLo.pc)]-eigsHi.pc[-1]
  min(which(eigsd<0))
}

find.npc.rulen<-function(x,rep=100,scaled=T,use.loaded=T){
  z<-prcomp(x,center = F)
  p <- ncol(x)
  n<-nrow(x)
  eigwerts<-z$sdev^2
  if(use.loaded==T & exists("cis") & exists("rulen.cits.ns") & exists("rulen.cits.ps") & n %in% rulen.cits.ns & p %in% rulen.cits.ps){
    eigsN<-cis[[n]][[p]]
  }else{
    eigsN<-ruleN.pc(n,p,reps=rep)
  }
  if(scaled==T){
    ruleNok<-eigwerts>eigsN
  }else{
    eigrel<-eigwerts/sum(eigwerts)
    eigspc<-eigrel*p
    ruleNok<- eigspc>eigsN
  }
  min(which(ruleNok==F))-1
}

get.npc.opt<-function(x,pc.opt,rep=100){
  if(pc.opt %in% c("north","rulen","empirical") == F) stop ("no valid method for pc truncation provided")
  if(pc.opt=="rulen"){
    npc<-find.npc.rulen(x,rep=rep,scaled=F)
    if(npc==0) npc<-1
  }
  if(pc.opt=="north"){
    npc<-find.npc.northrot(x)
  }
  if(pc.opt=="empirical"){
    npc<-find.npc.empirical(x,rep=rep)
  }
  npc
}

f.pc.opt<-function(x,pc.opt,pc.sample=0,rep=100,use.loaded=T){
  if(pc.opt %in% c("north","rulen","empirical") == F) stop ("no valid method for pc truncation provided")
  pc.z <- prcomp(x,center = F)
  eigwert <- pc.z$sdev^2
  p <- ncol(x)
  n<-nrow(x)
  if(pc.opt=="north"){
    eigsLo.pc <- eigwert * (1 - sqrt(2/n))
    eigsHi.pc <- eigwert * (1 + sqrt(2/n))
    eigsd<-eigsLo.pc[-length(eigsLo.pc)]-eigsHi.pc[-1]
    npc<-min(which(eigsd<0))
  }
  if(pc.opt=="rulen"){
    if(use.loaded==T & exists("cis") & exists("rulen.cits.ns") & exists("rulen.cits.ps") & n %in% rulen.cits.ns & p %in% rulen.cits.ps){
      eigsN<-cis[[n]][[p]]
    }else{
      eigsN<-ruleN.pc(n,p,reps=rep)
    }
    eigrel<-eigwert/sum(eigwert)
    eigspc<-eigrel*p
    ruleNok<- eigspc>eigsN
    npc<-min(which(ruleNok==F))-1
    if(npc==0) npc<-1
  }
  if(pc.opt=="empirical"){
    npc<-find.npc.empirical(x,rep=rep)
  }
  npc<-npc+sample((-1*pc.sample):(pc.sample),1)
  if(npc<1) npc<-1
  if(npc>dim(x)[2]) npc<-dim(x)[2]
  auswahl <- 1:npc
  h <- pc.z$x[,auswahl]
  list(a=pc.z$rotation[,auswahl], h=h,npc=npc)
}

ruleN.pc<-function(n,p,reps=100){
  xdat <- rnorm(n * p * reps)
  dim(xdat) <- c(n, p, reps)
  get.eigs <- function(x) prcomp(x,center = F)$sdev^2
  eigs <- apply(xdat, 3, get.eigs)
  q.95 <- apply(eigs, 1, quantile, probs = 0.95)
  round(q.95, 3)
}

find.npc.empirical<-function(z.mat,rep){
  eigwerts<-array(dim=c(dim(z.mat)[2],rep))
  for( i in 1:rep){
    sy<-sample(1:(dim(z.mat)[1]-dim(z.mat)[1]/10),1)
    z.mat.i<-z.mat[-(sy:(sy+dim(z.mat)[1]/10-1)),]
    pc.z <- prcomp(z.mat.i,center = F)
    eigwerts[1:length(pc.z$sdev^2),i] <- pc.z$sdev^2
  }
  eigwerts.95<-apply(eigwerts,1,function(x) quantile(x,probs=c(0.025,0.975),na.rm=T))
  lastpc<-dim(eigwerts.95)[2]
  for(i in 2:dim(eigwerts.95)[2]){
    if(eigwerts.95[2,i]>eigwerts.95[1,(i-1)]){
      lastpc<-i-1
      break
    }
  }
  lastpc
}

## CCA matrix helpers -----------------------------------------------------
## Functions shared by CCA fitting, prediction, and optional residual noise.

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

repmat<-function (a, n, m = n)
{
  if (length(a) == 0)
    return(c())
  if (!is.numeric(a) && !is.complex(a))
    stop("Argument 'a' must be a numeric or complex.")
  if (is.vector(a))
    a <- matrix(a, nrow = 1, ncol = length(a))
  if (!is.numeric(n) || !is.numeric(m) || length(n) != 1 ||
      length(m) != 1)
    stop("Arguments 'n' and 'm' must be single integers.")
  n <- max(floor(n), 0)
  m <- max(floor(m), 0)
  if (n <= 0 || m <= 0)
    return(matrix(0, nrow = n, ncol = m))
  matrix(1, n, m) %x% a
}

## Statistical and grid reshaping helpers ---------------------------------
## Functions for autocorrelation, covariance, and reshaping gridded output.

fastacf<-function(x){
  if(is.null(dim(x))){
  ac<-cor(x,c(x[-1],NA),use="complete.obs")
}else{
  ac<-apply(x,2,function(y) cor(y,c(y[-1],NA),use="complete.obs"))
}
ac
}

cov2 <- function(x){
  1/(NROW(x) -1) * crossprod(scale(x, TRUE , FALSE))
}

vecttomat<-function(x,dims,nona){
  neu<-array(dim=c(dim(x)[1],dims[1]*dims[2]))
  neu[,nona]<-x
  neu<-array(neu,dim=c(dim(x)[1],dims[1],dims[2]))
  neu<-aperm(neu,c(2,3,1))
}

vecttomat.1d<-function(x,dims,nona){
  neu<-vector(length=dims[1]*dims[2]);neu[]<-NA
  neu[nona]<-x
  neu<-array(neu,dim=c(dims[1],dims[2]))
}

vecttomat.quantiles<-function(x,nlevels,grid.dims,nyears,nona){
  neu<-array(dim=c(nlevels,nyears,grid.dims[1]*grid.dims[2]))
  neu[,,nona]<-x
  neu<-array(neu,dim=c(nlevels,nyears,grid.dims[1],grid.dims[2]))
  neu<-aperm(neu,c(3,4,1,2))
}





