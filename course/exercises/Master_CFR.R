######### Master for PC-reconstruction of fields adapted after Casty, and others.
### includes options for CPS and CCA and index recons. Some options and option-combinations are largely untested / unused for a long time so may not yet / no longer work
######################################################################################################

### "clean" version for Clim.recon Blockkurs; RN 2026/06/02

## This code runs reconstuctions of climate indices (single time series) or fields using either:
##  - PCR (à la Neukom / PAGES2k 2019 based on Luterbacher, 2002; Neukom, 2010, 2011, 2014)
##  - CPS (à la Neukom / PAGES2k 2019 based on Neukom 2010; 2014)
##  - CCA (fields only, à la Neukom 2019 based on Wang 2014, Smerdon 2010)

## This file is to define the settings and parameters and it will call the files to run the reconstructions
## The first part includes highest level required definitions
## The second & third parts contain settings that can also be run by default (on your own risk)

## run the recon by saving this script and running it: source("Master_GMST_clean.R") 

## The code will genearte a sub-folder in the current folder named like the "out.suffix" variable
## this folder will contain all output files, which habe the same suffix.

## brief explanations on each variable and option, see below.
## possible default values are provided for each variable in Part II/III

## !!! Important Notes: !!! 

## Ensemble reconstruction --------------
## The script currently produces a single-member reconstruction. For quantification of uncertainties,
## you can geberate an ensemble and resample reconstruction parameters and/or add noise to the reconstruction.
## This can all be set in Part II and III below:
## set 'nens' to the number of ensemble members you want to generat (and ideally turn parallel computing on)
## to sample recon parameters turn the following varialbes to TRUE and adapt the according parameters to your preferences.
## sample.proxies, sample.pcs (PCR only), sample.weights (PCR only), do.calibvar,
## sample.cps.weight (CPS only), sample.radius (CPS only), cca_options$sample_params (CCA only)
## To add noise to reconstruction ensemble members set  'add.arnoise' to TRUe and you will get an ensemble based on the reconstruction residuals.

## Output files -------------------------
## The script generates the following output files, where X ist the name defined in 'out.suffix'.
## Depending on the settings, not all of these files will be generated.
##  - output_X.nc: NetCDF file with the gridded reconstruction. Ensemble members are in levels.
##  - Recon_spatmean_ROSM_X.txt: Text file containing the index (or spatial mean of teh grid) reconstruction.
##    Years are in rows, Ensemble members in columns, first columns are the years, semicolon separated, no header.
##  - Analogous files to the above for verification statistics.
##  - Log_X.txt: A log file that contains the parameter settings used for this recon
##  - Target_field_mean_latw_X.txt: Latitude weighted field mean of the input grid
##  - Opt_PC_truncation_proxies_X.txt: Result of the PC truncation routine for the proxy matrix. One number is provided for each proxy nest.
##  - Opt_PC_truncation_instr_X.txt: Result of the PC truncation routine for the instrumental target grid. One number is provided for each proxy nest.
##  - foreach_out.txt: Output of the parallel routine for debugging.
##  - Calib-period_years_X.txt: Indicates, which years are used for calibration, using indices (1: first year of calibration period = 'calib.start'). One line per ensemble member. Only if do.calibvar == TRUE
##  - EOF-truncations_X.txt: EOF truncation that have been used for each ensemble members. numbers indicate the % of total variance explained by the selected PCs. Only if sample.pcs == TRUE
##  - Scale-factors_X.txt: Weight of each Proxy in the PCA. One line per proxy, one number per ensemble member. Only if sample.proxies == TRUE
##  - Folder "nests":
##      - Nest_years_X.txt: Each line stands for a unique proxy combination ("nest"). Numbers indicate, for which years of the reconstruciton period, this nest is representative (i.e. the years for which this proxy combination contributes to the final recon)
##      - Nest_proxies_X.txt: Each line stands for a unique proxy combination ("nest"). Numbers indicate the proxy records (columns in the input proxy matrix), which are used to generate the reconstruction of this nest.



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ###
### ### ### ###
### ### ### ### ### ### Part I: Overview And Daily-Use Settings
### ### ### ###
### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ######

### Paths And Experiment Name

# Define the code location and the experiment name used for the output folder and output file suffixes.

# In which directory are the reconstruction R files? 
recon.files.path <- "Recon_files/"

# # experiment name used for the output folder and as suffix for all output files
# out.suffix <- "validation.comparison.clean"

#add date-time-stamp to reconstruction to make each attempt unique
date_to_add <- format(Sys.time(), "%Y_%m_%d_%H%M%S")
out.suffix <- paste0(recon.name,"_",date_to_add)


### Core Reconstruction Setup

# Choose the reconstruction method, target type, and main time windows.

# Define Reconsruction method
# Default is PCR. Set do.cps OR do.cca to TRUE if you want to use one of these other methods.
#do.cps <- FALSE
do.cca <- FALSE  #note that CCA works only if do.field==T and do.index==F (i.e. only for field recons)


# Field reconstruction? 
#do.field <- TRUE

# make a (additional) reconstruction of the spatial mean?
# Or if the target file is a single time series, reconstruct  this index.
# In case of an index-only recon set TRUE here and FALSE above
#do.index <- TRUE


# Reconstruction Start and end years
# Provide a start year (how far should the reconstructions go back (same year or later than first year in proxy file))
# At least one proxy must be available in each year of the recon (two for CCA)
#startyear <- 1840
#end year (for the recon after the calib-overlap-period towards present
#endyear <- 2002


# Calibration periods

# start and end year ofcalibration period (no missing values allowed, automatic proxy infilling option provided below) 
# in the ensemble case, calibration / validation intervalls will be sampled within this period.
#calib.start <- 1881
#calib.end <- 2000

# Main verification switch: should verification be performed?
# Then calibration / verification will be done in blocks within calib.start:calib.end
# detailed settings below
do.verif <- FALSE # default: TRUE

# should additional early verification prior to the calibration/verification window be done
# fully independent from calibration also in the ensemble case.
# Ensemble score metrics will only be calculated over this period (but independt from the value of do.verif.early, see below)
do.verif.early <- FALSE
early.start <- 1851
early.end <- 1880



# Ensemble Size
# Set the number of reconstruction ensemble members.
# Note that for fields, the size of the output files may get large.(A grid with 1'000 cells yields about 4MB file size per 1000 years and each ensemble member)
nens <- 1

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ###
### ### ### ###
### ### ### ### ### ###  Part II: Input Files And General Workflow Settings
### ### ### ###
### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ######


### Target Input Data

# Define the instrumental target data file name and structure and how the target grid is described.

##does the instrumental target need to be read in or is it already available in the workspace
read.instr <- FALSE

targetdata <- reconstruction.target

# inputfile for the instrumental target grid
# can be txt file (no missing values, years in first column, then the cells, space or tab separated)
# or NetCDF
# For an index reconstruction, provide a txt file with two columns: year and climate data
#targetfile <- "HadCRUT4.3_GraphEM_SP80_18502014_Apr-Mar_corr.nc"

# options for the case that the instrumental data are read from file:
# is the grid a netcdf file? (annual resolution, all grid cells must have either zero missing values or only missing values)
ncfile <- TRUE
# if yes provide the variable names with the data, longitudes and latitudes
nc.var <- "tas"
nc.lon <- "lon"
nc.lat <- "lat"

# if yes provide the start year of the data in the netcdf-file
nc.startyear <- 1850


##if txt file: provide the column separator
targetfile.sep <- " "
#does the file have a header?
targetfile.head <- FALSE

#if txt file and a grid: provide the longitudes and latitueds covered by the grid and the nona-cells and the dimensions of the output grid
#(the nona-cells is kind of a land-sea mask, it lists the IDs of the cells of the original grid which are not NA ord std.dev<0.1)
#in case its an index recon just use 1 for lons and lats
lons <- 1
lats <- 1
nonafile <- "nona.txt" #can be left like this if it is not requiered (index recon ord NetCDF file)
grid.dims <- c(1,1) #can be left like this if it is not requiered (index recon ord NetCDF file)

#name and unit name of the output variable to be used in the output netcdf file
nc.varname <- "tas"
nc.varunit <- "deg_C"



### Proxy Input Data

# Define the proxy data file name and format.

# Are the proxies read in from a file or already in the workspace (needs to be a matrix called "proxy.table.input" with years in first row and proxies in all other rows)
read.proxies <- FALSE

# inputfile for the proxy data (years in first column, first line header, no missing values in calibration-overlap-period)
#proxyfile <- "Input_data/N-TREND2015_data_infilled_1880_2000.csv"

#column separator
proxyfile.sep <- ";"
#does the file have a header?
proxyfile.head <- TRUE


### Parallel Computing

# Control parallel execution across ensemble members.

# should the reconstruction be done in parallel mode (across members)
dopar <- FALSE
# maximum number of cores for parallel computing across ensemble members
# note that for index recon serial is often faster as the recons are fast and allocating takes more time.
maxcores <- 4
#can be tested as follows:
# library(parallel)
# maxcores <- detectCores()-1
# print(maxcores)


### Output And Performance

# Configure output writing, file size handling, and runtime-related options.

## should NetCDF files with all ensemble members (in levels) be written?
members.out <- TRUE # default: TRUE

## Is this a big job? (e.g. more than ca. 1000 members, 500 years and 1000 grid cells and members.out==T)?
# in this case the outputfiles will be split into 1 file per year to make the process faster
# only for recon data. use write.quantiles below to save space for verification data.
bigjob <- FALSE # default: FALSE

## write ensemble range (percentiles) instead of the full ensemble to these nc files? this will save a lot of disk space.
# write.quantiles.recon works only for bigjob==F
write.quantiles.recon <- FALSE # default: FALSE
write.quantiles.verif <- FALSE # default: FALSE
quantiles <- c(0.05,0.5,0.95) # default: c(0.05,0.5,0.95)

## should NetCDF files with the ensemble mean grid be written?
# works only if do.field==TRUE
do.ens.mean.nc <- FALSE # default: TRUE

#write reconstruction output for each nest?
write.nests <- FALSE # default: FALSE


### Uncertainty And Sampling

# Configure the main ensemble, calibration, and noise-sampling settings.

# First part defines the parameter-based ensemble.



# Proxy sampling options ### ### ###
sample.proxies <- FALSE # default: TRUE
# If yes: specify number of proxies to be retained from the proxy set for each ensemble member:
# Use negative numbers to sample a fraction of proxies to be removed. Eg -10 to remove 10% of proxies
nproxies.ret <- (-10) # default: (-10)

# sample proxies for each nest or for each member?
# if TRUE nproxies.ret proxies will be removed for each nest and member
# if FALSE nproxies.ret proxies will be removed from the full proxy set for each member and then used over all nests
# i.e. if false no proxies may be removed from early nests!
sample.proxies.nest <- FALSE # default: TRUE

# If TRUE: Should the proxy selections be written out to a file for each nest (in the "nests" folder)
write.proxy.selection <- FALSE # default: FALSE

# allow no proxy to be reomved as well
# works only if nproxies.ret==1
# if TRUE x is sampled betwen 0 and nproxies. x is the column-index of the  proxy that will be removed. if it's zero, none will be removed
allow.fullset <- FALSE # default: FALSE

# End proxy samling options ### ### ###



# Calibration and Verification period options ### ### ###

# sample over calibration/verification intervals?
# The periods are sampled ramdomly from the overlap period, with blocks of fixed or variable length
do.calibvar <- FALSE # default: TRUE

## if yes:  variable length of calibration and verificatoin period blocks?
calib.variable <- FALSE # default: TRUE

## leave this unchanged, this variable is needed below
calib.length <- calib.end-calib.start+1 

# If yes (do.calibvar == T): minimum length of Calibration period
mincl <- max(calib.length-30, 30) # default: max(calib.length-30, 30) (i.e. 30 years less than overlap period, but at least 30 years)
# maximum length of Calibration period (cannot be larger than calib.length)
maxcl <- calib.length-15 # default: calib.length-15

# if not (i.e. do.calibvar == F): length of calibration period
calib.sub.length <- calib.length-30 # default: calib.length-30

# how long should the verification blocks be (verif periods will be split up in nyears/calib.block.length years, with one shorter block if there is a remainder, the startyears of the blocks will be sampled, all remaining years will be used for claibration).
verif.block.length <- 10 # default: 10

##for fixed calibration years (do.calibvar==FALSE): define vector containing index years to be used for calibration (remainder will be used for verification, if applicable). use (1:calib.length) for using all years, i.e. no verification.
calib.years.fix <- (1:calib.length) # default: (1:calib.length)

# End Calibration and Verification period options ### ### ###



# Second part defines the ensemble based on regression residuals

# residual noise options ### ### ### 

# add residual-based AR noise to the recon to quantify uncertainties? 
# this creates an ensemble of reconstructions based on residual-based AR noise modeling
add.arnoise <- FALSE # default: TRUE

# two versions are possible
# "simple" : Simply adding residuals based AR(1) noise to the recon
# "gene" : Do the re-fitting after the noise calculation and addition as in Wahl & smerdon (2012). This will results in narrower uncertainty bands.
arnoise.version <- "simple" # default: "simple"

# For field recon only: If "simple": resample from residual covariance matrix to generate spatially consistent fields using rmvn (see Jianghaos papers)? Is slow.
# FALSE it adds to noise to each grid cell separately, i.e. no spatial consistency in the noise
add.arnoise.spat <- TRUE # default: TRUE
# if TRUE: should PCs of reisidual fields be taken before using rmvn. This makes it much faster and results are usually very similar
do.pc.arnoise.spat <- TRUE # default: TRUE


# how many noise ensemble members should be created? 
# Attention: If MCiterations is > 1: this will crate an ensemble of ONLY noise perturbations. nens will be set to MCiterations
# and all verification options and ensemble mean outputs will be turned off. dopar needs to be FALSE. run this with "optimal" reconstruction parameters for reasonable results because only one recn will be performed and the noise will be created based on the residuals of this single recon.
# for combined calibration/parameter uncertainties use MCiterations==1
# MCiterations==1 means that for each parameter ensemble member, noise will be added
MCiterations <- 1 # default: 1

#do the noise calculation parallel? (set T only if MCiterations>1 or dopar==F otherwise it probably messes up the parallel stuff)
dopar.addnoise <- FALSE # default: FALSE

# End residual noise options ### ### ### 




### Verification And Calibration Residual Output

# Choose which diagnostics, verification measures, and calibration residual outputs are written.

# This all applies if do.verif and/or do.verif.early are TRUE

# Calculate REs?
do.re <- FALSE # default: TRUE
# REs of spatial mean / index recon?
do.re.index <- FALSE # default: TRUE

# CEs?
do.ce <- FALSE # default: TRUE
# CEs of spatial mean / index recon?
do.ce.index <- FALSE # default: TRUE

# r-squared?
do.r2 <- FALSE # default: TRUE
# R2 of spatial mean / index recon?
do.r2.index <- FALSE # default: TRUE

# percentage of cells with positive REs? (only if do.re==TRUE and do.field==T)
do.proz.pos <- FALSE # default: TRUE

# calculate the RMSE of the ensemble mean for each year?
do.rmse <- FALSE # default: TRUE

# calculate the RMSE for each nest and member? (will only be done for the spatial mean)
do.rmse.all <- FALSE # default: TRUE

# calculate the RE/CE/R2 of the ensemble mean of all calib/verif years? (will only be done for the spatial mean)
# probably more robust than the ensemble mean of the metrics.
# it selects the years used for calibration for each member, generates the ensemble mean than calculate the verification statistics over the entire calibration / verification interval
# same for validation years. Details in Neukom et al. 2014 (SM)
do.res.calib.verif.mean <- FALSE # default: TRUE

# calculate ensemble score diagnostics in the early verification period? Based on Werner & Tingley 2015(CP)
# Note that nens needs to be >1 for this to work.
# It works also if do.verif.early==F, using the period early.start:early.end.
do.ens.scores <- FALSE
dopar.ens.scores <- TRUE # run the ensemble scores in parallel? 


# calculate and write the residuals and their AR1 coefficients?
# sd of residuals over all members will be written
do.residuals <- FALSE # default: FALSE
# median ar1 coefficient of ensemble members will be written. works only if do.residuals==T
do.residual.ar1 <- FALSE # default: FALSE

# If FALSE, write the residuals for each member and nest? (will only be done for the spatial mean)
write.residuals.all <- FALSE # default: FALSE



### General Advanced Settings

# Configure workflow options that apply across methods but are not usually changed for daily runs.


# Define the nests with a pre-defined block length reduce processing time
# recon period is divided in blocks of nest.blocks.l length. only proxies with no NA's within a block will be considered for the recon of these years.
# can reduce the number of nests drastically, but for some records a substantial amount of information may get truncated.
nest.blocks <- FALSE # default: FALSE

#define time period for nest blocks [years]
nest.blocks.l <- 200 # default: 200


# fill in missing data in proxy matrix during calibration period?
# In case there are missing values and this fill.in.proxies == F, there will be an error.
fill.in.proxies <- TRUE

# If yes, which method is to be used to infill the proxies?
# chose between:
#     - 'dineof': DINEOF (default, Taylor et al., 2013),
#     - 'cps' : CPS (preserves most variability),
#     - 'gap.eof': EOF based infilling from Scherrer & Appenzeller 2006. Can produce large outliers
#     - 'pc.recon': Similar to gap.eof but uses regression-based infilling not jsut SVD transformations
#     - 'lms.recon': Same as pc.recon but uses robust (lms) regression
infill.method <- "cps"


# Filtered data 
#Do some statistics for the filtered data (loess filter applied)?
do.filter <- FALSE # default: FALSE
# window length of the filter for the decadal time series
filterlength <- 31 # default: 31




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ###
### ### ### ###
### ### ### ### ### ###  Part III: Method-Specific Settings
### ### ### ###
### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ######


### PCR Settings

# Configure PCR-specific truncation, weighting, and variance-adjustment options.

# Adjust variance of nests? 
#adjust variance of reconstruction of each nest to the target? 
do.var.adj <- TRUE # default: TRUE


# PC Truncation
# sample number of PCs in ensembles members?
sample.pcs <- FALSE # default: TRUE

# should an optimization method for PC truncation be used?
do.pc.opt <- TRUE # default: FALSE
# 3 possible options:
# "north" : North rule of thumb (fastest)
# "rulen" : rule N (Overland and Preisendorfer 1982). Slow as PC confidence intervals are computed every time from random data (by default 100 times, can be changed in rep variable in f.pc.opt function)
# "empirical" : similar to north except that confidence intervalls are sampled based on the real data by bootstrapping (removing 10% of calibration years). Is also slow.
pc.opt <- "empirical" # default: "rulen"

# if yes, how many PCs around the "ideal" number from above should be used for sampling?
# (i.e. if this number is 1 and north yields 4 pcs then it is sampled between 3 and 5 proxies)
# use 0 for no sampling
# only implemented for proxy PCs currently (not instrumental PCs).
pc.sample <- 2 # default: 2

# pre-defined PC selection interval based on explained variance
# If sample.pcs==T and do.pc.opt==F provide PC truncation frame for proxy PCs (lower and upper bound in terms of explained variance)
proz.pc.min.s <- 0.40 # default: 0.4
proz.pc.max.s <- 0.90 # default: 0.9
# Provide PC truncation frame for instrumental PCs (lower and upper bound in terms of explained variance)
proz.pc.min.q <- 0.60 # default: 0.6
proz.pc.max.q <- 0.99 # default: 0.99

# If sample.pcs==F and do.pc.opt==F: specify the % of variance that should be covered by the pc's
# for proxies
proz.pc.s <- 0.7 # default: 0.7
# for instrumentals
proz.pc.q <- 0.95 # default: 0.95

# Or altenatively use a fixed number of PCs to be used
# for TRUE, sample.pcs must be FALSE
npc.fix <- FALSE # default: FALSE
#for proxies
npc.s <- 10 # default: 10
#for instruments
npc.q <- 10 # default: 10

# power of latitude weighting for instrumental field in the PC routine. use 1 to weigh by cosine, 0.5 to weigh by the square root of the cosine (see e.g. Wahl&Ammann 2007)
latweight.cos.power <- 0.5 # default: 0.5

# Sample Proxy weights in PC routine
# vary the weights for the proxies randomly in the ensemble members?
sample.weights <- FALSE # default: TRUE
## margins for the scaling of the predictors i.e. the maximum and minimum variance allowed in the scaled predictor matrix
minsc <- 0.67 # default: 0.67 (i.e. 1/1.5)
maxsc <- 1.5  # default: 1.5



### CPS Settings

# Configure CPS-specific weighting, distance, and search-radius options.

# This is all only relevant if do.cps == TRUE

# should the proxies be weighted (by correlation and potentially additional weight factors)
do.cps.cor.weight <- TRUE # default: TRUE

# sample weighting factor as in cook et al. 2010?
# if FALSE proxies will be directly weighted (multiplied) with their correlation with the target over the calib period
sample.cps.weight <- FALSE # default: TRUE

# If true, how should the weighting be done:
# "c": Do the weighting using r^p as in cook et al 2010
# "u" : uniform around the best estimate of the correlation?
# "c" gives more weight to higher correlated series but results will not be equally distributed around the best estimate based on r
# "u": multiplies the weightfactor by a factor drawn within [minsc.cps maxsc.cps], current default [1/1.5 1.5]
cps.weighting <- "u" # default: "u"
# define cps weight factors for "c" (cook et al. 2010)
cps.weightfactors <- c(0,0.5,0.67,1,1.5,2) # default: c(0,0.5,0.67,1,1.5,2)
# define cps weight factors for "u"
minsc.cps <- 0.67 # default: 0.67
maxsc.cps <- 1.5 # # default: 1.5


# weigh proxies by distance to each other? (to account for spatially clustered proxies)
# this will be combined with the correlation based weighting.
# coordinates need to be provided in file below
# experimental! Rarely tested
cps.weigh.distances <- FALSE # default: F

# options for field (Point-by-point) CPS
# file with the coordinates of the proxies
# order must be the same as in the proxy file
# format: rows: name, latitude, longitude; columns: proxy records
# attention: convention of coordinates must be the same as in the target grid file! (can be changed for lons, see next option)
coord.file <- "C:/Users/neukom/Documents/Database/PAGES2k_DB/Data/CFR_comparison/coords_true-annual_calib-selection.csv"

# change lon convention of proxies? (to match the convention in the target grid)
# "n": no change
# "a": change from -180 - +180 to 0-360
# "b": change from 0-360 to -180 - +180
changelon.proxies <- "a" # default: "a"

# should only proxies within a certain search radius be included?
do.radius <- FALSE # default: FALSE

# should this radius be sampled for each ensemble member?
sample.radius <- FALSE # default: FALSE
# If TRUE define range for the search radius (km)
min.radius <- 500 # default: 500
max.radius <- 10000 # default: 1000
# If FALSE define search radius (km)
search.radius <- 2000 # default: 2000
# minimum number of proxies required to do recon
minproxies.radius <- 20 # default: 20

# what if less than minproxies are available within the radius?
# TRUE: allocate NA to this cell
# FALSE: use the [minproxies.radius] closest proxies for recon
# note that if TRUE, early nests with nproxies<minproxies,radius will be globally NA
# this only becomes active if do.radius == TRUE
few.na <- FALSE # default: FALSE


# parallel computing of grid-wise CPs (should probably only be set to TRUE if dopar == F)
dopar.cps.field <- FALSE # default: F



### CCA Settings

# Configure CCA-specific truncation and parameter-file options.

# This is all only relevant if do.cca == TRUE

# Note that CCA works only if do.field==T and do.index==F (i.e. only for field recons)
# Details about the parameters see Smerdon 2010, or Wang 2014

cca_options <- list() # leave unchanged

#cca_options$method <- 'smerdon10'   # NE08 not implemented currently, do not change
cca_options$dp_max <- 3  # maximum dp # default: 10
cca_options$dt_max <- 3  # maximum dt # default: 10
cca_options$dcca_max <- 12  # maximum dcca (will be min(dt,dp,dcca_max)) # default: 10

cca_options$K <- 2        # user-defined, for k-fold Cross validation # default: 2, warning: calculation time for the truncation parameters multiplies with this variable

# Load/save truncation parameters?
# calculating the parameters takes most of the time, so this can be done separately.
# if loadparas==T the parameters will be loaded from a .mat file that has been saved in an earlier run using saveparas==T.
# parameter files are called "params_n.mat" where "n" is the nest number (called "patterns" in cca code).
# saveparas should be used with nens==1 so each nest writes a single parameter file; ensemble-style parameter resampling can then be done later when re-loading.
# if loadparas==F, parameters will be computed with dopar so define maxcores above. 
cca_options$loadparas <- FALSE 
# save one parameter file per model for future reconstructions?
cca_options$saveparas <- FALSE

#sample truncation parameters over the ensemble
cca_options$sample_params <- FALSE # default: T
# define sample range (parameters will be sample +- the sample variable from the best estimate)
cca_options$dcca_sample <- 2 # default: 2
cca_options$dt_sample <- 2  # default: 2
cca_options$dp_sample <- 2 # default: 2

#location of parameter .mat files in case of loadparas==T
cca_options$parafileroot <- "TEST/params_"


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ###
### ### ### ###
### ### ### ### ### ###  Run the reconstruction workflow
### ### ### ###
### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ######

recon.files.path.use <- if(nzchar(recon.files.path) && !grepl("[/\\\\]$", recon.files.path)) paste0(recon.files.path, "/") else recon.files.path
source(paste0(recon.files.path.use, "Recon_workflow_clean.r"))







