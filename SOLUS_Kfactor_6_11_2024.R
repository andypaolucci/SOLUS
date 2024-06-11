library(terra)
library(rpostgis)
library(aws.s3)

# Connect to projects db as channel
channel <-  dbConnect(
  RPostgreSQL::PostgreSQL(),
  #RPostgres::Postgres(),
  user = Sys.getenv("sagauserid"), 
  password = Sys.getenv("sagapwd"),
  host = Sys.getenv("prodDB"),
  port = Sys.getenv("prodport"),  
  dbname = Sys.getenv("projectsDB"))

# connect to S3 bucket
bucket <- get_bucket("s3-fpac-nrcs-dshub-prod")

# list files in bucket
bucket_files <- get_bucket_df(bucket = bucket, max = 20000) 
bucket_files[with(bucket_files, grepl("solus/", Key)), ]

# Kfactor calculation page 618 Handbook
  # Kf <- (2.1 * M^1.14 * 10^-4 * (12 - a) + 3.25 * (b - 2) + 2.5 * (c - 3)) / 100 
 
# parameter: Depths
depths <- c(0,5,15,30,60,100,150)

# Directory file structure 
# ~/SOLUS/Inputs/solus/ --  raw SOLUS rasters saved from AWS S3 bucket
# ~/SOLUS/Inputs/structure/ -- SSURGO structure tiles created in this script 
# ~/SSURGO/30m_gssurgo/conus_30m_tiles -- 30m gSSURGO saved from postgres Db
# ~/SOLUS/Inputs/solus_ksat/ -- SOLUS ksat rasters created in previous script and coded minimum ksat layer created in this script
#~/SOLUS/Outputs/kfactor/kf_raw/ -- SOLUS Kf factor values before reclassifying and setting frags >90% ==0.02
#~/SOLUS/Outputs/kfactor/kf_classifed/ -- SOLUS Kf factor values reclassified and soils with frags >90% set to 0.02
#~/SOLUS/Outputs/kfactor/kw_classifed/ -- SOLUS Kw factor values calculated using fragment volume, reclassified, and soils with frags >90% set to 0.02

#### a: SOM ####
# divide by 1000 to get raw-unadjusted SOC value
# convert from SOC to SOM (SOC * 1.72)
# assign >4% SOM the value 4 based on NASIS calculation 
for (i in depths){
  soc_n <- rast(paste0("~/SOLUS/Inputs/solus/soc_r_1000x_", i, "_cm_2D_QRFadj_bt.tif"))
  soc_n_adj <- soc_n/1000
  som_n <- soc_n_adj * 1.72
  som_n_adj <-  ifel(som_n > 4, 4, som_n)
  
  writeRaster(som_n_adj, paste0("~/SOLUS/Inputs/solus_som/som_1x_", i, "_cm.tif"), overwrite=TRUE)      
}

#### b: structure code ####
# (1 =vfgr/sg, 2=fgr, 3=mcgr, 4=blocky, platy, or massive)
# assumed no data values were blocky 
all.structures <-read.csv('~/SOLUS/Inputs/structure/str_data.csv')

# ssurgo tiles 
files <- list.files(path="~/SSURGO/30m_gssurgo/conus_30m_tiles", pattern="*.tif", recursive=FALSE)

# filter table and join to gSSURGO
for (i in depths){
  ratings <- all.structures[ , c('mukey', paste0('strcode_', i, 'cm'))]
  row.names(ratings)<-NULL
  
  # join table with spatial. This takes awhile but you can split up the tiles and run 5 sessions at once
  for (x in files) { 
    tile_file  <- paste0("~/SSURGO/30m_gssurgo/conus_30m_tiles/", x)
    # load raster
    ssurgo_tile <- rast(tile_file)
    # set names to mukey
    names(ssurgo_tile) <- "mukey"
    
    # Subset ratings table using unique of mukeys in tile 
    mukey.df <- data.frame(unique(ssurgo_tile$mukey))
    ratings_sub <- ratings[ratings$mukey %in% mukey.df$mukey,]
    
    # add table to tile 
    levels(ssurgo_tile) <-ratings_sub
    
    # catalyze each raster to set attribute to raster value 
    strcode <- catalyze(ssurgo_tile)
    
    strcode_ifel<- ifel(strcode >4 , 4, strcode)
    
    writeRaster(strcode_ifel, paste0("~/SOLUS/Inputs/structure/str_", i, "cm/str_", i,"cm_", x))
  }
}

# resample structure to SOLUS resolution
resample.str <- function(strdepth){
  # virtual raster from str tiles
  str <- vrt(list.files(paste0("~/SOLUS/Inputs/structure/str_", strdepth, "cm"), recursive=TRUE, full.names=TRUE))
  # SOLUS raster 
  SOLUS <- rast("~/SOLUS/Inputs/solus/fragvol_r_1x_60_cm_2D_QRFadj_bt.tif")
  # empty raster with dimensions of SOLUS 
  empty_raster <- rast(xmin=xmin(SOLUS), xmax=xmax(SOLUS), ymin=ymin(SOLUS), ymax=ymax(SOLUS), ncols=ncol(SOLUS), nrows=nrow(SOLUS), res=res(SOLUS))
  crs(empty_raster)<- crs(SOLUS)
  # resample SSURGO-> SOLUS resolution 
  str.resamp <- resample(str, empty_raster, threads=TRUE, filename =paste0("~/SOLUS/Inputs/structure/str_SOLUS_res/str_", strdepth, "cm_100m_SSURGO_12_27_2023.tif"), overwrite=TRUE)
  return(str.resamp)
}
depths <-c(15)
lapply(depths, resample.str)

#### c minimum ksat codce #### 
# ksat code (1-6): min ksat within permeability control section 
#The permeability control section is the zone from the top of the mineral soil
#layer being evaluated to a depth of 50 cm below the top of that soil layer but should not
#exceed a profile depth of 200 cm

# c permeability control sections used for each depth of SOLUS data 
#  0 cm : 0-50. use existing min ksat 0-50. calculated from 0,5,15,30,45 cm observations
#  5 cm : 5-55. calculate min ksat between 5, 15, 30, and 45 
# 15 cm : 15-65. calculate min ksat between 15, 30, 45, and 60
# 30 cm : 30-80. calculate min ksat between 30, 45, and 60
# 60 cm : 60-110 calculate min ksat between 60 and 100 
# 100 cm : 100-150 calculate min ksat between 100 and 150
# 150 cm : 150-200. use 150 cm ksat as its the only observation within this depth range

# c codes based on minimum ksat 
# minksat >30 then 1
# minksat >15 then 2
# minksat >4.8 then 3
# minksat >1.2 then 4
# minksat >0.3 then 5
# minksat <=0.3 then 6 

# function for getting ksat class code: c 
minksat.class <- function(minksat, reslayer, resdepth, fname){
  minksat[minksat > 30] <-1
  minksat[minksat > 15 & minksat <= 30] <-2
  minksat[minksat > 4.8 & minksat <= 15] <-3
  minksat[minksat > 1.2 & minksat <= 4.8] <-4
  minksat[minksat > 0.3 & minksat <= 1.2] <-5
  minksat[minksat <= 0.3] <-6
  # classify bedrock layer to flay 
  resksat <- ifel(reslayer < resdepth, 0, 7)
  # get adjusted min ksat for when bedrock is within control section 
  minksat_adj <- min(minksat, resksat, filename = fname)
  return(minksat_adj)
}

# depth to restrictive feature layer. this includes 
res.dept <- rast("~/SOLUS/Inputs/solus/resdept_r_1x_all_cm_2D_QRFadj.tif")

# 0 - 50 cm adjusted ksat 
# min ksat 0-50 cm 
minksat_0_50 <-  rast("~/SOLUS/Inputs/solus_ksat/minksat_1x_0_50cm.tif")

minksat.class(minksat_0_50, res.dept, 50, "~/SOLUS/Inputs/solus_ksat/minksat_adj_restrc_0_50cm.tif")

# 5 - 55 cm adjusted ksat 
# min ksat between 5, 15, 30, 45, and 60 cm depths
ksat_5 <-  rast("~/SOLUS/Inputs/solus_ksat/ksat_1x_5cm_saxton.tif")
ksat_15 <-  rast("~/SOLUS/Inputs/solus_ksat/ksat_1x_15cm_saxton.tif")
ksat_30 <-  rast("~/SOLUS/Inputs/solus_ksat/ksat_1x_30cm_saxton.tif")
ksat_45 <-  rast("~/SOLUS/Inputs/solus_ksat/ksat_1x_45cm_saxton.tif")
ksat_60 <-  rast("~/SOLUS/Inputs/solus_ksat/ksat_1x_60cm_saxton.tif")

# minimum kat 5-55
minksat_5_15 <- min(ksat_5, ksat_15)
minksat_5_30 <-min(minksat_5_15, ksat_30)
minksat_5_45 <-min(minksat_5_30, ksat_45)
minksat_5_55 <-min(minksat_5_45, ksat_60)

minksat.class(minksat_5_55, res.dept, 55, "~/SOLUS/Inputs/solus_ksat/minksat_adj_restrc_5_55cm.tif")

# 15 - 65 cm adjusted ksat 
# calculate min ksat between 15, 30, 45 and 60 cm depths
minksat_15_30 <-min(ksat_15, ksat_30)
minksat_15_45 <-min(minksat_15_30, ksat_45)
minksat_15_65 <-min(minksat_15_45, ksat_60)

minksat.class(minksat_15_65, res.dept, 65, "~/SOLUS/Inputs/solus_ksat/minksat_adj_restrc_15_65cm.tif")

# 30 - 80 cm adjusted ksat 
# calculate min ksat between 30, 45 and 60 cm depths
minksat_30_45 <-min(ksat_30, ksat_45)
minksat_30_80 <-min(minksat_30_45, ksat_60)

minksat.class(minksat_30_80, res.dept, 80, "~/SOLUS/Inputs/solus_ksat/minksat_adj_restrc_30_80cm.tif")

# 60 - 110 cm adjusted ksat 
# 100 cm ksat
ksat_100 <- rast("~/SOLUS/Inputs/solus_ksat/ksat_1x_100cm_saxton.tif")
# min ksat between 60 and 100 cm
minksat_60_110 <-min(ksat_60, ksat_100)

minksat.class(minksat_60_110, res.dept, 110, "~/SOLUS/Inputs/solus_ksat/minksat_adj_restrc_60_110cm.tif")

# 100 - 150 cm  
# 150 cm ksat
ksat_150 <- rast("~/SOLUS/Inputs/solus_ksat/ksat_1x_150cm_saxton.tif")
# min ksat between 100 and 150 cm
minksat_100_150 <-min(ksat_100, ksat_150)

minksat.class(minksat_100_150, res.dept, 150, "~/SOLUS/Inputs/solus_ksat/minksat_adj_restrc_100_150cm.tif")

# 150 - 200 cm adjusted ksat. NOTE SOLUS dataset goes to 150 cm so only one observation 
minksat.class(ksat_150, res.dept, 200, "~/SOLUS/Inputs/solus_ksat/minksat_adj_restrc_150_200cm.tif")

#### M:(silt + vfs) * (100 - clay) ####
# M function 
kfactor_m <- function(silt, vfs, clay){
  # Calculate "M" input for K Factor calculation
  # m = (silt + vfs) * (100 - clay)
  m = (silt + vfs) * (100 - clay)
  return(m)
}
for (i in depths){
  silt_n <- rast(paste0("~/SOLUS/Inputs/solus/silttotal_r_1x_", i, "_cm_2D_QRFadj.tif"))
  vfs_n <- rast(paste0("~/SOLUS/Inputs/solus/sandvf_r_1x_", i, "_cm_2D_QRFadj.tif"))
  clay_n <- rast(paste0("~/SOLUS/Inputs/solus/claytotal_r_1x_", i, "_cm_2D_QRFadj.tif"))
  
  m_n <- kfactor_m(silt = silt_n, 
                   vfs = vfs_n, 
                   clay = clay_n)
  
  writeRaster(m_n, paste0("~/SOLUS/Inputs/solus_m/m_1x_", i, "cm.tif"), overwrite=TRUE)      
}

#### K factor calculations ####
# load rasters for each depth and make tiles 
for (i in depths){
  
# a percent organic matter
  # assign >4% SOM the value 4
 a <- rast(paste0("~/SOLUS/Inputs/solus_som/som_1x_", i, "_cm.tif")) 
  
# b structure code (1 =vfgr/sg, 2=fgr, 3=mcgr, 4=blocky, platy, or massive)
 b <- rast(paste0("~/SOLUS/Inputs/structure/str_SOLUS_res/str_", i, "cm_100m_SSURGO_12_27_2023.tif"))

# c ksat code (1-6): min ksat within profile
  # assign bedrock ksat = 0
  # minksat >30 then 1
  # minksat >15 then 2
  # minksat >4.8 then 3
  # minksat >1.2 then 4
  # minksat >0.3 then 5
  # minksat <=0.3 then 6 
 c <- rast(paste0("~/SOLUS/Inputs/solus_ksat/minksat_", i, "cm_kfactor.tif"))

# M  (silt + vfs) * (100 - clay)
 M <- rast(paste0("~/SOLUS/Inputs/solus_m/m_1x_", i, "cm.tif"))
 plot(M)
# make stack
k.inputs<- c(a, b, c, M)
names(k.inputs) <-c("a", "b", "c", "M")
#writeRaster(k.inputs, filename = paste0("./Inputs/k_inputs", i, "_1_3_2024.tif"), overwrite=FALSE)

# create grid for tiles 
x <- rast(ncols=10, nrows=10, extent = ext(k.inputs), crs = crs(k.inputs))
plot(x)
# make tiles 
makeTiles(k.inputs, x, filename=paste0("~/SOLUS/Inputs/k_inputs/k_", i, "cm_tiles/k_tiles_", i, "cm_.tif"))
}

# Kf raw value function 
calcK <- function(a, b, c, M, Kf){
  Kf <- (2.1 * M^1.14 * 10^-4 * (12 - a) + 3.25 * (b - 2) + 2.5 * (c - 3)) / 100 
  return(Kf)
}

# Run Kf raw value function on tiles 
for (i  in depths){
# input tiles 
inputtiles <- list.files(paste0("~/SOLUS/Inputs/k_inputs/k_", i, "cm_tiles/"), full.names = TRUE)

# make empty raster and fill with output for each tile stack
for(j in 1:length(inputtiles)){
  # raster tile stack 
  y <- inputtiles[j]
  k.inputs <- rast(y)
  # create band for output
  Kf <- rast(xmin=xmin(k.inputs), xmax=xmax(k.inputs), ymin=ymin(k.inputs), ymax=ymax(k.inputs), ncols=ncol(k.inputs), nrows=nrow(k.inputs), res=100)
  crs(Kf) <- crs(k.inputs)
  # add to stack and set name to Kf
  k.inputs <- c(k.inputs, Kf)
  names(k.inputs)[5]<-"Kf"
  # run function
  lapp(c(k.inputs$a, k.inputs$b, k.inputs$c, k.inputs$M, k.inputs$Kf), 
       fun=calcK,
       filename = paste0("~/SOLUS/Outputs/kfactor/kf_raw/kfactor_", i, "cm/kfactor_", i, "cm_1_12_2024.tif"),
       overwrite = FALSE,
       wopt = list(datatype = "FLT4S"))
  }
}

# Kf reclassified value function 
classifyKf <- function(Kf, frags, Kf_class){
  # reclassify
  Kf_class[Kf > 0.595] <- 0.64
  Kf_class[Kf <= 0.595 & Kf > 0.520] <- 0.55
  Kf_class[Kf <= 0.520 & Kf > 0.460] <- 0.49
  Kf_class[Kf <= 0.460 & Kf > 0.400] <- 0.43
  Kf_class[Kf <= 0.400 & Kf > 0.345] <- 0.37
  Kf_class[Kf <= 0.345 & Kf > 0.300] <- 0.32
  Kf_class[Kf <= 0.300 & Kf > 0.260] <- 0.28
  Kf_class[Kf <= 0.260 & Kf > 0.220] <- 0.24
  Kf_class[Kf <= 0.220 & Kf > 0.185] <- 0.20
  Kf_class[Kf <= 0.185 & Kf > 0.160] <- 0.17
  Kf_class[Kf <= 0.160 & Kf > 0.125] <- 0.15
  Kf_class[Kf <= 0.125 & Kf > 0.075] <- 0.10
  Kf_class[Kf <= 0.075 & Kf > 0.035] <- 0.05
  Kf_class[Kf <= 0.035] <- 0.02
  # fragment content greater than 90% get 0.02 (NASIS script uses in-liu of texture)
  Kf_class[frags > 90] <- 0.02
  return(Kf_class)
}

# Kw reclassified value function
classifyKw <- function(Kf, frags, Kw_class){
  #IF fragrock >= 15 AND NOT ISNULL(fragrock) THEN kf*(1.0336209-0.048228992*POW(LOGN(fragrock),2)) ELSE kf.
  Kw_class[frags >=15] <- Kf[frags >= 15] * (1.0336209 - 0.048228992 * log(frags[frags >= 15])^2)
  Kw_class[frags < 15] <- Kf[frags < 15]
  # >90% fragments == 0.02
  Kw_class[frags > 90] <- 0.02
  # reclassify
  Kw_class[Kw_class > 0.595] <- 0.64
  Kw_class[Kw_class <= 0.595 & Kw_class > 0.520] <- 0.55
  Kw_class[Kw_class <= 0.520 & Kw_class > 0.460] <- 0.49
  Kw_class[Kw_class <= 0.460 & Kw_class > 0.400] <- 0.43
  Kw_class[Kw_class <= 0.400 & Kw_class > 0.345] <- 0.37
  Kw_class[Kw_class <= 0.345 & Kw_class > 0.300] <- 0.32
  Kw_class[Kw_class <= 0.300 & Kw_class > 0.260] <- 0.28
  Kw_class[Kw_class <= 0.260 & Kw_class > 0.220] <- 0.24
  Kw_class[Kw_class <= 0.220 & Kw_class > 0.185] <- 0.20
  Kw_class[Kw_class <= 0.185 & Kw_class > 0.160] <- 0.17
  Kw_class[Kw_class <= 0.160 & Kw_class > 0.125] <- 0.15
  Kw_class[Kw_class <= 0.125 & Kw_class > 0.075] <- 0.10
  Kw_class[Kw_class <= 0.075 & Kw_class > 0.035] <- 0.05
  Kw_class[Kw_class <= 0.035] <- 0.02
  return(Kw_class)
}

# classify Kf and calculate Kw
for (i  in depths){
  # input tiles 
  inputtiles <- list.files(paste0("./Outputs/kfactor/kf_raw/kfactor_", i, "cm/"), full.names = TRUE)
  # fragment layer
  frags<- rast(paste0('~/SOLUS/Inputs/solus/fragvol_r_1x_', i, '_cm_2D_QRFadj_bt.tif'))
  # make empty raster and fill with output for each tile stack
  for(k in 1:length(inputtiles)){
    # raster tile stack 
    z <- inputtiles[k]
    Kf <- rast(z)
    # crop/mask fragments layer and set na values to zero. 
    frags.crop <- crop(frags, Kf)
    frags.mask <- mask(frags.crop, Kf)
    frags[is.na(frags)] <- 0 
    #  empty raster for Kf output
    Kf_class <- rast(xmin=xmin(Kf), xmax=xmax(Kf), ymin=ymin(Kf), ymax=ymax(Kf), ncols=ncol(Kf), nrows=nrow(Kf), res=100)
    crs(Kf_class) <- crs(Kf)
    # classifyKf 
    lapp(c(Kf, frags, Kf_class), 
         fun=classifyKf,
         filename = paste0("~/SOLUS/Outputs/kfactor/kf_classified/kf_class_", i, "cm/kf_class_", i, "cm_", k, ".tif"),
         overwrite = TRUE,
         wopt = list(datatype = "FLT4S"))
  }
}
Kf <- rast("~/SOLUS/Outputs/kfactor/kf_raw/kfactor_15cm/kfactor_15cm_1_12_2024.tif")
i<-15 
# calculate Kw
for (i  in depths){
  # input tiles 
  inputtiles <- list.files(paste0("./Outputs/kfactor/kf_raw/kfactor_", i, "cm/"), full.names = TRUE)
  # fragment layer
  allfrags<- rast(paste0('~/SOLUS/Inputs/solus/fragvol_r_1x_', i, '_cm_2D_QRFadj_bt.tif'))
  # make empty raster and fill with output for each tile stack
  for(k in 1:length(inputtiles)){
    # raster tile stack 
    z <- inputtiles[k]
    Kf <- rast(z)
    # crop/mask fragments layer and set na values to zero. 
    frags.crop <- crop(allfrags, Kf)
    frags <- mask(frags.crop, Kf)
    frags[is.na(frags)] <- 0 
    # empty raster for Kw output
    Kw_class <- rast(xmin=xmin(Kf), xmax=xmax(Kf), ymin=ymin(Kf), ymax=ymax(Kf), ncols=ncol(Kf), nrows=nrow(Kf), res=100)
    crs(Kw_class) <- crs(Kf)
    # classifyKw 
    lapp(c(Kf, frags, Kw_class), 
         fun=classifyKw,
         filename = paste0("~/SOLUS/test2.tif"),
         overwrite = FALSE,
         wopt = list(datatype = "FLT4S"))
  }
}

# plot
output.rast <- vrt(list.files("~/SOLUS/Outputs/kfactor/kf_classified/kf_class_15cm/", recursive=TRUE, full.names=TRUE))
plot(output.rast)












