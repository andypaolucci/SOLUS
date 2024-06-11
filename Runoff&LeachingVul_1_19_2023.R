library(terra)
library(aws.s3)

# point to s3 bucket
bucket <- get_bucket("s3-fpac-nrcs-dshub-prod")

#### Load Raster Inputs ####
# taxcode
# 1 = histosols 
# 2 = non-histosols 
taxcode <- vrt(list.files("~/SOLUS/Inputs/taxcode/", recursive=TRUE, full.names=TRUE))
plot(taxcode, col=c('red', 'green'))

# wtcode
# 1 = apparent
# 2 = perched 
# 3 = not rated 
wtcode <- vrt(list.files("~/SOLUS/Inputs/wtcode/", recursive=TRUE, full.names=TRUE)) 
plot(wtcode, col=c('red','green','grey'))

# wt_24in_flag. CART uses 61 cm for runoff and 76 for leaching. we are just going to use 24" (61 cm) 
# 1 = water table <24" 
# 2 =  water table >24" 
wt24in <- vrt(list.files("~/SOLUS/Inputs/wt_24in_flag/", recursive=TRUE, full.names=TRUE)) 
plot(wt24in, col=c('red','green'))

# kwfactor 
# < 0.24 = 1
# 0.24 < 0.28 = 2
# 0.28 < 0.32 = 3
# >= 0.32 = 4
# original SSURGO kfactor layer. Thickest mineral layer between 0-15 cm. Doominant Component
#kwfact <- vrt(list.files("~/SOLUS/Inputs/kfactor/", recursive=TRUE, full.names=TRUE)) 
# 1/16/2024 New SOLUS Kfactor raster. Using 0 cm. TODO: consider aggregating 0, 5, and 15 cm depths 
  # load kfactor
  # kwfact <-vrt(list.files("~/SOLUS/Outputs/kfactor/kw_classified/kw_class_0cm/", recursive=TRUE, full.names=TRUE)) 
  # plot(kwfact)
  # # reclassify into codes 
  # kwfact[kwfact >= 0.32] = 4
  # kwfact[kwfact >= 0.28 & kwfact < 0.32] = 3
  # kwfact[kwfact >= 0.24 & kwfact < 0.28] = 2
  # kwfact[kwfact < 0.24] = 1
  # # write raster for later 
  # writeRaster(kwfact, "~/SOLUS/Outputs/kfactor/kw_classified/kwcode_0cm_1_16_2024.tif", overwrite=TRUE)
  kwfact <- rast("~/SOLUS/Outputs/kfactor/kw_classified/kwcode_0cm_1_16_2024.tif")
  # check 
  plot(kwfact, col=c("#ffffbf","#fee08b","#fc8d59","#d53e4f"))

# hsg 
hsg <- vrt(list.files("~/SOLUS/Outputs/hsgtiles/", recursive=TRUE, full.names=TRUE))
hsg<- rast("./Inputs/hsg_SOLUS_10_04_2023.vrt")
# 1 = A
# 2 = B
# 3 = C
# 4 = D
# 5 = A/D
# 6 = B/D
# 7 = C/D
plot(hsg, col=c("#ffffbf","#fee08b","#fc8d59","#d53e4f","#4d9221","#998ec3","#3288bd"))

# slope 
# not classified 
# load .tif
slope <- rast("./Inputs/solus/covariates_SLPNED6.tif")
plot(slope)

# # fragvol 
# # load frag vol tiffs
# frag0cm <- rast("./Inputs/solus/fragvol_r_1x_0_cm_2D_QRFadj_bt.tif")
# frag5cm <- rast("./Inputs/solus/fragvol_r_1x_5_cm_2D_QRFadj_bt.tif")
# frag15cm <- rast("./Inputs/solus/fragvol_r_1x_15_cm_2D_QRFadj_bt.tif")
# frag30cm <- rast("./Inputs/solus/fragvol_r_1x_30_cm_2D_QRFadj_bt.tif")
# frag60cm <- rast("./Inputs/solus/fragvol_r_1x_60_cm_2D_QRFadj_bt.tif")
# frag100cm <- rast("./Inputs/solus/fragvol_r_1x_100_cm_2D_QRFadj_bt.tif")
# frag150cm <- rast("./Inputs/solus/fragvol_r_1x_150_cm_2D_QRFadj_bt.tif")
# 
# # average across 7 depth observations. Not sure of any other way since these arent layers 
# average.frags <- (frag0cm + frag5cm + frag15cm + frag30cm + frag60cm + frag100cm + frag150cm)/7
# 
# writeRaster(average.frags, "./Inputs/solus/averagefragvol_solus_0_150cm.tif")
# 
# # classify for interpretation 
# # 1 = <= 10 
# # 2 = >10 <= 30
# # 3 = >30 
# 
# average.frags[average.frags <= 10] <- 1
# average.frags[average.frags > 10 & average.frags <= 30] <- 2
# average.frags[average.frags > 30] <- 3
# 
# plot(average.frags)
# 
# writeRaster(average.frags, "./Inputs/fragvol0_150_solus_10_4_2023.tif")
fragvol <- rast("./Inputs/fragvol0_150_solus_10_4_2023.tif")
plot(fragvol, col=c("green","yellow","red"))

#### Stack Inputs and Make Tiles ####

# make empty raster with dimensions of SOLUS 
empty_raster <- rast(xmin=xmin(hsg), xmax=xmax(hsg), ymin=ymin(hsg), ymax=ymax(hsg), ncols=ncol(hsg), nrows=nrow(hsg), res=res(hsg))
crs(empty_raster)<- crs(hsg)

# resample SSURGO rasters to SOLUS resolution 
taxcode.resamp <- resample(taxcode, empty_raster, threads=TRUE, filename ="./Inputs/taxcode_resamp.tif", overwrite=TRUE, method="near")
wtcode.resamp <- resample(wtcode, empty_raster, threads=TRUE, filename ="./Inputs/wtcode_resamp.tif", overwrite=TRUE, method="near")
wt24in.resamp <- resample(wt24in, empty_raster, threads=TRUE, filename ="./Inputs/wt24in_resamp.tif", overwrite=TRUE, method="near")
kwfact.resamp <- resample(kwfact, empty_raster, threads=TRUE, filename ="./Inputs/kwfact_resamp.tif", overwrite=TRUE, method="near")


# edit min and maxs of wt24in and 
taxcode.resamp <- rast("./Inputs/taxcode_resamp.tif")
wtcode.resamp <- rast("./Inputs/wtcode_resamp.tif")
wt24in.resamp <- rast("./Inputs/wt24in_resamp.tif")
#kwfact.resamp <- rast("./Inputs/kwfact_resamp.tif")

# reclassify any values above max to max 
wt24in.adj <- ifel(wt24in.resamp > 2, 2, wt24in.resamp, filename="./Inputs/wt24in_adj.tif", overwrite=TRUE)
#kwfact.adj <- ifel(kwfact.resamp > 4, 4, kwfact.resamp, filename="./Inputs/kwfact_adj.tif", overwrite=TRUE)

wt24in.adj<- rast("./Inputs/wt24in_adj.tif")
#ksfact.adj<- rast("./Inputs/kwfact_adj.tif")

# mask slope values that are outside SSURGO 
slope.mask <-terra::mask(slope, wtcode.resamp, filename="./Inputs/SOLUS_slope_mask.tif")
hsg.mask <-terra::mask(hsg, wtcode.resamp, filename="./Inputs/SOLUS_hsg_mask.tif")
fragvol.mask <-terra::mask(fragvol, wtcode.resamp, filename="./Inputs/SOLUS_fragvol_mask.tif")
kwfact.mask <-terra::mask(kwfact, wtcode.resamp, filename="./Inputs/SOLUS_kfactor0cm_mask.tif")


# load preprocessed inputs
kwfact.mask <- rast("./Inputs/SOLUS_kfactor0cm_mask.tif")
slope.mask <- rast("./Inputs/SOLUS_slope_mask.tif")
hsg.mask <- rast("./Inputs/SOLUS_hsg_mask.tif")
fragvol.mask <- rast("./Inputs/SOLUS_fragvol_mask.tif")
# make stack with inputs
vul.inputs <- c(taxcode.resamp, wtcode.resamp, wt24in.adj, kwfact.mask, hsg.mask, slope.mask, fragvol.mask)

names(vul.inputs) <-c("taxcode", "wtcode", "wt24in", "kwfact", "hsg", "slope", "fragvol")

writeRaster(vul.inputs, filename = "./Inputs/vulinputs_1_19_2024.tif", overwrite=TRUE)
vul.inputs <- rast("./Inputs/vulinputs.tif")
# create grid for tiles 
x <- rast(ncols=10, nrows=10, extent = ext(vul.inputs), crs = crs(vul.inputs))

# make tiles 
makeTiles(vul.inputs, x, filename="./Inputs/vultiles_1_19/vultile_.tif")

#### Runoff Vulnerability ####

# runoff criteria 
runoffvul <- function(wt24in, kwfact, slope, hsg, runoff){
  
  # percher or apparent water table within 61 cm or 24" is 3=high. Since its perched or aparrent I jsut used the 24" SHWT raster  
  runoff[wt24in == 1] <- 3 # High
  
  # Dual HSGs that are not drained are 3=high 
  runoff[wt24in > 1 & hsg >= 5] <- 3 # High
  
  # HSG: A 
  runoff[wt24in > 1 & hsg == 1] <- 0 # Low
  
  # HSG: B
  runoff[wt24in > 1 & hsg == 2 & slope < 4] <- 0 # Low
  runoff[wt24in > 1 & hsg == 2 & slope >=4 & slope <= 6 & kwfact <= 3] <- 1 # Mod
  runoff[wt24in > 1 & hsg == 2 & slope >=4 & slope <= 6 & kwfact == 4] <- 2 # ModHigh
  runoff[wt24in > 1 & hsg == 2 & slope > 6] <- 3 # High
  
  # HSG: C
  runoff[wt24in > 1 & hsg == 3 & slope < 2] <- 0 # Low
  runoff[wt24in > 1 & hsg == 3 & slope >=2 & slope <= 6 & kwfact <= 2] <- 1 # Mod
  runoff[wt24in > 1 & hsg == 3 & slope >=2 & slope <= 6 & kwfact >= 3] <- 2 # ModHigh
  runoff[wt24in > 1 & hsg == 3 & slope > 6] <- 3 # High
  
  # HSG: D
  runoff[wt24in > 1 & hsg == 4 & slope < 2 & kwfact <= 2] <- 0 # Low
  runoff[wt24in > 1 & hsg == 4 & slope < 2 & kwfact >= 3] <- 1 # Mod
  runoff[wt24in > 1 & hsg == 4 & slope >=2 & slope <= 4] <- 2 # ModHigh
  runoff[wt24in > 1 & hsg == 4 & slope > 4] <- 3 # High
  
  return(runoff)
}

# Calculate Runoff Vulnerability
# list input tiles 
inputtiles <- list.files("./Inputs/vultiles_1_19/", full.names = TRUE)

# run runoff function on tiles
# inputs (wtcode, wt24in, hsg, slope, kwfact, runoff)
for(i in 1:100){
  x <- inputtiles[i]
  vul.inputs <- rast(x)
  
  # create band for output
  runoff <- rast(xmin=xmin(vul.inputs), xmax=xmax(vul.inputs), ymin=ymin(vul.inputs), ymax=ymax(vul.inputs), ncols=ncol(vul.inputs), nrows=nrow(vul.inputs), res=100)
  crs(runoff) <- crs(vul.inputs)
  # add to stack and set name to hsg
  vul.inputs <- c(vul.inputs, runoff)
  names(vul.inputs)[7]<-"runoff"
  
  lapp(c(vul.inputs$wt24in, vul.inputs$kwfact, vul.inputs$slope, vul.inputs$hsg, vul.inputs$runoff), 
       fun=runoffvul,
       filename = paste0("./Outputs/runoffvul_1_19/runoffvultile_", i, ".tif"),
       overwrite = TRUE,
       wopt = list(datatype = "INT1U"))
}

#vrt
runoffvrt <- vrt(list.files("~/SOLUS/Outputs/runoffvul_1_19/", recursive=TRUE, full.names=TRUE), filename="~/SOLUS/Outputs/SOLUS_runoffvul_1_19_2024.vrt")

#plot
plot(runoffvrt, col=viridis::viridis(4))

#### leaching Vulnerability ####
leachingvul <- function(taxcode, wt24in, hsg, slope, kwfact, fragvol, leaching){
  
  # histosols is 3=high 
  leaching[taxcode == 1] <- 3 # High
  
  # percher or appart water table within 61 cm or 24" is 3=high. Note: CART Criteria uses 76cm for leaching but we are going to use 61 since at this scale 15cm isnt much 
  leaching[taxcode == 2 & wt24in == 1] <- 3 # High
  
  # Dual HSGs are 3=high 
  leaching[taxcode == 2 & wt24in > 1 & hsg >= 5] <- 3 # High
  
  # HSG: A 
  leaching[taxcode == 2 & wt24in > 1 & hsg == 1 & slope > 12] <- 2 # ModHigh
  leaching[taxcode == 2 & wt24in > 1 & hsg == 1 & slope < 12] <- 3 # High
  
  # HSG: B
  leaching[taxcode == 2 & wt24in > 1 & hsg == 2 & slope <= 12 & kwfact >= 2] <- 1 # Mod
  leaching[taxcode == 2 & wt24in > 1 & hsg == 2 & slope > 12] <- 1 # Mod
  leaching[taxcode == 2 & wt24in > 1 & hsg == 2 & slope >= 3 & slope <= 12 & kwfact == 1] <- 2 # ModHigh
  leaching[taxcode == 2 & wt24in > 1 & hsg == 2 & slope < 3 & kwfact == 1] <- 3 # High
  
  # HSG: C
  leaching[taxcode == 2 & wt24in > 1 & hsg == 3] <- 1 #Mod
  
  # HSG: D
  leaching[taxcode == 2 & wt24in > 1 & hsg == 4] <- 0 #Low
  
  # Coarse Fragment Corrections
  # > 30% fragments 
  leaching[leaching == 0 & fragvol == 3] <- 2 # 
  leaching[leaching == 1 & fragvol == 3] <- 3
  leaching[leaching == 2 & fragvol == 3] <- 3
  
  # 10-30% fragments 
  leaching[leaching == 0 & fragvol == 2] <- 1
  leaching[leaching == 1 & fragvol == 2] <- 2
  leaching[leaching == 2 & fragvol == 2] <- 3
  
  return(leaching)
}



# run leaching function on tiles
# inputs (taxcode, wtcode, wt24in, hsg, slope, kwfact, fragvol, leaching)
for(i in 1:100){
  x <- inputtiles[i]
  vul.inputs <- rast(x)
  
  # create band for output
  leaching <- rast(xmin=xmin(vul.inputs), xmax=xmax(vul.inputs), ymin=ymin(vul.inputs), ymax=ymax(vul.inputs), ncols=ncol(vul.inputs), nrows=nrow(vul.inputs), res=100)
  crs(leaching) <- crs(vul.inputs)
  # add to stack and set name to hsg
  vul.inputs <- c(vul.inputs, leaching)
  names(vul.inputs)[8]<-"leaching"
  
  lapp(c(vul.inputs$taxcode, vul.inputs$wt24in, vul.inputs$hsg, vul.inputs$slope, vul.inputs$kwfact, vul.inputs$fragvol, vul.inputs$leaching), 
       fun=leachingvul,
       filename = paste0("./Outputs/leachingvul_1_19/leachingvultile_", i, ".tif"),
       overwrite = TRUE,
       wopt = list(datatype = "INT1U"))
}

# .vrt 
leachingvrt <- vrt(list.files("~/SOLUS/Outputs/leachingvul_1_19/", recursive=TRUE, full.names=TRUE), filename="~/SOLUS/Outputs/SOLUS_leachingvul_1_19_2024.vrt", overwrite=TRUE)

#plot
plot(leachingvrt, col=viridis::viridis(4))

writeRaster(leachingvrt, filename="./Outputs/SOLUS_nslp_1_19_2024.tif")

# TODO: see why on NOTCOM areas area already masked out with leaching 


# load rast for masking 
#kwfact.adj <- rast("./Inputs/kwfact_adj.tif")

# mask
#leachingvrt.notcom<- terra::mask(leachingvrt, kwfact.adj)

# plot
#plot(leachingvrt.notcom, col=viridis::viridis(4))

