# Calculate ksat with SOLUS sand, clay, and bulk density and use that with SSURGO minimum water table, and SOLUS soil depth to calculate HSG
# references Part 630 National Engineering Handbook Table 7-1
# soil_ksat function from https://github.com/saraya209/soil_ksat/tree/master/Functions. See note about unit conversion

# set working directory 
setwd("~/SOLUS")

# load libraries
library(terra)
library(rpostgis)

#### Get SSURGO muaggatt min water table ####
# Connect to projects db as channel
channel <-  dbConnect(
  RPostgreSQL::PostgreSQL(),
  #RPostgres::Postgres(),
  user = Sys.getenv("sagauserid"), 
  password = Sys.getenv("sagapwd"),
  host = Sys.getenv("prodDB"),
  port = Sys.getenv("prodport"),  
  dbname = Sys.getenv("projectsDB"))

# get table with nutrient leaching and runoff ratings 
ratings <- dbGetQuery(channel, "SELECT mukey, wtdepannmin FROM saga.geo_water_quality_last_mukey")

# set na min wt depths to 202
ratings$wtdepannmin[is.na(ratings$wtdepannmin)] = 202

# classify min wt to align with HSG criteria
ratings$wtdepannmin[ratings$wtdepannmin < 60] = 1
ratings$wtdepannmin[ratings$wtdepannmin >= 60 & ratings$wtdepannmin <= 100] = 2
ratings$wtdepannmin[ratings$wtdepannmin > 100] = 3

# get list of tiles 
files <- list.files(path="~/DS_saga/SSURGO/30m_gssurgo/conus_30m_tiles", pattern="*.tif", recursive=FALSE)

# join table with spatial. This takes awhile but you can split up the tiles and run 5 sessions at once and
for (x in files) { 
  tile_file  <- paste0("~/DS_saga/SSURGO/30m_gssurgo/conus_30m_tiles/", x)
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
  wtdepannmin <- catalyze(ssurgo_tile)
  
  writeRaster(wtdepannmin, paste0("~/SOLUS/Inputs/wtdepannmin/wtmin_", x))
}

#### Calculate SOLUS ksat #### 
clay0cm_1x <- rast("./Inputs/solus/claytotal_r_1x_0_cm_2D_QRFadj.tif")
sand0cm_1x <- rast("./Inputs/solus/sandtotal_r_1x_0_cm_2D_QRFadj.tif")
bd0cm_1x <- rast("./Inputs/solus/dbovendry_r_1x_0_cm_2D_QRFadj.tif")

# Saxton function 
Saxton_ks <- function(sand, clay, bd, percent = TRUE){
  # Calculates saturated hydraulic conductivity 
  # using Saxton et al. (1986) model
  # 
  # Args:
  #   sand,clay,...: Ratio or percent of particles of size sand and clay
  #   bd: bulk density of soil (g/cm^3)
  #   percent: Whether the size fractions are in ratio (0-1)
  #             or in percentages (0-100).
  #
  # Returns: Saturated hydraulic conducivity [(um/s)]
  #
  # Convert ratio to percent
  if(percent == FALSE){
    sand = 100*sand
    clay = 100*clay
  }
  silt = 100-sand-clay  # Calculate silt percent
  # Calculate porosity from bulk density and particle density = 2.65
  phi <- 1 - (bd / 2.65)
  #
  ### Saxton et al. (1986) model (verified)
  # Ks (cm/hr)
  K8 = exp(12.012- 
             0.0755 * sand + 
             (-3.895 + 
                0.03671 * sand - 
                0.1103 * clay + 
                0.00087546 * (clay^2)) / 
             phi)
  # Convert (cm/hr) to (um/s); 1 cm = 10000um / 1 hr= 3600 sec. 10000/3600=2.777778
  K8 = K8 * 2.777778
  return(K8)
}

# Run function on 0cm dataset
ksat_0cm_umsec <- Saxton_ks(sand = sand0cm_1x, 
                  clay = clay0cm_1x, 
                  bd = bd0cm_1x, 
                  percent = TRUE)

writeRaster(ksat_0cm_umsec, "./Inputs/solus_ksat/ksat_1x_0cm_saxton.tif")

# plot 
#plot(ksat_0cm_umsec, 
#     breaks=c(0, 0.01, 0.1, 1, 10, 40, 99999), 
#     col=viridis::viridis(6), 
#     mar=c(3.1, 3.1, 2.1, 7.1),
#     main="Surface (0cm) Saturated Hydraulic Conductivity (um/s)",
#     legend=TRUE,
#     plg=list(title="ksat(um/s)"))


# apply functions to all other depths 
depths <- c(5,15,30,60,100,150)

for (i in depths){
  clay_1x <- rast(paste0("./Inputs/solus/claytotal_r_1x_", i, "_cm_2D_QRFadj.tif"))
  sand_1x <- rast(paste0("./Inputs/solus/sandtotal_r_1x_", i, "_cm_2D_QRFadj.tif"))
  bd_1x <- rast(paste0("./Inputs/solus/dbovendry_r_1x_", i, "_cm_2D_QRFadj.tif"))
  
  ksat_umsec <- Saxton_ks(sand = sand_1x, 
                              clay = clay_1x, 
                              bd = bd_1x, 
                              percent = TRUE)
  
  writeRaster(ksat_umsec, paste0("./Inputs/solus_ksat/ksat_1x_", i, "cm_saxton.tif"))      
}

#### Calculate Minimum Ksat ####
# load ksat rasters 
ksat_0cm <-rast("./Inputs/solus_ksat/ksat_1x_0cm_saxton.tif")
ksat_5cm <-rast("./Inputs/solus_ksat/ksat_1x_5cm_saxton.tif")
ksat_15cm <-rast("./Inputs/solus_ksat/ksat_1x_15cm_saxton.tif")
ksat_30cm <-rast("./Inputs/solus_ksat/ksat_1x_30cm_saxton.tif")
ksat_60cm <-rast("./Inputs/solus_ksat/ksat_1x_60cm_saxton.tif")
ksat_100cm <-rast("./Inputs/solus_ksat/ksat_1x_100cm_saxton.tif")

# average 30 + 60cm ksat to get 45 cm 
ksat_45cm <- ((ksat_30cm + ksat_60cm)/ 2) 

# check 
#plot(ksat_45cm, 
#     breaks=c(0, 0.01, 0.1, 1, 10, 40, 99999), 
#     col=viridis::viridis(6), 
#     mar=c(3.1, 3.1, 2.1, 7.1),
#     main="Surface (0cm) Saturated Hydraulic Conductivity (um/s)",
#     legend=TRUE,
#     plg=list(title="ksat(um/s)"))

writeRaster(ksat_45cm, "./Inputs/solus_ksat/ksat_1x_45cm_saxton.tif")

# get minimum ksat between 0-50 cm 
minksat0_50cm <- min(ksat_0cm, 
                     ksat_5cm, 
                     ksat_15cm, 
                     ksat_30cm,
                     ksat_45cm)

writeRaster(minksat0_50cm, "./Inputs/solus_ksat/minksat_1x_0_50cm.tif")

# get minimum ksat between 0-60 cm 
minksat0_60cm <- min(minksat0_50cm, 
                     ksat_60cm)

writeRaster(minksat0_60cm, "./Inputs/solus_ksat/minksat_1x_0_60cm.tif")

# get minimum ksat between 0-100 cm 
minksat0_100cm <- min(minksat0_60cm, 
                      ksat_100cm)

writeRaster(minksat0_100cm, "./Inputs/solus_ksat/minksat_1x_0_100cm.tif")

#### Reclassify inputs ####
# min ksat 0-50
ksat50 <- rast("./Inputs/solus_ksat/minksat_1x_0_50cm.tif")

ksat50[ksat50 <= 1] <-1
ksat50[ksat50 > 1 & ksat50 <= 10] <-2
ksat50[ksat50 > 10 & ksat50 <= 40] <-3
ksat50[ksat50 > 40] <-4

writeRaster(ksat50, "./Inputs/solus_ksat/ksat50_classified.tif")

# min ksat 0-60
ksat60 <- rast("./Inputs/solus_ksat/minksat_1x_0_60cm.tif")

ksat60[ksat60 <= 1] <-1
ksat60[ksat60 > 1 & ksat60 <= 10] <-2
ksat60[ksat60 > 10 & ksat60 <= 40] <-3
ksat60[ksat60 > 40] <-4

writeRaster(ksat60, "./Inputs/solus_ksat/ksat60_classified.tif")

# min ksat 0-100
ksat100 <- rast("./Inputs/solus_ksat/minksat_1x_0_100cm.tif")

ksat100[ksat100 <= 0.4] <-1
ksat100[ksat100 > 0.4 & ksat100 <= 4] <-2
ksat100[ksat100 > 4 & ksat100 <= 10] <-3
ksat100[ksat100 > 10] <-4

writeRaster(ksat100, "./Inputs/solus_ksat/ksat100_classified.tif")

# classify bedrock depth
# load SOLUS anylithic raster 
bedrock <- rast("./Inputs/solus/anylithicdpt_1x_all_cm_2D_QRF")
# create matrix of RKLS breaks and class values 
bedrockclassvalues <- c(0, 50, 1, 50, 100, 2, 100, 201, 3)
bedrockclmat <- matrix(bedrockclassvalues, ncol=3, byrow=TRUE)
# classify RKLS using breaks and values 1-4 
bedrock_class <- classify(bedrock, bedrockclmat, include.lowest=TRUE, filename = "./Inputs/SOLUS_anylithic_classified.tif")

# The minwt input data is already classified in min water table step above

#### Stack inputs and Make Tiles ####

# SOLUS soil depth from anylithicdpt_1x_all_cm_2D_QRFadj.tif 
# integer
#<50 = 1 
#50-100 = 2 
# >100 = 3
depth <- rast("./Inputs/SOLUS_anylithic_classified.tif")

# SSURGO minimum water table depth from muaggatt (NA set to 202), resampled using nearest neighbor to SOLUS dimensions
# integer
# <60 =  1
#60-100 = 2
# >100  = 3
minwt <- rast("./Inputs/wtdepannmin/CONUS_wtdepannmin_100m_adj.tif")

# SOLUS minimum ksat 0-30 cm calculated using 0, 5, 15, and 30 cm sand, clay, bd SOLUS layers and saxton PTF
# Note: Actually should be 0-50 cm per HSG criteria 
# integer
# <1 = 1 
# >= 1 < 10 = 2
# >=10 <=40 = 3
# > 40 = 4
minksat50cm <- rast("./Inputs/solus_ksat/ksat50_classified.tif")


# SOLUS minimum ksat 0-60 cm calculated using 0, 5, 15, 30, and 60 cm sand, clay, bd SOLUS layers and saxton PTF
# integer
# <1 = 1 
# >= 1 < 10 = 2
# >=10 <=40 = 3
# > 40 = 4
minksat60cm <- rast("./Inputs/solus_ksat/ksat60_classified.tif")


# SOLUS minimum ksat 0-100 cm calculated using 0, 5, 15, 30, 60, and 100 cm sand, clay, bd SOLUS layers and saxton PTF
# integer
# <0.4 = 1 
# >= 0.4 < 4 = 2
# >=4 <=10 = 3
# > 10 = 4
minksat100cm <- rast("./Inputs/solus_ksat/ksat100_classified.tif")

# make stack
hsg.inputs<- c(depth, minwt, minksat50cm, minksat60cm, minksat100cm)
names(hsg.inputs) <-c("depth", "minwt", "minksat50cm", "minksat60cm", "minksat100cm")
writeRaster(hsg.inputs, filename = "./Inputs/hsg.inputs.tif", overwrite=TRUE)

# create grid for tiles 
x <- rast(ncols=10, nrows=10, extent = ext(hsg.inputs), crs = crs(hsg.inputs))

# make tiles 
makeTiles(hsg.inputs, x, filename="./Inputs/tiles/inputtile_.tif")

#### Calculate HSG ####

classifyHSG<- function(depth, minwt, minksat50cm, minksat60cm, minksat100cm, hsg){
  
  # < 50 cm depth
  hsg[depth == 1] <- 4 # D
  
  # 50-100 cm depth, <60 cm minwt, ksat depth range is 0-60 cm
  hsg[depth == 2 & minwt == 1 & minksat60cm == 4] <- 5 # A/D
  hsg[depth == 2 & minwt == 1 & minksat60cm == 3] <- 6 # B/D
  hsg[depth == 2 & minwt == 1 & minksat60cm == 2] <- 7 # C/D
  hsg[depth == 2 & minwt == 1 & minksat60cm == 1] <- 4 # D
  
  # 50-100 cm depth, 60-100 cm minwt, ksat depth range is 0-50 cm 
  hsg[depth == 2 & minwt == 2 & minksat50cm == 4] <- 1 # A
  hsg[depth == 2 & minwt == 2 & minksat50cm == 3] <- 2 # B
  hsg[depth == 2 & minwt == 2 & minksat50cm == 2] <- 3 # C
  hsg[depth == 2 & minwt == 2 & minksat50cm == 1] <- 4 # D
  
  # 50-100 cm depth, >100 cm minwt, ksat depth range is 0-50 cm 
  hsg[depth == 2 & minwt == 3 & minksat50cm == 4] <- 1 # A
  hsg[depth == 2 & minwt == 3 & minksat50cm == 3] <- 2 # B
  hsg[depth == 2 & minwt == 3 & minksat50cm == 2] <- 3 # C
  hsg[depth == 2 & minwt == 3 & minksat50cm == 1] <- 4 # D
  
  # >100cm dm depth, <60 cm minwt, ksat depth range is 0-100 cm 
  hsg[depth == 3 & minwt == 1 & minksat100cm == 4] <- 5 # A/D
  hsg[depth == 3 & minwt == 1 & minksat100cm == 3] <- 6 # B/D
  hsg[depth == 3 & minwt == 1 & minksat100cm == 2] <- 7 # C/D
  hsg[depth == 3 & minwt == 1 & minksat100cm == 1] <- 4 # D
  
  # >100 cm depth, 60-100 cm minwt, ksat depth range is 0-50 cm 
  hsg[depth == 3 & minwt == 2 & minksat50cm == 4] <- 1 # A
  hsg[depth == 3 & minwt == 2 & minksat50cm == 3] <- 2 # B
  hsg[depth == 3 & minwt == 2 & minksat50cm == 2] <- 3 # C
  hsg[depth == 3 & minwt == 2 & minksat50cm == 1] <- 4 # D
  
  # > 100 cm depth, >100 cm minwt, ksat depth range is 0 to 100 cm
  hsg[depth == 3 & minwt == 3 & minksat100cm == 4] <- 1 # A
  hsg[depth == 3 & minwt == 3 & minksat100cm == 3] <- 2 # B
  hsg[depth == 3 & minwt == 3 & minksat100cm == 2] <- 3 # C
  hsg[depth == 3 & minwt == 3 & minksat100cm == 1] <- 4 # D
  
  return(hsg)
}

# list input tiles 
inputtiles <- list.files("./Inputs/tiles/", full.names = TRUE)

# run function on tiles. Splitting this up into 4-5 sessions speeds things up drastically 
for(i in 1:100){
  x <- inputtiles[i]
  hsg.inputs <- rast(x)
  
  # create band for output
  hsg <- rast(xmin=xmin(hsg.inputs), xmax=xmax(hsg.inputs), ymin=ymin(hsg.inputs), ymax=ymax(hsg.inputs), ncols=ncol(hsg.inputs), nrows=nrow(hsg.inputs), res=100)
  crs(hsg) <- crs(hsg.inputs)
  # add to stack and set name to hsg
  hsg.inputs <- c(hsg.inputs, hsg)
  names(hsg.inputs)[6]<-"hsg"
  
  lapp(c(hsg.inputs$depth, hsg.inputs$minwt, hsg.inputs$minksat50cm, hsg.inputs$minksat60cm, hsg.inputs$minksat100cm, hsg.inputs$hsg), 
       fun=classifyHSG,
       filename = paste0("./Outputs/tiles/hsgtile_", i, ".tif"),
       overwrite = TRUE,
       wopt = list(datatype = "INT1U"))
}

# hsg .vrt 
hsg <- vrt(list.files("~/SOLUS/Outputs/tiles/", recursive=TRUE, full.names=TRUE), filename="SOLUS_hsg_9_24_2023.vrt")

# make one big stack
hsg.outputs <- c(hsg.inputs, hsg)
names(hsg.outputs)[6]<-"hsg"


#### Plot outputs ####
# HSG Output
plot(hsg, type="classes", 
     levels=c("A","B","C","D","A/D","B/D","C/D"), 
     col=c("#ffffbf","#fee08b","#fc8d59","#d53e4f","#4d9221","#998ec3","#3288bd"),
     main="Hydrologic Soil Group",
     mar=c(5,1,1,1))

#  100 cm min ksat
minksat100cm <- rast("./Inputs/solus_ksat/minksat_1x_0_100cm.tif")
minksat100class <- rast("./Inputs/solus_ksat/ksat100_classified.tif")


plot(minksat100cm, 
     breaks=c(0, 0.01, 0.1, 1, 10, 100, 99999), 
     col=viridis::viridis(6), 
     main="Hydrologic Soil Group",
     legend=TRUE,
     plg=list(title="ksat(um/s)"))

plot(minksat100class, 
     #breaks=c(0, 0.01, 0.1, 1, 10, 100, 99999), 
     col=viridis::viridis(6), 
     mar=c(3.1, 3.1, 2.1, 7.1),
     main="Saturated Hydraulic Conductivity (um/s)",
     legend=TRUE,
     plg=list(title="ksat(um/s)"))

# Crop input/output to inspect area 
depth <- rast("./Inputs/SOLUS_anylithic_classified.tif")
box <- ext(-1500000,-1000000, 2000000, 2500000)
depth.crop <- crop(depth, box)
plot(depth.crop, main="depth", col=viridis::viridis(3))




