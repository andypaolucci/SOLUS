library(terra)

#### Load inputs #### 
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
  # Convert (cm/hr) to (um/s); 1 cm = 10000 um and 1 hr = 3600 seconds 
  K8 = K8 * (1000/3600)
  return(K8)
}

# Run function on 0cm dataset
ksat_0cm_umsec <- Saxton_ks(sand = sand0cm_1x, 
                  clay = clay0cm_1x, 
                  bd = bd0cm_1x, 
                  percent = TRUE)

# write raster. using compression and converting to INT would save space 
writeRaster(ksat_0cm_umsec, "./Inputs/solus_ksat/ksat_1x_0cm_saxton.tif", overwrite=TRUE)

ksat_0cm_umsec <- rast("./Inputs/solus_ksat/ksat_1x_100cm_saxton.tif")

# plot 
plot(ksat_0cm_umsec, 
     breaks=c(0, 0.01, 0.1, 1, 10, 40, 99999), 
     col=viridis::viridis(6), 
     mar=c(3.1, 3.1, 2.1, 7.1),
     main="Surface (0cm) Saturated Hydraulic Conductivity (um/s)",
     legend=TRUE,
     plg=list(title="ksat(um/s)"))


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
  
  writeRaster(ksat_umsec, paste0("./Inputs/solus_ksat/ksat_1x_", i, "cm_saxton.tif"), overwrite=TRUE)      
}





