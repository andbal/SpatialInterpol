### Could be needed
#library(gstat)
#library(raster)
#library(caret)

# variable names - MUST be coherent with names in "Masterfile"
# "Bor__B__im","Humus____","Kali__K_O_","Kalkbedarf","Karbonate_","Kupfer__Cu","Magnesium_","Mangan__Mn","Phosphat__","pH_Wert__i","Zink__Zn__"

### Needed for map creation (probably not called correctly by OrdKrig)
library("SpatialPosition")

### To install SpatialInterpol library
#    library("devtools")
#    install_github("JBrenn/SpatialInterpol") # ORIGINAL version
#    install_github("andbal/SpatialInterpol") # MOD version (very important fix applied)
library("SpatialInterpol")

### Set paths and variables
# working path (where output will be created)
wpath <- "/home/lv70864/abalotti/Krig/Simulation"
# Path to calibrated parameter
var_file <- read.csv("/home/lv70864/abalotti/Krig/Calibration/20160920_idw_calibration.csv", stringsAsFactors = F,header=TRUE )
krige <- FALSE	# FALSE for IDW, TRUE for Kriging
local <- TRUE	# FALSE=Global interp, TRUE=Local interp
res <- 100	# Output resolution (in meters) - try first 500 or 1000 for global

# Do the maps recursively (achieved with validation=False)
for (row in rownames(var_file))
{
  i <- as.integer(row)
  if (krige)
  {  
    if (local)
    {
    # Local Kriging
    OrdKrig(wpath = wpath, datafolder = "master", validation = FALSE, npix = res, local = local,
            variable=var_file$variable[i], model=var_file$model[i], cutoff=c("AdigeVenosta"=var_file$rad_fit[i]),
            psill=c("AdigeVenosta"=var_file$psill[i]), nugget=c("AdigeVenosta"=var_file$nugget[i]),
            nmax=c("AdigeVenosta"=var_file$nmax[i]), omax=c("AdigeVenosta"=var_file$omax[i]),
            nmin=c("AdigeVenosta"=var_file$nmin[i]) )
    } else {
    # Global Kriging
    OrdKrig(wpath = wpath, datafolder = "master", validation = FALSE, npix = res, local = local,
            variable=var_file$variable[i], nmax=c("AdigeVenosta"=var_file$nmax[i]),
            omax=c("AdigeVenosta"=var_file$omax[i]), nmin=c("AdigeVenosta"=var_file$nmin[i]) )
    }
  } else {
    if (local)
    {
    # Local IDW
    OrdKrig(wpath = wpath, datafolder = "master", validation = FALSE, npix = res, local = local,
            inverseDistWeigths = TRUE, variable=var_file$variable[i],idp=c("AdigeVenosta"=var_file$idp[i]),
            nmax=c("AdigeVenosta"=var_file$nmax[i]), omax=c("AdigeVenosta"=var_file$omax[i]),
            nmin=c("AdigeVenosta"=var_file$nmin[i]) )
    } else {
    # Global IDW
    OrdKrig(wpath = wpath, datafolder = "master", validation = FALSE, npix = res, local = local,
            inverseDistWeigths = TRUE, variable=var_file$variable[i], idp=c("AdigeVenosta"=var_file$idp[i])
            nmax=c("AdigeVenosta"=var_file$nmax[i]), omax=c("AdigeVenosta"=var_file$omax[i]),
            nmin=c("AdigeVenosta"=var_file$nmin[i]) )
    }
  }

}
