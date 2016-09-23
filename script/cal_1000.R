#library(sp)
#library(gstat)
#library(raster)
#library(caret)
#library(hydroPSO)

# Load library
#if (FALSE){
#install.packages("devtools",repos="http://cran.wu.ac.at")
#library("devtools")
#install_git("https://gitlab.inf.unibz.it/Samuel.Senoner/hydroPSO.git")
#install_git("https://github.com/andbal/SpatialInterpol.git")
#}
library("SpatialInterpol")
library("parallel")

# Set paths and variables
cores <- detectCores()
parts <- cores*2
iters <- 500		# iteration of hydroPSO
wpath <- "/home/lv70864/abalotti/Krig/Simulation"
#dataf <- "master/no_outmask/Masterfile_AdigeVenosta_no_outmask.txt"
dataf <- "master/Masterfile_AdigeVenosta.txt"
krige <- TRUE	# TRUE=calibration for kriging, FALSE=calibration for IDW
variables <- list("Bor__B__im","Humus____","Kali__K_O_","Kalkbedarf","Karbonate_","Kupfer__Cu","Magnesium_","Mangan__Mn","Phosphat__","pH_Wert__i","Zink__Zn__")

######## CALIBRATION #########
# Calibration module HydroPSO require the parameters to change during the calibration process at the beginning
# of the function that will be tested. But this parameters are different between IDW and Krige, so two different
# function must be created
##############################

# Call the calibration module giving: function, method to use, min-max range of value that will be tried
# for each parameters
# SELECT ONLY ONE AT TIMES BECAUSE REQUIRE A LOT OF TIME

# Remember to set MinMax option to 'max' if output is R2 and 'min' if RMSE!!!

for (var in variables)
{
	if (krige)
	{
  	print("--- Calibration for KRIGING function ---")
	} else {
	print("--- Calibration for IDW function ---")
	}
	print(paste("--- Performing calibration for variable",var,"---"))

	if (krige)
	{
    # create output path
	cal_path <- file.path(wpath,paste("PSO_krige",cores,iters,"mpi_FALSE",sep="_"),var)
	dir.create(cal_path, recursive = T)
    # run hydroPSO
	hydroPSO::hydroPSO(fn = "OrdKrig_optim_krige", method="spso2011",
                       wpath = wpath, datafile = dataf, variable = var,
                	   lower = c(200,4,1), upper = c(2000,120,30),
	                   control=list(drty.out = cal_path, npart=parts,
                       maxit=iters, MinMax="max",parallel="parallel",
	                   par.pkgs = c("SpatialInterpol","gstat","caret","sp"),
                       write2disk=TRUE,REPORT=100,digits=15) )
	} else {
    # create output path
	cal_path <- file.path(wpath,paste("PSO_idw",cores,iters,"mpi_FALSE",sep="_"),var)
	dir.create(cal_path, recursive = T)
    # run hydroPSO
	hydroPSO::hydroPSO(fn = "OrdKrig_optim_idw", method="spso2011",
	                   wpath = wpath, datafile = dataf, variable = var,
	                   lower = c(1,4,1), upper = c(16,120,30),
	                   control=list(drty.out = cal_path, npart=parts,
                       maxit=iters, MinMax="max", parallel="parallel",
	                   par.pkgs = c("SpatialInterpol","gstat","caret","sp"),
	                   write2disk=TRUE,REPORT=100,digits=15) )
	}
}
