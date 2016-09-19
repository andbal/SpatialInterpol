# Optimisation function for IDW

# choose variables to use from colnames of namefile (see below)
# "Bor__B__im","Humus____","Kali__K_O_","Kalkbedarf","Karbonate_","Kupfer__Cu","Magnesium_","Mangan__Mn","Phosphat__","pH_Wert__i","Zink__Zn__"

# wpath = "H:/Projekte/MONALISA/05_Arbeitsbereiche/BaA/05_Soil_Interpolation/02_additional_maps"
# datafile = "master/original_dataset/Masterfile_AdigeVenosta.txt"

# PAR       default     range
#idp        2           1-16            -> weighting power (idw interpolation)
#nmax       12          4-120           -> maximum number of points inside range (idw interpolation)
#omax       3           1-30            -> maximum number of points in each quadrant (idw interpolation)

# FIXED PAR
#radius     3000    -> range for idw interpolation
#nmin       1       -> minimum number of points inside range (idw interpolation)

### IMPORTANT NOTE
# "Model" need to be tested one at time!!!
# Every fold may use a different model, and make no sense
# do average on parameters related to different model.

library(gstat)
library(caret)
# library(hydroPSO)
library(sp)

OrdKrig_optim_idw <- function(par = c(idp = 2.0, nmax=12, omax=3),
                              par_intp = c(radius = 3000, nmin=1),
                              wpath = "H:/Projekte/MONALISA/05_Arbeitsbereiche/BaA/05_Soil_Interpolation/02_additional_maps",
                              datafile = "master/original_dataset/Masterfile_AdigeVenosta.txt",
                              variable = "Humus____", local = TRUE,
                              kfold=10)
{
    # # to comment for the function version
    # par = c(idp = 2.0, nmax=12, omax=3)
    # par_intp = c(radius = 3000, nmin=1)
    # wpath = "H:/Projekte/MONALISA/05_Arbeitsbereiche/BaA/05_Soil_Interpolation/02_additional_maps"
    # datafile = "master/original_dataset/Masterfile_AdigeVenosta.txt"
    # variable = "Humus____"
    # kfold=10
    # local=TRUE
    # ###
    
    # read table 
    worktab <- read.table(file = file.path(wpath, datafile), header = TRUE, sep = ",",dec = ".")
    worktab <- cbind(worktab$x_Coord, worktab$y_Coord, worktab[,variable])
  
    # matrix 2 data.frame
    worktab <- as.data.frame(worktab)
    # rename cols
    names(worktab) <- c("X","Y","VARIABLE")
    
    # zeros
    worktab[worktab$VARIABLE <= 0,"VARIABLE"] <- 0.001

    # get folds
    flds <- createFolds(y = worktab$VARIABLE, k = kfold, list = TRUE, returnTrain = FALSE)
    var_name = strsplit(x = variable,split = "_")[[1]][1]
    mydata_out <- data.frame()

    val_fold_df <- data.frame()
    val_out_df  <- data.frame()
    mydata_fold <- data.frame()
  
    for (i in 1:length(flds))
    {
        # training
        train_set <- worktab[-flds[[i]],]
        
        # validation set
        valid_set <- worktab[flds[[i]],]
        
        # IDW
        Xnew <- valid_set[,c("X", "Y")]
        Xnew <- SpatialPoints(Xnew)
        
        myloc <- data.frame("X" = train_set$X,"Y" = train_set$Y)
        myloc <- SpatialPoints(myloc)
        
        ord_krig <- idw(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew,
                        idp = par[1], nmax = par[2], nmin = par_intp[2], omax = par[3],
                        maxdist = par_intp[1] )
        
        names(ord_krig) <- c("predict")
        
        # frame with original data + predicted for ACTUAL fold
        val_fold_df <- data.frame(fold=i,valid_set, ord_krig$predict)
        
        # Statistics for ACTUAL fold
        rms_fold <- RMSE(pred = val_fold_df$ord_krig.predict, obs = val_fold_df$VARIABLE, na.rm = T)
        r2_fold <- lm(formula = ord_krig.predict~VARIABLE, data = val_fold_df)
        r2_fold <- c(summary(r2_fold)$r.squared,summary(r2_fold)$adj.r.squared)
        
        # frame with original data + frame with predicted for ALL folds
        val_out_df <- rbind(val_out_df, val_fold_df)
        mydata_fold <- rbind(mydata_fold, data.frame(i, rms_fold, r2_fold[1], r2_fold[2]) )
    }
    names(mydata_fold) <- c("fold","RMSE","R2","adj_R2")    
    
    # Create set of parameters/statistics as average of values got from 10 folds
    par_new <- as.numeric(apply(X = mydata_fold[-c(1)],MARGIN = 2,FUN = mean))
    par_new <- data.frame(t(par_new),
                          as.numeric(par["nmax"]), as.numeric(par_intp["nmin"]),
                          as.numeric(par["omax"]), stringsAsFactors = F )
    colnames(par_new) <- c("rmse","r2","adj_r2","nmax","nmin","omax")
    
    # statistics for ALL folds
    rms <- RMSE(pred = val_out_df$ord_krig.predict, obs = val_out_df$VARIABLE, na.rm = T)
    r2 <- lm(formula = ord_krig.predict~VARIABLE, data = val_out_df)
    r2 <- c(summary(r2)$r.squared,summary(r2)$adj.r.squared)
    
    # add statistics
    par_new <- data.frame(par_new,rms_alldata=rms,r2_alldata=r2[1],adj_r2_alldata=r2[2])
    mydata_out <- rbind(mydata_out, par_new)
    
  
    # return(RMSE(pred = val_out_df$ord_krig.predict, obs = val_out_df$VARIABLE, na.rm = T))
    return(r2[1])
  
}

# # keep care: trade of between search distance and number of NA estimations
# # the smaller the search radius, the better the estimation - but lot of NAs
# # How to solve?

# hydroPSO::hydroPSO(fn = OrdKrig_optim_idw, method="spso2011",
#                    lower = c(1,100,8,1), upper = c(16,1000,100,25),
#                    control=list(drty.out = "/home/jbre/R/OrdKrig/PSO_idw", npart=40, 
#                                 parallel="none", par.pkgs = c("gstat","caret","hydroGOF","sp")))
