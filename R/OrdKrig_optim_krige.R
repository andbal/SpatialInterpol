# Optimisation function for ordinary kriging

# choose variables to use from colnames of namefile (see below)
# "Humus____" "pH_Wert__i" "Karbonate_" "Kalkbedarf" "Phosphat__" "Kali__K_O_" "Magnesium_" "Bor__B__im" "Mangan__Mn" "Kupfer__Cu" "Zink__Zn__"

# wpath = "H:/Projekte/MONALISA/05_Arbeitsbereiche/BaA/05_Soil_Interpolation/02_additional_maps"
# datafile = "master/original_dataset/Masterfile_AdigeVenosta.txt"

# PAR       default     range
#cutoff     300         200-1600/2000   -> range for experimental variogram
#nmax       12          4-120           -> maximum number of points inside range (krigin interpolation)
#omax       3           1-30            -> maximum number of points in each quadrant (krigin interpolation)
# NOT USED IN CURRENT VERSION
#psill      0.9         0-10            -> partial sill for experimental variogram
#nugget     0.1         0-1             -> nugget for experimental variogram

# FIXED PAR
#radius     3000    -> range for krigin interpolation, not for variogram and model
#nmin       1       -> minimum number of points inside range (krigin interpolation)

### IMPORTANT NOTE
# "Model" need to be tested one at time!!!
# Every fold may use a different model, and make no sense
# do average on parameters related to different model.

library(gstat)
library(caret)
# library(hydroPSO)
library(sp)

OrdKrig_optim_krige <- function(par = c(cutoff=300, nmax=12, omax=3),
                                par_kri = c(radius = 3000, nmin=1),
                                wpath = "H:/Projekte/MONALISA/05_Arbeitsbereiche/BaA/05_Soil_Interpolation/02_additional_maps",
                                datafile = "master/original_dataset/Masterfile_AdigeVenosta.txt",
                                variable = "Humus____", local = TRUE,
                                model = c("Exp"), kfold = 10 )
{
    # # to comment for the function version
    # par = c(cutoff=300, nmax=12, omax=3, psill=0.9, nugget=0.1)
    # par_kri = c(radius = 3000, nmin=1)
    # wpath = "H:/Projekte/MONALISA/05_Arbeitsbereiche/BaA/05_Soil_Interpolation/02_additional_maps"
    # datafile = "master/original_dataset/Masterfile_AdigeVenosta.txt"
    # variable = "Humus____"
    # kfold=10
    # local=TRUE
    # # Need to be tested one at time!!!
    # # Every fold may use a different model, and at the enf make no sense
    # # do average on parameters related to different model
    # model=c("Exp")#,"Sph")
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
    
    # # variable range for experimental variogram
    # range1 <- seq(100,500,50)
    # range2 <- seq(300,2400,300)
    # range_vario <- c(range2)#,range2)
    # 
    # for (k in range_vario)
    # {
    ### OR ###
    k=par[1]
    ###
    
        val_fold_df <- data.frame()
        val_out_df  <- data.frame()
        mydata_fold <- data.frame()
  
        for (i in 1:length(flds))
        {
        # training set
        train_set <- worktab[-flds[[i]],]
        # validation set
        valid_set <- worktab[flds[[i]],]
        
        # gstatVariogram - Calculate experimental variogram and fit a model on it
        if (local){
            # local variogram 
            my_var <- variogram(log(VARIABLE)~1, data=train_set, locations = ~X+Y, cutoff = k)#, width=par[1]/50)
            # m <- vgm(psill = par[5], model = model, range = par[1], nugget = par[6])
            m <- vgm(model)
            my_var_fit <- fit.variogram(my_var, m)
            # plot(my_var,model=my_var_fit)
            # print(my_var_fit)
        } else {
            # global variogram
            my_var <- variogram(log(VARIABLE)~1, data=train_set, locations = ~X+Y )
            # m <- vgm(psill = par[5], model = model, nugget = par[6])
            m <- vgm(model)
            my_var_fit <- fit.variogram(my_var, m)
        }
        
        # Ordinary Kriging
        Xnew <- valid_set[,c("X", "Y")]
        Xnew <- SpatialPoints(Xnew)
        myloc <- data.frame("X" = train_set$X,"Y" = train_set$Y)
        myloc <- SpatialPoints(myloc)
        
        if (local){    
            ord_krig <- krige(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew,
                              model = my_var_fit, nmax = par[2], nmin = par_kri[2],
                              omax = par[3], maxdist = par_kri[1] ) #maxdist = my_var_fit$range[2]
        } else {
            ord_krig <- krige(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew, model = my_var_fit )
        }
        names(ord_krig) <- c("predict", "variance")
        
        # frame with original data + predicted for ACTUAL fold
        val_fold_df <- data.frame(fold=i,valid_set, ord_krig$predict, ord_krig$variance)
        
        # Statistics for ACTUAL fold
        rms_fold <- RMSE(pred = val_fold_df$ord_krig.predict, obs = val_fold_df$VARIABLE, na.rm = T)
        r2_fold <- lm(formula = ord_krig.predict~VARIABLE, data = val_fold_df)
        r2_fold <- c(summary(r2_fold)$r.squared,summary(r2_fold)$adj.r.squared)
        
        # frame with original data + frame with predicted for ALL folds
        val_out_df <- rbind(val_out_df, val_fold_df)
        mydata_fold <- rbind(mydata_fold, data.frame(i, my_var_fit$model[2],
                                                     my_var_fit$psill[2],
                                                     my_var_fit$psill[1],
                                                     my_var_fit$range[2], k,
                                                     rms_fold, r2_fold[1], r2_fold[2]) )
        }# end loop on folds
        names(mydata_fold) <- c("fold","model","psill","nugget","rad_fit","range_exp","RMSE","R2","adj_R2")    
        
        # Create set of parameters/statistics as average of values got from 10 folds
        par_new <- as.numeric(apply(X = mydata_fold[-c(1:2)],MARGIN = 2,FUN = mean))
        par_new <- data.frame(as.character(my_var_fit$model[2]), t(par_new),
                              as.numeric(par["nmax"]), as.numeric(par_kri["nmin"]),
                              as.numeric(par["omax"]), stringsAsFactors = F )
        colnames(par_new) <- c("model","psill","nugget","rad_fit","rad_mod","rmse","r2","adj_r2","nmax","nmin","omax")
        
        # statistics for ALL folds
        rms <- RMSE(pred = val_out_df$ord_krig.predict, obs = val_out_df$VARIABLE, na.rm = T)
        r2 <- lm(formula = ord_krig.predict~VARIABLE, data = val_out_df)
        r2 <- c(summary(r2)$r.squared,summary(r2)$adj.r.squared)

        # add statistics
        par_new <- data.frame(par_new,rms_alldata=rms,r2_alldata=r2[1],adj_r2_alldata=r2[2])
        mydata_out <- rbind(mydata_out, par_new)
        
        # # write output (result of each fold computation)
        # print(paste("Write summary for distance ",k))
        # write.csv(x = mydata_fold,file = paste(var_name,"_x",k,"_10fold_",model,"_tab.csv",sep = ""), row.names = F,quote = F)
    
    # #
    # }# end loop on variogram distance
    # ###
    
    # write output (parameters as average on folds)
    dir.create(file.path(getwd(), var_name), recursive = T)
    print("Write total summary")
    write.csv(x = mydata_out,
              file = file.path(getwd(), var_name,
                               paste(var_name,"_",par[1],"_",par[2],"_",par[3],"_",model,".csv",sep = "")),
              row.names = F,quote = F )
    
    # #
    # var_tot <- variogram(log(VARIABLE)~1,data = worktab,locations = ~X+Y,cutoff = par[1])#, width=par[1]/100)
    # mod_tot <- vgm(psill = par_new$psill,model = par_new$model,nugget = par_new$nugget,range = par_new$range)
    # var_tot <- my_var # NON serve, basta solo il modello per interp. kriging
    # mod_tot <- vgm(psill = par_new$psill, model = par_new$model,
    #                nugget = par_new$nugget, range = par_new$rad_fit) )
    # plot(var_tot,model=mod_tot,main =paste("RMSE: ",rms) )
    ###
    
    # return(RMSE(pred = val_out_df$ord_krig.predict, obs = val_out_df$VARIABLE, na.rm = T))
    return(r2[1])
}

## keep care: trade of between search distance and number of NA estimations
## the smaller the search radius, the better the estimation - but lot of NAs
## How to solve?
## Solved because radius of variogram and interpolation are different! Variogram
## radius explain correlation between points -> will vary for each variables.
## Interpolation radius define the search area in which point where selected
## and used to interpolate a single pixel (following the model defined by
## variogram) -> could be choosen as fixed value
#
# hydroPSO::hydroPSO(fn = OrdKrig_optim_krige, method="spso2011",
#                    lower = c(0,0,0.01,8,1,0), upper = c(1000,359,1,100,25,10),
#                    control=list(drty.out = "/home/jbre/R/OrdKrig/PSO_krige", npart=40,
#                                 parallel="none", par.pkgs = c("gstat","caret","hydroGOF","sp")))
