# Optimistation function for ordinary kriging

# choose variables to use from colnames of namefile (see below)
# "Bodenart__" "Humus____"  "pH_Wert__i" "Karbonate_" "Kalkbedarf" "Phosphat__" "Kali__K_O_" "Magnesium_" "Bor__B__im" "Mangan__Mn" "Kupfer__Cu" "Zink__Zn__" "NUM"       
# "ID_string"  "ID_suolo"   "Num_Soil"   "Soil_newCl" "Schluff"    "Tonig"      "Sand"

# wpath = "H:/Projekte/MONALISA/05_Arbeitsbereiche/BaA/05_Soil_Interpolation/02_additional_maps"
# datafile = "master/original_dataset/Masterfile_AdigeVenosta.txt"

library(gstat)
library(caret)
# library(hydroPSO)
library(sp)

# OrdKrig_optim_krige <- function(par = c(cutoff=3000, nmax=12, omax=3, psill=1, nugget=1),
#                                 wpath = "/home/jbre/R/OrdKrig", 
#                                 datafile = "master/Masterfile_AdigeVenosta.txt",
#                                 variable = "Humus____", local  =TRUE,
#                                 model = c("Sph","Exp"), kfold = 10 )
# {
    # to comment for the functiom
    par = c(cutoff=3000, nmax=12, nmin=1, omax=3, psill=1, nugget=1)
    wpath = "H:/Projekte/MONALISA/05_Arbeitsbereiche/BaA/05_Soil_Interpolation/02_additional_maps"
    datafile = "master/original_dataset/Masterfile_AdigeVenosta.txt"
    variable = "Kali__K_O_"
    model=c("Sph","Exp")
    kfold=10
    local=TRUE
    ###
    
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
  
  val_fold_df <- data.frame()
  val_out_df  <- data.frame()
  mydata <- data.frame()
  var_name = strsplit(x = variable,split = "_")[[1]][1]
  
  for (i in 1:length(flds))
  {
    # training set
    train_set <- worktab[-flds[[i]],]
    # validation set
    valid_set <- worktab[flds[[i]],]
    
    # gstatVariogram - Calculate experimental variogram and fit a model on it
    if (local){
        # local variogram 
        my_var <- variogram(log(VARIABLE)~1, data=train_set, locations = ~X+Y, cutoff = par[1], width=par[1]/100)
        # m <- vgm(psill = par[5], model = model, range = par[1], nugget = par[6])
        m <- vgm(model)
        my_var_fit <- fit.variogram(my_var, m)
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
        ord_krig <- krige(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew, model = my_var_fit,
                          nmax = par[2], nmin = par[3], omax = par[4], maxdist = par[1] ) #maxdist = my_var_fit$range[2]
    } else {
        ord_krig <- krige(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew, model = my_var_fit )
    }
    names(ord_krig) <- c("predict", "variance")
    # frame with original data + predicted for ACTUAL fold
    val_fold_df <- data.frame(fold=i,valid_set, ord_krig$predict, ord_krig$variance)
    # RMSE for ACTUAL fold
    rms_fold <- RMSE(pred = val_fold_df$ord_krig.predict, obs = val_fold_df$VARIABLE, na.rm = T)
    # frame with original data + predicted for ALL folds
    val_out_df <- rbind(val_out_df, val_fold_df)

    mydata <- rbind(mydata, data.frame(i, my_var_fit$model[2],
                                       my_var_fit$psill[2], my_var_fit$psill[1],
                                       my_var_fit$range[2], rms_fold) )
  }# end loop on folds
  
  names(mydata) <- c("fold","model","psill","nugget","range","RMSE")    
  # Create set of parameters as average of values got from 10 folds
  par_new <- as.numeric(apply(X = mydata[-c(1:2)],MARGIN = 2,FUN = mean))
  par_new <- data.frame(as.character(my_var_fit$model[2]), t(par_new),
                   as.numeric(par[2]), as.numeric(par[3]), as.numeric(par[4]),stringsAsFactors = F )
  colnames(par_new) <- c("model","psill","nugget","range","rmse","nmax","nmin","omax")
  write.csv(x = par_new,file = paste(var_name,"_10_fold.csv",sep = ""), row.names = F,quote = F)
  # RMSE for ALL folds
  rms <- RMSE(pred = val_out_df$ord_krig.predict, obs = val_out_df$VARIABLE, na.rm = T)
  
  var_tot <- variogram(log(VARIABLE)~1,data = worktab,locations = ~X+Y,cutoff = par[1])#, width=par[1]/100)
  mod_tot <- vgm(psill = par_new$psill,model = par_new$model,nugget = par_new$nugget,range = par_new$range)
  plot(var_tot,model=mod_tot,main =paste("RMSE: ",rms) )
  
#   return(RMSE(pred = val_out_df$ord_krig.predict, obs = val_out_df$VARIABLE, na.rm = T))
#   
# }

# # keep care: trade of between search distance and number of NA estimations
# # the smaller the search radius, the better the estimation - but lot of NAs
# # How to solve?
# 
# hydroPSO::hydroPSO(fn = OrdKrig_optim_krige, method="spso2011",
#                    lower = c(0,0,0.01,8,1,0), upper = c(1000,359,1,100,25,10),
#                    control=list(drty.out = "/home/jbre/R/OrdKrig/PSO_krige", npart=40, 
#                                 parallel="none", par.pkgs = c("gstat","caret","hydroGOF","sp")))
