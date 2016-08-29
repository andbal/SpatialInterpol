# Function containing tools for various secondary .
# According to 'use' parameter, the function output could be:
# 1) 'use' == 'fold'
#    exclude from the original dataset a user-specified percentual of data and
#    create a new "cutted" dataset that could be used for training/interpolation
# 2) 'use' == 'validation'
#    given the original (full) dataset, perform 

OrdKrig_tools <- function (wpath = "/home/jbre/R/OrdKrig",
                           datafolder = "master", rastermask = "mask/Mask_master.tif",
                           use = "validation", perc=0.8, kfold = 5  )
{ # main func start
    
    ###################################
    filesIN <- dir(file.path(wpath, datafolder))
    if (is.null(names(cutoff)))
    {
        # get namezones according to input files 
        namezones <- c()
        for (namefile in filesIN)
        {
            nr_char <- nchar(filesIN)
            names(nr_char) <- filesIN
            namezone <- strsplit(substr(namefile,1,nr_char[namefile]-4), "_")[[1]][2]
            namezones <- c(namezones, namezone)
            }
        
        # list of arguments
        args <- list("cutoff"=cutoff, "anis_deg"=anis_deg,"anis_ax"=anis_ax,"psill"=psill,
                     "nugget"=nugget,"nmax"=nmax,"nmin"=nmin,"idp"=idp)
        
        # assing argument names to parameter value (named vector)
        for (arg in 1:length(args))
        {
            x <- args[[arg]]
            if (length(x)==length(namezones)) {
                names(x) <- namezones
                assign(x = names(args)[arg], value = x) 
            } else {
                print(paste("argument",  names(args)[arg], "is not given with name of investigated zone and differs in length from input .csv files in the folder", datafolder, ". Please check!", sep=" " ))
                }
            }
        }
    
    val_list <- list()
    for (namefile in filesIN)
    {
        nr_char <- nchar(filesIN)
        names(nr_char) <- filesIN
        namezone <- strsplit(substr(namefile,1,nr_char[namefile]-4), "_")[[1]][2]
        print(paste("processing zone", namezone, "for variable", variable, sep=" "))
        # read table 
        worktab <- read.table(file = file.path(wpath, datafolder, namefile), header = TRUE, sep = ",",dec = ".")
        worktab <- cbind(worktab$x_Coord,worktab$y_Coord,worktab[,variable])
        # matrix 2 data.frame
        worktab <- as.data.frame(worktab)
        # rename cols
        names(worktab) <- c("X","Y","VARIABLE")
        # zeros
        worktab[worktab$VARIABLE <= 0,"VARIABLE"] <- 0.001    
    }
    
    #######################


# check the working mode
use_type=c("fold","validation","mask")
if (!use %in% use_type)
{
    stop("'use' parameter is incorrect. Choose between 'fold','validation' or 'mask'")
}
# check parameter for working mode 'fold'
if (use == use_type[1])
{
    if (!is.numeric(perc) | ( perc < 0 | perc > 1 )
    {
        stop("'perc' parameter is incorrect. Must be in the interval [0,1]")
        }
    
}
# check parameter for working mode 'validation'
if (use == use_type[2])
{
    }
# check parameter for working mode 'mask'
if (use == use_type[3])
{
    }

### 1: Fold
if (use == use_type[1])
{
    # Set user sample size
    size <- floor(perc * nrow(worktab))
    ## set the seed to make your partition reproductible
    #set.seed(123)
    train_ind <- sample(seq_len(nrow(worktab)), size = size)
    train <- worktab[train_ind, ]
    #test <- mtcars[-train_ind, ] # remaining part, not used
    return(train)
    }

### 2: Validation
if (use == use_type[2])
{
    flds <- createFolds(y = worktab$VARIABLE, k = kfold, list = TRUE, returnTrain = FALSE)
    val_out_var <- list()
    val_out_df  <- data.frame()
    
    for (i in 1:length(flds))
    {
        # training
        train_set <- worktab[-flds[[i]],]
        
        # train
        if (!inverseDistWeigths)
        {
            # gstatVariogram - Calculate Sample variogram 
            my_var <- variogram(log(VARIABLE)~1, data=train_set, locations = ~X+Y, cutoff = cutoff[namezone])
            # Fit a Variogram Model to a Sample Variogram  
            m <- vgm(psill = psill[namezone], model = model, range = cutoff[namezone], nugget = nugget[namezone], 
                    anis = c(anis_deg[namezone], anis_ax[namezone]))
            my_var_fit <- fit.variogram(my_var, m)
        }
        
        # validation set
        valid_set <- worktab[flds[[i]],]
        
        Xnew <- valid_set[,c("X", "Y")]
        Xnew <- SpatialPoints(Xnew)
        
        myloc <- data.frame("X" = train_set$X,"Y" = train_set$Y)
        myloc <- SpatialPoints(myloc)
        
        if (inverseDistWeigths) {
            
            # Inverse Distance Weighting (local only)
            ord_krig <- gstat::idw(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew, idp = idp[namezone],
                                  nmax = nmax[namezone], nmin = nmin[namezone], omax = omax[namezone], maxdist = cutoff[namezone])
            
            names(ord_krig) <- c("predict")
            
            val_out_df <- rbind(val_out_df, data.frame(fold=i ,valid_set, ord_krig$predict))
    
        } else {
            
            # Ordinary Kriging (local only)
            ord_krig <- gstat::krige(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew, model = m,
                                    nmax = nmax[namezone], nmin = nmin[namezone], omax = omax[namezone], maxdist = m$range[2])
            
            names(ord_krig) <- c("predict", "variance")
            
            val_out_var[[paste("fold", i, sep="")]] <- list(vario = my_var, vario_fit = my_var_fit)
            val_out_df <- rbind(val_out_df, data.frame(fold=i ,valid_set, ord_krig$predict, ord_krig$variance))
        }
    }
    val_list[[namezone]] <- list(variogram = val_out_var, df = val_out_df)
    return(val_list[[namezone]])
}

### 3: Mask
if (use == use_type[2])
{
    }

} #end main func