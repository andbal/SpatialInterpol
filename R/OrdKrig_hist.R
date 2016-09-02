################################################################################
# Name:         check_data
# Author:       Andrea Balotti - balotti.and@gmail.com
# Date:         08/2016 - creation
# Copyright:    Andrea Balotti (c) 2016

### Description:
# Script to visually check if data are suitable for geostatistical method:
# 1) normally distribuited
# 2) stationary (mean and variance do not vary significantly in space)

### Reference:
# [1] G.Bohling, INTRODUCTION TO GEOSTATISTICS And VARIOGRAM ANALYSIS, C&PE 940, 17 October 2005
#     http://people.ku.edu/~gbohling/cpe940/Variograms.pdf
# [2] G.Bohling, KRIGING, C&PE 940, 19 October 2005
#     http://people.ku.edu/~gbohling/cpe940/Kriging.pdf
################################################################################

# Input path and file
# wpath = "H:/Projekte/MONALISA/05_Arbeitsbereiche/BaA/05_Soil_Interpolation/02_additional_maps"
# datafile = "master/original_dataset/Masterfile_AdigeVenosta.txt"

OrdKrig_hist <- function (wpath = "tmp/path",
                          datafile = "master/file.txt" )
{
    vars = list("Humus____","pH_Wert__i","Karbonate_","Kalkbedarf","Phosphat__",
                "Kali__K_O_","Magnesium_","Bor__B__im","Mangan__Mn","Kupfer__Cu",
                "Zink__Zn__")
    
    for (var in vars){
    
        # define input data table 
        worktab <- read.table(file = file.path(wpath, datafile), header = TRUE, sep = ",",dec = ".")
        worktab <- cbind(worktab$x_Coord, worktab$y_Coord, worktab[,var])
        worktab <- as.data.frame(worktab)
        names(worktab) <- c("X","Y","VARIABLE")
        
        # zeros
        zero <- which(worktab$VARIABLE <= 0)
        if (length(zero) > 0){
            print(paste(length(zero),"Values less or equal to ZERO found for",var,". Will be reassigned as 0.001") )
            worktab[worktab$VARIABLE <= 0,"VARIABLE"] <- 0.001
        }
        
        # Input data shortcut
        t=worktab$VARIABLE
        
        ### 2 plot in a frame
        par(oma=c(0,0,1,0))
        par(mfrow=c(1,2))
        # plot histogram (play with breaks value to obtain nice plot)
        hist(t,breaks = nrow(worktab)/500,freq = F,
             col = "gray", main = "Histogram vs. Normal Distribution",
             xlab = "Variable value") #,xlim=c(min(t),max(t)))
        
        # plot gaussian
        m <- mean(t,na.rm = T)      # mean
        std <- sd(t,na.rm = T)      # std.dev
        mingau=max(min(t),m-3*std)  # inf limit as min data or 3 std.dev. from mean
        maxgau=min(max(t),m+3*std)  # sup limit as max data or 3 std.dev. from mean
        curve(dnorm(t,m,std),xname = "t",from = mingau,to = maxgau,add=T,col="red",lwd=2)
        
        # plot QQ
        qqnorm(t,main = "Normal Q-Q plot")
        qqline(t)
        title(var, outer=T)
    }
}