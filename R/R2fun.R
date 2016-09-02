### R2 for variogram
# Description: Function to compute correlation coefficient after fit.variogram
# Copyright: code found on web

### Reference:
# http://r-sig-geo.2731867.n2.nabble.com/R2-values-from-SSErr-fit-variogram-attribute-td6167483.html
# https://stat.ethz.ch/pipermail/r-sig-geo/2011-March/011108.html
# http://grokbase.com/t/r/r-sig-geo/109242k2dm/how-to-calculate-r2-coefficient-of-determination-for-variogram

### INPUT:
# vario = variogram output
# fit = fit.variogram output
# fit.method = choose from [1,2,6,7] the fit.method used also for fit.variogram
#   1 ->    weights : $N_j$
#               with $N_j$ the number of point pairs.
#   2 ->    weights : $N_j/{γ(h_j)}^2$
#               with $N_j$ the number of point pairs and
#               $γ(h_j)$ the variance of the group of point-pairs.
#   6 ->    this method uses no weights
#   7 ->    weights : $N_j/h_j^2$ [DEFAULT]
#               with $N_j$ the number of point pairs and $h_j$ the distance.

R2fun <- function(vario, fit, fit.method=7)
{
    if(!is(vario, "gstatVariogram")) stop("Fit must be a gstatVariogram!\n")
    if(!is(fit, "variogramModel")) stop("Fit must be a variogramModel!\n")
    if(length(fit.method)!=1) stop("One fit.method must be selected!\n")
    if(!any(fit.method==c(1,2,6,7))) stop("The selected fit.method is not managed!\n")
    SSErr<-attr(fit,"SSErr")
    if(fit.method==1){
        weig<-vario$np
        SSTot<- sum( (weig* (vario$gamma-mean(vario$gamma)))^2 )
    }
    if(fit.method==2){
        weig<-vario$np/vario$gamma^2
        SSTot<- sum( (weig * (vario$gamma-mean(vario$gamma)))^2 )
    }
    if(fit.method==6){
        SSTot<- sum((vario$gamma-mean(vario$gamma))^2 )
    }
    if(fit.method==7){
        weig<-vario$np/vario$dist^2
        SSTot<- sum( (weig * (vario$gamma-mean(vario$gamma)))^2 )
    }
    SSReg <- SSTot-SSErr
    R2 <- 1-SSErr/SSTot
    # cat("\nFind2R values -> SSErr: ", SSErr, " SSReg: ", SSReg, " SSTot: ", SSTot, " R2: ", R2, "\n")
    return(R2)
}
