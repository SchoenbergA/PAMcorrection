#' estPAMcorr - Function for estimation of PAM correction.
#' @description visualizes the correction of PAM values by correction value.
#' @param ctr data.frame column - with Control values.
#' @param pam data.frame column - with PAM values.
#' @param cf numeric - desired correction factor.
#' @param al numeric - a vector of numeric values on the x axes to draw lines. Default=0
#' @param modef bolean - If TRUE uses input 'cf' to multiply PAM, if FALSE uses
#' substracts 'cf' from PAM. Default= TRUE

#' @return returns a plot containing the original PAM values, control values and
#' corrected PAM values along the absolute difference between corrected PAM and
#' control values (green). Further returns mean and sum of this difference for
#' quality estimation of the result.
#' @details
#'  * cf - it is recommended to calculate the mean of the diffenrence
#' between the original PAM and control values. Either use devision Control/PAM
#' for correction by factor and abs(PAM)-abs(Control) for correction by absolute values.
#'  * al - is used to visualise the borders of differnet classes for the values
#' @note For 'cf' it is recommended to calculate the mean of the diffenrence
#' between the original PAM and control values. Either use devision Control/PAM
#' for factor and abs(PAM)-abs(Control) for correction by absolute values.
#' @author Andreas Schönberg
#' @export estPAMcorr
#' @aliases estPAMcorr
#' @examples
#' # load data
#' dat <- read.csv(system.file("extdata","exp_PAM.csv",package = "PAMcorrection"))
#' head(dat)
#'
#' # visualize Control and PAM values
#' estPAMcorr(ctr = dat$CTR,pam = dat$PAM)
#'
#' # view borders of dialect
#' estPAMcorr(ctr = dat$CTR,pam = dat$PAM,al = c(9,15,21,36))
#'
#' # use a correction factor
#' cf <- mean(dat$CTR/dat$PAM)
#'
#' estPAMcorr(ctr = dat$CTR,pam = dat$PAM,al = c(9,15,21,36),cf = cf)
#'
#' # use absolute value for correction
#' cfa <- abs(mean(dat$CTR-dat$PAM))
#'
#' estPAMcorr(ctr = dat$CTR,pam = dat$PAM,al = c(9,15,21,36),cf = cfa,modef = FALSE)



estPAMcorr <- function(ctr,pam,cf=NULL,al=-1,modef=T){

  # check input
  #mean(ctr/pam)
  #if(modef==T & cf>)

  # abs or factor
  if(modef==T){
    if(missing(cf)){
      cf=1
    }
    cat("using factor for correction",sep="\n")
    plot((ctr),col="blue",
         ylim=c(0,2),
         main=paste0("PAM Correction factor ",cf),
         #sub =paste0("empty yet"),
         xlab="ID",
         ylab="Value")

    points((pam*cf),col="red")
    points(pam,col="orange")
    lines(ctr,col="blue")

    lines(pam*cf,col="red")
    lines(pam,col="orange")


    div <- ctr-(pam*cf)

    lines(abs(div),col="green")
    points(abs(div),col="green")
    abline(v=al)
    # labels
    mtext(paste0("absDiffSums: ",round(sum(abs(div)),2)),line = -18, adj = 0.01)
    mtext(paste0("absDiffMean: ",round(mean(abs(div)),2)),line = -19, adj = 0.01)
    # legend
    legend("topleft", legend=c("PAM", "PAM_corrected","Control","Difference"),
           col=c("orange","red","blue", "green"), lty=1, cex=0.8)

  }
  if (modef==F){
    if(missing(cf)){
      cf=0
    }
    cat("using absolute value for correction",sep="\n")
    plot((ctr),col="blue",
         ylim=c(0,2),
         main=paste0("PAM Correction absolut -",cf),
         #sub =paste0("empty yet"),
         xlab="ID",
         ylab="Value")

    points((pam-cf),col="red")
    points(pam,col="orange")
    lines(ctr,col="blue")

    lines(pam-cf,col="red")
    lines(pam,col="orange")


    div <- ctr-(pam-cf)

    lines(abs(div),col="green")
    points(abs(div),col="green")
    abline(v=al)
    # labels
    mtext(paste0("absDiffSums: ",round(sum(abs(div)),2)),line = -18, adj = 0.01)
    mtext(paste0("absDiffMean: ",round(mean(abs(div)),2)),line = -19, adj = 0.01)
    # legend
    legend("topleft", legend=c("PAM", "PAM_corrected","Control","Difference"),
           col=c("orange","red","blue", "green"), lty=1, cex=0.8)
  }
}# end of function
