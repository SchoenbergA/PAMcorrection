#' estPAMcorr2 - advanced Function for estimation of PAM correction.
#' @description visualizes the correction of PAM values by correction value.
#' @param ctr data.frame column - with Control values.
#' @param pam data.frame column - with PAM values.
#' @param cf numeric - desired correction factor.
#' @param al numeric - a vector of numeric values on the x axes to draw lines. Default=0
#' @param modef bolean - If TRUE uses input 'cf' to multiply PAM, if FALSE uses
#' substracts 'cf' from PAM. Default= TRUE
#' @param yl numeric - ylim adjustment. Default=2

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
#' @export estPAMcorr2
#' @aliases estPAMcorr2
#' @examples
#' # load data
#' dat <- read.csv(system.file("extdata","exp_PAM.csv",package = "PAMcorrection"))
#' head(dat)
#'
#' # plot control and PAM
#' estPAMcorr2(dat)
#'
#' # adjust y limit
#' estPAMcorr2(dat,yl = 4)
#'
#' # use vertical lines
#' abl <- seq(10,40,10)
#' estPAMcorr2(dat,yl = 4,al = abl)
#'
#' # order values by PAM
#' estPAMcorr2(dat,sortby = "PAM")
#'
#' # get correction value by mean
#' mf <-mean(dat$CTR/dat$PAM) # factor
#' ma <-mean(dat$PAM-dat$CTR) # absolute
#'
#' # plot correction
#' estPAMcorr2(dat,sortby = "PAM",cf = mf)# for factor
#' # for absolute use 'modef=F'
#' estPAMcorr2(dat,sortby = "PAM",cf = ma,modef = F) # for absolute




estPAMcorr2 <- function(df,cf=NULL,al=-1,modef=T,yl=2,sortby=NULL){

  # check input
  #mean(ctr/pam)
  #if(modef==T & cf>)



  # handle sort

  if(is.null(sortby)==F){
    posby <-which(colnames(df)==sortby)
    df <- df[order(df[,posby]),]

    cat(paste0("dataframe sorted by ",sortby),sep="\n")

    al= -1
    cat("set 'al' to -1 due to sorting data",sep="\n")

  }

  # get data
  ctr <- df$CTR
  pam <- df$PAM

  #  using factor
  if(modef==T){

    # without correction
    if(missing(cf)){
      cat("using original input PAM (not corrected by function)",sep="\n")
      plot((ctr),col="blue",
           ylim=c(0,yl),
           main=paste0("PAM"),
           #sub =paste0("empty yet"),
           xlab="ID",
           ylab="Value")

      #points((pam*cf),col="red")
      points(pam,col="orange")
      lines(ctr,col="blue")

      #lines(pam*cf,col="red")
      lines(pam,col="orange")


      div <- ctr-(pam)

      lines(abs(div),col="green")
      points(abs(div),col="green")
      abline(v=al)
      # labels
      mtext(paste0("absDiffSums: ",round(sum(abs(div)),2)),line = -18, adj = 0.01)
      mtext(paste0("absDiffMean: ",round(mean(abs(div)),2)),line = -19, adj = 0.01)
      mtext(paste0("absDiffSD: ",round(sd(abs(div)),2)),line = -20, adj = 0.01)
      # legend
      legend("topleft", legend=c("PAM", "Control","Difference"),
             col=c("orange","blue", "green"), lty=1, cex=0.8)

    } else {
      # if cf is given
      cat("using factor for correction",sep="\n")
      plot((ctr),col="blue",
           ylim=c(0,yl),
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
      mtext(paste0("absDiffSD: ",round(sd(abs(div)),2)),line = -20, adj = 0.01)
      # legend
      legend("topleft", legend=c("PAM", "PAM_corrected","Control","Difference"),
             col=c("orange","red","blue", "green"), lty=1, cex=0.8)
    }
  }
  # using absolut values
  if (modef==F){
    if(missing(cf)){# without correction

      cat("using original input PAM (not corrected by function)",sep="\n")
      plot((ctr),col="blue",
           ylim=c(0,yl),
           main=paste0("PAM"),
           #sub =paste0("empty yet"),
           xlab="ID",
           ylab="Value")

      #points((pam-cf),col="red")
      points(pam,col="orange")
      lines(ctr,col="blue")

      #lines(pam-cf,col="red")
      lines(pam,col="orange")


      div <- ctr-(pam)

      lines(abs(div),col="green")
      points(abs(div),col="green")
      abline(v=al)
      # labels
      mtext(paste0("absDiffSums: ",round(sum(abs(div)),2)),line = -18, adj = 0.01)
      mtext(paste0("absDiffMean: ",round(mean(abs(div)),2)),line = -19, adj = 0.01)
      mtext(paste0("absDiffSD: ",round(sd(abs(div)),2)),line = -20, adj = 0.01)
      # legend
      legend("topleft", legend=c("PAM", "Control","Difference"),
             col=c("orange","blue", "green"), lty=1, cex=0.8)


    }else {
      cat("using absolute value for correction",sep="\n")
      plot((ctr),col="blue",
           ylim=c(0,yl),
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
      mtext(paste0("absDiffSD: ",round(sd(abs(div)),2)),line = -20, adj = 0.01)
      # legend
      legend("topleft", legend=c("PAM", "PAM_corrected","Control","Difference"),
             col=c("orange","red","blue", "green"), lty=1, cex=0.8)
    }
  }
}# end of function
