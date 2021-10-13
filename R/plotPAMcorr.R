#' plotPAMcorr - Function for plotting of PAM correction.
#' @description visualizes the correction of PAM values by correction value.
#' @param df data.frame - with column "CTR" for control values, "PAM" for pam values
#' and optional "PAM_corr" for corrected PAM values.
#' @param al numeric - a vector of numeric values on the x axes to draw lines. Default=0
#' @param yl numeric - ylim adjustment. Default=2
#' @param sortby - character - column name to reorder the inoput datafram 'df'. IF set the argument 'al' will not be used.
#' @param  titel - character - desired maintitel for the plot . Default="PAMcorrection".
#' @return returns a plot containing the original PAM (orange) and Control (blue) values along with the difference (grey).
#' If corrected PAM values are available those will be plotted in red along with the differnece between the corrected PAM and COntrol values.
#' Additionally prints the total Difference along with mean and standart divation for PAM or if given PAM corrected and Control
#' @details
#' @note For 'cf' it is recommended to calculate the mean of the difference
#' between the original PAM and control values. Either use devision Control/PAM
#' for factor and abs(PAM)-abs(Control) for correction by absolute values.
#' @author Andreas Schönberg
#' @export plotPAMcorr
#' @aliases plotPAMcorr
#' @examples
#' # load data
#' dat <- read.csv(system.file("extdata","exp_PAM.csv",package = "PAMcorrection"))
#' head(dat)
#'
#' # plot control and PAM
#' plotPAMcorr(dat)
#'
#' # adjust y limit
#' plotPAMcorr(dat,yl = 4)
#'
#' # use vertical lines
#' abl <- seq(10,40,10)
#' plotPAMcorr(dat,yl = 4,al = abl)
#'
#' # order values by PAM
#' plotPAMcorr(dat,sortby = "PAM")
#'
#' # get correction value by mean
#' datf <-dat
#' datf$PAM_corr <-dat$PAM * mean(dat$CTR/dat$PAM) # factor
#' data <-dat
#' data$PAM_corr <-dat$PAM * mean(dat$PAM-dat$CTR) # factor
#'
#' # plot correction and use titel
#' plotPAMcorr(datf,sortby = "PAM",titel="PAM corr factor")# for factor
#' plotPAMcorr(data,sortby = "PAM",titel="PAM corr absolute") # for absolute




plotPAMcorr <- function(df,al=-1,yl=2,sortby=NULL,titel="PAM Correction"){

  if(is.null(sortby)==F){
    posby <-which(colnames(df)==sortby)
    df <- df[order(df[,posby]),]
    # plot
    cat(paste0("dataframe sorted by ",sortby),sep="\n")

    al= -1
    cat("set 'al' to -1 due to sorting data",sep="\n")
  }


  # get data
  ctr <- df$CTR
  pam <- df$PAM
  pamc<- df$PAM_corr
  if(is.null(pamc)){
    cat("no corrected PAM detected")
    titel<-paste0(titel," uncorrected")
    }
    # plot
    cat("using original input PAM (not corrected by function)",sep="\n")
    plot(ctr,col="blue",
         ylim=c(0,yl),
         main=paste0(titel),
         #sub =paste0("empty yet"),

         xlab=if(is.null(sortby)==F){paste0("sorted by ",sortby)}else{"ID"},
         ylab="Value")

    points((pamc),col="red")
    points(pam,col="orange")
    lines(ctr,col="blue")
    lines(pamc,col="red")
    lines(pam,col="orange")

    # div between control and corrected PAM
    div <- ctr-(pamc)
    lines(abs(div),col="green")
    points(abs(div),col="green")
    # div between control and org PAM
    divorg <- ctr-(pam)
    lines(abs(divorg),col="grey")
    points(abs(divorg),col="grey")
    # set ablines
    abline(h=0)
    abline(v=al)
    # labels
    if(is.null(pamc)){
      mtext(paste0("orgDiffSums: ",round(sum(abs(divorg)),2)),line = -18, adj = 0.01)
      mtext(paste0("orgDiffMean: ",round(mean(abs(divorg)),2)),line = -19, adj = 0.01)
      mtext(paste0("orgDiffSD: ",round(sd(abs(divorg)),2)),line = -20, adj = 0.01)
    } else {
    mtext(paste0("absDiffSums: ",round(sum(abs(div)),2)),line = -18, adj = 0.01)
    mtext(paste0("absDiffMean: ",round(mean(abs(div)),2)),line = -19, adj = 0.01)
    mtext(paste0("absDiffSD: ",round(sd(abs(div)),2)),line = -20, adj = 0.01)
    }
    # legend
    if(is.null(pamc)){
    legend("topleft", legend=c("PAM", "Control","org Difference"),
             col=c("orange","blue","grey"), lty=1, cex=0.8)
    } else {
    legend("topleft", legend=c("PAM","PAM corrected", "Control","Difference","org Difference"),
           col=c("orange","red","blue", "green","grey"), lty=1, cex=0.8)
    }


}# end of function

