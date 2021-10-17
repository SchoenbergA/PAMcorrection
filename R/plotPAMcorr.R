#' plotPAMcorr - Function for plotting of PAM correction.
#' @description visualizes the correction of PAM values by correction value.
#' @param df data.frame - with column "CTR" for control values, "PAM" for pam values
#' and optional "PAM_corr" for corrected PAM values.
#' @param al numeric - a vector of numeric values on the x axes to draw lines. Default=0
#' @param yl numeric - ylim adjustment. Default=c(0,2)
#' @param sortby character - column name to reorder the inoput datafram 'df'. IF set the argument 'al' will not be used.
#' @param abs_dif boolean - IF TRUE plots the differences as absolute values (no negative values). Default=TRUE.
#' @param titel character - desired maintitel for the plot . Default="PAMcorrection".
#' @return returns a plot containing the original PAM (orange) and Control (blue) values along with the difference (grey).
#' If corrected PAM values are available those will be plotted in red along with the differnece between the corrected PAM and COntrol values.
#' Additionally prints the total Difference along with mean and standart divation for PAM or if given PAM corrected and Control
#' @note For visualization purposes the difference between the corrected PAM and Control is used as absolute values (negative values are plotted as positive).
#' Further the mean and sd is calculated for the absolute values.
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
#' plotPAMcorr(dat,yl = c(-1,4))
#'
#' # use vertical lines
#' abl <- seq(10,40,10)
#' plotPAMcorr(dat,yl = c(0,4),al = abl)
#'
#' # order values by PAM
#' plotPAMcorr(dat,sortby = "PAM")
#'
#' # get correction value by mean
#' datf <-dat
#' datf$PAM_corr <-dat$PAM * mean(dat$CTR/dat$PAM) # factor
#'
#' # plot correction and use titel
#' plotPAMcorr(datf,sortby = "PAM",titel="PAM corr factor")# for factor
#'
#' # show difference in both poisitve and negative directions
#' plotPAMcorr(datf,sortby = "PAM",abs_dif=F,yl=c(-1,2))
#' # Note: The total and mean diffenrece is still calculcated in absolute values.






plotPAMcorr <- function(df,al=-1,yl=c(0,2),sortby=NULL,abs_dif=T,titel="PAM Correction"){

  if(is.null(sortby)==F){
    posby <-which(colnames(df)==sortby)
    df <- df[order(df[,posby]),]
    # plot
    cat(paste0("dataframe sorted by ",sortby),sep="\n")

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
    plot(ctr,col="blue",
         ylim=yl,
         main=paste0(titel),
         #sub =paste0("empty yet"),

         xlab=if(is.null(sortby)==F){paste0("sorted by ",sortby)}else{"Index"},

         ylab="Value")

    points((pamc),col="red")
    points(pam,col="orange")
    lines(ctr,col="blue")
    lines(pamc,col="red")
    lines(pam,col="orange")

    if(abs_dif==T){
      # div between control and corrected PAM
      div <- ctr-(pamc)
      lines(abs(div),col="green")
      points(abs(div),col="green")
      abline(h=mean(abs(div)),col="green")

    } else if (abs_dif==F){
      # div between control and corrected PAM
      div <- ctr-(pamc)
      lines((div),col="green")
      points((div),col="green")
      cat("plotting difference (green) including negative values",sep="\n")
      cat("plotting no mean fordifference (green)",sep="\n")
    }
      # div between control and org PAM
      divorg <- ctr-(pam)
      lines(abs(divorg),col="grey")
      points(abs(divorg),col="grey")
      abline(h=mean(abs(divorg)),col="grey")


    # set global ablines
    abline(h=0)
    abline(v=al)
    # labels
    if(is.null(pamc)){
      mtext(paste0("orgDiffSums (abs): ",round(sum(abs(divorg)),2)),line = -18, adj = 0.01)
      mtext(paste0("orgDiffMean (abs): ",round(mean(abs(divorg)),2)),line = -19, adj = 0.01)
      mtext(paste0("orgDiffSD (abs): ",round(sd(abs(divorg)),2)),line = -20, adj = 0.01)
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

