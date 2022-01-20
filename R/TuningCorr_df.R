#' TuningCorr_df - Automated correction of PAM using tuned mean for specified attributes.
#' @description uses the mean based on the data clipped by 'Whisker' and 'Quantil' for specified subgroups to correct the PAM values.
#' @param df data.frame - containing Control and PAM values with colnames "CTR" and "PAM"
#' @param att1 character - Column name of data frame which should be used to correct by mean. If missing, will
#' use the global mean to correct PAM.
#' @param att2 character - Additional column name of data frame which should be used to correct by mean.
#' Only necessary if the correction should use subgroups for the data in 'pos1'
#' @param tuning character - method for clipping the data used to calculate the mean. Methods can be "whisker", "quantil" or "threshold". If 'NULL' will perform no tuning. Default= NULL.
#' @param threshold numeric - specifies the threshold to clip values for tuning. Only needed if 'tuning' == "threshold".
#' @return returns a data.frame with an additional column "PAM_corr" with the corrected PAM values
#' and prints the respective mean values for each class or combination of classes. Uses 'tuned' data to calculate the mean values.
#' @note If used without any attributes will use the global mean over all values to correct PAM. IF 'att1' is set
#' will calculate the mean for each group in 'att1' to correct PAM. IF 'att2' is given will use the mean for each group
#' with the combination for each attribute.
#' @details * tuning - available modes are "whisker" (which clips all values </> lower / upper whisker threshold in boxplot) and "quantil" (which clips all values </> 25 and 75 quantils).
#' "threshold" will clip all values > a specified threshold.
#' * If a specific class or combination of attributes has no values left after tuning the mean will be calculated based on the original data. This will lead to a warning.
#' The class or combination is marked as 'noData' in the result table.
#' * If a class or combination has a 'NaN' the class or combination is not existing in the dataframe. This is possible when using att1 and att2 an if not every att1 has all att2 classes.
#' @author Andreas Schönberg
#' @export TuningCorr_df
#' @aliases TuningCorr_df
#' @examples
#' # load data
#' dat <- read.csv(system.file("extdata","exp_PAM.csv",package = "PAMcorrection"))
#' head(dat)
#'
#' # how the tuning works ("whisker" and "quantil")
#' bp <-boxplot(dat$PAM-dat$CTR)
#' bp$stats[1,1] # lower whisker
#' bp$stats[2,1] # lower quantil
#' bp$stats[3,1] # median
#' bp$stats[4,1] # upper quantil
#' bp$stats[5,1] # upper whisker
#'
#' # correct df by global means
#' corrected <-TuningCorr_df(dat,tuning = "whisker")
#' corrected <-TuningCorr_df(dat,tuning = "quantil")
#'
#' # use one class (no tuning)
#' corrected <-TuningCorr_df(dat,att1 = "dilect")
#' # use one class with threshold tuning
#' corrected <-TuningCorr_df(dat,att1 = "dilect",tuning="threshold",threshold=0.8)
#' # use two classes (means for each generation depending on dilect) tuned by 'whisker'
#' corrected <-TuningCorr_df(dat,att1 = "generation",att2 = "dilect",tuning = "whisker")



TuningCorr_df <-function(df,att1=NULL,att2=NULL,tuning=NULL,threshold=NULL){

  # check input
  if(is.null(att1) & is.null(att2)==F){
    stop("missing 'att1' but given 'att2'")
  }

  ##############################################################################
  # 1. Tuning
  if (is.null(tuning)){
    # get copy with no changes
    dft <-df
    cat("No tuning - using original input data",sep="\n")
    } else if(tuning=="whisker"){
    # copy df
    dft <- df
    # get absolut difference
    dft$diffa <- abs(dft$PAM-dft$CTR)
    # get boxplot
    bp <-boxplot(dft$diffa)
    # clip values by whisker
    boxplot(dft$diffa)
    #print(dft$diffa<bp$stats[1,1] | dft$diffa>bp$stats[5,1])
    dft <-dft[!(dft$diffa<bp$stats[1,1] | dft$diffa>bp$stats[5,1]),] # clip out of whisker
    cat(paste0("Tuning data by ",tuning),sep="\n")
    cat(paste0("clipping ",nrow(df)-nrow(dft), " entries",sep="\n"))
  } else if(tuning=="quantil"){
    # copy df
    dft <- df
    # get absolut difference
    dft$diffa <- abs(dft$PAM-dft$CTR)
    # get boxplot
    bp <-boxplot(dft$diffa)
    # clip values by whisker
    boxplot(dft$diffa)
    #print(dft$diffa<bp$stats[2,1] | dft$diffa>bp$stats[4,1])
    dft <-dft[!(dft$diffa<bp$stats[2,1] | dft$diffa>bp$stats[4,1]),] # clip out of quantil
    cat(paste0("Tuning data by ",tuning),sep="\n")
    cat(paste0("clipping ",nrow(df)-nrow(dft), " entries",sep="\n"))
  } else if(tuning=="threshold"){
    # copy df
    dft <- df
    # get absolut difference
    dft$diffa <- abs(dft$PAM-dft$CTR)
    dft <-dft[!(dft$diffa>threshold),] # clip out of whisker
    cat(paste0("Tuning data by threshold ",threshold),sep="\n")
    cat(paste0("clipping ",nrow(df)-nrow(dft), " entries",sep="\n"))
  }


  ##############################################################################


  # 2. Global mean #############################################################
  # if no attribute is given, calculate global mean
  if(is.null(att1)){

      df$PAM_corr <-df$PAM * mean(dft$CTR / dft$PAM)
      cat(paste0("global mean (factor) "),mean(dft$CTR / dft$PAM))
      return(df)

  } # end global mean
  ##############################################################################

  ##############################################################################
  # Prepare col position of attributes and init dataframes

  # get col position of attributes
  if(is.null(att2)){
  pos1 <- which(colnames(df)==att1)
  } else {
    pos1 <- which(colnames(df)==att1)
    pos2 <- which(colnames(df)==att2)
  }

  # init dataframe to store results
  ck <-data.frame() # to check for missing data
  re <-data.frame() # store resulting mean vlaues for class
  ##############################################################################


  # 3. Mean factor #############################################################

  # for factor correction
  if(is.null(att2)){
  # correct by one attribute
  cat("using att 1 only",sep="\n")
  # get uniques for each pos
  u1 <- unique(df[,pos1])
  # inner loop for pos1
  for (i in 1:length(u1)) {
    # get results to check for missing results
    ck <- rbind(ck,c(length(dft$PAM[ dft[,pos1]==u1[i] ]),u1[i],round(mean(dft$CTR[ dft[,pos1]==u1[i] ] / dft$PAM[ dft[,pos1]==u1[i] ]),digits = 4)))
    if(ck[nrow(ck),3]=="NaN"){
      #cat(paste0("Class '",re[nrow(ck),2],"' has no values left after tuning, using original data"),sep="\n")
      # use original data
      df$PAM_corr[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ] * mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ])
      # store in df
      re <- rbind(re,c(length(df$PAM[ df[,pos1]==u1[i] ]),u1[i],round(mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ]),digits = 4),"noData"))
    } else {
    df$PAM_corr[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ]    * mean(dft$CTR[ dft[,pos1]==u1[i] ] / dft$PAM[ dft[,pos1]==u1[i] ])
    re <- rbind(re,c(length(dft$PAM[ dft[,pos1]==u1[i] ]),u1[i],round(mean(dft$CTR[ dft[,pos1]==u1[i] ] / dft$PAM[ dft[,pos1]==u1[i] ]),digits = 4),"fine"))
    } # end fork for checking
     }# end i loop

  # return
  #cat("Platzhalter1",sep="\n")
  colnames(re) <- c("n_obj","type","mean_f","tuning")
  re$tuning[re$mean_f=="NaN" & re$tuning=="noData"] <- "missing"
  if(any(re$tuning=="noData")){
    warning("Some combinations have no values left after tuning, using original data")
  }
  print(re)
  #ls <-list(df,re)
  #return(ls)
  return(df)
  } else {
    cat("using att 1 and att 2",sep="\n")
  # correct by two attributes
  # get uniques for each pos
  u1 <- unique(df[,pos1])
  u2 <- unique(df[,pos2])

  # outer loop over pos2
  for (j in 1: length(u2)) {

    # inner loop for pos1
    for (i in 1:length(u1)) {
      ck <-rbind(ck,c(length(dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(dft$CTR[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]] / dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),digits = 4)))
      # get results to check for missing results
      if(ck[nrow(ck),3]=="NaN"){
        #cat(paste0("Class '",re[nrow(re),2],"' has no values left after tuning, using original data"),sep="\n")
        df$PAM_corr[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] * mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]])
        # store in df
        re <-rbind(re,c(length(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),digits = 4),"noData"))
        } else {
        df$PAM_corr[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]     * mean(dft$CTR[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]] / dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]])
        re <-rbind(re,c(length(dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(dft$CTR[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]] / dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),digits = 4),"fine"))

        }# end fork for checking
        }# end i loop
  }# end j loop

  # return
  #cat("Platzhalter2",sep="\n")
  colnames(re) <- c("n_obj","type","mean_f","tuning")
  re$tuning[re$mean_f=="NaN" & re$tuning=="noData"] <- "missing"
  if(any(re$tuning=="noData")){
    warning("Some combinations have no values left after tuning, using original data")
  }
  print(re)
  #ls <-list(df=df,re=re)
  #return(ls)
  return(df)
  }


  ##############################################################################

} # end of function

