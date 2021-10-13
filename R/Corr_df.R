#' Corr_df - Automated correction of PAM using mean for specified attributes.
#' @description uses the mean for specified subgroups to correct the PAM values.
#' @param df data.frame - containing Control and PAM values with colnames "CTR" and "PAM"
#' @param att1 character - Column name of data frame which should be used to correct by mean. If missing, will
#' use the global mean to correct PAM.
#' @param att2 character - Additional column name of data frame which should be used to correct by mean.
#' Only necessary if the correction should use subgroups for the data in 'pos1'
#' @param modef - boolean - if TRUE uses factor correction, if FALSE uses absolute values. Default=TRUE.
#' @return returns a data.frame with an additional column "PAM_corr" with the corrected PAM values
#' and prints the respective mean values for each class or combination of classes.
#' @note If used without any attributes will unse the global mean over all values to correct PAM. IF 'att1' is set
#' will calculate the mean for each group in 'att1' to correct PAM. IF 'att2' is given will use the mean for each group
#' with the combination for each attribute.
#' @author Andreas Schönberg
#' @export Corr_df
#' @aliases Corr_df
#' @examples
#' # load data
#' dat <- read.csv(system.file("extdata","exp_PAM.csv",package = "PAMcorrection"))
#' head(dat)
#'
#' # correct df by means for "generation"
#' corrected <-Corr_df(dat,att1 = "generation")
#'
#' # use two classes (means for each generation depending on dilect)
#' corrected <-Corr_df(dat,att1 = "generation",att2 = "dilect")
#'
#' # use absolute values
#' corrected <-Corr_df(dat,att1 = "generation",modef = F)
#' corrected <-Corr_df(dat,att1 = "generation",att2 = "dilect",modef = F)



Corr_df <-function(df,att1=NULL,att2=NULL,modef=T){

  # check input
  if(is.null(att1) & is.null(att2)==F){
    stop("missing 'att1' but given 'att2'")
  }

  # no attribute is given, calculate global mean
  if(is.null(att1)){
    if(modef==T){
      df$PAM_corr <-df$PAM * mean(df$CTR / df$PAM)
      cat(paste0("global mean (factor) "),mean(df$CTR / df$PAM))
      return(df)
    }
    if(modef==F){
      df$PAM_corr <-df$PAM * mean(df$PAM - df$PAM)
      cat(paste0("global mean (factor) "),mean(df$PAM - df$PAM))
      return(df)
    }

  }

  # get col position of attributes
  if(is.null(att2)){
  pos1 <- which(colnames(df)==att1)
  } else {
    pos1 <- which(colnames(df)==att1)
    pos2 <- which(colnames(df)==att2)
  }


  # init dataframe to store results
  re <-data.frame()
  # for factor correction
  if(modef==T){
  if(is.null(att2)){
  # correct by one attribute
  cat("using att 1 only",sep="\n")
  # get uniques for each pos
  u1 <- unique(df[,pos1])
  # inner loop for pos1
  for (i in 1:length(u1)) {
    #print(paste0(u1[i]," correction value ",mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ])))
    re <- rbind(re,c(length(df$PAM[ df[,pos1]==u1[i] ]),u1[i],round(mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ]),digits = 4)))
    df$PAM_corr[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ] * mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ])

     }# end i loop
  colnames(re) <- c("n_obj","type","mean_f")
  print(re)
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

      #print(paste0(u1[i]," ",u2[j]," correction value ",mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]])))
      re <-rbind(re,c(length(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),digits = 4)))

      df$PAM_corr[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] * mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]])
    }# end i loop
  }# end j loop
  colnames(re) <- c("n_obj","type","mean_f")
  print(re)
  return(df)
  }
  }# end if T
  if(modef==F){
    if(is.null(att2)){
    # correct by one attribute
    cat("using att 1 only",sep="\n")
    # get uniques for each pos
    u1 <- unique(df[,pos1])
    # inner loop for pos1
    for (i in 1:length(u1)) {
      #print(paste0(u1[i]," correction value ",mean(df$PAM[ df[,pos1]==u1[i] ] - df$CTR[ df[,pos1]==u1[i] ])))
      re <- rbind(re,c(length(df$PAM[ df[,pos1]==u1[i] ]),u1[i],round(mean(df$PAM[ df[,pos1]==u1[i] ] - df$CTR[ df[,pos1]==u1[i] ]),digits = 4)))
      df$PAM_corr[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ] - mean(df$PAM[ df[,pos1]==u1[i] ] - df$CTR[ df[,pos1]==u1[i] ])

    }# end i loop
    colnames(re) <- c("n_obj","type","mean_abs")
    print(re)
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

        #print(paste0(u1[i]," ",u2[j]," correction value ",mean(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]])))
        re <-rbind(re,c(length(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),digits = 4)))

        df$PAM_corr[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - mean(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]])
      }# end i loop
    }# end j loop
    colnames(re) <- c("n_obj","type","mean_abs")
    print(re)
    return(df)
  }
  }# end if F
} # end of function

