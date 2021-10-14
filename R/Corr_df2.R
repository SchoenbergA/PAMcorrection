#' Corr_df2 - Automated correction of PAM using tuned mean for specified attributes.
#' @description uses the mean based on the data clipped by 'Whisker' and 'Quantil' for specified subgroups to correct the PAM values.
#' @param df data.frame - containing Control and PAM values with colnames "CTR" and "PAM"
#' @param att1 character - Column name of data frame which should be used to correct by mean. If missing, will
#' use the global mean to correct PAM.
#' @param att2 character - Additional column name of data frame which should be used to correct by mean.
#' Only necessary if the correction should use subgroups for the data in 'pos1'
#' @param tuning - character - method for clipping the data used to calculate the mean. Either "whisker" or "quantil". If 'NULL' will perform no tuning. Default= NULL.
#' @param modef - boolean - if TRUE uses factor correction, if FALSE uses absolute values. Default=TRUE.
#' @return returns a data.frame with an additional column "PAM_corr" with the corrected PAM values
#' and prints the respective mean values for each class or combination of classes. Uses 'tuned' data to calculate the mean values.
#' @note If used without any attributes will use the global mean over all values to correct PAM. IF 'att1' is set
#' will calculate the mean for each group in 'att1' to correct PAM. IF 'att2' is given will use the mean for each group
#' with the combination for each attribute.
#' @details If a specific class or combination of attributes has no values left after tuning the mean will be calculated based on the original data. This will lead to a warning.
#' The class or combination is marked as 'keep' in the result table.
#' @author Andreas Schönberg
#' @export Corr_df2
#' @aliases Corr_df2
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


Corr_df2 <-function(df,att1=NULL,att2=NULL,tuning=NULL,modef=T){

  # check input
  if(is.null(att1) & is.null(att2)==F){
    stop("missing 'att1' but given 'att2'")
  }

  ##############################################################################

  if (is.null(tuning)){
    # get copy with no changes
    dft <-df
    cat("No tuning - using original input data",sep="\n")
    warning("No tuning was used - but result table will show 'tuned' for each row. This is NO error.")
  } else if(tuning=="whisker"){
    # copy df
    dft <- df
    # get absolut difference
    dft$diffa <- dft$PAM-dft$CTR
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
    dft$diffa <- dft$PAM-dft$CTR
    # get boxplot
    bp <-boxplot(dft$diffa)
    # clip values by whisker
    boxplot(dft$diffa)
    #print(dft$diffa<bp$stats[2,1] | dft$diffa>bp$stats[4,1])
    dft <-dft[!(dft$diffa<bp$stats[2,1] | dft$diffa>bp$stats[4,1]),] # clip out of quantil
    cat(paste0("Tuning data by ",tuning),sep="\n")
    cat(paste0("clipping ",nrow(df)-nrow(dft), " entries",sep="\n"))
  }


  ##############################################################################


  # 1. Global mean #############################################################
  # if no attribute is given, calculate global mean
  if(is.null(att1)){
    if(modef==T){
      df$PAM_corr <-df$PAM * mean(dft$CTR / dft$PAM)
      cat(paste0("global mean (factor) "),mean(dft$CTR / dft$PAM))
      return(df)
    }
    if(modef==F){
      df$PAM_corr <-df$PAM * mean(dft$PAM - dft$CTR)
      cat(paste0("global mean (absolut) "),mean(dft$PAM - dft$CTR))
      return(df)
    }
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


  # 2. Mean factor #############################################################

  # for factor correction
  if(modef==T){
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
      cat(paste0("Class '",re[nrow(ck),2],"' has no values left after tuning, using original data"),sep="\n")
      warning(paste0("Class '",re[nrow(ck),2],"' has no values left after tuning, using original data"))
      # use original data
      df$PAM_corr[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ] * mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ])
      # store in df
      re <- rbind(re,c(length(df$PAM[ df[,pos1]==u1[i] ]),u1[i],round(mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ]),digits = 4),"keep"))
    } else {
    df$PAM_corr[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ]    * mean(dft$CTR[ dft[,pos1]==u1[i] ] / dft$PAM[ dft[,pos1]==u1[i] ])
    re <- rbind(re,c(length(dft$PAM[ dft[,pos1]==u1[i] ]),u1[i],round(mean(dft$CTR[ dft[,pos1]==u1[i] ] / dft$PAM[ dft[,pos1]==u1[i] ]),digits = 4),"tuned"))
    } # end fork for checking
     }# end i loop

  # return
  cat("Platzhalter1",sep="\n")
  colnames(re) <- c("n_obj","type","mean_f","tuning")
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
      ck <-rbind(ck,c(length(dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(dft$CTR[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]] / dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),digits = 4)))
      # get results to check for missing results
      if(ck[nrow(ck),3]=="NaN"){
        cat(paste0("Class '",re[nrow(re),2],"' has no values left after tuning, using original data"),sep="\n")
        warning(paste0("Class '",re[nrow(re),2],"' has no values left after tuning, using original data"))
        df$PAM_corr[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] * mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]])
        # store in df
        re <-rbind(re,c(length(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),digits = 4),"keep"))
        } else {
        df$PAM_corr[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]     * mean(dft$CTR[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]] / dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]])
        re <-rbind(re,c(length(dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(dft$CTR[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]] / dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),digits = 4),"tuned"))

        }# end fork for checking
        }# end i loop
  }# end j loop

  # return
  cat("Platzhalter2",sep="\n")
  colnames(re) <- c("n_obj","type","mean_f","tuning")
  print(re)
  return(df)
  }
  }# end if modef==T

  ##############################################################################


  # 3. Mean absolute #############################################################
  if(modef==F){
    if(is.null(att2)){
    # correct by one attribute
    cat("using att 1 only",sep="\n")
    # get uniques for each pos
    u1 <- unique(df[,pos1])
    # inner loop for pos1
    for (i in 1:length(u1)) {
      # get results to check for missing results
      ck <- rbind(ck,c(length(dft$PAM[ dft[,pos1]==u1[i] ]),u1[i],round(mean(dft$PAM[ dft[,pos1]==u1[i] ] - dft$CTR[ dft[,pos1]==u1[i] ]),digits = 4)))
      if(ck[nrow(ck),3]=="NaN"){
      cat(paste0("Class '",re[nrow(ck),2],"' has no values left after tuning, using original data"),sep="\n")
      warning(paste0("Class '",re[nrow(ck),2],"' has no values left after tuning, using original data"))
      # use original data
      df$PAM_corr[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ] - mean(df$PAM[ df[,pos1]==u1[i] ] - df$CTR[ df[,pos1]==u1[i] ])
      # store in df
      re <- rbind(re,c(length(df$PAM[ df[,pos1]==u1[i] ]),u1[i],round(mean(df$PAM[ df[,pos1]==u1[i] ] - df$CTR[ df[,pos1]==u1[i] ]),digits = 4),"keep"))
      } else {
      df$PAM_corr[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ] -    mean(dft$PAM[ dft[,pos1]==u1[i] ] - dft$CTR[ dft[,pos1]==u1[i] ])
      re <- rbind(re,c(length(dft$PAM[ dft[,pos1]==u1[i] ]),u1[i],round(mean(dft$PAM[ dft[,pos1]==u1[i] ] - dft$CTR[ dft[,pos1]==u1[i] ]),digits = 4),"tuned"))
      } # end fork for checking
    }# end i loop

    # return
    cat("Platzhalter1",sep="\n")
    colnames(re) <- c("n_obj","type","mean_abs","tuning")
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
        ck <-rbind(ck,c(length(dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]] - dft$CTR[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),digits = 4)))
        # get results to check for missing results
        if(ck[nrow(ck),3]=="NaN"){
          cat(paste0("Class '",re[nrow(re),2],"' has no values left after tuning, using original data"),sep="\n")
          warning(paste0("Class '",re[nrow(re),2],"' has no values left after tuning, using original data"))
          df$PAM_corr[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] * mean(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]])
          # store in df
          re <-rbind(re,c(length(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),digits = 4),"keep"))
        } else {
          df$PAM_corr[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]     * mean(dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]] - dft$CTR[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]])
          re <-rbind(re,c(length(dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(dft$PAM[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]] - dft$CTR[dft[,pos1]==u1[i] & dft[,pos2]==u2[j]]),digits = 4),"tuned"))

        }# end fork for checking
      }# end i loop
    }# end j loop

    # return
    cat("Platzhalter1",sep="\n")
    colnames(re) <- c("n_obj","type","mean_abs","tuning")
    print(re)
    return(df)
  }
  }# end if F
} # end of function

