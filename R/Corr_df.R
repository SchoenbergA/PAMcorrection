#' Corr_df - Automated correction of PAM using mean for specifiedy classes.
#' @description visualizes the correction of PAM values by correction value.
#' @param df data.frame - containing Control and PAM values with colnames "CTR" and "PAM"
#' @param pos1 numeric - Column of data frame which should be used to correct by mean.
#' @param pos2 numeric - Additional column of data frame which should be used to correct by mean.
#' Only necessary if the correction should use subgroups for the data in 'pos1'
#' @param modef - bolean - if TRUE uses factor correction, if FALSE uses absolute values. Default=TRUE.
#' @return returns a data.frame with the corrected PAM values and prints the respective mean vlaues for each class or combination of classes.
#' @author Andreas Schönberg
#' @export Corr_df
#' @aliases Corr_df
#' @examples
#' # load data
#' dat <- read.csv(system.file("extdata","exp_PAM.csv",package = "PAMcorrection"))
#' head(dat)
#'
#' # correct df by means for "generation"
#' which(colnames(dat)=="generation")
#' corrected <-Corr_df(dat,pos1 = 3)
#'
#' # use two classes (means for each generation depending on dilect)
#' which(colnames(dat)=="generation")
#' which(colnames(dat)=="dilect")
#' corrected <-Corr_df(dat,pos1 = 3,pos2 = 2)
#'
#' # use absolute values
#' corrected <-Corr_df(dat,pos1 = 3,modef = F)
#' corrected <-Corr_df(dat,pos1 = 3,pos2 = 2,modef = F)



Corr_df <-function(df,pos1,pos2=NULL,modef=T){
  # init dataframe to store results
  re <-data.frame()
  # for factor correction
  if(modef==T){
  if(is.null(pos2)){
  # correct by one attribute
  cat("using pos 1 only",sep="\n")
  # get uniques for each pos
  u1 <- unique(df[,pos1])
  # inner loop for pos1
  for (i in 1:length(u1)) {
    print(paste0(u1[i]," correction value ",mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ])))
    re <- rbind(re,c(length(df$PAM[ df[,pos1]==u1[i] ]),u1[i],round(mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ]),digits = 4)))
    df$PAM[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ] * mean(df$CTR[ df[,pos1]==u1[i] ] / df$PAM[ df[,pos1]==u1[i] ])

     }# end i loop
  colnames(re) <- c("n_obj","type","mean")
  print(re)
  return(df)
  } else {
    cat("using pos 1 and pos2",sep="\n")
  # correct by two attributes
  # get uniques for each pos
  u1 <- unique(df[,pos1])
  u2 <- unique(df[,pos2])

  # outer loop over pos2
  for (j in 1: length(u2)) {

    # inner loop for pos1
    for (i in 1:length(u1)) {

      print(paste0(u1[i]," ",u2[j]," correction value ",mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]])))
      re <-rbind(re,c(length(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),digits = 4)))

      df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] * mean(df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]] / df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]])
    }# end i loop
  }# end j loop
  colnames(re) <- c("n_obj","type","mean")
  print(re)
  return(df)
  }
  }# end if T
  if(modef==F){
    if(is.null(pos2)){
    # correct by one attribute
    cat("using pos 1 only",sep="\n")
    # get uniques for each pos
    u1 <- unique(df[,pos1])
    # inner loop for pos1
    for (i in 1:length(u1)) {
      print(paste0(u1[i]," correction value ",mean(df$PAM[ df[,pos1]==u1[i] ] - df$CTR[ df[,pos1]==u1[i] ])))
      re <- rbind(re,c(length(df$PAM[ df[,pos1]==u1[i] ]),u1[i],round(mean(df$PAM[ df[,pos1]==u1[i] ] - df$CTR[ df[,pos1]==u1[i] ]),digits = 4)))
      df$PAM[ df[,pos1]==u1[i] ] <-df$PAM[ df[,pos1]==u1[i] ] - mean(df$PAM[ df[,pos1]==u1[i] ] - df$CTR[ df[,pos1]==u1[i] ])

    }# end i loop
    colnames(re) <- c("n_obj","type","mean")
    print(re)
    return(df)
  } else {
    cat("using pos 1 and pos2",sep="\n")
    # correct by two attributes
    # get uniques for each pos
    u1 <- unique(df[,pos1])
    u2 <- unique(df[,pos2])

    # outer loop over pos2
    for (j in 1: length(u2)) {

      # inner loop for pos1
      for (i in 1:length(u1)) {

        print(paste0(u1[i]," ",u2[j]," correction value ",mean(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]])))
        re <-rbind(re,c(length(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),paste0(u1[i]," ",u2[j]),round(mean(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]]),digits = 4)))

        df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] <-df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - mean(df$PAM[df[,pos1]==u1[i] & df[,pos2]==u2[j]] - df$CTR[df[,pos1]==u1[i] & df[,pos2]==u2[j]])
      }# end i loop
    }# end j loop
    colnames(re) <- c("n_obj","type","mean")
    print(re)
    return(df)
  }
  }# end if F
} # end of function

