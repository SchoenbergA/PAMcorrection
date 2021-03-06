% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TuningCorr_df.R
\name{TuningCorr_df}
\alias{TuningCorr_df}
\title{TuningCorr_df - Automated correction of PAM using tuned mean for specified attributes.}
\usage{
TuningCorr_df(df, att1 = NULL, att2 = NULL, tuning = NULL, threshold = NULL)
}
\arguments{
\item{df}{data.frame - containing Control and PAM values with colnames "CTR" and "PAM"}

\item{att1}{character - Column name of data frame which should be used to correct by mean. If missing, will
use the global mean to correct PAM.}

\item{att2}{character - Additional column name of data frame which should be used to correct by mean.
Only necessary if the correction should use subgroups for the data in 'pos1'}

\item{tuning}{character - method for clipping the data used to calculate the mean. Methods can be "whisker", "quantil" or "threshold". If 'NULL' will perform no tuning. Default= NULL.}

\item{threshold}{numeric - specifies the threshold to clip values for tuning. Only needed if 'tuning' == "threshold".}
}
\value{
returns a data.frame with an additional column "PAM_corr" with the corrected PAM values
and prints the respective mean values for each class or combination of classes. Uses 'tuned' data to calculate the mean values.
}
\description{
uses the mean based on the data clipped by 'Whisker' and 'Quantil' for specified subgroups to correct the PAM values.
}
\details{
\itemize{
\item tuning - available modes are "whisker" (which clips all values </> lower / upper whisker threshold in boxplot) and "quantil" (which clips all values </> 25 and 75 quantils).
"threshold" will clip all values > a specified threshold.
\item If a specific class or combination of attributes has no values left after tuning the mean will be calculated based on the original data. This will lead to a warning.
The class or combination is marked as 'noData' in the result table.
\item If a class or combination has a 'NaN' the class or combination is not existing in the dataframe. This is possible when using att1 and att2 an if not every att1 has all att2 classes.
}
}
\note{
If used without any attributes will use the global mean over all values to correct PAM. IF 'att1' is set
will calculate the mean for each group in 'att1' to correct PAM. IF 'att2' is given will use the mean for each group
with the combination for each attribute.
}
\examples{
# load data
dat <- read.csv(system.file("extdata","exp_PAM.csv",package = "PAMcorrection"))
head(dat)

# how the tuning works ("whisker" and "quantil")
bp <-boxplot(dat$PAM-dat$CTR)
bp$stats[1,1] # lower whisker
bp$stats[2,1] # lower quantil
bp$stats[3,1] # median
bp$stats[4,1] # upper quantil
bp$stats[5,1] # upper whisker

# correct df by global means
corrected <-TuningCorr_df(dat,tuning = "whisker")
corrected <-TuningCorr_df(dat,tuning = "quantil")

# use one class (no tuning)
corrected <-TuningCorr_df(dat,att1 = "dilect")
# use one class with threshold tuning
corrected <-TuningCorr_df(dat,att1 = "dilect",tuning="threshold",threshold=0.8)
# use two classes (means for each generation depending on dilect) tuned by 'whisker'
corrected <-TuningCorr_df(dat,att1 = "generation",att2 = "dilect",tuning = "whisker")
}
\author{
Andreas Schönberg
}
