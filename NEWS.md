# PAMcorrection 0.0.0.98
Develop version Major Fixes

* TuningCorr_df - now uses abs diffenrece for boxplot
                - delete modus for mean by absolute diffenrece
* plotPAMcorr - add option to print green difference as real differnece (including negative values)
- ylim can now be adjusted in both directions

# PAMcorrection 0.0.0.98
* Bugfixes and fixes in documentation

# PAMcorrection 0.0.0.98
Develop version

* add new function 'TuningCorr_df' uses tresholds by whisker or quantil to crop data for calculating the means and uses them on the original data.
* some fixes for examples

# PAMcorrection 0.0.0.97
Develop version

* change input for 'Corr_df': from position to colname.
* 'Corr_df' now saves result in new column named 'PAM_corr' (which is used in 'plotPAMcorrection')
* add new function ' plotPAMcorrection'. Alternative version of 'estPAMcorr2' without in-function correction and advanced visualisation.
* improve printed results and add sigh for either factor or absolute values

# PAMcorrection 0.0.0.96
Develop version

* add stats for 'Corr_df'. Now prints amount of class objects and the mean value used to correct in a dataframe (in function printed not output)

# PAMcorrection 0.0.0.95
Develop version

* add absolute value mode to 'Corr_df'

# PAMcorrection 0.0.0.94
Develop version

* add Corr_df - function to calculate means by selected class or class and subclass
* add exmaples for functions

# PAMcorrection 0.0.0.93

* add estPAMcorr2 - now can sort data and has simplified input data but need correct format in df (cannot select pam and ctr is alternative names)

# PAMcorrection 0.0.0.91

Develop version

* add estPAMcorr - function for visualization of PAM and estimating correction values
* add example data

# PAMcorrection 0.0.0.9000

Initial version

* add Package script
* add NEWS

