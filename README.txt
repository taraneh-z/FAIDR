FAIDR ("Feature Analysis of Intrinsically Disordered Regions")

Implementation of a probabilistic model to classify IDR function based on bulk molecular features. Implemented by Alan M. Moses (alan.moses@utoronto.ca) with contribution from Taraneh Zarin (taraneh.zarin@mail.utoronto.ca). 

Instructions:
1. Install R
2. In a directory, you should have the following files (provided):  
	i. "evolsig_unconstrained_10orthol_50perc_non_na.txt" (this is a table with predicted IDRs in the yeast proteome in each row, with its corresponding vector of z-scores representing the evolutionary signature for 164 molecular features [columns] (described in Zarin et al., eLife, 2019). 
	ii. "targets_unconstrained_evolsig.txt" (this is a table with corresponding rows of IDRs in the yeast proteome (as above, 2-1), but with a 1 value for each column if the protein is annotated for the corresponding function or phenotype [columns], and a 0 if it is not. 
	iii. FAIDR_script.R (this is the script that runs the model and outputs the results)
3. Set your directory to that containing the above 3 files. You can run
the script from the commandline with Rscript ./FAIDR_script.R, or from within R, you can use the command source("./FAIDR_script.R"). Note that this script depends on the R package "glmnet". The script should automatically install glmnet if you do not already  have it in your library.
4. The script will give you the following output tables: 
	i. "PostIDR_df.csv" (Posterior probability that a given IDR [rows] is associated with a given function or phenotype [columns])
	ii. "Tstats_df.csv" (T-statistic representing L1 regularized coefficients of the model; Molecular features [rows] by function or phenotype [columns]).
	iii. "Pred_df.csv" (Prediction that a gene (row) is associated with a function (column)
