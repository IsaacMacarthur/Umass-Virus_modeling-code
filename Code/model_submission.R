# script to output the parquet file for submission, run model_fitting and data cleaning before this file.
set.seed(1)
submission_df <- prediction_sampler(stan = model,N=100, given_date  = target_date, splines = T)
file_name <- paste0(paste(target_date, "UMass", "HMLR", sep = "-"), ".parquet") 
write_parquet(submission_df,file_name )
