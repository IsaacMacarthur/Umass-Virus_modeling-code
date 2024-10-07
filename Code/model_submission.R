# script to output the parquet file for submission, run model_fitting and data cleaning before this file.
submission_df <- prediction_sampler(stan = model, given_date  = target_date)
file_name <- paste(paste(target_date, "UMass", "HMLR", sep = "-"), ".parquet") 
write_parquet(submission_df,file_name )
