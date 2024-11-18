# run data_cleaning before this file
model <- stan_maker_splines(data = data_used, target_date = target_date, num_seq =  1, stan_file = "./Code/Heir_MLR_Splines_v2_ncp.stan", warmup = 7000, interations = 15000) # I don't know exactly what we are doing 
#about states with a low number of sequences, I set the threshold at 20 seq in the last 60 to be included, but even this 
#does not get all the states.
