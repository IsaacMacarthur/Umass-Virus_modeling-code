# run data_cleaning before this file
model <- stan_maker(data = data_used, target_date = target_date, num_seq =  1) # I don't know exactly what we are doing 
#about states with a low number of sequences, I set the threshold at 20 seq in the last 60 to be included, but even this 
#does not get all the states.
