# master function list 
library(tidyverse)
library(ggplot2)
library(vcdExtra)
library(rstan)
library(nnet)
library(zoo)
library("scoringRules")
library(dplyr)
loc_finder <- function(data, num_seq = 90, target_date = Sys.Date()){
  # takes in a virus dataset and returns all locations with at least num_seq cases in the last 60 days
  target_lo <- c()
  i = 1
  for(place in unique(data$location)){ # finding locations with over num_seq in last 60 days
    Virus_counts_state <- filter(data, date >= as.Date(target_date) - 60, location == place)
    # we are modeling from the first of June
    if(sum(Virus_counts_state$sequences) >= num_seq ){
      target_lo[i] <- place
      i = i + 1
    }
  }
  return(target_lo)
}
stan_maker <- function(data, num_seq = 90, target_date = Sys.Date(), num_days = 150, target_loc = NULL, interations = 3000, warmup = 1000, stan_file = "Weighted_Heir_MLR.stan" ){
  # data is a dataframe containing the columns sequences, location and date. 
  # num_seq is the number of seq to use as cutoff point, a numeric integer 
  # target date is the last day of the data to use, needs to be a date object
  # num_days is the number of days before the target date to use, a numeric interger
  # target_loc allows you to choose what locations to use, a string vector
  # returns a list containing the stan object, the number of locations, the number of clades, and the locations
  if( is.null(target_loc)){
    target_lo <- loc_finder(data = data, num_seq = num_seq, target_date = target_date)
  } else{
    target_lo = target_loc
  }
  data_case <- filter(data, location %in% target_lo, date >= as.Date(target_date) - num_days, date <= as.Date(target_date)) # dating from April 1st
  data_case$mlr <- c(1:length(data_case$clade)) # need numeric levels for the mlr
  j = 1
  for( k in unique(data_case$clade)){ # giving each clade a numeric level
    data_case$mlr <- ifelse(data_case$clade == k,j,data_case$mlr )
    j = j + 1
  }
  data_case$ll <- c(1:length(data_case$clade)) # need numeric levels for the locations
  j = 1
  for( k in unique(data_case$location)){
    data_case$ll <- ifelse(data_case$location == k,j,data_case$ll )
    j = j + 1
  }
  data_case$days <- as.numeric(as_date(data_case$date)) - as.numeric(as_date(as.Date(target_date) - num_days)) # days from start of dataset
  L = length(unique(data_case$ll))
  K = length(unique(data_case$clade))
  mlr_data <- list(
    weights = data_case$sequences, 
    L = L, # number of locations
    ll = data_case$ll, # where each case was located
    y = data_case$mlr, # the clade of each case 
    x = data_case$days, # the days the cases happened (first day of dataset treated as 0)
    N = length(data_case$location), # how many cases we had
    K = K # number of different clades
  )
  mlr_fit <- stan( # fitting the model
    file = stan_file,  
    data = mlr_data,    
    chains = 1,             
    warmup = warmup,          
    iter = interations,            
    cores = 8,
    refresh = 500,
  )
  return(list(mlr_fit = mlr_fit, L = L, K = K, target_lo = target_lo ))
}
stan_maker_Mech <- function(data, num_seq = 90, target_date = Sys.Date(), num_days = 150, target_loc = NULL, inital = NULL, alt_file = NULL, interations = 3000, warmup = 1000){
  # data is a dataframe containing the columns sequences, location and date. 
  # num_seq is the number of seq to use as cutoff point, a numeric integer 
  # target date is the last day of the data to use, needs to be a date object
  # num_days is the number of days before the target date to use, a numeric interger
  # target_loc allows you to choose what locations to use, a string vector
  # returns a list containing the stan object, the number of locations, the number of clades, and the locations
  if( is.null(target_loc)){
    target_lo <- loc_finder(data = data, num_seq = num_seq, target_date = target_date)
  } else{
    target_lo = target_loc
  }
  s_v = c(rep(-1, length(unique(data$clade)) - 1))
  j = 1
  for(clades in unique(data$clade)){
    if(clades == "other"){
      next
    }
    s_v[j] <- num_days - as.numeric(target_date - min(filter(data, clade == clades)$date))
    j = j + 1
  }
  data_case <- filter(data, location %in% target_lo, date >= as.Date(target_date) - num_days, date <= as.Date(target_date)) # dating from start date
  data_case$mlr <- c(1:length(data_case$clade)) # need numeric levels for the mlr
  j = 1
  for( k in unique(data_case$clade)){ # giving each clade a numeric level
    data_case$mlr <- ifelse(data_case$clade == k,j,data_case$mlr )
    j = j + 1
  }
  data_case$ll <- c(1:length(data_case$clade)) # need numeric levels for the locations
  j = 1
  for( k in unique(data_case$location)){
    data_case$ll <- ifelse(data_case$location == k,j,data_case$ll )
    j = j + 1
  }
  data_case$days <- as.numeric(as_date(data_case$date)) - as.numeric(as_date(as.Date(target_date) - num_days)) + 1 # days from start of dataset
  L = length(unique(data_case$ll))
  V = length(unique(data_case$clade))
  if(!is.null(inital)){
    I = list(I = inital)
    mlr_data <- list(
      weights = data_case$sequences, 
      L = L, # number of locations
      ll = data_case$ll, # where each case was located
      y = data_case$mlr, # the clade of each case 
      t = data_case$days, # the days the cases happened (first day of dataset treated as 0)
      N = length(data_case$location), # how many cases we had
      V = V, # number of different clades
      s_v = s_v,
      I = I
    )
    mlr_fit <- stan( # fitting the model
      file = "Mech_intial_v2.stan",  
      data = mlr_data,
      init = "random",
      init_r = 0.004,
      chains = 1,             
      warmup = warmup,          
      iter = interations,            
      cores = 8,
      refresh = 150,
    )
    return(list(mlr_fit = mlr_fit, L = L, V = V, target_lo = target_lo, start_times = s_v, I = I ))
  } else if(!is.null(alt_file)){
    mlr_data <- list(
      weights = data_case$sequences, 
      L = L, # number of locations
      ll = data_case$ll, # where each case was located
      y = data_case$mlr, # the clade of each case 
      t = data_case$days, # the days the cases happened (first day of dataset treated as 0)
      N = length(data_case$location), # how many cases we had
      V = V, # number of different clades
      s_v = s_v
    )
    mlr_fit <- stan( # fitting the model
      file = alt_file,  
      data = mlr_data,
      control = list(max_treedepth = 10, adapt_delta = 0.9),
      init = "random",
      init_r = 0.004,
      chains = 1,
      warmup = warmup,          
      iter = interations,            
      cores = 8,
      refresh = 150,
    )
    return(list(mlr_fit = mlr_fit, L = L, V = V, target_lo = target_lo, start_times = s_v ))
  } else {
    mlr_data <- list(
      weights = data_case$sequences, 
      L = L, # number of locations
      ll = data_case$ll, # where each case was located
      y = data_case$mlr, # the clade of each case 
      t = data_case$days, # the days the cases happened (first day of dataset treated as 0)
      N = length(data_case$location), # how many cases we had
      V = V, # number of different clades
      s_v = s_v
    )
    mlr_fit <- stan( # fitting the model
      file = "Mech_model_constant_nosv.stan",  
      data = mlr_data,    
      chains = 4,
      control = list(max_treedepth = 10, adapt_delta = 0.9),
      init = "random",
      init_r = 0.004,
      warmup = warmup,          
      iter = interations,            
      cores = 8,
      refresh = 150,
    )
  }
}
combine_clades <- function(data, target_date = Sys.Date, num_days = 100, num_seq = 50, recombinant = T){
  # combines minor clades into other
  # data is the virus count data frame
  # target_date is the day you want to predict from
  # num_days is the amount of days before target date you want to consider 
  # num_seq is the cutoff for the minor clades
  # recombinant is do you want recombinamt in minor clade, even if it exceeds the cutoff, default True
  minor_clades <- c() # the clades to be combined together
  j = 1
  for( clades in unique(data$clade )){ # finding clades with less than 50 variants
    if(sum(filter(data, clade == clades, date >= target_date - num_days, date <= target_date)$sequences) < num_seq ){
      minor_clades[j] <- clades
      j = j + 1
    }
  }
  if(recombinant){
    minor_clades[j] <- "recombinant"
  }
  
  data$clade <- fct_collapse(data$clade, other = minor_clades)
  return(data$clade)
}

mlr_fitter <- function(data, locations, target_date, days_before = 150){
  # data is the virus counts
  # locations is a vector of locations you want fitted
  # target_date is the date to the model from
  # days before is how many days you want to use to model, default 150
  # returns a named list over mlr models one for each location 
  mlrs <- list()
  if(is.null(data$days)){
    data$days <- as.numeric(as_date(data$date)) - as.numeric(as_date(as.Date(target_date) - days_before))
  }
  for( locat in locations){
    new_data <- filter(data, location == locat, date <= target_date, date >= target_date - days_before) 
    model <- multinom(clade ~ days, data = new_data, weights = sequences, Hess = T )
    mlrs[[locat]] <- model
  }
  return(mlrs)
}
get_energy_scores <- function(stan, old_data, new_data, target_date, dates = c(120:160), nowcast_length = 30, num_draws = 2000, N = 1, mlr_basic = NULL, skipped = NULL, shifted = F){
  # stan takes a list returned by stan maker
  # old_new is the data the stan was fit to
  # new_data is the data you want to compute the energy scores over
  # target date is the date the model was fit on
  # dates gives the dates from beginning of model to fit
  # nowcast_length gives the length of the nowcast
  # num_draws is how many non-warmup draws the stan model has
  # N is the number of times each posterior sample is used to generate samples  
  # mlr_basic can take a mlr model to compute energy scores over
  # skipped can take a vector of string names of locations to skip
  # shifted means the Mech model with s_v not equal to the zero vector is being used
  # returns a list of energy scores over locations
  energy_scores <- list() # the list of energy scores by location
  if( is.null(stan$K)){
    K <- stan$V 
  } else{
    K <- stan$K
  }
  L <- stan$L
  clades <- levels(old_data$clade)
  clades <- clades[-1]
  clades[K] <- "other" # have to put the clades in the right order
  draws <- extract(stan$mlr_fit)
  predicted_mat <- matrix(nrow = K, ncol = 100*N)
  random_draws <- matrix(nrow = K - 1, ncol = 2 )
  if(!is.null(mlr_basic)){
    energy_scores_mlr <- list()
  }
  for(l in 1:L){
    if( !is.null(skipped)){
      if(any(stan$target_lo[l] == skipped)){
        print(stan$target_lo[l])
        next
      }
    }
    old_data_loc <- filter(old_data, location == stan$target_lo[l])
    new_data_loc <- filter(new_data, location == stan$target_lo[l])
    if(!is.null(mlr_basic)){
      model <- mlr_basic[[l]]
      energy_vector_mlr <- c()
    }
    energy_vector <- c() # the vector of energy scores one per day
    observed_data <- c(rep(0, length(clades))) # the data we observed for the days
    for(i in 1:length(dates)){
      if(sum(filter(old_data_loc, date == target_date - nowcast_length + i )$sequences) != 0 || sum(filter(new_data_loc, date == target_date - nowcast_length + i )$sequences) == 0){
        next 
      } else{
        for(j in 1:length(clades)){
          observed_data[j] <-sum(filter(new_data_loc, date == target_date - nowcast_length + i, clade == clades[j] )$sequences) 
        } # the matrix of X_i
        for( m in 1:100){ # the matrix of random draws
          for(q in 1:(K-1)){
            random_draws[q, ] <- c(draws$raw_alpha[ceiling(runif(1, min = 0, max = num_draws)),l,q], draws$raw_beta[ceiling(runif(1, min = 0, max = num_draws)),l,q] ) # getting the random draws 
          }
          if(!shifted){
            days <- exp(random_draws[, 1] + random_draws[, 2]*dates[i])/(sum(exp(random_draws[, 1] + random_draws[, 2]*dates[i]))+1) # softmaxing
          } else{
            days <- exp(random_draws[, 1] + random_draws[, 2]*(dates[i] - stan$start_times))
            for( you in 1:(K-1)){
              if(stan$start_times[you] >= i ){
                days[you] <- 0.001
              }
            }
            days <- days/(sum(days) + 1)
          }
          days[K] <- 1 - sum(days)
          predicted_mat[, (1 + (m-1)*(N)):(m*N)] <- rmultinom(N, sum(observed_data),days)
        }
        
        energy_vector[length(energy_vector) + 1] <- es_sample(observed_data, predicted_mat)
        if(!is.null(mlr_basic)){
          test_days <- data.frame(days = dates) # the days we are predicting on
          prediction <- predict(model, type = "probs", newdata = test_days)
          prediction <- t(prediction)
          predicted_mat <- rmultinom(100*N, sum(observed_data),c(prediction[2:K, i], prediction[1, i]))
          energy_vector_mlr[length(energy_vector_mlr) + 1] <- es_sample(observed_data, predicted_mat)
        }
      }
    }
    energy_scores[[stan$target_lo[l]]] <- energy_vector
    if(!is.null(mlr_basic)){
      energy_scores_mlr[[stan$target_lo[l]]] <- energy_vector_mlr
    }
  }
  if(is.null(mlr_basic)){
    return(energy_scores)
  } else{
    return(list(energy_scores = energy_scores, energy_scores_mlr = energy_scores_mlr))
  }
}

mlr_probs <- function(stan, num_days, shifted = F){
  # returns a list of probs from day 0 to the given day
  # takes a stan object, returned by stan_maker and the number of days probs wanted, if shifted is
  # true, then uses time since the variant was introduced instead of days from beginning of dataset, used for
  # Mech model non-zero s_v
  full_probs <- list()
  means <- get_posterior_mean(stan$mlr_fit, pars = c("raw_alpha", "raw_beta")) # the alpha and beta
  L = stan$L
  if(!is.null(stan$K)){
    K = stan$K
  } else{
    K = stan$V
  }
  for(l in 1:L){
    intercepts <- means[(1 + (K-1)*(l-1)):(K - 1 + (K-1)*(l-1))] # the alphas
    coef <- means[(1 + (K-1)*L + (K-1)*(l-1)):((K-1 + (K-1)*L + (K-1)*(l-1))) ] # the beta's
    probs <- matrix( nrow = K, ncol = num_days) # the matrix of probabilities for each day
    dates <- c(1:num_days)
    if(!(shifted)){
      for( i in dates){ # getting the probabilities 
        days <- exp(intercepts + coef*dates[i])/(sum(exp(intercepts + coef*dates[i]))+1)
        days[K] <- 1 - sum(days)
        probs[, i] <- days
      } 
    } else {
      for( i in dates){ # getting the probabilities 
        days <- exp(intercepts + coef*(dates[i] - stan$start_times))
        for( m in 1:(K-1)){
          if(stan$start_times[m] >= i ){
            days[m] <- 0.001
          }
        }
        days <- days/(sum(days) + 1)
        days[K] <- 1 - sum(days)
        probs[, i] <- days
      }
    }
    full_probs[[stan$target_lo[l]]] <- probs
  }
  return(full_probs)
}
CI_maker <- function(stan, num_days, CI_level = 0.9, shifted = F){
  # takes in a stan object returned by stan maker, the number of days of probablities wanted, and the CI level, with default of 
  # 0.9, returns a named list of CI upper and lower bounds for each location, set shifted to True, if using Mech model
  # with non-zero s_v
  L <- stan$L
  CIs <- list()
  if( is.null(stan$K)){
    K <- stan$V
  } else{
    K <- stan$K
  }
  for( l in 1:L){
    raw_alpha <- extract(stan$mlr_fit)$raw_alpha[, l, ]
    raw_beta <- extract(stan$mlr_fit)$raw_beta[, l, ]
    upper <- matrix(nrow = K, ncol = num_days)
    lower <- matrix(nrow = K, ncol = num_days)
    sample_probs <- matrix(nrow = K, ncol = length(raw_alpha[, 1]))
    for(i in 1:num_days){
      if(!shifted){
        sample_probs[1:(K-1),] <- t(exp( raw_alpha + raw_beta*i)/(rowSums(exp( raw_alpha + raw_beta*i)) + 1))
      } else{
        temp <- matrix( nrow = K-1, ncol = length(raw_alpha[, 1]))
        temp <- exp( t(raw_alpha) + t(raw_beta)*(i - stan$start_times))
        temp <- t(temp)
        for( j in 1:(K-1)){
          if(stan$start_times[j] >= i){
            temp[, j ] <- c(rep(0.001, length(temp[,j ])))
          }
        }
        temp <- temp/(rowSums(temp) + 1)
        sample_probs[1:(K-1), ] <- t(temp)
      }
      sample_probs[K, ] <- 1 - colSums(sample_probs[1:(K-1),])
      for( k in 1:K){
        upper[k, i] <- quantile(sample_probs[k, ], probs = CI_level)
        lower[k, i] <- quantile(sample_probs[k, ], probs = 1- CI_level)
      }
    }
    name <- 
      CIs[[paste(stan$target_lo[l], "upper")]] <- upper
    CIs[[paste(stan$target_lo[l], "lower")]] <- lower
  }
  return(CIs)
}
plot_data <- function(stan, data, colors = c("black","blue", "red", "green", "yellow"), target_date = Sys.Date(), num_days = 150, shifted = F, CI = F, CI_level = 0.9, other_probs = NULL){
  # stan is the object returned by stan maker or stan maker mech
  # data is the virus data sets which points are to be ploted
  # colors is the list of colors to be used for plotting
  # target date is the last date you want to be plotted
  # num_days is the amount of days to be plotted
  # if shifted is true that means this is an s_v model, default False
  # If CI is true, CI lines will be plotted default False
  # CI level gives the CI level, does nothing if CI = False
  # other_probs can take in a list of other probabilities to be plotted. The defualt is null
  L <- stan$L
  if(!(is.null(stan$K))){
    K <- stan$K
  } else {
    K <- stan$V
  }
  model_probs <- mlr_probs(stan = stan, num_days = num_days, shifted = shifted)
  if(CI){
    CI_probs <- CI_maker(stan = stan, num_days = num_days, CI_level = CI_level, shifted = shifted)
  }
  for(l in 1:L){ 
    dates <- c(1:num_days)
    clades <- levels(data$clade)
    clades <- clades[-1]
    clades[K] <- "other" # have to put the clades in the right order
    sample_probs <- matrix( nrow = K, ncol = num_days) # the matrix of observed probabilities
    weights <- c(rep(0, num_days)) # weights for sizes of sample points
    for( i in dates){ # getting the probabilities 
      total_seq <- filter(data, date == (target_date - num_days)  + i, location == stan$target_lo[l]  )
      weights[i] <- sum(total_seq$sequences) # counting the number of seq per day
      days <- (rep(0, K))
      if(sum(total_seq$sequences) == 0){
        sample_probs[, i] <- c(rep(0, K)) # if no seq on a day, report 0
      } else{
        for( p in 1:K){
          days[p] <- sum(filter(total_seq, clade == clades[p])$sequences)/sum(total_seq$sequences) # the sample probabilities per day
        }
        sample_probs[, i] <- days
      }
    }
    weights <- 4*weights/max(weights)
    main = paste(stan$target_lo[l], "Virus Probabilities(observed and predicted)") # creating the plots
    plot(c(1:num_days), model_probs[[stan$target_lo[l]]][1, ], type = 'l', ylim = c(0,1), ylab = "Probability", xlab = paste("time from", target_date - num_days), col = colors[1], main = main) # Heir_MLR probs
    for( num in 1:K){
      points(c(1:num_days), sample_probs[num, ],col = colors[num], cex = weights ) # sample probs
      if(CI){
        lines(c(1:num_days), CI_probs[[paste(stan$target_lo[l], "upper")]][num, ], col = colors[num], lty = 3)
        lines(c(1:num_days), CI_probs[[paste(stan$target_lo[l], "lower")]][num, ], col = colors[num], lty = 3)
      }
      if(!(is.null(other_probs))){
        lines(c(1:num_days), other_probs[[stan$target_lo[l]]][num, ], col = colors[num], lty = 4)
      }
    }
    for( num in 2:K){
      lines(c(1:num_days), model_probs[[stan$target_lo[l]]][num, ], col = colors[num])
    }
  }
}
#Function to delete a percentage of observations between a start and end date, chat GPT created code
delete_percentage <- function(df, start_date, end_date, percentage) {
  # Ensure the percentage is between 0 and 1
  if (percentage < 0 | percentage > 1) {
    stop("Percentage must be between 0 and 1.")
  }
  
  # Convert start and end dates to Date class if they aren't already
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  # Filter data between start_date and end_date
  filtered_data <- df %>%
    filter(date >= start_date & date <= end_date)
  
  # Determine the number of rows to delete based on the percentage
  num_to_delete <- round(nrow(filtered_data) * percentage)
  
  # Randomly sample the rows to delete
  rows_to_delete <- sample(1:nrow(filtered_data), size = num_to_delete)
  
  # Remove the selected rows from the filtered data
  filtered_data_deleted <- filtered_data[-rows_to_delete, ]
  
  # Combine back the rows outside the date range with the remaining rows
  final_df <- df %>%
    filter(!(date >= start_date & date <= end_date)) %>%
    bind_rows(filtered_data_deleted)
  
  return(final_df)
}
