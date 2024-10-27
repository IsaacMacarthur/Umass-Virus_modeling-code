
   library(splines)
stan_maker_splines <- function(data, B = 3, num_seq = 90, target_date = Sys.Date(), num_days = 150, target_loc = NULL, interations = 3000, warmup = 1000, stan_file = "Heir_MLR_Splines.stan" ){
  # data is a dataframe containing the columns sequences, location and date.
  # B is the degree of the spline
  # num_seq is the number of seq to use as cutoff point, a numeric integer 
  # target date is the last day of the data to use, needs to be a date object
  # num_days is the number of days before the target date to use, a numeric interger
  # target_loc allows you to choose what locations to use, a string vector
  # returns a list containing the stan object, the number of locations, the number of clades, the locations, and the 
  # number of basis functions
  if( is.null(target_loc)){
    target_lo <- loc_finder(data = data, num_seq = num_seq, target_date = target_date)
  } else{
    target_lo = target_loc
  }
  data_case <- filter(data, location %in% target_lo, date >= as.Date(target_date) - num_days, date <= as.Date(target_date)) # dating from April 1st
  data_case$mlr <- c(1:length(data_case$clade)) # need numeric levels for the mlr
  j = 1
  for( k in unique(data$clade)){ # giving each clade a numeric level
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
  clades <- unique(data$clade)
  temp <- clades[1]
  clades <- clades[-1]
  clades[K] <- temp
  mlr_data <- list(
    weights = data_case$sequences, # the number of each trio
    L = L, # number of locations
    ll = data_case$ll, # where each case was located
    y = data_case$mlr, # the clade of each case 
    x = t(bs(data_case$days, degree =  B)), # the days the cases happened (first day of dataset treated as 0)
    N = length(data_case$location), # how many cases we had
    K = K, # number of different clades
    B = B
  )
  mlr_fit <- stan( # fitting the model
    file = stan_file,  
    data = mlr_data,    
    chains = 1,             
    warmup = warmup,          
    iter = interations,            
    cores = 8,
    refresh = 150,
  )
  return(list(mlr_fit = mlr_fit, L = L, K = K, target_lo = target_lo, B = B, clades = clades ))
}
mlr_probs_splines <- function(stan, num_days){
  # returns a list of probs from day 0 to the given day
  # takes a stan object from the spline model, returned by stan_maker and the number of days probs wanted,
  full_probs <- list()
  means <- extract(stan$mlr_fit, pars = c("raw_alpha", "raw_beta")) # the alpha and beta 
  L = stan$L
  B = stan$B
  if(!is.null(stan$K)){
    K = stan$K
  } else{
    K = stan$V
  }
  for(l in 1:L){
    intercepts <- means$raw_alpha[, l, ] # the alphas
    coef <- means$raw_beta[ , l , , ] # the beta's
    probs <- array( dim = c(K, num_days, length(intercepts[1,]))) # the matrix of probabilities for each day
    mean_probs <- matrix( nrow = K, ncol = num_days)
    dates <- bs(c(1:num_days), degree =  B)
    days <- c(rep(0, K))
    for(j in 1:length(intercepts[1, ])){
      for( i in 1:num_days){ # getting the probabilities 
        for(k in 1:K-1){
          days[k] <- exp(intercepts[j, k] + sum(coef[j,k,]*dates[i, ])) 
        }
        days[1:K-1] <- days[1:K-1]/(sum(days[1:K-1]) +1)
        days[K] <- 1 - sum(days[1:(K-1)])
        probs[  , i, j] <- days
      }
    }
    for(i in 1:num_days){
      for(j in 1:K){
        mean_probs[j, i] <- mean(probs[j, i, ])
      }
    }
    full_probs[[stan$target_lo[l]]] <- mean_probs
  }
  for(lo in stan$target_lo){
    row.names(full_probs[[lo]]) <- stan$clades
  }
  return(full_probs)
}
get_energy_scores_splines <- function(stan, old_data, new_data, target_date, dates = c(119:160), nowcast_length = 31, N = 1, num_draws = 2000, mlr_basic = NULL){
  # stan takes a list returned by stan maker, here you are using the spline model
  # old_new is the data the stan was fit to
  # new_data is the data you want to compute the energy scores over
  # target date is the date the model was fit on
  # dates gives the dates from beginning of model to fit
  # N is the number of times to sample for each posterior sample 
  # nowcast_length gives the length of the nowcast
  # num_draws is how many non-warmup draws the stan model has
  # mlr_basic can take a mlr model to compute energy scores over
  # returns a list of named energy scores over locations
  energy_scores <- list() # the list of energy scores by location
  if( is.null(stan$K)){
    K <- stan$V 
  } else{
    K <- stan$K
  }
  L <- stan$L
  B <- stan$B
  clades <- levels(old_data$clade)
  clades <- clades[-1]
  clades[K] <- "other" # have to put the clades in the right order
  draws <- extract(stan$mlr_fit)
  splines <- t(bs(c(1:max(dates)), degree =  B))
  predicted_mat <- matrix(nrow = K, ncol = 100*N) # the matrix of X_i
  random_draws <- matrix(nrow = K - 1, ncol = B + 1 ) # the matrix of random draws
  num_draws <- length(draws$raw_alpha[,1,1])
  if(!is.null(mlr_basic)){
    energy_scores_mlr <- list()
  }
  for(l in 1:L){
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
        } 
        for( m in 1:100){
          days <- c(rep(0, K))
          rn <- ceiling(runif(1, min = 0, max = num_draws))
          for(q in 1:(K-1)){
            random_draws[q, 1] <- draws$raw_alpha[rn,l,q]
            for( b in 1:B){
              random_draws[q, b + 1] <- draws$raw_beta[rn,l,q,b] # getting the random draws 
            }
            days[q] <- exp(random_draws[q, 1] + sum(random_draws[q, 2:(B+1)]*splines[, dates[i]]))
          }
          days[1:(K-1)] <- days[1:(K-1)]/(sum(days[1:(K-1)]) + 1) # softmaxing
          days[K] <- 1 - sum(days)
          predicted_mat[, (1 + (m-1)*N):(m*N)] <- rmultinom(N, sum(observed_data),days)
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
