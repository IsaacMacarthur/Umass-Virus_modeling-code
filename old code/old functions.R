# this is old functions that are probably not useful
Mech_energy_scores <- function(stan, old_data, new_data, target_date, dates = c(120:160), nowcast_length = 30, mlr_basic = NULL, skipped = NULL){
  # used for the linear Mech model
  energy_scores <- list() # the list of energy scores by location
  V <- stan$V
  L <- stan$L
  clades <- unique(old_data$clade)
  clades <- clades[-1]
  clades[V] <- "other" # have to put the clades in the right order
  draws <- extract(stan$mlr_fit)
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
        }
        predicted_mat <- matrix(nrow = V, ncol = 100) # the matrix of X_i
        for( m in 1:100){
          random_draws <- matrix(nrow = V - 1, ncol = 2 ) # the matrix of random draws
          for(q in 1:(V-1)){
            random_draws[q, ] <- c(draws$I[ceiling(runif(1, min = 0, max = 3000)),l,q], draws$r[ceiling(runif(1, min = 0, max = 3000)),l,q] ) # getting the random draws 
          }
          days <- c(rep(0, V-1))
          for(w in 1:(V-1)){ # calculating the probabilities
            days[w] <- random_draws[w,1]
            if( dates[i] <= stan$start_times[w]){ # skipping clades before their starting time
              days[w] = 0
            } else{
              for(s in 1:(dates[i] - stan$start_times[w] - 1)){
                days[w] = days[w]*(exp(s*random_draws[w,2]))
              }
              days[w] = days[w]*expm1(random_draws[w,2]*(dates[i] - stan$start_times[w])) 
            }
          }
          days = days/(sum(days) + 1)
          days[V] <- 1 - sum(days)
          predicted_mat[, m] <- rmultinom(1, sum(observed_data),days)
        }
        energy_vector[length(energy_vector) + 1] <- es_sample(observed_data, predicted_mat)
        if(!is.null(mlr_basic)){
          test_days <- data.frame(days = dates) # the days we are predicting on
          prediction <- predict(model, type = "probs", newdata = test_days)
          prediction <- t(prediction)
          predicted_mat <- rmultinom(100, sum(observed_data),c(prediction[2:V, i], prediction[1, i]))
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
stan_maker_shifted <- function(data, num_seq = 90, target_date = Sys.Date(), num_days = 150, target_loc = NULL){
  # data is a dataframe containing the columns sequences, location and date.
  # file_name is the name of the stan file you want to use.
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
  data_case$days <- as.numeric(as_date(data_case$date)) - as.numeric(as_date(as.Date(target_date) - num_days)) # days from start of dataset
  L = length(unique(data_case$ll))
  V = length(unique(data_case$clade))
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
    file = "shifted_time_model_v2.stan",  
    data = mlr_data,    
    chains = 1,
    control = list(max_treedepth = 10, adapt_delta = 0.95),
    init = "random",
    init_r = 0.004,
    warmup = 2000,          
    iter = 5000,            
    cores = 8,
    refresh = 500,
  )
  return(list(mlr_fit = mlr_fit, L = L, V = V, target_lo = target_lo, start_times = s_v ))
}
plotter <- function(stan, data, colors, target_date = Sys.Date(), days_before = 150, days_ahead = 10, plot_rolling = T, plot_data = T, mlr_basic = NULL){
  # stan is the object returned by stan_maker
  # data is the Virus_count data-frame
  # colors is a vector of color names of length K
  # target_date is the date you are predicting from
  # days before is how many days before the target date you are using as data
  # days_ahead is how many days you want to predict past the target date
  # returns a  named list containing plots a named list of plots, mAE_rolling_hier a matrix of MAE for the heir model for each location, and mAE_rolling_MLR a matrix of MAE for the MLR for each location, and MLR_model and named list of MLR
  #models
  num_days <- days_before + days_ahead # the number of days to plot
  L = stan$L
  K = stan$K
  plots <- list() # the plots
  models <- list() # the mlr_models for each location
  mAE_rolling_hier <- matrix(nrow = L, ncol = 2*days_ahead ) # the rolling averages
  mAE_rolling_MLR <- matrix(nrow = L, ncol = 2*days_ahead) 
  for(l in 1:L){ # one plot for each location
    means <- get_posterior_mean(stan$mlr_fit, pars = c("raw_alpha", "raw_beta")) # the alpha and beta
    intercepts <- means[(1 + (K-1)*(l-1)):(K - 1 + (K-1)*(l-1))] # the alphas
    coef <- means[(1 + (K-1)*L + (K-1)*(l-1)):((K-1 + (K-1)*L + (K-1)*(l-1))) ] # the beta's
    probs <- matrix( nrow = K, ncol = num_days) # the matrix of probabilities for each day
    dates <- c(1:num_days)
    clades <- unique(data$clade)
    clades <- clades[-1]
    clades[K] <- "other" # have to put the clades in the right order
    sample_probs <- matrix( nrow = K, ncol = num_days) # the matrix of observed probabilities
    weights <- c(rep(0, num_days)) # weights for sizes of sample points
    for( i in dates){ # getting the probabilities 
      days <- exp(intercepts + coef*dates[i])/(sum(exp(intercepts + coef*dates[i]))+1)
      days[K] <- 1 - sum(days)
      probs[, i] <- days
      total_seq <- filter(data, date == (target_date - days_before)  + i, location == stan$target_lo[l]  )
      weights[i] <- sum(total_seq$sequences) # counting the number of seq per day
      if(sum(total_seq$sequences) == 0){
        sample_probs[, i] <- c(rep(0, K)) # if no seq on a day, report 0
      } else{
        for( p in 1:K){
          days[p] <- sum(filter(total_seq, clade == clades[p])$sequences)/sum(total_seq$sequences) # the sample probabilities per day
        }
        sample_probs[, i] <- days
      }
    }
    rolldates <- matrix(nrow = K, ncol = 6) # the extra dates needed for the rolling average
    j = 1
    for( i in -5:0){
      total_seq <- filter(data, date == (target_date - days_before)  + i, location == stan$target_lo[l]  )
      if(sum(total_seq$sequences) == 0){
        rolldates[, j] <- c(rep(0, K))
      } else{
        for( p in 1:K){
          days[p] <- sum(filter(total_seq, clade == clades[p])$sequences)/sum(total_seq$sequences)
        }
      }
      rolldates[ , j] <- days
      j = j + 1
    }
    weights <- 4*weights/max(weights)
    rolling <- cbind(rolldates, sample_probs)
    main = paste(stan$target_lo[l], "Virus Probabilities(observed and predicted)") # creating the plots
    plot(c(1:num_days), probs[ 1, ], type = 'l', ylim = c(0,1), ylab = "Probability", xlab = paste("time from", target_date - days_before), col = colors[1], main = main) # Heir_MLR probs
    if( all(c(plot_data, plot_rolling))){
      for( num in 1:K){
        points(c(1:num_days), sample_probs[num, ],col = colors[num], cex = weights ) # sample probs, weighted by how many seq on that day, max size is normalized to 4
        lines(c(1:num_days),rollmean(rolling[num,], k = 7, align = "left"), col = colors[num], lty = "dashed" ) # seven day rolling averages
      }
    } else if(plot_data){
      for( num in 1:K){
        points(c(1:num_days), sample_probs[num, ],col = colors[num], cex = weights ) # sample probs
      }
    } else if( plot_rolling){
      for( num in 1:K){
        lines(c(1:num_days),rollmean(rolling[num,], k = 7, align = "left"), col = colors[num], lty = "dashed" )
      }
    }
    lines(c(days_before, days_before),c(0,1), col = "black" ) # the end of the generated data 
    lines(c(days_before -30, days_before - 30), c(0,1), col = "gray")
    for( num in 2:K){
      lines(c(1:num_days), probs[ num, ], col = colors[num])
    }
    word <- stan$target_lo[l]
    if(all(rowSums(sample_probs[, 1:days_before]) != 0)){ # if any of the locations don't have at least one seq of each clade we skip the mlr, because it won't work
      if(is.null(data$days)){
        data$days <- as.numeric(as_date(data$date)) - as.numeric(as_date(as.Date(target_date) - days_before))
      }
      if(is.null(mlr_basic)){
        new_data <- filter(data, location == stan$target_lo[l], date <= target_date, date >= target_date - days_before) 
        model <- multinom(clade ~ days, data = new_data, weights = sequences, Hess = T ) # mlr model 
      } else{
        model <- mlr_basic[[stan$target_lo[l]]]
      }
      
      test_days <- data.frame(days = c(1:num_days)) # the days we are predicting on
      prediction <- predict(model, type = "probs", newdata = test_days)
      lines(c(1:num_days),prediction[ , 1], col = colors[K], lty = "dotted" ) # mlr fits other
      f <- t(prediction)
      error <- abs(f[1, (days_before - days_ahead+1):(num_days)] - rollmean(rolling[K, ], k = 7, align = "left")[(days_before+1 - days_ahead):(num_days)] )
      models[[word]] <- model
      for( num in 1:(K-1)){
        lines(c(1:num_days),prediction[ , num + 1], col = colors[num], lty = "dotted" )
        error = error + abs(f[num + 1, (days_before - days_ahead+1):(num_days)] - rollmean(rolling[num, ], k = 7, align = "left")[(days_before+1 - days_ahead):(num_days)] ) # have to adjust for the shift in varients
      }
      mAE_rolling_MLR[l, ] <- (1/K)*error
    } else{
      mAE_rolling_MLR[l, ] <- rep(-1, 2*days_ahead ) # report -1
    }
    error2 <- 0
    for( num in 1:K){
      error2 = error2 + abs(probs[num,(days_before - days_ahead+1):(num_days) ] - rollmean(rolling[num, ], k = 7, align = "left")[(days_before+1 - days_ahead):(num_days)] ) # the error
    }
    mAE_rolling_hier[l, ] <- (1/K)*error2
    plotsave <- recordPlot()
    plots[[word]] <- plotsave # the plots
  }
  return(list( plots = plots, mAE_rolling_MLR = mAE_rolling_MLR, mAE_rolling_hier = mAE_rolling_hier, MLR_models = models))
}
logit_plotter <- function(stan, mlr, data, colors, start_date, number_days = 180, num_predicted = 30){
  # takes in the object for stan_maker a saved named list of MlR models, returned by plotter, a list of colors, a start date, and number_days the number of days that should be plotted, and the number of days predicted into the future  
  K <- stan$K
  L <- stan$L
  plots <- list()
  for(l in 1:stan$L){
    means <- get_posterior_mean(stan$mlr_fit, pars = c("raw_alpha", "raw_beta")) # the alpha and beta
    intercepts <- means[(1 + (K-1)*(l-1)):(K - 1 + (K-1)*(l-1))] # the alphas
    coef <- means[(1 + (K-1)*L + (K-1)*(l-1)):((K-1 + (K-1)*L + (K-1)*(l-1))) ] # the beta's
    num_days <- c(0:number_days)
    vals <- matrix(nrow = K-1, ncol = length(num_days))
    sample_log_prop <- matrix(nrow = K-1, ncol = length(num_days))
    clades <- unique(data$clade)
    clades <- clades[-1]
    clades[K] <- "other" # putting the clade names in the right order
    data_loc <- filter(data, location == stan$target_lo[l]) # pre-filter to make code run faster
    weights <- matrix(nrow = K-1, ncol = length(num_days)) # the weights for the logits
    for( i in 1:(K-1)){
      vals[i, ] <- intercepts[i] + coef[i]*num_days
      for(W in 0:num_days[length(num_days)]){
        clade_sum <- sum(filter(data_loc, clade == clades[i], date == start_date + W)$sequences) # sum of seq of target clade on target date
        other_sum <- sum(filter(data_loc, clade == "other", date == start_date + W)$sequences) # sum of seq of other on target day
        if(clade_sum != 0 & other_sum != 0){
          sample_log_prop[i, W + 1] <- log(clade_sum/other_sum)
          weights[i, W + 1] <- clade_sum + other_sum
        } else{
          sample_log_prop[i, W + 1 ] <- 70000
          weights[i, W + 1] <- 0
        }
        
      }
      weights[i, ] <- 4*weights[i, ]/max(weights[i, ]) 
    }
    
    main = paste(stan$target_lo[l], "Logits initial vs mean") # creating the plots
    plot(num_days, vals[ 1, ], type = 'l', ylim = c(min(vals), max(vals)), ylab = "Logit (with respect to other)", xlab = paste("time from", start_date), col = colors[1], main = main)
    points(sample_log_prop[1,], col = colors[1], cex = weights[1, ])
    lines(c(number_days - num_predicted, number_days - num_predicted),c(min(vals), max(vals))) # the date we are predicting from
    for(i in 2:(K-1)){
      lines(num_days, vals[ i, ], col = colors[i])
      points(sample_log_prop[i, ], col = colors[i], cex = weights[i, ])
    }
    if(!(is.null(mlr[[stan$target_lo[l]]]))){ # ploting the mlr if given
      mlr_coef <- coef(mlr[[stan$target_lo[l]]])
      for(i in 1:(K-1)){
        lines(num_days, (mlr_coef[i, 1] + mlr_coef[i, 2]*num_days), col = colors[i], lty = "dashed")
      }
    }
    word <- stan$target_lo[l]
    plotsave <- recordPlot()
    plots[[word]] <- plotsave
  }
  return(plots)
}
logit_clade_plotter <- function(stan, data, start_date, colors, number_days = 180){
  # takes a stan object returned by stan maker, a start date and an end date
  # returns a list of plots each containing all logits of a clade over all locations
  K <- stan$K
  L <- stan$L
  plots <- list()
  clades <- unique(data$clade)
  clades <- clades[-1]
  clades[K] <- "other" # putting the clade names in the right order
  num_days <- c(0:number_days)
  line_types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  combos <- list()
  for(k in 1:(stan$K - 1)){
    means <- get_posterior_mean(stan$mlr_fit, pars = c("raw_alpha", "raw_beta")) # the alpha and beta
    intercepts <- means[seq(k,length(means)/2,(K-1) )] # the alphas of the clade we want
    coef <- means[seq(k + (length(means)/2), length(means), (K-1)) ] # the beta's of the clade we want
    vals <- matrix(nrow = L, ncol = length(num_days))
    for( i in 1:L){
      vals[i, ] <- intercepts[i] + coef[i]*num_days
    }
    main = paste(clades[k], "Logits for each location") # creating the plots
    plot(num_days, vals[ 1, ], type = 'l', ylim = c(min(vals), max(vals)), ylab = "Logit (with respect to other)", xlab = paste("time from", start_date), col = colors[1], main = main)
    word <- paste(line_types[1], colors[1] )
    combos[[word]] <- stan$target_lo[1]
    for(i in 2:L){
      lines(num_days, vals[ i, ], col = colors[(i %% 7) + 1] , lty = line_types[(i %% 6) + 1])
      word <- paste(line_types[(i %% 6) + 1], colors[(i %% 7) + 1] )
      combos[[word]] <- stan$target_lo[i]
    }
    word <- clades[k]
    plotsave <- recordPlot()
    plots[[word]] <- plotsave
  }
  return( list(plots = plots, combos = combos))
}
