data {
  int<lower=0> N; // the number of samples
  int<lower=1> L; // number of locations
  int<lower=1> K; // number of clades
  array[N] int<lower=1, upper=K> y; // the clades
  array[N] int<lower=1, upper=L> ll; // locations
  array[N] real x; // the days of the samples
  row_vector[N] weights; // number of seq per location 
}
parameters {
  vector<lower=0>[K-1] bsd; // prior sd for betas
  vector<lower=0>[K-1] asd; // prior sd for the alphas
  array[K-1] real bloc; // prior means for betas
  array[K-1] real aloc; // prior means for the alphas
  array[L] vector[K-1] alpha_noncentered; // non-centered alpha parameters
  array[L] vector[K-1] beta_noncentered; // non-centered beta parameters
}
transformed parameters {
  // Centered alpha and beta parameters

  // Calculated alpha_raw and beta_raw
  array[L] vector[K-1] raw_alpha;
  array[L] vector[K-1] raw_beta;

  // Apply non-centered transformation and back-calculate raw values
  for (l in 1:L) {
    for( k in 1:(K-1)){
      raw_alpha[l, k] = aloc[k] + asd[k]*alpha_noncentered[l,k];
      raw_beta[l,k] = bloc[k] + bsd[k]*beta_noncentered[l,k];
    }
  }
}
model {
  // Priors for hyperparameters
  bsd ~ normal(1, 1); // prior for beta standard deviation
  asd ~ normal(1, 1); // prior for alpha standard deviation
  bloc ~ normal(0, 1); // prior for beta mean
  aloc ~ normal(0, 1); // prior for alpha mean

  // Standard normal priors for non-centered parameters
  for (l in 1:L) {
    alpha_noncentered[l] ~ normal(0, 1);
    beta_noncentered[l] ~ normal(0, 1);
  }

  {
    // Temporary variables to avoid repeated calculation
    array[L] vector[K-1] alpha;
    array[L] vector[K-1] beta;
    for (l in 1:L) {
      alpha[l] = append_row(0, raw_alpha[l]);
      beta[l] = append_row(0, raw_beta[l]);
    }

    // Vectorized likelihood calculation
    vector[N] log_prob;
    for (n in 1:N) {
      log_prob[n] = categorical_logit_lpmf(y[n] | alpha[ll[n]] + beta[ll[n]] * x[n]);
    }
    target += weights * log_prob;
  }
}
