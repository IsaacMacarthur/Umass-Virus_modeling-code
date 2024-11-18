data {
  int<lower=0> N; // the number of samples
  int<lower=1> L; // number of locations
  int<lower=1> K; // number of clades
  int<lower=1> B; // the number of basis functions
  array[N] int<lower=1, upper=K> y; // the clades
  array[N] int<lower=1, upper=L> ll; // locations
  matrix[B, N] x; // the spline over the days of the samples
  row_vector[N] weights; // number of seq per location 
}

parameters {
  vector<lower=0>[K-1] bsd; // prior sd for betas
  vector<lower=0>[K-1] asd; // prior sd for alphas
  array[K-1] real aloc; // prior means for the alpha's
  array[K-1] real bloc; // prior means for betas
  
  // Non-centered parameterization
  array[L] vector[K-1] alpha_nc; // raw alpha's
  array[L] matrix[K-1, B] beta_nc; // raw betas
}

transformed parameters {
  // Apply non-centered parameterization
  array[L] vector[K-1] raw_alpha;
  array[L] matrix[K-1, B] raw_beta;
  for (l in 1:L) {
    for (k in 1:(K-1)) {
      raw_alpha[l, k] = aloc[k] + asd[k] .* alpha_nc[l, k]; // Reparameterized alpha
      for(b in 1:(B)){
        raw_beta[l, k, b] = bloc[k] + bsd[k] * beta_nc[l, k,b]; // Reparameterized beta
      }
    }
  }
}

model {
  // Priors
  bsd ~ normal(1, 10000); // prior for the sd
  bloc ~ normal(0, 10000); // priors for the betas 
  asd ~ normal(1, 10000); // prior for the alpha sd
  aloc ~ normal(0, 10000); // prior for the alpha means

  // Priors on raw parameters (standard normal)
  for (l in 1:L) {
    alpha_nc[l] ~ normal(0, 1); // standard normal for alpha_nc
    for (k in 1:(K-1)) {
      beta_nc[l, k] ~ normal(0, 1); // standard normal for beta_nc
    }
  }

  // Likelihood
  {
    // Temporary variables to avoid repeated calculation
    array[L] vector[K] alpha;
    array[L] matrix[K, B] beta;
    for (l in 1:L) {
      alpha[l] = append_row(0, raw_alpha[l]);
      for (b in 1:B) {
        beta[l, :, b] = append_row(0, raw_beta[l, :, b]);
      }
    }
    // Vectorized likelihood calculation
    vector[N] log_prob;
    for (n in 1:N) {
      log_prob[n] = categorical_logit_lpmf(y[n] | alpha[ll[n]] + beta[ll[n], :, :] * x[:, n]);
    }
    target += weights * log_prob;
  }
}


