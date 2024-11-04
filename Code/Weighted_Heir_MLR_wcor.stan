data {
  int<lower=0> N; // the number of samples
  int<lower=1> L; // number of locations
  int<lower=1> K; // number of clades
  array[N] int<lower=1, upper=K> y; // the clades
  array[N] int<lower=1, upper=L> ll; // locations
  array[N] real x; // the days of the samples
  row_vector [N]  weights; // number of seq per location 
}
parameters {
  array[L*(K-1)] corr_matrix[2] simga; // prior correlation for the variants 
  array[K-1] real bloc; // prior means for betas
  array[K-1] real aloc; // prior means for the alpha's
  array[L] vector[K-1] raw_alpha; // alpha's without the reference variant
  array[L] vector[K-1] raw_beta; // betas without the reference variant
}
model {
  bloc ~ normal(0, 10000); // priors for the betas 
  aloc ~ normal(0, 10000); // prior for the alpha means
  for(i in 1:(L)*(K-1)){
    simga[i] ~ lkj_corr(2);
  }
  array[L*(K-1)] vector[2] gamma; // temporary variables to store the alpha_beta pairs 
  array[L*(K-1)] vector[2] means;
  for (k in 1:(K-1)) {
    for(l in 1:L){
      gamma[l + L*(k-1)] = to_vector({raw_alpha[l,k], raw_beta[l,k]});
      means[l + L*(k-1)] = to_vector({aloc[k], bloc[k]});
      gamma[l + L*(k-1)] ~ multi_normal(means[l + L*(k-1)], simga[l + L*(k-1)]);
    }
  }
  {
    // Temporary variables to avoid repeated calculation
    array[L] vector[K] alpha;
    array[L] vector[K] beta;
    for (l in 1:L) {
      alpha[l] = append_row(0, raw_alpha[l]);
      beta[l] = append_row(0, raw_beta[l]);
    }

    // Vectorized likelihood calculation
    vector[N] log_prob;
    for (n in 1:N) {
      log_prob[n] = categorical_logit_lpmf(y[n] | alpha[ll[n]] + beta[ll[n]] * x[n]);
    }
    target += weights*log_prob;
  }
}

