# this is the constant mechanistic model with s_v equal to the zero vector 
data {
  int<lower=1> L; // the number of locations
  int<lower=1> V; // the number of different clades
  int<lower=0> N; // the number of location, variant, times, trios
  vector[V - 1] s_v; // the vector of starting times for each variant
  array[N] int<lower=1, upper=V> y; // the clades
  array[N] int<lower=1, upper=L> ll; // the locations
  array[N] int t; // the days of the samples
  row_vector[N] weights; // number of seq per location
}
transformed data{
  vector[V] s_v_true = append_row(0, s_v);
}
parameters {
  array[L] vector<lower=0>[V-1] raw_beta; // the growth rates for the different variant-location pair
  array[L] vector<lower=0>[V-1] I; // the I_(s_v)'s
  real<lower=0> bsd; // the prior standard deviations for the variants' growth rates
  array[V-1] real bloc; // the prior means for variants' growth rates 
}
transformed parameters{
 array[L] vector[V-1] raw_alpha; # the intercepts
 for(v in 1:(V-1)){
   for( l in 1:L){
    raw_alpha[l,v] = log(expm1(raw_beta[l, v])) - raw_beta[l,v] + log(I[l,v]); 
   }
 }
}
model{
for (v in 1:(V-1)) {
    raw_beta[:, v] ~ normal(bloc[v], bsd);
    I[:, v] ~ normal(0, 3);
  }

    // Temporary variables to avoid repeated calculation
    array[L] vector[V] alpha;
    array[L] vector[V] beta;
    for (l in 1:L) {
      alpha[l] = append_row(0, raw_alpha[l]);
      beta[l] = append_row(0, raw_beta[l]);
    }

    // Vectorized likelihood calculation
    vector[N] log_prob;
    for (n in 1:N) {
      vector[V] time_coef;
      for(v in 1:V){
        time_coef[v] = beta[ll[n], v]*(t[n] - s_v_true[v]); # calculating the growth since t = s_v
      }
      log_prob[n] = categorical_logit_lpmf(y[n] | alpha[ll[n]] + time_coef); # here we let s_v vary
    }
    target += weights*log_prob;
  }

