// SIR Model for Bayesian Inference with CmdStanR
// Frequency-dependent transmission matching BEAST2 simulations
// Poisson likelihood on daily new recoveries using diff(R) for consistency with MLE

functions {
  // SIR ODE system with frequency-dependent transmission
  // Matches BEAST2 and MLE baseline implementation
  // Parameters passed as array for Stan ODE interface
  vector sir_ode(real t, vector y, array[] real theta) {
    real beta = theta[1];
    real gamma = theta[2];
    
    vector[3] dydt;
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real N = S + I + R;  // Total population (constant)
    
    // Frequency-dependent transmission: dS/dt = -(beta/N)*S*I
    dydt[1] = -(beta / N) * S * I;
    dydt[2] = (beta / N) * S * I - gamma * I;
    dydt[3] = gamma * I;
    
    return dydt;
  }
}

data {
  int<lower=1> n_obs;                    // Number of observed time points (days)
  array[n_obs] int<lower=0> y;          // Daily new recoveries (counts)
  array[n_obs] real<lower=0> ts;        // Day values (0, 1, 2, ...)
  real<lower=0> t0;                      // Initial time (0)
  
  // Initial conditions
  real<lower=0> S0_fixed;                // Fixed susceptible population
  real<lower=0> I0_fixed;                // Fixed initial infected (always 1)
  
  // Prior bounds (from simulation_params.yaml, unbuffered)
  real<lower=0> R0_min;
  real<lower=0> R0_max;
  real<lower=0> recovery_time_min;
  real<lower=0> recovery_time_max;
  real<lower=0, upper=1> prior_buffer;   // Buffer fraction (e.g., 0.1 = 10%)
}

transformed data {
  // Apply buffer to priors (same strategy as MLE for fair comparison)
  // This widens priors slightly to avoid edge effects in MCMC sampling
  real R0_prior_min = R0_min * (1.0 - prior_buffer);
  real R0_prior_max = R0_max * (1.0 + prior_buffer);
  real rt_prior_min = recovery_time_min * (1.0 - prior_buffer);
  real rt_prior_max = recovery_time_max * (1.0 + prior_buffer);
  real R0_init_fixed = 0.0;              // Initial recovered (always 0)
  
  // FIXED: Extend time grid to get proper diff(R) alignment
  // Day 0 corresponds to recoveries from t=0 to t=1, need R(0) and R(1)
  // Day n corresponds to recoveries from t=n to t=n+1, need R(n) and R(n+1)
  // IMPORTANT: Stan ODE requires t0 < ts[1], so we use t0 = -1e-10
  int n_times = n_obs + 1;
  array[n_times] real ts_extended;
  for (i in 1:n_obs) {
    ts_extended[i] = ts[i];
  }
  ts_extended[n_times] = ts[n_obs] + 1.0;  // Add one more time point
  
  // Stan ODE requires t0 < first observation time
  // If ts[1] == 0, use a tiny negative t0
  real t0_adjusted = ts[1] > 0 ? t0 : -1e-10;
}


parameters {
  real<lower=R0_prior_min, upper=R0_prior_max> R0;
  real<lower=rt_prior_min, upper=rt_prior_max> recovery_time;
}

transformed parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  
  // Convert from R0 and recovery_time to beta and gamma
  gamma = 1.0 / recovery_time;
  beta = R0 * gamma;
}

model {
  // Uniform priors (implicit from parameter constraints)
  
  // Initial conditions
  vector[3] y0;
  y0[1] = S0_fixed;
  y0[2] = I0_fixed;
  y0[3] = R0_init_fixed;
  
  // Solve ODE at extended time points to get diff(R) alignment
  array[n_times] vector[3] y_hat;
  array[2] real theta = {beta, gamma};
  y_hat = ode_rk45(sir_ode, y0, t0_adjusted, ts_extended, theta);
  
  // Likelihood: Poisson for daily new recoveries using diff(R)
  // diff(R)[i] = R(ts[i+1]) - R(ts[i]) = expected cases for day i
  for (i in 1:n_obs) {
    real expected_cases = y_hat[i+1, 3] - y_hat[i, 3];  // R(t+1) - R(t)
    expected_cases = fmax(expected_cases, 1e-6);  // Avoid zero
    y[i] ~ poisson(expected_cases);
  }
}

generated quantities {
  // Log-likelihood for model comparison
  real log_lik = 0;
  
  // Initial conditions
  vector[3] y0;
  y0[1] = S0_fixed;
  y0[2] = I0_fixed;
  y0[3] = R0_init_fixed;
  
  // Re-solve ODE (needed for generated quantities block)
  array[n_times] vector[3] y_hat;
  array[2] real theta = {beta, gamma};
  y_hat = ode_rk45(sir_ode, y0, t0_adjusted, ts_extended, theta);
  
  for (i in 1:n_obs) {
    real expected_cases = y_hat[i+1, 3] - y_hat[i, 3];
    expected_cases = fmax(expected_cases, 1e-6);
    log_lik += poisson_lpmf(y[i] | expected_cases);
  }
}

