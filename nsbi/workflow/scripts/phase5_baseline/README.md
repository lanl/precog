# Phase 5: Baseline Methods

This phase implements standard epidemiological methods for comparison with the neural SBI approach.

## SIR MLE Baseline

Maximum likelihood estimation of SIR model parameters using 1-day binned case counts (recoveries).

### Method

1. **Model**: Standard SIR ODE with frequency-dependent transmission (matching BEAST2)
   - dS/dt = -(β/N) S I
   - dI/dt = (β/N) S I - γ I
   - dR/dt = γ I
   - Where N = S + I + R (constant)

2. **Parameters Estimated**:
   - β (transmission rate)
   - γ (recovery rate)
   - Derived: R0 = β/γ, recovery time = 1/γ

3. **Data**: 
   - **New recoveries** (R increases) binned at 1-day resolution (consistent with NSBI)
   - This matches the phylogenetic tree construction (sampling from R compartment)
   - NOT new infections (S decreases)

4. **Likelihood**: Poisson likelihood on binned new recoveries
   - Expected recoveries[t] = R[t] - R[t-1] from ODE solution
   - Observed recoveries ~ Poisson(Expected recoveries)

5. **Optimization**: L-BFGS-B method with parameter bounds
   - Bounds derived dynamically from simulation priors with 10% buffer:
     - γ ∈ [(1/rt_max)×0.9, (1/rt_min)×1.1] (e.g., [0.064, 1.1] for rt ∈ [1, 14])
     - β ∈ [(R0_min/rt_max)×0.9, (R0_max/rt_min)×1.1] (e.g., [0.064, 11.0] for R0 ∈ [1, 10])
   - Buffer allows for stochastic variation and model misspecification
   - Prevents unrealistic parameter combinations while avoiding boundary constraints

6. **Starting Values**: 
   - R0: geometric mean of prior bounds (log scale)
   - Recovery time: arithmetic mean of prior bounds
   - Uses prior knowledge, not true values

7. **Uncertainty Quantification**:
   - 95% confidence intervals via Monte Carlo sampling from Hessian
   - Accounts for β-γ correlation
   - 10,000 samples from multivariate normal

### Known Limitations

1. **Temporal Binning Loss**: 
   - For fast outbreaks (short recovery time ~1-2 days), 1-day binning loses significant temporal information
   - Can lead to identifiability issues (multiple parameter sets produce similar binned data)
   - Example: High β + high γ can mimic low β + low γ when only observing recoveries

2. **Model Mismatch**:
   - MLE fits deterministic ODE to stochastic simulation output
   - Poisson observation model may underestimate uncertainty from epidemic stochasticity

3. **Coverage**: 
   - 95% CIs may have lower than nominal coverage due to (1) and (2)

### Implementation

- **Batching**: 500 simulations per batch for parallel processing
- **Inputs**: 
  - Case counts: `results/summaries/{sim_id}_case_counts.csv`
  - True parameters: `results/parameters/test_design.csv`
- **Outputs**:
  - `test_sir_mle_results.csv`: Consolidated results with 34 columns:
    - Columns 1-26: Per-simulation metadata, true parameters, estimates, confidence intervals, and fit statistics
    - Columns 27-34: Diagnostic metrics (error, absolute error, coverage, CI width for R0 and recovery_time)
    - **Note**: Diagnostic columns (27-34) contain NA for non-converged fits (where `converged = FALSE`)
  - `test_sir_mle_summary.yaml`: Overall performance metrics computed from converged fits only

### Performance Metrics

The summary includes:
- Convergence rate
- MAE (mean absolute error)
- RMSE (root mean squared error)
- Bias
- Coverage (% of CIs containing true value)
- Mean CI width

### Running

```bash
# Run only baseline
snakemake --cores N --use-conda all_baseline

# Baseline runs by default with full pipeline
snakemake --cores N --use-conda all
```

### Comparison Visualization

When both MLE and Bayesian baselines are enabled, the workflow automatically:
1. Identifies simulations that converged in both methods
2. Randomly selects N examples (default: 3, configurable in `config/bayesian_params.yaml`)
3. Re-fits selected simulations to extract full covariance matrices (MLE) and posterior draws (Bayes)
4. Generates comparison plot: `reports/mle_vs_bayes_comparison.pdf`

The comparison uses random sampling with the workflow seed for reproducibility.
To change the number of examples, edit `bayesian_baseline.comparison.n_examples` in config.

**Outputs:**
- `results/phase5_baseline/comparison/converged_sims.txt` - List of all converged sim_ids
- `results/phase5_baseline/comparison/selected_sims.txt` - Randomly selected subset
- `results/phase5_baseline/comparison/{sim_id}_*` - Detailed comparison data per simulation
- `reports/mle_vs_bayes_comparison.pdf` - Multi-panel visualization

### SIR Fit Visualization

To validate whether the MLE and Bayesian methods produce reasonable epidemic dynamics, use the `plot_sir_fits.py` script:

```bash
cd /Users/zenalapp/Documents/research/forecasting_dr/neural_sbi
python workflow/scripts/phase5_baseline/plot_sir_fits.py
```

This generates `reports/sir_fits_comparison.pdf` showing:
- **Left column**: I compartment (infected prevalence) over time
- **Right column**: R compartment (cumulative recoveries) vs observed case counts

Each plot compares three parameter sets:
- **Green lines**: True parameters from simulation
- **Red dashed lines**: MLE estimates
- **Blue dotted lines**: Bayes mean estimates
- **Black points**: Observed case count data (right column only)

**Purpose**: This visualization helps determine if poor parameter estimates are due to:
1. **Implementation issues**: If fits don't match data, methods need revisiting
2. **Model identifiability**: If fits match data despite wrong parameters, multiple parameter sets can produce similar epidemic curves (expected for binned data)

The script automatically processes all simulations in `results/phase5_baseline/comparison/selected_sims.txt`.


### Dependencies

R packages (conda environment `r-baseline.yaml`):
- deSolve (ODE solver)
- bbmle (MLE fitting)
- mvtnorm (confidence intervals)
- yaml (output formatting)

### Notes

- Only runs on TEST simulations (not training data)
- S (population size) is fixed and known (not estimated)
- Frequency-dependent transmission: R0 = β/γ (independent of N)
- Convergence issues are logged but don't block pipeline

