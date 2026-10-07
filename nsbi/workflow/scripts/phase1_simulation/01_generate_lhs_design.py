"""
Generate Latin Hypercube Sampling (LHS) design for SIR parameters.
This script can generate either train or test parameter sets.
"""

import numpy as np
import pandas as pd
from pyDOE2 import lhs

# Import centralized seed utility
import sys
import os
# Add scripts root to path (for imports from subdirectories)
scripts_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if scripts_root not in sys.path:
    sys.path.insert(0, scripts_root)
from utils.seed_utils import get_seed

# Access snakemake variables
dataset = snakemake.params.dataset  # "train" or "test"
n_samples = snakemake.params.n_samples
n_reps = snakemake.params.n_reps
S_min = snakemake.params.S_min
S_max = snakemake.params.S_max
S_fixed = snakemake.params.get('S_fixed', None)  # NEW: Optional fixed S value
R0_min = snakemake.params.R0_min
R0_max = snakemake.params.R0_max
rt_min = snakemake.params.rt_min
rt_max = snakemake.params.rt_max
seed = snakemake.params.seed

output_file = snakemake.output[0]

def generate_lhs_parameters(n_samples, seed, S_min, S_max, R0_min, R0_max, 
                           rt_min, rt_max, n_replicates, dataset_prefix, S_fixed=None):
    """
    Generate LHS samples for SIR parameters.
    
    Parameters:
    - n_samples: Number of unique parameter combinations
    - seed: Random seed for reproducibility
    - S_min, S_max: Bounds for susceptible population (ignored if S_fixed is set)
    - R0_min, R0_max: Bounds for basic reproductive number (log scale)
    - rt_min, rt_max: Bounds for recovery time (1/gamma)
    - n_replicates: Number of replicates per parameter combination (NOT USED - kept for compatibility)
    - dataset_prefix: Prefix for dataset (e.g., "train" or "test") (NOT USED - kept for compatibility)
    - S_fixed: If provided, fix S at this value (only sample R0 and recovery_time)
    
    Returns:
    - DataFrame with columns: param_id, S, R0, recovery_time, gamma, beta
      (One row per unique parameter set)
    """
    np.random.seed(seed)
    
    # Determine number of parameters to sample
    if S_fixed is not None:
        # Only sample R0 and recovery_time (2D LHS)
        n_params = 2
        print(f"  S fixed at {S_fixed} - sampling only R0 and recovery_time")
    else:
        # Sample all 3 parameters (3D LHS)
        n_params = 3
        print(f"  Sampling all 3 parameters: S, R0, recovery_time")
    
    # Generate LHS design (n_samples × n_params)
    # For n_samples=1, use random sampling instead of LHS (maximin fails)
    if n_samples == 1:
        lhs_design = np.random.rand(1, n_params)
    else:
        lhs_design = lhs(n_params, samples=n_samples, criterion='maximin')
    
    # Scale to parameter bounds
    if S_fixed is not None:
        # S is fixed, LHS columns are [R0, recovery_time]
        S_samples = np.full(n_samples, S_fixed)
        
        # R0: Log scale (geometric spacing)
        log_R0_min = np.log10(R0_min)
        log_R0_max = np.log10(R0_max)
        R0_samples = 10 ** (log_R0_min + (log_R0_max - log_R0_min) * lhs_design[:, 0])
        
        # Recovery time: Linear scale
        rt_samples = rt_min + (rt_max - rt_min) * lhs_design[:, 1]
    else:
        # All 3 parameters sampled, LHS columns are [S, R0, recovery_time]
        # S: Linear scale
        S_samples = S_min + (S_max - S_min) * lhs_design[:, 0]
        
        # R0: Log scale (geometric spacing)
        log_R0_min = np.log10(R0_min)
        log_R0_max = np.log10(R0_max)
        R0_samples = 10 ** (log_R0_min + (log_R0_max - log_R0_min) * lhs_design[:, 1])
        
        # Recovery time: Linear scale
        rt_samples = rt_min + (rt_max - rt_min) * lhs_design[:, 2]
    
    # Calculate derived parameters
    gamma_samples = 1.0 / rt_samples  # Recovery rate
    beta_samples = R0_samples * gamma_samples  # Infection rate
    
    # Create DataFrame - ONE ROW PER UNIQUE PARAMETER SET
    # Replicates will use the same parameters (BEAST stochasticity provides variation)
    records = []
    for param_id in range(n_samples):
        records.append({
            'param_id': param_id,
            'S': S_samples[param_id],
            'R0': R0_samples[param_id],
            'recovery_time': rt_samples[param_id],
            'gamma': gamma_samples[param_id],
            'beta': beta_samples[param_id]
        })
    
    df = pd.DataFrame(records)
    
    # Round for readability
    df['S'] = df['S'].round(0).astype(int)
    df['R0'] = df['R0'].round(4)
    df['recovery_time'] = df['recovery_time'].round(4)
    df['gamma'] = df['gamma'].round(6)
    df['beta'] = df['beta'].round(6)
    
    return df

# Generate design
print(f"Generating {dataset} design:")
print(f"  {n_samples} unique parameter sets")
print(f"  {n_reps} replicates per set will reuse the same XML (BEAST stochasticity gives variation)")
print(f"  Total simulations: {n_samples * n_reps}")

# Handle random seed generation
actual_seed = get_seed(seed)

if S_fixed is not None:
    print(f"  ⚠️  S is FIXED at {S_fixed} (not estimated by model)")
else:
    print(f"  S will be sampled from [{S_min}, {S_max}]")

# if n_samples > 0:
df = generate_lhs_parameters(n_samples, actual_seed, S_min, S_max, 
                             R0_min, R0_max, rt_min, rt_max, 
                             n_reps, dataset, S_fixed=S_fixed)
    
# Save to CSV
df.to_csv(output_file, index=False)
    
print(f"\n{dataset.capitalize()} design saved to: {output_file}")
print(f"\nParameter ranges:")
print(df[['S', 'R0', 'recovery_time', 'gamma', 'beta']].describe())
# else:
#     # Create empty dataframe with proper columns (NEW: no sim_id or replicate columns)
#     df = pd.DataFrame(columns=['param_id', 'S', 'R0', 'recovery_time', 'gamma', 'beta'])
#     df.to_csv(output_file, index=False)
    
#     print(f"  (Creating empty {dataset} set - n_samples = 0)")
#     print(f"\nEmpty {dataset} design saved to: {output_file}")

print(f"\n✓ {dataset.capitalize()} LHS design complete!")
