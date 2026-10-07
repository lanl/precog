#!/usr/bin/env python3
"""
Identify simulations that converged in both MLE and Bayesian methods,
then randomly select N for comparison visualization.

This script is designed for Snakemake integration and uses the snakemake object.
"""
import pandas as pd
import numpy as np

# Parse snakemake inputs
mle_df = pd.read_csv(snakemake.input.mle)
bayes_df = pd.read_csv(snakemake.input.bayes)
n_examples = snakemake.params.n_examples
seed = snakemake.params.seed

print("\n" + "="*70)
print("IDENTIFYING CONVERGED SIMULATIONS FOR COMPARISON")
print("="*70)

# Find converged simulations in each method
mle_conv = set(mle_df[mle_df['converged'] == True]['sim_id'])
bayes_conv = set(bayes_df[bayes_df['converged'] == True]['sim_id'])
common = sorted(mle_conv & bayes_conv)

print(f"\nMLE converged: {len(mle_conv)}")
print(f"Bayes converged: {len(bayes_conv)}")
print(f"Both converged: {len(common)}")

if len(common) == 0:
    print("\n⚠️  WARNING: No simulations converged in both methods!")
    print("Comparison visualization will not be generated.")
    # Create empty files to satisfy Snakemake dependencies
    with open(snakemake.output.converged, 'w') as f:
        pass
    with open(snakemake.output.selected, 'w') as f:
        pass
else:
    # Write all converged to file
    with open(snakemake.output.converged, 'w') as f:
        f.write('\n'.join(common))
    
    # Random selection for visualization
    np.random.seed(seed)
    
    # Handle n_examples = -1 or null (use all)
    if n_examples is None or n_examples < 0:
        n_to_select = len(common)
    else:
        n_to_select = min(n_examples, len(common))
    
    if n_to_select < len(common):
        selected = sorted(np.random.choice(common, n_to_select, replace=False))
        print(f"\nRandomly selected {n_to_select} simulations for visualization (seed={seed})")
    else:
        selected = common
        print(f"\nUsing all {len(common)} converged simulations for visualization")
    
    # Write selected to file
    with open(snakemake.output.selected, 'w') as f:
        f.write('\n'.join(selected))
    
    print("\nSelected simulations:")
    for sim_id in selected:
        mle_row = mle_df[mle_df['sim_id'] == sim_id].iloc[0]
        bayes_row = bayes_df[bayes_df['sim_id'] == sim_id].iloc[0]
        print(f"  {sim_id}: R0={mle_row['R0_true']:.2f}, D={mle_row['D_true']:.2f}, "
              f"Bayes ESS_R0={bayes_row['ess_bulk_R0']:.0f}")

print("\n" + "="*70)
print("✅ Identification complete!")
print(f"   Converged list: {snakemake.output.converged}")
print(f"   Selected list: {snakemake.output.selected}")
print("="*70)
