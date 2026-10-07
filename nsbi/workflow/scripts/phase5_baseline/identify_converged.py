#!/usr/bin/env python3
"""
Identify simulations that converged in both MLE and Bayesian methods.
"""
import pandas as pd
import sys

# Load results
mle_results = pd.read_csv('results/phase5_baseline/sir_mle/results.csv')
bayes_results = pd.read_csv('results/phase5_baseline/sir_bayes/results.csv')

# Filter for converged only
mle_converged = mle_results[mle_results['converged'] == True]['sim_id'].tolist()
bayes_converged = bayes_results[bayes_results['converged'] == True]['sim_id'].tolist()

# Find common converged simulations
common_converged = sorted(set(mle_converged) & set(bayes_converged))

print(f"\nMLE converged: {len(mle_converged)}")
print(f"Bayes converged: {len(bayes_converged)}")
print(f"Both converged: {len(common_converged)}")
print(f"\nCommon converged simulations:")
for sim_id in common_converged:
    mle_row = mle_results[mle_results['sim_id'] == sim_id].iloc[0]
    bayes_row = bayes_results[bayes_results['sim_id'] == sim_id].iloc[0]
    print(f"  {sim_id}: R0_true={mle_row['R0_true']:.2f}, D_true={mle_row['D_true']:.2f}, "
          f"Bayes ESS_R0={bayes_row['ess_bulk_R0']:.0f}, ESS_D={bayes_row['ess_bulk_recovery_time']:.0f}")

# Select 3 representative examples
if len(common_converged) >= 3:
    # Pick diverse parameter values
    selected = [common_converged[0], common_converged[len(common_converged)//2], common_converged[-1]]
    print(f"\nSelected examples for comparison: {selected}")
    print(" ".join(selected))
else:
    print(f"\nAll {len(common_converged)} examples selected")
    print(" ".join(common_converged))
