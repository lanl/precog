#!/usr/bin/env python3
"""
Compare MLE vs Bayesian inference by visualizing posterior draws with MLE contours.

This script works both as a standalone script and within Snakemake.
- Snakemake mode: Reads selected_sims.txt from snakemake.input.selected
- Standalone mode: Reads from command-line args or auto-detects from comparison dir
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import yaml
import sys
import os

# Determine execution mode
try:
    # Snakemake mode
    selected_file = snakemake.input.selected
    output_file = snakemake.output.plot
    comparison_dir = "results/phase5_baseline/comparison"
    
    # Read selected simulations from file
    with open(selected_file, 'r') as f:
        sim_ids = [line.strip() for line in f if line.strip()]
    
    print(f"Running in Snakemake mode")
    print(f"Selected simulations file: {selected_file}")
    
except NameError:
    # Standalone mode
    comparison_dir = "results/phase5_baseline/comparison"
    output_file = "reports/mle_vs_bayes_comparison.pdf"
    
    if len(sys.argv) > 1:
        # Use command-line arguments
        sim_ids = sys.argv[1:]
    else:
        # Auto-detect from comparison directory
        files = os.listdir(comparison_dir)
        sim_ids = sorted(set([f.replace("_mle_cov.csv", "").replace("_bayes_draws.csv", "")
                              .replace("_mle_summary.yaml", "").replace("_bayes_summary.yaml", "") 
                              for f in files if "test_p" in f]))
        sim_ids = [s for s in sim_ids if os.path.exists(f"{comparison_dir}/{s}_mle_cov.csv") and 
                                           os.path.exists(f"{comparison_dir}/{s}_bayes_draws.csv")]
    
    print(f"Running in standalone mode")

os.makedirs(os.path.dirname(output_file), exist_ok=True)

if len(sim_ids) == 0:
    print("\n⚠️  WARNING: No simulations found for comparison!")
    print("   This may be because no simulations converged in both methods.")
    print("   Skipping plot generation.")
    # Create empty output file to satisfy Snakemake
    with open(output_file, 'w') as f:
        pass
    sys.exit(0)

print(f"\nGenerating comparison plots for {len(sim_ids)} simulations:")
for sim_id in sim_ids:
    print(f"  - {sim_id}")

def confidence_ellipse(mean, cov, ax, n_std=1.96, facecolor='none', edgecolor='red', linewidth=2, label=None):
    """Draw confidence ellipse from mean and covariance matrix."""
    eigvals, eigvecs = np.linalg.eigh(cov)
    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
    width, height = 2 * n_std * np.sqrt(eigvals)
    ellip = Ellipse(xy=mean, width=width, height=height, angle=angle,
                   facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth,
                   label=label, alpha=0.7, zorder=10)
    ax.add_patch(ellip)
    return ellip

n_sims = len(sim_ids)

# Calculate grid dimensions for roughly square layout
if n_sims == 1:
    nrows, ncols = 1, 1
elif n_sims == 2:
    nrows, ncols = 1, 2
else:
    ncols = int(np.ceil(np.sqrt(n_sims)))
    nrows = int(np.ceil(n_sims / ncols))

print(f"\nCreating {nrows}x{ncols} grid for {n_sims} plots")

fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))

# Flatten axes array for easier iteration
if n_sims == 1:
    axes = [axes]
else:
    axes = axes.flatten()

for idx, sim_id in enumerate(sim_ids):
    ax = axes[idx]
    bayes_draws = pd.read_csv(f"{comparison_dir}/{sim_id}_bayes_draws.csv")
    mle_cov = pd.read_csv(f"{comparison_dir}/{sim_id}_mle_cov.csv", index_col=0)
    with open(f"{comparison_dir}/{sim_id}_mle_summary.yaml", 'r') as f:
        mle_meta = yaml.safe_load(f)
    with open(f"{comparison_dir}/{sim_id}_bayes_summary.yaml", 'r') as f:
        bayes_meta = yaml.safe_load(f)
    
    R0_true = mle_meta['R0_true']
    D_true = mle_meta['D_true']
    R0_mle = mle_meta['R0_mle']
    D_mle = mle_meta['D_mle']
    R0_bayes = bayes_meta['R0_mean']
    D_bayes = bayes_meta['D_mean']
    cov_matrix = mle_cov.values
    mle_mean = np.array([R0_mle, D_mle])
    
    ax.scatter(bayes_draws['R0'], bayes_draws['D'], 
              alpha=0.3, s=10, c='blue', label='Bayesian posterior', zorder=1)
    confidence_ellipse(mle_mean, cov_matrix, ax, n_std=1.96, 
                      edgecolor='red', linewidth=2.5, label='MLE 95% CI')
    ax.plot(R0_true, D_true, 'go', markersize=12, markeredgewidth=2, 
           markerfacecolor='lightgreen', label='True values', zorder=20)
    ax.plot(R0_mle, D_mle, 'rs', markersize=10, label='MLE estimate', zorder=15)
    ax.plot(R0_bayes, D_bayes, 'b^', markersize=10, label='Bayes mean', zorder=15)
    
    ax.set_xlabel('$R_0$', fontsize=12)
    ax.set_ylabel('Recovery time $D$ (days)', fontsize=12)
    title_text = f'{sim_id}\n$R_0$={R0_true:.2f}, $D$={D_true:.2f}'
    ax.set_title(title_text, fontsize=11)
    ax.legend(loc='best', fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    
    R0_min_plot = min(bayes_draws['R0'].min(), R0_true, R0_mle, R0_bayes) - 0.5
    R0_max_plot = max(bayes_draws['R0'].max(), R0_true, R0_mle, R0_bayes) + 0.5
    D_min_plot = min(bayes_draws['D'].min(), D_true, D_mle, D_bayes) - 1.0
    D_max_plot = max(bayes_draws['D'].max(), D_true, D_mle, D_bayes) + 1.0
    ax.set_xlim(R0_min_plot, R0_max_plot)
    ax.set_ylim(D_min_plot, D_max_plot)
    
    rhat_r0 = bayes_meta['rhat_R0']
    ess_r0 = int(bayes_meta['ess_bulk_R0'])
    textstr = f'Bayes: Rhat={rhat_r0:.3f}, ESS={ess_r0}\nMLE NLL={mle_meta["nll"]:.1f}'
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=8,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Hide unused subplots
for idx in range(n_sims, len(axes)):
    axes[idx].set_visible(False)

plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nPlot saved to: {output_file}")

print("\n" + "="*70)
print("Summary Comparison")
print("="*70)
for sim_id in sim_ids:
    with open(f"{comparison_dir}/{sim_id}_mle_summary.yaml", 'r') as f:
        mle_meta = yaml.safe_load(f)
    with open(f"{comparison_dir}/{sim_id}_bayes_summary.yaml", 'r') as f:
        bayes_meta = yaml.safe_load(f)
    
    print(f"\n{sim_id}:")
    print(f"  True: R0={mle_meta['R0_true']:.3f}, D={mle_meta['D_true']:.3f}")
    print(f"  MLE:  R0={mle_meta['R0_mle']:.3f}, D={mle_meta['D_mle']:.3f}")
    print(f"  Bayes: R0={bayes_meta['R0_mean']:.3f}, D={bayes_meta['D_mean']:.3f}")
    print(f"  MLE error:   R0={abs(mle_meta['R0_mle']-mle_meta['R0_true']):.3f}, D={abs(mle_meta['D_mle']-mle_meta['D_true']):.3f}")
    print(f"  Bayes error: R0={abs(bayes_meta['R0_mean']-mle_meta['R0_true']):.3f}, D={abs(bayes_meta['D_mean']-mle_meta['D_true']):.3f}")
    print(f"  Bayes diagnostics: Rhat_R0={bayes_meta['rhat_R0']:.3f}, ESS_R0={int(bayes_meta['ess_bulk_R0'])}")

print("\n" + "="*70)
print("✅ Comparison complete!")
print("="*70)




