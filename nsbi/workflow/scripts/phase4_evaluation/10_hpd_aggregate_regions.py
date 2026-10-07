#!/usr/bin/env python3
"""
Aggregate HPD credible region batches and generate summary.

Takes all batch CSV files, concatenates them, and produces:
1. Final hpd_confidence_regions.csv (all observations)
2. Summary text file with comprehensive statistics

Inputs (Snakemake):
    - batches: List of batch CSV files
    - test_data: Test data (for metadata)

Outputs:
    - regions: Final aggregated CSV
    - summary: Human-readable summary text file

Usage:
    Invoked by Snakemake rule aggregate_hpd_regions
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.hpd_calibration import HPDCalibrator

# Snakemake I/O
batch_files = snakemake.input.batches
test_data_file = snakemake.input.test_data
regions_output = snakemake.output.regions
summary_output = snakemake.output.summary
log_file = snakemake.log[0]

# Setup logging
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
log = open(log_file, 'w')

def write_log(msg):
    """Write to both log file and flush immediately"""
    log.write(msg + '\n')
    log.flush()
    print(msg)  # Also to stdout

write_log(f"{'='*70}")
write_log(f"HPD Credible Regions - Aggregation")
write_log(f"{'='*70}\n")

# Create output directory
Path(regions_output).parent.mkdir(parents=True, exist_ok=True)

# Load test data (for metadata)
write_log(f"Loading test data: {test_data_file}")
test_data = np.load(test_data_file, allow_pickle=True)
write_log(f"  Test samples: {len(test_data['sim_ids'])}\n")

# Load and concatenate all batch files
write_log(f"Loading {len(batch_files)} batch files...")
batch_dfs = []
for i, batch_file in enumerate(batch_files, 1):
    df = pd.read_csv(batch_file)
    batch_dfs.append(df)
    if i % 10 == 0 or i == 1 or i == len(batch_files):
        write_log(f"  [{i}/{len(batch_files)}] Loaded {Path(batch_file).name}: {len(df)} rows")

write_log(f"\nConcatenating batches...")
results_df = pd.concat(batch_dfs, ignore_index=True)
write_log(f"  Total rows: {len(results_df)}")
write_log(f"  Unique observations: {results_df['sim_id'].nunique()}\n")

# Sort by sim_id for consistency
results_df = results_df.sort_values('sim_id').reset_index(drop=True)

# Save final CSV
write_log(f"Saving final HPD credible regions CSV...")
results_df.to_csv(regions_output, index=False)
write_log(f"  ✓ Saved: {regions_output}")
write_log(f"  Rows: {len(results_df)}")
write_log(f"  Columns: {len(results_df.columns)}\n")

# Generate summary statistics
write_log(f"{'='*70}")
write_log(f"Generating summary statistics...")
write_log(f"{'='*70}\n")

summary_lines = []
summary_lines.append("="*70)
summary_lines.append("HPD Credible Regions Summary")
summary_lines.append("="*70)
summary_lines.append("")
summary_lines.append("Method: Highest Posterior Density (HPD)")
summary_lines.append("  - Bayesian credible regions from neural ratio estimator")
summary_lines.append("  - Assumes flat prior")
summary_lines.append("  - No calibration step (unlike LF2I)")
summary_lines.append("")
summary_lines.append(f"Test Set: {len(results_df)} observations")
summary_lines.append(f"  R0 range: [{results_df['R0_true'].min():.2f}, {results_df['R0_true'].max():.2f}]")
summary_lines.append(f"  Recovery time (D) range: [{results_df['D_true'].min():.2f}, {results_df['D_true'].max():.2f}]")
summary_lines.append("")
summary_lines.append(f"Credible Level: {(1-results_df['alpha'].iloc[0])*100:.0f}% (alpha={results_df['alpha'].iloc[0]})")
summary_lines.append("")
summary_lines.append("Empirical Coverage (Test Set):")
summary_lines.append(f"  2D coverage: {results_df['coverage_2d'].mean():.1%} ({results_df['coverage_2d'].sum()}/{len(results_df)})")
summary_lines.append(f"  R0 marginal: {results_df['R0_coverage'].mean():.1%} ({results_df['R0_coverage'].sum()}/{len(results_df)})")
summary_lines.append(f"  D marginal: {results_df['D_coverage'].mean():.1%} ({results_df['D_coverage'].sum()}/{len(results_df)})")
summary_lines.append("")
summary_lines.append("Credible Set Sizes:")
summary_lines.append(f"  Mean: {results_df['set_size'].mean():.1f} grid points")
summary_lines.append(f"  Median: {results_df['set_size'].median():.0f} grid points")
summary_lines.append(f"  Range: [{results_df['set_size'].min()}, {results_df['set_size'].max()}]")
summary_lines.append("")
summary_lines.append("Credible Set Widths:")
summary_lines.append(f"  R0 mean width: {results_df['R0_interval_width'].mean():.3f}")
summary_lines.append(f"  R0 median width: {results_df['R0_interval_width'].median():.3f}")
summary_lines.append(f"  D mean width: {results_df['D_interval_width'].mean():.3f}")
summary_lines.append(f"  D median width: {results_df['D_interval_width'].median():.3f}")
summary_lines.append("")
summary_lines.append("Point Estimate Errors (MAE):")
summary_lines.append(f"  R0 MAE: {results_df['R0_abs_error'].mean():.3f}")
summary_lines.append(f"  D MAE: {results_df['D_abs_error'].mean():.3f}")
summary_lines.append("")

# Add probabilistic scores if present
if 'R0_crps' in results_df.columns:
    summary_lines.append("Probabilistic Scores (Test Set):")
    summary_lines.append(f"  CRPS (R0): {results_df['R0_crps'].mean():.4f} ± {results_df['R0_crps'].std():.4f}")
    summary_lines.append(f"  CRPS (D): {results_df['D_crps'].mean():.4f} ± {results_df['D_crps'].std():.4f}")
    summary_lines.append(f"  Energy Score: {results_df['energy_score'].mean():.4f} ± {results_df['energy_score'].std():.4f}")
    summary_lines.append(f"  Interval Score (R0): {results_df['R0_interval_score'].mean():.4f} ± {results_df['R0_interval_score'].std():.4f}")
    summary_lines.append(f"  Interval Score (D): {results_df['D_interval_score'].mean():.4f} ± {results_df['D_interval_score'].std():.4f}")
    summary_lines.append(f"  Log Score: {results_df['log_score'].mean():.2f} ± {results_df['log_score'].std():.2f}")
    summary_lines.append("")
    summary_lines.append("  Notes:")
    summary_lines.append("    - Lower CRPS, Energy Score, and Interval Scores are better")
    summary_lines.append("    - Higher Log Score is better")
    summary_lines.append("")

# Check if coverage deviates from nominal
nominal_coverage = 1 - results_df['alpha'].iloc[0]
coverage_2d = results_df['coverage_2d'].mean()
if abs(coverage_2d - nominal_coverage) > 0.10:
    summary_lines.append("⚠ NOTES:")
    summary_lines.append(f"  2D coverage ({coverage_2d:.1%}) deviates from nominal ({nominal_coverage:.1%})")
    summary_lines.append("  This is expected for HPD without calibration.")
    summary_lines.append("  The neural estimator may have systematic biases.")
    summary_lines.append("  Compare with LF2I (calibrated) results for reference.")
    summary_lines.append("")

summary_lines.append("="*70)
summary_lines.append("Comparison with LF2I:")
summary_lines.append("  - HPD: Bayesian, interpretable, no calibration")
summary_lines.append("  - LF2I: Frequentist, calibrated for coverage guarantee")
summary_lines.append("  - Both use same neural ratio estimator")
summary_lines.append("  - HPD typically has narrower intervals (may be optimistic)")
summary_lines.append("="*70)
summary_lines.append("")
summary_lines.append("Outputs:")
summary_lines.append(f"  CSV: {regions_output}")
summary_lines.append("="*70)

summary_text = "\n".join(summary_lines)

# Write summary
with open(summary_output, 'w') as f:
    f.write(summary_text)

write_log(summary_text)
write_log(f"\n✓ Summary written to: {summary_output}\n")

write_log(f"{'='*70}")
write_log(f"✓ Aggregation complete!")
write_log(f"{'='*70}")

log.close()

# Print summary to stdout for Snakemake log
print("\n" + summary_text)

