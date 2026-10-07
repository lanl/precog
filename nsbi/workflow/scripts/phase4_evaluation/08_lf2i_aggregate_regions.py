#!/usr/bin/env python3
"""
Aggregate LF2I confidence region batches and generate summary.

Takes all batch CSV files, concatenates them, and produces:
1. Final confidence_regions.csv (all observations)
2. Summary text file with comprehensive statistics

Inputs (Snakemake):
    - batches: List of batch CSV files
    - calibration: Fitted calibrator (for metadata)
    - test_data: Test data (for metadata)

Outputs:
    - regions: Final aggregated CSV
    - summary: Human-readable summary text file

Usage:
    Invoked by Snakemake rule aggregate_lf2i_regions
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Add scripts directory to path (must match import pattern in calibration script)
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.lf2i_calibration import LF2ICalibrator

# Snakemake I/O
batch_files = snakemake.input.batches
calibration_file = snakemake.input.calibration
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
write_log(f"LF2I Confidence Regions - Aggregation")
write_log(f"{'='*70}\n")

# Create output directory
Path(regions_output).parent.mkdir(parents=True, exist_ok=True)

# Load calibrator (for metadata)
write_log(f"Loading calibrator: {calibration_file}")
calibrator = LF2ICalibrator.load(calibration_file)
write_log(f"  Alpha: {calibrator.alpha}")
write_log(f"  Quantile: {calibrator.quantile}")
write_log(f"  Calibration samples: {len(calibrator.calibration_data)}\n")

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
write_log(f"Saving final confidence regions CSV...")
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
summary_lines.append("LF2I Calibration Summary")
summary_lines.append("="*70)
summary_lines.append("")
summary_lines.append("Calibration Set:")
summary_lines.append(f"  Samples: {len(calibrator.calibration_data)}")
summary_lines.append(f"  R0 range: [{calibrator.calibration_data['R0'].min():.2f}, {calibrator.calibration_data['R0'].max():.2f}]")
summary_lines.append(f"  Recovery time range: [{calibrator.calibration_data['recovery_time'].min():.2f}, {calibrator.calibration_data['recovery_time'].max():.2f}]")
summary_lines.append("")
summary_lines.append(f"Confidence Level: {(1-calibrator.alpha)*100:.0f}% (alpha={calibrator.alpha})")
summary_lines.append("")
summary_lines.append("Empirical Coverage (Test Set):")
summary_lines.append(f"  2D coverage: {results_df['coverage_2d'].mean():.1%} ({results_df['coverage_2d'].sum()}/{len(results_df)})")
summary_lines.append(f"  R0 marginal: {results_df['R0_coverage'].mean():.1%} ({results_df['R0_coverage'].sum()}/{len(results_df)})")
summary_lines.append(f"  Recovery time marginal: {results_df['recovery_time_coverage'].mean():.1%} ({results_df['recovery_time_coverage'].sum()}/{len(results_df)})")
summary_lines.append("")
summary_lines.append("Confidence Set Sizes:")
summary_lines.append(f"  Mean: {results_df['set_size'].mean():.1f} grid points")
summary_lines.append(f"  Median: {results_df['set_size'].median():.0f} grid points")
summary_lines.append(f"  Range: [{results_df['set_size'].min()}, {results_df['set_size'].max()}]")
summary_lines.append("")
summary_lines.append("Confidence Set Widths:")
summary_lines.append(f"  R0 mean width: {results_df['R0_width'].mean():.3f}")
summary_lines.append(f"  R0 median width: {results_df['R0_width'].median():.3f}")
summary_lines.append(f"  Recovery time mean width: {results_df['recovery_time_width'].mean():.3f}")
summary_lines.append(f"  Recovery time median width: {results_df['recovery_time_width'].median():.3f}")
summary_lines.append("")

# Interval Scores (if available)
if 'interval_score_R0' in results_df.columns and 'interval_score_recovery_time' in results_df.columns:
    summary_lines.append("Interval Scores (Test Set):")
    summary_lines.append(f"  R0: {results_df['interval_score_R0'].mean():.3f} ± {results_df['interval_score_R0'].std():.3f}")
    summary_lines.append(f"  Recovery time: {results_df['interval_score_recovery_time'].mean():.3f} ± {results_df['interval_score_recovery_time'].std():.3f}")
    summary_lines.append("")
    summary_lines.append("  Note: Interval Score = width + coverage penalties (lower is better)")
    summary_lines.append("        Balances sharpness (narrow intervals) vs coverage")
    summary_lines.append("")

summary_lines.append("MLE Statistics:")
summary_lines.append(f"  R0 MAE: {np.abs(results_df['mle_R0'] - results_df['true_R0']).mean():.3f}")
summary_lines.append(f"  Recovery time MAE: {np.abs(results_df['mle_recovery_time'] - results_df['true_recovery_time']).mean():.3f}")
summary_lines.append("")

# Warnings
if len(calibrator.calibration_data) < 50:
    summary_lines.append("⚠ WARNINGS:")
    summary_lines.append(f"  Small calibration set (n={len(calibrator.calibration_data)})")
    summary_lines.append("  Coverage estimates have high uncertainty")
    summary_lines.append("  Recommend n≥50 for reliable calibration")
    summary_lines.append("")

# Check if coverage is far from nominal
coverage_2d = results_df['coverage_2d'].mean()
if abs(coverage_2d - calibrator.quantile) > 0.15:
    if '⚠ WARNINGS:' not in summary_lines:
        summary_lines.append("⚠ WARNINGS:")
    summary_lines.append(f"  2D coverage ({coverage_2d:.1%}) deviates from nominal ({calibrator.quantile:.1%})")
    summary_lines.append("  This may indicate:")
    summary_lines.append("    - Small sample size effects")
    summary_lines.append("    - Model misspecification")
    summary_lines.append("    - Need for larger grid resolution")
    summary_lines.append("")

summary_lines.append("="*70)
summary_lines.append("Outputs:")
summary_lines.append(f"  CSV: {regions_output}")
summary_lines.append(f"  Calibration: {calibration_file}")
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
