"""
Aggregate Parameter Scan Batch Results

Combines individual batch CSV files from parameter scan evaluation into a single
output file and validates completeness. Also validates that all batches used
the same parameter grid for consistency.
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

# Access snakemake variables
batch_files = snakemake.input.batch_files
grid_file = snakemake.input.grid_file
test_data_file = snakemake.input.test_data
output_file = snakemake.output[0]
n_r0_points = snakemake.params.n_r0_points
n_rt_points = snakemake.params.n_rt_points
log_file = snakemake.log[0]

# Setup logging
log_dir = Path(log_file).parent
log_dir.mkdir(parents=True, exist_ok=True)

with open(log_file, 'w') as log:
    log.write("="*70 + "\n")
    log.write("Aggregating Parameter Scan Batches\n")
    log.write("="*70 + "\n")
    log.write(f"Input batches: {len(batch_files)}\n")
    log.write(f"Grid file: {grid_file}\n")
    log.write(f"Test data file: {test_data_file}\n\n")
    log.flush()
    
    # Load test data to get actual number of test samples (after NaN filtering)
    log.write("Loading test data to determine expected row count...\n")
    test_data = np.load(test_data_file, allow_pickle=True)
    n_test_samples = len(test_data['sim_ids'])
    expected_rows = n_test_samples * n_r0_points * n_rt_points
    log.write(f"  Test samples (after filtering): {n_test_samples}\n")
    log.write(f"  R0 grid points: {n_r0_points}\n")
    log.write(f"  RT grid points: {n_rt_points}\n")
    log.write(f"  Expected total rows: {expected_rows}\n\n")
    log.flush()
    
    # Load grid for validation
    log.write("Loading reference grid...\n")
    grid_data = np.load(grid_file, allow_pickle=True)
    expected_r0 = sorted(grid_data['r0_grid'])
    expected_rt = sorted(grid_data['rt_grid'])
    grid_metadata = grid_data['metadata'].item()
    log.write(f"  R0 grid points: {len(expected_r0)}\n")
    log.write(f"  RT grid points: {len(expected_rt)}\n")
    log.write(f"  Grid metadata: {grid_metadata}\n\n")
    log.flush()
    
    # Load all batch files
    batch_dfs = []
    for batch_file in batch_files:
        log.write(f"Loading {batch_file}... ")
        log.flush()
        
        df = pd.read_csv(batch_file)
        batch_dfs.append(df)
        log.write(f"✓ ({len(df)} rows)\n")
        log.flush()
    
    # Concatenate all batches
    combined_df = pd.concat(batch_dfs, ignore_index=True)
    log.write(f"\nCombined: {len(combined_df)} total rows\n")
    log.flush()
    
    # Validate completeness
    if len(combined_df) != expected_rows:
        missing = expected_rows - len(combined_df)
        pct_missing = (missing / expected_rows) * 100 if expected_rows > 0 else 0
        error_msg = (f"ERROR: Expected {expected_rows} rows, "
                    f"but got {len(combined_df)}\n"
                    f"  Missing: {missing} rows ({pct_missing:.2f}%)\n"
                    f"  Based on {n_test_samples} test samples (after NaN filtering)\n"
                    f"  Grid: {n_r0_points} R0 points × {n_rt_points} RT points\n"
                    f"  This suggests some parameter scan batches are incomplete or missing.")
        log.write(f"\n{error_msg}\n")
        raise ValueError(error_msg)
    
    # Check for duplicates
    dup_mask = combined_df.duplicated(subset=['sim_id', 'R0', 'recovery_time'], keep=False)
    duplicates = combined_df[dup_mask]
    if len(duplicates) > 0:
        error_msg = f"ERROR: Found {len(duplicates)} duplicate rows"
        log.write(f"\n{error_msg}\n")
        log.write(f"Duplicates:\n{duplicates.head()}\n")
        raise ValueError(error_msg)
    
    # Validate grid consistency
    log.write("\nValidating grid consistency...\n")
    actual_r0 = sorted(combined_df['R0'].unique())
    actual_rt = sorted(combined_df['recovery_time'].unique())
    
    r0_match = np.allclose(actual_r0, expected_r0, rtol=1e-10)
    rt_match = np.allclose(actual_rt, expected_rt, rtol=1e-10)
    
    if not r0_match:
        error_msg = (f"ERROR: R0 values in data don't match expected grid!\n"
                    f"Expected {len(expected_r0)} unique R0 values, "
                    f"got {len(actual_r0)}")
        log.write(f"\n{error_msg}\n")
        log.write(f"Expected R0 (first 5): {expected_r0[:5]}\n")
        log.write(f"Actual R0 (first 5): {actual_r0[:5]}\n")
        raise ValueError(error_msg)
    
    if not rt_match:
        error_msg = (f"ERROR: recovery_time values in data don't match expected grid!\n"
                    f"Expected {len(expected_rt)} unique RT values, "
                    f"got {len(actual_rt)}")
        log.write(f"\n{error_msg}\n")
        log.write(f"Expected RT (first 5): {expected_rt[:5]}\n")
        log.write(f"Actual RT (first 5): {actual_rt[:5]}\n")
        raise ValueError(error_msg)
    
    log.write(f"  ✓ R0 grid matches: {len(actual_r0)} unique values\n")
    log.write(f"  ✓ RT grid matches: {len(actual_rt)} unique values\n")
    log.write(f"  ✓ All batches used consistent grid!\n")
    log.flush()
    
    # Sort by sim_id, R0, recovery_time for consistent output
    combined_df = combined_df.sort_values(['sim_id', 'R0', 'recovery_time']).reset_index(drop=True)
    
    # Save aggregated results
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_df.to_csv(output_file, index=False)
    
    log.write("\n" + "="*70 + "\n")
    log.write("AGGREGATION SUMMARY\n")
    log.write("="*70 + "\n")
    log.write(f"Total rows: {len(combined_df)}\n")
    log.write(f"Unique test cases: {combined_df['sim_id'].nunique()}\n")
    log.write(f"Grid points per case: {len(combined_df) // combined_df['sim_id'].nunique()}\n")
    log.write(f"Output: {output_file}\n")
    log.write(f"\nPrediction Statistics:\n")
    log.write(f"  Mean: {combined_df['predicted_probability'].mean():.4f}\n")
    log.write(f"  Std: {combined_df['predicted_probability'].std():.4f}\n")
    log.write(f"  Min: {combined_df['predicted_probability'].min():.4f}\n")
    log.write(f"  Max: {combined_df['predicted_probability'].max():.4f}\n")
    log.write(f"\n✓ Aggregation complete!\n")

print(f"Aggregated {len(combined_df)} rows from {len(batch_files)} batches -> {output_file}")
