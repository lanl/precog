#!/usr/bin/env python3
"""
Compute LF2I confidence regions for a BATCH of test observations.

This is the parallelized version - each batch processes a subset of test samples.
Outputs a partial CSV that will be aggregated later.

Inputs (Snakemake):
    - calibration: Fitted LF2I calibrator (pickle)
    - param_scan: Full parameter scan CSV
    - test_data: Test data (for extracting batch sim_ids)
    - params.sim_ids: List of sim_ids to process in this batch

Output:
    - Batch CSV with confidence region information (subset of observations)

Usage:
    Invoked by Snakemake rule compute_lf2i_regions_batch
"""

import sys
import traceback
import numpy as np
import pandas as pd
from pathlib import Path

# Snakemake I/O (get these first before any imports that might fail)
try:
    calibration_file = snakemake.input.calibration
    param_scan_file = snakemake.input.param_scan
    test_data_file = snakemake.input.test_data
    batch_sim_ids = snakemake.params.sim_ids
    batch_output = snakemake.output[0]
    log_file = snakemake.log[0]
except NameError as e:
    print(f"ERROR: snakemake object not available: {e}")
    print("This script must be run via Snakemake")
    sys.exit(1)

# Setup logging early
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
log = open(log_file, 'w')

def write_log(msg):
    """Write to both log file and flush immediately"""
    log.write(msg + '\n')
    log.flush()
    print(msg)  # Also print to stdout for Snakemake

# Wrap everything in try/except to catch any errors
try:
    write_log(f"{'='*70}")
    write_log(f"LF2I Confidence Regions - Batch Processing")
    write_log(f"{'='*70}\n")
    
    # Add scripts directory to path (not utils directly)
    # This must match the import pattern used in 01_calibrate_lf2i.py
    # so that pickle can find the utils.lf2i_calibration module
    scripts_path = str(Path(__file__).parent.parent)
    if scripts_path not in sys.path:
        sys.path.insert(0, scripts_path)
    write_log(f"Added to path: {scripts_path}")
    
    # Import with utils. prefix (must match how calibrator was pickled)
    from utils.lf2i_calibration import LF2ICalibrator
    write_log(f"✓ Imported LF2ICalibrator\n")

    # Create output directory
    Path(batch_output).parent.mkdir(parents=True, exist_ok=True)

    # Load calibrator
    write_log(f"Loading calibrator from: {calibration_file}")
    calibrator = LF2ICalibrator.load(calibration_file)
    write_log(f"  ✓ Loaded (alpha={calibrator.alpha}, quantile={calibrator.quantile})\n")

    # Load parameter scan (full file, we'll filter it)
    write_log(f"Loading parameter scan: {param_scan_file}")
    param_scan_df = pd.read_csv(param_scan_file)
    write_log(f"  Total rows: {len(param_scan_df)}")
    write_log(f"  Total observations: {param_scan_df['sim_id'].nunique()}\n")

    # Filter to batch sim_ids only
    write_log(f"Filtering to batch sim_ids...")
    write_log(f"  Batch size: {len(batch_sim_ids)} observations")
    param_scan_batch = param_scan_df[param_scan_df['sim_id'].isin(batch_sim_ids)].copy()
    write_log(f"  Filtered rows: {len(param_scan_batch)}")
    write_log(f"  Observations in batch: {param_scan_batch['sim_id'].nunique()}\n")

    # Validate all batch sim_ids are present
    missing_ids = set(batch_sim_ids) - set(param_scan_batch['sim_id'].unique())
    if missing_ids:
        write_log(f"⚠ WARNING: {len(missing_ids)} sim_ids missing from param_scan:")
        for mid in list(missing_ids)[:5]:
            write_log(f"    {mid}")
        if len(missing_ids) > 5:
            write_log(f"    ... and {len(missing_ids)-5} more")
        write_log("")

    # Load test data (for true parameters)
    write_log(f"Loading test data: {test_data_file}")
    test_data = np.load(test_data_file, allow_pickle=True)
    write_log(f"  Test samples: {len(test_data['sim_ids'])}\n")

    # Compute confidence sets for each observation in batch
    write_log(f"{'='*70}")
    write_log(f"Computing confidence regions for {len(batch_sim_ids)} observations...")
    write_log(f"{'='*70}\n")

    results = []
    for i, sim_id in enumerate(batch_sim_ids, 1):
        # Skip if not in param_scan (already warned above)
        if sim_id not in param_scan_batch['sim_id'].values:
            continue
        
        if i % 10 == 0 or i == 1 or i == len(batch_sim_ids):
            write_log(f"  [{i}/{len(batch_sim_ids)}] Processing {sim_id}...")
        
        # Compute confidence set (passing test_data for true parameters and coverage)
        cs_summary = calibrator.compute_confidence_set(sim_id, param_scan_batch, test_data=test_data)
        results.append(cs_summary)

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Save batch CSV
    write_log(f"\n{'='*70}")
    write_log(f"Saving batch results...")
    write_log(f"{'='*70}\n")

    results_df.to_csv(batch_output, index=False)
    write_log(f"✓ Batch CSV saved: {batch_output}")
    write_log(f"  Rows: {len(results_df)}")
    write_log(f"  Columns: {len(results_df.columns)}\n")

    # Quick stats for this batch
    write_log(f"Batch Statistics:")
    write_log(f"  2D coverage: {results_df['coverage_2d'].mean():.1%} ({results_df['coverage_2d'].sum()}/{len(results_df)})")
    write_log(f"  R0 coverage: {results_df['R0_coverage'].mean():.1%}")
    write_log(f"  RT coverage: {results_df['recovery_time_coverage'].mean():.1%}")
    write_log(f"  Mean set size: {results_df['set_size'].mean():.1f} points")
    write_log(f"  Mean R0 width: {results_df['R0_width'].mean():.3f}")
    write_log(f"  Mean RT width: {results_df['recovery_time_width'].mean():.3f}\n")

    write_log(f"{'='*70}")
    write_log(f"✓ Batch processing complete!")
    write_log(f"{'='*70}")

except Exception as e:
    # Log any errors that occur
    write_log(f"\n{'='*70}")
    write_log(f"ERROR: Batch processing failed!")
    write_log(f"{'='*70}\n")
    write_log(f"Exception: {type(e).__name__}")
    write_log(f"Message: {str(e)}\n")
    write_log(f"Traceback:")
    write_log(traceback.format_exc())
    log.close()
    sys.exit(1)

log.close()
