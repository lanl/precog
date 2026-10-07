"""
Extract BINNED case counts from BEAST trajectory files for tsfeatures calculation.

This script bins SAMPLING EVENTS (I→R transitions) into daily intervals.
Each output CSV contains only the actual outbreak duration (no zero-padding).

Key features:
- Bins RECOVERIES (sampling events) into 1-day intervals using np.histogram()
- Tracks I→R transitions (when individuals are sampled)
- No padding - each file contains only valid bins
- Parallel processing using multiprocessing.Pool
- Output: CSV files with binned daily sampling counts

Note: In phylodynamic inference, sampling occurs at recovery (I→R transition),
      so this correctly matches the tree statistic 'number_of_lineages'.
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
from multiprocessing import Pool
from functools import partial


def bin_case_counts(times, vals, bin_width=1.0):
    """Bin case counts into daily intervals - NO PADDING."""
    if len(times) == 0:
        return np.array([0.0])
    
    max_time = times.max()
    n_bins = int(np.ceil(max_time / bin_width)) + 1
    bins = np.arange(0, (n_bins + 1) * bin_width, bin_width)
    binned, _ = np.histogram(times, bins=bins, weights=vals)
    return binned[:n_bins]


def extract_binned_case_counts_one_file(sim_id, bin_width=1.0, verbose=False):
    """Extract binned case counts for a single simulation."""
    try:
        traj_file = f"results/phase1_simulation/simulations/{sim_id}.traj"
        output_file = f"results/phase2_processing/case_counts/{sim_id}_binned_case_counts.csv"
        
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Read trajectory file
        df = pd.read_csv(traj_file, sep='\t')
        df = df[df['Sample'] == 0]
        
        # Pivot to wide format
        df_wide = df.pivot_table(index='t', columns='population', values='value', aggfunc='first')
        df_wide = df_wide.reset_index().sort_values('t').reset_index(drop=True)
        
        # Calculate new recoveries (sampling events)
        # In phylodynamic inference, sampling occurs when individuals recover (I→R transition)
        # This matches the tree statistic 'number_of_lineages' (total sampled individuals)
        if 'R' in df_wide.columns:
            R = df_wide['R'].values
            new_recoveries = np.diff(R, prepend=R[0])
            new_recoveries[0] = 0
        else:
            new_recoveries = np.zeros(len(df_wide))
        
        # Bin the sampling events by day
        times = df_wide['t'].values
        binned_case_counts = bin_case_counts(times, new_recoveries, bin_width)
        
        # Save to CSV (simple format for R)
        output_df = pd.DataFrame({
            'day': np.arange(len(binned_case_counts)),
            'case_count': binned_case_counts
        })
        output_df.to_csv(output_file, index=False)
        
        if verbose:
            print(f"✓ {sim_id}: {len(binned_case_counts)} days")
        
        return (sim_id, True, None)
        
    except Exception as e:
        error_msg = f"{type(e).__name__}: {str(e)}"
        if verbose:
            print(f"✗ {sim_id}: {error_msg}")
        return (sim_id, False, error_msg)


# Main execution
if __name__ == "__main__":
    # Access snakemake variables
    sim_ids = snakemake.params.sim_ids
    batch_id = snakemake.wildcards.batch_id
    bin_width = snakemake.params.bin_width
    log_file = snakemake.log[0]
    n_threads = snakemake.threads
    
    # Setup logging
    log_dir = Path(log_file).parent
    log_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Extracting binned case counts - Batch {batch_id}")
    print(f"Processing {len(sim_ids)} files with {n_threads} cores")
    print(f"Bin width: {bin_width} day(s)")
    
    # Process files in parallel
    with Pool(processes=n_threads) as pool:
        results = pool.map(
            partial(extract_binned_case_counts_one_file, bin_width=bin_width, verbose=False), 
            sim_ids
        )
    
    # Analyze results
    successful = [sim_id for sim_id, success, _ in results if success]
    failed = [(sim_id, error) for sim_id, success, error in results if not success]
    
    # Write detailed log
    with open(log_file, 'w') as log:
        log.write(f"Binned Case Counts Extraction - Batch {batch_id}\n")
        log.write(f"=" * 80 + "\n")
        log.write(f"Bin width: {bin_width} day(s)\n")
        log.write(f"Output: Daily binned sampling events (I→R, no padding)\n")
        log.write(f"Total files: {len(sim_ids)}\n")
        log.write(f"Successful: {len(successful)}\n")
        log.write(f"Failed: {len(failed)}\n")
        log.write(f"Threads used: {n_threads}\n\n")
        
        if failed:
            log.write(f"Failed extractions:\n")
            for sim_id, error in failed:
                log.write(f"  - {sim_id}: {error}\n")
            
            # Write failures to separate file
            failure_file = log_dir / f"batch_{batch_id}_binned_failures.txt"
            with open(failure_file, 'w') as f:
                f.write('\n'.join([sim_id for sim_id, _ in failed]))
            log.write(f"\nFailure list saved to: {failure_file}\n")
    
    # Print summary
    print(f"Batch {batch_id} complete:")
    print(f"  Successful: {len(successful)}")
    print(f"  Failed: {len(failed)}")
    
    if failed:
        print(f"  See log for details: {log_file}")
