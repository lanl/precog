"""
Aggregate tsfeatures from individual CSV files into final dataset file.

Loads individual tsfeatures.csv files, merges with parameters, and creates
a single aggregated dataset file.
"""

import pandas as pd
import glob
import re
import os
import sys
from pathlib import Path

# Access snakemake variables
params_file = snakemake.input.params
output_file = snakemake.output[0]
dataset = snakemake.params.dataset
expected_sims = snakemake.params.expected_sims
log_file = snakemake.log[0]

# Setup logging
log_dir = Path(log_file).parent
log_dir.mkdir(parents=True, exist_ok=True)

with open(log_file, 'w') as log:
    log.write(f"Aggregating {dataset} tsfeatures\n")
    log.write(f"=" * 80 + "\n")
    log.write(f"Expected simulations: {expected_sims}\n\n")
    log.flush()
    
    # Load parameter file
    log.write(f"Loading parameter file: {params_file}\n")
    params_df = pd.read_csv(params_file)
    log.write(f"  Parameters: {len(params_df)} parameter sets\n\n")
    log.flush()
    
    # Find all tsfeatures files for this dataset
    tsfeatures_dir = "results/phase2_processing/features/temp/tsfeatures"
    pattern = f"{tsfeatures_dir}/{dataset}_p*_r*_tsfeatures.csv"
    tsfeatures_files = sorted(glob.glob(pattern))
    
    log.write(f"Found {len(tsfeatures_files)} tsfeatures files\n")
    log.flush()
    
    if len(tsfeatures_files) == 0:
        error_msg = f"ERROR: No tsfeatures files found: {pattern}"
        log.write(f"\n{error_msg}\n")
        raise FileNotFoundError(error_msg)
    
    # Load all tsfeatures files
    log.write(f"\nLoading tsfeatures files...\n")
    all_data = []
    
    for i, tsfeatures_file in enumerate(tsfeatures_files):
        if (i + 1) % 100 == 0:
            log.write(f"  Loaded {i + 1}/{len(tsfeatures_files)} files...\n")
            log.flush()
        
        # Parse sim_id from filename
        filename = os.path.basename(tsfeatures_file)
        match = re.match(r'(train|test)_p(\d{5})_r(\d+)_tsfeatures\.csv', filename)
        
        if not match:
            log.write(f"  WARNING: Could not parse filename: {filename}\n")
            continue
        
        dataset_name = match.group(1)
        param_id = int(match.group(2))
        replicate = int(match.group(3))
        sim_id = f"{dataset_name}_p{param_id:05d}_r{replicate}"
        
        # Load tsfeatures
        try:
            tsfeatures_df = pd.read_csv(tsfeatures_file)
            
            if len(tsfeatures_df) == 0:
                log.write(f"  WARNING: Empty file: {filename}\n")
                continue
            
            stats_row = tsfeatures_df.iloc[0].to_dict()
            param_row = params_df[params_df['param_id'] == param_id]
            
            if len(param_row) == 0:
                log.write(f"  WARNING: No parameters for param_id {param_id}\n")
                continue
            
            param_row = param_row.iloc[0]
            
            # Combine metadata, parameters, and tsfeatures
            row_data = {
                'sim_id': sim_id,
                'param_id': param_id,
                'replicate': replicate,
                'R0': param_row['R0'],
                'recovery_time': param_row['recovery_time']
            }
            row_data.update(stats_row)
            all_data.append(row_data)
            
        except Exception as e:
            log.write(f"  ERROR loading {filename}: {str(e)}\n")
            continue
    
    log.write(f"\nSuccessfully loaded {len(all_data)} simulations\n")
    log.flush()
    
    # Create final DataFrame
    final_df = pd.DataFrame(all_data)
    
    # Check if we got all expected simulations
    if len(final_df) < expected_sims:
        log.write(f"\nWARNING: Expected {expected_sims} simulations but got {len(final_df)}\n")
        log.write(f"Missing: {expected_sims - len(final_df)} simulations\n")
    elif len(final_df) > expected_sims:
        log.write(f"\nWARNING: Expected {expected_sims} simulations but got {len(final_df)}\n")
        log.write(f"Extra: {len(final_df) - expected_sims} simulations\n")
    else:
        log.write(f"\n✓ All {expected_sims} simulations accounted for\n")
    
    # Get tsfeature column names (all except metadata)
    metadata_cols = ['sim_id', 'param_id', 'replicate', 'R0', 'recovery_time']
    tsfeature_cols = [col for col in final_df.columns if col not in metadata_cols]
    
    log.write(f"\ntsfeatures summary:\n")
    log.write(f"  Total tsfeatures: {len(tsfeature_cols)}\n")
    log.write(f"  First 10 features: {tsfeature_cols[:10]}\n")
    log.write(f"  Last 10 features: {tsfeature_cols[-10:]}\n\n")
    
    # Check for missing values
    n_missing = final_df[tsfeature_cols].isna().sum().sum()
    if n_missing > 0:
        log.write(f"WARNING: {n_missing} missing values detected in tsfeatures\n")
        missing_per_feature = final_df[tsfeature_cols].isna().sum()
        missing_features = missing_per_feature[missing_per_feature > 0]
        log.write(f"Features with missing values:\n")
        for feat, count in missing_features.items():
            log.write(f"  {feat}: {count} missing\n")
    else:
        log.write(f"✓ No missing values in tsfeatures\n")
    
    # Save aggregated data
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    final_df.to_csv(output_file, index=False)
    log.write(f"\n✓ Aggregated data saved to: {output_file}\n")
    log.write(f"  Shape: {final_df.shape}\n")
    log.write(f"  Columns: {list(final_df.columns)[:10]}... (showing first 10)\n")
    
    log.write(f"\n{'='*80}\n")
    log.write(f"Aggregation complete!\n")
    log.flush()

print(f"✓ {dataset} tsfeatures aggregated: {len(all_data)} simulations")
print(f"  Output: {output_file}")
