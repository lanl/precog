"""
Aggregate tree features from individual CSV files into final dataset file.

Loads individual treefeatures.csv files, merges with parameters, and creates
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
    log.write(f"Aggregating {dataset} Tree Features\n")
    log.write(f"=" * 80 + "\n")
    log.write(f"Expected simulations: {expected_sims}\n\n")
    log.flush()
    
    # Load parameter file
    log.write(f"Loading parameter file: {params_file}\n")
    params_df = pd.read_csv(params_file)
    log.write(f"  Parameters: {len(params_df)} parameter sets\n\n")
    log.flush()
    
    # Find all tree features files for this dataset
    treefeatures_dir = "results/phase2_processing/features/temp/treefeatures"
    pattern = f"{treefeatures_dir}/{dataset}_p*_r*_treefeatures.csv"
    treefeatures_files = sorted(glob.glob(pattern))
    
    log.write(f"Found {len(treefeatures_files)} tree features files\n")
    log.flush()
    
    if len(treefeatures_files) == 0:
        error_msg = f"ERROR: No tree features files found: {pattern}"
        log.write(f"\n{error_msg}\n")
        raise FileNotFoundError(error_msg)
    
    # Load all tree features files
    log.write(f"\nLoading tree features files...\n")
    all_data = []
    
    for i, tree_file in enumerate(treefeatures_files):
        if (i + 1) % 100 == 0:
            log.write(f"  Loaded {i + 1}/{len(treefeatures_files)} files...\n")
            log.flush()
        
        # Parse sim_id from filename
        filename = os.path.basename(tree_file)
        match = re.match(r'(train|test)_p(\d{5})_r(\d+)_treefeatures\.csv', filename)
        
        if not match:
            log.write(f"  WARNING: Could not parse filename: {filename}\n")
            continue
        
        dataset_name = match.group(1)
        param_id = int(match.group(2))
        replicate = int(match.group(3))
        sim_id = f"{dataset_name}_p{param_id:05d}_r{replicate}"
        
        # Load tree features
        try:
            treefeatures_df = pd.read_csv(tree_file)
            
            if len(treefeatures_df) == 0:
                log.write(f"  WARNING: Empty file: {filename}\n")
                continue
            
            stats_row = treefeatures_df.iloc[0].to_dict()
            param_row = params_df[params_df['param_id'] == param_id]
            
            if len(param_row) == 0:
                log.write(f"  WARNING: No parameters for param_id {param_id}\n")
                continue
            
            param_row = param_row.iloc[0]
            
            # Combine metadata, parameters, and tree features
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
    
    log.write(f"\n✓ Loaded {len(all_data)} simulation results\n")
    
    # Create DataFrame
    combined_df = pd.DataFrame(all_data)
    log.write(f"\nCombined: {len(combined_df)} total rows\n")
    log.flush()
    
    # Validate completeness
    if len(combined_df) != expected_sims:
        error_msg = (f"ERROR: Expected {expected_sims} simulations, "
                    f"but got {len(combined_df)}")
        log.write(f"\n{error_msg}\n")
        raise ValueError(error_msg)
    
    # Check for duplicates
    duplicates = combined_df[combined_df.duplicated(subset=['sim_id'], keep=False)]
    if len(duplicates) > 0:
        error_msg = f"ERROR: Found {len(duplicates)} duplicate sim_ids"
        log.write(f"\n{error_msg}\n")
        log.write(f"Duplicates:\n{duplicates}\n")
        raise ValueError(error_msg)
    
    # Sort by param_id and replicate for consistent output
    combined_df = combined_df.sort_values(['param_id', 'replicate']).reset_index(drop=True)
    
    # Get tree feature columns
    metadata_cols = ['sim_id', 'param_id', 'replicate', 'R0', 'recovery_time']
    treefeature_cols = [col for col in combined_df.columns if col not in metadata_cols]
    
    # Save aggregated results
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_df.to_csv(output_file, index=False)
    
    log.write(f"\n" + "=" * 80 + "\n")
    log.write(f"AGGREGATION SUMMARY\n")
    log.write(f"=" * 80 + "\n")
    log.write(f"Dataset: {dataset}\n")
    log.write(f"Total simulations: {len(combined_df)}\n")
    log.write(f"Tree features: {len(treefeature_cols)} features\n")
    log.write(f"Output: {output_file}\n")
    log.write(f"\n✓ Aggregation complete!\n")
    
    # Remove intermediate CSV files
    log.write(f"\n" + "=" * 80 + "\n")
    log.write(f"Cleaning up intermediate files...\n")
    removed_count = 0
    for tree_file in treefeatures_files:
        try:
            os.remove(tree_file)
            removed_count += 1
        except Exception as e:
            log.write(f"  WARNING: Could not remove {tree_file}: {str(e)}\n")
    
    log.write(f"✓ Removed {removed_count} intermediate CSV files\n")

print(f"Aggregated {dataset}: {len(combined_df)} simulations with {len(treefeature_cols)} tree features -> {output_file}")
