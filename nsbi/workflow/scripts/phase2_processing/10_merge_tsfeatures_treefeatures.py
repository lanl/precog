"""
Merge tsfeatures with tree features into combined dataset.

Loads both tsfeatures and tree features CSVs, merges them on sim_id,
and creates a single combined dataset with all features.

The combined dataset includes:
- Time series features (tsfeatures): 56 features
- Tree features: 51 features
- Total: 107 features + metadata
"""

import pandas as pd
import sys
from pathlib import Path

# Access snakemake variables
tsfeatures_file = snakemake.input.tsfeatures
treefeatures_file = snakemake.input.treefeatures
output_file = snakemake.output[0]
log_file = snakemake.log[0]

# Setup logging
log_dir = Path(log_file).parent
log_dir.mkdir(parents=True, exist_ok=True)

with open(log_file, 'w') as log:
    log.write("=" * 80 + "\n")
    log.write("Merging TSFeatures + Tree Features\n")
    log.write("=" * 80 + "\n\n")
    
    # Load tsfeatures
    log.write(f"Loading tsfeatures: {tsfeatures_file}\n")
    tsfeatures_df = pd.read_csv(tsfeatures_file)
    log.write(f"  TSFeatures: {len(tsfeatures_df)} rows, {len(tsfeatures_df.columns)} columns\n")
    
    # Get tsfeature columns (exclude metadata)
    metadata_cols = ['sim_id', 'param_id', 'replicate', 'R0', 'recovery_time']
    tsfeature_cols = [col for col in tsfeatures_df.columns if col not in metadata_cols]
    log.write(f"  TSFeature columns: {len(tsfeature_cols)}\n")
    log.write(f"  First 10 tsfeatures: {tsfeature_cols[:10]}\n\n")
    
    # Load tree features
    log.write(f"Loading tree features: {treefeatures_file}\n")
    tree_df = pd.read_csv(treefeatures_file)
    log.write(f"  Tree features: {len(tree_df)} rows, {len(tree_df.columns)} columns\n")
    
    # Get tree feature columns (exclude metadata)
    treefeature_cols = [col for col in tree_df.columns if col not in metadata_cols]
    log.write(f"  Tree feature columns: {len(treefeature_cols)} columns\n")
    log.write(f"  First 10 tree features: {treefeature_cols[:10]}\n\n")
    
    # Check that both have the same sim_ids
    tsfeatures_sims = set(tsfeatures_df['sim_id'])
    tree_sims = set(tree_df['sim_id'])
    
    if tsfeatures_sims != tree_sims:
        missing_in_tree = tsfeatures_sims - tree_sims
        missing_in_tsfeatures = tree_sims - tsfeatures_sims
        
        error_msg = "ERROR: Mismatch in sim_ids between tsfeatures and tree features\n"
        if missing_in_tree:
            error_msg += f"  Missing in tree features: {sorted(list(missing_in_tree))[:10]}...\n"
        if missing_in_tsfeatures:
            error_msg += f"  Missing in tsfeatures: {sorted(list(missing_in_tsfeatures))[:10]}...\n"
        
        log.write(error_msg)
        raise ValueError(error_msg)
    
    log.write(f"✓ Both datasets have matching sim_ids: {len(tsfeatures_sims)} simulations\n\n")
    
    # Merge on sim_id (inner join - should have all rows)
    log.write("Merging datasets on sim_id...\n")
    
    # Extract tsfeatures
    tsfeatures_features = tsfeatures_df[['sim_id'] + tsfeature_cols]
    
    # Extract tree features  
    tree_features = tree_df[['sim_id'] + treefeature_cols]
    
    # Merge features
    combined_features = pd.merge(tsfeatures_features, tree_features, on='sim_id', how='inner')
    
    # Add metadata back from tsfeatures_df
    metadata = tsfeatures_df[['sim_id', 'param_id', 'replicate', 'R0', 'recovery_time']]
    combined_df = pd.merge(metadata, combined_features, on='sim_id', how='inner')
    
    # Reorder columns: metadata first, then tsfeatures, then tree features
    final_cols = ['sim_id', 'param_id', 'replicate', 'R0', 'recovery_time'] + tsfeature_cols + treefeature_cols
    combined_df = combined_df[final_cols]
    
    log.write(f"✓ Merge complete\n")
    log.write(f"  Combined dataset: {len(combined_df)} rows, {len(combined_df.columns)} columns\n")
    log.write(f"  Metadata columns: 5\n")
    log.write(f"  TSFeature columns: {len(tsfeature_cols)}\n")
    log.write(f"  Tree feature columns: {len(treefeature_cols)}\n")
    log.write(f"  Total features: {len(tsfeature_cols) + len(treefeature_cols)}\n\n")
    
    # Validate no missing data
    if combined_df.isnull().any().any():
        null_counts = combined_df.isnull().sum()
        null_cols = null_counts[null_counts > 0]
        log.write(f"WARNING: Missing values detected:\n{null_cols}\n\n")
    else:
        log.write("✓ No missing values\n\n")
    
    # Sort by param_id and replicate for consistency
    combined_df = combined_df.sort_values(['param_id', 'replicate']).reset_index(drop=True)
    
    # Save combined dataset
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_df.to_csv(output_file, index=False)
    
    log.write("=" * 80 + "\n")
    log.write("MERGE SUMMARY\n")
    log.write("=" * 80 + "\n")
    log.write(f"Output: {output_file}\n")
    log.write(f"Total simulations: {len(combined_df)}\n")
    log.write(f"Total features: {len(tsfeature_cols) + len(treefeature_cols)}\n")
    log.write(f"  - TSFeatures: {len(tsfeature_cols)}\n")
    log.write(f"  - Tree features: {len(treefeature_cols)}\n")
    log.write("\n✓ Merge complete!\n")

print(f"Combined tsfeatures + tree features created: {len(combined_df)} simulations, "
      f"{len(tsfeature_cols) + len(treefeature_cols)} features -> {output_file}")
