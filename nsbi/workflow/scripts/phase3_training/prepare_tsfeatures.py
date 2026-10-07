"""
Prepare training and test datasets for NRE model using tsfeatures with z-score normalization.

Loads tsfeatures and parameters, applies z-score normalization,
and creates numpy arrays suitable for training the Neural Ratio Estimation classifier.

Key points:
- Filters out rows with NaN or Inf values before normalization
- Normalizer is FIT on training data only
- Same normalization applied to test data (no data leakage)
- R0 normalized in log10 space (since it's sampled log-scale)
- Normalization stats saved for inference-time use

Inf values handling:
- Inf values corrupt Z-score statistics (mean/std become Inf/NaN)
- Both +Inf and -Inf are detected and filtered
- Rows with Inf in features OR parameters are removed
"""

import numpy as np
import pandas as pd
import sys
import os

# Add utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.normalization import ZScoreNormalizer, ParameterNormalizer, save_normalization_stats, create_feature_masks

# Access snakemake variables
train_stats_file = snakemake.input.train_input
test_stats_file = snakemake.input.test_input
train_output = snakemake.output.train_data
test_output = snakemake.output.test_data

def load_and_prepare_data(stats_file):
    """
    Load tsfeatures and prepare feature and parameter arrays.
    
    Args:
        stats_file: Path to tsfeatures CSV
        
    Returns:
        features: numpy array of shape (N, n_tsfeatures) with time series features
        parameters: numpy array of shape (N, 2) with R0 and recovery_time
        param_ids: numpy array of shape (N,) with parameter IDs
        sim_ids: array of simulation IDs
    """
    df = pd.read_csv(stats_file)
    
    # Extract metadata columns
    metadata_cols = ['sim_id', 'param_id', 'replicate', 'R0', 'recovery_time']
    
    # Extract tsfeature columns (all columns except metadata)
    tsfeature_cols = [col for col in df.columns if col not in metadata_cols]
    
    # Extract features (all tsfeatures)
    features = df[tsfeature_cols].values
    
    # Extract parameters (2 parameters)
    parameters = df[['R0', 'recovery_time']].values
    
    # Extract parameter IDs (critical for safe negative sampling)
    param_ids = df['param_id'].values
    
    # Extract simulation IDs for reference
    sim_ids = df['sim_id'].values
    
    return features, parameters, param_ids, sim_ids, tsfeature_cols

print("="*70)
print("Preparing NRE Training Data (tsfeatures) with Z-Score Normalization")
print("="*70)

# Load training data
print(f"\nLoading training data from {train_stats_file}...")
X_train, theta_train, param_ids_train, sim_ids_train, tsfeature_cols = load_and_prepare_data(train_stats_file)
print(f"  Training samples: {len(X_train)}")
print(f"  Features shape: {X_train.shape}")
print(f"  tsfeatures: {len(tsfeature_cols)}")
print(f"  Parameters shape: {theta_train.shape}")
print(f"  Unique param_ids: {len(np.unique(param_ids_train))}")

# ============================================================================
# DATA VALIDATION: Check and filter NaN/Inf values
# ============================================================================
print("\n" + "="*70)
print("DATA VALIDATION: Checking and Filtering NaN/Inf values")
print("="*70)

# Filter training data
train_report_path = train_output.replace('.npz', '_masking_report.txt')
X_train, theta_train, param_ids_train, sim_ids_train, masks_train, train_report = create_feature_masks(
    X_train, theta_train, param_ids_train, sim_ids_train,
    feature_names=tsfeature_cols,
    dataset_name="Training",
    output_report_path=train_report_path
)

# Print masking summary
if train_report['masking']['n_features_masked'] > 0:
    print(f"\nℹ️  Masked {train_report['masking']['n_features_masked']} feature values in training data")
    print(f"   Samples with masked features: {train_report['masking']['n_samples_with_masked_features']} "
          f"({train_report['masking']['pct_samples_with_masked_features']:.2f}%)")
    
    # Show masked features
    if train_report['masking']['features_with_masking']:
        print(f"   Masked features:")
        for feat in train_report['masking']['features_with_masking'][:10]:  # Show first 10
            print(f"     - {feat['name']}: {feat['n_masked']} values masked ({feat['pct_masked']:.2f}%)")
        if len(train_report['masking']['features_with_masking']) > 10:
            print(f"     ... and {len(train_report['masking']['features_with_masking']) - 10} more")
else:
    print(f"\n✓ Training data: No NaN or Inf values detected")

print(f"✓ Training data after filtering: {len(X_train)} samples")
print(f"✓ Training parameters: No NaN or Inf values detected")

# Load test data
print(f"\nLoading test data from {test_stats_file}...")
X_test, theta_test, param_ids_test, sim_ids_test, _ = load_and_prepare_data(test_stats_file)
print(f"  Test samples: {len(X_test)}")
print(f"  Features shape: {X_test.shape}")
print(f"  Parameters shape: {theta_test.shape}")
print(f"  Unique param_ids: {len(np.unique(param_ids_test))}")

# Filter test data
test_report_path = test_output.replace('.npz', '_masking_report.txt')
X_test, theta_test, param_ids_test, sim_ids_test, masks_test, test_report = create_feature_masks(
    X_test, theta_test, param_ids_test, sim_ids_test,
    feature_names=tsfeature_cols,
    dataset_name="Test",
    output_report_path=test_report_path
)

# Print masking summary
if test_report['masking']['n_features_masked'] > 0:
    print(f"\nℹ️  Masked {test_report['masking']['n_features_masked']} feature values in test data")
    print(f"   Samples with masked features: {test_report['masking']['n_samples_with_masked_features']} "
          f"({test_report['masking']['pct_samples_with_masked_features']:.2f}%)")
    
    # Show masked features
    if test_report['masking']['features_with_masking']:
        print(f"   Masked features:")
        for feat in test_report['masking']['features_with_masking'][:10]:  # Show first 10
            print(f"     - {feat['name']}: {feat['n_masked']} values masked ({feat['pct_masked']:.2f}%)")
        if len(test_report['masking']['features_with_masking']) > 10:
            print(f"     ... and {len(test_report['masking']['features_with_masking']) - 10} more")
else:
    print(f"\n✓ Test data: No masked features (all complete)")

print(f"✓ Test data: {len(X_test)} samples (100% retained)")



# Print raw statistics
print("\n" + "="*70)
print("RAW Data Statistics (before normalization)")
print("="*70)
print(f"\nFeatures (train) - showing first 5 tsfeatures:")
print(f"  Mean: {X_train.mean(axis=0)[:5]}")
print(f"  Std:  {X_train.std(axis=0)[:5]}")
print(f"\nParameters (train):")
print(f"  R0            - mean: {theta_train[:, 0].mean():.4f}, std: {theta_train[:, 0].std():.4f}")
print(f"  recovery_time - mean: {theta_train[:, 1].mean():.4f}, std: {theta_train[:, 1].std():.4f}")
print(f"  log10(R0)     - mean: {np.log10(theta_train[:, 0]).mean():.4f}, std: {np.log10(theta_train[:, 0]).std():.4f}")

# Initialize normalizers
print("\n" + "="*70)
print("Fitting Normalizers on Training Data")
print("="*70)

# Feature normalizer (standard z-score) - fit with masks
feature_normalizer = ZScoreNormalizer()
feature_normalizer.fit(X_train, mask=masks_train)
print(f"\n✓ Feature normalizer fitted")
print(f"  Feature means (first 5): {feature_normalizer.mean_[:5]}")
print(f"  Feature stds (first 5):  {feature_normalizer.std_[:5]}")

# Parameter normalizer (R0 in log10 space, recovery_time in linear space)
param_normalizer = ParameterNormalizer(is_log_scale=np.array([True, False]))
param_normalizer.fit(theta_train)
print(f"\n✓ Parameter normalizer fitted")
print(f"  Param means (in transform space): {param_normalizer.mean_}")
print(f"  Param stds (in transform space):  {param_normalizer.std_}")
print(f"  Transform space: [log10(R0), recovery_time]")

# Apply normalization to training data
print("\n" + "="*70)
print("Applying Normalization")
print("="*70)

X_train_normalized = feature_normalizer.transform(X_train)
theta_train_normalized = param_normalizer.transform(theta_train)
print(f"\n✓ Training data normalized")

# Apply same normalization to test data (critical: use training stats!)
X_test_normalized = feature_normalizer.transform(X_test)
theta_test_normalized = param_normalizer.transform(theta_test)
print(f"✓ Test data normalized (using training statistics)")

# Verify normalization
print("\n" + "="*70)
print("NORMALIZED Data Statistics (should be ~N(0,1) for training)")
print("="*70)
print(f"\nFeatures (train) - showing first 5:")
print(f"  Mean: {X_train_normalized.mean(axis=0)[:5]}")
print(f"  Std:  {X_train_normalized.std(axis=0)[:5]}")
print(f"\nParameters (train):")
print(f"  Mean: {theta_train_normalized.mean(axis=0)}")
print(f"  Std:  {theta_train_normalized.std(axis=0)}")
print(f"\nFeatures (test) - showing first 5:")
print(f"  Mean: {X_test_normalized.mean(axis=0)[:5]}")
print(f"  Std:  {X_test_normalized.std(axis=0)[:5]}")
print(f"\nParameters (test):")
print(f"  Mean: {theta_test_normalized.mean(axis=0)}")
print(f"  Std:  {theta_test_normalized.std(axis=0)}")

# Save training data with normalization stats
np.savez(
    train_output,
    features=X_train_normalized,
    parameters=theta_train_normalized,
    param_ids=param_ids_train,
    sim_ids=sim_ids_train,
    feature_masks=masks_train,  # Add feature masks
    feature_mean=feature_normalizer.mean_,
    feature_std=feature_normalizer.std_,
    param_mean=param_normalizer.mean_,
    param_std=param_normalizer.std_,
    param_is_log=param_normalizer.is_log_scale
)
print(f"\n✓ Training data saved to: {train_output}")

# Save test data (with same normalization stats for reference)
np.savez(
    test_output,
    features=X_test_normalized,
    parameters=theta_test_normalized,
    param_ids=param_ids_test,
    sim_ids=sim_ids_test,
    feature_masks=masks_test,  # Add feature masks
    feature_mean=feature_normalizer.mean_,
    feature_std=feature_normalizer.std_,
    param_mean=param_normalizer.mean_,
    param_std=param_normalizer.std_,
    param_is_log=param_normalizer.is_log_scale
)
print(f"✓ Test data saved to: {test_output}")

print("\n" + "="*70)
print("Data Preparation Complete")
print("="*70)
print("\nNormalization Details:")
print("  - Features: Standard z-score normalization")
print("  - R0: Z-score normalization in log10 space")
print("  - recovery_time: Standard z-score normalization")
print("  - Statistics fitted on training data ONLY")
print("  - Same statistics applied to test data")
print(f"\ntsfeatures: {tsfeature_cols[:10]}... (showing first 10)")
print("\nAll data is now ready for NRE training!")
