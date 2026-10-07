"""
Prepare training and test datasets for NRE model using combined tsfeatures + tree features.

Loads combined tsfeatures+treefeatures and parameters, applies z-score normalization,
and creates numpy arrays suitable for training the Neural Ratio Estimation classifier.

Key points:
- Normalizer is FIT on training data only
- Same normalization applied to test data (no data leakage)
- R0 normalized in log10 space (since it's sampled log-scale)
- Normalization stats saved for inference-time use
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
    """Load combined tsfeatures+treefeatures and prepare arrays."""
    df = pd.read_csv(stats_file)
    metadata_cols = ['sim_id', 'param_id', 'replicate', 'R0', 'recovery_time']
    feature_cols = [col for col in df.columns if col not in metadata_cols]
    features = df[feature_cols].values
    R0 = df['R0'].values
    recovery_time = df['recovery_time'].values
    parameters = np.column_stack([R0, recovery_time])
    param_ids = df['param_id'].values
    sim_ids = df['sim_id'].values
    return features, parameters, param_ids, sim_ids, feature_cols

print("="*70)
print("Preparing Combined TSFeatures + Tree Features Training Data")
print("="*70)
print()

# Load data
print(f"Loading training data: {train_stats_file}")
X_train, theta_train, param_ids_train, sim_ids_train, feature_cols = load_and_prepare_data(train_stats_file)
print(f"  Training samples: {len(X_train)}, Combined features: {X_train.shape[1]}")

print(f"Loading test data: {test_stats_file}")
X_test, theta_test, param_ids_test, sim_ids_test, _ = load_and_prepare_data(test_stats_file)
print(f"  Test samples: {len(X_test)}, Combined features: {X_test.shape[1]}")
print()

# Check for NaN/Inf and filter
print("Checking for invalid values and filtering...")

train_report_path = train_output.replace('.npz', '_masking_report.txt')
X_train, theta_train, param_ids_train, sim_ids_train, masks_train, train_report = create_feature_masks(
    X_train, theta_train, param_ids_train, sim_ids_train,
    feature_names=feature_cols,
    dataset_name="Training",
    output_report_path=train_report_path
)

test_report_path = test_output.replace('.npz', '_masking_report.txt')
X_test, theta_test, param_ids_test, sim_ids_test, masks_test, test_report = create_feature_masks(
    X_test, theta_test, param_ids_test, sim_ids_test,
    feature_names=feature_cols,
    dataset_name="Test",
    output_report_path=test_report_path
)

# Report masking results
if train_report['masking']['n_features_masked'] > 0:
    print(f"  ℹ️  Masked {train_report['masking']['n_features_masked']} feature values in training data")
if test_report['masking']['n_features_masked'] > 0:
    print(f"  ℹ️  Masked {test_report['masking']['n_features_masked']} feature values in test data")

print(f"  ✓ Training data: {len(X_train)} samples (100% retained)")
print(f"  ✓ Test data: {len(X_test)} samples (100% retained)")
print()

# Normalize
print("Normalizing features (z-score)...")
feature_normalizer = ZScoreNormalizer()
X_train_norm = feature_normalizer.fit_transform(X_train, mask=masks_train)
X_test_norm = feature_normalizer.transform(X_test)
print(f"  Train: mean={X_train_norm.mean():.6f}, std={X_train_norm.std():.6f}")
print(f"  Test: mean={X_test_norm.mean():.6f}, std={X_test_norm.std():.6f}")
print()

print("Normalizing parameters...")
# R0 is in log-scale (True), recovery_time is linear (False)
param_normalizer = ParameterNormalizer(is_log_scale=np.array([True, False]))
theta_train_norm = param_normalizer.fit_transform(theta_train)
theta_test_norm = param_normalizer.transform(theta_test)
print(f"  Train: mean={theta_train_norm.mean():.6f}, std={theta_train_norm.std():.6f}")
print(f"  Test: mean={theta_test_norm.mean():.6f}, std={theta_test_norm.std():.6f}")
print()

# Save training data with normalization statistics
print(f"Saving training data: {train_output}")
os.makedirs(os.path.dirname(train_output), exist_ok=True)
np.savez(
    train_output,
    features=X_train_norm,
    parameters=theta_train_norm,
    param_ids=param_ids_train,
    sim_ids=sim_ids_train,
    feature_masks=masks_train,  # Add feature masks
    # Save normalization statistics
    feature_mean=feature_normalizer.mean_,
    feature_std=feature_normalizer.std_,
    param_mean=param_normalizer.mean_,
    param_std=param_normalizer.std_,
    param_is_log=param_normalizer.is_log_scale
)

print(f"Saving test data: {test_output}")
# Save test data (with same normalization stats for reference)
np.savez(
    test_output,
    features=X_test_norm,
    parameters=theta_test_norm,
    param_ids=param_ids_test,
    sim_ids=sim_ids_test,
    feature_masks=masks_test,  # Add feature masks
    # Include normalization stats (same as training)
    feature_mean=feature_normalizer.mean_,
    feature_std=feature_normalizer.std_,
    param_mean=param_normalizer.mean_,
    param_std=param_normalizer.std_,
    param_is_log=param_normalizer.is_log_scale
)

print()
print("="*70)
print("Data Preparation Complete")
print("="*70)
print(f"Combined features: {X_train_norm.shape[1]} (tsfeatures + tree features)")
print(f"Total input dimension: {X_train_norm.shape[1] + theta_train_norm.shape[1]}")
