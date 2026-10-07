# Masking utilities - functions to create masks for NA/Inf values
import numpy as np


def print_mask_statistics(masks, dataset_name="Dataset"):
    """
    Print summary statistics about mask coverage.
    
    Args:
        masks: Binary mask array of shape (n_samples, n_features)
               where 1=valid, 0=masked
        dataset_name: Name to display in output
    """
    n_samples, n_features = masks.shape
    
    # Overall statistics
    total_values = n_samples * n_features
    valid_values = int(masks.sum())
    masked_values = total_values - valid_values
    
    print(f"\n📊 {dataset_name} Mask Statistics:")
    print(f"  Total values: {total_values:,}")
    print(f"  Valid values: {valid_values:,} ({valid_values/total_values*100:.2f}%)")
    print(f"  Masked values: {masked_values:,} ({masked_values/total_values*100:.2f}%)")
    
    # Per-sample statistics
    valid_per_sample = masks.sum(axis=1)
    print(f"  Valid features per sample: {valid_per_sample.mean():.1f} ± {valid_per_sample.std():.1f} (max={n_features})")
    
    # Samples with any masking
    samples_with_masking = int((masks.sum(axis=1) < n_features).sum())
    print(f"  Samples with masked features: {samples_with_masking} ({samples_with_masking/n_samples*100:.2f}%)")
    
    # Fully masked samples (if any)
    fully_masked = int((masks.sum(axis=1) == 0).sum())
    if fully_masked > 0:
        print(f"  ⚠️  Fully masked samples: {fully_masked} (all features invalid)")


def create_feature_masks(features, parameters, param_ids, sim_ids, feature_names=None, dataset_name="data", output_report_path=None):
    """Create binary masks for features with NaN or Inf values instead of removing rows.
    
    Args:
        features: Feature array (may contain NaN/Inf)
        parameters: Parameter array (must NOT contain NaN/Inf - will raise error)
        param_ids: Parameter IDs
        sim_ids: Simulation IDs
        feature_names: Optional list of feature names for reporting
        dataset_name: Name for reporting
        output_report_path: Optional path to save report
    
    Returns:
        features: Cleaned features (invalid values replaced with 0)
        parameters: Original parameters (unchanged)
        param_ids: Original param_ids (unchanged)
        sim_ids: Original sim_ids (unchanged)
        masks: Binary masks (1=valid, 0=invalid)
        report: Dictionary with masking statistics
    """
    from datetime import datetime
    import numpy as np
    
    n_samples = len(features)
    if len(parameters) != n_samples or len(param_ids) != n_samples or len(sim_ids) != n_samples:
        raise ValueError("Array length mismatch")
    
    n_features = features.shape[1] if features.ndim > 1 else 1
    n_params = parameters.shape[1] if parameters.ndim > 1 else 1
    n_unique_param_ids = len(np.unique(param_ids))
    
    # Check for invalid parameters - this should NEVER happen
    nan_in_params = np.isnan(parameters).any()
    inf_in_params = np.isinf(parameters).any()
    
    if nan_in_params or inf_in_params:
        n_nan = np.isnan(parameters).sum()
        n_inf = np.isinf(parameters).sum()
        affected_rows = np.where(np.isnan(parameters).any(axis=1) | np.isinf(parameters).any(axis=1))[0]
        raise ValueError(
            f"Invalid parameters detected in {dataset_name}! Parameters should never contain NA/Inf.\\n"
            f"  NaN count: {n_nan}\\n"
            f"  Inf count: {n_inf}\\n"
            f"  Affected rows: {len(affected_rows)} (first 10: {affected_rows[:10].tolist()})\\n"
            f"This indicates a problem in the simulation/data generation pipeline."
        )
    
    # Copy features to avoid modifying input
    features = features.copy()
    
    # Detect invalid values (NaN or Inf) in features
    invalid_mask = np.isnan(features) | np.isinf(features)
    
    # Create masks (1 = valid, 0 = invalid)
    feature_masks = (~invalid_mask).astype(np.float32)
    
    # Fill invalid positions with 0
    features = np.where(invalid_mask, 0.0, features)
    
    # Compute masking statistics
    n_features_masked = (feature_masks == 0).sum()
    n_samples_with_masked_features = (feature_masks.sum(axis=1) < n_features).sum()
    masked_per_feature = (feature_masks == 0).sum(axis=0)
    
    # Per-feature masking details
    features_with_masking = []
    for feat_idx in range(n_features):
        n_masked = int(masked_per_feature[feat_idx])
        if n_masked > 0:
            feat_name = feature_names[feat_idx] if feature_names is not None else f"feature_{feat_idx}"
            pct_masked = (n_masked / n_samples * 100) if n_samples > 0 else 0
            features_with_masking.append({
                'index': int(feat_idx),
                'name': feat_name,
                'n_masked': n_masked,
                'pct_masked': pct_masked
            })
    
    # Build report
    report = {
        'dataset_name': dataset_name,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'data': {
            'n_samples': n_samples,
            'n_features': n_features,
            'n_params': n_params,
            'n_unique_param_ids': n_unique_param_ids
        },
        'masking': {
            'n_features_masked': int(n_features_masked),
            'n_samples_with_masked_features': int(n_samples_with_masked_features),
            'pct_samples_with_masked_features': (n_samples_with_masked_features / n_samples * 100) if n_samples > 0 else 0,
            'features_with_masking': features_with_masking
        },
        'output_shapes': {
            'features': tuple(features.shape),
            'parameters': tuple(parameters.shape),
            'masks': tuple(feature_masks.shape)
        }
    }
    
    # Save report if requested
    if output_report_path:
        save_masking_report(report, output_report_path)
    
    return features, parameters, param_ids, sim_ids, feature_masks, report



def save_masking_report(report, filepath):
    """Save feature masking report to a text file."""
    import os
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    with open(filepath, 'w') as f:
        f.write("=" * 70 + "\\n")
        f.write(f"FEATURE MASKING REPORT - {report['dataset_name']}\\n")
        f.write("=" * 70 + "\\n")
        f.write(f"Generated: {report['timestamp']}\\n\\n")
        f.write("STRATEGY: Mask features with NA/Inf (kept in dataset)\\n")
        f.write("          Error if parameters contain NA/Inf (should never happen)\\n\\n")
        
        f.write("-" * 70 + "\\n")
        f.write("DATA\\n")
        f.write("-" * 70 + "\\n")
        data = report['data']
        f.write(f"Samples: {data['n_samples']}\\n")
        f.write(f"Features: {data['n_features']}\\n")
        f.write(f"Parameters: {data['n_params']}\\n")
        f.write(f"Unique param_ids: {data['n_unique_param_ids']}\\n\\n")
        
        f.write("-" * 70 + "\\n")
        f.write("FEATURE MASKING (NA/Inf values in features)\\n")
        f.write("-" * 70 + "\\n")
        masking = report['masking']
        f.write(f"Total feature values masked: {masking['n_features_masked']}\\n")
        f.write(f"Samples with ≥1 masked feature: {masking['n_samples_with_masked_features']} ({masking['pct_samples_with_masked_features']:.2f}%)\\n")
        f.write(f"Features affected: {len(masking['features_with_masking'])}\\n\\n")
        
        if masking['features_with_masking']:
            f.write("Features with masked values:\\n")
            for feat_info in masking['features_with_masking']:
                f.write(f"  [{feat_info['index']}] {feat_info['name']}: {feat_info['n_masked']} masked ({feat_info['pct_masked']:.2f}%)\\n")
            f.write("\\n")
        else:
            f.write("No masked features - all features are complete!\\n\\n")
        
        f.write("-" * 70 + "\\n")
        f.write("OUTPUT SHAPES\\n")
        f.write("-" * 70 + "\\n")
        shapes = report['output_shapes']
        f.write(f"Features: {shapes['features']}\\n")
        f.write(f"Parameters: {shapes['parameters']}\\n")
        f.write(f"Masks: {shapes['masks']}\\n\\n")
        
        f.write("=" * 70 + "\\n")
        f.write("SUMMARY\\n")
        f.write("=" * 70 + "\\n")
        f.write(f"Approach: Masking (keep all simulations, mask invalid features)\\n")
        f.write(f"Features masked: {masking['n_features_masked']} values across {len(masking['features_with_masking'])} features\\n")
        f.write(f"Simulations retained: {data['n_samples']} (100%)\\n")
        f.write(f"Result: All simulations retained with masks for invalid feature values\\n")
        f.write("=" * 70 + "\\n")
    
    print(f"✓ Masking report saved to: {filepath}")
