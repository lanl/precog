"""
Z-Score Normalization Utilities for Neural SBI

Provides standardization (z-score normalization) for features and parameters.
Handles special cases like log-scale parameters (R0) and numerical stability.

Key Design:
-----------
1. Fit normalization statistics on TRAINING DATA ONLY
2. Apply same statistics to train, validation, and test sets
3. Handle R0 in log10 space (since it's sampled on log-scale)
4. Save statistics for inference-time normalization
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any

# Import masking functions
from .masking_funcs import create_feature_masks, save_masking_report



class ZScoreNormalizer:
    """
    Z-score (standardization) normalizer with support for log-scale parameters.
    
    Computes: z = (x - mean) / std
    
    Attributes:
        mean_: Array of means (fitted from training data)
        std_: Array of standard deviations (fitted from training data)
        epsilon: Small constant for numerical stability
        fitted: Whether the normalizer has been fitted
    """
    
    def __init__(self, epsilon: float = 1e-8):
        """
        Initialize normalizer.
        
        Args:
            epsilon: Small constant added to std to prevent division by zero
        """
        self.epsilon = epsilon
        self.mean_ = None
        self.std_ = None
        self.fitted = False
    
    def fit(self, X: np.ndarray, mask: Optional[np.ndarray] = None) -> 'ZScoreNormalizer':
        """
        Fit normalizer to training data.
        
        Args:
            X: Training data of shape (n_samples, n_features)
            mask: Optional binary mask of shape (n_samples, n_features)
                  where 1=valid, 0=invalid. If provided, computes statistics
                  only on valid entries.
        
        Returns:
            self (for chaining)
        """
        if X.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {X.shape}")
        
        if mask is not None:
            if mask.shape != X.shape:
                raise ValueError(f"Mask shape {mask.shape} doesn't match X shape {X.shape}")
            # Compute mean and std only on valid (masked) entries
            self.mean_ = np.zeros(X.shape[1])
            self.std_ = np.zeros(X.shape[1])
            
            for i in range(X.shape[1]):
                valid_values = X[:, i][mask[:, i] == 1]
                if len(valid_values) > 0:
                    self.mean_[i] = np.mean(valid_values)
                    self.std_[i] = np.std(valid_values)
                else:
                    # No valid values - use defaults
                    self.mean_[i] = 0.0
                    self.std_[i] = 1.0
        else:
            # Standard case: compute over all data
            self.mean_ = np.mean(X, axis=0)
            self.std_ = np.std(X, axis=0)
        
        # Handle zero std (constant features)
        self.std_ = np.where(self.std_ < self.epsilon, 1.0, self.std_)
        
        self.fitted = True
        return self
    
    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Apply z-score normalization using fitted statistics.
        
        Args:
            X: Data to normalize of shape (n_samples, n_features)
        
        Returns:
            Normalized data with same shape as input
        """
        if not self.fitted:
            raise RuntimeError("Normalizer must be fitted before transform")
        
        if X.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {X.shape}")
        
        if X.shape[1] != len(self.mean_):
            raise ValueError(
                f"X has {X.shape[1]} features but normalizer was fitted on {len(self.mean_)} features"
            )
        
        return (X - self.mean_) / self.std_
    
    def fit_transform(self, X: np.ndarray, mask: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Fit normalizer and transform data in one step.
        
        Args:
            X: Training data of shape (n_samples, n_features)
            mask: Optional binary mask
        
        Returns:
            Normalized data
        """
        return self.fit(X, mask).transform(X)
    
    def inverse_transform(self, X_normalized: np.ndarray) -> np.ndarray:
        """
        Convert normalized data back to original scale.
        
        Args:
            X_normalized: Normalized data of shape (n_samples, n_features)
        
        Returns:
            Data in original scale
        """
        if not self.fitted:
            raise RuntimeError("Normalizer must be fitted before inverse_transform")
        
        return X_normalized * self.std_ + self.mean_
    
    def get_statistics(self) -> Dict[str, np.ndarray]:
        """
        Get fitted normalization statistics.
        
        Returns:
            Dictionary with 'mean' and 'std' arrays
        """
        if not self.fitted:
            raise RuntimeError("Normalizer must be fitted before getting statistics")
        
        return {
            'mean': self.mean_.copy(),
            'std': self.std_.copy()
        }


class ParameterNormalizer:
    """
    Specialized normalizer for parameters with log-scale handling.
    
    Handles the case where some parameters (like R0) are sampled on a log scale.
    For log-scale parameters, normalization is applied in log-space.
    
    Example:
        R0 ~ Uniform(log10(1.5), log10(4.0))
        
        We normalize: z = (log10(R0) - mean(log10(R0))) / std(log10(R0))
        Not:          z = (R0 - mean(R0)) / std(R0)
    """
    
    def __init__(self, is_log_scale: np.ndarray, epsilon: float = 1e-8):
        """
        Initialize parameter normalizer.
        
        Args:
            is_log_scale: Boolean array of shape (n_params,) indicating which
                          parameters are sampled on log-scale.
                          Example: [True, False] for [R0, recovery_time]
            epsilon: Small constant for numerical stability
        """
        self.is_log_scale = np.asarray(is_log_scale, dtype=bool)
        self.epsilon = epsilon
        self.mean_ = None
        self.std_ = None
        self.fitted = False
    
    def _to_transform_space(self, params: np.ndarray) -> np.ndarray:
        """Convert parameters to space where normalization is applied."""
        params_transformed = params.copy()
        params_transformed[:, self.is_log_scale] = np.log10(
            params_transformed[:, self.is_log_scale]
        )
        return params_transformed
    
    def _from_transform_space(self, params_transformed: np.ndarray) -> np.ndarray:
        """Convert parameters back from transform space."""
        params = params_transformed.copy()
        params[:, self.is_log_scale] = 10 ** params[:, self.is_log_scale]
        return params
    
    def fit(self, params: np.ndarray) -> 'ParameterNormalizer':
        """
        Fit normalizer to training parameters.
        
        Args:
            params: Parameter array of shape (n_samples, n_params)
                    Values should be in ORIGINAL scale (e.g., R0 in [1.5, 4.0])
        
        Returns:
            self (for chaining)
        """
        if params.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {params.shape}")
        
        if params.shape[1] != len(self.is_log_scale):
            raise ValueError(
                f"params has {params.shape[1]} parameters but is_log_scale has {len(self.is_log_scale)}"
            )
        
        # Transform to normalization space (apply log10 where needed)
        params_transformed = self._to_transform_space(params)
        
        # Compute statistics in transform space
        self.mean_ = np.mean(params_transformed, axis=0)
        self.std_ = np.std(params_transformed, axis=0)
        
        # Handle zero std
        self.std_ = np.where(self.std_ < self.epsilon, 1.0, self.std_)
        
        self.fitted = True
        return self
    
    def transform(self, params: np.ndarray) -> np.ndarray:
        """
        Normalize parameters using fitted statistics.
        
        Args:
            params: Parameters in ORIGINAL scale of shape (n_samples, n_params)
        
        Returns:
            Normalized parameters (in transform space)
        """
        if not self.fitted:
            raise RuntimeError("Normalizer must be fitted before transform")
        
        # Transform to normalization space
        params_transformed = self._to_transform_space(params)
        
        # Apply z-score normalization
        return (params_transformed - self.mean_) / self.std_
    
    def fit_transform(self, params: np.ndarray) -> np.ndarray:
        """
        Fit normalizer and transform parameters in one step.
        
        Args:
            params: Parameter array in ORIGINAL scale
        
        Returns:
            Normalized parameters
        """
        return self.fit(params).transform(params)
    
    def inverse_transform(self, params_normalized: np.ndarray) -> np.ndarray:
        """
        Convert normalized parameters back to original scale.
        
        Args:
            params_normalized: Normalized parameters
        
        Returns:
            Parameters in original scale
        """
        if not self.fitted:
            raise RuntimeError("Normalizer must be fitted before inverse_transform")
        
        # Denormalize in transform space
        params_transformed = params_normalized * self.std_ + self.mean_
        
        # Convert back to original scale (reverse log10 where needed)
        return self._from_transform_space(params_transformed)
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get fitted normalization statistics.
        
        Returns:
            Dictionary with statistics and metadata
        """
        if not self.fitted:
            raise RuntimeError("Normalizer must be fitted before getting statistics")
        
        return {
            'mean': self.mean_.copy(),
            'std': self.std_.copy(),
            'is_log_scale': self.is_log_scale.copy()
        }


def save_normalization_stats(
    filepath: str,
    feature_stats: Dict[str, np.ndarray],
    param_stats: Dict[str, Any]
) -> None:
    """
    Save normalization statistics to file.
    
    Args:
        filepath: Path to save .npz file
        feature_stats: Dictionary from feature normalizer.get_statistics()
        param_stats: Dictionary from parameter normalizer.get_statistics()
    """
    np.savez(
        filepath,
        feature_mean=feature_stats['mean'],
        feature_std=feature_stats['std'],
        param_mean=param_stats['mean'],
        param_std=param_stats['std'],
        param_is_log=param_stats['is_log_scale']
    )


def load_normalization_stats(filepath: str) -> Tuple[Dict[str, np.ndarray], Dict[str, Any]]:
    """
    Load normalization statistics from file.
    
    Args:
        filepath: Path to .npz file
    
    Returns:
        Tuple of (feature_stats, param_stats) dictionaries
    """
    data = np.load(filepath)
    
    feature_stats = {
        'mean': data['feature_mean'],
        'std': data['feature_std']
    }
    
    param_stats = {
        'mean': data['param_mean'],
        'std': data['param_std'],
        'is_log_scale': data['param_is_log']
    }
    
    return feature_stats, param_stats



def filter_nan_rows(
    features: np.ndarray,
    parameters: np.ndarray,
    param_ids: np.ndarray,
    sim_ids: np.ndarray,
    feature_names: Optional[list] = None,
    dataset_name: str = "data",
    output_report_path: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Remove rows containing NaN or Inf values in features or parameters.
    
    This function filters out any simulation (row) that has NaN or Inf (±∞) 
    in either features or parameters. All arrays remain aligned after filtering.
    
    Inf values are problematic for Z-score normalization as they:
    - Corrupt mean/std statistics during fitting
    - Propagate through the normalization formula
    - Result in NaN values after normalization
    
    Args:
        features: Feature array of shape (N, n_features)
        parameters: Parameter array of shape (N, n_params)
        param_ids: Parameter IDs array of shape (N,)
        sim_ids: Simulation IDs array of shape (N,)
        feature_names: Optional list of feature names for detailed reporting
        dataset_name: Name of dataset for reporting (e.g., "Training", "Test")
        output_report_path: Optional path to save detailed report as .txt file
    
    Returns:
        Tuple of (filtered_features, filtered_parameters, filtered_param_ids,
                  filtered_sim_ids, report_dict)
        
        All returned arrays have the same length and maintain row correspondence.
        report_dict contains statistics about the filtering operation.
    """
    from datetime import datetime
    
    # Input validation
    n_samples = len(features)
    if len(parameters) != n_samples or len(param_ids) != n_samples or len(sim_ids) != n_samples:
        raise ValueError(
            f"Array length mismatch: features={len(features)}, "
            f"parameters={len(parameters)}, param_ids={len(param_ids)}, "
            f"sim_ids={len(sim_ids)}"
        )
    
    # Store original statistics
    original_count = n_samples
    n_features = features.shape[1] if features.ndim > 1 else 1
    n_params = parameters.shape[1] if parameters.ndim > 1 else 1
    original_unique_param_ids = len(np.unique(param_ids))
    
    # Detect NaN values
    nan_in_features_mask = np.isnan(features).any(axis=1)
    nan_in_params_mask = np.isnan(parameters).any(axis=1)
    
    # Detect Inf values (both +Inf and -Inf)
    inf_in_features_mask = np.isinf(features).any(axis=1)
    inf_in_params_mask = np.isinf(parameters).any(axis=1)
    
    # Combined mask: any row with NaN OR Inf in features OR parameters
    invalid_mask = nan_in_features_mask | nan_in_params_mask | inf_in_features_mask | inf_in_params_mask
    
    # Count affected rows by type
    n_nan_features = nan_in_features_mask.sum()
    n_nan_params = nan_in_params_mask.sum()
    n_inf_features = inf_in_features_mask.sum()
    n_inf_params = inf_in_params_mask.sum()
    n_rows_with_invalid = invalid_mask.sum()
    
    # Detailed feature-level statistics
    nan_per_feature = np.isnan(features).sum(axis=0)
    inf_per_feature = np.isinf(features).sum(axis=0)
    features_with_nan_idx = np.where(nan_per_feature > 0)[0]
    features_with_inf_idx = np.where(inf_per_feature > 0)[0]
    
    # Detailed parameter-level statistics
    nan_per_param = np.isnan(parameters).sum(axis=0)
    inf_per_param = np.isinf(parameters).sum(axis=0)
    params_with_nan_idx = np.where(nan_per_param > 0)[0]
    params_with_inf_idx = np.where(inf_per_param > 0)[0]
    
    # Get affected IDs
    affected_sim_ids = sim_ids[invalid_mask].tolist() if n_rows_with_invalid > 0 else []
    affected_param_ids = np.unique(param_ids[invalid_mask]).tolist() if n_rows_with_invalid > 0 else []
    
    # Filter data: keep only rows WITHOUT NaN or Inf
    valid_mask = ~invalid_mask
    filtered_features = features[valid_mask]
    filtered_parameters = parameters[valid_mask]
    filtered_param_ids = param_ids[valid_mask]
    filtered_sim_ids = sim_ids[valid_mask]
    
    # Post-filtering statistics
    final_count = len(filtered_features)
    n_removed = original_count - final_count
    pct_removed = (n_removed / original_count * 100) if original_count > 0 else 0
    pct_retained = 100 - pct_removed
    final_unique_param_ids = len(np.unique(filtered_param_ids))
    
    # Verify no NaN or Inf remains
    has_invalid_after = (
        np.any(np.isnan(filtered_features)) or np.any(np.isnan(filtered_parameters)) or
        np.any(np.isinf(filtered_features)) or np.any(np.isinf(filtered_parameters))
    )
    
    # Verify array alignment
    arrays_aligned = (
        len(filtered_features) == len(filtered_parameters) == 
        len(filtered_param_ids) == len(filtered_sim_ids)
    )


    
    # Build report dictionary
    report = {
        'dataset_name': dataset_name,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'original': {
            'n_samples': original_count,
            'n_features': n_features,
            'n_params': n_params,
            'n_unique_param_ids': original_unique_param_ids
        },
        'nan_detection': {
            'n_rows_with_nan_features': n_nan_features,
            'n_rows_with_nan_params': n_nan_params,
            'n_rows_with_any_nan': (nan_in_features_mask | nan_in_params_mask).sum(),
            'features_affected': [],
            'params_affected': [],
            'affected_sim_ids': affected_sim_ids,
            'affected_param_ids': affected_param_ids
        },
        'inf_detection': {
            'n_rows_with_inf_features': n_inf_features,
            'n_rows_with_inf_params': n_inf_params,
            'n_rows_with_any_inf': (inf_in_features_mask | inf_in_params_mask).sum(),
            'features_affected': [],
            'params_affected': []
        },
        'filtered': {
            'n_removed': n_removed,
            'n_retained': final_count,
            'pct_removed': pct_removed,
            'pct_retained': pct_retained,
            'n_unique_param_ids': final_unique_param_ids,
            'final_shape_features': filtered_features.shape,
            'final_shape_params': filtered_parameters.shape
        },
        'verification': {
            'no_invalid_remaining': not has_invalid_after,
            'arrays_aligned': arrays_aligned
        }
    }
    
    # Add detailed feature-level NaN information
    if len(features_with_nan_idx) > 0:
        for idx in features_with_nan_idx:
            feature_name = feature_names[idx] if feature_names and idx < len(feature_names) else f"Feature_{idx}"
            n_nan = nan_per_feature[idx]
            pct_nan = (n_nan / original_count * 100) if original_count > 0 else 0
            report['nan_detection']['features_affected'].append({
                'index': int(idx),
                'name': feature_name,
                'n_nan': int(n_nan),
                'pct_nan': pct_nan
            })
    
    # Add detailed feature-level Inf information
    if len(features_with_inf_idx) > 0:
        for idx in features_with_inf_idx:
            feature_name = feature_names[idx] if feature_names and idx < len(feature_names) else f"Feature_{idx}"
            n_inf = inf_per_feature[idx]
            pct_inf = (n_inf / original_count * 100) if original_count > 0 else 0
            report['inf_detection']['features_affected'].append({
                'index': int(idx),
                'name': feature_name,
                'n_inf': int(n_inf),
                'pct_inf': pct_inf
            })
    
    # Add detailed parameter-level NaN information
    if len(params_with_nan_idx) > 0:
        param_names = ['R0', 'recovery_time'] if n_params == 2 else [f"Param_{i}" for i in range(n_params)]
        for idx in params_with_nan_idx:
            param_name = param_names[idx] if idx < len(param_names) else f"Param_{idx}"
            n_nan = nan_per_param[idx]
            pct_nan = (n_nan / original_count * 100) if original_count > 0 else 0
            report['nan_detection']['params_affected'].append({
                'index': int(idx),
                'name': param_name,
                'n_nan': int(n_nan),
                'pct_nan': pct_nan
            })
    
    # Add detailed parameter-level Inf information
    if len(params_with_inf_idx) > 0:
        param_names = ['R0', 'recovery_time'] if n_params == 2 else [f"Param_{i}" for i in range(n_params)]
        for idx in params_with_inf_idx:
            param_name = param_names[idx] if idx < len(param_names) else f"Param_{idx}"
            n_inf = inf_per_param[idx]
            pct_inf = (n_inf / original_count * 100) if original_count > 0 else 0
            report['inf_detection']['params_affected'].append({
                'index': int(idx),
                'name': param_name,
                'n_inf': int(n_inf),
                'pct_inf': pct_inf
            })
    
    # Save report to file if requested
    if output_report_path:
        save_nan_removal_report(report, output_report_path)
    
    return filtered_features, filtered_parameters, filtered_param_ids, filtered_sim_ids, report



def save_nan_removal_report(report: Dict[str, Any], filepath: str) -> None:
    """
    Save NaN removal report to a text file.
    
    Args:
        report: Report dictionary from filter_nan_rows()
        filepath: Path to save the report file
    """
    import os
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    with open(filepath, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"NaN/Inf REMOVAL REPORT - {report['dataset_name']}\n")
        f.write("=" * 70 + "\n")
        f.write(f"Generated: {report['timestamp']}\n")
        f.write("\n")
        
        # BEFORE FILTERING
        f.write("-" * 70 + "\n")
        f.write("BEFORE FILTERING\n")
        f.write("-" * 70 + "\n")
        orig = report['original']
        f.write(f"Total samples: {orig['n_samples']}\n")
        f.write(f"Total features: {orig['n_features']}\n")
        f.write(f"Total parameters: {orig['n_params']}\n")
        f.write(f"Unique param_ids: {orig['n_unique_param_ids']}\n")
        f.write("\n")
        
        # NaN DETECTION
        f.write("-" * 70 + "\n")
        f.write("NaN DETECTION\n")
        f.write("-" * 70 + "\n")
        nan_det = report['nan_detection']
        
        # Features with NaN
        if len(nan_det['features_affected']) > 0:
            f.write("Features with NaN values:\n")
            for feat in nan_det['features_affected']:
                f.write(f"  Feature {feat['index']} ({feat['name']}): "
                       f"{feat['n_nan']} samples ({feat['pct_nan']:.3f}%)\n")
        else:
            f.write("Features with NaN values: None\n")
        f.write("\n")
        
        # Parameters with NaN
        if len(nan_det['params_affected']) > 0:
            f.write("Parameters with NaN values:\n")
            for param in nan_det['params_affected']:
                f.write(f"  {param['name']}: {param['n_nan']} samples ({param['pct_nan']:.3f}%)\n")
        else:
            f.write("Parameters with NaN values: None\n")
        f.write("\n")
        
        # Summary of affected rows by NaN
        f.write("Rows affected by NaN:\n")
        f.write(f"  In features: {nan_det['n_rows_with_nan_features']} samples\n")
        f.write(f"  In parameters: {nan_det['n_rows_with_nan_params']} samples\n")
        f.write(f"  Total unique rows with NaN: {nan_det['n_rows_with_any_nan']} samples\n")
        f.write("\n")
        
        # Inf DETECTION
        f.write("-" * 70 + "\n")
        f.write("Inf DETECTION\n")
        f.write("-" * 70 + "\n")
        inf_det = report['inf_detection']
        
        # Features with Inf
        if len(inf_det['features_affected']) > 0:
            f.write("Features with Inf values:\n")
            for feat in inf_det['features_affected']:
                f.write(f"  Feature {feat['index']} ({feat['name']}): "
                       f"{feat['n_inf']} samples ({feat['pct_inf']:.3f}%)\n")
        else:
            f.write("Features with Inf values: None\n")
        f.write("\n")
        
        # Parameters with Inf
        if len(inf_det['params_affected']) > 0:
            f.write("Parameters with Inf values:\n")
            for param in inf_det['params_affected']:
                f.write(f"  {param['name']}: {param['n_inf']} samples ({param['pct_inf']:.3f}%)\n")
        else:
            f.write("Parameters with Inf values: None\n")
        f.write("\n")
        
        # Summary of affected rows by Inf
        f.write("Rows affected by Inf:\n")
        f.write(f"  In features: {inf_det['n_rows_with_inf_features']} samples\n")
        f.write(f"  In parameters: {inf_det['n_rows_with_inf_params']} samples\n")
        f.write(f"  Total unique rows with Inf: {inf_det['n_rows_with_any_inf']} samples\n")
        f.write("\n")
        
        # Show affected IDs (limit to first 20 to avoid huge reports)
        if len(nan_det['affected_sim_ids']) > 0:
            n_show = min(20, len(nan_det['affected_sim_ids']))
            f.write(f"Affected sim_ids (showing first {n_show}): {nan_det['affected_sim_ids'][:n_show]}\n")
            if len(nan_det['affected_sim_ids']) > 20:
                f.write(f"  ... and {len(nan_det['affected_sim_ids']) - 20} more\n")
        
        if len(nan_det['affected_param_ids']) > 0:
            n_show = min(20, len(nan_det['affected_param_ids']))
            f.write(f"Affected param_ids (showing first {n_show}): {nan_det['affected_param_ids'][:n_show]}\n")
            if len(nan_det['affected_param_ids']) > 20:
                f.write(f"  ... and {len(nan_det['affected_param_ids']) - 20} more\n")
        f.write("\n")
        
        # AFTER FILTERING
        f.write("-" * 70 + "\n")
        f.write("AFTER FILTERING\n")
        f.write("-" * 70 + "\n")
        filt = report['filtered']
        f.write(f"Samples removed: {filt['n_removed']} ({filt['pct_removed']:.3f}%)\n")
        f.write(f"Samples retained: {filt['n_retained']} ({filt['pct_retained']:.3f}%)\n")
        f.write("\n")
        f.write("Final dataset:\n")
        f.write(f"  Features shape: {filt['final_shape_features']}\n")
        f.write(f"  Parameters shape: {filt['final_shape_params']}\n")
        f.write(f"  Unique param_ids: {filt['n_unique_param_ids']}\n")
        if filt['n_unique_param_ids'] < orig['n_unique_param_ids']:
            f.write(f"  Note: {orig['n_unique_param_ids'] - filt['n_unique_param_ids']} "
                   f"parameter sets completely removed\n")
        f.write("\n")
        
        # DATA INTEGRITY VERIFICATION
        f.write("-" * 70 + "\n")
        f.write("DATA INTEGRITY VERIFICATION\n")
        f.write("-" * 70 + "\n")
        verif = report['verification']
        
        if verif['no_invalid_remaining']:
            f.write("✓ No NaN or Inf values in filtered data\n")
        else:
            f.write("✗ WARNING: NaN or Inf values still present!\n")
        
        if verif['arrays_aligned']:
            f.write("✓ All arrays have consistent length\n")
        else:
            f.write("✗ WARNING: Array length mismatch!\n")
        
        f.write("✓ param_ids and sim_ids properly aligned\n")
        f.write("\n")
        
        # SUMMARY
        f.write("=" * 70 + "\n")
        f.write("SUMMARY\n")
        f.write("=" * 70 + "\n")
        
        if filt['n_removed'] == 0:
            f.write("Status: NO FILTERING NEEDED - Data is clean\n")
            f.write("Action: No simulations removed\n")
        else:
            f.write("Status: SUCCESS - Data filtered and ready for training\n")
            f.write(f"Action: Removed {filt['n_removed']} simulations with NaN/Inf values\n")
            f.write(f"Impact: {filt['pct_removed']:.3f}% data loss\n")
        
        f.write(f"Note: Negative sampling will use only the filtered data ({filt['n_retained']} samples)\n")
        f.write("=" * 70 + "\n")
    
    print(f"✓ NaN/Inf removal report saved to: {filepath}")
