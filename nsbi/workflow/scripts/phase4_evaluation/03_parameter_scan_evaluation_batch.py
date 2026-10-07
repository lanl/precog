"""
Parameter Scan Evaluation for NRE Model

For each positive test case, scans across a uniform grid of (R0, recovery_time)
values and obtains model predictions. Outputs results in long format CSV for
visualization and analysis of the likelihood surface.

This script loads a pre-generated parameter grid to ensure all batches use
identical grid values.

Output CSV columns:
- param_id: Parameter set ID from test set
- replicate: Replicate number for this parameter set
- sim_id: Simulation identifier (param_id_replicate)
- R0: Candidate R0 value (grid point)
- recovery_time: Candidate recovery time value (grid point)
- predicted_probability: Model prediction [0-1] for this (statistics, parameters) pair
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import yaml
import sys
from pathlib import Path

# Add utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, 'workflow/scripts')
from utils.nre_model import NREClassifier
from utils.normalization import ParameterNormalizer

# Access snakemake variables
test_data_file = snakemake.input.test_data
model_path = snakemake.input.model
config_file = snakemake.input.config_file
sim_params_file = snakemake.input.sim_params
grid_file = snakemake.input.grid_file  # Pre-generated parameter grid
batch_output = snakemake.output[0]
sim_ids_batch = snakemake.params.sim_ids  # Batch-specific test case IDs
model_size = snakemake.params.model_size  # e.g., 'small', 'medium', 'large'
batch_id = snakemake.wildcards.batch_id
log_file = snakemake.log[0]
n_r0_points = snakemake.params.n_r0_points
n_rt_points = snakemake.params.n_rt_points

def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def load_model(checkpoint_path, config, device, feature_dim, model_size):
    """Load trained model from checkpoint.
    
    Args:
        checkpoint_path: Path to model checkpoint
        config: Configuration dict
        device: torch device
        feature_dim: Number of features in test data (determines input_dim with +2 for parameters)
        model_size: Size variant (e.g., 'small', 'medium', 'large')
    """
    use_masking = config['model'].get('use_masking', False)
    n_params = 2  # R0, recovery_time
    
    # Adjust input_dim based on masking configuration (must match training)
    if use_masking:
        # When masking enabled, masks are concatenated as additional features
        input_dim = feature_dim + feature_dim + n_params
    else:
        input_dim = feature_dim + n_params
    
    # Get architecture from model size configuration
    if model_size not in config['model_sizes']:
        raise ValueError(f"Unknown model size '{model_size}'. Available: {list(config['model_sizes'].keys())}")
    
    size_config = config['model_sizes'][model_size]
    hidden_dims = size_config['hidden_dims']
    dropout_rate = config['model']['dropout_rate']
    
    model = NREClassifier(
        input_dim=input_dim,
        hidden_dims=hidden_dims,
        dropout_rate=dropout_rate,
        use_batch_norm=config['model']['use_batch_norm'],
        use_masking=use_masking,
        n_features=feature_dim  # Required for masking to work correctly
    ).to(device)
    
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    return model, checkpoint['epoch'], checkpoint['loss']

def load_parameter_normalizer(test_data_file):
    """
    Load parameter normalization statistics from test data file.
    
    Args:
        test_data_file: Path to test data NPZ file with normalization stats
    
    Returns:
        ParameterNormalizer instance configured with training statistics
    """
    data = np.load(test_data_file, allow_pickle=True)
    
    # Check if normalization statistics are present
    if 'param_mean' not in data or 'param_std' not in data or 'param_is_log' not in data:
        raise ValueError(
            f"Test data file {test_data_file} is missing normalization statistics. "
            "Expected keys: param_mean, param_std, param_is_log"
        )
    
    param_mean = data['param_mean']
    param_std = data['param_std']
    param_is_log = data['param_is_log']
    
    # Create and configure normalizer
    normalizer = ParameterNormalizer(is_log_scale=param_is_log)
    normalizer.mean_ = param_mean
    normalizer.std_ = param_std
    normalizer.fitted = True
    
    return normalizer

def load_parameter_grid(grid_file, expected_n_r0, expected_n_rt):
    """
    Load pre-generated parameter grid from file and validate dimensions.
    
    Args:
        grid_file: Path to NPZ file containing grid
        expected_n_r0: Expected number of R0 points
        expected_n_rt: Expected number of recovery_time points
    
    Returns:
        r0_grid: (n_r0_points,) array
        rt_grid: (n_rt_points,) array
        metadata: Dict with grid parameters
    """
    grid_data = np.load(grid_file, allow_pickle=True)
    r0_grid = grid_data['r0_grid']
    rt_grid = grid_data['rt_grid']
    metadata = grid_data['metadata'].item()
    
    # Validate dimensions
    if len(r0_grid) != expected_n_r0:
        raise ValueError(
            f"Grid dimension mismatch: expected {expected_n_r0} R0 points, "
            f"but grid file has {len(r0_grid)}"
        )
    if len(rt_grid) != expected_n_rt:
        raise ValueError(
            f"Grid dimension mismatch: expected {expected_n_rt} RT points, "
            f"but grid file has {len(rt_grid)}"
        )
    
    return r0_grid, rt_grid, metadata

# Setup logging
log_dir = Path(log_file).parent
log_dir.mkdir(parents=True, exist_ok=True)

with open(log_file, 'w') as log:
    log.write("="*70 + "\n")
    log.write(f"Parameter Scan - Batch {batch_id}\n")
    log.write("="*70 + "\n")
    log.flush()
    
    config = load_config(config_file)
    log.write(f"Loaded NRE config from: {config_file}\n")
    log.flush()
    
    log.write(f"Loading prior config from: {sim_params_file}\n")
    prior_config = load_config(sim_params_file)
    prior_bounds = {
        'R0_min': prior_config['R0_min'],
        'R0_max': prior_config['R0_max'],
        'recovery_time_min': prior_config['recovery_time_min'],
        'recovery_time_max': prior_config['recovery_time_max']
    }
    log.write(f"Prior: R0[{prior_bounds['R0_min']}-{prior_bounds['R0_max']}], RT[{prior_bounds['recovery_time_min']}-{prior_bounds['recovery_time_max']}]\n")
    log.flush()
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    log.write(f"Using device: {device}\n")
    log.flush()
    
    log.write(f"Loading test data from: {test_data_file}\n")
    data = np.load(test_data_file, allow_pickle=True)
    X_test_all = data['features']
    sim_ids_all = data['sim_ids']
    feature_dim = X_test_all.shape[1]
    
    # Load masks if available
    use_masking = config['model'].get('use_masking', False)
    if 'feature_masks' in data and use_masking:
        masks_all = data['feature_masks']
        log.write(f"  Loaded masks: {masks_all.shape}\n")
        log.write(f"  Avg valid bins: {masks_all.sum(axis=1).mean():.1f}/{masks_all.shape[1]}\n")
    else:
        masks_all = None
        if use_masking:
            log.write("  WARNING: use_masking=True but no masks found in data file. Disabling masking.\n")
            use_masking = False
    
    log.write(f"  Total test samples: {len(X_test_all)}\n")
    log.write(f"  Feature dimension: {feature_dim}\n")
    log.write(f"  Masking enabled: {use_masking}\n")
    log.flush()
    
    # Load parameter normalization statistics
    log.write(f"\\nLoading parameter normalization statistics from test data...\\n")
    param_normalizer = load_parameter_normalizer(test_data_file)
    log.write(f"  Parameter normalization loaded successfully\\n")
    log.write(f"  param_mean: {param_normalizer.mean_}\\n")
    log.write(f"  param_std: {param_normalizer.std_}\\n")
    log.write(f"  param_is_log: {param_normalizer.is_log_scale}\\n")
    log.flush()
    
    # Validate normalization with a test case
    log.write(f"\\nValidating parameter normalization...\\n")
    test_params = np.array([[5.0, 7.0]])  # Mid-range values
    test_params_norm = param_normalizer.transform(test_params)
    log.write(f"  Test: [R0=5.0, RT=7.0] -> normalized: {test_params_norm[0]}\\n")
    test_params_back = param_normalizer.inverse_transform(test_params_norm)
    log.write(f"  Inverse test: {test_params_norm[0]} -> {test_params_back[0]}\\n")
    log.flush()

    
    log.write(f"Loading model from: {model_path}\n")
    model, trained_epochs, final_loss = load_model(model_path, config, device, feature_dim, model_size)
    log.write(f"  Trained for {trained_epochs} epochs, loss={final_loss:.6f}\n")
    
    # Log correct input dimension based on masking
    if use_masking:
        model_input_dim = feature_dim + feature_dim + 2  # features + masks + params
        log.write(f"  Model input dimension: {model_input_dim} (features={feature_dim} + masks={feature_dim} + params=2)\n")
    else:
        model_input_dim = feature_dim + 2
        log.write(f"  Model input dimension: {model_input_dim} (features={feature_dim} + params=2)\n")
    log.flush()
    
    log.write(f"Filtering to batch {batch_id} ({len(sim_ids_batch)} IDs)\n")
    batch_indices = [i for i, sid in enumerate(sim_ids_all) if sid in sim_ids_batch]
    X_test = X_test_all[batch_indices]
    sim_ids = [sim_ids_all[i] for i in batch_indices]
    
    # Filter masks if available
    if masks_all is not None:
        masks_test = masks_all[batch_indices]
    else:
        masks_test = None
    
    log.write(f"  Filtered to {len(sim_ids)} test cases\n")
    log.flush()
    
    log.write(f"Loading pre-generated grid from: {grid_file}\n")
    r0_grid, rt_grid, grid_metadata = load_parameter_grid(grid_file, n_r0_points, n_rt_points)
    log.write(f"  Grid dimensions validated: {n_r0_points} R0 × {n_rt_points} RT = {n_r0_points * n_rt_points} points per case\n")
    log.write(f"  R0 range: [{r0_grid.min():.4f}, {r0_grid.max():.4f}]\n")
    log.write(f"  RT range: [{rt_grid.min():.4f}, {rt_grid.max():.4f}]\n")
    log.write(f"  Grid metadata: R0 spacing={grid_metadata['r0_spacing']}, RT spacing={grid_metadata['rt_spacing']}\n")
    log.flush()
    
    log.write(f"Scanning {len(sim_ids)} test cases...\n")
    log.flush()
    
    first_sim_probs = []  # Track probabilities for first simulation for diagnostics

    results = []
    for i, sim_id in enumerate(sim_ids):
        if (i + 1) % 10 == 0 or i == 0 or i == len(sim_ids) - 1:
            log.write(f"  [{i+1}/{len(sim_ids)}] {sim_id}\n")
            log.flush()
        
        parts = sim_id.split('_')
        param_id = int(parts[1][1:])
        replicate = int(parts[2][1:])
        x_obs = X_test[i]
        
        # Get mask for this observation (if available)
        mask_obs = masks_test[i] if masks_test is not None else None
        
        for r0_val in r0_grid:
            for rt_val in rt_grid:
                # Create candidate parameter (raw/unnormalized)
                theta_candidate = np.array([r0_val, rt_val])
                
                # CRITICAL FIX: Normalize parameters before feeding to model
                # Model was trained on normalized parameters, so grid values must be normalized too
                theta_candidate_normalized = param_normalizer.transform(theta_candidate.reshape(1, -1))[0]
                
                # Build input vector based on masking configuration
                # MASKING FIX: Include masks in input when masking is enabled
                if mask_obs is not None and use_masking:
                    # Include masks as features when masking is enabled
                    # Input: [features, masks, parameters]
                    input_vec = np.concatenate([x_obs, mask_obs, theta_candidate_normalized])
                else:
                    # Standard input without masking
                    # Input: [features, parameters]
                    input_vec = np.concatenate([x_obs, theta_candidate_normalized])
                
                input_tensor = torch.FloatTensor(input_vec).unsqueeze(0).to(device)
                
                with torch.no_grad():
                    logit = model(input_tensor)  # No mask argument - already concatenated
                    prob = torch.sigmoid(logit).cpu().numpy()[0, 0]
                
                # Track first simulation probabilities for diagnostics
                if i == 0:
                    first_sim_probs.append(prob)
                
                results.append({
                    'param_id': param_id,
                    'replicate': replicate,
                    'sim_id': sim_id,
                    'R0': r0_val,
                    'recovery_time': rt_val,
                    'predicted_probability': prob
                })
        
        # Log diagnostics after first simulation
        if i == 0 and len(first_sim_probs) > 0:
            log.write(f"  First case diagnostics:\\n")
            log.write(f"    Prob range: [{min(first_sim_probs):.4e}, {max(first_sim_probs):.4e}]\\n")
            log.write(f"    Mean prob: {np.mean(first_sim_probs):.4f}\\n")
            log.flush()
    
    results_df = pd.DataFrame(results)
    log.write(f"Created DataFrame with {len(results_df)} predictions\n")
    log.flush()
    
    output_dir = Path(batch_output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(batch_output, index=False)
    log.write(f"Saved to: {batch_output}\n")
    
    log.write(f"Stats: Mean={results_df['predicted_probability'].mean():.4f}, ")
    log.write(f"Std={results_df['predicted_probability'].std():.4f}\n")
    log.write(f"Batch {batch_id} complete!\n")

print(f"Batch {batch_id}: {len(results_df)} predictions")
print("="*70)
print("Parameter Scan Evaluation")
print("="*70)

# Load config
config = load_config(config_file)
print(f"\nLoaded NRE config from: {config_file}")
