"""
Probabilistic Scoring Rules for Posterior Evaluation

Implements proper scoring rules: CRPS, Energy Score, Interval Score, Log Score
These metrics complement coverage-based evaluation.

References:
- Gneiting & Raftery (2007). JASA.
- Gneiting & Ranjan (2011). JBES.
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

try:
    import properscoring as ps
    HAS_PROPERSCORING = True
except ImportError:
    HAS_PROPERSCORING = False
    import warnings
    warnings.warn("properscoring not available. Falling back to custom implementations.")


def sample_from_grid_posterior(grid_df, param_cols, n_samples=1000, random_state=None):
    """
    Draw samples from a discrete posterior distribution defined on a grid.
    
    Args:
        grid_df: DataFrame with parameter columns and 'posterior' column
        param_cols: List of parameter column names
        n_samples: Number of samples to draw (default: 1000)
        random_state: Random seed for reproducibility
        
    Returns:
        samples: Array of shape (n_samples, n_params)
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    grid_points = grid_df[param_cols].values
    posterior_probs = grid_df['posterior'].values
    posterior_probs = posterior_probs / posterior_probs.sum()
    
    sampled_indices = np.random.choice(
        len(grid_points), 
        size=n_samples, 
        replace=True, 
        p=posterior_probs
    )
    
    return grid_points[sampled_indices]


def compute_crps_from_samples(samples, true_value):
    """Compute CRPS from samples. Lower is better."""
    if HAS_PROPERSCORING:
        return ps.crps_ensemble(true_value, samples)
    else:
        samples = np.asarray(samples).flatten()
        term1 = np.mean(np.abs(samples - true_value))
        pairwise_diffs = np.abs(samples[:, None] - samples[None, :])
        term2 = 0.5 * np.mean(pairwise_diffs)
        return term1 - term2


def compute_energy_score(samples, true_params):
    """Compute Energy Score (multivariate CRPS). Lower is better."""
    # Properscoring doesn't include energy_score, so use custom implementation
    samples = np.asarray(samples)
    true_params = np.asarray(true_params).flatten()
    if samples.ndim == 1:
        samples = samples.reshape(-1, 1)
    
    distances_to_true = np.linalg.norm(samples - true_params, axis=1)
    term1 = np.mean(distances_to_true)
    
    pairwise_distances = cdist(samples, samples, metric='euclidean')
    triu_indices = np.triu_indices(len(samples), k=1)
    term2 = 0.5 * np.mean(pairwise_distances[triu_indices]) * 2
    
    return term1 - term2


def compute_interval_score(lower, upper, true_value, alpha=0.05):
    """Compute Interval Score for credible intervals. Lower is better."""
    width = upper - lower
    penalty_lower = (2.0 / alpha) * (lower - true_value) * (true_value < lower)
    penalty_upper = (2.0 / alpha) * (true_value - upper) * (true_value > upper)
    return width + penalty_lower + penalty_upper


def compute_log_score_nearest_neighbor(grid_df, param_cols, true_params):
    """Compute log posterior density at true parameters. Higher is better."""
    grid_points = grid_df[param_cols].values
    posterior_probs = grid_df['posterior'].values
    true_params = np.asarray(true_params).flatten()
    
    distances = np.linalg.norm(grid_points - true_params, axis=1)
    nearest_idx = np.argmin(distances)
    posterior_at_nearest = posterior_probs[nearest_idx]
    
    eps = 1e-300
    return np.log(posterior_at_nearest + eps)


def compute_all_scores(
    grid_df, 
    true_R0, 
    true_rt,
    R0_lower, 
    R0_upper, 
    rt_lower, 
    rt_upper,
    alpha=0.05,
    n_samples=1000,
    random_state=None
):
    """
    Compute all probabilistic scores for a single posterior.
    
    Returns dict with: crps_R0, crps_recovery_time, energy_score,
                      interval_score_R0, interval_score_recovery_time, log_score
    """
    param_cols = ['R0', 'recovery_time']
    true_params = np.array([true_R0, true_rt])
    
    # Draw samples from grid posterior
    samples = sample_from_grid_posterior(grid_df, param_cols, n_samples, random_state)
    
    # Compute scores
    crps_R0 = compute_crps_from_samples(samples[:, 0], true_R0)
    crps_rt = compute_crps_from_samples(samples[:, 1], true_rt)
    energy_score = compute_energy_score(samples, true_params)
    interval_score_R0 = compute_interval_score(R0_lower, R0_upper, true_R0, alpha)
    interval_score_rt = compute_interval_score(rt_lower, rt_upper, true_rt, alpha)
    log_score = compute_log_score_nearest_neighbor(grid_df, param_cols, true_params)
    
    return {
        'crps_R0': float(crps_R0),
        'crps_recovery_time': float(crps_rt),
        'energy_score': float(energy_score),
        'interval_score_R0': float(interval_score_R0),
        'interval_score_recovery_time': float(interval_score_rt),
        'log_score': float(log_score),
        'n_posterior_samples': n_samples
    }


def aggregate_scores(scores_df):
    """Compute summary statistics for probabilistic scores."""
    score_cols = [
        'crps_R0', 'crps_recovery_time', 'energy_score',
        'interval_score_R0', 'interval_score_recovery_time', 'log_score'
    ]
    
    summary = {}
    for col in score_cols:
        if col in scores_df.columns:
            summary[f'{col}_mean'] = float(scores_df[col].mean())
            summary[f'{col}_median'] = float(scores_df[col].median())
            summary[f'{col}_std'] = float(scores_df[col].std())
            summary[f'{col}_min'] = float(scores_df[col].min())
            summary[f'{col}_max'] = float(scores_df[col].max())
    
    return summary

