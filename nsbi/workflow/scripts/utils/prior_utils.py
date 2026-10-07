"""
Prior sampling utilities for Neural SBI.

Provides consistent parameter sampling from prior distributions,
matching the LHS sampling scheme used in simulation design.
"""

import numpy as np


def sample_from_prior(n_samples, prior_config):
    """
    Sample parameters from the prior distribution.
    
    Matches the LHS sampling scheme:
    - R0: Log10-uniform (geometric spacing)
    - recovery_time: Linear uniform
    - S: Linear uniform (if not fixed)
    
    Args:
        n_samples: Number of samples to generate
        prior_config: Dict with keys:
            - 'R0_min', 'R0_max': R0 bounds
            - 'recovery_time_min', 'recovery_time_max': Recovery time bounds
            - 'S_min', 'S_max': S bounds (if S_fixed is None)
            - 'S_fixed': Fixed S value or None
    
    Returns:
        If S_fixed: (n_samples, 2) array [R0, recovery_time]
        If S_fixed=None: (n_samples, 3) array [S, R0, recovery_time]
    """
    # Sample R0 on log10 scale (geometric spacing, matches LHS)
    log_R0_min = np.log10(prior_config['R0_min'])
    log_R0_max = np.log10(prior_config['R0_max'])
    u_R0 = np.random.uniform(0, 1, n_samples)
    R0_samples = 10 ** (log_R0_min + u_R0 * (log_R0_max - log_R0_min))
    
    # Sample recovery_time on linear scale (matches LHS)
    recovery_samples = np.random.uniform(
        prior_config['recovery_time_min'],
        prior_config['recovery_time_max'],
        n_samples
    )
    
    S_fixed = prior_config.get('S_fixed', None)
    
    if S_fixed is not None:
        # Return [R0, recovery_time] only
        return np.column_stack([R0_samples, recovery_samples])
    else:
        # Sample S on linear scale (matches LHS)
        S_samples = np.random.uniform(
            prior_config['S_min'],
            prior_config['S_max'],
            n_samples
        )
        # Return [S, R0, recovery_time]
        return np.column_stack([S_samples, R0_samples, recovery_samples])
