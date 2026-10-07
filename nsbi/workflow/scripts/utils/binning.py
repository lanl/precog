"""
Shared utility functions for binning time series data.

This module provides consistent binning methods used across the pipeline
to ensure Neural SBI and baseline methods process data identically.
"""

import numpy as np


def bin_time_series(times, values, bin_width, max_time=None):
    """
    Bin time series data into fixed-width bins using histogram method.
    
    This is the canonical binning function used throughout the pipeline.
    It matches the implementation in 04b_agg_single.py.
    
    Args:
        times: Array of event times
        values: Array of values (e.g., case counts) corresponding to each time
        bin_width: Width of each time bin
        max_time: Optional maximum time (if None, uses times.max())
    
    Returns:
        bin_times: Array of bin center times
        binned_values: Array of summed values in each bin
    
    Example:
        >>> times = np.array([0.5, 1.2, 1.8, 2.3, 3.1])
        >>> cases = np.array([1, 2, 1, 3, 2])
        >>> bin_times, binned = bin_time_series(times, cases, bin_width=1.0)
        >>> print(binned)  # [1, 3, 3, 2, 0]
    """
    if len(times) == 0:
        return np.array([]), np.array([])
    
    # Determine maximum time and number of bins
    mt = max_time if max_time is not None else times.max()
    n_bins = int(np.ceil(mt / bin_width)) + 1
    
    # Create bin edges
    bins = np.arange(0, (n_bins + 1) * bin_width, bin_width)
    
    # Bin the data using histogram with weights
    binned_values, _ = np.histogram(times, bins=bins, weights=values)
    
    # Create bin center times
    bin_times = np.arange(0, n_bins * bin_width, bin_width)
    
    return bin_times, binned_values[:n_bins]


def load_case_counts(filepath):
    """
    Load case counts from CSV file.
    
    Args:
        filepath: Path to case counts CSV file
        
    Returns:
        times: Array of event times
        cases: Array of new infection counts
    """
    import pandas as pd
    df = pd.read_csv(filepath)
    return df['t'].values, df['new_infections'].values
