"""
Centralized seed management utilities for Neural SBI pipeline.

This module provides consistent random seed handling across all scripts
to ensure reproducibility and proper random number generation.

Author: Neural SBI Pipeline Refactoring
Date: 2026-08-19
Updated: 2026-09-09 - Added workflow seed system
"""

import numpy as np


def get_seed(config_seed):
    """
    Generate random seed if config_seed is None, otherwise return config_seed.
    
    This function handles null/None seeds consistently across the pipeline,
    generating a random seed when needed and logging the result.
    
    Args:
        config_seed: Seed value from configuration (int or None)
                     If None, generates a random seed
    
    Returns:
        int: The seed to use (either config_seed or newly generated)
    
    Example:
        >>> seed = get_seed(42)  # Returns 42
        >>> seed = get_seed(None)  # Generates and returns random seed
    """
    if config_seed is None:
        actual_seed = np.random.randint(0, 2**31 - 1)
        print(f"  No seed specified - generated random seed: {actual_seed}")
        return actual_seed
    print(f"  Using specified seed: {config_seed}")
    return config_seed


def derive_seed(master_seed, offset):
    """
    Derive a new seed from a master seed using an offset.
    
    This ensures different random streams for different purposes while
    maintaining reproducibility from a single master seed.
    
    Args:
        master_seed: The master/base seed
        offset: Integer offset to create derived seed
    
    Returns:
        int: Derived seed (master_seed + offset)
    
    Example:
        >>> training_seed = derive_seed(42, 0)      # 42
        >>> validation_seed = derive_seed(42, 1000)  # 1042
        >>> eval_seed = derive_seed(42, 2000)        # 2042
    """
    return master_seed + offset


def print_workflow_seed(seed, is_random=False):
    """
    Print the workflow seed prominently for reproducibility.
    
    Args:
        seed: The workflow seed being used
        is_random: Whether the seed was randomly generated (vs. from config)
    """
    print("\n" + "="*70)
    if is_random:
        print(f"🎲 WORKFLOW SEED: {seed} (randomly generated)")
        print("   Use this seed to reproduce these exact results:")
        print(f"   Set 'workflow_seed: {seed}' in config/nre_config.yaml")
    else:
        print(f"🎲 WORKFLOW SEED: {seed} (from config)")
        print("   Results are reproducible with this seed")
    print("="*70 + "\n")


def set_random_seeds(seed, verbose=True):
    """
    Set random seeds for all relevant libraries (numpy, torch, etc.).
    
    Args:
        seed: Random seed value
        verbose: Whether to print confirmation message
    
    Returns:
        seed: The seed that was set (for confirmation)
    """
    import numpy as np
    
    np.random.seed(seed)
    
    # Set torch seed if available
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass
    
    if verbose:
        print(f"  Random seeds set to: {seed}")
    
    return seed
