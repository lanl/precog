"""
Negative sampling utilities for NRE training.

Provides different strategies for generating negative samples:
- Within-batch: Permute parameters within current batch (sbi approach)
- Cross-batch: Sample parameters from full dataset (more diversity)
- Prior random: Sample from prior distribution (easy negatives)

All strategies ensure no self-matches: negative pairs (x_i, theta_j) always
have param_id[i] != param_id[j], preventing replicates from matching with
their own parameter values.
"""

import numpy as np


def create_within_batch_negatives(theta_batch, param_ids_batch, 
                                   n_negatives_per_positive=1, seed=None):
    """
    Create negative samples by permuting parameters within the current batch only.
    
    This is the approach used by the `sbi` package (within-batch permutation):
    - Fast: O(B) per batch
    - Perfect marginal matching: negatives come from exact same distribution as batch
    - Simple: just reorder the batch parameters
    - Works well for batch_size >= 32
    
    Ensures no self-matches: For each sample i, the negative parameter comes
    from a different param_id to prevent replicates from being paired together.
    
    Args:
        theta_batch: (B, d) array of parameters for current batch
        param_ids_batch: (B,) array of parameter IDs for current batch
        n_negatives_per_positive: Number of negative samples per positive (default=1)
        seed: Random seed for reproducibility
        
    Returns:
        theta_neg: (B * n_negatives_per_positive, d) array of permuted parameters
        collision_info: Dict with statistics about collision fixing
        
    Raises:
        ValueError: If batch contains only one unique param_id (cannot create valid negatives)
        Warning: If batch_size < 32 (may have insufficient diversity)
    """
    if seed is not None:
        np.random.seed(seed)
    
    B = len(theta_batch)
    d = theta_batch.shape[1]
    
    # Check if we have enough unique param_ids
    unique_param_ids = np.unique(param_ids_batch)
    if len(unique_param_ids) == 1:
        raise ValueError(
            f"Cannot create valid negatives: batch contains only one unique param_id "
            f"({unique_param_ids[0]}). All samples would self-match."
        )
    
    # Storage for all negatives
    all_negatives = []
    total_initial_collisions = 0
    total_fixed_collisions = 0
    total_iterations = 0
    
    # Generate n_negatives_per_positive permutations
    for neg_idx in range(n_negatives_per_positive):
        # Create initial random permutation
        perm = np.random.permutation(B)
        
        # Check for self-matches (where param_id[i] == param_id[perm[i]])
        self_matches = (param_ids_batch == param_ids_batch[perm])
        n_initial_collisions = self_matches.sum()
        total_initial_collisions += n_initial_collisions
        
        # Fix self-matches by iterative swapping
        max_iterations = 100
        iteration = 0
        
        while self_matches.any() and iteration < max_iterations:
            match_idxs = np.where(self_matches)[0]
            non_match_idxs = np.where(~self_matches)[0]
            
            if len(non_match_idxs) == 0:
                # All positions are self-matches - should not happen if we have >1 unique param_id
                raise RuntimeError(
                    f"Failed to fix self-matches: all positions match. "
                    f"This should not happen. Batch param_ids: {np.unique(param_ids_batch)}"
                )
            
            # Swap each self-match with a random non-self-match
            for match_idx in match_idxs:
                swap_idx = np.random.choice(non_match_idxs)
                perm[match_idx], perm[swap_idx] = perm[swap_idx], perm[match_idx]
            
            # Re-compute self-matches after swapping
            self_matches = (param_ids_batch == param_ids_batch[perm])
            iteration += 1
        
        total_iterations += iteration
        
        if self_matches.any():
            raise RuntimeError(
                f"Failed to fix self-matches after {max_iterations} iterations. "
                f"Remaining collisions: {self_matches.sum()}"
            )
        
        total_fixed_collisions += n_initial_collisions
        
        # Extract permuted parameters
        theta_neg = theta_batch[perm]
        all_negatives.append(theta_neg)
    
    # Combine all negatives
    theta_neg_combined = np.concatenate(all_negatives, axis=0)
    
    collision_info = {
        'method': 'within_batch',
        'total_initial_collisions': total_initial_collisions,
        'total_fixed_collisions': total_fixed_collisions,
        'total_iterations': total_iterations,
        'n_negatives_per_positive': n_negatives_per_positive,
        'batch_size': B,
        'unique_param_ids': len(unique_param_ids)
    }
    
    return theta_neg_combined, collision_info


def create_permuted_negatives(theta, param_ids, n_negatives_per_positive=1, method='random', seed=None):
    """
    Create negative samples by permuting parameters while avoiding self-matches.
    
    This ensures that for each sample i with param_id[i], the negative sample(s)
    have different param_ids, preventing replicates from matching with their own
    parameter values.
    
    Args:
        theta: (N, d) array of parameters
        param_ids: (N,) array of parameter IDs for each sample
        n_negatives_per_positive: Number of negative samples per positive (default=1)
        method: 'random' for random permutation, 'circular' for deterministic shift
        seed: Random seed for reproducibility (only used if method='random')
        
    Returns:
        theta_neg: (N * n_negatives_per_positive, d) array of permuted parameters
        collision_info: Dict with statistics about collision fixing
    """
    if seed is not None and method == 'random':
        np.random.seed(seed)
    
    N = len(theta)
    d = theta.shape[1]
    
    # Storage for all negatives
    all_negatives = []
    total_initial_collisions = 0
    total_fixed_collisions = 0
    
    # Generate n_negatives_per_positive permutations
    for neg_idx in range(n_negatives_per_positive):
        # Create initial permutation
        if method == 'circular':
            # Deterministic circular shift
            shift = neg_idx + 1
            perm_indices = np.roll(np.arange(N), shift)
        elif method == 'random':
            # Random permutation
            perm_indices = np.random.permutation(N)
        else:
            raise ValueError(f"Unknown permutation method: {method}")
        
        # Get permuted param_ids
        param_ids_permuted = param_ids[perm_indices]
        
        # Find collisions (where param_id matches after permutation)
        collisions = (param_ids == param_ids_permuted)
        n_collisions = collisions.sum()
        total_initial_collisions += n_collisions
        
        if n_collisions > 0:
            # Fix collisions by swapping
            collision_indices = np.where(collisions)[0]
            
            for col_idx in collision_indices:
                # Find a non-collision index to swap with
                current_param = param_ids[col_idx]
                current_perm_param = param_ids[perm_indices[col_idx]]
                
                # Try to find a valid swap partner
                found_swap = False
                for candidate in np.random.permutation(N):
                    if candidate == col_idx:
                        continue
                    
                    candidate_param = param_ids[candidate]
                    candidate_perm_param = param_ids[perm_indices[candidate]]
                    
                    # Check if swap would fix collision without creating new ones
                    would_fix_col = (current_param != candidate_perm_param)
                    would_fix_candidate = (candidate_param != current_perm_param)
                    
                    if would_fix_col and would_fix_candidate:
                        # Perform swap
                        perm_indices[col_idx], perm_indices[candidate] = \
                            perm_indices[candidate], perm_indices[col_idx]
                        found_swap = True
                        total_fixed_collisions += 1
                        break
                
                if not found_swap:
                    # Fallback: use random different param_id
                    candidates = np.where(param_ids != current_param)[0]
                    if len(candidates) > 0:
                        swap_target = np.random.choice(candidates)
                        perm_indices[col_idx], perm_indices[swap_target] = \
                            perm_indices[swap_target], perm_indices[col_idx]
                        total_fixed_collisions += 1
        
        # Verify no self-matches remain
        param_ids_permuted_final = param_ids[perm_indices]
        remaining_collisions = (param_ids == param_ids_permuted_final).sum()
        
        if remaining_collisions > 0:
            print(f"  ⚠️  Warning: {remaining_collisions} self-matches remain after fixing")
        
        # Get permuted parameters
        theta_neg = theta[perm_indices]
        all_negatives.append(theta_neg)
    
    # Combine all negatives
    theta_neg_combined = np.concatenate(all_negatives, axis=0)
    
    collision_info = {
        'total_initial_collisions': total_initial_collisions,
        'total_fixed_collisions': total_fixed_collisions,
        'n_permutations': n_negatives_per_positive,
        'samples_per_permutation': N
    }
    
    return theta_neg_combined, collision_info


def create_cross_batch_negatives(theta_batch, param_ids_batch, 
                                  theta_full, param_ids_full,
                                  n_negatives_per_positive=1, seed=None):
    """
    Create negative samples by drawing from the full dataset (cross-batch sampling).
    
    For each sample in the batch, randomly selects a negative parameter from the
    entire dataset, ensuring different param_id to avoid self-matches.
    
    Advantages over within-batch:
    - Higher diversity: access to all N parameters in dataset
    - Better for small batches: not limited by batch size
    - Better for imbalanced data: can always find rare parameters
    
    Disadvantages:
    - Slower: O(B × N) operations per batch
    - Sampling with replacement: can duplicate negatives within batch
    - Requires passing full dataset to training loop
    
    Args:
        theta_batch: (B, d) parameters for current batch
        param_ids_batch: (B,) param_ids for current batch
        theta_full: (N, d) all parameters in dataset (to draw negatives from)
        param_ids_full: (N,) all param_ids in dataset
        n_negatives_per_positive: Number of negatives per positive
        seed: Random seed
        
    Returns:
        theta_neg: (B * n_negatives_per_positive, d) negative parameters
        collision_info: Dict with statistics
    """
    if seed is not None:
        np.random.seed(seed)
    
    B = len(theta_batch)
    all_negatives = []
    n_collisions_fixed = 0
    
    for neg_idx in range(n_negatives_per_positive):
        theta_neg_batch = np.zeros_like(theta_batch)
        
        for i, param_id in enumerate(param_ids_batch):
            # Find all samples in full dataset with different param_id
            valid_indices = np.where(param_ids_full != param_id)[0]
            
            if len(valid_indices) == 0:
                # Fallback: if somehow all samples have same param_id, use random
                valid_indices = np.arange(len(param_ids_full))
                n_collisions_fixed += 1
            
            # Randomly select one
            selected_idx = np.random.choice(valid_indices)
            theta_neg_batch[i] = theta_full[selected_idx]
        
        all_negatives.append(theta_neg_batch)
    
    theta_neg_combined = np.concatenate(all_negatives, axis=0)
    
    collision_info = {
        'method': 'cross_batch',
        'n_collisions_fixed': n_collisions_fixed,
        'batch_size': B,
        'n_negatives_per_positive': n_negatives_per_positive,
        'dataset_size': len(theta_full)
    }
    
    return theta_neg_combined, collision_info

