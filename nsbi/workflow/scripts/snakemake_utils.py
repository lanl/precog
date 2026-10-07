"""
Utility functions for Snakemake workflow.
Provides ID parsing, batch creation, and other helper functions.
"""

import re
from typing import List, Tuple, Optional


# ============================================================================
# Simulation ID Parsing
# ============================================================================

def get_param_id_from_sim_id(sim_id: str) -> Optional[int]:
    """Extract param_id from sim_id like 'train_p00042_r2' -> 42"""
    match = re.search(r'_p(\d+)_r\d+', sim_id)
    return int(match.group(1)) if match else None


def get_replicate_from_sim_id(sim_id: str) -> Optional[int]:
    """Extract replicate from sim_id like 'train_p00042_r2' -> 2"""
    match = re.search(r'_r(\d+)', sim_id)
    return int(match.group(1)) if match else None


def get_dataset_from_sim_id(sim_id: str) -> str:
    """Extract dataset from sim_id like 'train_p00042_r2' -> 'train'"""
    return sim_id.split('_')[0]


# ============================================================================
# Batch Creation
# ============================================================================

def get_batches(ids: List[str], batch_size: int) -> List[Tuple[int, List[str]]]:
    """
    Split IDs into batches.
    
    Args:
        ids: List of ID strings to batch
        batch_size: Number of IDs per batch
        
    Returns:
        List of (batch_id, batch_ids) tuples
    """
    batches = []
    for i in range(0, len(ids), batch_size):
        batch_id = i // batch_size
        batch_items = ids[i:i+batch_size]
        batches.append((batch_id, batch_items))
    return batches


def get_batch_id_for_item(item_id: str, all_batches: List[Tuple[int, List[str]]]) -> int:
    """
    Get the batch_id for a given item_id.
    
    Args:
        item_id: The ID to find
        all_batches: List of (batch_id, batch_items) tuples
        
    Returns:
        batch_id containing the item
        
    Raises:
        ValueError: If item not found in any batch
    """
    for batch_id, batch_items in all_batches:
        if item_id in batch_items:
            return batch_id
    raise ValueError(f"Item {item_id} not found in any batch")


# ============================================================================
# Baseline Batching

# ============================================================================
# Configuration Helpers
# ============================================================================

def calculate_n_train_batches(n_train_xmls: int, batch_size: int) -> int:
    """Calculate number of training batches needed"""
    return (n_train_xmls + batch_size - 1) // batch_size


def calculate_n_test_batches(n_test_xmls: int, batch_size: int) -> int:
    """Calculate number of test batches needed"""
    return (n_test_xmls + batch_size - 1) // batch_size
