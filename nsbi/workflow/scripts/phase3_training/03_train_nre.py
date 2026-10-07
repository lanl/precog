"""
Train Neural Ratio Estimation (NRE) model with early stopping.

Trains a binary classifier to distinguish true parameter-statistic pairs
from randomly shuffled pairs. Uses early stopping based on validation loss
convergence to automatically determine optimal number of epochs.
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import yaml
import sys

# Add utils to path and import model
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, 'workflow/scripts')
from utils.prior_utils import sample_from_prior
from utils.nre_model import NREClassifier

# Access snakemake variables
train_data_file = snakemake.input.train_data
config_file = snakemake.input.config_file
model_final_path = snakemake.output.model_final
model_best_path = snakemake.output.model_best
training_log_path = snakemake.output.training_log
training_seed = snakemake.params.training_seed
validation_seed = snakemake.params.validation_seed
model_size = snakemake.params.model_size  # e.g., 'small', 'medium', 'large'

class EarlyStopping:
    """Early stopping based on validation loss convergence."""
    
    def __init__(self, patience=20, min_delta=0.0001, mode='min'):
        """
        Args:
            patience: Number of epochs to wait for improvement
            min_delta: Minimum change to qualify as improvement
            mode: 'min' for loss (lower is better)
        """
        self.patience = patience
        self.min_delta = min_delta
        self.mode = mode
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.best_epoch = 0
        
    def __call__(self, current_score, epoch):
        """Check if training should stop."""
        if self.best_score is None:
            self.best_score = current_score
            self.best_epoch = epoch
            return False
        
        improved = current_score < (self.best_score - self.min_delta)
        
        if improved:
            self.best_score = current_score
            self.best_epoch = epoch
            self.counter = 0
            return False
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
                return True
            return False

def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def set_seed(seed):
    """Set random seeds for reproducibility."""
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

def create_batches(X, theta, masks, batch_size):
    """Create mini-batches from data with optional masks.
    
    Returns:
        Generator yielding (X_batch, theta_batch, masks_batch, batch_indices)
        where batch_indices are the actual indices into the original arrays.
    """
    n_samples = len(X)
    indices = np.random.permutation(n_samples)
    
    for start_idx in range(0, n_samples, batch_size):
        end_idx = min(start_idx + batch_size, n_samples)
        batch_indices = indices[start_idx:end_idx]
        if masks is not None:
            yield X[batch_indices], theta[batch_indices], masks[batch_indices], batch_indices
        else:
            yield X[batch_indices], theta[batch_indices], None, batch_indices

def train_epoch(model, optimizer, X, theta, param_ids, masks, prior_config, 
                negative_config, batch_size, device, gradient_clip_norm=None):
    """Train for one epoch with configurable negative sampling.
    
    Supports three negative sampling strategies:
    1. 'within_batch': Permute parameters within current batch (sbi approach)
       - Fast, simple, perfect marginal matching
       - Recommended for batch_size >= 32
    2. 'cross_batch': Sample parameters from full dataset
       - More diversity, better for small batches or imbalanced data
       - Slower but more flexible
    3. 'prior_random': Sample from prior distribution
       - Creates easy negatives (may be out-of-distribution)
       - Use with higher negative_ratio (5-10)
    
    All strategies prevent self-matches: negative pairs (x_i, theta_j) always
    have param_id[i] != param_id[j], ensuring different parameter values even
    with multiple replicates.
    
    Args:
        model: NRE model
        optimizer: Optimizer
        X: Features (N, d_features)
        theta: Parameters (N, d_params)
        param_ids: Parameter IDs (N,) - for permuted negatives
        masks: Masks (N, d_features) or None
        prior_config: Prior bounds for random sampling
        negative_config: Dict with 'strategy', 'method', 'ratio'
        batch_size: Batch size
        device: torch device
        gradient_clip_norm: If not None, clip gradients to this max norm for stability
    
    Returns:
        Average training loss
    """
    from utils.negative_sampling import (
        create_within_batch_negatives,
        create_cross_batch_negatives,
    )
    
    model.train()
    total_loss = 0.0
    n_batches = 0
    
    criterion = nn.BCEWithLogitsLoss()
    
    # Get negative sampling configuration
    negative_strategy = negative_config.get('strategy', 'within_batch')
    negative_ratio = int(negative_config.get('ratio', 1))
    
    for batch_data in create_batches(X, theta, masks, batch_size):
        X_batch, theta_batch, masks_batch, batch_indices = batch_data
        batch_size_actual = len(X_batch)
        
        # Get param_ids for this batch using actual batch indices
        param_ids_batch = param_ids[batch_indices]
        
        # Create negative samples based on strategy
        if negative_strategy == 'within_batch':
            # NEW: sbi-style within-batch permutation (fast, simple)
            # Permutes parameters within the current batch only
            theta_neg, _ = create_within_batch_negatives(
                theta_batch,
                param_ids_batch,
                n_negatives_per_positive=negative_ratio,
                seed=None  # Different each batch for variety
            )
            
        elif negative_strategy == 'cross_batch':
            # Sample from full dataset (more diversity, slower)
            theta_neg, _ = create_cross_batch_negatives(
                theta_batch,
                param_ids_batch,
                theta,  # Full dataset
                param_ids,  # Full param_ids
                n_negatives_per_positive=negative_ratio,
                seed=None  # Different each batch for variety
            )
            
        elif negative_strategy == 'prior_random':
            # Sample from prior distribution (easy negatives)
            n_negatives_total = batch_size_actual * negative_ratio
            theta_neg = sample_from_prior(n_negatives_total, prior_config)
            
        else:
            raise ValueError(
                f"Unknown negative_strategy: '{negative_strategy}'. "
                f"Valid options: 'within_batch', 'cross_batch', 'prior_random'"
            )
        
        # Positive samples: true pairs (label=1)
        X_pos = torch.FloatTensor(X_batch).to(device)
        theta_pos = torch.FloatTensor(theta_batch).to(device)
        
        # Include masks in input if available
        if masks_batch is not None:
            masks_pos = torch.FloatTensor(masks_batch).to(device)
            input_pos = torch.cat([X_pos, masks_pos, theta_pos], dim=1)
        else:
            input_pos = torch.cat([X_pos, theta_pos], dim=1)
        
        labels_pos = torch.ones(batch_size_actual, 1).to(device)
        
        # Negative samples: features with shuffled/sampled parameters (label=0)
        # Repeat features to match number of negatives
        X_neg = np.tile(X_batch, (negative_ratio, 1))
        X_neg_t = torch.FloatTensor(X_neg).to(device)
        theta_neg_t = torch.FloatTensor(theta_neg).to(device)
        
        # Include masks in input if available (repeat masks for negatives)
        if masks_batch is not None:
            masks_neg = np.tile(masks_batch, (negative_ratio, 1))
            masks_neg_t = torch.FloatTensor(masks_neg).to(device)
            input_neg = torch.cat([X_neg_t, masks_neg_t, theta_neg_t], dim=1)
        else:
            input_neg = torch.cat([X_neg_t, theta_neg_t], dim=1)
        
        labels_neg = torch.zeros(len(theta_neg), 1).to(device)
        
        # Combine positive and negative samples
        inputs = torch.cat([input_pos, input_neg], dim=0)
        labels = torch.cat([labels_pos, labels_neg], dim=0)
        
        # Shuffle combined batch to prevent positional bias
        perm = torch.randperm(len(inputs))
        inputs = inputs[perm]
        labels = labels[perm]
        
        # Forward pass
        optimizer.zero_grad()
        outputs = model(inputs)  # No mask argument - already concatenated in inputs
        loss = criterion(outputs, labels)
        
        # Backward pass
        loss.backward()
        
        # Gradient clipping for stability
        if gradient_clip_norm is not None:
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=gradient_clip_norm)
        
        optimizer.step()
        
        total_loss += loss.item()
        n_batches += 1
    
    return total_loss / n_batches


def validate_epoch(model, X, theta, theta_neg_fixed, masks, batch_size, device, negative_ratio=1):
    """Validate with pre-generated fixed negative samples.
    
    Args:
        model: The NRE model
        X: Validation features
        theta: Validation parameters (positive samples)
        theta_neg_fixed: Pre-generated negative parameters (length = len(theta) * negative_ratio)
        masks: Optional masks for features
        batch_size: Batch size for validation
        device: torch device
        negative_ratio: Number of negatives per positive (must match theta_neg_fixed size)
    
    Returns:
        Average validation loss across all batches
    """
    model.eval()
    total_loss = 0.0
    n_batches = 0
    
    criterion = nn.BCEWithLogitsLoss()
    
    with torch.no_grad():
        # Use sequential ordering (not random) for consistent batches
        n_samples = len(X)
        indices = np.arange(n_samples)  # Sequential: [0, 1, 2, ...]
        
        for start_idx in range(0, n_samples, batch_size):
            end_idx = min(start_idx + batch_size, n_samples)
            batch_indices = indices[start_idx:end_idx]
            batch_size_actual = len(batch_indices)
            
            X_batch = X[batch_indices]
            theta_batch = theta[batch_indices]
            
            # Get corresponding negatives (accounting for negative_ratio)
            neg_start_idx = start_idx * negative_ratio
            neg_end_idx = end_idx * negative_ratio
            theta_neg_batch = theta_neg_fixed[neg_start_idx:neg_end_idx]
            
            masks_batch = masks[batch_indices] if masks is not None else None
            
            # Positive samples
            X_pos = torch.FloatTensor(X_batch).to(device)
            theta_pos = torch.FloatTensor(theta_batch).to(device)
            
            # Include masks in input if available
            if masks_batch is not None:
                masks_pos = torch.FloatTensor(masks_batch).to(device)
                input_pos = torch.cat([X_pos, masks_pos, theta_pos], dim=1)
            else:
                input_pos = torch.cat([X_pos, theta_pos], dim=1)
            
            labels_pos = torch.ones(batch_size_actual, 1).to(device)
            
            # Negative samples (features repeated for each negative)
            X_neg = np.tile(X_batch, (negative_ratio, 1))
            X_neg_t = torch.FloatTensor(X_neg).to(device)
            theta_neg_t = torch.FloatTensor(theta_neg_batch).to(device)
            
            # Include masks in input if available (repeat masks for negatives)
            if masks_batch is not None:
                masks_neg = np.tile(masks_batch, (negative_ratio, 1))
                masks_neg_t = torch.FloatTensor(masks_neg).to(device)
                input_neg = torch.cat([X_neg_t, masks_neg_t, theta_neg_t], dim=1)
            else:
                input_neg = torch.cat([X_neg_t, theta_neg_t], dim=1)
            
            labels_neg = torch.zeros(len(theta_neg_batch), 1).to(device)
            
            inputs = torch.cat([input_pos, input_neg], dim=0)
            labels = torch.cat([labels_pos, labels_neg], dim=0)
            
            # Shuffle within batch to prevent positional bias
            perm = torch.randperm(len(inputs))
            inputs = inputs[perm]
            labels = labels[perm]
            
            outputs = model(inputs)  # No mask argument - already concatenated in inputs
            loss = criterion(outputs, labels)
            
            total_loss += loss.item()
            n_batches += 1
    
    return total_loss / n_batches

def save_checkpoint(model, optimizer, epoch, loss, path):
    """Save model checkpoint."""
    torch.save({
        'epoch': epoch,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'loss': loss,
    }, path)

# Main training
print("="*70)
print("Training Neural Ratio Estimation Model")
print("="*70)

set_seed(training_seed)
print(f"\nTraining seed: {training_seed}")
print(f"Validation seed: {validation_seed}")
print(f"Model size: {model_size}")

config = load_config(config_file)
print(f"Loaded config from: {config_file}")

# Load prior configuration for negative sampling
print(f"\nLoading prior config from: config/simulation_params.yaml")
prior_config = load_config("config/simulation_params.yaml")
prior_bounds = {
    'R0_min': prior_config['R0_min'],
    'R0_max': prior_config['R0_max'],
    'recovery_time_min': prior_config['recovery_time_min'],
    'recovery_time_max': prior_config['recovery_time_max'],
    'S_fixed': prior_config.get('S_fixed', None)
}
if prior_bounds['S_fixed'] is None:
    prior_bounds['S_min'] = prior_config['S_min']
    prior_bounds['S_max'] = prior_config['S_max']
    print(f"  S prior: Uniform({prior_bounds['S_min']}, {prior_bounds['S_max']})")
else:
    print(f"  S fixed at: {prior_bounds['S_fixed']}")
print(f"  R0 prior: Log-Uniform({prior_bounds['R0_min']}, {prior_bounds['R0_max']})")
print(f"  recovery_time prior: Uniform({prior_bounds['recovery_time_min']}, {prior_bounds['recovery_time_max']})")

print(f"\nLoading training data from: {train_data_file}")
data = np.load(train_data_file)
X = data['features']
theta = data['parameters']

# Load param_ids (required for proper train/val splitting)
if 'param_ids' not in data:
    raise ValueError("ERROR: param_ids not found in training data file! This is required for proper train/val splitting.")
param_ids = data['param_ids']

# Load masks if available
use_masking = config['model'].get('use_masking', True)  # Default to True
if 'feature_masks' in data and use_masking:
    # Standard feature masks (for most models)
    masks = data['feature_masks']
    print(f"  Loaded feature_masks: {masks.shape}")
    print(f"  Avg valid features: {masks.sum(axis=1).mean():.1f}/{masks.shape[1]}")
elif 'masks' in data and use_masking:
    # Legacy bin masks (for case_counts model)
    masks = data['masks']
    print(f"  Loaded bin masks: {masks.shape}")
    print(f"  Avg valid bins: {masks.sum(axis=1).mean():.1f}/{masks.shape[1]}")
else:
    masks = None
    if use_masking:
        print("  WARNING: use_masking=True but no masks found in data file. Disabling masking.")
        use_masking = False

print(f"  Total samples: {len(X)}, Feature dim: {X.shape[1]}, Param dim: {theta.shape[1]}")
print(f"  Unique param_ids: {len(np.unique(param_ids))}")
print(f"  Masking enabled: {use_masking}")

# Split train/validation at PARAMETER level (not simulation level)
# This prevents data leakage where replicates of the same parameter are split across train/val
print("\n" + "="*70)
print("Splitting Train/Validation at Parameter Level")
print("="*70)

val_split = config['training']['validation_split']

# Get unique parameter IDs
unique_param_ids = np.unique(param_ids)
n_unique_params = len(unique_param_ids)
print(f"\nTotal unique parameters: {n_unique_params}")

# Split at parameter level (not simulation level)
n_val_params = int(n_unique_params * val_split)
n_train_params = n_unique_params - n_val_params

# Shuffle parameter IDs (not simulation IDs)
shuffled_param_ids = np.random.permutation(unique_param_ids)
val_param_ids_set = set(shuffled_param_ids[:n_val_params])
train_param_ids_set = set(shuffled_param_ids[n_val_params:])

print(f"Train parameters: {n_train_params}")
print(f"Validation parameters: {n_val_params}")

# Create train/val indices based on param_ids
# This ensures all replicates of a parameter stay together
train_indices = np.array([i for i, pid in enumerate(param_ids) if pid in train_param_ids_set])
val_indices = np.array([i for i, pid in enumerate(param_ids) if pid in val_param_ids_set])

print(f"\nTrain simulations: {len(train_indices)} (from {n_train_params} unique parameters)")
print(f"Validation simulations: {len(val_indices)} (from {n_val_params} unique parameters)")

# Verify no parameter overlap between train and val
train_params_in_split = set(param_ids[train_indices])
val_params_in_split = set(param_ids[val_indices])
overlap = train_params_in_split.intersection(val_params_in_split)
if len(overlap) > 0:
    raise ValueError(f"ERROR: Found {len(overlap)} parameters in both train and val sets! This should never happen.")
print(f"✓ Verified: No parameter overlap between train and validation")

# Sanity checks
assert len(train_indices) + len(val_indices) == len(X), "Train + val indices should equal total samples"
assert len(set(train_indices).intersection(set(val_indices))) == 0, "Train and val indices should not overlap"

# Check average replicates per parameter
avg_reps_train = len(train_indices) / len(train_params_in_split) if len(train_params_in_split) > 0 else 0
avg_reps_val = len(val_indices) / len(val_params_in_split) if len(val_params_in_split) > 0 else 0
print(f"Average replicates per parameter - Train: {avg_reps_train:.1f}, Val: {avg_reps_val:.1f}")

# Extract features, parameters, and param_ids for train/val
X_train, theta_train = X[train_indices], theta[train_indices]
X_val, theta_val = X[val_indices], theta[val_indices]
param_ids_train = param_ids[train_indices]
param_ids_val = param_ids[val_indices]

# Split masks if available
if masks is not None:
    masks_train = masks[train_indices]
    masks_val = masks[val_indices]
else:
    masks_train = None
    masks_val = None

print(f"\n{'='*70}")

# Get negative sampling configuration
negative_config = config.get('negative_sampling', {})
negative_strategy = negative_config.get('strategy', 'within_batch')
negative_ratio = int(config['training'].get('negative_ratio', 1))

print(f"\nNegative Sampling Configuration:")
print(f"  Strategy: {negative_strategy}")
print(f"  Negative ratio: {negative_ratio}:1 (negatives per positive)")

# Pre-generate fixed validation negatives for consistent evaluation
print(f"\nGenerating fixed validation negatives (seed={validation_seed})...")
np.random.seed(validation_seed)
torch.manual_seed(validation_seed)

if negative_strategy in ['within_batch', 'cross_batch']:
    # For both within_batch and cross_batch strategies, use cross_batch for validation
    # (validation uses full dataset, not batches)
    # Import here to avoid circular import issues
    from utils.negative_sampling import create_cross_batch_negatives
    
    theta_val_neg, val_collision_info = create_cross_batch_negatives(
        theta_val,
        param_ids_val,
        theta_val,  # Use validation set as full dataset
        param_ids_val,
        n_negatives_per_positive=negative_ratio,
        seed=validation_seed
    )
    print(f"  Generated {len(theta_val_neg)} permuted negative samples")
    print(f"  Fixed {val_collision_info['n_collisions_fixed']} self-match collisions")
elif negative_strategy == 'prior_random':
    n_val_negatives = len(X_val) * negative_ratio
    theta_val_neg = sample_from_prior(n_val_negatives, prior_bounds)
    print(f"  Generated {len(theta_val_neg)} random negative samples from prior")
else:
    raise ValueError(
        f"Unknown negative_strategy: '{negative_strategy}'. "
        f"Valid options: 'within_batch', 'cross_batch', 'prior_random'"
    )

print(f"  → Validation will use the same negatives every epoch for low-variance metrics")



device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# Initialize model with feature-specific architecture
n_features = X.shape[1]
n_params = theta.shape[1]

# Adjust input_dim based on masking configuration
if use_masking:
    # When masking enabled, masks are concatenated as additional features
    # Input becomes: [features, mask, parameters]
    input_dim = n_features + n_features + n_params  # features + masks + params
    print(f"\n⚠️  Masking enabled: masks will be concatenated as {n_features} additional features")
    print(f"   Input dim: {n_features} (features) + {n_features} (masks) + {n_params} (params) = {input_dim}")
else:
    input_dim = n_features + n_params
    print(f"\nInput dim: {n_features} (features) + {n_params} (params) = {input_dim}")

# Get architecture from model size configuration
if model_size is None:
    raise ValueError("model_size parameter is required (e.g., 'small', 'medium', 'large')")

if model_size not in config['model_sizes']:
    raise ValueError(f"Unknown model size '{model_size}'. Available: {list(config['model_sizes'].keys())}")

size_config = config['model_sizes'][model_size]
hidden_dims = size_config['hidden_dims']
dropout_rate = config['model']['dropout_rate']

print(f"\n→ Model size: {model_size}")
print(f"  Description: {size_config.get('description', 'N/A')}")
print(f"  Hidden dims: {hidden_dims}")
print(f"  Dropout: {dropout_rate}")

model = NREClassifier(
    input_dim=input_dim,
    hidden_dims=hidden_dims,
    dropout_rate=dropout_rate,
    use_batch_norm=config['model']['use_batch_norm'],
    use_masking=use_masking,
    n_features=n_features  # Required for masking to work correctly
).to(device)
print(f"\nModel: Input={input_dim}, Hidden={hidden_dims}, Dropout={dropout_rate}, Masking={use_masking}")

optimizer = optim.Adam(
    model.parameters(),
    lr=config['training']['learning_rate'],
    weight_decay=config['training']['weight_decay']
)

# Learning rate scheduler (optional)
scheduler = None
if config['training'].get('lr_scheduler', {}).get('enabled', False):
    scheduler_config = config['training']['lr_scheduler']
    if scheduler_config['type'] == 'ReduceLROnPlateau':
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            mode='min',
            factor=scheduler_config.get('factor', 0.5),
            patience=scheduler_config.get('patience', 5),
            min_lr=scheduler_config.get('min_lr', 1e-6),
        )
        print(f"Learning rate scheduler: ReduceLROnPlateau (factor={scheduler_config.get('factor', 0.5)}, patience={scheduler_config.get('patience', 5)})")

early_stopping = EarlyStopping(
    patience=config['training']['early_stopping']['patience'],
    min_delta=config['training']['early_stopping']['min_delta'],
    mode='min'
)
print(f"Early stopping: patience={early_stopping.patience}, min_delta={early_stopping.min_delta}")

# Gradient clipping config
gradient_clip_norm = None
if config['training'].get('gradient_clipping', {}).get('enabled', False):
    gradient_clip_norm = config['training']['gradient_clipping'].get('max_norm', 1.0)
    print(f"Gradient clipping: max_norm={gradient_clip_norm}")

# Training loop
best_val_loss = float('inf')
training_history = []
batch_size = config['training']['batch_size']
max_epochs = config['training']['max_epochs']
log_interval = config['training'].get('log_interval', 10)

# Warn once if configured batch_size is small for within_batch strategy
if negative_strategy == 'within_batch' and batch_size < 32:
    print(f"\n⚠️  Warning: within_batch strategy with batch_size={batch_size} (< 32).")
    print(f"   Consider increasing batch_size to ≥32 or using 'cross_batch' strategy for more diversity.")

print(f"\n{'='*70}")
print(f"Starting Training (max {max_epochs} epochs)")
print(f"{'='*70}\n")

for epoch in range(max_epochs):
    # Train with configurable negative sampling
    neg_config_dict = {
        'strategy': negative_strategy,
        'ratio': negative_ratio
    }
    train_loss = train_epoch(
        model, optimizer, X_train, theta_train, param_ids_train, 
        masks_train, prior_bounds, neg_config_dict, batch_size, device,
        gradient_clip_norm=gradient_clip_norm
    )
    
    # Validate with pre-generated fixed negatives
    val_loss = validate_epoch(
        model, X_val, theta_val, theta_val_neg, masks_val, batch_size, device,
        negative_ratio=negative_ratio
    )
    
    # Update learning rate scheduler
    if scheduler is not None:
        scheduler.step(val_loss)
    
    # Log
    if epoch % log_interval == 0 or epoch < 10:
        current_lr = optimizer.param_groups[0]['lr']
        print(f"Epoch {epoch:4d}: train_loss={train_loss:.6f}, val_loss={val_loss:.6f}, lr={current_lr:.6f}")
    
    # Save history
    training_history.append({
        'epoch': epoch,
        'train_loss': train_loss,
        'val_loss': val_loss
    })
    
    # Save best model
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        save_checkpoint(model, optimizer, epoch, val_loss, model_best_path)
        if epoch % log_interval == 0 or epoch < 10:
            print(f"  → New best model (val_loss={val_loss:.6f})")
    
    # Check early stopping
    if early_stopping(val_loss, epoch):
        print(f"\n{'='*70}")
        print(f"Early stopping at epoch {epoch}")
        print(f"Best val_loss: {early_stopping.best_score:.6f} at epoch {early_stopping.best_epoch}")
        print(f"No improvement for {early_stopping.patience} epochs")
        print(f"{'='*70}\n")
        break

# Save final model
save_checkpoint(model, optimizer, epoch, val_loss, model_final_path)

# Save training history
history_df = pd.DataFrame(training_history)
history_df.to_csv(training_log_path, index=False)

print(f"\n{'='*70}")
print(f"Training Complete")
print(f"{'='*70}")
print(f"Total epochs: {epoch + 1}")
print(f"Best val_loss: {best_val_loss:.6f}")
print(f"Final val_loss: {val_loss:.6f}")
print(f"Best model: {model_best_path}")
print(f"Final model: {model_final_path}")
print(f"Training log: {training_log_path}")



