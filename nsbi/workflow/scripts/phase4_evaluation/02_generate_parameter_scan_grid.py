"""
Generate Parameter Scan Grid

Creates a uniform grid of parameter values for parameter scan evaluation.
This grid is generated once and shared across all batches to ensure consistency.

The grid is saved with metadata for reproducibility and validation.

Output file (NPZ format) contains:
- r0_grid: Array of R0 values (log-uniform spacing)
- rt_grid: Array of recovery_time values (linear spacing)
- metadata: Dict with grid parameters and bounds
"""

import os
import numpy as np
import yaml
from pathlib import Path

# Access snakemake variables
config_file = snakemake.input.config_file
sim_params_file = snakemake.input.sim_params
output_file = snakemake.output[0]
log_file = snakemake.log[0]
n_r0_points = snakemake.params.n_r0_points
n_rt_points = snakemake.params.n_rt_points

def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def create_parameter_grid(prior_bounds, n_r0_points, n_rt_points):
    """
    Create uniform grid of parameter values.
    
    R0: Log10-uniform spacing (geometric)
    recovery_time: Linear uniform spacing
    
    Returns:
        r0_grid: (n_r0_points,) array
        rt_grid: (n_rt_points,) array
    """
    # R0 on log10 scale
    log_r0_min = np.log10(prior_bounds['R0_min'])
    log_r0_max = np.log10(prior_bounds['R0_max'])
    log_r0_grid = np.linspace(log_r0_min, log_r0_max, n_r0_points)
    r0_grid = 10 ** log_r0_grid
    
    # Recovery time on linear scale
    rt_grid = np.linspace(
        prior_bounds['recovery_time_min'],
        prior_bounds['recovery_time_max'],
        n_rt_points
    )
    
    return r0_grid, rt_grid

# Setup logging
log_dir = Path(log_file).parent
log_dir.mkdir(parents=True, exist_ok=True)

with open(log_file, 'w') as log:
    log.write("="*70 + "\n")
    log.write("Generate Parameter Scan Grid\n")
    log.write("="*70 + "\n")
    log.flush()
    
    log.write(f"Loading simulation params from: {sim_params_file}\n")
    sim_params = load_config(sim_params_file)
    log.flush()
    
    # Extract prior bounds
    prior_bounds = {
        'R0_min': sim_params['R0_min'],
        'R0_max': sim_params['R0_max'],
        'recovery_time_min': sim_params['recovery_time_min'],
        'recovery_time_max': sim_params['recovery_time_max']
    }
    log.write(f"Prior bounds:\n")
    log.write(f"  R0: [{prior_bounds['R0_min']}, {prior_bounds['R0_max']}]\n")
    log.write(f"  recovery_time: [{prior_bounds['recovery_time_min']}, {prior_bounds['recovery_time_max']}]\n")
    log.flush()
    
    log.write(f"\nGrid dimensions:\n")
    log.write(f"  n_r0_points: {n_r0_points}\n")
    log.write(f"  n_rt_points: {n_rt_points}\n")
    log.write(f"  Total grid points: {n_r0_points * n_rt_points}\n")
    log.flush()
    
    log.write(f"\nCreating parameter grid...\n")
    r0_grid, rt_grid = create_parameter_grid(prior_bounds, n_r0_points, n_rt_points)
    log.flush()
    
    # Prepare metadata
    metadata = {
        'n_r0_points': n_r0_points,
        'n_rt_points': n_rt_points,
        'R0_min': prior_bounds['R0_min'],
        'R0_max': prior_bounds['R0_max'],
        'recovery_time_min': prior_bounds['recovery_time_min'],
        'recovery_time_max': prior_bounds['recovery_time_max'],
        'r0_spacing': 'log-uniform',
        'rt_spacing': 'linear'
    }
    
    log.write(f"\nGrid statistics:\n")
    log.write(f"  R0 grid:\n")
    log.write(f"    Min: {r0_grid.min():.6f}\n")
    log.write(f"    Max: {r0_grid.max():.6f}\n")
    log.write(f"    First 5: {r0_grid[:5]}\n")
    log.write(f"    Last 5: {r0_grid[-5:]}\n")
    log.write(f"  RT grid:\n")
    log.write(f"    Min: {rt_grid.min():.6f}\n")
    log.write(f"    Max: {rt_grid.max():.6f}\n")
    log.write(f"    First 5: {rt_grid[:5]}\n")
    log.write(f"    Last 5: {rt_grid[-5:]}\n")
    log.flush()
    
    # Save grid with metadata
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    np.savez(
        output_file,
        r0_grid=r0_grid,
        rt_grid=rt_grid,
        metadata=metadata
    )
    log.write(f"\nSaved parameter grid to: {output_file}\n")
    log.write(f"✓ Grid generation complete!\n")

print(f"Generated parameter scan grid: {n_r0_points} × {n_rt_points} = {n_r0_points * n_rt_points} points")
print(f"Saved to: {output_file}")
