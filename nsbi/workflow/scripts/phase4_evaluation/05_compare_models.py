"""Compare NRE models using HPD inference results (pre-computed point estimates and scores)."""

import yaml
import pandas as pd
import numpy as np

# Access snakemake variables
hpd_region_files = snakemake.input.hpd_regions
model_ids = snakemake.params.model_ids
comparison_yaml = snakemake.output.comparison_yaml
comparison_csv = snakemake.output.comparison_csv

print("="*70)
print("Comparing NRE Models (HPD-Based)")
print("="*70)
print("\nMetrics: Point estimate errors, log scores, coverage\n")

def compute_model_metrics(hpd_df):
    """Extract comparison metrics from HPD confidence regions."""
    
    # Point estimate errors (using posterior mode)
    point_errors_r0 = np.abs(hpd_df['R0_point_estimate'] - hpd_df['R0_true'])
    point_errors_d = np.abs(hpd_df['D_point_estimate'] - hpd_df['D_true'])
    
    metrics = {
        'n_observations': len(hpd_df),
        'mean_point_error_r0': float(point_errors_r0.mean()),
        'median_point_error_r0': float(point_errors_r0.median()),
        'mean_point_error_d': float(point_errors_d.mean()),
        'median_point_error_d': float(point_errors_d.median()),
        'mean_log_score': float(hpd_df['log_score'].mean()),
        'mean_crps_r0': float(hpd_df['R0_crps'].mean()),
        'mean_crps_d': float(hpd_df['D_crps'].mean()),
        'mean_energy_score': float(hpd_df['energy_score'].mean()),
        'coverage_2d': float(hpd_df['coverage_2d'].mean()),
        'coverage_r0': float(hpd_df['R0_coverage'].mean()),
        'coverage_d': float(hpd_df['D_coverage'].mean()),
        'mean_r0_width': float(hpd_df['R0_interval_width'].mean()),
        'mean_d_width': float(hpd_df['D_interval_width'].mean())
    }
    
    return metrics

# Compute metrics for each model
print(f"Loading HPD results for {len(model_ids)} models...\n")
metrics_by_model = {}

for model_id, hpd_file in zip(model_ids, hpd_region_files):
    print(f"Processing {model_id}...")
    print(f"  Loading: {hpd_file}")
    
    hpd_df = pd.read_csv(hpd_file)
    metrics = compute_model_metrics(hpd_df)
    metrics_by_model[model_id] = metrics
    
    print(f"  Observations: {metrics['n_observations']}")
    print(f"  Median point errors - R0: {metrics['median_point_error_r0']:.4f}, D: {metrics['median_point_error_d']:.4f}")
    print(f"  Mean log score: {metrics['mean_log_score']:.2f}")
    print(f"  2D coverage: {metrics['coverage_2d']:.1%}\n")

# Display comparison
print("="*70)
print("MODEL COMPARISON SUMMARY")
print("="*70)

for model_id, metrics in metrics_by_model.items():
    print(f"\n{model_id.replace('_', ' ').title()}:")
    print(f"  Point Estimate Errors (Median):")
    print(f"    R0: {metrics['median_point_error_r0']:.4f}")
    print(f"    D: {metrics['median_point_error_d']:.4f}")
    print(f"  Probabilistic Scores (Mean):")
    print(f"    Log Score: {metrics['mean_log_score']:.2f} (higher is better)")
    print(f"    CRPS R0: {metrics['mean_crps_r0']:.4f} (lower is better)")
    print(f"    CRPS D: {metrics['mean_crps_d']:.4f} (lower is better)")
    print(f"    Energy Score: {metrics['mean_energy_score']:.4f} (lower is better)")
    print(f"  Coverage (95% credible regions):")
    print(f"    2D: {metrics['coverage_2d']:.1%}")
    print(f"    R0: {metrics['coverage_r0']:.1%}")
    print(f"    D: {metrics['coverage_d']:.1%}")
    print(f"  Interval Widths (Mean):")
    print(f"    R0: {metrics['mean_r0_width']:.3f}")
    print(f"    D: {metrics['mean_d_width']:.3f}")

# Save results
print(f"\n{'='*70}")
print("Saving comparison results...")
print("="*70 + "\n")

comparison = {f'{mid}_model': m for mid, m in metrics_by_model.items()}

with open(comparison_yaml, 'w') as f:
    yaml.dump(comparison, f, default_flow_style=False)
print(f"✓ YAML saved to: {comparison_yaml}")

# CSV with key metrics for quick comparison
df_comparison = pd.DataFrame({
    'model': list(metrics_by_model.keys()),
    'median_point_error_r0': [m['median_point_error_r0'] for m in metrics_by_model.values()],
    'median_point_error_d': [m['median_point_error_d'] for m in metrics_by_model.values()],
    'mean_log_score': [m['mean_log_score'] for m in metrics_by_model.values()],
    'mean_energy_score': [m['mean_energy_score'] for m in metrics_by_model.values()],
    'coverage_2d': [m['coverage_2d'] for m in metrics_by_model.values()]
})

# Sort by point estimate error (primary) and log score (secondary)
df_comparison = df_comparison.sort_values(['median_point_error_r0', 'mean_log_score'], 
                                         ascending=[True, False])
df_comparison.to_csv(comparison_csv, index=False)
print(f"✓ CSV saved to: {comparison_csv}")

print("\n" + "="*70)
print("✓ Model Comparison Complete")
print("="*70)
print(f"\nBest model (by R0 point error): {df_comparison.iloc[0]['model']}")

