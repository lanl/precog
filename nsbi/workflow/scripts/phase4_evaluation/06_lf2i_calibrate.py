"""
Calibrate LF2I quantile regressor using existing parameter scan.

INPUTS (all pre-existing):
- param_scan.csv: Pre-computed grid evaluations from Phase 4
- test_data.npz: True parameters from test set

OUTPUTS:
- lf2i_calibration.pkl: Fitted quantile regressor (for inference and plotting)
- lf2i_diagnostics.yaml: Cross-validation results

NO NEW SIMULATIONS - Pure post-processing of existing results.
"""

import sys
import yaml
import numpy as np
import pandas as pd
from pathlib import Path

# Add scripts directory to path (use absolute path for robustness)
scripts_path = str(Path(__file__).parent.parent)
sys.path.insert(0, scripts_path)
from utils.lf2i_calibration import LF2ICalibrator

# Snakemake inputs
param_scan_file = snakemake.input.param_scan
test_data_file = snakemake.input.test_data
calibration_output = snakemake.output.calibration
diagnostics_output = snakemake.output.diagnostics
alpha = snakemake.params.alpha
k_folds = snakemake.params.k_folds
log_file = snakemake.log[0]

# Setup logging
Path(log_file).parent.mkdir(parents=True, exist_ok=True)

with open(log_file, 'w') as log:
    log.write("="*70 + "\n")
    log.write("LF2I Calibration - Fit Quantile Regressor\n")
    log.write("="*70 + "\n\n")
    log.flush()
    
    # Load data
    log.write(f"Loading parameter scan: {param_scan_file}\n")
    param_scan_df = pd.read_csv(param_scan_file)
    log.write(f"  Rows: {len(param_scan_df)}\n")
    log.write(f"  Unique sim_ids: {param_scan_df['sim_id'].nunique()}\n")
    log.write(f"  Grid points per obs: {len(param_scan_df) // param_scan_df['sim_id'].nunique()}\n\n")
    log.flush()
    
    log.write(f"Loading test data: {test_data_file}\n")
    test_data = np.load(test_data_file, allow_pickle=True)
    log.write(f"  Test samples: {len(test_data['sim_ids'])}\n")
    log.write(f"  Parameters shape: {test_data['parameters'].shape}\n\n")
    log.flush()
    
    # Initialize calibrator
    log.write(f"Initializing LF2I calibrator:\n")
    log.write(f"  Alpha: {alpha}\n")
    log.write(f"  Nominal coverage: {(1-alpha):.1%}\n\n")
    log.flush()
    
    calibrator = LF2ICalibrator(alpha=alpha)
    
    # Fit
    log.write("="*70 + "\n")
    log.write("FITTING QUANTILE REGRESSOR\n")
    log.write("="*70 + "\n\n")
    log.flush()
    
    # Redirect print to log
    import io
    from contextlib import redirect_stdout
    
    string_buffer = io.StringIO()
    with redirect_stdout(string_buffer):
        calibrator.fit(param_scan_df, test_data)
    
    log.write(string_buffer.getvalue())
    log.write("\n")
    log.flush()
    
    log.write(f"✓ Calibration complete\n")
    log.write(f"  Calibration samples: {len(calibrator.calibration_data)}\n")
    log.write(f"  Regressor type: {type(calibrator.quantile_regressor).__name__}\n\n")
    log.flush()
    
    # Cross-validation
    log.write("="*70 + "\n")
    log.write(f"CROSS-VALIDATION ({k_folds}-fold)\n")
    log.write("="*70 + "\n\n")
    log.flush()
    
    cv_results = calibrator.cross_validate_coverage(k_folds=k_folds)
    
    log.write(f"Fold coverages:\n")
    for i, cov in enumerate(cv_results['fold_coverages']):
        log.write(f"  Fold {i+1}: {cov:.3f}\n")
    log.write(f"\n")
    log.write(f"Summary:\n")
    log.write(f"  Mean coverage: {cv_results['mean_coverage']:.3f}\n")
    log.write(f"  Std dev: {cv_results['std_coverage']:.3f}\n")
    log.write(f"  Std error: {cv_results['se_coverage']:.3f}\n")
    log.write(f"  Nominal coverage: {cv_results['nominal_coverage']:.3f}\n")
    log.write(f"  Deviation: {(cv_results['mean_coverage'] - cv_results['nominal_coverage']):.3f}\n\n")
    log.flush()
    
    # Check if coverage is reasonable
    coverage_deviation = abs(cv_results['mean_coverage'] - cv_results['nominal_coverage'])
    expected_se = np.sqrt(cv_results['nominal_coverage'] * (1 - cv_results['nominal_coverage']) / 
                          cv_results['n_calib_samples'])
    
    if coverage_deviation > 3 * expected_se:
        log.write("⚠ WARNING: Coverage significantly deviates from nominal level\n")
        log.write(f"  Deviation: {coverage_deviation:.3f}\n")
        log.write(f"  Expected SE: {expected_se:.3f}\n")
        log.write("  This may indicate model misspecification or small sample size\n\n")
    else:
        log.write(f"✓ Coverage within expected range (3 SE = {3*expected_se:.3f})\n\n")
    log.flush()
    
    # Save
    log.write("="*70 + "\n")
    log.write("SAVING RESULTS\n")
    log.write("="*70 + "\n\n")
    log.flush()
    
    Path(calibration_output).parent.mkdir(parents=True, exist_ok=True)
    calibrator.save(calibration_output)
    log.write(f"✓ Saved calibration to: {calibration_output}\n")
    log.flush()
    
    with open(diagnostics_output, 'w') as f:
        yaml.dump(cv_results, f, default_flow_style=False)
    log.write(f"✓ Saved diagnostics to: {diagnostics_output}\n\n")
    log.flush()
    
    log.write("="*70 + "\n")
    log.write("CALIBRATION COMPLETE!\n")
    log.write("="*70 + "\n")
    log.flush()

print(f"✓ LF2I calibration complete: {calibration_output}")
print(f"  Calibration samples: {len(calibrator.calibration_data)}")
print(f"  Mean CV coverage: {cv_results['mean_coverage']:.3f} (nominal: {cv_results['nominal_coverage']:.3f})")
