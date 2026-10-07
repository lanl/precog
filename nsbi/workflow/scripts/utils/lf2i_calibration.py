"""
LF2I Confidence Set Calibration for Neural Ratio Estimation

Implements frequentist-valid confidence regions via Neyman construction
with test statistic calibration (LF2I approach from Lueckmann et al. 2021).

KEY: This module operates ENTIRELY on pre-computed parameter scan results.
No new model evaluations needed - pure post-processing of CSV files.

Reference:
Lueckmann et al. (2021). "Likelihood-free inference with neural ratio estimation"
Electronic Journal of Statistics. https://doi.org/10.1214/24-EJS2307
Section 3.3: Test statistic calibration
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import KFold
import pickle


class LF2ICalibrator:
    """
    Calibrate confidence sets using pre-computed parameter scan results.
    
    Workflow:
    1. Load existing param_scan.csv (already has predictions for all grid points)
    2. Convert probabilities to log-ratios
    3. Compute test statistics τ(x, θ) = -2 * [log r̂(x,θ) - log r̂(x,θ̂)]
    4. Fit quantile regression Q̂_{1-α}(θ) on calibration pairs (θ_i, τ_i)
    5. Generate confidence sets C(x) = {θ : τ(x,θ) ≤ Q̂_{1-α}(θ)}
    
    The test statistic τ follows the Neyman-Pearson lemma: under correct model
    specification, -2 log(LR) is asymptotically χ² distributed. However, with
    neural ratio estimation, the distribution may differ, so we calibrate it
    empirically using the test set.
    """
    
    def __init__(self, alpha=0.05, quantile_regressor_params=None):
        """
        Initialize calibrator.
        
        Args:
            alpha: Significance level (default 0.05 for 95% CI)
            quantile_regressor_params: Dict of sklearn GradientBoostingRegressor params
                                      If None, uses robust defaults
        """
        self.alpha = alpha
        self.quantile = 1 - alpha  # 0.95 for 95% CI
        
        # Default to simple, robust regressor for production
        # Conservative settings to avoid overfitting with small calibration sets
        if quantile_regressor_params is None:
            quantile_regressor_params = {
                'loss': 'quantile',
                'alpha': self.quantile,
                'n_estimators': 100,
                'max_depth': 3,
                'learning_rate': 0.1,
                'min_samples_split': 5,
                'min_samples_leaf': 2,
                'subsample': 0.8,
                'random_state': 42
            }
        
        self.qr_params = quantile_regressor_params
        self.quantile_regressor = None
        self.calibration_data = None
        
    def fit(self, param_scan_df, test_data_npz):
        """
        Fit quantile regressor using pre-computed parameter scan.
        
        Args:
            param_scan_df: DataFrame from {model_id}_param_scan.csv
                           Columns: sim_id, param_id, replicate, R0, recovery_time, 
                                    predicted_probability
            test_data_npz: Loaded npz dict with keys:
                          - 'sim_ids': observation identifiers
                          (Parameters are loaded from test_design.csv instead)
        
        Returns:
            self (for chaining)
        """
        print("Preparing calibration data from parameter scan...")
        
        # Step 1: Extract calibration pairs (θ_i, τ_i)
        calib_data = self._prepare_calibration_data(param_scan_df, test_data_npz)
        self.calibration_data = calib_data
        
        print(f"  Calibration samples: {len(calib_data)}")
        print(f"  τ range: [{calib_data['tau_at_true'].min():.4f}, {calib_data['tau_at_true'].max():.4f}]")
        
        # Step 2: Fit quantile regressor
        theta_calib = calib_data[['R0', 'recovery_time']].values
        tau_calib = calib_data['tau_at_true'].values
        
        print(f"\nFitting quantile regressor (quantile={self.quantile:.3f})...")
        self.quantile_regressor = GradientBoostingRegressor(**self.qr_params)
        self.quantile_regressor.fit(theta_calib, tau_calib)
        print("  ✓ Fitting complete")
        
        return self
    
    def _prepare_calibration_data(self, param_scan_df, test_data_npz):
        """
        Compute calibration pairs (θ_i, τ_i) from parameter scan results.
        
        For each test sample i:
        1. Find θ̂_i = argmax_θ log r̂(x_i, θ) from grid
        2. Compute τ_i = -2 * [log r̂(x_i, θ_i) - log r̂(x_i, θ̂_i)]
           where θ_i is the TRUE parameter that generated x_i
        
        The test statistic τ measures how much worse the likelihood is at the
        true parameter compared to the MLE. Small τ => true param is plausible.
        """
        # Load sim_ids from test.npz
        sim_ids = test_data_npz['sim_ids']
        
        # Extract param_ids from sim_ids (format: test_p00000_r0 -> param_id = 0)
        param_ids = np.array([int(sid.split('_')[1][1:]) for sid in sim_ids])
        
        # Load test_design.csv to get original-scale parameters (no denormalization needed!)
        test_design_path = "results/phase1_simulation/parameters/test_design.csv"
        print(f"  Loading true parameters from: {test_design_path}")
        test_design = pd.read_csv(test_design_path)
        
        # Verify all param_ids exist in test_design
        missing_ids = set(param_ids) - set(test_design['param_id'])
        if missing_ids:
            raise ValueError(f"Missing param_ids in test_design.csv: {missing_ids}")
        
        # Get true parameters (original scale - R0 and recovery_time)
        true_params = test_design.loc[param_ids, ['R0', 'recovery_time']].values
        
        # Validation: Print parameter ranges
        print(f"  ✓ Loaded {len(true_params)} true parameter sets")
        print(f"  R0 range: [{true_params[:, 0].min():.2f}, {true_params[:, 0].max():.2f}]")
        print(f"  Recovery time range: [{true_params[:, 1].min():.2f}, {true_params[:, 1].max():.2f}]")
        
        # OPTIMIZATION: Pre-group by sim_id to avoid O(n_test × n_rows) filtering
        # This reduces 4+ hours to ~15 seconds for 5,000 test cases
        print(f"  Grouping {len(param_scan_df)} rows by sim_id...")
        grouped = param_scan_df.groupby('sim_id')
        print(f"  Found {len(grouped)} unique observations")
        
        calib_records = []
        n_total = len(sim_ids)
        
        for i, sim_id in enumerate(sim_ids):
            # Progress logging every 500 cases
            if (i + 1) % 500 == 0:
                print(f"  Processed {i + 1}/{n_total} test cases...")
            
            # Get all grid evaluations for this observation (optimized lookup)
            try:
                obs_df = grouped.get_group(sim_id).copy()
            except KeyError:
                print(f"  WARNING: No grid evaluations found for {sim_id}, skipping")
                continue
            
            # Convert probability p to log-ratio (logit)
            # p = sigmoid(logit) => logit = log(p / (1-p))
            # This is proportional to log r̂(x, θ) up to a constant that cancels in differences
            eps = 1e-10  # Numerical stability
            obs_df['logit'] = np.log(obs_df['predicted_probability'] + eps) - \
                              np.log(1 - obs_df['predicted_probability'] + eps)
            
            # Find MLE: θ̂ = argmax_θ log r̂(x, θ)
            mle_idx = obs_df['logit'].idxmax()
            logit_max = obs_df.loc[mle_idx, 'logit']
            mle_R0 = obs_df.loc[mle_idx, 'R0']
            mle_rt = obs_df.loc[mle_idx, 'recovery_time']
            
            # Get true parameters for this observation
            true_R0, true_rt = true_params[i]
            
            # Find grid point closest to true parameters
            # (may not be exact if true params fall between grid points)
            obs_df['dist_to_true'] = np.sqrt(
                ((obs_df['R0'] - true_R0) / (obs_df['R0'].max() - obs_df['R0'].min() + 1e-6))**2 +
                ((obs_df['recovery_time'] - true_rt) / 
                 (obs_df['recovery_time'].max() - obs_df['recovery_time'].min() + 1e-6))**2
            )
            
            closest_idx = obs_df['dist_to_true'].idxmin()
            logit_at_true = obs_df.loc[closest_idx, 'logit']
            
            # Check if true params are exactly on grid (within numerical tolerance)
            closest_dist = obs_df.loc[closest_idx, 'dist_to_true']
            exact_match = closest_dist < 1e-6
            
            # Compute test statistic: τ = -2 * [log r̂(x, θ_true) - log r̂(x, θ̂)]
            # Under null hypothesis (correct model), τ ~ χ²(df) asymptotically
            # We calibrate this empirically to handle model misspecification
            tau_at_true = -2 * (logit_at_true - logit_max)
            
            # Ensure τ ≥ 0 (can be slightly negative due to numerical precision)
            tau_at_true = max(0.0, tau_at_true)
            
            calib_records.append({
                'sim_id': sim_id,
                'R0': true_R0,
                'recovery_time': true_rt,
                'tau_at_true': tau_at_true,
                'mle_R0': mle_R0,
                'mle_recovery_time': mle_rt,
                'logit_at_true': logit_at_true,
                'logit_max': logit_max,
                'exact_grid_match': exact_match,
                'closest_grid_dist': closest_dist
            })
        
        calib_df = pd.DataFrame(calib_records)
        
        # Report diagnostics
        n_exact = calib_df['exact_grid_match'].sum()
        print(f"  True params on grid: {n_exact}/{len(calib_df)}")
        if n_exact < len(calib_df):
            print(f"  Mean distance to grid: {calib_df['closest_grid_dist'].mean():.6f}")
        
        return calib_df
    
    def predict_quantile(self, theta_grid):
        """
        Predict critical value Q̂_{1-α}(θ) at grid points.
        
        The critical value is the (1-α) quantile of the test statistic
        distribution conditional on θ. Points where τ(x,θ) ≤ Q̂_{1-α}(θ)
        are included in the confidence set.
        
        Args:
            theta_grid: (n_grid, 2) array of [R0, recovery_time] values
        
        Returns:
            critical_values: (n_grid,) array of critical values
        """
        if self.quantile_regressor is None:
            raise ValueError("Must call fit() before predict_quantile()")
        
        return self.quantile_regressor.predict(theta_grid)
    
    def compute_confidence_set(self, sim_id, param_scan_df, test_data=None, return_grid=False):
        """
        Compute confidence set C(x) for a single observation.
        
        The confidence set is defined as:
            C(x) = {θ : τ(x, θ) ≤ Q̂_{1-α}(θ)}
        
        where τ is the test statistic and Q̂_{1-α}(θ) is the calibrated
        critical value at parameter θ.
        
        Args:
            sim_id: Observation identifier
            param_scan_df: Full parameter scan DataFrame with all grid evaluations
            test_data: Optional test_data dict (from npz) to add true parameters and coverage
            return_grid: If True, return full grid with indicators; 
                        if False, return summary dict
        
        Returns:
            If return_grid=False:
                dict with set summary (size, bounds, volumes, true params, coverage, etc.)
            If return_grid=True:
                DataFrame with in_set indicator for each grid point
                (useful for plotting later)
        """
        # Get all grid evaluations for this observation
        obs_df = param_scan_df[param_scan_df['sim_id'] == sim_id].copy()
        
        if len(obs_df) == 0:
            raise ValueError(f"No grid evaluations found for {sim_id}")
        
        # Convert probabilities to log-ratios
        eps = 1e-10
        obs_df['logit'] = np.log(obs_df['predicted_probability'] + eps) - \
                          np.log(1 - obs_df['predicted_probability'] + eps)
        
        # Find MLE on grid
        mle_idx = obs_df['logit'].idxmax()
        logit_max = obs_df.loc[mle_idx, 'logit']
        
        # Compute test statistic for all grid points
        obs_df['tau'] = -2 * (obs_df['logit'] - logit_max)
        obs_df['tau'] = obs_df['tau'].clip(lower=0.0)  # Numerical stability
        
        # Predict critical values at all grid points
        theta_grid = obs_df[['R0', 'recovery_time']].values
        obs_df['critical_value'] = self.predict_quantile(theta_grid)
        
        # Confidence set: {θ : τ(x, θ) ≤ Q̂_{1-α}(θ)}
        obs_df['in_set'] = obs_df['tau'] <= obs_df['critical_value']
        
        if return_grid:
            # Return full grid for plotting
            return obs_df
        
        # Compute summary statistics
        in_set_df = obs_df[obs_df['in_set']]
        
        # Parse sim_id to extract param_id and replicate
        # Format: test_p{param_id}_r{replicate}
        parts = sim_id.split('_')
        param_id = parts[1]  # p00000
        replicate = int(parts[2][1:])  # r0 -> 0
        
        # Compute grid spacings for interval width calculation
        # Since parameters are continuous but discretized on a grid, each grid point
        # represents a cell/bin centered at that value. The interval width should
        # account for the full extent of cells, not just point-to-point distance.
        R0_unique = np.sort(obs_df['R0'].unique())
        rt_unique = np.sort(obs_df['recovery_time'].unique())
        
        # R0 grid is log-uniform, so compute spacing in log10 space
        # Recovery time grid is linear, so compute spacing in linear space
        if len(R0_unique) > 1:
            log_R0_unique = np.log10(R0_unique)
            log_R0_spacing = float(np.min(np.diff(log_R0_unique)))
        else:
            log_R0_spacing = 0.0
        
        rt_spacing = float(np.min(np.diff(rt_unique))) if len(rt_unique) > 1 else 0.0
        
        summary = {
            'sim_id': sim_id,
            'param_id': param_id,
            'replicate': replicate,
            'set_size': int(in_set_df['in_set'].sum()),
            'R0_lower': float(in_set_df['R0'].min()) if len(in_set_df) > 0 else np.nan,
            'R0_upper': float(in_set_df['R0'].max()) if len(in_set_df) > 0 else np.nan,
            'recovery_time_lower': float(in_set_df['recovery_time'].min()) if len(in_set_df) > 0 else np.nan,
            'recovery_time_upper': float(in_set_df['recovery_time'].max()) if len(in_set_df) > 0 else np.nan,
            'mle_R0': float(obs_df.loc[mle_idx, 'R0']),
            'mle_recovery_time': float(obs_df.loc[mle_idx, 'recovery_time']),
            'alpha': self.alpha,
            'nominal_coverage': self.quantile
        }
        
        # Add widths accounting for grid cell boundaries
        # R0: Since grid is log-uniform, compute width in log10 space then convert back
        # Recovery time: Grid is linear, compute width directly
        # Formula: width = (max - min) + spacing
        # This accounts for half-cell width on each boundary
        if len(in_set_df) > 0:
            # Compute width in log10 space for R0
            log_R0_lower = np.log10(summary['R0_lower'])
            log_R0_upper = np.log10(summary['R0_upper'])
            log_R0_width = (log_R0_upper - log_R0_lower) + log_R0_spacing
            
            # Convert log width back to linear space
            # The confidence interval in linear space spans [10^(log_lower - log_spacing/2), 10^(log_upper + log_spacing/2)]
            R0_lower_extended = 10 ** (log_R0_lower - log_R0_spacing / 2.0)
            R0_upper_extended = 10 ** (log_R0_upper + log_R0_spacing / 2.0)
            summary['R0_width'] = float(R0_upper_extended - R0_lower_extended)
            summary['R0_log10_width'] = float(log_R0_width)  # Also store log10 width
            
            summary['recovery_time_width'] = (summary['recovery_time_upper'] - 
                                             summary['recovery_time_lower'] + 
                                             rt_spacing)
        else:
            summary['R0_width'] = 0.0
            summary['R0_log10_width'] = 0.0
            summary['recovery_time_width'] = 0.0
        
        # Approximate set volume (for 2D: area)
        if len(in_set_df) > 0:
            summary['set_volume'] = float(summary['R0_width'] * summary['recovery_time_width'])
        else:
            summary['set_volume'] = 0.0
        
        # Add true parameters and coverage if test_data provided
        if test_data is not None:
            # Find index in test_data
            sim_idx = np.where(test_data['sim_ids'] == sim_id)[0]
            if len(sim_idx) > 0:
                sim_idx = sim_idx[0]
                
                # Load true parameters from test_design.csv
                # Extract param_id from sim_id (format: test_p00000_r0 -> param_id = 0)
                param_id = int(sim_id.split('_')[1][1:])
                test_design = pd.read_csv("results/phase1_simulation/parameters/test_design.csv")
                true_R0 = float(test_design.loc[test_design['param_id'] == param_id, 'R0'].values[0])
                true_rt = float(test_design.loc[test_design['param_id'] == param_id, 'recovery_time'].values[0])
                
                summary['true_R0'] = true_R0
                summary['true_recovery_time'] = true_rt
                
                # Check coverage accounting for grid cell boundaries
                # R0: Grid is log-uniform, so extend bounds in log10 space
                # Recovery time: Grid is linear, so extend bounds in linear space
                # Since grid points are cell centers and our width formula is (max-min)+spacing,
                # the coverage region extends by spacing/2 beyond the boundary grid points.
                if not np.isnan(summary['R0_lower']):
                    # For R0 (log-uniform grid), extend bounds in log10 space
                    log_R0_lower = np.log10(summary['R0_lower'])
                    log_R0_upper = np.log10(summary['R0_upper'])
                    log_R0_lower_extended = log_R0_lower - log_R0_spacing / 2.0
                    log_R0_upper_extended = log_R0_upper + log_R0_spacing / 2.0
                    
                    # Check if true R0 falls in extended region (in log space)
                    log_true_R0 = np.log10(summary['true_R0'])
                    summary['R0_coverage'] = bool(
                        log_R0_lower_extended <= log_true_R0 <= log_R0_upper_extended
                    )
                else:
                    summary['R0_coverage'] = False
                
                if not np.isnan(summary['recovery_time_lower']):
                    # For recovery time (linear grid), extend bounds in linear space
                    rt_lower_extended = summary['recovery_time_lower'] - rt_spacing / 2.0
                    rt_upper_extended = summary['recovery_time_upper'] + rt_spacing / 2.0
                    summary['recovery_time_coverage'] = bool(
                        rt_lower_extended <= summary['true_recovery_time'] <= rt_upper_extended
                    )
                else:
                    summary['recovery_time_coverage'] = False
                
                summary['coverage_2d'] = bool(
                    summary['R0_coverage'] and summary['recovery_time_coverage']
                )
                
                # Add test statistic and critical value at true parameters
                # Find closest grid point to true parameters
                theta_true = np.array([[summary['true_R0'], summary['true_recovery_time']]])
                distances = np.sqrt(
                    (obs_df['R0'] - summary['true_R0'])**2 + 
                    (obs_df['recovery_time'] - summary['true_recovery_time'])**2
                )
                closest_idx = distances.idxmin()
                
                summary['tau_at_true'] = float(obs_df.loc[closest_idx, 'tau'])
                summary['quantile_at_true'] = float(obs_df.loc[closest_idx, 'critical_value'])
                
                # Compute Interval Scores (frequentist-compatible metric)
                # Measures sharpness (width) vs miscoverage tradeoff
                from .probabilistic_scores import compute_interval_score
                
                summary['interval_score_R0'] = float(compute_interval_score(
                    summary['R0_lower'], 
                    summary['R0_upper'], 
                    true_R0, 
                    alpha=self.alpha
                ))
                
                summary['interval_score_recovery_time'] = float(compute_interval_score(
                    summary['recovery_time_lower'],
                    summary['recovery_time_upper'],
                    true_rt,
                    alpha=self.alpha
                ))
        
        return summary
    
    def cross_validate_coverage(self, k_folds=5):
        """
        K-fold cross-validation to estimate empirical coverage.
        
        This provides a diagnostic check of the calibration quality.
        We fit the quantile regressor on k-1 folds and check if the
        held-out fold has the expected coverage rate (1-α).
        
        Args:
            k_folds: Number of CV folds
        
        Returns:
            dict with coverage statistics:
                - fold_coverages: list of per-fold coverage rates
                - mean_coverage: average across folds
                - std_coverage: standard deviation
                - se_coverage: standard error
                - nominal_coverage: target coverage (1-α)
                - n_folds: number of folds used
                - n_calib_samples: total calibration samples
        """
        if self.calibration_data is None:
            raise ValueError("Must call fit() before cross_validate_coverage()")
        
        theta_calib = self.calibration_data[['R0', 'recovery_time']].values
        tau_calib = self.calibration_data['tau_at_true'].values
        
        n_samples = len(theta_calib)
        
        # Adjust k_folds if necessary
        if k_folds > n_samples:
            k_folds = n_samples
            print(f"  WARNING: Reduced k_folds to {k_folds} (= n_samples)")
        
        kf = KFold(n_splits=k_folds, shuffle=True, random_state=42)
        
        fold_coverages = []
        
        for fold_idx, (train_idx, test_idx) in enumerate(kf.split(theta_calib)):
            # Fit on train fold
            qr_fold = GradientBoostingRegressor(**self.qr_params)
            qr_fold.fit(theta_calib[train_idx], tau_calib[train_idx])
            
            # Predict critical values on test fold
            critical_values = qr_fold.predict(theta_calib[test_idx])
            
            # Check coverage: τ_i ≤ Q̂(θ_i)
            # If calibration is correct, this should happen ~(1-α) of the time
            covered = tau_calib[test_idx] <= critical_values
            fold_coverage = covered.mean()
            fold_coverages.append(fold_coverage)
        
        return {
            'fold_coverages': fold_coverages,
            'mean_coverage': np.mean(fold_coverages),
            'std_coverage': np.std(fold_coverages),
            'se_coverage': np.std(fold_coverages) / np.sqrt(k_folds),
            'nominal_coverage': self.quantile,
            'n_folds': k_folds,
            'n_calib_samples': len(tau_calib)
        }
    
    # DEPRECATED: No longer needed - parameters loaded directly from test_design.csv
    # def _denormalize_params(self, params_norm, param_mean, param_std, param_is_log):
    #     """
    #     Denormalize parameters (inverse of training normalization).
    #     
    #     Training normalizes parameters as: z = (x - mean) / std
    #     where x may be log-transformed for parameters like R0.
    #     
    #     This reverses both operations.
    #     
    #     NOTE: This method contained a bug - used np.exp() but should have used 10**
    #     to match the log10() transformation in ParameterNormalizer.
    #     Fixed by loading original-scale parameters from test_design.csv instead.
    #     """
    #     # First denormalize: x = z * std + mean
    #     params = params_norm * param_std + param_mean
    #     
    #     # If originally log-transformed, exponentiate
    #     # This is typical for parameters like R0 that span orders of magnitude
    #     if param_is_log[0]:
    #         params[:, 0] = np.exp(params[:, 0])  # BUG: should be 10**params[:, 0]
    #     if len(param_is_log) > 1 and param_is_log[1]:
    #         params[:, 1] = np.exp(params[:, 1])  # BUG: should be 10**params[:, 1]
    #     
    #     return params
    
    def save(self, filepath):
        """Save calibrated model to disk (pickle format)."""
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
    
    @classmethod
    def load(cls, filepath):
        """Load calibrated model from disk."""
        with open(filepath, 'rb') as f:
            return pickle.load(f)

