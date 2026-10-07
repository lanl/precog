"""
Highest Posterior Density (HPD) Credible Regions for Neural Ratio Estimation

Implements Bayesian credible intervals/regions using the posterior distribution
implied by the neural ratio estimator (assuming flat prior).

KEY DIFFERENCES FROM LF2I:
- Bayesian credible regions (not frequentist confidence sets)
- No calibration step - directly uses posterior from neural estimator
- Simpler and more interpretable, but coverage depends on estimator quality

WORKFLOW:
1. Convert predicted probabilities to likelihood ratios: r = p/(1-p)
2. Normalize to posterior (assuming flat prior): π(θ|x) ∝ r(x,θ)
3. For marginals: integrate (sum) over other parameters
4. Select HPD: smallest region containing (1-α) posterior mass
5. Check coverage: true parameter in HPD region?

Reference:
Chen, M.-H., & Shao, Q.-M. (1999). Monte Carlo Estimation of Bayesian 
Credible and HPD Intervals. Journal of Computational and Graphical Statistics.
"""

import numpy as np
import pandas as pd
import pickle


class HPDCalibrator:
    """
    Compute HPD credible regions from pre-computed parameter scan.
    
    Unlike LF2I, this is a pure post-processing method that doesn't require
    calibration - it directly uses the posterior distribution from the neural
    ratio estimator.
    
    Supports both:
    - Marginal 1D credible intervals (for each parameter separately)
    - Joint 2D credible regions (both parameters together)
    """
    
    def __init__(self, alpha=0.05):
        """
        Initialize HPD calibrator.
        
        Args:
            alpha: Significance level (default 0.05 for 95% credible intervals)
        """
        self.alpha = alpha
        self.credible_level = 1 - alpha  # 0.95 for 95% credible interval
    
    def compute_posterior_from_grid(self, obs_df):
        """
        Convert predicted probabilities to normalized posterior distribution.
        
        Assumes flat (uniform) prior, so posterior ∝ likelihood ratio.
        
        Args:
            obs_df: DataFrame with 'predicted_probability' column
        
        Returns:
            posterior: Normalized posterior probabilities (sums to 1)
        """
        eps = 1e-10  # Numerical stability
        
        # Convert probability to likelihood ratio: r = p/(1-p)
        probs = obs_df['predicted_probability'].values
        likelihood_ratio = probs / (1 - probs + eps)
        
        # Normalize to get posterior (assuming flat prior)
        posterior = likelihood_ratio / likelihood_ratio.sum()
        
        return posterior
    
    def compute_marginal_hpd_1d(self, obs_df, param_name):
        """
        Compute 1D marginal HPD credible interval for a single parameter.
        
        Marginalizes over the other parameter by summing posterior mass.
        
        Args:
            obs_df: DataFrame with columns [param_name, 'posterior']
            param_name: Name of parameter ('R0' or 'recovery_time')
        
        Returns:
            hpd_values: Array of parameter values in the HPD region
            marginal_posterior: Series with marginal posterior for each unique value
        """
        # Marginalize: sum posterior over other parameter
        marginal = obs_df.groupby(param_name)['posterior'].sum()
        
        # Sort by posterior density (descending)
        sorted_marginal = marginal.sort_values(ascending=False)
        
        # Compute cumulative mass
        cumsum = sorted_marginal.cumsum()
        
        # HPD region: smallest set containing credible_level mass
        in_hpd = cumsum <= self.credible_level
        hpd_values = sorted_marginal[in_hpd].index.values
        
        return hpd_values, marginal
    
    def compute_joint_hpd_2d(self, obs_df):
        """
        Compute 2D joint HPD credible region.
        
        Selects grid points with highest posterior density until
        credible_level mass is accumulated.
        
        Args:
            obs_df: DataFrame with columns ['R0', 'recovery_time', 'posterior']
        
        Returns:
            hpd_df: DataFrame of grid points in the HPD region
        """
        # Sort all grid points by posterior density (descending)
        sorted_df = obs_df.sort_values('posterior', ascending=False).copy()
        
        # Compute cumulative mass
        sorted_df['cumsum'] = sorted_df['posterior'].cumsum()
        
        # HPD region: points until cumulative mass reaches credible_level
        sorted_df['in_hpd'] = sorted_df['cumsum'] <= self.credible_level
        
        # Return points in HPD
        hpd_df = sorted_df[sorted_df['in_hpd']].copy()
        
        return hpd_df
    
    def compute_joint_mode(self, obs_df):
        """
        Compute joint mode (MAP estimate).
        
        Returns (R0, rt) at maximum joint posterior density.
        
        Args:
            obs_df: DataFrame with 'posterior', 'R0', 'recovery_time' columns
        
        Returns:
            tuple: (R0_mode, rt_mode)
        """
        mode_idx = obs_df['posterior'].idxmax()
        return (
            float(obs_df.loc[mode_idx, 'R0']),
            float(obs_df.loc[mode_idx, 'recovery_time'])
        )
    
    def compute_marginal_mode(self, marginal_posterior):
        """
        Compute marginal mode from marginal posterior.
        
        Args:
            marginal_posterior: Series with parameter values as index, 
                               posterior probabilities as values
        
        Returns:
            float: Parameter value with highest marginal posterior
        """
        return float(marginal_posterior.idxmax())
    
    def compute_marginal_mean(self, marginal_posterior):
        """
        Compute marginal mean (expected value).
        
        Optimal for minimizing component-wise squared error.
        
        Args:
            marginal_posterior: Series with parameter values as index,
                               posterior probabilities as values
        
        Returns:
            float: E[θ|x] = Σ θ · p(θ|x)
        """
        param_values = marginal_posterior.index.values
        probabilities = marginal_posterior.values
        return float(np.sum(param_values * probabilities))
    
    def compute_marginal_median(self, marginal_posterior):
        """
        Compute marginal median (50th percentile).
        
        Optimal for minimizing component-wise absolute error (MAE).
        
        Args:
            marginal_posterior: Series with parameter values as index,
                               posterior probabilities as values
        
        Returns:
            float: Value where cumulative probability = 0.5
        """
        # Sort by parameter value (not by density)
        sorted_marginal = marginal_posterior.sort_index()
        
        # Compute cumulative distribution function
        cdf = sorted_marginal.cumsum()
        
        # Find first value where CDF >= 0.5
        median_idx = (cdf >= 0.5).idxmax()
        
        return float(median_idx)
    
    def compute_all_point_estimates(self, obs_df, marginal_R0, marginal_rt):
        """
        Compute all point estimate types.
        
        Args:
            obs_df: DataFrame with joint posterior
            marginal_R0: Series with marginal posterior for R0
            marginal_rt: Series with marginal posterior for recovery_time
        
        Returns:
            dict with keys:
                - joint_mode_R0, joint_mode_rt
                - marginal_mode_R0, marginal_mode_rt
                - marginal_mean_R0, marginal_mean_rt (PRIMARY)
                - marginal_median_R0, marginal_median_rt
        """
        # Joint mode
        joint_mode_R0, joint_mode_rt = self.compute_joint_mode(obs_df)
        
        # Marginal modes
        marginal_mode_R0 = self.compute_marginal_mode(marginal_R0)
        marginal_mode_rt = self.compute_marginal_mode(marginal_rt)
        
        # Marginal means (PRIMARY ESTIMATES - optimal for MSE)
        marginal_mean_R0 = self.compute_marginal_mean(marginal_R0)
        marginal_mean_rt = self.compute_marginal_mean(marginal_rt)
        
        # Marginal medians (optimal for MAE)
        marginal_median_R0 = self.compute_marginal_median(marginal_R0)
        marginal_median_rt = self.compute_marginal_median(marginal_rt)
        
        return {
            'joint_mode_R0': joint_mode_R0,
            'joint_mode_rt': joint_mode_rt,
            'marginal_mode_R0': marginal_mode_R0,
            'marginal_mode_rt': marginal_mode_rt,
            'marginal_mean_R0': marginal_mean_R0,
            'marginal_mean_rt': marginal_mean_rt,
            'marginal_median_R0': marginal_median_R0,
            'marginal_median_rt': marginal_median_rt,
        }
    
    def check_coverage_1d(self, true_value, all_values, hpd_values):
        """
        Check if true parameter value is covered by 1D HPD region.
        
        Uses nearest neighbor approach: if the closest grid point to the
        true value is in the HPD, we count it as covered.
        
        Args:
            true_value: True parameter value (scalar)
            all_values: Array of ALL parameter values on the grid
            hpd_values: Array of parameter values in HPD region
        
        Returns:
            bool: True if covered, False otherwise
        """
        if len(all_values) == 0:
            return False
        
        if len(hpd_values) == 0:
            return False
        
        # Find nearest grid point to true value (from ALL grid points)
        distances = np.abs(all_values - true_value)
        nearest_value = all_values[distances.argmin()]
        
        # Check if that nearest grid point is in the HPD region
        is_covered = nearest_value in hpd_values
        
        return bool(is_covered)
    
    def check_coverage_2d(self, true_R0, true_rt, all_grid_df, hpd_df):
        """
        Check if true parameter vector is covered by 2D HPD region.
        
        Uses nearest neighbor in Euclidean distance.
        
        Args:
            true_R0: True R0 value
            true_rt: True recovery time value
            all_grid_df: DataFrame of ALL grid points (columns: R0, recovery_time)
            hpd_df: DataFrame of grid points in HPD (columns: R0, recovery_time)
        
        Returns:
            bool: True if covered, False otherwise
        """
        if len(all_grid_df) == 0:
            return False
        
        if len(hpd_df) == 0:
            return False
        
        # Find nearest grid point to true parameters (from ALL grid points)
        distances = np.sqrt(
            (all_grid_df['R0'].values - true_R0)**2 + 
            (all_grid_df['recovery_time'].values - true_rt)**2
        )
        nearest_idx = distances.argmin()
        nearest_R0 = all_grid_df.iloc[nearest_idx]['R0']
        nearest_rt = all_grid_df.iloc[nearest_idx]['recovery_time']
        
        # Check if that nearest grid point is in the HPD region
        is_covered = ((hpd_df['R0'] == nearest_R0) & 
                      (hpd_df['recovery_time'] == nearest_rt)).any()
        
        return bool(is_covered)
    
    def compute_confidence_set(self, sim_id, param_scan_df, test_data=None, return_grid=False, 
                              compute_scores=True, n_samples=1000):
        """
        Compute HPD credible regions for a single observation.
        
        Computes:
        1. Marginal 1D HPD intervals for R0 and recovery_time
        2. Joint 2D HPD region
        3. Coverage statistics if test_data provided
        4. Probabilistic scores if compute_scores=True
        
        Args:
            sim_id: Observation identifier
            param_scan_df: Full parameter scan DataFrame with all grid evaluations
            test_data: Optional test_data dict (from npz) to add true parameters and coverage
            return_grid: If True, return full grid with indicators; 
                        if False, return summary dict
            compute_scores: If True, compute CRPS, energy score, interval score, log score
            n_samples: Number of posterior samples for score computation (default: 1000)
        
        Returns:
            If return_grid=False:
                dict with set summary (size, bounds, coverage, scores, etc.)
            If return_grid=True:
                DataFrame with in_set indicator for each grid point
        """
        # Get all grid evaluations for this observation
        obs_df = param_scan_df[param_scan_df['sim_id'] == sim_id].copy()
        
        if len(obs_df) == 0:
            raise ValueError(f"No grid evaluations found for {sim_id}")
        
        # Compute posterior distribution
        obs_df['posterior'] = self.compute_posterior_from_grid(obs_df)
        
        # Compute marginal HPD for R0
        hpd_R0_values, marginal_R0 = self.compute_marginal_hpd_1d(obs_df, 'R0')
        
        # Compute marginal HPD for recovery_time
        hpd_rt_values, marginal_rt = self.compute_marginal_hpd_1d(obs_df, 'recovery_time')
        
        # Compute joint 2D HPD
        hpd_2d_df = self.compute_joint_hpd_2d(obs_df)
        
        # Compute all point estimates
        point_estimates = self.compute_all_point_estimates(obs_df, marginal_R0, marginal_rt)
        
        # Mark points in marginal HPDs and joint HPD
        obs_df['in_R0_hpd'] = obs_df['R0'].isin(hpd_R0_values)
        obs_df['in_rt_hpd'] = obs_df['recovery_time'].isin(hpd_rt_values)
        obs_df['in_2d_hpd'] = obs_df.index.isin(hpd_2d_df.index)
        
        if return_grid:
            # Return full grid for plotting
            return obs_df
        
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
        
        # Compute summary statistics with standardized column names
        summary = {
            'sim_id': sim_id,
            'param_id': param_id,
            'replicate': replicate,
            'method': 'hpd',
            'alpha': self.alpha,
            'nominal_coverage': self.credible_level,
            'set_size': int(len(hpd_2d_df)),  # Number of grid points in 2D HPD
            'R0_lower_bound': float(hpd_R0_values.min()) if len(hpd_R0_values) > 0 else np.nan,
            'R0_upper_bound': float(hpd_R0_values.max()) if len(hpd_R0_values) > 0 else np.nan,
            'D_lower_bound': float(hpd_rt_values.min()) if len(hpd_rt_values) > 0 else np.nan,
            'D_upper_bound': float(hpd_rt_values.max()) if len(hpd_rt_values) > 0 else np.nan,
            
            # PRIMARY POINT ESTIMATES (marginal mean - optimal for MSE)
            'R0_point_estimate': point_estimates['marginal_mean_R0'],
            'D_point_estimate': point_estimates['marginal_mean_rt'],
            
            # JOINT ESTIMATORS
            'R0_joint_mode': point_estimates['joint_mode_R0'],
            'D_joint_mode': point_estimates['joint_mode_rt'],
            
            # MARGINAL ESTIMATORS (explicit for completeness and comparison)
            'R0_marginal_mode': point_estimates['marginal_mode_R0'],
            'D_marginal_mode': point_estimates['marginal_mode_rt'],
            'R0_marginal_mean': point_estimates['marginal_mean_R0'],
            'D_marginal_mean': point_estimates['marginal_mean_rt'],
            'R0_marginal_median': point_estimates['marginal_median_R0'],
            'D_marginal_median': point_estimates['marginal_median_rt'],
        }
        
        # Add widths accounting for grid cell boundaries
        # R0: Since grid is log-uniform, compute width in log10 space then convert back
        # Recovery time: Grid is linear, compute width directly
        # Formula: width = (max - min) + spacing
        # This accounts for half-cell width on each boundary
        if len(hpd_R0_values) > 0:
            # Compute width in log10 space for R0
            log_R0_lower = np.log10(summary['R0_lower_bound'])
            log_R0_upper = np.log10(summary['R0_upper_bound'])
            log_R0_width = (log_R0_upper - log_R0_lower) + log_R0_spacing
            
            # Convert log width back to linear space
            # The credible interval in linear space spans [10^(log_lower - log_spacing/2), 10^(log_upper + log_spacing/2)]
            R0_lower_extended = 10 ** (log_R0_lower - log_R0_spacing / 2.0)
            R0_upper_extended = 10 ** (log_R0_upper + log_R0_spacing / 2.0)
            summary['R0_interval_width'] = float(R0_upper_extended - R0_lower_extended)
            summary['R0_log10_interval_width'] = float(log_R0_width)  # Also store log10 width
        else:
            summary['R0_interval_width'] = 0.0
            summary['R0_log10_interval_width'] = 0.0
        
        if len(hpd_rt_values) > 0:
            summary['D_interval_width'] = (summary['D_upper_bound'] - 
                                          summary['D_lower_bound'] + 
                                          rt_spacing)
        else:
            summary['D_interval_width'] = 0.0
        
        # Approximate set volume (for 2D: area)
        if summary['set_size'] > 0:
            summary['set_volume'] = float(summary['R0_interval_width'] * summary['D_interval_width'])
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
                param_id_int = int(sim_id.split('_')[1][1:])
                test_design = pd.read_csv("results/phase1_simulation/parameters/test_design.csv")
                true_R0 = float(test_design.loc[test_design['param_id'] == param_id_int, 'R0'].values[0])
                true_rt = float(test_design.loc[test_design['param_id'] == param_id_int, 'recovery_time'].values[0])
                
                summary['R0_true'] = true_R0
                summary['D_true'] = true_rt
                
                # Compute errors for PRIMARY estimate (marginal mean)
                summary['R0_error'] = float(summary['R0_point_estimate'] - true_R0)
                summary['D_error'] = float(summary['D_point_estimate'] - true_rt)
                summary['R0_abs_error'] = float(abs(summary['R0_error']))
                summary['D_abs_error'] = float(abs(summary['D_error']))
                
                # Compute errors for ALL estimator types (for comparison)
                for est_type in ['joint_mode', 'marginal_mode', 'marginal_median']:
                    summary[f'R0_error_{est_type}'] = float(summary[f'R0_{est_type}'] - true_R0)
                    summary[f'D_error_{est_type}'] = float(summary[f'D_{est_type}'] - true_rt)
                    summary[f'R0_abs_error_{est_type}'] = float(abs(summary[f'R0_error_{est_type}']))
                    summary[f'D_abs_error_{est_type}'] = float(abs(summary[f'D_error_{est_type}']))
                
                # Check coverage using nearest neighbor approach
                # Need to pass all grid values to find nearest neighbor properly
                all_R0_values = obs_df['R0'].unique()
                all_rt_values = obs_df['recovery_time'].unique()
                
                summary['R0_coverage'] = bool(self.check_coverage_1d(true_R0, all_R0_values, hpd_R0_values))
                summary['D_coverage'] = bool(self.check_coverage_1d(true_rt, all_rt_values, hpd_rt_values))
                summary['coverage_2d'] = bool(self.check_coverage_2d(true_R0, true_rt, obs_df[['R0', 'recovery_time']], hpd_2d_df))
                
                # Compute probabilistic scores if requested
                if compute_scores:
                    from .probabilistic_scores import compute_all_scores
                    
                    scores = compute_all_scores(
                        grid_df=obs_df,
                        true_R0=true_R0,
                        true_rt=true_rt,
                        R0_lower=summary['R0_lower_bound'],
                        R0_upper=summary['R0_upper_bound'],
                        rt_lower=summary['D_lower_bound'],
                        rt_upper=summary['D_upper_bound'],
                        alpha=self.alpha,
                        n_samples=n_samples,
                        random_state=42
                    )
                    
                    # Rename score keys to standardized names
                    score_mapping = {
                        'crps_R0': 'R0_crps',
                        'crps_recovery_time': 'D_crps',
                        'interval_score_R0': 'R0_interval_score',
                        'interval_score_recovery_time': 'D_interval_score'
                    }
                    
                    renamed_scores = {}
                    for old_key, new_key in score_mapping.items():
                        if old_key in scores:
                            renamed_scores[new_key] = scores[old_key]
                    
                    # Add other scores unchanged
                    for key in ['energy_score', 'log_score']:
                        if key in scores:
                            renamed_scores[key] = scores[key]
                    
                    # Add scores to summary
                    summary.update(renamed_scores)
        
        return summary
    
    def assess_coverage(self, param_scan_df, test_data_npz, verbose=True, 
                       compute_scores=True, n_samples=1000):
        """
        Assess empirical coverage over entire test set.
        
        Args:
            param_scan_df: Full parameter scan DataFrame
            test_data_npz: Loaded test data npz file
            verbose: Print progress
            compute_scores: If True, compute probabilistic scores for each observation
            n_samples: Number of posterior samples for score computation
        
        Returns:
            dict with coverage statistics and scores:
                - mean_coverage_2d: empirical 2D coverage
                - mean_coverage_R0: empirical R0 marginal coverage
                - mean_coverage_rt: empirical recovery time marginal coverage
                - results_df: Full results DataFrame (with scores if compute_scores=True)
        """
        sim_ids = test_data_npz['sim_ids']
        results = []
        
        if verbose:
            print(f"Assessing coverage for {len(sim_ids)} test observations...")
        
        for i, sim_id in enumerate(sim_ids):
            if verbose and (i % 50 == 0 or i == len(sim_ids) - 1):
                print(f"  [{i+1}/{len(sim_ids)}] Processing {sim_id}")
            
            summary = self.compute_confidence_set(
                sim_id, param_scan_df, test_data=test_data_npz,
                compute_scores=compute_scores, n_samples=n_samples
            )
            results.append(summary)
        
        results_df = pd.DataFrame(results)
        
        # Compute coverage statistics
        coverage_stats = {
            'mean_coverage_2d': results_df['coverage_2d'].mean(),
            'mean_coverage_R0': results_df['R0_coverage'].mean(),
            'mean_coverage_D': results_df['D_coverage'].mean(),
            'mean_R0_width': results_df['R0_interval_width'].mean(),
            'mean_D_width': results_df['D_interval_width'].mean(),
            'mean_set_size': results_df['set_size'].mean(),
            'n_observations': len(results_df),
            'results_df': results_df
        }
        
        # Add probabilistic score summaries if computed
        if compute_scores and 'R0_crps' in results_df.columns:
            from .probabilistic_scores import aggregate_scores
            score_summary = aggregate_scores(results_df)
            coverage_stats.update(score_summary)
        
        if verbose:
            print(f"\nCoverage Statistics:")
            print(f"  2D coverage: {coverage_stats['mean_coverage_2d']:.1%}")
            print(f"  R0 coverage: {coverage_stats['mean_coverage_R0']:.1%}")
            print(f"  D coverage: {coverage_stats['mean_coverage_D']:.1%}")
            print(f"  Mean R0 width: {coverage_stats['mean_R0_width']:.3f}")
            print(f"  Mean D width: {coverage_stats['mean_D_width']:.3f}")
            print(f"  Mean set size: {coverage_stats['mean_set_size']:.1f} grid points")
            
            if compute_scores and 'R0_crps_mean' in coverage_stats:
                print(f"\nProbabilistic Scores:")
                print(f"  Mean CRPS (R0): {coverage_stats['R0_crps_mean']:.4f}")
                print(f"  Mean CRPS (D): {coverage_stats['D_crps_mean']:.4f}")
                print(f"  Mean Energy Score: {coverage_stats['energy_score_mean']:.4f}")
                print(f"  Mean Interval Score (R0): {coverage_stats['R0_interval_score_mean']:.4f}")
                print(f"  Mean Interval Score (D): {coverage_stats['D_interval_score_mean']:.4f}")
                print(f"  Mean Log Score: {coverage_stats['log_score_mean']:.2f}")
        
        return coverage_stats
    
    def save(self, filepath):
        """Save calibrator to disk (pickle format)."""
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
    
    @classmethod
    def load(cls, filepath):
        """Load calibrator from disk."""
        with open(filepath, 'rb') as f:
            return pickle.load(f)


