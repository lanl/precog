# tsfeatures Benchmarking

This sub-workflow systematically benchmarks time series features from the R `tsfeatures` package to identify fast, fixed-length features for case count time series analysis.

## Overview

**Problem**: The `tsfeatures` package provides many time series features, but:
1. Some are slow for large-scale analysis (processing 10,000+ series)
2. Some return variable-length outputs that can't be used directly in neural networks
3. Not all features are appropriate for epidemic case count curves

**Solution**: Benchmark all features on case count-like time series, identify fast fixed-length ones, and provide an optimized function that computes only these features.

## Directory Structure

```
benchmark/tsfeatures/
├── README.md                   # This file
├── workflow/
│   ├── Snakefile              # Benchmark workflow
│   └── scripts/
│       └── benchmark_tsfeatures.R
├── config/
│   └── config.yaml            # Benchmark configuration
└── resources/                  # Documentation and results
```

## Benchmarking Methodology

We systematically benchmark all tsfeatures functions on progressively longer time series:
- **Series lengths tested**: 100, 500, and 1000 time points
- **Timeout protection**: Features taking >30 seconds are terminated
- **Variable-length detection**: Features that return different output lengths for different input lengths are excluded
- **Smart skipping**: Features classified as slow (>1 second) on shorter series are skipped on longer series
- **Replication**: Fast features are tested 3 times for reliable measurements

### Key Innovation: Fixed-Length Filtering

Unlike tree statistics, time series features can return:
- **Fixed-length output**: Same number of features regardless of series length (✓ Good for ML)
- **Variable-length output**: Output size depends on input length (✗ Can't use directly)

We test each feature on 3 different lengths (50, 100, 200) and exclude any with variable output.

## Usage

### Running the Benchmark

```bash
cd benchmark/tsfeatures
Rscript workflow/scripts/benchmark_tsfeatures.R
```

**Output**:
- `results/benchmark_results.csv` - Timing data for all features
- `results/fast_features.txt` - List of fast, fixed-length features
- `results/variable_length_features.txt` - Features with variable output (excluded)
- `results/fast_features.RData` - R data object with fast feature names

### Configuration

Edit `config/config.yaml` to modify:
```yaml
series_lengths: [100, 500, 1000]  # Series lengths to test
timeout: 30                        # Timeout in seconds
n_replicates: 3                    # Replicates for fast features
frequency: 1                       # Time series frequency (daily)
```

## Expected Results

### Feature Categories

| Category | Criteria | Usage |
|----------|----------|-------|
| **Fast, Fixed-Length** | < 0.1s AND fixed output | ✅ Use in workflow |
| **Moderate, Fixed-Length** | 0.1-1.0s AND fixed output | ⚠️ Consider if needed |
| **Slow** | > 1.0s | ❌ Avoid |
| **Variable-Length** | Output size varies | ❌ Can't use directly |
| **Failed/Timeout** | Errors or >30s | ❌ Avoid |

### Features Tested

Based on tsfeatures package (28 functions tested):

**Fast features (<0.1s, 25 functions → 56 scalar features)**:
- `acf_features` (6 features) - Autocorrelation features
- `arch_stat` (1) - ARCH LM statistic
- `autocorr_features` (7) - Additional autocorrelation features
- `crossing_points` (1) - Crossing points
- `dist_features` (2) - Distribution features
- `entropy` (1) - Spectral entropy
- `firstzero_ac` (1) - First zero crossing of ACF
- `flat_spots` (1) - Flat spot detection
- `heterogeneity` (4) - Heterogeneity measures
- `holt_parameters` (2) - Holt's linear trend
- `hurst` (1) - Hurst exponent
- `hw_parameters` (3) - Holt-Winters parameters
- `lumpiness` (1) - Variance of variances
- `max_kl_shift` (2) - Maximum KL shift
- `max_level_shift` (2) - Maximum level shift
- `max_var_shift` (2) - Maximum variance shift
- `nonlinearity` (1) - Nonlinearity measure
- `pacf_features` (3) - Partial autocorrelation features
- `station_features` (3) - Stationarity features
- `stl_features` (8) - STL decomposition features
- `stability` (1) - Stability measure
- `unitroot_kpss` (1) - KPSS unit root test
- `unitroot_pp` (1) - Phillips-Perron test
- `zero_proportion` (1) - Proportion of zeros

**Excluded features (3 functions)**:
- `binarize_mean` - Variable-length output
- `compengine` (1.9s) - Too slow, 16 features
- `pred_features` (1.3s) - Too slow, 3 features
- `scal_features` (0.5s) - Moderate speed, 1 feature

## Integration with Main Workflow

This benchmark workflow is **standalone** but informs the main workflow:

1. **Benchmark identifies** fast, fixed-length features
2. **Main workflow uses** these features in `workflow/scripts/utils/calc_fast_tsfeatures.R`
3. **Phase 2** extracts these features from case count trajectories
4. **Phase 3** trains NRE model using these features

## Synthetic Data Generation

The benchmark uses synthetic epidemic-like time series:
- Gamma-distributed curves (mimicking epidemic trajectories)
- Variable peak times and magnitudes
- Realistic noise
- Non-negative integer case counts (including zeros)

This ensures features are tested on data similar to actual epidemic case counts.

## Performance Expectations

For ~10,000 case count series:
- **All features**: Could take hours (if many are slow)
- **Fast features only**: ~10-30 minutes (estimated)
- **Speedup**: 10-50x faster by using only fast features

## Maintenance

**When to re-run**:
- tsfeatures package updates
- Need to test additional features
- Change in time series characteristics (e.g., longer series)

**Update frequency**: Quarterly or after major tsfeatures releases

**Dependencies**:
- R >= 4.0
- R packages: tsfeatures, forecast (tsfeatures dependency)

## Files Generated

After running the benchmark:
- `results/benchmark_results.csv` - Full timing data
- `results/fast_features.txt` - List of fast feature names (one per line)
- `results/fast_features.RData` - R object for easy loading
- `results/variable_length_features.txt` - Excluded features

## Next Steps

After benchmarking:
1. Review `results/fast_features.txt` 
2. Implement `calc_fast_tsfeatures()` function using these features
3. Integrate into main workflow (Phase 2: feature extraction)
4. Train NRE model using tsfeatures

## Historical Context

- **Benchmarking date**: September 2026
- **tsfeatures version**: Latest from CRAN
- **Platform**: macOS (results similar on Linux/Windows)
- **Use case**: Epidemic case count time series from SIR simulations
