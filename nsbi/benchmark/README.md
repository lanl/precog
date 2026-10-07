# TreeStats Benchmarking

This sub-workflow systematically benchmarks phylogenetic tree statistics from the R `treestats` package to identify fast, reliable statistics for large-scale analysis.

## Overview

**Problem**: The `treestats::calc_all_stats()` function computes 70+ statistics, but some are extremely slow for large trees (1,000-10,000 tips), making them impractical when processing tens of thousands of trees.

**Solution**: Benchmark all statistics on progressively larger trees, identify fast ones, and provide an optimized function that computes only the fast statistics.

## Directory Structure

```
benchmark/treestats/
├── README.md                   # This file
├── workflow/
│   ├── Snakefile              # Benchmark workflow
│   └── scripts/
│       └── benchmark_treestats.R
├── config/
│   └── config.yaml            # Benchmark configuration
├── resources/                  # Documentation and reference data
│   ├── SIZE_AGNOSTIC_STATS.md
│   ├── LTT_NTT_NORMALIZATION.md
│   ├── ALL_STATISTICS_REFERENCE.md
│   ├── benchmark_results.csv
│   └── size_agnostic_statistics_list.txt
└── results/                    # Generated outputs (git-ignored)
```

## Benchmarking Methodology

We systematically benchmarked all 70 statistics from the `treestats` package on progressively larger trees:
- **Tree sizes tested**: 1,000, 5,000, and 10,000 tips
- **Timeout protection**: Statistics taking >30 seconds were terminated
- **Smart skipping**: Statistics classified as slow (>1 second) on smaller trees were skipped on larger trees
- **Replication**: Fast statistics were tested 3 times for reliable measurements

## Usage

### For Users: Using Fast Statistics in Your Analysis

The main workflow uses an optimized function that computes only fast statistics:

```r
# The main workflow already uses this in:
# workflow/scripts/phase2_processing/05_extract_tree_stats_batch.R

# Load your tree
library(ape)
tree <- read.nexus("my_tree.nex")

# Compute only fast statistics (automatically done by main workflow)
# See resources/size_agnostic_statistics_list.txt for the curated list
```

The main workflow automatically uses 14 size-agnostic statistics selected from the 43 fast statistics. See `resources/SIZE_AGNOSTIC_STATS.md` for details.

### For Maintainers: Re-running the Benchmark

If the treestats package updates or you need to test new statistics:

```bash
cd benchmark/treestats
snakemake --use-conda --cores 1
```

**Output**:
- `results/benchmark_results.csv` - Timing data for all statistics
- `results/fast_statistics_list.txt` - List of fast statistics
- `results/benchmark_summary.md` - Human-readable summary

**Configuration**: Edit `config/config.yaml` to modify:
```yaml
tree_sizes: [1000, 5000, 10000]  # Tree sizes to test
timeout: 30                       # Timeout in seconds
n_replicates: 3                   # Replicates for fast stats
```

**Force re-run**:
```bash
snakemake --forceall --use-conda --cores 1
```

## Key Results

### For 10,000 Tip Trees:

| Category | Count | Definition |
|----------|-------|------------|
| **Fast** | 43 | < 0.1 seconds |
| **Moderate** | 3 | 0.1 - 1.0 seconds |
| **Slow** | 2 | > 1.0 seconds |
| **Failed/Unavailable** | 22 | Errors or dependencies missing |

### Performance Comparison

**Single tree (10,000 tips):**
- All successful statistics: **~5.03 seconds**
- Fast statistics only: **~0.02 seconds**
- **Speedup: 279x faster**

### Time Projections for Large Datasets

| Number of Trees | All Stats | Fast Stats Only | Time Saved |
|-----------------|-----------|-----------------|------------|
| 10,000 | 14.0 hours | 0.1 hours | 13.9 hours |
| 50,000 | 69.8 hours | 0.3 hours | 69.5 hours |
| 100,000 | 139.6 hours | 0.5 hours | 139.1 hours |

## Fast Statistics (43 total)

These statistics complete in < 0.1 seconds for 10,000 tip trees:

### Fastest (< 0.001 sec):
1. `number_of_lineages` - 0.0000 sec
2. `mean_branch_length` - 0.0000 sec
3. `var_branch_length` - 0.0001 sec
4. `treeness` - 0.0002 sec
5. `max_width` - 0.0002 sec

### Branch Length Statistics:
- `mean_branch_length`
- `mean_branch_length_int`
- `mean_branch_length_ext`
- `var_branch_length`
- `var_branch_length_int`
- `var_branch_length_ext`

### Balance/Imbalance Metrics:
- `colless` - Classic Colless index
- `colless_corr` - Corrected Colless
- `colless_quad` - Quadratic Colless
- `ew_colless` - Equal weights Colless
- `sackin` - Sackin index
- `b1`, `b2` - Beta statistics
- `blum` - Blum index
- `rogers` - Rogers J index
- `root_imbalance` - Root imbalance

### Topology Metrics:
- `cherries` - Number of cherry pairs
- `double_cherries` - Number of double cherries
- `four_prong` - Four-prong configurations
- `pitchforks` - Pitchfork configurations
- `stairs`, `stairs2` - Staircase configurations
- `avg_ladder`, `max_ladder` - Ladder metrics

### Tree Shape:
- `diameter` - Tree diameter
- `max_depth` - Maximum depth
- `max_width` - Maximum width
- `average_leaf_depth` - Average leaf depth
- `avg_vert_depth` - Average vertical depth
- `max_del_width` - Maximum delta width

### Distance Metrics:
- `tot_coph` - Total cophenetic distance
- `tot_internal_path` - Total internal path length
- `area_per_pair` - Area per pair
- `wiener` - Wiener index

### Other:
- `i_stat` - I statistic
- `j_one` - J1 statistic
- `rquartet` - Quartet index
- `psv` - Phylogenetic species variability
- `max_betweenness` - Maximum betweenness
- `max_closeness` - Maximum closeness
- `mw_over_md` - Max width over max depth ratio

## Moderate Statistics (3 total)

These take 0.1-1.0 seconds for 10,000 tips - **excluded from fast set**:
- `crown_age` - 0.63 sec
- `tree_height` - 0.63 sec
- `pigot_rho` - 0.88 sec

## Slow Statistics (2 total)

These take > 1 second for 10,000 tips - **must avoid**:
- `mntd` (Mean Nearest Taxon Distance) - 2.86 sec
- `eigen_centrality` - 10.6 sec (skipped on large trees)

## Failed/Unavailable Statistics (22 total)

These either failed due to errors or have missing dependencies:
- `beta`, `gamma` - May require specific tree properties
- `il_number`, `j_stat`, `imbalance_steps` - Implementation issues
- `laplace_spectrum_*` (4 variants) - Spectral methods
- `max_adj`, `min_adj`, `max_laplace`, `min_laplace` - Graph metrics
- `eigen_centralityW`, `max_closenessW` - Weighted versions
- `mpd`, `phylogenetic_div`, `vpd` - Phylogenetic diversity metrics
- `nltt_base` - NLTT base
- `symmetry_nodes` - Symmetry metric
- `tot_path` - Total path length
- `var_depth` - Depth variance

## Files in This Directory

### Core Files:
- **`benchmark_treestats.R`** - Main benchmarking script with timeout protection
- **`calc_fast_stats.R`** - Optimized function to compute only fast statistics
- **`benchmark_results.csv`** - Raw timing data for all statistics
- **`fast_statistics.txt`** - List of 43 fast statistic names
- **`fast_statistics.RData`** - R data object with fast statistic names
- **`benchmark_output.log`** - Full console output from benchmarking run
- **`README.md`** - This file

## Usage

### Using the Optimized Function

```r
# Load the optimized function
source("treestats_benchmark/calc_fast_stats.R")

# Load your tree
library(ape)
tree <- read.nexus("my_tree.nex")

# Compute only fast statistics
stats <- calc_fast_stats(tree)

# Convert to data frame
stats_df <- as.data.frame(stats)

# Save results
write.csv(stats_df, "tree_stats_fast.csv", row.names = FALSE)
```

### Integrating into Existing Pipeline

Replace this:
```r
stats_list <- treestats::calc_all_stats(tree)
```

With this:
```r
source("treestats_benchmark/calc_fast_stats.R")
stats_list <- calc_fast_stats(tree)
```

### Re-running the Benchmark

```bash
cd treestats_benchmark
Rscript benchmark_treestats.R
```

The script will:
1. Test each statistic on 1,000, 5,000, and 10,000 tip trees
2. Apply 30-second timeout to prevent hanging
3. Skip slow statistics on larger trees
4. Generate updated results and recommendations

## Recommendations

### For Large-Scale Analyses (10,000+ trees):
✅ **Use `calc_fast_stats()`** - Computes 43 informative statistics in ~0.02 sec per tree

❌ **Avoid these slow statistics:**
- `eigen_centrality` - 10+ seconds
- `mntd` - ~3 seconds
- `crown_age`, `tree_height`, `pigot_rho` - 0.6-0.9 seconds

### For Small Datasets (< 1,000 trees):
- You can use `treestats::calc_all_stats()` if time permits
- Still avoid `eigen_centrality` unless absolutely necessary

### For Maximum Information:
- The 43 fast statistics cover:
  - Balance/imbalance (10 metrics)
  - Topology (10 metrics)
  - Branch lengths (6 metrics)
  - Tree shape (7 metrics)
  - Distance metrics (4 metrics)
  - Other useful metrics (6 metrics)

## Integration with Main Workflow

This benchmark workflow is **standalone** but informs the main workflow:

1. **Benchmark identifies** fast statistics (43 total)
2. **Curates** size-agnostic subset (14 statistics)
3. **Main workflow uses** the curated list in `workflow/scripts/phase2_processing/05_extract_tree_stats_batch.R`
4. **If treestats updates**, re-run benchmark and potentially update main workflow

See `resources/SIZE_AGNOSTIC_STATS.md` for how the 14 size-agnostic statistics were selected from the 43 fast ones.

## Maintenance

**When to re-run**:
- treestats package updates
- Need to support larger trees
- Want to test new statistics

**Update frequency**: Quarterly or after major treestats releases

**Dependencies**:
- R >= 4.0
- R packages: ape, treestats, yaml
- Snakemake >= 7.0
- Uses same conda environment as main workflow (`workflow/envs/r-treestats.yaml`)

## Historical Context

- **Benchmarking date**: July 2026
- **treestats version**: Latest from CRAN
- **Platform**: macOS (results similar on Linux/Windows)
- **Pre-computed results**: See `resources/benchmark_results.csv`
