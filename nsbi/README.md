# Neural Simulation-Based Inference for Epidemic Forecasting

A Snakemake-based pipeline for training neural ratio estimation (NRE) models to perform inference on epidemic parameters from phylogenetic trees and case count data.

The full_output directory contains the parameters, simulation outputs, training and testing data, models, and evaluation metrics. The reports directory contains the quarto files to generate the paper tables, figures, and summaries.

## Overview

This pipeline simulates SIR (Susceptible-Infectious-Recovered) epidemic dynamics using BEAST2, extracts summary statistics from phylogenetic trees and time series data, and trains neural networks to estimate posterior distributions of epidemic parameters.

**Key Features:**
- Automated end-to-end workflow from simulation to inference
- Frequency-dependent transmission for population-independent R₀
- Tree size-agnostic statistics (14 inherently size-independent metrics)
- Cross-platform conda environments (macOS/Linux)
- Auto-installation of R packages from CRAN on first run
- Fixed-width time binning for temporal alignment
- Multiple model comparison: tree statistics, time series features (tsfeatures), and combined models
- Support for local and HPC cluster execution

## Quick Start

```bash
# 1. Clone repository
git clone <repo-url>
cd neural_sbi

# 2. Create workflow environment
conda env create -f workflow/envs/snakemake.yaml
conda activate neural_sbi_workflow

# 3. Configure BEAST2 path
nano config/simulation_params.yaml
# Set: beast_path: "/path/to/beast"

# 4. Run pipeline (creates conda envs automatically)
snakemake --use-conda --cores 4

# 5. First run: R packages auto-install from CRAN (~5-10 min compilation)
# 6. Subsequent runs: Fast (packages cached)

# 7. View results
ls -lh results/evaluation/
```

## Project Structure

Following [Snakemake best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html):

```
neural_sbi/
├── README.md                     # This file
├── CHANGELOG.md                  # Complete change history
├── workflow/                     # Snakemake workflow (STANDARD STRUCTURE)
│   ├── Snakefile                # Main workflow entrypoint
│   ├── rules/                   # Modular workflow rules
│   │   ├── common.smk              # Configuration and global variables
│   │   ├── phase1_simulation.smk   # LHS, XML generation, BEAST
│   │   ├── phase2_processing.smk   # Extraction and aggregation
│   │   ├── phase3_training.smk     # Training data and models
│   │   ├── phase4_evaluation.smk   # Evaluation and visualization
│   │   └── phase5_baselines.smk    # Baseline methods
│   ├── envs/                    # Fine-grained conda environments
│   │   ├── python-core.yaml        # Python analysis tools
│   │   ├── r-treestats.yaml        # R tree statistics
│   │   ├── r-epiestim.yaml         # R baseline methods
│   │   ├── r-quarto.yaml           # Notebook rendering
│   │   └── snakemake.yaml          # Workflow management
│   ├── scripts/                 # All workflow scripts (organized by phase)
│   │   ├── phase1_simulation/      # LHS design, XML, BEAST
│   │   ├── phase2_processing/      # Case counts, tree stats, aggregation
│   │   ├── phase3_training/        # Feature selection, training data, NRE
│   │   ├── phase4_evaluation/      # Model evaluation and comparison
│   │   ├── phase5_baselines/       # Baseline methods and comparison
│   │   ├── utils/                  # Shared utility functions
│   │   └── snakemake_utils.py      # Snakemake helper functions
│   └── notebooks/               # Quarto/Jupyter notebooks
│       └── eda.qmd                 # Exploratory data analysis
├── config/                      # Configuration files (KEPT SEPARATE)
│   ├── simulation_params.yaml      # Phase 1-2: Simulations
│   ├── processing_params.yaml      # Phase 3: Processing
│   ├── training_params.yaml        # Phase 4-5: Training
│   ├── evaluation_params.yaml      # Phase 6+: Evaluation
│   ├── system_params.yaml          # System settings
│   ├── cluster_config.yaml         # SLURM resource allocation
│   └── slurm_submit.sh             # SLURM submission script
├── resources/                   # Static input files
│   └── sim_sir_trees.xml           # BEAST2 XML template
├── benchmark/                   # Benchmarking sub-workflows
│   └── treestats/                  # Tree statistics performance benchmarking
│       ├── workflow/
│       │   ├── Snakefile
│       │   └── scripts/
│       ├── config/
│       ├── resources/
│       └── README.md
├── docs/                        # Documentation
│   ├── QUICKSTART.md               # Quick start guide
│   ├── INSTALL_NOTES.md            # Installation instructions
│   └── CLUSTER_GUIDE.md            # Complete HPC cluster guide
├── results/                     # Pipeline outputs (generated, git-ignored)
└── logs/                        # Execution logs (generated, git-ignored)
```

## Pipeline Steps

1. **LHS Parameter Sampling** → Design training/test parameter sets
2. **XML Generation** → Create BEAST2 input files with frequency-dependent transmission
3. **BEAST2 Simulation** → Run epidemic simulations with phylogenetic trees
4. **Case Count Extraction** → Extract I(t) time series from trajectories
5. **Tree Statistics** → Compute 70+ phylogenetic statistics using R
6. **Data Aggregation** → Combine time series and tree stats with fixed-width binning
7. **Training Data Preparation** → Create positive/negative pairs for NRE
8. **Model Training** → Train two models: (θ,y) and (θ,y,z)
9. **Model Evaluation** → Compare posterior inference performance

## Key Configuration Parameters

Configuration is split across multiple files for better organization. See `docs/CONFIG_QUICK_REFERENCE.md` for details.

**`config/simulation_params.yaml`** - Simulation settings:
```yaml
# Number of simulations
n_train_params: 20       # Start with 20, scale to 25,000+ for production
n_test_params: 5
n_replicates: 2

# Parameter ranges (uniform priors)
S_min: 100               # Susceptible population
S_max: 10000
S_fixed: 5000           # Optional: fix S, estimate only R0 & recovery_time
R0_min: 1.0             # Basic reproduction number
R0_max: 10.0
recovery_time_min: 1.0  # Days to recovery
recovery_time_max: 14.0
```

**`config/processing_params.yaml`** - Processing settings:
```yaml
# Time binning
bin_width: 1.0          # Fixed 1-day bins

# Tree statistics
tree_stats:
  size_agnostic: true   # Use size-normalized statistics
  normalization: "both" # Yule and PDA normalizations
  normalize_time_series: true

# Feature selection
feature_selection:
  enabled: false
  redundancy_threshold: 0.95
```

**`config/training_params.yaml`** - Neural network settings:
```yaml
nre_hidden_features: 50
nre_num_transforms: 5
nre_learning_rate: 0.0005
nre_training_batch_size: 512
early_stopping_patience: 50
```

**`config/system_params.yaml`** - Paths and resources:
```yaml
beast_path: "/Applications/BEAST 2.7.7/bin/beast"
threads: 8
```

# Random seeds (for reproducibility)
lhs_train_seed: 42
lhs_test_seed: 123
training_seed: 42
```

## Critical Implementation Details

### Frequency-Dependent Transmission
The pipeline uses **frequency-dependent transmission** to ensure R₀ is independent of population size:

```
beta = R₀ × gamma / N
```

This prevents effective R₀ from scaling with population size (critical for valid inference).

### Fixed-Width Time Binning
Time series are binned using **fixed 1-day bins** (not fixed bin count) with zero-padding to ensure:
- Temporal alignment across simulations
- Consistent vector dimensions
- Interpretable time units

### Comprehensive Tree Statistics (Strict Size-Agnostic by Default)
The pipeline uses R's `treestats` package with **strict size-agnostic selection** (UPDATED):

**Default behavior** (`size_agnostic=true`):
- Computes **14 inherently size-agnostic statistics** (vs 66 previously)
- No Yule/PDA normalization (those methods are not truly tree size-agnostic)
- Only includes statistics that are independent of tree size by their nature
- ~70% faster computation than normalized approach

**Why strict agnostic mode?**
- Trees are generated from SIR coalescent process (backward-time)
- Yule (forward birth) and PDA (uniform topology) assume specific tree models
- Normalization introduces model assumptions that may not hold for epidemic trees
- Cleaner: use only inherently independent statistics

**Inherently size-agnostic statistics (14 total):**
- **Branch length measures (7)**: mean, variance (internal/external/all), treeness
- **Balance ratios (3)**: beta, equal-weights Colless, root imbalance  
- **Ratio metric (1)**: max width / max depth
- **Information-theoretic (3)**: I-statistic, J-one, J-statistic
- **Lineages Through Time (LTT)** - normalized by n_tips (proportion of final diversity)
- **Nodes Through Time (NTT)** - normalized by (n_tips - 1) (probability distribution)

**Time series normalization:**
- **LTT**: Divided by n_tips → values range from 0.0 to 1.0
- **NTT**: Divided by (n_tips - 1) → sums to 1.0
- Makes time series independent of tree size for better generalization

See `treestats_benchmark/SIZE_AGNOSTIC_STATS.md` for detailed explanation.

**To use all 51 original statistics:**
```yaml
tree_stats:
  size_agnostic: false  # Computes 51 statistics (all original, no normalization)
  normalize_time_series: false  # Use raw LTT/NTT counts
```

## Documentation

- **[QUICKSTART.md](docs/QUICKSTART.md)** - Get started in 5 minutes
- **[INSTALL_NOTES.md](docs/INSTALL_NOTES.md)** - Detailed installation guide
- **[CLUSTER_GUIDE.md](docs/CLUSTER_GUIDE.md)** - Complete HPC cluster deployment guide
- **[TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md)** - Common issues and solutions
- **[ARCHITECTURE.md](docs/ARCHITECTURE.md)** - Technical design and implementation details
- **[CHANGELOG.md](CHANGELOG.md)** - Complete change history and version notes
- **[DEV_GUIDELINES.md](DEV_GUIDELINES.md)** - Development guidelines and best practices

### For Contributors
**Please read [DEV_GUIDELINES.md](DEV_GUIDELINES.md)** before contributing. Key rules:
- ❌ **Zero documentation for small fixes** (typos, minor changes)
- ✏️ **Concise documentation for major changes** (1 paragraph in CHANGELOG.md)
- 🗑️ **Delete test files immediately after use**
- 📦 **Archive old documentation immediately**
- 🚫 **Never create** `BUGFIX_*.md`, `IMPLEMENTATION_*.md`, `SUMMARY_*.md`

## Requirements

### Software
- **Python 3.10+** (managed by conda)
- **R 4.0+** (managed by conda)
- **BEAST2** (tested with 2.7.7) - requires separate installation
- **R treestats package** (manual install after conda env setup)

### Hardware
- ~16GB RAM for training (scales with simulation count)
- GPU optional (CPU training works fine for small datasets)

### Installation
All Python and R dependencies are managed by conda through fine-grained environment files:
```bash
# Create workflow environment
conda env create -f workflow/envs/snakemake.yaml
conda activate neural_sbi_workflow

# Install R treestats package (not available via conda)
conda activate neural_sbi_treestats
Rscript -e 'install.packages("treestats", repos="http://cran.r-project.org")'
```

**Note:** The R treestats package cannot be installed via conda-forge and must be installed from CRAN manually.

## Citation

If you use this pipeline, please cite:
- BEAST2: Bouckaert et al. (2019)
- sbi package: Tejero-Cantero et al. (2020)
- treestats: van der Bijl (2023)

## License

MIT License - See LICENSE file for details

## Contact

For questions or issues, please open a GitHub issue or contact the maintainers.

---

**Status**: Pipeline validated and production-ready ✅  
**Last Updated**: August 2026  
**Version**: Production v1.0

