"""
Common configuration, global variables, and wildcard constraints.
"""

import sys
import os

# Add workflow scripts to path
sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from snakemake_utils import get_batches


# ============================================================================
# Configuration Loading
# ============================================================================

configfile: "config/simulation_params.yaml"
configfile: "config/system_params.yaml"
configfile: "config/nre_config.yaml"
configfile: "config/evaluation_params.yaml"
configfile: "config/bayesian_params.yaml"  # Bayesian baseline MCMC parameters
configfile: "config/models_config.yaml"


# ============================================================================
# Workflow Seed Management
# ============================================================================

import numpy as np

def get_workflow_seed():
    """
    Get or generate the master workflow seed for reproducibility.
    
    This seed controls all random operations in the workflow (except LHS sampling,
    which has separate frozen seeds). All other seeds are derived from this master seed.
    
    Returns:
        int: The master workflow seed
    """
    workflow_seed = config.get("workflow_seed", None)
    
    if workflow_seed is None:
        # Generate random seed
        workflow_seed = np.random.randint(0, 2**31 - 1)
        print("\n" + "="*70)
        print("🎲 WORKFLOW SEED: {} (randomly generated)".format(workflow_seed))
        print("   Use this seed to reproduce these exact results:")
        print("   Set 'workflow_seed: {}' in config/nre_config.yaml".format(workflow_seed))
        print("="*70 + "\n")
    else:
        print("\n" + "="*70)
        print("🎲 WORKFLOW SEED: {} (from config)".format(workflow_seed))
        print("   Results are reproducible with this seed")
        print("="*70 + "\n")
    
    return workflow_seed

# Generate master seed
WORKFLOW_SEED = get_workflow_seed()

# Derive all other seeds deterministically from master seed
# These ensure different random streams for different purposes
TRAINING_SEED = WORKFLOW_SEED
VALIDATION_SEED = WORKFLOW_SEED + 1000
NEGATIVE_SAMPLING_SEED = WORKFLOW_SEED + 2000
EVALUATION_SEED = WORKFLOW_SEED + 3000

print("Derived seeds:")
print(f"  Training seed: {TRAINING_SEED}")
print(f"  Validation seed: {VALIDATION_SEED}")
print(f"  Negative sampling seed: {NEGATIVE_SAMPLING_SEED}")
print(f"  Evaluation seed: {EVALUATION_SEED}")
print()


# ============================================================================
# Model Configuration Helpers
# ============================================================================

def get_enabled_base_models():
    """Get list of enabled base model IDs (without size suffix)."""
    return [model_id for model_id, cfg in config["models"].items() 
            if cfg.get("enabled", False)]

def get_model_sizes():
    """Get list of available model sizes from nre_config."""
    return list(config["model_sizes"].keys())

def expand_models_with_sizes():
    """Expand base models into full model IDs with size variants.
    
    Example: 'summary_stats' -> ['summary_stats_small', 'summary_stats_medium', 'summary_stats_large']
    """
    base_models = get_enabled_base_models()
    sizes = get_model_sizes()
    
    expanded = []
    for base_model in base_models:
        for size in sizes:
            expanded.append(f"{base_model}_{size}")
    
    return expanded

def get_enabled_uncertainty_methods():
    """Get list of enabled uncertainty quantification methods.
    
    Returns:
        list: Enabled methods (e.g., ['hpd'], ['lf2i'], or ['hpd', 'lf2i'])
    """
    methods = config.get("uncertainty_methods", ["hpd"])
    # Ensure it's a list
    if isinstance(methods, str):
        methods = [methods]
    # Validate methods
    valid_methods = {"hpd", "lf2i"}
    for method in methods:
        if method not in valid_methods:
            raise ValueError(f"Invalid uncertainty method '{method}'. Must be one of: {valid_methods}")
    return methods

def get_base_model_config(model_id):
    """Get configuration for the base model (strips size suffix if present).
    
    Args:
        model_id: Full model ID like 'summary_stats_small' or base model 'summary_stats'
    
    Returns:
        Base model configuration from models_config.yaml
    """
    # Check if model_id has a size suffix
    sizes = get_model_sizes()
    for size in sizes:
        if model_id.endswith(f"_{size}"):
            base_model = model_id[:-len(f"_{size}")]
            return config["models"][base_model]
    
    # No size suffix, return as-is
    return config["models"][model_id]

def get_model_size(model_id):
    """Extract size from model_id (e.g., 'summary_stats_small' -> 'small').
    
    Returns None if no size suffix found.
    """
    sizes = get_model_sizes()
    for size in sizes:
        if model_id.endswith(f"_{size}"):
            return size
    return None

def get_base_model_from_full_id(model_id):
    """Extract base model name from full model_id by removing size suffix.
    
    Args:
        model_id: Full model ID like 'summary_stats_small' or 'treefeatures_medium'
    
    Returns:
        Base model name like 'summary_stats' or 'treefeatures'
        Returns model_id unchanged if no size suffix found.
    """
    sizes = get_model_sizes()
    for size in sizes:
        if model_id.endswith(f"_{size}"):
            return model_id[:-len(f"_{size}")]
    return model_id

# List of enabled base models (without size suffix)
ENABLED_BASE_MODELS = get_enabled_base_models()

# List of all model variants (base_model × size)
ENABLED_MODELS = expand_models_with_sizes()

# List of available sizes
MODEL_SIZES = get_model_sizes()


# ============================================================================
# Simulation Counts
# ============================================================================

N_TRAIN_SIMS = config["n_train_params"] * config["n_replicates"]
N_TEST_SIMS = config["n_test_params"] * config["n_replicates"]

# ============================================================================
# Simulation IDs
# ============================================================================
# Format: {dataset}_p{param_id:05d}_r{rep}
# Replicates share the same parameters/XML

TRAIN_SIM_IDS = [
    f"train_p{i:05d}_r{r}"
    for i in range(config["n_train_params"])
    for r in range(config["n_replicates"])
]

TEST_SIM_IDS = [
    f"test_p{i:05d}_r{r}"
    for i in range(config["n_test_params"])
    for r in range(config["n_replicates"])
]

ALL_SIM_IDS = TRAIN_SIM_IDS + TEST_SIM_IDS

# ============================================================================
# XML Generation (One per parameter set, shared across replicates)
# ============================================================================

TRAIN_PARAM_IDS = list(range(config["n_train_params"]))
TEST_PARAM_IDS = list(range(config["n_test_params"]))

TRAIN_XML_IDS = [f"train_p{i:05d}" for i in TRAIN_PARAM_IDS]
TEST_XML_IDS = [f"test_p{i:05d}" for i in TEST_PARAM_IDS]
ALL_XML_IDS = TRAIN_XML_IDS + TEST_XML_IDS

# ============================================================================
# Batch Configuration
# ============================================================================

XML_BATCH_SIZE = 1000
BEAST_BATCH_SIZE = 100
PHASE2_BATCH_SIZE = 500
CASE_COUNTS_BATCH_SIZE = 1000
TREEFEATURES_BATCH_SIZE = 1000
BASELINE_BATCH_SIZE = 50  # Reduced from 500 for MCMC stability (50 sims × 3 chains = 150 processes)
EPIESTIM_BATCH_SIZE = 500

# ============================================================================
# XML Batching
# ============================================================================

TRAIN_XML_BATCHES = get_batches(TRAIN_XML_IDS, XML_BATCH_SIZE)
n_train_xml_batches = len(TRAIN_XML_BATCHES)

TEST_XML_BATCHES_RAW = get_batches(TEST_XML_IDS, XML_BATCH_SIZE)
# Offset test batch IDs to be globally unique
TEST_XML_BATCHES = [
    (batch_id + n_train_xml_batches, batch_xmls)
    for batch_id, batch_xmls in TEST_XML_BATCHES_RAW
]

ALL_XML_BATCHES = TRAIN_XML_BATCHES + TEST_XML_BATCHES
N_XML_BATCHES = len(ALL_XML_BATCHES)

# ============================================================================
# BEAST Batching
# ============================================================================

TRAIN_BEAST_BATCHES = get_batches(TRAIN_SIM_IDS, BEAST_BATCH_SIZE)
n_train_beast_batches = len(TRAIN_BEAST_BATCHES)

TEST_BEAST_BATCHES_RAW = get_batches(TEST_SIM_IDS, BEAST_BATCH_SIZE)
# Offset test batch IDs
TEST_BEAST_BATCHES = [
    (batch_id + n_train_beast_batches, batch_sims)
    for batch_id, batch_sims in TEST_BEAST_BATCHES_RAW
]

ALL_BEAST_BATCHES = TRAIN_BEAST_BATCHES + TEST_BEAST_BATCHES
N_BEAST_BATCHES = len(ALL_BEAST_BATCHES)

# ============================================================================
# Phase2 Processing Batching
# ============================================================================

TRAIN_PHASE2_BATCHES = get_batches(TRAIN_SIM_IDS, PHASE2_BATCH_SIZE)
n_train_phase2_batches = len(TRAIN_PHASE2_BATCHES)

TEST_PHASE2_BATCHES_RAW = get_batches(TEST_SIM_IDS, PHASE2_BATCH_SIZE)
# Offset test batch IDs to be globally unique
TEST_PHASE2_BATCHES = [
    (batch_id + n_train_phase2_batches, batch_sims)
    for batch_id, batch_sims in TEST_PHASE2_BATCHES_RAW
]

ALL_PHASE2_BATCHES = TRAIN_PHASE2_BATCHES + TEST_PHASE2_BATCHES
N_PHASE2_BATCHES = len(ALL_PHASE2_BATCHES)

# ============================================================================
# Parameter Scan Evaluation Batching
# ============================================================================

PARAM_SCAN_BATCH_SIZE = config.get("param_scan_batch_size", 500)

PARAM_SCAN_TEST_IDS = TEST_SIM_IDS  # Reuse existing test IDs
PARAM_SCAN_BATCHES = get_batches(PARAM_SCAN_TEST_IDS, PARAM_SCAN_BATCH_SIZE)
N_PARAM_SCAN_BATCHES = len(PARAM_SCAN_BATCHES)

# ============================================================================
# Baseline Method Batching (TEST simulations only)
# ============================================================================

BASELINE_TEST_IDS = TEST_SIM_IDS
BASELINE_BATCHES = get_batches(BASELINE_TEST_IDS, BASELINE_BATCH_SIZE)
N_BASELINE_BATCHES = len(BASELINE_BATCHES)

# ============================================================================
# Bayesian Baseline Batching (Stratified subset of TEST simulations)
# ============================================================================

def get_bayesian_test_ids(n_param_combos, include_all_reps=True):
    """
    Select representative test simulations spanning the parameter space.
    Uses stratified sampling on (R0, recovery_time) to ensure coverage.
    
    Args:
        n_param_combos: Number of parameter combinations to select (not total sims)
        include_all_reps: If True, include all replicates; if False, only first replicate
    
    Returns:
        List of sim_ids for Bayesian baseline
    """
    if n_param_combos is None or n_param_combos <= 0 or n_param_combos < 0:
        # Run all test simulations
        return TEST_SIM_IDS
    
    import pandas as pd
    import numpy as np
    
    # Load test parameters
    params_file = "results/phase1_simulation/parameters/test_design.csv"
    if not os.path.exists(params_file):
        # Fallback: use first n_param_combos if file doesn't exist yet
        print(f"⚠️  {params_file} not found, using first {n_param_combos} parameter combos")
        n_reps = config['n_replicates']
        if include_all_reps:
            return TEST_SIM_IDS[:n_param_combos * n_reps]
        else:
            return [f"test_p{i:05d}_r0" for i in range(n_param_combos)]
    
    # Read test parameters (all rows in test_design.csv are test params)
    params_df = pd.read_csv(params_file)
    
    # Get unique parameter combinations (before replication)
    unique_params = params_df.drop_duplicates(subset=['param_id']).copy()
    
    # If requesting more than available, use all
    if n_param_combos >= len(unique_params):
        return TEST_SIM_IDS
    
    # Stratified sampling: Create 2D bins for (R0, recovery_time)
    n_bins = max(2, int(np.sqrt(n_param_combos)))  # At least 2x2 grid
    
    unique_params['R0_bin'] = pd.qcut(unique_params['R0'], q=n_bins, 
                                       labels=False, duplicates='drop')
    unique_params['rt_bin'] = pd.qcut(unique_params['recovery_time'], q=n_bins,
                                        labels=False, duplicates='drop')
    
    # Sample one parameter set from each bin
    selected_params = (unique_params
                       .groupby(['R0_bin', 'rt_bin'], dropna=False)
                       .sample(n=1, random_state=WORKFLOW_SEED)
                       .head(n_param_combos))
    
    selected_param_ids = selected_params['param_id'].tolist()
    
    # Get sim_ids for selected parameters
    if include_all_reps:
        # All replicates for selected parameters
        bayes_sim_ids = [
            sim_id for sim_id in TEST_SIM_IDS
            if any(f"_p{pid:05d}_" in sim_id for pid in selected_param_ids)
        ]
    else:
        # First replicate only
        bayes_sim_ids = [f"test_p{pid:05d}_r0" for pid in selected_param_ids]
    
    n_reps = config['n_replicates'] if include_all_reps else 1
    print(f"\n📊 Bayesian baseline: Selected {len(selected_param_ids)} parameter combinations")
    print(f"   × {n_reps} replicate(s) = {len(bayes_sim_ids)} total MCMC runs")
    print(f"   R0 range: [{selected_params['R0'].min():.2f}, {selected_params['R0'].max():.2f}]")
    print(f"   Recovery time range: [{selected_params['recovery_time'].min():.2f}, "
          f"{selected_params['recovery_time'].max():.2f}]\n")
    
    return bayes_sim_ids

# Apply Bayesian baseline selection
bayes_config = config.get("bayesian_baseline", {})
if bayes_config.get("enabled", False):
    n_param_combos = bayes_config.get("n_param_combos", 4)
    include_all_reps = bayes_config.get("include_all_replicates", True)
    bayes_batch_size = bayes_config.get("batch_size", 100)
    
    BAYES_TEST_IDS = get_bayesian_test_ids(n_param_combos, include_all_reps)
    BAYES_BATCHES = get_batches(BAYES_TEST_IDS, bayes_batch_size)
    N_BAYES_BATCHES = len(BAYES_BATCHES)
else:
    BAYES_TEST_IDS = []
    BAYES_BATCHES = []
    N_BAYES_BATCHES = 0
    print("📊 Bayesian baseline: Disabled\n")

# ============================================================================
# Wildcard Constraints
# ============================================================================

wildcard_constraints:
    sim_id="(train|test)_p\\d{5}_r\\d+",
    xml_id="(train|test)_p\\d{5}",
    dataset="train|test",
    batch_id="\\d+",
    base_model="[a-z_]+",  # e.g., treefeatures, tsfeatures, combined_stats
    model_id="[a-z_]+",  # e.g., summary_stats_small, case_counts_medium, treefeatures_large
    model_name="[a-zA-Z0-9_]+",
    feat_choice="[a-zA-Z0-9_]+",
