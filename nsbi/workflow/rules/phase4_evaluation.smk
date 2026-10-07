"""
Phase 4: Evaluate trained NRE models

Evaluation pipeline with three components:
1. Parameter scan: Evaluate model on parameter grid to generate likelihood surface
2. Model comparison: Compare performance across models using parameter scan results
3. Uncertainty quantification: LF2I and HPD confidence/credible regions from likelihood surface

Uses config/models_config.yaml to determine which models to evaluate.
"""

# ============================================================================
# Parameter Scan: Likelihood Surface Evaluation
# ============================================================================

rule generate_parameter_scan_grid:
    """Generate grid of parameter values for likelihood surface evaluation."""
    input:
        config_file="config/nre_config.yaml",
        sim_params="config/simulation_params.yaml"
    output:
        "results/phase4_evaluation/parameter_grid/parameter_scan_grid.npz"
    params:
        n_r0_points=config.get("param_scan_r0_points", 30),
        n_rt_points=config.get("param_scan_recovery_time_points", 30)
    log:
        "logs/generate_parameter_scan_grid.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase4_evaluation/02_generate_parameter_scan_grid.py"

rule parameter_scan_batch:
    """Evaluate model on parameter grid for a batch of test simulations."""
    input:
        test_data=lambda w: f"results/phase3_training/training_data/{get_base_model_from_full_id(w.model_id)}_test.npz",
        model="results/phase3_training/trained_models/checkpoints/{model_id}_best.pt",
        config_file="config/nre_config.yaml",
        sim_params="config/simulation_params.yaml",
        grid_file="results/phase4_evaluation/parameter_grid/parameter_scan_grid.npz"
    output:
        temp("results/phase4_evaluation/parameter_scans/{model_id}_param_scan_batch_{batch_id}.csv")
    params:
        sim_ids=lambda wildcards: PARAM_SCAN_BATCHES[int(wildcards.batch_id)][1],
        n_r0_points=config.get("param_scan_r0_points", 30),
        n_rt_points=config.get("param_scan_recovery_time_points", 30),
        model_size=lambda w: get_model_size(w.model_id)
    log:
        "logs/parameter_scan_{model_id}/batch_{batch_id}.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase4_evaluation/03_parameter_scan_evaluation_batch.py"

rule aggregate_parameter_scan:
    """Aggregate all parameter scan batches for a model."""
    input:
        batch_files=lambda w: expand(
            "results/phase4_evaluation/parameter_scans/{model_id}_param_scan_batch_{batch_id}.csv",
            model_id=w.model_id,
            batch_id=range(N_PARAM_SCAN_BATCHES)
        ),
        grid_file="results/phase4_evaluation/parameter_grid/parameter_scan_grid.npz",
        test_data=lambda w: f"results/phase3_training/training_data/{get_base_model_from_full_id(w.model_id)}_test.npz"
    output:
        "results/phase4_evaluation/parameter_scans/{model_id}_param_scan.csv"
    params:
        n_r0_points=config.get("param_scan_r0_points", 30),
        n_rt_points=config.get("param_scan_recovery_time_points", 30)
    log:
        "logs/aggregate_parameter_scan_{model_id}.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase4_evaluation/04_aggregate_parameter_scan.py"

# ============================================================================
# Model Comparison
# ============================================================================

rule compare_models:
    """Compare all enabled models using HPD inference results (point estimates, coverage, scores)."""
    input:
        hpd_regions=lambda w: [f"results/phase4_evaluation/inference/{model_id}/hpd_confidence_regions.csv" 
                              for model_id in ENABLED_MODELS]
    output:
        comparison_yaml="results/phase4_evaluation/model_comparison.yaml",
        comparison_csv="results/phase4_evaluation/model_comparison.csv"
    params:
        model_ids=ENABLED_MODELS
    log:
        "logs/compare_models.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase4_evaluation/05_compare_models.py"


# ============================================================================
# Uncertainty Quantification: LF2I Confidence Regions
# ============================================================================

rule lf2i_calibrate:
    """
    Calibrate LF2I confidence sets using existing parameter scan results.
    
    LF2I (Likelihood-Free Inference with neural ratio estimation):
    - Frequentist-valid confidence regions
    - Requires calibration step using test data
    
    Reference:
    Lueckmann et al. (2021). "Likelihood-free inference with neural ratio estimation"
    Electronic Journal of Statistics. Section 3.3: Test statistic calibration.
    """
    input:
        param_scan="results/phase4_evaluation/parameter_scans/{model_id}_param_scan.csv",
        test_data=lambda w: f"results/phase3_training/training_data/{get_base_model_from_full_id(w.model_id)}_test.npz"
    output:
        calibration="results/phase4_evaluation/inference/{model_id}/lf2i_calibration.pkl",
        diagnostics="results/phase4_evaluation/inference/{model_id}/lf2i_diagnostics.yaml"
    params:
        alpha=0.05,  # 95% confidence intervals
        k_folds=5
    log:
        "logs/calibration/lf2i_calibrate_{model_id}.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase4_evaluation/06_lf2i_calibrate.py"

rule lf2i_compute_regions_batch:
    """Compute LF2I confidence regions for a batch of test observations."""
    input:
        calibration="results/phase4_evaluation/inference/{model_id}/lf2i_calibration.pkl",
        param_scan="results/phase4_evaluation/parameter_scans/{model_id}_param_scan.csv",
        test_data=lambda w: f"results/phase3_training/training_data/{get_base_model_from_full_id(w.model_id)}_test.npz"
    output:
        temp("results/phase4_evaluation/inference/{model_id}/batches/lf2i_batch_{batch_id}.csv")
    params:
        sim_ids=lambda wildcards: PARAM_SCAN_BATCHES[int(wildcards.batch_id)][1]
    log:
        "logs/calibration/lf2i_regions_{model_id}/batch_{batch_id}.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase4_evaluation/07_lf2i_compute_regions_batch.py"

rule lf2i_aggregate_regions:
    """Aggregate all LF2I confidence region batches and generate summary."""
    input:
        batches=lambda w: expand(
            "results/phase4_evaluation/inference/{{model_id}}/batches/lf2i_batch_{batch_id}.csv",
            batch_id=range(N_PARAM_SCAN_BATCHES)
        ),
        calibration="results/phase4_evaluation/inference/{model_id}/lf2i_calibration.pkl",
        test_data=lambda w: f"results/phase3_training/training_data/{get_base_model_from_full_id(w.model_id)}_test.npz"
    output:
        regions="results/phase4_evaluation/inference/{model_id}/lf2i_confidence_regions.csv",
        summary="results/phase4_evaluation/inference/{model_id}/lf2i_summary.txt"
    log:
        "logs/calibration/lf2i_aggregate_{model_id}.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase4_evaluation/08_lf2i_aggregate_regions.py"

# ============================================================================
# Uncertainty Quantification: HPD Credible Regions
# ============================================================================

rule hpd_compute_regions_batch:
    """
    Compute HPD credible regions for a batch of test observations.
    
    HPD (Highest Posterior Density) method:
    - Bayesian credible regions from neural ratio estimator
    - No calibration step required (unlike LF2I)
    - Assumes flat prior
    - Simpler and more interpretable than LF2I
    """
    input:
        param_scan="results/phase4_evaluation/parameter_scans/{model_id}_param_scan.csv",
        test_data=lambda w: f"results/phase3_training/training_data/{get_base_model_from_full_id(w.model_id)}_test.npz"
    output:
        temp("results/phase4_evaluation/inference/{model_id}/batches/hpd_batch_{batch_id}.csv")
    params:
        sim_ids=lambda wildcards: PARAM_SCAN_BATCHES[int(wildcards.batch_id)][1]
    resources:
        mem_mb=2000,  # 2GB should be enough with chunked reading
        runtime=60    # 60 minutes max per batch
    log:
        "logs/calibration/hpd_regions_{model_id}/batch_{batch_id}.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase4_evaluation/09_hpd_compute_regions_batch.py"

rule hpd_aggregate_regions:
    """Aggregate all HPD credible region batches and generate summary."""
    input:
        batches=lambda w: expand(
            "results/phase4_evaluation/inference/{{model_id}}/batches/hpd_batch_{batch_id}.csv",
            batch_id=range(N_PARAM_SCAN_BATCHES)
        ),
        test_data=lambda w: f"results/phase3_training/training_data/{get_base_model_from_full_id(w.model_id)}_test.npz"
    output:
        regions="results/phase4_evaluation/inference/{model_id}/hpd_confidence_regions.csv",
        summary="results/phase4_evaluation/inference/{model_id}/hpd_summary.txt"
    log:
        "logs/calibration/hpd_aggregate_{model_id}.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase4_evaluation/10_hpd_aggregate_regions.py"

# ============================================================================
# Aggregate Rules
# ============================================================================

rule all_lf2i:
    """Generate LF2I confidence regions for all enabled models."""
    input:
        expand(
            "results/phase4_evaluation/inference/{model_id}/lf2i_confidence_regions.csv",
            model_id=ENABLED_MODELS
        )

rule all_hpd:
    """Generate HPD credible regions for all enabled models."""
    input:
        expand(
            "results/phase4_evaluation/inference/{model_id}/hpd_confidence_regions.csv",
            model_id=ENABLED_MODELS
        )

rule all_uncertainty_quantification:
    """
    Generate uncertainty quantification results based on config settings.
    
    Default: HPD credible regions only (faster, no calibration step)
    To enable LF2I: set uncertainty_methods: ['hpd', 'lf2i'] in config/evaluation_params.yaml
    To run only LF2I: set uncertainty_methods: ['lf2i']
    """
    input:
        lambda w: (
            expand("results/phase4_evaluation/inference/{model_id}/lf2i_confidence_regions.csv", 
                   model_id=ENABLED_MODELS) if 'lf2i' in get_enabled_uncertainty_methods() else []
        ) + (
            expand("results/phase4_evaluation/inference/{model_id}/hpd_confidence_regions.csv", 
                   model_id=ENABLED_MODELS) if 'hpd' in get_enabled_uncertainty_methods() else []
        )

