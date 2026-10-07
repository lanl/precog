"""
Phase 3: Prepare training data and train NRE models

Uses config/models_config.yaml to determine which models to train.
Data preparation happens once per base model (feature set).
Training happens for each size variant using the shared data.
"""

# ============================================================================
# Generic Rules (work for any model via wildcards)
# ============================================================================

rule prepare_base_model_data:
    """Prepare training data for a base model (feature set).
    
    This rule runs once per base model (e.g., 'treefeatures', 'tsfeatures').
    The prepared data is shared across all size variants (small, medium, large).
    """
    input:
        train_input=lambda w: config["models"][w.base_model]["input_data"]["train"],
        test_input=lambda w: config["models"][w.base_model]["input_data"]["test"]
    output:
        train_data="results/phase3_training/training_data/{base_model}_train.npz",
        test_data="results/phase3_training/training_data/{base_model}_test.npz"
    params:
        script=lambda w: config["models"][w.base_model]["prep_script"]
    log:
        "logs/prepare_{base_model}_data.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/{params.script}"

rule train_model:
    """Train NRE model for a specific model with specific size.
    
    Model ID format: {base_model}_{size} (e.g., summary_stats_small)
    The size determines the hidden_dims used for the architecture.
    All size variants share the same training data (base model).
    """
    input:
        train_data=lambda w: f"results/phase3_training/training_data/{get_base_model_from_full_id(w.model_id)}_train.npz",
        config_file="config/nre_config.yaml"
    output:
        model_final="results/phase3_training/trained_models/checkpoints/{model_id}_final.pt",
        model_best="results/phase3_training/trained_models/checkpoints/{model_id}_best.pt",
        training_log="results/phase3_training/trained_models/training_logs/{model_id}_training.csv"
    params:
        training_seed=TRAINING_SEED,
        validation_seed=VALIDATION_SEED,
        model_size=lambda w: get_model_size(w.model_id)
    log:
        "logs/train_{model_id}.log"
    threads: 8
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase3_training/03_train_nre.py"
