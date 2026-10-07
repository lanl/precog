"""
Phase 2: Extract summary statistics from simulations (BATCHED)
"""

# ============================================================================
# Binned Case Counts Extraction (for tsfeatures)
# ============================================================================

rule extract_binned_case_counts_batch:
    """Extract binned (daily) case counts for tsfeatures calculation."""
    input:
        params_train="results/phase1_simulation/parameters/train_design.csv",
        params_test="results/phase1_simulation/parameters/test_design.csv",
        simulations=expand(
            "results/batch_markers/beast_batch_{batch_id}.done",
            batch_id=range(N_BEAST_BATCHES)
        )
    output:
        touch("results/batch_markers/binned_case_counts_batch_{batch_id}.done")
    params:
        sim_ids=lambda wildcards: ALL_PHASE2_BATCHES[int(wildcards.batch_id)][1],
        bin_width=1.0  # 1 day bins
    log:
        "logs/extract_binned_case_counts/batch_{batch_id}.log"
    threads: 8
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase2_processing/04b_extract_binned_case_counts_batch.py"

# ============================================================================
# Tree Statistics Setup and Extraction
# ============================================================================

rule setup_r_treestats:
    """
    One-time setup: Install treestats and remaining CRAN dependencies.
    
    This rule runs ONCE before any parallel batch jobs to avoid race conditions.
    It installs 5 packages from CRAN that aren't available in conda:
      - pso, subplex, treebalance, DDD, treestats
    
    All other dependencies (31 packages) are pre-installed from conda.
    """
    output:
        marker=touch("results/batch_markers/r_treestats_setup.done")
    log:
        "logs/setup_r_treestats.log"
    conda:
        "../envs/r-treestats.yaml"
    script:
        "../scripts/phase2_processing/00_setup_r_treestats.R"

rule extract_treefeatures_batch:
    input:
        params_train="results/phase1_simulation/parameters/train_design.csv",
        params_test="results/phase1_simulation/parameters/test_design.csv",
        simulations=expand(
            "results/batch_markers/beast_batch_{batch_id}.done",
            batch_id=range(N_BEAST_BATCHES)
        ),
        setup_marker="results/batch_markers/r_treestats_setup.done"
    output:
        touch("results/batch_markers/treefeatures_batch_{batch_id}.done")
    params:
        sim_ids=lambda wildcards: ALL_PHASE2_BATCHES[int(wildcards.batch_id)][1]
    log:
        "logs/extract_treefeatures/batch_{batch_id}.log"
    threads: 8
    conda:
        "../envs/r-treestats.yaml"
    script:
        "../scripts/phase2_processing/05_extract_treefeatures_batch.R"

rule aggregate_train_treefeatures:
    input:
        expand(
            "results/batch_markers/treefeatures_batch_{batch_id}.done",
            batch_id=[batch_id for batch_id, _ in TRAIN_PHASE2_BATCHES]
        ),
        params="results/phase1_simulation/parameters/train_design.csv"
    output:
        "results/phase2_processing/features/treefeatures/train_treefeatures.csv"
    params:
        dataset="train",
        expected_sims=N_TRAIN_SIMS
    log:
        "logs/aggregate_treefeatures_train.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase2_processing/06_aggregate_treefeatures.py"

rule aggregate_test_treefeatures:
    input:
        expand(
            "results/batch_markers/treefeatures_batch_{batch_id}.done",
            batch_id=[batch_id for batch_id, _ in TEST_PHASE2_BATCHES]
        ),
        params="results/phase1_simulation/parameters/test_design.csv"
    output:
        "results/phase2_processing/features/treefeatures/test_treefeatures.csv"
    params:
        dataset="test",
        expected_sims=N_TEST_SIMS
    log:
        "logs/aggregate_treefeatures_test.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase2_processing/06_aggregate_treefeatures.py"


# ============================================================================
# tsfeatures Extraction and Aggregation
# ============================================================================

rule extract_tsfeatures_batch:
    """Extract time series features from BINNED case count trajectories in batch."""
    input:
        params_train="results/phase1_simulation/parameters/train_design.csv",
        params_test="results/phase1_simulation/parameters/test_design.csv",
        binned_case_counts=expand(
            "results/batch_markers/binned_case_counts_batch_{batch_id}.done",
            batch_id=range(N_PHASE2_BATCHES)
        )
    output:
        touch("results/batch_markers/tsfeatures_batch_{batch_id}.done")
    params:
        sim_ids=lambda wildcards: ALL_PHASE2_BATCHES[int(wildcards.batch_id)][1]
    log:
        "logs/extract_tsfeatures/batch_{batch_id}.log"
    threads: 8
    conda:
        "../envs/r-tsfeatures.yaml"
    script:
        "../scripts/phase2_processing/08_extract_tsfeatures_batch.R"

rule aggregate_train_tsfeatures:
    """Aggregate training tsfeatures into single CSV."""
    input:
        expand(
            "results/batch_markers/tsfeatures_batch_{batch_id}.done",
            batch_id=[batch_id for batch_id, _ in TRAIN_PHASE2_BATCHES]
        ),
        params="results/phase1_simulation/parameters/train_design.csv"
    output:
        "results/phase2_processing/features/tsfeatures/train_tsfeatures.csv"
    params:
        dataset="train",
        expected_sims=N_TRAIN_SIMS
    log:
        "logs/aggregate_tsfeatures_train.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase2_processing/09_aggregate_tsfeatures.py"

rule aggregate_test_tsfeatures:
    """Aggregate test tsfeatures into single CSV."""
    input:
        expand(
            "results/batch_markers/tsfeatures_batch_{batch_id}.done",
            batch_id=[batch_id for batch_id, _ in TEST_PHASE2_BATCHES]
        ),
        params="results/phase1_simulation/parameters/test_design.csv"
    output:
        "results/phase2_processing/features/tsfeatures/test_tsfeatures.csv"
    params:
        dataset="test",
        expected_sims=N_TEST_SIMS
    log:
        "logs/aggregate_tsfeatures_test.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase2_processing/09_aggregate_tsfeatures.py"


# ============================================================================
# Combined TSFeatures + Tree Statistics
# ============================================================================

rule merge_train_tsfeatures_treefeatures:
    """Merge tsfeatures with tree statistics for training data."""
    input:
        tsfeatures="results/phase2_processing/features/tsfeatures/train_tsfeatures.csv",
        treefeatures="results/phase2_processing/features/treefeatures/train_treefeatures.csv"
    output:
        "results/phase2_processing/features/tsfeatures_treefeatures/train_tsfeatures_treefeatures.csv"
    log:
        "logs/merge_tsfeatures_treefeatures_train.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase2_processing/10_merge_tsfeatures_treefeatures.py"

rule merge_test_tsfeatures_treefeatures:
    """Merge tsfeatures with tree statistics for test data."""
    input:
        tsfeatures="results/phase2_processing/features/tsfeatures/test_tsfeatures.csv",
        treefeatures="results/phase2_processing/features/treefeatures/test_treefeatures.csv"
    output:
        "results/phase2_processing/features/tsfeatures_treefeatures/test_tsfeatures_treefeatures.csv"
    log:
        "logs/merge_tsfeatures_treefeatures_test.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase2_processing/10_merge_tsfeatures_treefeatures.py"


# ============================================================================
# Cleanup Individual Feature Files (Post-Aggregation)
# ============================================================================

rule cleanup_treefeatures:
    """Remove individual treefeatures CSV files after aggregation completes."""
    input:
        train="results/phase2_processing/features/treefeatures/train_treefeatures.csv",
        test="results/phase2_processing/features/treefeatures/test_treefeatures.csv"
    output:
        touch("results/batch_markers/cleanup_treefeatures.done")
    log:
        "logs/cleanup_treefeatures.log"
    shell:
        """
        echo "Cleaning up individual treefeatures files..." > {log}
        
        # Remove individual feature files
        rm -f results/phase2_processing/features/temp/treefeatures/*_treefeatures.csv 2>> {log}
        
        # Remove temp directory if empty
        rmdir results/phase2_processing/features/temp/treefeatures 2>> {log} || true
        rmdir results/phase2_processing/features/temp 2>> {log} || true
        
        echo "Cleanup complete" >> {log}
        """

rule cleanup_tsfeatures:
    """Remove individual tsfeatures CSV files after aggregation completes."""
    input:
        train="results/phase2_processing/features/tsfeatures/train_tsfeatures.csv",
        test="results/phase2_processing/features/tsfeatures/test_tsfeatures.csv"
    output:
        touch("results/batch_markers/cleanup_tsfeatures.done")
    log:
        "logs/cleanup_tsfeatures.log"
    shell:
        """
        echo "Cleaning up individual tsfeatures files..." > {log}
        
        # Remove individual feature files
        rm -f results/phase2_processing/features/temp/tsfeatures/*_tsfeatures.csv 2>> {log}
        
        # Remove temp directory if empty
        rmdir results/phase2_processing/features/temp/tsfeatures 2>> {log} || true
        rmdir results/phase2_processing/features/temp 2>> {log} || true
        
        echo "Cleanup complete" >> {log}
        """

rule cleanup_all_features:
    """Cleanup all individual feature files after aggregation."""
    input:
        "results/batch_markers/cleanup_treefeatures.done",
        "results/batch_markers/cleanup_tsfeatures.done"
    output:
        touch("results/batch_markers/cleanup_all_features.done")
    log:
        "logs/cleanup_all_features.log"
    shell:
        """
        echo "All feature cleanup complete" > {log}
        echo "Individual feature files removed" >> {log}
        echo "Aggregated features retained in:" >> {log}
        echo "  - results/phase2_processing/features/treefeatures/" >> {log}
        echo "  - results/phase2_processing/features/tsfeatures/" >> {log}
        """

