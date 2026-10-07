"""
Phase 1: LHS Parameter Design, XML Generation, and BEAST2 Simulation
"""

# ============================================================================
# LHS Parameter Design
# ============================================================================

rule generate_train_lhs:
    output:
        "results/phase1_simulation/parameters/train_design.csv"
    params:
        dataset="train",
        n_samples=config["n_train_params"],
        n_reps=config["n_replicates"],
        S_min=config["S_min"],
        S_max=config["S_max"],
        S_fixed=config.get("S_fixed", None),
        R0_min=config["R0_min"],
        R0_max=config["R0_max"],
        rt_min=config["recovery_time_min"],
        rt_max=config["recovery_time_max"],
        seed=config["lhs_train_seed"]
    log:
        "logs/lhs_design_train.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase1_simulation/01_generate_lhs_design.py"

rule generate_test_lhs:
    output:
        "results/phase1_simulation/parameters/test_design.csv"
    params:
        dataset="test",
        n_samples=config["n_test_params"],
        n_reps=config["n_replicates"],
        S_min=config["S_min"],
        S_max=config["S_max"],
        S_fixed=config.get("S_fixed", None),
        R0_min=config["R0_min"],
        R0_max=config["R0_max"],
        rt_min=config["recovery_time_min"],
        rt_max=config["recovery_time_max"],
        seed=config["lhs_test_seed"]
    log:
        "logs/lhs_design_test.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase1_simulation/01_generate_lhs_design.py"

# ============================================================================
# XML Generation (Batched)
# ============================================================================

rule generate_xml_batch:
    input:
        template="resources/sim_sir_trees.xml",
        params_train="results/phase1_simulation/parameters/train_design.csv",
        params_test="results/phase1_simulation/parameters/test_design.csv"
    output:
        marker=touch("results/batch_markers/xml_batch_{batch_id}.done")
    params:
        xml_ids=lambda wildcards: [
            xml_id for batch_id, batch_xmls in ALL_XML_BATCHES
            for xml_id in batch_xmls
            if batch_id == int(wildcards.batch_id)
        ],
        initial_infected=config["initial_infected"],
        min_recovered=config["min_recovered"]
    log:
        "logs/generate_xml_batch/batch_{batch_id}.log"
    threads: 8
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase1_simulation/02_generate_xml_batch.py"

# ============================================================================
# BEAST2 Simulation (Batched)
# ============================================================================

rule run_beast_batch:
    input:
        xml_batches=expand(
            "results/batch_markers/xml_batch_{batch_id}.done",
            batch_id=range(N_XML_BATCHES)
        )
    output:
        marker=touch("results/batch_markers/beast_batch_{batch_id}.done")
    params:
        sim_ids=lambda wildcards: ALL_BEAST_BATCHES[int(wildcards.batch_id)][1],
        beast_path=config["beast_path"]
    log:
        "logs/beast_batch/batch_{batch_id}.log"
    threads: 1
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase1_simulation/03_run_beast_batch.py"
