#!/bin/bash

set -euo pipefail

# Clean embeddings once at the very beginning, before any Slurm jobs are submitted
EMB_DIR="data/embeddings_gam_real"
echo "Cleaning ${EMB_DIR} ..."
rm -rf "${EMB_DIR:?}"
mkdir "${EMB_DIR:?}"


rm -f \
data/fail_msg.RData \
data/inputs_and_outputs_sirMappings_byRegime_Crash.RData \
data/inputs_and_outputs_sirMappings_byRegime_Dec.RData \
data/inputs_and_outputs_sirMappings_byRegime_Inc.RData \
data/inputs_and_outputs_sirMappings_byRegime_NearPeak.RData \
data/inputs_and_outputs_sirMappings_byRegime_Surge.RData \
data/inputs_sirOutputFilling_byRegime_Crash.RData \
data/inputs_sirOutputFilling_byRegime_Dec.RData \
data/inputs_sirOutputFilling_byRegime_Inc.RData \
data/inputs_sirOutputFilling_byRegime_NearPeak.RData \
data/inputs_sirOutputFilling_byRegime_Surge.RData \
data/outputs_sirOutputFilling_byRegime_Crash.RData \
data/outputs_sirOutputFilling_byRegime_Dec.RData \
data/outputs_sirOutputFilling_byRegime_Inc.RData \
data/outputs_sirOutputFilling_byRegime_NearPeak.RData \
data/outputs_sirOutputFilling_byRegime_Surge.RData \
data/synthetic_embeddings.RData

rm -f \
viz/sir_by_regime_experiment.png \
viz/smoa_performances.png \
viz/average_distances_to_synthetic.png \
viz/performance_boxplots.png \
viz/performance_rw.png \
viz/performance_boxplots_byregime.png \
viz/pca.png

echo "Clearning synthetic logfiles..."
rm -rf data/stochastic_sir
rm -rf SLURMOUT
rm -rf logfiles

mkdir data/stochastic_sir
mkdir SLURMOUT
mkdir logfiles

echo "Submitting first array job: regime sampling..."
job1=$(sbatch --parsable regime_outputs_run.sh)
echo "Submitted regime sampling array as job ${job1}"

echo "Submitting embedding preparation job, dependent on completion of regime sampling..."
job2=$(sbatch --parsable --dependency=afterok:${job1} prep_embeddings.sh)
echo "Submitted embedding preparation job as ${job2}"

echo "Submitting work manifest creation job, dependent on completion of embedding preparation..."
job3=$(sbatch --parsable --dependency=afterok:${job2} create_manifest.sh)
echo "Submitted work manifest creation job as ${job3}"

echo "Submitting fourth array job: SMOA overall, dependent on completion of manifest creation..."
job4=$(sbatch --parsable --dependency=afterok:${job3} perform_smoa_allData.sh)
echo "Submitted SMOA overall array as job ${job4}"

echo "Submitting final plotting job, dependent on entire SMOA array..."
job5=$(sbatch --parsable --dependency=afterok:${job4} submit_final_plots.sh)
echo "Submitted final plotting job as ${job5}"

echo "Done."
echo "First array job:         ${job1} (regime sampling)"
echo "Second job:              ${job2} (embedding preparation)"
echo "Third job:               ${job3} (work manifest creation)"
echo "Fourth array job:        ${job4} (SMOA overall)"
echo "Final plot job:          ${job5}"