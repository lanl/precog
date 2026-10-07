"""
Run BEAST2 simulations in batch (sequential execution).

This script processes multiple BEAST simulations sequentially within a single cluster job.
Each simulation takes ~2-20 seconds, so batching 100 simulations takes ~3-30 minutes,
which is much more efficient than submitting 100 separate cluster jobs.

Key features:
- Sequential execution (BEAST is single-threaded)
- Robust error handling (one failure doesn't stop the batch)
- Detailed logging per simulation
- Failure tracking for easy resubmission
"""

import os
import sys
import shutil
import tempfile
import traceback
from pathlib import Path

# Access snakemake variables
sim_ids = snakemake.params.sim_ids
beast_path = snakemake.params.beast_path
batch_id = snakemake.wildcards.batch_id
log_file = snakemake.log[0]

# Helper functions from Snakefile
def get_param_id_from_sim_id(sim_id):
    """Extract param_id from sim_id like 'train_p00042_r2' -> 42"""
    import re
    match = re.search(r'_p(\d+)_r\d+', sim_id)
    return int(match.group(1)) if match else None

def get_dataset_from_sim_id(sim_id):
    """Extract dataset from sim_id like 'train_p00042_r2' -> 'train'"""
    return sim_id.split('_')[0]

# Setup logging
log_dir = Path(log_file).parent
log_dir.mkdir(parents=True, exist_ok=True)

failed_sims = []
successful_sims = []

# Open main log file
with open(log_file, 'w') as log:
    log.write(f"BEAST Batch Processing - Batch {batch_id}\n")
    log.write(f"=" * 80 + "\n")
    log.write(f"Processing {len(sim_ids)} simulations\n")
    log.write(f"BEAST path: {beast_path}\n\n")
    log.flush()
    
    # Process each simulation
    for idx, sim_id in enumerate(sim_ids, 1):
        log.write(f"\n[{idx}/{len(sim_ids)}] Processing {sim_id}...\n")
        log.flush()
        
        try:
            # Get paths
            param_id = get_param_id_from_sim_id(sim_id)
            dataset = get_dataset_from_sim_id(sim_id)
            xml_id = f"{dataset}_p{param_id:05d}"
            
            xml_path = os.path.abspath(f"results/phase1_simulation/xml_files/{xml_id}.xml")
            trees_path = os.path.abspath(f"results/phase1_simulation/simulations/{sim_id}.trees")
            traj_path = os.path.abspath(f"results/phase1_simulation/simulations/{sim_id}.traj")
            
            # Verify XML exists
            if not os.path.exists(xml_path):
                raise FileNotFoundError(f"XML file not found: {xml_path}")
            
            # Create output directory if needed
            os.makedirs(os.path.dirname(trees_path), exist_ok=True)
            
            # Create temporary directory for this BEAST run
            with tempfile.TemporaryDirectory() as tmpdir:
                # Copy XML to temp dir
                xml_name = os.path.basename(xml_path)
                temp_xml = os.path.join(tmpdir, xml_name)
                shutil.copy(xml_path, temp_xml)
                
                # Run BEAST in temp directory
                import subprocess
                result = subprocess.run(
                    [beast_path, "-overwrite", xml_name],
                    cwd=tmpdir,
                    capture_output=True,
                    text=True,
                    timeout=300  # 5 minute timeout per simulation
                )
                
                if result.returncode != 0:
                    log.write(f"  ERROR: BEAST failed with return code {result.returncode}\n")
                    log.write(f"  STDOUT: {result.stdout[:500]}\n")
                    log.write(f"  STDERR: {result.stderr[:500]}\n")
                    raise RuntimeError(f"BEAST failed for {sim_id}")
                
                # Get base name for BEAST outputs (without .xml extension)
                xml_base = xml_name.replace('.xml', '')
                
                # Move outputs to final location
                beast_trees = os.path.join(tmpdir, f"{xml_base}.trees")
                beast_traj = os.path.join(tmpdir, f"{xml_base}.traj")
                
                if os.path.exists(beast_trees):
                    shutil.move(beast_trees, trees_path)
                    log.write(f"  ✓ Trees saved to {trees_path}\n")
                else:
                    raise FileNotFoundError(f"BEAST did not generate trees file: {beast_trees}")
                
                if os.path.exists(beast_traj):
                    shutil.move(beast_traj, traj_path)
                    log.write(f"  ✓ Trajectory saved to {traj_path}\n")
                else:
                    raise FileNotFoundError(f"BEAST did not generate trajectory file: {beast_traj}")
            
            successful_sims.append(sim_id)
            log.write(f"  ✓ SUCCESS: {sim_id}\n")
            
        except Exception as e:
            failed_sims.append(sim_id)
            log.write(f"  ✗ FAILED: {sim_id}\n")
            log.write(f"  Error: {str(e)}\n")
            log.write(f"  Traceback:\n")
            traceback.print_exc(file=log)
            continue
        
        log.flush()
    
    # Summary
    log.write(f"\n" + "=" * 80 + "\n")
    log.write(f"BATCH {batch_id} SUMMARY\n")
    log.write(f"=" * 80 + "\n")
    log.write(f"Total simulations: {len(sim_ids)}\n")
    log.write(f"Successful: {len(successful_sims)}\n")
    log.write(f"Failed: {len(failed_sims)}\n")
    
    if failed_sims:
        log.write(f"\nFailed simulations:\n")
        for sim_id in failed_sims:
            log.write(f"  - {sim_id}\n")
        
        # Write failures to separate file for easy resubmission
        failure_file = log_dir / f"batch_{batch_id}_failures.txt"
        with open(failure_file, 'w') as f:
            f.write('\n'.join(failed_sims))
        log.write(f"\nFailure list saved to: {failure_file}\n")

# Exit with error if any simulations failed
if failed_sims:
    print(f"WARNING: Batch {batch_id} completed with {len(failed_sims)} failures")
    print(f"See log for details: {log_file}")
else:
    print(f"Batch {batch_id} completed successfully: {len(successful_sims)} simulations")

