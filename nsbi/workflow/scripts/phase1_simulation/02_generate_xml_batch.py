"""
Generate BEAST2 XML files in batch with parallel processing.
NEW: Generates ONE XML per unique parameter set (not per replicate).
Replicates will reuse the same XML - BEAST stochasticity provides variation.
"""

import pandas as pd
import os
from multiprocessing import Pool
from functools import partial

def generate_one_xml(xml_id, params_train_df, params_test_df, initial_infected, min_recovered):
    """
    Generate a single XML file for one parameter set.
    
    Args:
        xml_id: XML ID (e.g., "train_p00042" for param_id=42)
        params_train_df: DataFrame with training parameters
        params_test_df: DataFrame with test parameters
        initial_infected: Number of initially infected individuals
        min_recovered: Minimum recovered required for valid outbreak
    
    Returns:
        Path to generated XML file
    """
    # Determine if this is train or test and extract param_id
    if xml_id.startswith('train_'):
        param_id = int(xml_id.replace('train_p', ''))
        params_df = params_train_df
    elif xml_id.startswith('test_'):
        param_id = int(xml_id.replace('test_p', ''))
        params_df = params_test_df
    else:
        raise ValueError(f"Unknown xml_id format: {xml_id}")
    
    # Get parameters for this param_id
    params = params_df.iloc[param_id]
    
    S = params['S']
    R0 = params['R0']
    recovery_time = params['recovery_time']
    gamma = params['gamma']
    beta = params['beta']
    
    # Calculate total population size (constant in closed SIR model)
    N = S + initial_infected  # S + I + R at time 0 (R=0 initially)
    
    # Calculate frequency-dependent transmission rate
    # This ensures R0 is independent of population size
    beta_freq = beta / N
    
    # Create XML content (generic, no replicate-specific info)
    # Output filenames will be specified at BEAST runtime
    xml_output = f"""<!-- Stochastic SIR trajectory simulation for parameter set {param_id} -->
<!-- Parameters: S={S}, R0={R0:.4f}, recovery_time={recovery_time:.4f}, beta={beta:.6f}, gamma={gamma:.6f} -->
<!-- N={N}, beta/N={beta_freq:.8f} (frequency-dependent transmission) -->
<!-- NOTE: This XML will be reused for multiple replicates. BEAST's stochastic nature provides variation. -->

<beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">

    <run spec="Simulator" nSims="1">
        <simulate spec="SimulatedTree" id="tree">
            <trajectory spec="StochasticTrajectory" id="SIRTrajectory" mustHave="R&gt;{min_recovered}">
                <population spec="RealParameter" id="S" value="{S}"/>
                <population spec="RealParameter" id="I" value="{initial_infected}"/>
                <samplePopulation spec="RealParameter" id="R" value="0"/>

                <reaction spec="Reaction" rate="{beta_freq}"> S + I -> 2I </reaction>
                <reaction spec="Reaction" rate="{gamma}"> I -> R </reaction>
            </trajectory>
        </simulate>

        <logger spec="Logger" fileName="{xml_id}.traj">
            <log idref="SIRTrajectory"/>
        </logger>

        <logger spec="Logger" mode="tree" fileName="{xml_id}.trees">
          <log spec="TypedTreeLogger" tree="@tree" />
        </logger>

    </run>
</beast>
"""
    
    # Write XML file
    output_xml = f"results/phase1_simulation/xml_files/{xml_id}.xml"
    os.makedirs(os.path.dirname(output_xml), exist_ok=True)
    with open(output_xml, 'w') as f:
        f.write(xml_output)
    
    return output_xml


# Main batch processing logic
if __name__ == '__main__':
    # Access snakemake variables
    template_file = snakemake.input.template
    params_train_file = snakemake.input.params_train
    params_test_file = snakemake.input.params_test
    xml_ids = snakemake.params.xml_ids
    initial_infected = snakemake.params.initial_infected
    min_recovered = snakemake.params.min_recovered
    n_threads = snakemake.threads
    
    print(f"Generating {len(xml_ids)} XML files in parallel using {n_threads} threads")
    print(f"XML IDs: {xml_ids[0]} ... {xml_ids[-1]}")
    
    # Load parameter DataFrames
    params_train_df = pd.read_csv(params_train_file)
    params_test_df = pd.read_csv(params_test_file)
    
    print(f"Loaded parameters:")
    print(f"  Train: {len(params_train_df)} unique parameter sets")
    print(f"  Test: {len(params_test_df)} unique parameter sets")
    
    # Create partial function with fixed arguments
    generate_func = partial(
        generate_one_xml,
        params_train_df=params_train_df,
        params_test_df=params_test_df,
        initial_infected=initial_infected,
        min_recovered=min_recovered
    )
    
    # Generate XMLs in parallel
    with Pool(n_threads) as pool:
        results = pool.map(generate_func, xml_ids)
    
    print(f"✓ Generated {len(results)} XML files")
    print(f"  These XMLs will be reused for multiple BEAST runs (replicates)")
