/*
 * Module: PP (Preprocessing)
 * 
 * Implements bin/preprocessing.py for scRNA-seq data preprocessing
 * Supports both single parameter and parameter grid execution
 */

process PP {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}/preprocessing", mode: 'copy'
    
    // Container for reproducibility (optional)
    // container 'scanpy/scanpy:latest'
    
    input:
    // Exercise 1: Single parameter input
    // TODO: Define input tuple for Exercise 1
    
    // Exercise 2: Parameter grid input  
    // TODO: Define input tuple for Exercise 2 with parameter combinations
    
    output:
    // TODO: Define output tuple - processed data and QC plots
    
    script:
    // Exercise 1: Single parameter execution
    // TODO: Construct preprocessing.py command with single parameters
    
    // Exercise 2: Parameter grid execution
    // TODO: Construct preprocessing.py command with parameter grid values
    // TODO: Include parameter combination in output naming
    
    stub:
    // TODO: Create stub outputs for testing
}