#!/usr/bin/env nextflow

/*
 * Pipeline 1: Basic Nextflow Concepts
 * 
 * This pipeline demonstrates fundamental Nextflow concepts:
 * - Process definitions
 * - Input/output declarations  
 * - Channel operations
 * - Workflow composition
 * - Parameter handling
 */

nextflow.enable.dsl = 2

// Pipeline parameters with default values
params.input = "$baseDir/data/sample.fa"
params.outdir = "results"

/*
 * PROCESS: Split sequences
 * 
 * Takes a FASTA file and splits it into individual sequence files.
 * Each sequence gets its own file named seq_1, seq_2, etc.
 */
process splitSequences {
    // Process directive to specify output directory
    publishDir "${params.outdir}/split", mode: 'copy'

    input:
    path input_file

    output:
    path 'seq_*'

    script:
    """
    echo "Splitting sequences from ${input_file}..."
    awk '/^>/{f="seq_"++d} {print > f}' < ${input_file}
    """
}

/*
 * PROCESS: Reverse sequences
 * 
 * Takes individual sequence files and reverses their content.
 * Demonstrates process chaining and parallel execution.
 */
process reverseSequences {
    // Process directive for output
    publishDir "${params.outdir}/reversed", mode: 'copy'
    
    input:
    path sequence_file

    output:
    path "${sequence_file.baseName}_reversed.txt"

    script:
    """
    echo "Processing ${sequence_file}..."
    cat ${sequence_file} | rev > ${sequence_file.baseName}_reversed.txt
    """
}

/*
 * WORKFLOW: Main workflow definition
 * 
 * Defines the execution flow:
 * 1. Split input FASTA into individual sequences
 * 2. Reverse each sequence in parallel
 * 3. Display results
 */
workflow {
    // Create input channel from parameter
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Execute splitSequences process
    split_ch = splitSequences(input_ch)
    
    // Flatten channel to process each file individually
    sequences_ch = split_ch.flatten()
    
    // Execute reverseSequences for each split file
    reversed_ch = reverseSequences(sequences_ch)
    
    // Display results (optional - for learning purposes)
    reversed_ch.view { file -> "Processed: ${file}" }
}

/*
 * WORKFLOW SUMMARY:
 * 
 * This pipeline demonstrates:
 * - Parameter definition and usage
 * - Process input/output declarations
 * - Channel operations (fromPath, flatten, view)
 * - Process chaining with channels
 * - Automatic parallelization
 * - Output publishing with publishDir
 */