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
params.input = "${baseDir}/data/sample.fa"
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
    path ('seq_*'), emit: sequences

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
    path "${sequence_file.baseName}_reversed.txt", emit: reverse_equences

    script:
    """
    echo "Processing ${sequence_file}..."
    awk '/^>/ {print} !/^>/ {for(i=length;i>=1;i--) printf("%c", substr(\$0,i,1)); print ""}' ${sequence_file} > ${sequence_file.baseName}_reversed.txt
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
    input_ch = Channel.fromPath(params.input)

    // Split
    splitSequences(input_ch)

    // Reverse
    reverseSequences(splitSequences.out.sequences.flatten())
}
