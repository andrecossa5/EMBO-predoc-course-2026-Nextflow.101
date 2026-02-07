#!/usr/bin/env nextflow

/*
 * Pipeline 3: Advanced scRNA-seq Preprocessing with Parameter Grids
 * 
 * This pipeline demonstrates advanced Nextflow concepts:
 * - Parameter grids and hyperparameter exploration
 * - Matrix market (mtx) file input handling
 * - Advanced subworkflow implementation
 * - Channel operations for parameter combinations
 */

nextflow.enable.dsl = 2

// Include subworkflows
include { PREPROCESS } from './subworkflows/preprocess'

// Pipeline parameters
params.manifest = "$baseDir/data/manifest.csv"
params.outdir = "results"
params.help = false

// Single-value parameters (Exercise 1)
// TODO: Add default values for all preprocessing parameters

// Multi-value parameters (Exercise 2) 
// TODO: Add parameter arrays for hyperparameter grid

// Print help message
def helpMessage() {
    log.info"""
    ===================================
    Pipeline 3: Advanced scRNA-seq Course
    ===================================
    
    Usage:
      nextflow run main.nf [options]
    
    Options:
      --manifest        Input manifest CSV file (default: data/manifest.csv)
      --outdir          Output directory (default: results)
      --help            Show this help message
    
    Input manifest format:
      sample,mtx_dir
      sample1,/path/to/sample1_mtx_dir/
      sample2,/path/to/sample2_mtx_dir/
      
    The pipeline will:
      1. Read all MTX directories from the manifest
      2. Preprocess each sample with specified parameters
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * MAIN WORKFLOW
 */
workflow {
    // Print pipeline information
    log.info """
    ===================================
    Pipeline 3: Advanced scRNA-seq Course
    ===================================
    manifest    : ${params.manifest}
    outdir      : ${params.outdir}
    ===================================
    """.stripIndent()

    // Call preprocessing subworkflow
    PREPROCESS(params.manifest)
    
    // Display results
    PREPROCESS.out.concatenated_data
        .view { "Concatenated dataset created: ${it}" }
}