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
params.manifest = "${baseDir}/data/manifest.csv"
params.outdir = "results"

/*
 * MAIN WORKFLOW
 */
workflow {

  // Call preprocessing subworkflow
  PREPROCESS(params.manifest)

  // Display results
  PREPROCESS.out.preprocessed_output.view()
}
