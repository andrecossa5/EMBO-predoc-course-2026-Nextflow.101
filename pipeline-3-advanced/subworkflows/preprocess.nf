/*
 * Subworkflow: PREPROCESS
 * 
 * This subworkflow demonstrates:
 * - Integration with custom Python scripts
 * - Parameter passing to external processes
 * - Handling of different input formats (mtx vs h5ad)
 * - Sample concatenation and batch integration
 */

include { PP } from '../modules/pp'
include { CONCATENATE } from '../modules/concatenate'

workflow PREPROCESS {
    take:
    manifest   // Channel: manifest.csv file
    
    // Exercise 2: Parameter grid inputs (for future extension)
    // TODO: Define input channels for Exercise 2 with parameter combinations

    main:
    
    // Direct concatenation from manifest
    // This reads all MTX files from manifest and concatenates them
    concatenated_data = CONCATENATE(manifest)
    
    // TODO: Add individual sample preprocessing (PP module) for advanced workflows
    // Call preprocessing module for individual samples (Exercise 2)
    // TODO: Call PP module with parameter grid for individual processing

    emit:
    // Output channels
    concatenated_data = concatenated_data  // Channel: concatenated.h5ad
    // individual_data = // Channel: [sample, processed.h5ad] from PP (for Exercise 2)
}