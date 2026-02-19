/*
 * Subworkflow: PREPROCESS
 *
 * This subworkflow demonstrates:
 * - Integration with custom Python scripts
 * - Parameter passing to external processes
 * - Handling of different input formats (mtx vs h5ad)
 * - Sample concatenation and batch integration
 */

include { CONCATENATE } from '../modules/concatenate'
// include { PP } from '../modules/pp'

workflow PREPROCESS {
    take:
    manifest // Channel: manifest.csv file

    main:

    // Direct concatenation from manifest
    CONCATENATE(manifest)

    emit:
    preprocessed_output = CONCATENATE.out.adata
}
