/*
 * Module: CONCATENATE
 *
 * Concatenates MTX files from manifest using bin/concatenate.py
 * Reads manifest directly and outputs single concatenated h5ad file
 */

process CONCATENATE {
    input:
    // Takes manifest file directly
    path manifest

    output:
    // Outputs single concatenated h5ad file
    path "concatenated.h5ad", emit: adata

    script:
    """
    python ${baseDir}/bin/concatenate.py \\
        --manifest ${manifest} \\
        --output concatenated.h5ad
    """

    stub:
    """
    touch concatenated.h5ad
    """
}
