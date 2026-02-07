#!/usr/bin/env nextflow

/*
 * Pipeline 2: Variant Calling with Subworkflows and Conditionals
 * 
 * This pipeline demonstrates advanced Nextflow concepts:
 * - Subworkflow implementation
 * - Conditional process execution 
 * - Module organization
 * - Manifest-driven input handling
 * - Multiple variant calling tools
 */

nextflow.enable.dsl = 2

// Include modules and subworkflows
include { PROCESS_BAM } from './subworkflows/process_bam'

// Pipeline parameters
params.manifest = "$baseDir/data/manifest.csv"
params.reference = "$baseDir/data/reference.fa"  
params.outdir = "results"
params.caller = "gatk"  // Options: gatk, freebayes, bcftools
params.help = false

// Print help message
def helpMessage() {
    log.info"""
    ===================================
    Pipeline 2: Variant Calling Course
    ===================================
    
    Usage:
      nextflow run main.nf [options]
    
    Options:
      --manifest        Input manifest CSV file (default: data/manifest.csv)
      --reference       Reference genome FASTA (default: data/reference.fa)
      --caller          Variant caller to use: gatk, freebayes, bcftools (default: gatk)
      --outdir          Output directory (default: results)
      --help            Show this help message
    
    Input manifest format:
      sample,bam_file,bam_index
      sample1,/path/to/sample1.bam,/path/to/sample1.bam.bai
      sample2,/path/to/sample2.bam,/path/to/sample2.bam.bai
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
    Pipeline 2: Variant Calling Course
    ===================================
    manifest    : ${params.manifest}
    reference   : ${params.reference}
    caller      : ${params.caller}
    outdir      : ${params.outdir}
    ===================================
    """.stripIndent()

    // Create channels from manifest
    Channel
        .fromPath(params.manifest)
        .splitCsv(header: true)
        .map { row -> 
            [
                row.sample, 
                file(row.bam_file), 
                file(row.bam_index)
            ]
        }
        .set { bam_ch }

    // Create reference channel
    reference_ch = Channel.fromPath(params.reference)

    // Execute variant calling subworkflow for each sample
    PROCESS_BAM(bam_ch, reference_ch)

    // Display results summary
    PROCESS_BAM.out.vcf_files
    .view { sample, vcf -> 
        "Sample ${sample}: Generated VCF file: ${vcf.name}"
    }
}