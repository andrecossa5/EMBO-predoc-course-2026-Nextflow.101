/*
 * Subworkflow: PROCESS_BAM
 * 
 * This subworkflow demonstrates:
 * - Conditional process execution based on parameters
 * - Single variant calling tool selection
 * - Simplified input/output handling
 */

include { VARIANT_CALLING_GATK } from '../modules/gatk'
include { VARIANT_CALLING_FREEBAYES } from '../modules/freebayes'  
include { VARIANT_CALLING_BCFTOOLS } from '../modules/bcftools'

workflow PROCESS_BAM {
    take:
    bam_ch        // Channel: [sample, bam_file, bam_index]
    reference_ch  // Channel: reference.fa

    main:
    
    // Select single caller based on params.caller
    if (params.caller == "gatk") {
        vcf_results = VARIANT_CALLING_GATK(bam_ch, reference_ch)
    } 
    else if (params.caller == "freebayes") {
        vcf_results = VARIANT_CALLING_FREEBAYES(bam_ch, reference_ch)
    }
    else if (params.caller == "bcftools") {
        vcf_results = VARIANT_CALLING_BCFTOOLS(bam_ch, reference_ch)
    }
    else {
        error "Invalid caller specified: ${params.caller}. Must be one of: gatk, freebayes, bcftools"
    }

    emit:
    vcf_files = vcf_results  // Channel: [sample, vcf]
}