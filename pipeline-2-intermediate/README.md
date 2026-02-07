# Pipeline 2: Variant Calling with Subworkflows and Conditionals

This pipeline demonstrates intermediate Nextflow concepts through a variant calling workflow that supports multiple tools and conditional execution.

## Learning Objectives

After completing this pipeline, you will understand:

- ✅ **Subworkflow organization** and composition
- ✅ **Conditional process execution** based on parameters
- ✅ **Module system** for reusable components
- ✅ **Manifest-driven input** handling
- ✅ **Multiple tool integration** and comparison
- ✅ **Advanced channel operations** and data flow

## Pipeline Overview

```
Manifest CSV → BAM files → Subworkflow (PROCESS_BAM) → VCF outputs
                ↓                     ↓
          sample info            Conditional execution:
                                 - GATK HaplotypeCaller
                                 - FreeBayes  
                                 - BCFtools mpileup/call
```

## Architecture

```
main.nf                          # Main workflow entry point
├── subworkflows/
│   └── process_bam.nf          # BAM processing subworkflow
├── modules/                    # Individual tool modules
│   ├── gatk.nf                # GATK HaplotypeCaller
│   ├── freebayes.nf           # FreeBayes variant caller
│   └── bcftools.nf            # BCFtools variant caller
└── data/                      # Example input files
    ├── manifest.csv           # Sample manifest
    ├── reference.fa           # Reference genome
    └── *.bam + *.bam.bai     # Sample BAM files
```

## Quick Start

1. **Run all variant callers** (default):
   ```bash
   nextflow run main.nf
   ```

2. **Run specific caller**:
   ```bash
   nextflow run main.nf --variant_caller gatk
   nextflow run main.nf --variant_caller freebayes
   nextflow run main.nf --variant_caller bcftools
   ```

3. **Use custom manifest**:
   ```bash
   nextflow run main.nf --manifest /path/to/your/manifest.csv
   ```

## Input Files

### Manifest Format (CSV)
```csv
sample,bam_file,bam_index
sample_A,data/sample_A.bam,data/sample_A.bam.bai
sample_B,data/sample_B.bam,data/sample_B.bam.bai
```

### Required Files
- **manifest.csv**: Sample information and file paths
- **reference.fa**: Reference genome FASTA
- **\*.bam**: Aligned sequence data (BAM format)
- **\*.bam.bai**: BAM index files

## Output Structure

```
results/
├── sample_A/
│   ├── gatk/
│   │   └── sample_A.gatk.vcf.gz
│   ├── freebayes/
│   │   └── sample_A.freebayes.vcf.gz
│   └── bcftools/
│       └── sample_A.bcftools.vcf.gz
├── sample_B/
│   └── ... (same structure)
├── timeline.html         # Execution timeline
├── report.html          # Execution report
├── trace.txt            # Execution trace
└── dag.svg              # Workflow DAG visualization
```

## Key Concepts Demonstrated

### 1. Subworkflow Definition
```nextflow
// subworkflows/process_bam.nf
workflow PROCESS_BAM {
    take:
    bam_ch        // Input channel
    reference_ch  // Reference genome
    caller        // Caller selection parameter

    main:
    // Conditional execution logic
    if (caller == "gatk" || caller == "all") {
        gatk_vcf = VARIANT_CALLING_GATK(bam_ch, reference_ch)
    }
    // ... more conditions

    emit:
    vcf_files = grouped_vcf  // Output channel
}
```

### 2. Conditional Process Execution
```nextflow
// Execute different callers based on parameter
if (caller == "gatk" || caller == "all") {
    gatk_vcf = VARIANT_CALLING_GATK(bam_ch, reference_ch)
}

if (caller == "freebayes" || caller == "all") {
    freebayes_vcf = VARIANT_CALLING_FREEBAYES(bam_ch, reference_ch)
}
```

### 3. Module Organization
```nextflow
// modules/gatk.nf
process VARIANT_CALLING_GATK {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}/gatk", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai)
    path reference
    
    output:
    tuple val(sample), path("${sample}.gatk.vcf.gz")
    
    script:
    // GATK-specific variant calling logic
}
```

### 4. Manifest-Driven Input
```nextflow
Channel
    .fromPath(params.manifest)
    .splitCsv(header: true)
    .map { row -> 
        [row.sample, file(row.bam_file), file(row.bam_index)]
    }
    .set { bam_ch }
```

## Advanced Features

### Execution Profiles
- **standard**: Local execution
- **docker**: Containerized execution with specific tool containers
- **cluster**: SLURM cluster execution

### Resource Management
```nextflow
process {
    withName: 'VARIANT_CALLING_GATK' {
        cpus = 2
        memory = '8 GB' 
        time = '4h'
    }
}
```

### Container Integration
```nextflow
withName: 'VARIANT_CALLING_GATK' {
    container = 'broadinstitute/gatk:4.4.0.0'
}
```

## Exercises

1. **Add a new variant caller**: Create a module for VarScan2
2. **Implement quality filtering**: Add VCF filtering steps
3. **Compare outputs**: Add a process to compare variant calls between tools
4. **Add annotations**: Integrate variant annotation tools

## Common Patterns Demonstrated

### Channel Mixing
```nextflow
vcf_results = vcf_results.mix(
    gatk_vcf.map { sample, vcf -> [sample, [vcf], "GATK"] }
)
```

### Channel Grouping
```nextflow
grouped_vcf = vcf_results
    .groupTuple(by: 0)
    .map { sample, vcf_lists, callers ->
        [sample, vcf_lists.flatten(), callers]
    }
```

### Tag Usage for Logging
```nextflow
process VARIANT_CALLING_GATK {
    tag "${sample}"  // Shows sample name in logs
    // ...
}
```

## Exercise

1. Add a process that merge each all vcf files into a single one
2. Modify the pipeline to scatter variant calling across chromosomes (advanced)

## Next Steps
Move to **Pipeline 3**.

## Resources

- [Nextflow Modules](https://nextflow.io/docs/latest/module.html)
- [Subworkflows](https://nextflow.io/docs/latest/workflow.html)
- [Conditional Execution](https://nextflow.io/docs/latest/process.html#when)
- [nf-core Modules](https://nf-co.re/modules)