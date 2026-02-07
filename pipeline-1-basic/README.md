# Pipeline 1: Basic Nextflow Concepts

This pipeline introduces fundamental Nextflow concepts through a simple two-process workflow that splits FASTA sequences and reverses their content.

## Pipeline Overview

```
Input FASTA → splitSequences → reverseSequences → Results
    ↓               ↓                 ↓
sample.fa      seq_1, seq_2...   reversed files
```

## Quick Start

1. **Run with default parameters**:
   ```bash
   nextflow run main.nf
   ```

2. **Run with custom input**:
   ```bash
   nextflow run main.nf --input /path/to/your/sequences.fa
   ```

3. **Specify output directory**:
   ```bash
   nextflow run main.nf --outdir my_results
   ```

## Input Files

- **sample.fa**: Example FASTA file with 4 sequences (provided)
- You can substitute with any FASTA file

## Output Structure

```
results/
├── split/           # Individual sequence files
│   ├── seq_1
│   ├── seq_2
│   └── ...
├── reversed/        # Reversed sequence files  
│   ├── seq_1_reversed.txt
│   ├── seq_2_reversed.txt
│   └── ...
├── timeline.html    # Execution timeline
├── report.html      # Execution report
└── trace.txt        # Execution trace
```

## Key Concepts Demonstrated

### 1. Process Definition
```nextflow
process splitSequences {
    publishDir "${params.outdir}/split", mode: 'copy'
    
    input:
    path input_file
    
    output:
    path 'seq_*'
    
    script:
    """
    awk '/^>/{f="seq_"++d} {print > f}' < ${input_file}
    """
}
```

### 2. Channel Operations
```nextflow
// Create channel from file parameter
input_ch = Channel.fromPath(params.input, checkIfExists: true)

// Flatten channel for parallel processing
sequences_ch = split_ch.flatten()

// Display results
reversed_ch.view { file -> "Processed: ${file}" }
```

### 3. Workflow Composition
```nextflow
workflow {
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    split_ch = splitSequences(input_ch)
    sequences_ch = split_ch.flatten()
    reversed_ch = reverseSequences(sequences_ch)
    reversed_ch.view { file -> "Processed: ${file}" }
}
```

## Exercises

1. **Modify the script**: Change the `reverseSequences` process to count nucleotides instead
2. **Add parameters**: Add a parameter to control the output file naming
3. **Add validation**: Add input file validation in the workflow
4. **Extend functionality**: Add a third process that counts sequence lengths

## Next Steps

Once comfortable with this pipeline:
1. Move to **Pipeline 2** for subworkflows and conditionals
2. Experiment with different input files
3. Try running with different execution profiles
4. Explore the generated reports (timeline.html, report.html)

## Resources

- [Nextflow Processes](https://nextflow.io/docs/latest/process.html)
- [Channel Operations](https://nextflow.io/docs/latest/channel.html)
- [Configuration Files](https://nextflow.io/docs/latest/config.html)