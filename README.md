# EMBO Nextflow Course 2026

A comprehensive course covering Nextflow workflow development from basic concepts to advanced patterns, specifically designed for bioinformatics applications.

## Course Overview

This repository contains three progressive pipelines that teach Nextflow development:

### 📚 [Pipeline 1: Basic Concepts](pipeline-1-basic/)
- **Objective**: Learn fundamental Nextflow syntax and concepts
- **Topics**: Processes, channels, operators, basic workflows
- **Use Case**: Simple FASTA file processing and analysis
- **Key Skills**: DSL2 syntax, channel operations, process definition

### 🔬 [Pipeline 2: Intermediate Patterns](pipeline-2-intermediate/)
- **Objective**: Master subworkflows, modules, and conditional execution
- **Topics**: Modular design, variant calling workflows, parameter handling
- **Use Case**: Multi-tool variant calling pipeline (GATK, FreeBayes, BCFtools)
- **Key Skills**: Subworkflows, modules, conditional processes, manifest handling

### 🚀 [Pipeline 3: Advanced Integration](pipeline-3-advanced/)
- **Objective**: Implement advanced data handling and integration patterns
- **Topics**: Multi-sample integration, MTX file handling, data concatenation
- **Use Case**: Single-cell RNA-seq data concatenation from multiple samples
## Learning Path

```
Pipeline 1 (Basic) → Pipeline 2 (Intermediate) → Pipeline 3 (Advanced)
      ↓                      ↓                        ↓
  Fundamentals         Modular Design          Integration Patterns
  Processes           Subworkflows            Multi-sample Handling
  Channels            Modules                 External Scripts
  Operators           Conditionals            Data Aggregation
```

## Prerequisites

- Basic command line experience
- Understanding of bioinformatics file formats
- Basic Python knowledge (for Pipeline 3)
- Nextflow installed (version ≥ 1.0)

## Getting Started

1. **Clone this repository**:
   ```bash
   git clone <repository-url>
   cd nextflow_course_2026
   ```

2. **Start with Pipeline 1**:
   ```bash
   cd pipeline-1-basic
   cat README.md  # Read the instructions
   nextflow run main.nf
   ```

3. **Progress through each pipeline**:
   - Complete exercises in order
   - Review README files for detailed instructions
   - Test your understanding with the provided examples

## Course Structure

```
nextflow_course_2026/
├── README.md                 # This overview
├── .gitignore               # Git ignore patterns
├── pipeline-1-basic/        # Basic Nextflow concepts
│   ├── main.nf
│   ├── nextflow.config
│   ├── README.md
│   └── data/
├── pipeline-2-intermediate/  # Subworkflows and modules
│   ├── main.nf
│   ├── nextflow.config
│   ├── README.md
│   ├── modules/
│   ├── subworkflows/
│   └── data/
└── pipeline-3-advanced/     # Advanced integration patterns
    ├── main.nf
    ├── nextflow.config
    ├── README.md
    ├── bin/
    ├── modules/
    ├── subworkflows/
    └── data/
```

## Key Learning Outcomes

By completing this course, you will be able to:

✅ **Design and implement** complete Nextflow workflows  
✅ **Structure projects** using modules and subworkflows  
✅ **Handle complex data flows** with channels and operators  
✅ **Integrate external tools** and custom scripts  
✅ **Manage parameters** and configuration effectively  
✅ **Debug and troubleshoot** workflow issues  
✅ **Follow best practices** for reproducible workflows  
✅ **Scale workflows** for different execution environments  

## Support and Resources

### Course Resources:
- 📖 **README files** in each pipeline directory
- 💡 **Inline comments** explaining key concepts
- 🔧 **Working examples** with sample data
- 🐛 **Troubleshooting guides** for common issues

### External Resources:
- [Official Nextflow Documentation](https://nextflow.io/docs/latest/)
- [Nextflow Patterns](https://nextflow-io.github.io/patterns/)
- [nf-core Guidelines](https://nf-co.re/developers/guidelines)
- [Nextflow Community](https://community.nextflow.io/)

## License

This course is released under the MIT License.

## Acknowledgments

- **EMBO** for supporting bioinformatics training
- **Nextflow team** for creating an amazing workflow system  
- **nf-core community** for establishing best practices
- **Course contributors** and beta testers

---

**Happy Learning!** 🧬🔬

*Start your Nextflow journey with Pipeline 1 and work your way through to become a workflow expert.*