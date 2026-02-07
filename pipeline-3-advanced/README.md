# Pipeline 3: Advanced scRNA-seq Data Concatenation

This pipeline demonstrates advanced Nextflow concepts for single-cell RNA-seq data handling, focusing on reading multiple MTX files from a manifest and concatenating them into a single analysis-ready dataset.

## Overview

This pipeline processes single-cell RNA-seq data in Matrix Market (mtx) format using a custom Python concatenation script. It reads sample information from a manifest file, loads each sample's MTX data, and combines all samples into a single h5ad file with proper sample annotations.

## Current Functionality

The pipeline implements a streamlined workflow that:
1. **Reads a manifest CSV** containing sample names and MTX directory paths
2. **Loads MTX data** for each sample using scanpy
3. **Adds sample covariates** to distinguish cells from different samples
4. **Concatenates all samples** into a single AnnData object
5. **Outputs a unified h5ad file** ready for downstream analysis

## Input Data

- **Format**: Matrix Market (mtx) directories containing:
  - `matrix.mtx`: Count matrix (Market Matrix format)
  - `barcodes.tsv`: Cell barcodes  
  - `features.tsv` or `genes.tsv`: Gene features/names
- **Manifest**: CSV file with columns:
  - `sample`: Sample identifier
  - `mtx_dir`: Path to MTX directory for that sample

## Pipeline Structure

### Current Implementation

**Objective**: Create a functional pipeline that reads MTX files from manifest and outputs concatenated data.

**Components**:
1. **Main workflow** ([main.nf](main.nf)):
   - Creates manifest channel from CSV file
   - Calls PREPROCESS subworkflow
   - Displays results summary

2. **PREPROCESS subworkflow** ([subworkflows/preprocess.nf](subworkflows/preprocess.nf)):
   - Takes manifest as input
   - Calls CONCATENATE module
   - Returns concatenated dataset

3. **CONCATENATE module** ([modules/concatenate.nf](modules/concatenate.nf)):
   - Executes `bin/concatenate.py` script
   - Takes manifest file as input
   - Outputs single concatenated h5ad file

4. **Concatenation script** ([bin/concatenate.py](bin/concatenate.py)):
   - Reads manifest CSV file
   - Loads each MTX directory using scanpy
   - Adds sample information as `.obs['sample']`
   - Concatenates all samples using `anndata.concat()`
   - Saves unified h5ad file

**Expected Outputs**:
- Single concatenated h5ad file containing all samples
- Sample information stored in `.obs['sample']` column
- Unified gene expression matrix across all samples
- Concatenation metadata in `.uns['concatenation_info']`

### Future Extensions (Exercise 2)

The pipeline structure supports future extensions for:
- Individual sample preprocessing before concatenation
- Parameter grid exploration for preprocessing parameters
- Multiple preprocessing approaches comparison
- Advanced quality control and filtering

## Key Learning Objectives

1. **Manifest-based Input Handling**:
   - Reading CSV manifests in Nextflow
   - Channel creation from file paths
   - Batch processing coordination

2. **MTX Data Integration**:
   - Loading Matrix Market format data
   - Sample annotation and metadata handling
   - Data concatenation strategies

3. **Python-Nextflow Integration**:
   - Custom script integration
   - Parameter passing between Nextflow and Python
   - File I/O coordination

4. **Scalable Pipeline Design**:
   - Modular workflow construction
   - Subworkflow implementation
   - Future-proof architecture

## Usage

### Basic Execution:
```bash
nextflow run main.nf --manifest data/manifest.csv --outdir results
```

### Parameters:
- `--manifest`: Path to manifest CSV file (default: data/manifest.csv)
- `--outdir`: Output directory (default: results)

### Manifest Format:
```csv
sample,mtx_dir
sample_A,data/sample_A_mtx
sample_B,data/sample_B_mtx
sample_C,data/sample_C_mtx
```

### Expected Directory Structure:
```
data/
├── manifest.csv
├── sample_A_mtx/
│   ├── matrix.mtx
│   ├── barcodes.tsv
│   └── features.tsv
├── sample_B_mtx/
│   ├── matrix.mtx
│   ├── barcodes.tsv
│   └── features.tsv
└── sample_C_mtx/
    ├── matrix.mtx
    ├── barcodes.tsv
    └── features.tsv
```

## Script Features

The `concatenate.py` script provides:
- **Comprehensive validation**: Checks manifest format and file existence
- **Flexible gene naming**: Supports gene symbols or gene IDs
- **Join strategies**: Inner or outer joins for different gene sets
- **Unique cell identifiers**: Prevents barcode conflicts across samples
- **Detailed logging**: Progress tracking and summary statistics
- **Error handling**: Clear error messages for troubleshooting

## Output Structure

The concatenated h5ad file contains:
- **X**: Combined expression matrix (cells × genes)
- **obs**: Cell metadata including `sample` column
- **var**: Gene metadata unified across samples
- **uns**: Concatenation metadata and parameters

## File Structure
```
pipeline-3-advanced/
├── main.nf                    # Main pipeline workflow
├── nextflow.config            # Pipeline configuration
├── README.md                  # This documentation
├── bin/
│   ├── preprocessing.py       # Python preprocessing script (for future use)
│   └── concatenate.py         # Python concatenation script
├── subworkflows/
│   └── preprocess.nf         # Preprocessing subworkflow
├── modules/
│   ├── pp.nf                 # Preprocessing module (for future use)
│   └── concatenate.nf        # Concatenation module
└── data/
    ├── manifest.csv          # Sample manifest
    └── sample_*_mtx/         # MTX data directories
```
   - Create sample mtx data in `data/` directory
   - Write manifest.csv with sample entries
   - Run pipeline with default parameters
   - Verify output generation

**Expected Outputs**:
- Individual processed h5ad files with normalized, filtered data
- Concatenated h5ad file with all samples and sample covariates
- QC plots (PCA and UMAP visualizations)  
- Analysis reports

### Exercise 2: Parameter Grid Exploration

**Objective**: Extend the pipeline to systematically explore hyperparameter combinations.

**Tasks**:
1. **Create parameter grid channels**:
   - Convert multi-value parameters to channels
   - Use `combine()` to create parameter combinations
   - Handle channel operations for grid creation

2. **Modify the PP module**:
   - Accept parameter combination inputs
   - Include parameter values in output naming
   - Generate unique outputs for each parameter set

3. **Update configuration**:
   - Define parameter arrays in `nextflow.config`
   - Set resource requirements for grid execution
   - Configure output organization

4. **Implement concatenation for parameter grids**:
   - Group processed samples by parameter combinations
   - Concatenate samples within each parameter set
   - Generate parameter-specific concatenated datasets

5. **Implement grid execution**:
   - Combine sample data with parameter grids
   - Execute preprocessing for all combinations
   - Organize results by parameter sets

**Parameter Grid Variables**:
- `n_neighbors`: [10, 15, 20, 30]
- `n_pcs`: [20, 30, 40, 50] 
- `n_hvgs`: [1000, 2000, 3000, 4000]
- `resolution`: [0.1, 0.3, 0.5, 0.8, 1.0]

**Expected Outputs**:
- Multiple individual processed datasets per sample (one per parameter combination)
- Multiple concatenated datasets (one per parameter combination) with sample covariates
- Systematic directory structure organized by parameters
- Comparative analysis results
- Parameter performance metrics

## Key Learning Objectives

1. **Advanced Channel Operations**:
   - Parameter grid generation
   - Channel combinations and transformations
   - Complex data flow patterns

2. **Parameterization Strategies**:
   - Single vs. multi-value parameter handling
   - Hyperparameter optimization workflows
   - Result organization and comparison

3. **Scalable Pipeline Design**:
   - Efficient resource utilization
   - Parallel parameter exploration
   - Modular workflow construction

4. **Integration with External Scripts**:
   - Python script parameter passing
   - File format handling (mtx to h5ad)
   - Custom tool integration

## Usage

### Exercise 1 (Single Parameters):
```bash
nextflow run main.nf --manifest data/manifest.csv --outdir results_single
```

### Exercise 2 (Parameter Grid):
```bash
nextflow run main.nf --manifest data/manifest.csv --outdir results_grid -profile grid
```

## File Structure
```
pipeline-3-advanced/
├── main.nf                    # Main pipeline workflow
├── nextflow.config            # Pipeline configuration
├── README.md                  # This documentation
├── bin/
│   ├── preprocessing.py       # Python preprocessing script
│   └── concatenate.py         # Python concatenation script
├── subworkflows/
│   └── preprocess.nf         # Preprocessing subworkflow
├── modules/
│   ├── pp.nf                 # Preprocessing module
│   └── concatenate.nf        # Concatenation module
└── data/
    ├── manifest.csv          # Sample manifest (to be created)
    └── sample_*_mtx/         # MTX data directories
```

## Development Notes

This pipeline currently implements the core concatenation functionality. Future development can extend it with:

### Planned Extensions:
1. **Individual Sample Preprocessing**:
   - Add PP module implementation for quality control
   - Include filtering, normalization, and dimensionality reduction
   - Generate per-sample QC reports

2. **Parameter Grid Exploration**:
   - Implement parameter combinations for preprocessing
   - Support hyperparameter optimization workflows
   - Compare different preprocessing strategies

3. **Advanced Features**:
   - Batch effect correction methods
   - Integration algorithms (Harmony, Scanorama, etc.)
   - Automated quality assessment

### Technical Considerations:
- **Memory Management**: Current configuration allocates 32GB for concatenation
- **Scalability**: Designed to handle multiple large samples efficiently
- **Error Recovery**: Comprehensive validation and error reporting
- **Reproducibility**: Container support and parameter logging

## Troubleshooting

### Common Issues:

1. **MTX Directory Not Found**:
   - Check paths in manifest.csv are correct
   - Ensure MTX directories contain required files (`matrix.mtx`, `barcodes.tsv`, `features.tsv`)

2. **Memory Issues**:
   - Increase memory allocation in nextflow.config
   - Consider processing fewer samples per run

3. **Gene Name Conflicts**:
   - Use `--gene_names gene_ids` for consistent gene identifiers
   - Check that all samples have compatible gene annotations

4. **File Format Issues**:
   - Ensure MTX files are properly formatted
   - Check that barcodes.tsv and features.tsv match matrix dimensions

### Debug Mode:
```bash
# Run with detailed logging
nextflow run main.nf --manifest data/manifest.csv -with-report -with-timeline

# Check concatenate script directly
python bin/concatenate.py --manifest data/manifest.csv --output test.h5ad
```

## Testing

### Quick Test:
1. Create test MTX directories in `data/`
2. Update manifest.csv with correct paths
3. Run pipeline: `nextflow run main.nf`
4. Check output in `results/concatenated/`

### Expected Results:
- Single h5ad file with all samples
- Sample annotations in `.obs['sample']`
- Preserved gene information in `.var`
- Concatenation metadata in `.uns`

## Contributing

This pipeline serves as a foundation for advanced scRNA-seq workflows. Contributions welcome for:
- Additional preprocessing modules
- Integration with other single-cell tools  
- Performance optimizations
- Documentation improvements
- Test data and examples

## Contact

For questions about this pipeline or the EMBO Nextflow course, please contact the course organizers.