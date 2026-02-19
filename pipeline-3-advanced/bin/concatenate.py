#!/usr/bin/env python

"""
Concatenation script for scRNA-seq data.
Reads MTX files from manifest and concatenates into single h5ad with sample covariates.
"""

import os
import argparse
import pandas as pd
import scanpy as sc
import anndata as ad


##


def parse_arguments():
    parser = argparse.ArgumentParser(description='Concatenate scRNA-seq samples from manifest')

    # Required arguments
    parser.add_argument('--manifest', '-m', required=True,
                        help='Input manifest CSV file with sample,mtx_dir columns')
    parser.add_argument('--output', '-o', required=True,
                        help='Output concatenated h5ad file')

    # Optional arguments
    parser.add_argument('--batch_key', default='sample',
                        help='Key name for sample batch information (default: sample)')
    parser.add_argument('--join', default='outer',
                        choices=['inner', 'outer'],
                        help='How to join matrices (default: outer)')

    return parser.parse_args()


def main():
    args = parse_arguments()

    # Read manifest file
    print(f"Reading manifest: {args.manifest}")
    try:
        manifest_df = pd.read_csv(args.manifest)
    except Exception as e:
        raise ValueError(f"Error reading manifest file: {e}")

    # Validate manifest format
    required_columns = ['sample', 'adata']
    if not all(col in manifest_df.columns for col in required_columns):
        raise ValueError(f"Manifest must contain columns: {required_columns}")

    print(f"Found {len(manifest_df)} samples in manifest")

    # Load all AnnData objects from MTX directories
    adatas = []
    for i, row in manifest_df.iterrows():
        sample_name = row['sample']
        adata_file = row['adata']

        print(f"Loading sample {i+1}/{len(manifest_df)}: {sample_name} from {adata_file}")

        # Check if AnnData file exists
        if not os.path.exists(adata_file):
            raise FileNotFoundError(f"AnnData file not found: {adata_file}")

        # Load AnnData file using scanpy
        try:
            adata = sc.read(adata_file)
            adata.obs[args.batch_key] = sample_name
            adatas.append(adata)

            print(f"  Loaded {adata.n_obs} cells, {adata.n_vars} genes")

        except Exception as e:
            raise ValueError(f"Error loading adata from {adata_file}: {e}")

    # Concatenate all samples
    print(f"\nConcatenating {len(adatas)} samples...")
    adata_concat = ad.concat(
        adatas,
        join=args.join,
        fill_value=0
    )


    # Print summary statistics
    print(f"\nConcatenation summary:")
    print(f"  Total cells: {adata_concat.n_obs}")
    print(f"  Total genes: {adata_concat.n_vars}")
    print(f"  Samples: {len(manifest_df)}")
    print(f"  Cells per sample:")
    for sample in adata_concat.obs[args.batch_key].unique():
        n_cells = (adata_concat.obs[args.batch_key] == sample).sum()
        print(f"    {sample}: {n_cells}")

    # Save concatenated data
    print(f"\nSaving concatenated data to: {args.output}")
    adata_concat.write(args.output)

    print("Concatenation completed successfully!")


if __name__ == '__main__':
    main()
