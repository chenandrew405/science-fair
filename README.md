

## Project Overview

This is a genomics science fair project that analyzes 23andMe genetic data against pharmacogenomics databases to compute risk scores for diseases and drug interactions. The pipeline converts raw genetic data, matches it against external databases (ClinPGx), and uses machine learning to generate risk scores.

## Development Setup

### Virtual Environment

The project uses a Python virtual environment located at `.venv/`:

```bash
# Activate virtual environment
source .venv/bin/activate

# Deactivate when done
deactivate
```

### Key Dependencies

- os
- pandas
- numpy
- re
- sklearn.ensemble.RandomForestRegressor
- sklearn.metrics.mean_squared_error
- sklearn.metrics.r2_score
- warnings
- logging
- sys
- time
- argparse
- fnmatch
- signal
- pathlib.Path
- typing.Optional

## Running the Pipeline

### 1. Data Conversion (src/convert_txt_to_csv.py)

Converts 23andMe .txt files to CSV format:

```bash
# From project root
python src/convert_txt_to_csv.py
```

- Input: `.txt` files in `data/` folder (tab-separated 23andMe genome files)
- Output: `.csv` files in `data/` folder
- Expected columns: rsid, chromosome, position, genotype

### 2. Processing Algorithm (src/analyse_data.py)

Runs the main analysis pipeline:

```bash
# From project root
python src/analyse_data
```

- Input: Converted genome CSV in `data/` folder
- External data: TSV files in `external_data/` (ClinPGx database)
- Output: `data_visualization/final_scores.csv` with risk scores

## Code Architecture

### Core Components

**src/convert_txt_to_csv.py**: Simple conversion utility
- `convert_txt_to_csv()`: Converts single file with column normalization
- `convert_all_txt()`: Batch processes all .txt files
- Handles missing columns gracefully with warnings

**src/analyse_data.py**: Main analysis engine
- `load_user_file()`: Loads converted genome CSV, normalizes genotypes
- `load_external()`: Loads ClinPGx TSV files with error handling
- `compute_score()`: Core ML scoring function using RandomForestClassifier
  - Matches user SNPs (rsid) against external database
  - Encodes genotypes with LabelEncoder
  - Trains RF classifier on matched variants
  - Returns probability-based score weighted by match count
- `main()`: Orchestrates the pipeline
  - Loads user genome and 6 external databases
  - Computes scores for 3 diseases and 3 drugs
  - Exports combined results to CSV

### External Data Files (external_data/)

ClinPGx pharmacogenomics database (CC BY-SA 4.0 license):
- `clinicalVariants.tsv`: Clinical variant annotations
- `var_drug_ann.tsv`: Variant-drug associations
- `var_pheno_ann.tsv`: Variant-phenotype associations
- `var_fa_ann.tsv`: Functional annotation
- `drugLabels.tsv`: Drug label information
- `drugLabels.byGene.tsv`: Gene-specific drug labels

All external TSV files use lowercase column names after loading.

### Key Implementation Details

**Genotype Normalization**: Both user and external data use "NN" for missing/no-call genotypes (converted from "--")

**Scoring Algorithm**:
1. Filter external database for specific disease/drug
2. Match filtered variants to user's genome by rsid
3. Encode all genotypes numerically
4. Train RandomForest on matched variants (positive examples)
5. Predict on user's entire genome
6. Average prediction probabilities
7. Weight by match ratio: `score * (1 + matches/total_snps)`
8. Clamp to [0, 1]

**Column Mappings**:
- User data: rsid, genotype
- External data: "variant annotation id" → rsid, "alleles" → genotype
- Disease filtering: "phenotype" column (case-insensitive substring match)
- Drug filtering: "drug(s)" column (case-insensitive substring match)

## Important Notes

- The project processes personal genomic data - never commit raw .txt or .csv files to git
- External data is licensed under CC BY-SA 4.0 (see external_data/LICENSE.txt)
- The ML scoring approach is experimental/educational - not for clinical use
- Missing external files are handled gracefully with warnings and zero scores
