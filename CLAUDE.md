
# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

- pandas 2.3.3 - Data manipulation
- numpy 2.3.4 - Numerical operations
- scikit-learn 1.7.2 - Machine learning (RandomForestClassifier)
- torch 2.9.1 - Deep learning (may be used for future extensions)
- matplotlib 3.10.7 - Visualization

## Running the Pipeline

### 1. Data Conversion (data/data.py)

Converts 23andMe .txt files to CSV format:

```bash
# From project root
python data/data.py
```

- Input: `.txt` files in `data/` folder (tab-separated 23andMe genome files)
- Output: `.csv` files in `processing_algorithm/` folder
- Expected columns: rsid, chromosome, position, genotype

### 2. Processing Algorithm (processing_algorithm/processing.py)

Runs the main analysis pipeline:

```bash
# From project root
python processing_algorithm/processing.py
```

- Input: Converted genome CSV in `processing_algorithm/` folder
- External data: TSV files in `external_data/` (ClinPGx database)
- Output: `data_visualization/final_scores.csv` with risk scores

## Code Architecture

### Data Flow

```
data/*.txt → data/data.py → processing_algorithm/*.csv →
processing_algorithm/processing.py → data_visualization/final_scores.csv
```

### Core Components

**data/data.py**: Simple conversion utility
- `convert_txt_to_csv()`: Converts single file with column normalization
- `convert_all_txt()`: Batch processes all .txt files
- Handles missing columns gracefully with warnings

**processing_algorithm/processing.py**: Main analysis engine
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
