import os
import pandas as pd
import numpy as np
from pathlib import Path

def load_all_csvs(data_dir='data'):
    """Load all aggregated CSV files from the data directory."""
    csv_files = list(Path(data_dir).glob('aggregated_*.csv'))
    print(f"Found {len(csv_files)} CSV files")

    all_data = []
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        df['source_file'] = csv_file.name
        all_data.append(df)

    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"Combined data shape: {combined_df.shape}")
    return combined_df

def min_max_normalize(series):
    """Apply min-max normalization to scale values between 0 and 1."""
    min_val = series.min()
    max_val = series.max()
    if max_val - min_val == 0:
        return pd.Series([0] * len(series))
    return (series - min_val) / (max_val - min_val)

def quantile_normalize(df, columns):
    """Apply quantile normalization across specified columns."""
    # Extract numeric columns
    data = df[columns].values

    # Sort each column
    sorted_data = np.sort(data, axis=0)

    # Calculate mean of each row (quantile)
    mean_quantiles = sorted_data.mean(axis=1)

    # Get ranks for each column
    ranks = data.argsort(axis=0).argsort(axis=0)

    # Replace values with mean quantile values
    normalized_data = np.zeros_like(data)
    for i in range(data.shape[1]):
        normalized_data[:, i] = mean_quantiles[ranks[:, i]]

    return pd.DataFrame(normalized_data, columns=[f'{col}_quantile_norm' for col in columns])

def log_transform(series, method='log1p'):
    """Apply log transformation to handle values including zeros."""
    if method == 'log1p':
        # log(1 + x) - handles zeros well
        return np.log1p(series)
    elif method == 'log2':
        # Add small constant to avoid log(0)
        return np.log2(series + 1e-10)
    elif method == 'log10':
        return np.log10(series + 1e-10)
    else:
        return np.log(series + 1e-10)

def calculate_statistics(df, numeric_cols):
    """Calculate comprehensive statistics for numeric columns."""
    stats_results = []

    for col in numeric_cols:
        data = df[col].dropna()

        stats_dict = {
            'column': col,
            'mean': data.mean(),
            'median': data.median(),
            'mode': data.mode()[0] if len(data.mode()) > 0 else np.nan,
            'std_dev': data.std(),
            'variance': data.var(),
            'min': data.min(),
            'max': data.max(),
            'range': data.max() - data.min(),
            'q1': data.quantile(0.25),
            'q3': data.quantile(0.75),
            'iqr': data.quantile(0.75) - data.quantile(0.25),
            'skewness': data.skew(),
            'kurtosis': data.kurtosis(),
            'count': len(data),
            'sum': data.sum()
        }

        stats_results.append(stats_dict)

    return pd.DataFrame(stats_results)

def main():
    """Main analysis pipeline."""
    print("=" * 80)
    print("COMPREHENSIVE STATISTICAL ANALYSIS PIPELINE")
    print("=" * 80)

    # Ensure output directory exists
    os.makedirs('data_visualization', exist_ok=True)

    # 1. Load all CSV files
    print("\n1. Loading CSV files...")
    df = load_all_csvs('data')
    print(f"Columns: {df.columns.tolist()}")
    print(f"\nFirst few rows:\n{df.head()}")

    # Identify numeric columns (excluding source_file)
    numeric_cols = ['avg_score', 'avg_matched_snps']

    # 2. Apply normalizations and transformations to the full dataset
    print("\n" + "=" * 80)
    print("2. Applying Normalizations and Transformations")
    print("=" * 80)

    # Min-Max Normalization
    for col in numeric_cols:
        df[f'{col}_minmax'] = min_max_normalize(df[col])
        print(f"Applied Min-Max Normalization to {col}")

    # Quantile Normalization (only on avg_score since it's probabilistic)
    quantile_norm_df = quantile_normalize(df, ['avg_score'])
    df = pd.concat([df, quantile_norm_df], axis=1)
    # Also normalize avg_matched_snps separately but we won't use it in final output
    quantile_norm_snps = quantile_normalize(df, ['avg_matched_snps'])
    df = pd.concat([df, quantile_norm_snps], axis=1)
    print(f"Applied Quantile Normalization")

    # Log Transformations
    for col in numeric_cols:
        df[f'{col}_log1p'] = log_transform(df[col], 'log1p')
        df[f'{col}_log2'] = log_transform(df[col], 'log2')
        df[f'{col}_log10'] = log_transform(df[col], 'log10')
        print(f"Applied Log Transformations to {col}")

    # 3. Create final single CSV - ONE ROW PER DISEASE/DRUG with 12 columns
    print("\n" + "=" * 80)
    print("3. Generating Final Single CSV (6 rows, 12 columns)")
    print("=" * 80)

    final_results = []

    for disease in sorted(df['disease_or_drug'].unique()):
        disease_data = df[df['disease_or_drug'] == disease]

        # Calculate statistics for avg_score (primary metric)
        score_data = disease_data['avg_score']

        result_row = {
            'Disease or Drug': disease,
            'Average Risk Score': score_data.mean(),
            'Average Matched Single Nucleotide Polymorphisms': disease_data['avg_matched_snps'].mean(),
            'Median': score_data.median(),
            'Mode': score_data.mode()[0] if len(score_data.mode()) > 0 else 0.0,
            'Minimum': score_data.min(),
            'Maximum': score_data.max(),
            'Range': score_data.max() - score_data.min(),
            'Quantile Normalization': disease_data['avg_score_quantile_norm'].mean(),
            'Log Normalization': disease_data['avg_score_log1p'].mean(),
            'Min-Max Normalization': disease_data['avg_score_minmax'].mean(),
            'Skewness': score_data.skew() if not pd.isna(score_data.skew()) else 0.0
        }

        final_results.append(result_row)

    final_df = pd.DataFrame(final_results)

    # Save final single CSV
    output_file = 'data_visualization/final_analysis.csv'
    final_df.to_csv(output_file, index=False)

    print(f"\nFinal analysis saved to: {output_file}")
    print(f"Shape: {final_df.shape}")
    print(f"Rows: {len(final_df)} (one per disease/drug)")
    print(f"Columns: {len(final_df.columns)}")

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\nGenerated single file: {output_file}")
    print("\n12 Columns:")
    print("  1. Disease or Drug")
    print("  2. Average Risk Score")
    print("  3. Average Matched Single Nucleotide Polymorphisms")
    print("  4. Median")
    print("  5. Mode")
    print("  6. Minimum")
    print("  7. Maximum")
    print("  8. Range")
    print("  9. Quantile Normalization")
    print("  10. Log Normalization")
    print("  11. Min-Max Normalization")
    print("  12. Skewness")

if __name__ == '__main__':
    main()
