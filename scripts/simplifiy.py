import pandas as pd
import os

# Get the script directory and navigate to project root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_FOLDER = os.path.join(PROJECT_ROOT, "data")
INPUT_FOLDER = os.path.join(DATA_FOLDER, "visualization")
OUTPUT_FOLDER = os.path.join(DATA_FOLDER, "results")
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Ensure input folder exists
if not os.path.exists(INPUT_FOLDER):
    raise FileNotFoundError(f"[ERROR] Input folder does not exist: {INPUT_FOLDER}")

# Collect all final_scores_*.csv files in the visualization folder
data_files = sorted([f for f in os.listdir(INPUT_FOLDER) if f.startswith("final_scores_") and f.endswith('.csv')])

if not data_files:
    print(f"[WARNING] No final_scores_*.csv files found in {INPUT_FOLDER}")
    print(f"  Available files: {os.listdir(INPUT_FOLDER)}")
else:
    print(f"[INFO] Found {len(data_files)} final_scores files to process")

simplified_rows = []

for filename in data_files:
    file_path = os.path.join(INPUT_FOLDER, filename)
    try:
        df = pd.read_csv(file_path)
        print(f"[INFO] Processing {filename}: {len(df)} rows, columns: {df.columns.tolist()}")
        # Only keep rows where 'score' and 'matched_snps' columns exist
        if {'sample', 'name', 'score', 'matched_snps'}.issubset(df.columns):
            for _, row in df.iterrows():
                simplified_rows.append({
                    'sample': row['sample'],
                    'disease': row['name'],
                    'score': row['score'],
                    'matched_snps': row['matched_snps']
                })
            
            # Generate per-dataset aggregated scores (averages for this specific CSV)
            dataset_df = df[['name', 'score', 'matched_snps']].copy()
            aggregated = dataset_df.groupby('name').agg({
                'score': 'mean',
                'matched_snps': 'mean'
            }).round(6)
            
            aggregated.columns = ['avg_score', 'avg_matched_snps']
            aggregated = aggregated.reset_index()
            aggregated = aggregated.rename(columns={'name': 'disease_or_drug'})
            aggregated = aggregated.sort_values('disease_or_drug').reset_index(drop=True)
            
            # Save per-dataset aggregated scores
            base_name = filename.replace('final_scores_', 'aggregated_').replace('.csv', '')
            aggregated_path = os.path.join(OUTPUT_FOLDER, f"{base_name}.csv")
            aggregated.to_csv(aggregated_path, index=False)
            print(f"  [INFO] Dataset aggregated scores written to {os.path.basename(aggregated_path)}")
        else:
            print(f"[WARNING] {filename} missing required columns")
    except Exception as e:
        print(f"[WARNING] Could not process {filename}: {e}")

# Create DataFrame and save to results/simplified_scores.csv
if simplified_rows:
    simplified_df = pd.DataFrame(simplified_rows)
    output_path = os.path.join(OUTPUT_FOLDER, "simplified_scores.csv")
    simplified_df.to_csv(output_path, index=False)
    print(f"\n[INFO] Simplified scores written to {output_path}")
    print(f"[INFO] Total rows: {len(simplified_df)}")
else:
    print("[ERROR] No data rows to write!")