import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

# --------------------------
# Paths
# --------------------------
BASE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "./data")
USER_FILE = os.path.join(BASE_PATH, "23andme_genome_kat_suricata_v5_full_20171221130201.csv")
DATA_VISUALIZATION = os.path.join(BASE_PATH, "visualization")
DATA_EXTERNAL = os.path.join(BASE_PATH, "external")

EXTERNAL_FILES = {
    "clinvar": "clinicalVariants.tsv",
    "drug_ann": "var_drug_ann.tsv",
    "pheno_ann": "var_pheno_ann.tsv",
    "fa_ann": "var_fa_ann.tsv",
    "drug_labels": "drugLabels.tsv",
    "drug_labels_gene": "drugLabels.byGene.tsv"
}

# --------------------------
# Load user file
# --------------------------
def load_user_file(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"[ERROR] User file not found: {filepath}")

    # 23andMe columns: rsid, chromosome, position, genotype
    df = pd.read_csv(filepath, sep=",", header=None, dtype=str, low_memory=False)
    if df.shape[1] < 4:
        raise ValueError(f"[ERROR] User file has too few columns: {df.shape[1]}")
    df.columns = ["rsid", "chromosome", "position", "genotype"]
    
    # Keep only rsid and genotype, uppercase genotypes
    user_df = df[["rsid", "genotype"]].copy()
    user_df["genotype"] = user_df["genotype"].str.upper().replace("--", "NN").fillna("NN")
    print(f"[INFO] Loaded user file with {len(user_df)} rows, columns: {user_df.columns.tolist()}")
    return user_df

# --------------------------
# Load external TSV
# --------------------------
def load_external(file_key):
    filename = EXTERNAL_FILES[file_key]
    path = os.path.join(DATA_EXTERNAL, filename)
    if not os.path.exists(path):
        print(f"[WARNING] Missing external file (ignored): {filename}")
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df.columns = [c.lower() for c in df.columns]  # normalize columns to lowercase
    print(f"[INFO] Loaded {filename}: {len(df)} rows, columns: {df.columns.tolist()}")
    return df

# --------------------------
# Compute score using RF
# --------------------------
def compute_score(user_df, external_df, rsid_col="rsid", alleles_col="alleles"):
    if external_df.empty:
        return 0.0

    rsid_col = rsid_col.lower()
    alleles_col = alleles_col.lower()
    if rsid_col not in external_df.columns:
        return 0.0
    if alleles_col not in external_df.columns:
        external_df[alleles_col] = "NN"

    # Clean alleles
    external_df[alleles_col] = external_df[alleles_col].str.upper().replace("--", "NN").fillna("NN")
    
    # Match user SNPs
    matches = external_df[external_df[rsid_col].isin(user_df["rsid"])]
    if matches.empty:
        return 0.0

    # Encode genotypes
    le = LabelEncoder()
    combined_genotypes = pd.concat([user_df["genotype"], matches[alleles_col]]).astype(str)
    le.fit(combined_genotypes)

    X_train = le.transform(matches[alleles_col]).reshape(-1, 1)
    y_train = np.ones(len(X_train))  # positive examples
    clf = RandomForestClassifier(n_estimators=50, random_state=42)
    clf.fit(X_train, y_train)

    X_user = le.transform(user_df["genotype"]).reshape(-1, 1)
    try:
        base_score = clf.predict_proba(X_user)[:, 1].mean()  # average probability
    except IndexError:
        base_score = clf.predict_proba(X_user)[:, 0].mean()

    # Reward proportional to number of matches
    score = base_score * (1 + len(matches)/len(user_df))
    return min(score, 1.0)  # clamp to [0,1]

# --------------------------
# Main function
# --------------------------
def main():
    user_df = load_user_file(USER_FILE)

    # Load external data
    clinvar = load_external("clinvar")
    drug_ann = load_external("drug_ann")
    pheno_ann = load_external("pheno_ann")
    fa_ann = load_external("fa_ann")
    drug_labels = load_external("drug_labels")
    drug_labels_gene = load_external("drug_labels_gene")

    # Normalize column names
    pheno_ann.columns = [c.lower() for c in pheno_ann.columns]
    drug_ann.columns = [c.lower() for c in drug_ann.columns]

    # Diseases and drugs
    diseases = ["skin cancer", "diabetes type 2", "lung cancer"]
    drugs = ["carbamazepine", "abacavir", "phenytoin"]

    # Compute disease scores
    disease_scores = {}
    for disease in diseases:
        if "phenotype" in pheno_ann.columns:
            subset = pheno_ann[pheno_ann["phenotype"].str.contains(disease, case=False, na=False)]
            disease_scores[disease] = compute_score(user_df, subset, rsid_col="variant annotation id", alleles_col="alleles")
        else:
            disease_scores[disease] = 0.0

    # Compute drug scores
    drug_scores = {}
    for drug in drugs:
        drug_col = "drug(s)"
        if drug_col in drug_ann.columns:
            subset = drug_ann[drug_ann[drug_col].str.contains(drug, case=False, na=False)]
            drug_scores[drug] = compute_score(user_df, subset, rsid_col="variant annotation id", alleles_col="alleles")
        else:
            drug_scores[drug] = 0.0

    # Combine and export
    final_scores = pd.DataFrame([disease_scores | drug_scores])
    if not os.path.exists(DATA_VISUALIZATION):
        os.makedirs(DATA_VISUALIZATION)
    output_file = os.path.join(DATA_VISUALIZATION, "final_scores.csv")
    final_scores.to_csv(output_file, index=False)
    print(f"[SUCCESS] Exported score table → {output_file}")

if __name__ == "__main__":
    main()
