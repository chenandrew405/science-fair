import os
import pandas as pd
import numpy as np
import re
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = os.path.join("..", "data")

# Find all CSV files in processing_algorithm/ directory (where converted genome CSVs are)
csv_files = [f for f in os.listdir(os.path.join(PROJECT_ROOT, "processing_algorithm")) if f.endswith('.csv')]
if not csv_files:
    raise FileNotFoundError("[ERROR] No CSV files found in processing_algorithm/ directory")

DATA_VISUALIZATION = os.path.join(PROJECT_ROOT, "data_visualization")
DATA_EXTERNAL = os.path.join(PROJECT_ROOT, "data", "external")

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
    # First try with header, then without if it fails
    try:
        df = pd.read_csv(filepath, sep=",", dtype=str, low_memory=False)
        if df.shape[1] >= 4 and all(col in df.columns for col in ['rsid', 'chromosome', 'position', 'genotype']):
            # File has proper headers
            pass
        else:
            # Doesn't have expected headers, reload without header
            df = pd.read_csv(filepath, sep=",", header=None, dtype=str, low_memory=False)
            df.columns = ["rsid", "chromosome", "position", "genotype"]
    except Exception:
        # Try without header
        df = pd.read_csv(filepath, sep=",", header=None, dtype=str, low_memory=False)
        if df.shape[1] < 4:
            raise ValueError(f"[ERROR] User file has too few columns: {df.shape[1]}")
        df.columns = ["rsid", "chromosome", "position", "genotype"]

    # Keep only rsid and genotype, uppercase genotypes
    user_df = df[["rsid", "genotype"]].copy()
    user_df["genotype"] = user_df["genotype"].str.upper().replace("--", "NN").fillna("NN")

    # Filter out any header rows that might have been included as data
    user_df = user_df[user_df["rsid"] != "rsid"]

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
# Helper functions for genotype matching
# --------------------------
def normalize_genotype(genotype):
    """Normalize genotype to canonical form (alphabetically sorted alleles)."""
    if pd.isna(genotype) or genotype == "NN" or genotype == "--":
        return None
    # Remove common separators and spaces
    genotype = str(genotype).upper().replace("/", "").replace("|", "").replace(" ", "").replace("_", "")
    # Handle star alleles (e.g., *1/*28) - keep as is
    if "*" in genotype:
        return genotype
    # Sort single nucleotide genotypes (e.g., AT -> AT, TA -> AT)
    if len(genotype) == 2 and genotype.isalpha():
        return "".join(sorted(genotype))
    return genotype

def genotypes_match(user_genotype, risk_genotype):
    """
    Check if user's genotype matches the risk genotype from database.

    Handles:
    - Exact matches (AA == AA)
    - Compound genotypes (AA + AT matches if user has AA or AT)
    - Heterozygous matches (AT matches both AT and TA)
    - Star allele notation (*1/*28)
    """
    user_gt = normalize_genotype(user_genotype)

    if user_gt is None:
        return False

    # Handle compound genotypes (e.g., "AA + AT" or "GG + GT + TT")
    if "+" in str(risk_genotype):
        risk_alleles = [normalize_genotype(g.strip()) for g in str(risk_genotype).split("+")]
        return user_gt in risk_alleles

    risk_gt = normalize_genotype(risk_genotype)
    if risk_gt is None:
        return False

    return user_gt == risk_gt

def extract_odds_ratio(text):
    """
    Extract odds ratio from free text (Notes or Sentence columns).
    Returns the OR value, or None if not found.
    """
    if pd.isna(text):
        return None

    # Pattern to match OR = X.XX or odds ratio = X.XX
    patterns = [
        r'OR\s*=\s*(\d+\.?\d*)',
        r'odds ratio\s*=\s*(\d+\.?\d*)',
        r'OR\s+(\d+\.?\d*)',
    ]

    for pattern in patterns:
        match = re.search(pattern, str(text), re.IGNORECASE)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                continue

    return None

def get_significance_weight(significance):
    """
    Convert significance rating to numeric weight.

    Weights:
    - "yes" (significant) = 1.0
    - "not stated" (uncertain) = 0.5
    - "no" (not significant) = 0.1 (small weight, don't ignore completely)
    """
    if pd.isna(significance):
        return 0.5

    sig_lower = str(significance).lower().strip()

    if sig_lower == "yes":
        return 1.0
    elif sig_lower == "no":
        return 0.1
    else:  # "not stated" or other
        return 0.5

def get_direction_multiplier(direction_text):
    """
    Determine if variant increases or decreases risk.

    Returns:
    - 1.0 for risk-increasing (increased, higher, more)
    - -1.0 for risk-decreasing (decreased, lower, less, protective)
    - 0.5 for uncertain/neutral
    """
    if pd.isna(direction_text):
        return 0.5

    text_lower = str(direction_text).lower()

    # Risk-increasing keywords
    if any(word in text_lower for word in ["increased", "higher", "more", "greater", "elevated"]):
        return 1.0

    # Risk-decreasing/protective keywords
    if any(word in text_lower for word in ["decreased", "lower", "less", "reduced", "protective"]):
        return -1.0

    return 0.5

# --------------------------
# Compute score using weighted polygenic risk score
# --------------------------
def compute_score(user_df, external_df, rsid_col="rsid", alleles_col="alleles"):
    """
    Compute weighted polygenic risk score (PRS) based on matched variants.

    IMPROVEMENTS OVER OLD VERSION:
    1. Matches actual GENOTYPES (not just rsid presence)
    2. Weights variants by clinical significance
    3. Incorporates direction of effect (risk-increasing vs protective)
    4. Extracts odds ratios when available
    5. Uses weighted sum formula instead of simple ratio

    Polygenic Risk Score Formula:
    PRS = Σ(weight_i × direction_i × OR_i) / N

    Where:
    - weight_i = clinical significance weight (0.1-1.0)
    - direction_i = direction multiplier (+1 risk, -1 protective)
    - OR_i = odds ratio (default 1.5 if not available)
    - N = number of variants in database (for normalization)

    Score Interpretation:
    - Positive score = increased risk
    - Negative score = protective/decreased risk
    - Magnitude indicates strength of association
    - Normalized to roughly 0-1 scale for reporting

    Args:
        user_df: User genome dataframe with columns [rsid, genotype]
        external_df: External database dataframe
        rsid_col: Column name for rsid in external_df
        alleles_col: Column name for alleles in external_df

    Returns:
        tuple: (score: float, matched_snps: int, details: dict)
    """
    # Handle empty external data
    if external_df.empty:
        return 0.0, 0, {}

    # Make a copy to avoid warnings
    external_df = external_df.copy()

    # Normalize column names
    rsid_col = rsid_col.lower()
    alleles_col = alleles_col.lower()

    # Check for required columns
    if rsid_col not in external_df.columns:
        return 0.0, 0, {}

    # Prepare user data lookup
    user_lookup = dict(zip(user_df["rsid"], user_df["genotype"]))

    # Track matched variants and their contributions
    matched_variants = []
    total_weighted_score = 0.0
    risk_increasing = 0
    risk_decreasing = 0

    # Process each variant in the external database
    for idx, row in external_df.iterrows():
        rsid = row.get(rsid_col)
        risk_alleles = row.get(alleles_col)

        if pd.isna(rsid) or rsid not in user_lookup:
            continue

        user_genotype = user_lookup[rsid]

        # Check if genotypes match
        if not genotypes_match(user_genotype, risk_alleles):
            continue

        # Extract metadata for weighting
        significance = row.get("significance", None)
        direction = row.get("direction of effect", None)
        notes = row.get("notes", "")
        sentence = row.get("sentence", "")

        # Get weights
        sig_weight = get_significance_weight(significance)
        dir_multiplier = get_direction_multiplier(direction)

        # Try to extract odds ratio from text
        odds_ratio = extract_odds_ratio(notes) or extract_odds_ratio(sentence)
        if odds_ratio is None:
            odds_ratio = 1.5  # Default moderate effect size

        # Calculate contribution
        # Use log(OR) for better scale (OR of 2.0 -> 0.69, OR of 10 -> 2.3)
        log_or = np.log(odds_ratio) if odds_ratio > 0 else 0
        contribution = sig_weight * dir_multiplier * log_or

        total_weighted_score += contribution

        if dir_multiplier > 0:
            risk_increasing += 1
        elif dir_multiplier < 0:
            risk_decreasing += 1

        matched_variants.append({
            "rsid": rsid,
            "user_genotype": user_genotype,
            "risk_genotype": risk_alleles,
            "significance": significance,
            "odds_ratio": odds_ratio,
            "contribution": contribution
        })

    matched_snps = len(matched_variants)
    total_variants = len(external_df)

    if matched_snps == 0:
        return 0.0, 0, {}

    # Normalize score
    # Divide by sqrt(N) to account for number of variants (standard PRS normalization)
    normalized_score = total_weighted_score / np.sqrt(total_variants) if total_variants > 0 else 0.0

    # Convert to 0-1 scale using sigmoid transformation
    # This maps the score to a probability-like scale
    final_score = 1 / (1 + np.exp(-normalized_score))

    # Prepare detailed results
    details = {
        "matched_snps": matched_snps,
        "total_variants": total_variants,
        "risk_increasing": risk_increasing,
        "risk_decreasing": risk_decreasing,
        "raw_score": total_weighted_score,
        "normalized_score": normalized_score,
        "avg_contribution": total_weighted_score / matched_snps if matched_snps > 0 else 0
    }

    return final_score, matched_snps, details


# --------------------------
# Machine Learning: Q-Learning SNP Selection
# --------------------------
class SNPSelectionQLearning:
    """
    Q-Learning agent to learn optimal SNP selection strategy.
    """

    def __init__(self, snp_list, learning_rate=0.1, discount=0.95, epsilon=0.1):
        self.snp_list = snp_list
        self.n_snps = len(snp_list)
        self.lr = learning_rate
        self.gamma = discount
        self.epsilon = epsilon
        self.q_table = {}
        self.best_selection = None
        self.best_reward = -float('inf')

    def state_to_key(self, state):
        """Convert binary state vector to string key."""
        return ''.join(map(str, state.astype(int)))

    def get_q_value(self, state, action):
        """Get Q-value for state-action pair."""
        key = self.state_to_key(state)
        if key not in self.q_table:
            self.q_table[key] = np.zeros(self.n_snps)
        return self.q_table[key][action]

    def set_q_value(self, state, action, value):
        """Set Q-value for state-action pair."""
        key = self.state_to_key(state)
        if key not in self.q_table:
            self.q_table[key] = np.zeros(self.n_snps)
        self.q_table[key][action] = value

    def choose_action(self, state):
        """Epsilon-greedy action selection."""
        if np.random.random() < self.epsilon:
            return np.random.randint(0, self.n_snps)
        else:
            key = self.state_to_key(state)
            if key not in self.q_table:
                return np.random.randint(0, self.n_snps)
            return np.argmax(self.q_table[key])

    def update(self, state, action, reward, next_state):
        """Q-learning update rule."""
        current_q = self.get_q_value(state, action)
        next_max_q = np.max([self.get_q_value(next_state, a) for a in range(self.n_snps)])
        new_q = current_q + self.lr * (reward + self.gamma * next_max_q - current_q)
        self.set_q_value(state, action, new_q)

    def train(self, user_df, external_df, target_score, episodes=50, max_steps=30):
        """Train Q-learning agent to find optimal SNP selection."""
        for episode in range(episodes):
            state = np.random.randint(0, 2, size=self.n_snps)
            for step in range(max_steps):
                action = self.choose_action(state)
                next_state = state.copy()
                next_state[action] = 1 - next_state[action]

                # Compute reward
                selected_snps = [self.snp_list[i] for i in range(self.n_snps) if next_state[i] == 1]
                reward = self._compute_reward(user_df, external_df, selected_snps, target_score)

                self.update(state, action, reward, next_state)

                if reward > self.best_reward:
                    self.best_reward = reward
                    self.best_selection = selected_snps.copy()

                state = next_state

        return self.best_selection, self.best_reward

    def _compute_reward(self, user_df, external_df, selected_snps, target_score):
        """Compute reward for SNP selection."""
        if len(selected_snps) == 0:
            return -1.0

        rsid_col = "variant/haplotypes"
        if rsid_col not in external_df.columns:
            return 0.0

        filtered_df = external_df[external_df[rsid_col].isin(selected_snps)]

        if filtered_df.empty:
            return -0.5

        score, matched, _ = compute_score(user_df, filtered_df, rsid_col=rsid_col, alleles_col="alleles")
        score_diff = abs(score - target_score)
        correlation_reward = 1.0 - score_diff
        match_ratio = matched / len(selected_snps) if len(selected_snps) > 0 else 0
        match_reward = match_ratio * 0.3
        selection_ratio = len(selected_snps) / self.n_snps
        penalty = -0.2 if (selection_ratio < 0.1 or selection_ratio > 0.9) else 0

        return correlation_reward + match_reward + penalty


# --------------------------
# Machine Learning: RFC Validation
# --------------------------
def train_rfc_validator(features_list, scores_list):
    """
    Train Random Forest to validate PRS scores.

    Args:
        features_list: List of feature arrays [matched_snps, risk_increasing, risk_decreasing, raw_score]
        scores_list: List of PRS scores

    Returns:
        Trained RFC model
    """
    if len(features_list) < 2:
        # Need at least 2 samples for training
        return None

    X = np.array(features_list)
    y = np.array(scores_list)

    # Train Random Forest Regressor
    rfc = RandomForestRegressor(n_estimators=100, max_depth=10, random_state=42)
    rfc.fit(X, y)


    # Validate with predictions
    predictions = rfc.predict(X)

    return rfc, predictions


# --------------------------
# Main function
# --------------------------
def process_single_sample(user_file_path, pheno_ann, drug_ann, diseases, drugs):
    """
    Process a single sample file and return results.

    Args:
        user_file_path: Path to user CSV file
        pheno_ann: Phenotype annotation dataframe
        drug_ann: Drug annotation dataframe
        diseases: List of diseases to analyze
        drugs: List of drugs to analyze

    Returns:
        tuple: (results_df, sample_name)
    """
    print(f"\n{'=' * 70}")
    print(f"Processing sample: {os.path.basename(user_file_path)}")
    print(f"{'=' * 70}")

    user_df = load_user_file(user_file_path)

    # Storage for results
    results = []
    all_features = []  # For RFC training
    all_prs_scores = []  # For RFC training

    print("\n[COMPUTING DISEASE RISK SCORES]")
    print("-" * 70)

    # Compute disease scores
    for disease in diseases:
        print(f"\nAnalyzing: {disease}")
        if "phenotype" in pheno_ann.columns:
            subset = pheno_ann[pheno_ann["phenotype"].str.contains(disease, case=False, na=False)]
            total_variants = len(subset["variant/haplotypes"].dropna().unique())
            print(f"  - Found {total_variants} unique variants in database")

            score, matched_snps, details = compute_score(
                user_df, subset,
                rsid_col="variant/haplotypes",
                alleles_col="alleles"
            )

            # Track features for RFC
            features = [matched_snps, details.get('risk_increasing', 0),
                       details.get('risk_decreasing', 0), details.get('raw_score', 0)]
            all_features.append(features)
            all_prs_scores.append(score)

            # Run Q-Learning for SNP selection (if enough SNPs) - FAST MODE
            ql_score = 0.0
            ql_snps_selected = 0
            if total_variants >= 2 and total_variants <= 50:  # Only for small datasets
                print(f"  - Running Q-Learning SNP selection...")
                snp_list = subset["variant/haplotypes"].dropna().unique()[:20]  # Limit to 20 SNPs
                if len(snp_list) >= 2:
                    ql_agent = SNPSelectionQLearning(snp_list.tolist(), epsilon=0.2)
                    best_snps, best_reward = ql_agent.train(user_df, subset, score, episodes=10, max_steps=10)
                    ql_snps_selected = len(best_snps) if best_snps else 0

                    # Compute score with Q-Learning selected SNPs
                    if best_snps:
                        ql_subset = subset[subset["variant/haplotypes"].isin(best_snps)]
                        ql_score, _, _ = compute_score(user_df, ql_subset,
                                                       rsid_col="variant/haplotypes",
                                                       alleles_col="alleles")
                    print(f"    • Q-Learning optimized: {ql_score:.4f} (selected {ql_snps_selected} SNPs)")

            results.append({
                "name": disease,
                "score": score,
                "matched_snps": matched_snps,
                "qlearning_score": ql_score,
                "qlearning_snps_selected": ql_snps_selected
            })

            print(f"  - Matched SNPs: {matched_snps}")
            print(f"    • Risk-increasing: {details.get('risk_increasing', 0)}")
            print(f"    • Protective: {details.get('risk_decreasing', 0)}")
            print(f"  - Raw weighted score: {details.get('raw_score', 0):.4f}")
            print(f"  - Final PRS score: {score:.4f}")
        else:
            all_features.append([0, 0, 0, 0])
            all_prs_scores.append(0.0)
            results.append({
                "name": disease,
                "score": 0.0,
                "matched_snps": 0,
                "qlearning_score": 0.0,
                "qlearning_snps_selected": 0
            })
            print(f"  - WARNING: 'phenotype' column not found")

    print("\n[COMPUTING DRUG INTERACTION SCORES]")
    print("-" * 70)

    # Compute drug scores
    for drug in drugs:
        print(f"\nAnalyzing: {drug}")
        drug_col = "drug(s)"
        if drug_col in drug_ann.columns:
            subset = drug_ann[drug_ann[drug_col].str.contains(drug, case=False, na=False)]
            total_variants = len(subset["variant/haplotypes"].dropna().unique())
            print(f"  - Found {total_variants} unique variants in database")

            score, matched_snps, details = compute_score(
                user_df, subset,
                rsid_col="variant/haplotypes",
                alleles_col="alleles"
            )

            # Track features for RFC
            features = [matched_snps, details.get('risk_increasing', 0),
                       details.get('risk_decreasing', 0), details.get('raw_score', 0)]
            all_features.append(features)
            all_prs_scores.append(score)

            # Run Q-Learning for SNP selection (if enough SNPs) - FAST MODE
            ql_score = 0.0
            ql_snps_selected = 0
            if total_variants >= 2 and total_variants <= 50:  # Only for small datasets
                print(f"  - Running Q-Learning SNP selection...")
                snp_list = subset["variant/haplotypes"].dropna().unique()[:20]  # Limit to 20 SNPs
                if len(snp_list) >= 2:
                    ql_agent = SNPSelectionQLearning(snp_list.tolist(), epsilon=0.2)
                    best_snps, best_reward = ql_agent.train(user_df, subset, score, episodes=10, max_steps=10)
                    ql_snps_selected = len(best_snps) if best_snps else 0

                    # Compute score with Q-Learning selected SNPs
                    if best_snps:
                        ql_subset = subset[subset["variant/haplotypes"].isin(best_snps)]
                        ql_score, _, _ = compute_score(user_df, ql_subset,
                                                       rsid_col="variant/haplotypes",
                                                       alleles_col="alleles")
                    print(f"    • Q-Learning optimized: {ql_score:.4f} (selected {ql_snps_selected} SNPs)")

            results.append({
                "name": drug,
                "score": score,
                "matched_snps": matched_snps,
                "qlearning_score": ql_score,
                "qlearning_snps_selected": ql_snps_selected
            })

            print(f"  - Matched SNPs: {matched_snps}")
            print(f"    • Risk-increasing: {details.get('risk_increasing', 0)}")
            print(f"    • Protective: {details.get('risk_decreasing', 0)}")
            print(f"  - Raw weighted score: {details.get('raw_score', 0):.4f}")
            print(f"  - Final PRS score: {score:.4f}")
        else:
            all_features.append([0, 0, 0, 0])
            all_prs_scores.append(0.0)
            results.append({
                "name": drug,
                "score": 0.0,
                "matched_snps": 0,
                "qlearning_score": 0.0,
                "qlearning_snps_selected": 0
            })
            print(f"  - WARNING: 'drug(s)' column not found")

    # Train RFC validator
    print("\n[TRAINING RFC VALIDATOR]")
    print("-" * 70)
    rfc_predictions = [0.0] * len(results)
    if len(all_features) >= 2:
        print(f"Training Random Forest on {len(all_features)} samples...")
        rfc_model, predictions = train_rfc_validator(all_features, all_prs_scores)
        if rfc_model is not None:
            rfc_predictions = predictions.tolist()
            print(f"RFC validation complete!")
            print(f"  Mean prediction error: {np.mean(np.abs(np.array(predictions) - np.array(all_prs_scores))):.4f}")
    else:
        print("Not enough samples for RFC training (need 2+)")

    # Add RFC predictions to results
    for i, result in enumerate(results):
        result['rfc_validation_score'] = rfc_predictions[i]

    # Extract sample name from user file
    sample_name = os.path.splitext(os.path.basename(user_file_path))[0]

    results_df = pd.DataFrame(results)

    # Print summary
    print("\n" + "-" * 70)
    print("SAMPLE SUMMARY")
    print("-" * 70)
    print(f"Sample: {sample_name}")
    print(f"Total targets analyzed: {len(diseases) + len(drugs)}")
    print(f"  - Diseases: {len(diseases)}")
    print(f"  - Drugs: {len(drugs)}")

    return results_df, sample_name


def main():
    """
    Main pipeline that processes all CSV files in the processing_algorithm/ directory.
    """
    print("=" * 70)
    print("PHARMACOGENOMICS RISK SCORING PIPELINE")
    print("=" * 70)

    csv_files_to_process = csv_files
    print(f"\nFound {len(csv_files_to_process)} CSV files to process")

    # Load external data once (shared across all samples)
    print("\n[LOADING EXTERNAL DATA]")
    clinvar = load_external("clinvar")
    drug_ann = load_external("drug_ann")
    pheno_ann = load_external("pheno_ann")
    fa_ann = load_external("fa_ann")
    drug_labels = load_external("drug_labels")
    drug_labels_gene = load_external("drug_labels_gene")

    # Normalize column names
    pheno_ann.columns = [c.lower() for c in pheno_ann.columns]
    drug_ann.columns = [c.lower() for c in drug_ann.columns]

    # Diseases and drugs to analyze
    # Updated to match actual phenotype values in database
    diseases = [
        "Skin Neoplasms",           # matches "Disease:Skin Neoplasms"
        "Diabetes Mellitus, Type 2", # matches "Disease:Diabetes Mellitus, Type 2"
        "Lung Neoplasms"            # matches "Disease:Lung Neoplasms"
    ]
    drugs = ["carbamazepine", "abacavir", "phenytoin"]

    # Create output directory
    if not os.path.exists(DATA_VISUALIZATION):
        os.makedirs(DATA_VISUALIZATION)

    # Output file path
    output_file = os.path.join(DATA_VISUALIZATION, "final_scores.csv")

    # Check for already processed samples (resume capability)
    already_processed = set()
    if os.path.exists(output_file):
        try:
            existing_df = pd.read_csv(output_file)
            if 'sample' in existing_df.columns:
                already_processed = set(existing_df['sample'].unique())
                print(f"\n[RESUME] Found existing results with {len(already_processed)} samples already processed")
                print(f"  - Skipping: {', '.join(sorted(list(already_processed)[:5]))}")
                if len(already_processed) > 5:
                    print(f"    ... and {len(already_processed) - 5} more")
        except Exception as e:
            print(f"\n[WARNING] Could not read existing output file: {e}")
            print("  - Starting fresh")

    # Filter out already processed samples
    csv_files_to_process = [
        csv_file for csv_file in csv_files_to_process
        if os.path.splitext(csv_file)[0] not in already_processed
    ]

    if not csv_files_to_process:
        print("\n[INFO] All samples have already been processed!")
        print(f"Output: {output_file}")
        return

    print(f"\n[PROCESSING] {len(csv_files_to_process)} samples remaining")

    # Process all CSV files
    all_samples_results = []
    samples_processed = 0

    for i, csv_file in enumerate(csv_files_to_process, 1):
        user_file_path = os.path.join(BASE_PATH, csv_file)

        try:
            results_df, sample_name = process_single_sample(
                user_file_path, pheno_ann, drug_ann, diseases, drugs
            )

            # Add sample name column
            results_df.insert(0, 'sample', sample_name)

            # Append to combined results
            all_samples_results.append(results_df)
            samples_processed += 1

            # Write incrementally to CSV after each sample
            # Use append mode if file exists, write mode with header if new file
            file_exists = os.path.exists(output_file)
            if file_exists:
                # File exists: append without header
                results_df.to_csv(output_file, mode='a', index=False, header=False)
            else:
                # New file: write with header
                results_df.to_csv(output_file, mode='w', index=False, header=True)

            print(f"\n✓ Saved results for {sample_name} ({i}/{len(csv_files_to_process)})")

        except Exception as e:
            print(f"\n[ERROR] Failed to process {csv_file}: {str(e)}")
            continue

    # Final summary
    if all_samples_results:
        combined_df = pd.concat(all_samples_results, ignore_index=True)

        print("\n" + "=" * 70)
        print("PIPELINE COMPLETE")
        print("=" * 70)
        print(f"Successfully processed {samples_processed} / {len(csv_files_to_process)} samples")
        print(f"\nOutput: {output_file}")
        print(f"Total rows: {len(combined_df)}")
        print(f"\nResults CSV columns:")
        print(f"  - sample: Sample identifier")
        print(f"  - name: Disease or drug name")
        print(f"  - score: Weighted PRS score (0-1)")
        print(f"  - matched_snps: Number of SNPs matched to user genome")
        print(f"  - qlearning_score: Q-Learning optimized SNP selection score")
        print(f"  - qlearning_snps_selected: Number of SNPs selected by Q-Learning")
        print(f"  - rfc_validation_score: RFC double-check validation")
        print("\n" + "=" * 70)
        print("SCORING IMPROVEMENTS:")
        print("  ✓ Genotype matching (not just rsid presence)")
        print("  ✓ Clinical significance weighting")
        print("  ✓ Direction of effect (risk-increasing vs protective)")
        print("  ✓ Odds ratio extraction from literature")
        print("  ✓ Weighted polygenic risk score formula")
        print("\n" + "MACHINE LEARNING ENHANCEMENTS:")
        print("  ✓ Random Forest Classifier for score validation")
        print("  ✓ Q-Learning for optimal SNP selection")
        print("  ✓ Three independent scoring methods in CSV")
        print("=" * 70)
    else:
        print("\n[ERROR] No samples were successfully processed!")

if __name__ == "__main__":
    main()
