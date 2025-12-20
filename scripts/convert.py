import os
import pandas as pd
import logging

PROJECT_ROOT = os.path.join("..", "data")
DATA_FOLDER = os.path.join(PROJECT_ROOT, "pgp_downloads")

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')


def convert_file_to_csv(file_path):
    """Convert a single data file (tab-separated) to .csv and export to processing_algorithm folder."""
    try:
        # Check if file is binary/compressed (common extensions)
        filename = os.path.basename(file_path).lower()
        compressed_patterns = {'.zip', '.gz', '.bz2', '.xz', '.tar', '.rar', '.7z'}

        # Check both single extension and multi-part extensions (e.g., .vcf.gz)
        if any(filename.endswith(pattern) for pattern in compressed_patterns):
            logging.info(f"Skipping compressed/binary file: {file_path}")
            return None

        # Read file and find the header line (starts with # but contains column names)
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        # Find header line (last comment line that looks like column names)
        header_line = None
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith('#'):
                # Check if this looks like a header (contains rsid, chromosome, etc.)
                if 'rsid' in line.lower():
                    header_line = line.lstrip('#').strip()
                    data_start = i + 1
            else:
                break

        if header_line:
            # Read with explicit header
            df = pd.read_csv(file_path, sep="\t", skiprows=data_start, names=header_line.split('\t'), dtype=str, low_memory=False)
        else:
            # Fallback to original method
            df = pd.read_csv(file_path, sep="\t", comment="#", dtype=str, low_memory=False)

        # Try to standardize required columns
        # We'll assume the most common ones: rsid, chromosome, position, genotype
        required_cols = ["rsid", "chromosome", "position", "genotype"]
        col_map = {}
        for col in required_cols:
            matches = [c for c in df.columns if c.lower() == col]
            if matches:
                col_map[matches[0]] = col
            else:
                logging.warning(f"Column '{col}' not found in {file_path}, will skip this column")

        df.rename(columns=col_map, inplace=True)

        # Clean genotype if exists
        if "genotype" in df.columns:
            df["genotype"] = df["genotype"].replace({"--": None})

        # Output CSV path - save to processing_algorithm folder
        csv_file = os.path.splitext(os.path.basename(file_path))[0] + ".csv"
        processing_folder = os.path.join(PROJECT_ROOT, "processing_algorithm")
        os.makedirs(processing_folder, exist_ok=True)
        csv_path = os.path.join(processing_folder, csv_file)
        df.to_csv(csv_path, index=False)
        logging.info(f"Converted {file_path} -> {csv_path}")
        return csv_path

    except Exception as e:
        logging.error(f"Failed converting {file_path}: {e}")
        return None


def convert_txt_to_csv(file_path):
    """Deprecated: Use convert_file_to_csv() instead. Kept for backwards compatibility."""
    logging.warning("convert_txt_to_csv() is deprecated, using convert_file_to_csv() instead")
    return convert_file_to_csv(file_path)


def convert_all_data_files():
    """Convert all data files (non-CSV, non-code) in pgp_downloads folder to .csv in processing_algorithm."""
    # File extensions to exclude (code, config, already converted)
    excluded_extensions = {'.csv', '.py', '.pyc', '.pyo', '.pyw', '.sh', '.bash',
                          '.json', '.yaml', '.yml', '.toml', '.md', '.gitignore', '.vcf', '.zip'}

    # Common non-data filenames to exclude (case-insensitive)
    excluded_names = {'license', 'readme', 'copying', 'changelog', 'makefile'}

    # Get all files in data folder
    all_files = []
    for f in os.listdir(DATA_FOLDER):
        file_path = os.path.join(DATA_FOLDER, f)

        # Skip directories
        if os.path.isdir(file_path):
            continue

        # Get file extension and base name
        _, ext = os.path.splitext(f)
        base_name = os.path.basename(f).lower()

        # Skip if extension is in excluded list
        if ext.lower() in excluded_extensions:
            continue

        # Skip if filename (without extension) is in excluded names
        name_without_ext = os.path.splitext(base_name)[0]
        if name_without_ext in excluded_names or base_name in excluded_names:
            continue

        all_files.append(f)

    if not all_files:
        logging.warning("No data files found to convert in pgp_downloads folder.")
        return []

    logging.info(f"Found {len(all_files)} file(s) to convert: {', '.join(all_files)}")

    csv_paths = []
    for data_file in all_files:
        file_path = os.path.join(DATA_FOLDER, data_file)
        csv_path = convert_file_to_csv(file_path)
        if csv_path:
            csv_paths.append(csv_path)
    return csv_paths


def convert_all_txt():
    """Deprecated: Use convert_all_data_files() instead. Kept for backwards compatibility."""
    logging.warning("convert_all_txt() is deprecated, using convert_all_data_files() instead")
    return convert_all_data_files()


if __name__ == "__main__":
    convert_all_data_files()
