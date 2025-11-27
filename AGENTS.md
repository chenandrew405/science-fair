# Repository Guidelines

## Project Structure & Module Organization
Keep raw genomes inside `data/` together with `data.py`, which standardizes rsid/chromosome/position/genotype fields before exporting CSVs into `processing_algorithm/`. `processing_algorithm/processing.py` merges the cleaned genome with ClinVar/PharmGKB TSVs shipped in `external_data/` and owns the scoring logic. Derived artifacts, including `data_visualization/final_scores.csv`, should be treated as disposable build outputs; delete and regenerate them rather than editing by hand.

## Build, Test, and Development Commands
- `python -m venv .venv && source .venv/bin/activate` – optional but keeps pandas/numpy/scikit-learn isolated.
- `pip install pandas numpy scikit-learn` – minimal runtime stack; pin versions in `requirements.txt` if you add more libs.
- `python data/data.py` – converts every `.txt` genome in `data/` to CSV in `processing_algorithm/`, logging missing columns.
- `python processing_algorithm/processing.py` – runs the scoring pipeline and writes `data_visualization/final_scores.csv`; rerun whenever external TSVs change.

## Coding Style & Naming Conventions
Follow PEP 8 with 4-space indents, snake_case functions, and UPPER_SNAKE constants (see `DATA_FOLDER`, `EXTERNAL_FILES`). Favor descriptive DataFrame names (`drug_ann`, `disease_scores`) over single letters, and keep logging messages bracketed `[INFO]` / `[WARNING]` for grep-friendly output. Document non-obvious transformations (e.g., genotype normalization) inline.

## Testing Guidelines
Add pytest-based coverage under `tests/`, mirroring module names (`tests/test_processing.py`). Use lightweight fixtures created from small slices of the TSVs to keep runs fast. Minimum expectation is to unit-test `convert_txt_to_csv`, `load_external`, and `compute_score`, stubbing I/O with tmp paths. Run `pytest -q` before opening a PR and attach coverage deltas if you introduce new modules.

## Commit & Pull Request Guidelines
History shows short, imperative commits (“Add main.py”), so keep messages under 72 chars with a concise body when needed. Every PR should describe the change, list validation commands, and link the relevant issue or dataset ticket. Include screenshots or CSV diffs only if they illustrate an algorithmic shift; otherwise rely on reproducible commands.

## Security & Data Handling
Genomic files are personally identifiable; never commit raw `.txt` genomes or TSVs sourced from partners. Prefer `.gitignore` entries for temporary exports inside `data_visualization/` and scrub metadata from sample files before sharing. When debugging, truncate data to the first 50 rows and store under `data/samples/` so collaborators can reproduce issues without accessing sensitive records.
