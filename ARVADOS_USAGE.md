# Arvados API Usage Guide

This guide explains how to use the enhanced PGP scraper with direct Arvados API authentication to download genetic files without HTML interference.

## Prerequisites

1. **Install dependencies**:
   ```bash
   source .venv/bin/activate
   pip install -r requirements.txt
   ```

2. **Get your Arvados API token**:
   - Log into your Arvados instance (e.g., https://workbench.su92l.arvadosapi.com)
   - Navigate to: User menu → Account Settings → Current token
   - Copy your API token (format: `v2/xxxxx-xxxxx-xxxxxxxxxxxxx/...`)

## Usage Modes

### Mode 1: Traditional HTML Scraping (with authentication)

Uses HTML scraping with added authentication headers:

```bash
# Using environment variable (recommended)
export ARVADOS_API_TOKEN="v2/xxxxx-xxxxx-xxxxxxxxxxxxx/..."
python data/pgp_scraper.py

# Or using command-line argument
python data/pgp_scraper.py --arvados-token "v2/xxxxx-xxxxx-xxxxxxxxxxxxx/..."
```

### Mode 2: Direct Arvados API (recommended)

Downloads files directly via Arvados API, bypassing HTML:

```bash
# Set your token and host
export ARVADOS_API_TOKEN="v2/xxxxx-xxxxx-xxxxxxxxxxxxx/..."

# Download from all PGP collections with 23andMe data
python data/pgp_scraper.py \
  --use-arvados-api \
  --arvados-host "su92l.arvadosapi.com" \
  --limit 10

# Download from a specific collection
python data/pgp_scraper.py \
  --use-arvados-api \
  --arvados-host "su92l.arvadosapi.com" \
  --collection-uuid "su92l-4zz18-xxxxxxxxxxxxx" \
  --file-pattern "*.txt"
```

## Command-Line Options

### Arvados-Specific Options

- `--use-arvados-api`: Enable direct Arvados API mode (bypasses HTML scraping)
- `--arvados-token TOKEN`: API token (or set `ARVADOS_API_TOKEN` env var)
- `--arvados-host HOST`: Arvados API host (e.g., `su92l.arvadosapi.com`)
- `--collection-uuid UUID`: Download from specific collection
- `--file-pattern PATTERN`: File pattern to match (default: `*.txt`)

### General Options

- `--output-dir DIR`: Output directory (default: `data/pgp_downloads`)
- `--limit N`: Maximum number of files to download
- `--overwrite`: Overwrite existing files
- `--delay SECONDS`: Delay between downloads (default: 1.0)

## Examples

### Example 1: Download specific collection

```bash
export ARVADOS_API_TOKEN="your-token-here"

python data/pgp_scraper.py \
  --use-arvados-api \
  --arvados-host "su92l.arvadosapi.com" \
  --collection-uuid "su92l-4zz18-123456789abcdef" \
  --output-dir "./genetic_data"
```

### Example 2: Search and download 23andMe files

```bash
export ARVADOS_API_TOKEN="your-token-here"

python data/pgp_scraper.py \
  --use-arvados-api \
  --arvados-host "su92l.arvadosapi.com" \
  --file-pattern "*.txt" \
  --limit 5 \
  --output-dir "./pgp_data"
```

### Example 3: Download all file types from a collection

```bash
export ARVADOS_API_TOKEN="your-token-here"

python data/pgp_scraper.py \
  --use-arvados-api \
  --arvados-host "su92l.arvadosapi.com" \
  --collection-uuid "su92l-4zz18-123456789abcdef" \
  --file-pattern "*"
```

## Finding Collection UUIDs

You can find collection UUIDs in several ways:

1. **Workbench UI**: Browse collections in the Arvados Workbench, copy UUID from URL
2. **API directly**: Use `arv` CLI tool or Python API to list collections
3. **From HTML links**: Look at download URLs in PGP portal (contains collection UUID)

Example collection UUID format: `su92l-4zz18-xxxxxxxxxxxxx`

## Troubleshooting

### "Arvados API mode requires arvados-python-client"
```bash
source .venv/bin/activate
pip install arvados-python-client
```

### "Arvados API mode requires an API token"
Set your token:
```bash
export ARVADOS_API_TOKEN="your-token-here"
```

### Authentication errors
- Verify token is correct and not expired
- Ensure `--arvados-host` matches your instance
- Check token has read permissions for collections

## Security Notes

- Never commit API tokens to git
- Use environment variables for tokens
- Tokens grant access to your Arvados account - keep them secure
- Consider using temporary tokens for automated downloads
