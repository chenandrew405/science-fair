#!/usr/bin/env python3
"""
PGP Genetic Data Scraper with Arvados API Support

Downloads 23andMe genetic data from Personal Genome Project with two modes:
1. Traditional HTML scraping with authentication headers
2. Direct Arvados API access (recommended for bulk downloads)
"""

import os
import sys
import time
import argparse
import fnmatch
import signal
from pathlib import Path
from typing import Optional

try:
    import requests
    from bs4 import BeautifulSoup
except ImportError:
    print("Error: Required packages not installed. Run: pip install requests beautifulsoup4")
    sys.exit(1)


class TimeoutException(Exception):
    pass


def timeout_handler(signum, frame):
    raise TimeoutException("Operation timed out")


def download_with_html_scraping(
    url: str,
    output_dir: Path,
    arvados_token: Optional[str] = None,
    limit: Optional[int] = None,
    overwrite: bool = False,
    delay: float = 1.0,
    file_timeout: int = 300
) -> None:
    """Download files using HTML scraping with optional authentication."""
    from urllib.parse import urljoin

    print(f"Fetching PGP data from: {url}")

    headers = {}
    if arvados_token:
        headers['Authorization'] = f'Bearer {arvados_token}'
        print("Using Arvados authentication token")

    response = requests.get(url, headers=headers)
    response.raise_for_status()

    soup = BeautifulSoup(response.text, 'html.parser')

    # Find all /user_file/download/ links in the PGP table
    download_ids = []
    for link in soup.find_all('a', href=True):
        href = link['href']
        if '/user_file/download/' in href:
            # Extract just the ID or keep full relative URL
            full_url = urljoin('https://my.pgp-hms.org', href)
            download_ids.append(full_url)

    print(f"Found {len(download_ids)} download links on PGP portal")

    if limit:
        download_ids = download_ids[:limit]
        print(f"Limited to {limit} files")

    if not download_ids:
        print("No download links found. Check the page structure.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # Skip problematic collections that cause hangs
    skip_collections = ['/user_file/download/683', '/user_file/download/679', '/user_file/download/675']

    downloaded = 0
    for i, pgp_url in enumerate(download_ids, 1):
        # Skip known problematic collections
        if any(skip_id in pgp_url for skip_id in skip_collections):
            print(f"\n[{i}/{len(download_ids)}] Skipping problematic collection: {pgp_url}")
            continue

        print(f"\n[{i}/{len(download_ids)}] Fetching collection page: {pgp_url}")

        # Set a 60-second alarm for processing this entire collection
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(60)

        try:
            # Step 1: Get the Arvados collection landing page
            collection_response = requests.get(pgp_url, headers=headers, allow_redirects=True, timeout=30)
            collection_response.raise_for_status()

            # Step 2: Parse the collection page to find actual file links
            collection_soup = BeautifulSoup(collection_response.text, 'html.parser')

            # Find all file links in the collection (all genetic data files)
            file_links = []
            for link in collection_soup.find_all('a', href=True, class_='item'):
                href = link['href']
                # Download all files (txt, zip, vcf, etc.)
                file_url = urljoin(collection_response.url, href)
                filename = link.text.strip() or href.split('/')[-1]
                file_links.append((file_url, filename))

            if not file_links:
                print(f"  No files found in collection")
                continue

            # Step 3: Download each file in the collection
            for file_url, filename in file_links:
                output_path = output_dir / filename

                if output_path.exists() and not overwrite:
                    print(f"  Skipping {filename} (already exists)")
                    continue

                print(f"  Downloading {filename}...")

                try:
                    # Add timeout to prevent hanging on large files
                    file_response = requests.get(file_url, headers=headers, stream=True, timeout=file_timeout)
                    file_response.raise_for_status()

                    # Download with progress for large files
                    total_size = int(file_response.headers.get('content-length', 0))
                    downloaded_size = 0
                    start_time = time.time()

                    with open(output_path, 'wb') as f:
                        for chunk in file_response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                downloaded_size += len(chunk)
                                # Show progress for files > 50MB every 10MB
                                if total_size > 50_000_000 and downloaded_size % 10_000_000 < 8192:
                                    progress = (downloaded_size / total_size * 100) if total_size else 0
                                    elapsed = time.time() - start_time
                                    speed = downloaded_size / elapsed / 1_000_000  # MB/s
                                    print(f"    Progress: {downloaded_size:,} / {total_size:,} bytes ({progress:.1f}%) - {speed:.2f} MB/s")

                    size = output_path.stat().st_size
                    elapsed = time.time() - start_time
                    print(f"    ✓ Saved: {output_path} ({size:,} bytes in {elapsed:.1f}s)")
                    downloaded += 1

                except requests.exceptions.Timeout:
                    print(f"    ✗ Timeout downloading {filename} (>{file_timeout}s) - skipping")
                    if output_path.exists():
                        output_path.unlink()  # Remove partial download
                except Exception as e:
                    print(f"    ✗ Error downloading {filename}: {e}")
                    if output_path.exists():
                        output_path.unlink()  # Remove partial download

            if i < len(download_ids):
                time.sleep(delay)

            # Cancel the alarm if we completed successfully
            signal.alarm(0)

        except TimeoutException:
            signal.alarm(0)  # Cancel the alarm
            print(f"  ✗ Timeout processing collection (>60s) - skipping")
        except Exception as e:
            signal.alarm(0)  # Cancel the alarm
            print(f"  Error processing {pgp_url}: {e}")

    print(f"\n{'='*60}")
    print(f"Download complete! Total files downloaded: {downloaded}")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}")


def download_with_arvados_api(
    arvados_host: str,
    arvados_token: str,
    output_dir: Path,
    collection_uuid: Optional[str] = None,
    file_pattern: str = "*",
    limit: Optional[int] = None,
    overwrite: bool = False,
    delay: float = 1.0
) -> None:
    """Download files directly via Arvados API."""
    try:
        import arvados
    except ImportError:
        print("Error: Arvados API mode requires arvados-python-client")
        print("Install with: pip install arvados-python-client")
        sys.exit(1)

    if not arvados_token:
        print("Error: Arvados API mode requires an API token")
        print("Set ARVADOS_API_TOKEN environment variable or use --arvados-token")
        sys.exit(1)

    print(f"Connecting to Arvados at: {arvados_host}")

    # Initialize Arvados API client
    api = arvados.api(
        'v1',
        host=arvados_host,
        token=arvados_token,
        insecure=False
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    collections_to_process = []

    if collection_uuid:
        # Download from specific collection
        print(f"Downloading from collection: {collection_uuid}")
        collections_to_process.append(collection_uuid)
    else:
        # Search for collections with 23andMe data
        print("Searching for collections with 23andMe genetic data...")

        # Search collections by name containing "23andme"
        filters = [
            ['name', 'ilike', '%23andme%']
        ]

        collections = api.collections().list(
            filters=filters,
            limit=100
        ).execute()

        collection_items = collections.get('items', [])
        print(f"Found {len(collection_items)} collections matching '23andme'")

        if not collection_items:
            print("No collections found. Try specifying --collection-uuid directly.")
            return

        for coll in collection_items:
            collections_to_process.append(coll['uuid'])
            print(f"  - {coll['uuid']}: {coll.get('name', 'Unnamed')}")

    total_downloaded = 0

    for coll_uuid in collections_to_process:
        if limit and total_downloaded >= limit:
            print(f"\nReached download limit of {limit} files")
            break

        try:
            print(f"\nProcessing collection: {coll_uuid}")

            # Get collection
            collection = arvados.collection.CollectionReader(
                coll_uuid,
                api_client=api
            )

            # Find matching files
            matching_files = []
            for stream_path, _, file_list in collection.walk():
                for filename in file_list:
                    full_path = os.path.join(stream_path, filename).lstrip('./')
                    if fnmatch.fnmatch(filename, file_pattern):
                        matching_files.append(full_path)

            print(f"Found {len(matching_files)} files matching pattern '{file_pattern}'")

            # Download files
            for file_path in matching_files:
                if limit and total_downloaded >= limit:
                    break

                filename = os.path.basename(file_path)
                output_path = output_dir / filename

                if output_path.exists() and not overwrite:
                    print(f"  Skipping {filename} (already exists)")
                    continue

                print(f"  Downloading {filename}...")

                try:
                    with collection.open(file_path, 'rb') as src:
                        data = src.read()
                        output_path.write_bytes(data)

                    print(f"    Saved to: {output_path} ({len(data)} bytes)")
                    total_downloaded += 1

                    time.sleep(delay)

                except Exception as e:
                    print(f"    Error: {e}")

        except Exception as e:
            print(f"  Error processing collection {coll_uuid}: {e}")

    print(f"\nDownload complete! Total files downloaded: {total_downloaded}")


def main():
    parser = argparse.ArgumentParser(
        description='Download genetic data from Personal Genome Project',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # HTML scraping with authentication
  python %(prog)s --arvados-token "v2/..."

  # Direct Arvados API (recommended)
  python %(prog)s --use-arvados-api --arvados-host jutro.arvadosapi.com --limit 10

  # Download specific collection
  python %(prog)s --use-arvados-api --arvados-host jutro.arvadosapi.com \\
    --collection-uuid su92l-4zz18-xxxxxxxxxxxxx
        """
    )

    # Mode selection
    parser.add_argument(
        '--use-arvados-api',
        action='store_true',
        help='Use direct Arvados API instead of HTML scraping'
    )

    # Arvados options
    parser.add_argument(
        '--arvados-token',
        help='Arvados API token (or set ARVADOS_API_TOKEN env var)'
    )
    parser.add_argument(
        '--arvados-host',
        default='su92l.arvadosapi.com',
        help='Arvados API host (default: su92l.arvadosapi.com)'
    )
    parser.add_argument(
        '--collection-uuid',
        help='Specific collection UUID to download from'
    )
    parser.add_argument(
        '--file-pattern',
        default='*',
        help='File pattern to match for Arvados API mode (default: *)'
    )

    # HTML scraping options
    parser.add_argument(
        '--url',
        default='https://my.pgp-hms.org/public_genetic_data?utf8=%E2%9C%93&data_type=23andMe&commit=Search',
        help='PGP URL to scrape (for HTML mode)'
    )
    parser.add_argument(
        '--file-timeout',
        type=int,
        default=300,
        help='Timeout for downloading individual files in seconds (default: 300)'
    )

    # General options
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/pgp_downloads'),
        help='Output directory (default: data/pgp_downloads)'
    )
    parser.add_argument(
        '--limit',
        type=int,
        help='Maximum number of files to download'
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite existing files'
    )
    parser.add_argument(
        '--delay',
        type=float,
        default=1.0,
        help='Delay between downloads in seconds (default: 1.0)'
    )

    args = parser.parse_args()

    # Get API token from args or environment
    arvados_token = args.arvados_token or os.environ.get('ARVADOS_API_TOKEN')

    if args.use_arvados_api:
        download_with_arvados_api(
            arvados_host=args.arvados_host,
            arvados_token=arvados_token,
            output_dir=args.output_dir,
            collection_uuid=args.collection_uuid,
            file_pattern=args.file_pattern,
            limit=args.limit,
            overwrite=args.overwrite,
            delay=args.delay
        )
    else:
        download_with_html_scraping(
            url=args.url,
            output_dir=args.output_dir,
            arvados_token=arvados_token,
            limit=args.limit,
            overwrite=args.overwrite,
            delay=args.delay,
            file_timeout=args.file_timeout
        )


if __name__ == '__main__':
    main()
