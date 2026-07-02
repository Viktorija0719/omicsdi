# OmicsDI Multi-omics Dataset Discovery and Download Workflow

## 1. Environment setup

Create and activate a conda environment:

```bash
conda create -n omicsdi python=3.11 -y
conda activate omicsdi
```

Install required packages and command-line tools:

```bash
conda install -c conda-forge requests pandas openpyxl wget curl -y
pip install ddipy
```

On the HPC cluster, activate the environment with:

```bash
module load anaconda/conda-24.9.2
source activate omicsdi
```

Check that the environment is ready:

```bash
python --version
which python
which wget
which curl
```

## 2. About `ddipy`

[`ddipy`](https://omicsdi.readthedocs.io/en/latest/ddipy.html) is the Python client package for accessing OmicsDI through the OmicsDI REST web services.

OmicsDI is mainly a dataset discovery/indexing resource. It helps find datasets, metadata, related accessions, and file links, but full data downloads often still happen from the original repositories, such as GEO, ENA, PRIDE, MetaboLights, MassIVE, or other source archives.

Minimal `ddipy` example:

```python
from ddipy.dataset_client import DatasetClient

client = DatasetClient()

results = client.search(
    "cardiovascular multi-omics",
    "publication_date",
    "ascending"
)

print(results)
```

Get details for one dataset:

```python
from ddipy.dataset_client import DatasetClient

client = DatasetClient()

details = client.get_dataset_details("repository_name", "ACCESSION", False)
print(details)
```

Useful `DatasetClient` methods include:

| Method | Purpose |
|---|---|
| `search()` | Search OmicsDI datasets |
| `get_dataset_details()` | Retrieve normalized metadata for one dataset |
| `get_dataset_files()` | Retrieve file information for one dataset |
| `get_file_links()` | Retrieve downloadable file links when available |
| `get_similar()` | Retrieve similar or related datasets |
| `batch()` | Retrieve metadata for several datasets |
| `latest()` | Retrieve latest datasets |

## 3. General workflow

The workflow has three stages:

1. **Search OmicsDI** for candidate datasets.
2. **Inspect and filter candidates** to choose relevant multi-omics studies.
3. **Download metadata, processed files, and optionally raw files** from the original repositories.

Recommended order:

```bash
python scripts/search_cardio_omicsdi.py
bash scripts/run_cardio_downloads.sh dry
bash scripts/run_cardio_downloads.sh metadata
bash scripts/run_cardio_downloads.sh processed
```

Use full raw-data download only after checking disk space:

```bash
df -h .
bash scripts/run_cardio_downloads.sh all
```

## 4. How the scripts work

### `scripts/search_cardio_omicsdi.py`

This script searches OmicsDI using predefined cardiovascular and multi-omics search queries. It calls the OmicsDI REST API, retrieves dataset summaries, fetches detailed JSON metadata where possible, extracts URLs and related accessions, and writes output tables.

Main outputs:

```text
metadata/raw_search_results.json
metadata/cardio_omicsdi_enriched.json
metadata/candidate_datasets.tsv
links/all_extracted_urls.txt
```

The candidate table includes fields such as source repository, accession, title, publication date, likely multi-omics status, number of extracted URLs, and associated accessions.

### `scripts/download_cardio_omicsdi_selected.py`

This is the main flexible downloader for selected rows from a candidate table. It reads a TSV/CSV selection file, saves OmicsDI metadata as provenance, and downloads files according to repository type.

Supported source patterns include:

| Source type | Download behavior |
|---|---|
| GEO | Downloads SOFT record, matrix files, MINiML metadata, and supplementary files |
| ENA / BioProject | Creates ENA run reports and FASTQ/submitted-file URL lists; downloads raw files only in `all` mode |
| PRIDE | Queries PRIDE file metadata and downloads selected file URLs |
| MetaboLights | Downloads metadata/processed files or the full public study folder |
| MassIVE | Attempts public FTP download and saves source-page information |

Basic usage:

```bash
python scripts/download_cardio_omicsdi_selected.py \
  --input cardio_omicsdi_priority.csv \
  --outdir downloads_selected \
  --mode processed
```

To actually download files, add `--execute`:

```bash
python scripts/download_cardio_omicsdi_selected.py \
  --input cardio_omicsdi_priority.csv \
  --outdir downloads_selected \
  --mode processed \
  --execute
```

Download only selected accessions:

```bash
python scripts/download_cardio_omicsdi_selected.py \
  --input cardio_omicsdi_priority.csv \
  --outdir downloads_selected \
  --mode processed \
  --execute \
  --only ACCESSION1 ACCESSION2
```

Modes:

| Mode | Meaning |
|---|---|
| `metadata` | Download source pages, OmicsDI JSON, and metadata reports |
| `processed` | Download metadata plus processed/support files where available |
| `all` | Download raw files too; may require large storage |

### `scripts/run_cardio_downloads.sh`

This is a convenience wrapper around `download_cardio_omicsdi_selected.py`. It uses:

```text
INPUT=cardio_omicsdi_priority.csv
OUTDIR=downloads_selected
SCRIPT=scripts/download_cardio_omicsdi_selected.py
```

Run modes:

```bash
bash scripts/run_cardio_downloads.sh dry
bash scripts/run_cardio_downloads.sh metadata
bash scripts/run_cardio_downloads.sh processed
bash scripts/run_cardio_downloads.sh all
```

`dry` prints commands without downloading. `all` asks for confirmation because raw files can be very large.

### Project-specific combined downloader script

A project-specific shell script is also included for downloading one preselected multi-repository dataset combination. It has four modes:

```bash
bash scripts/<project_specific_downloader>.sh dry
bash scripts/<project_specific_downloader>.sh metadata
bash scripts/<project_specific_downloader>.sh processed
bash scripts/<project_specific_downloader>.sh all
```

Its logic is:

1. Download OmicsDI JSON metadata for each repository accession.
2. Download the source repository metadata page or record.
3. Extract URLs from OmicsDI JSON into manifest files.
4. Split manifests into processed/support files and raw-data files.
5. Download processed files first.
6. Download raw files only in `all` mode after user confirmation.

Typical output structure:

```text
downloads_<dataset_name>/
├── omicsdi_json/
├── manifests/
├── repository_1/
├── repository_2/
├── repository_3/
└── logs/
```

## 5. Recommended safety checks

Before downloading raw files:

```bash
df -h .
du -sh downloads_*/ 2>/dev/null
```

Create a file inventory:

```bash
find downloads_* -type f | sort > file_inventory.txt
```

Save folder sizes:

```bash
du -sh downloads_*/* > folder_sizes.txt
```

For PRIDE or other archives that provide checksum files, verify downloads when possible:

```bash
sha1sum -c checksum.txt 2>/dev/null || md5sum -c checksum.txt 2>/dev/null || true
```

## 6. Notes for reproducibility

Keep these files with every downloaded dataset:

- OmicsDI JSON metadata
- URL manifests
- source repository records or source pages
- download logs
- file inventory
- checksum files, when available

This makes it possible to reconstruct how the data were found, which source URLs were used, and which files were downloaded.
