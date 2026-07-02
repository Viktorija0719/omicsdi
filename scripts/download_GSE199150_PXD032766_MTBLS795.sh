#!/usr/bin/env bash
set -euo pipefail

# Download all data associated with:
#   GSE199150   transcriptomics (GEO)
#   PXD032766   proteomics/phosphoproteomics (PRIDE)
#   MTBLS795    metabolomics (MetaboLights)
# Project: Integrated multi-omics analysis of adverse cardiac remodeling and metabolic inflexibility upon ErbB2 and ERRα deficiency
#
# Usage:
#   bash download_GSE199150_PXD032766_MTBLS795.sh dry
#   bash download_GSE199150_PXD032766_MTBLS795.sh metadata
#   bash download_GSE199150_PXD032766_MTBLS795.sh processed
#   bash download_GSE199150_PXD032766_MTBLS795.sh all
#
# Modes:
#   dry        create metadata and manifests, print what would be downloaded
#   metadata   download OmicsDI JSON + GEO text record + manifests only
#   processed  download processed/support files, not huge raw .raw/.mzML files
#   all        download all files found in OmicsDI/GEO/PRIDE/MetaboLights manifests

MODE="${1:-dry}"
OUTDIR="${2:-downloads_ErbB2_ERRa_multiomics}"

if [[ ! "$MODE" =~ ^(dry|metadata|processed|all)$ ]]; then
  echo "ERROR: mode must be one of: dry, metadata, processed, all" >&2
  exit 1
fi

mkdir -p "$OUTDIR"/{omicsdi_json,manifests,geo_GSE199150,pride_PXD032766,metabolights_MTBLS795,logs}
LOG="$OUTDIR/logs/download_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command: $1" >&2; exit 1; }
}
need_cmd curl
need_cmd wget
need_cmd python3

run_cmd() {
  echo
  echo "\$ $*"
  if [[ "$MODE" != "dry" ]]; then
    "$@"
  fi
}

fetch_url() {
  local url="$1"
  local out="$2"
  echo
  echo "Fetching: $url"
  echo "Saving:   $out"
  if [[ "$MODE" != "dry" ]]; then
    curl -L --fail --retry 3 --connect-timeout 30 --max-time 300 "$url" -o "$out"
  fi
}

echo "Mode: $MODE"
echo "Output: $OUTDIR"
echo "Log: $LOG"

# 1) Save source records and OmicsDI normalized metadata.
# OmicsDI is a discovery layer; final files are downloaded from source archives.
fetch_url "https://www.omicsdi.org/ws/dataset/geo/GSE199150.json" \
  "$OUTDIR/omicsdi_json/geo_GSE199150.json"
fetch_url "https://www.omicsdi.org/ws/dataset/pride/PXD032766.json" \
  "$OUTDIR/omicsdi_json/pride_PXD032766.json"

# OmicsDI has used metabolights_dataset as the database slug for this record.
if [[ "$MODE" == "dry" ]]; then
  echo
  echo "Would try OmicsDI MetaboLights endpoint: metabolights_dataset/MTBLS795.json"
else
  if ! curl -L --fail --retry 3 --connect-timeout 30 --max-time 300 \
      "https://www.omicsdi.org/ws/dataset/metabolights_dataset/MTBLS795.json" \
      -o "$OUTDIR/omicsdi_json/metabolights_dataset_MTBLS795.json"; then
    echo "Primary MetaboLights OmicsDI endpoint failed; trying alternative slug metabolights."
    curl -L --fail --retry 3 --connect-timeout 30 --max-time 300 \
      "https://www.omicsdi.org/ws/dataset/metabolights/MTBLS795.json" \
      -o "$OUTDIR/omicsdi_json/metabolights_dataset_MTBLS795.json"
  fi
fi

fetch_url "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199150&targ=self&view=full&form=text" \
  "$OUTDIR/geo_GSE199150/GSE199150_full_record.soft.txt"

# 2) Build URL manifests from OmicsDI JSON.
# Converts EBI FTP URLs to HTTPS, and URL-encodes spaces in paths.
cat > "$OUTDIR/manifests/_extract_urls.py" <<'PY'
import json
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit, quote

outdir = Path(__import__('sys').argv[1])
json_dir = outdir / "omicsdi_json"
manifest_dir = outdir / "manifests"
manifest_dir.mkdir(parents=True, exist_ok=True)

RAW_SUFFIXES = (
    ".raw", ".mzml", ".mzxml", ".mgf", ".wiff", ".scan", ".d", ".fastq", ".fq", ".bam", ".cram"
)
PROCESSED_SUFFIXES = (
    ".txt", ".tsv", ".csv", ".xlsx", ".xls", ".json", ".xml", ".mzid", ".mztab", ".sdrf", ".pdf", ".html", ".zip", ".gz", ".tar"
)

def normalize_url(url: str) -> str:
    url = url.strip()
    if not url:
        return url
    # Prefer HTTPS over FTP where EBI/NCBI public FTP is mirrored over HTTPS.
    if url.startswith("ftp://ftp.pride.ebi.ac.uk/"):
        url = "https://ftp.pride.ebi.ac.uk/" + url[len("ftp://ftp.pride.ebi.ac.uk/"):]
    elif url.startswith("ftp://ftp.ebi.ac.uk/"):
        url = "https://ftp.ebi.ac.uk/" + url[len("ftp://ftp.ebi.ac.uk/"):]
    elif url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        url = "https://ftp.ncbi.nlm.nih.gov/" + url[len("ftp://ftp.ncbi.nlm.nih.gov/"):]

    parts = urlsplit(url)
    # Encode path spaces and special characters but keep slashes and existing percent escapes.
    path = quote(parts.path, safe="/%:@&=+$,;~*'()!-._")
    return urlunsplit((parts.scheme, parts.netloc, path, parts.query, parts.fragment))

def extract_urls(obj):
    urls = []
    if isinstance(obj, dict):
        for v in obj.values():
            urls.extend(extract_urls(v))
    elif isinstance(obj, list):
        for v in obj:
            urls.extend(extract_urls(v))
    elif isinstance(obj, str) and obj.startswith(("http://", "https://", "ftp://")):
        urls.append(normalize_url(obj))
    return urls

def write_list(path: Path, urls):
    urls = sorted(set(u for u in urls if u))
    path.write_text("\n".join(urls) + ("\n" if urls else ""), encoding="utf-8")
    print(f"{path}: {len(urls)} URL(s)")
    return urls

all_urls = []
for json_path in sorted(json_dir.glob("*.json")):
    try:
        obj = json.loads(json_path.read_text(encoding="utf-8"))
    except Exception as e:
        print(f"Could not read {json_path}: {e}")
        continue
    urls = extract_urls(obj)
    all_urls.extend(urls)
    write_list(manifest_dir / f"{json_path.stem}_all_urls.txt", urls)

# Main per-accession manifests.
pride = [u for u in all_urls if "PXD032766" in u or "pride/archive/projects/PXD032766" in u]
mtbls = [u for u in all_urls if "MTBLS795" in u]
geo = [u for u in all_urls if "GSE199150" in u]

# Add canonical GEO family directories explicitly.
geo_base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199150/"
geo += [
    geo_base,
    geo_base + "matrix/",
    geo_base + "miniml/",
    geo_base + "suppl/",
]

# Split processed/support from raw MS files.
def is_raw(u):
    low = u.lower()
    return low.endswith(RAW_SUFFIXES) or "/FILES/Method%20" in u and low.endswith(".mzml")

def is_processed(u):
    low = u.lower()
    return low.endswith(PROCESSED_SUFFIXES) and not is_raw(u)

write_list(manifest_dir / "GSE199150_geo_urls.txt", geo)
write_list(manifest_dir / "PXD032766_all_urls.txt", pride)
write_list(manifest_dir / "PXD032766_processed_urls.txt", [u for u in pride if is_processed(u)])
write_list(manifest_dir / "PXD032766_raw_urls.txt", [u for u in pride if is_raw(u)])
write_list(manifest_dir / "MTBLS795_all_urls.txt", mtbls)
write_list(manifest_dir / "MTBLS795_processed_urls.txt", [u for u in mtbls if is_processed(u)])
write_list(manifest_dir / "MTBLS795_raw_urls.txt", [u for u in mtbls if is_raw(u)])
write_list(manifest_dir / "ALL_project_urls.txt", all_urls)
PY

if [[ "$MODE" == "dry" ]]; then
  echo
  echo "Dry mode: metadata files are not downloaded, so URL manifests cannot be populated yet."
  echo "Run metadata first, then processed or all."
else
  python3 "$OUTDIR/manifests/_extract_urls.py" "$OUTDIR"
fi

if [[ "$MODE" == "metadata" || "$MODE" == "dry" ]]; then
  echo
  echo "Finished metadata/dry mode."
  exit 0
fi

# 3) Download GEO family directory.
# Keep folders because GEO has matrix/miniml/suppl subdirectories.
echo
if [[ "$MODE" == "processed" || "$MODE" == "all" ]]; then
  echo "Downloading GEO GSE199150 family folder: matrix, miniml, supplementary files."
  run_cmd wget -r -np -nH --cut-dirs=4 -R 'index.html*' --tries=5 --timeout=60 \
    -P "$OUTDIR/geo_GSE199150" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199150/"
fi

# 4) Download PRIDE PXD032766.
if [[ "$MODE" == "processed" ]]; then
  echo
  echo "Downloading PRIDE processed/support files only."
  run_cmd wget -c --tries=5 --timeout=60 --content-disposition \
    -i "$OUTDIR/manifests/PXD032766_processed_urls.txt" \
    -P "$OUTDIR/pride_PXD032766"
elif [[ "$MODE" == "all" ]]; then
  echo
  echo "ALL mode includes PRIDE raw .raw files. This can be large."
  read -r -p "Type YES to download all PXD032766 PRIDE files: " answer
  if [[ "$answer" == "YES" ]]; then
    run_cmd wget -c --tries=5 --timeout=60 --content-disposition \
      -i "$OUTDIR/manifests/PXD032766_all_urls.txt" \
      -P "$OUTDIR/pride_PXD032766"
  else
    echo "Skipping PRIDE all-file download."
  fi
fi

# 5) Download MetaboLights MTBLS795.
# Use -x/-nH/--cut-dirs=6 to preserve FILES/Method 1/... folder structure and avoid filename collisions.
if [[ "$MODE" == "processed" ]]; then
  echo
  echo "Downloading MetaboLights processed/support files only."
  run_cmd wget -c -x -nH --cut-dirs=6 --tries=5 --timeout=60 --content-disposition \
    -i "$OUTDIR/manifests/MTBLS795_processed_urls.txt" \
    -P "$OUTDIR/metabolights_MTBLS795"
elif [[ "$MODE" == "all" ]]; then
  echo
  echo "ALL mode includes MetaboLights raw .mzML files. This can be large."
  read -r -p "Type YES to download all MTBLS795 MetaboLights files: " answer
  if [[ "$answer" == "YES" ]]; then
    run_cmd wget -c -x -nH --cut-dirs=6 --tries=5 --timeout=60 --content-disposition \
      -i "$OUTDIR/manifests/MTBLS795_all_urls.txt" \
      -P "$OUTDIR/metabolights_MTBLS795"
  else
    echo "Skipping MetaboLights all-file download."
  fi
fi

echo
cat > "$OUTDIR/README_download_summary.txt" <<EOF
Project: Integrated multi-omics analysis of adverse cardiac remodeling and metabolic inflexibility upon ErbB2 and ERRa deficiency

Downloaded/attempted accessions:
  GSE199150  GEO transcriptomics
  PXD032766  PRIDE proteomics/phosphoproteomics
  MTBLS795   MetaboLights metabolomics

Important folders:
  omicsdi_json/                 normalized OmicsDI metadata
  manifests/                    URL lists extracted from OmicsDI JSON
  geo_GSE199150/                GEO record + matrix/miniml/supplementary files
  pride_PXD032766/              PRIDE files
  metabolights_MTBLS795/        MetaboLights files
  logs/                         command logs

To check files:
  find $OUTDIR -type f | sort > $OUTDIR/file_inventory.txt
  du -sh $OUTDIR/*

To verify PRIDE checksums if checksum.txt was downloaded:
  cd $OUTDIR/pride_PXD032766
  sha1sum -c checksum.txt 2>/dev/null || md5sum -c checksum.txt 2>/dev/null || true
EOF

find "$OUTDIR" -type f | sort > "$OUTDIR/file_inventory.txt"
du -sh "$OUTDIR"/* || true

echo
if [[ "$MODE" == "processed" ]]; then
  echo "Finished processed mode. Raw .raw/.mzML files were not fully downloaded."
  echo "Run: bash $0 all   when you are ready for complete raw-data download."
else
  echo "Finished all mode."
fi
