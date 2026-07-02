#!/usr/bin/env python3
"""
Download selected cardio multi-omics datasets discovered via OmicsDI.

Input:
  - cardio_omicsdi_priority.csv  (your file is tab-separated even if extension is .csv)

Default behavior:
  - dry-run: prints what would be downloaded
  - downloads metadata/processed files only when --execute is used
  - raw FASTQ/proteomics/metabolomics directories are downloaded only with --mode all

Usage examples:
  python scripts/download_cardio_omicsdi_selected.py --input cardio_omicsdi_priority.csv
  python scripts/download_cardio_omicsdi_selected.py --input cardio_omicsdi_priority.csv --execute --mode processed
  python scripts/download_cardio_omicsdi_selected.py --input cardio_omicsdi_priority.csv --execute --mode all
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable, Any
from urllib.parse import quote_plus

try:
    import requests
except ImportError:
    requests = None


OMICS_DI_BASE = "https://www.omicsdi.org/ws"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# Extensions likely to be metadata or processed result files.
# Raw files can still sneak in if submitters zipped them as supplementary files.
PROCESSED_EXTENSIONS = (
    ".txt", ".tsv", ".csv", ".json", ".xml", ".xlsx", ".xls",
    ".pdf", ".html", ".htm", ".soft.gz", ".tar.gz", ".tgz",
    ".rdata", ".rds", ".mzid", ".mztab", ".sdrf"
)

RAW_EXTENSIONS = (
    ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".cram", ".sam",
    ".raw", ".wiff", ".wiff.scan", ".d", ".mzml", ".mzxml", ".mgf"
)


def eprint(*args: Any) -> None:
    print(*args, file=sys.stderr)


def read_selected_table(path: Path) -> list[dict[str, str]]:
    """
    Read selected dataset table. It auto-detects tab vs comma.
    Your current cardio_omicsdi_priority.csv is tab-separated.
    """
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    first_line = text.splitlines()[0]
    delimiter = "\t" if first_line.count("\t") >= first_line.count(",") else ","

    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row in reader:
            rows.append({k: (v or "").strip() for k, v in row.items()})
    return rows


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def geo_series_folder(gse: str) -> str:
    """
    GEO FTP series folders use the pattern:
      GSE108157 -> GSE108nnn
      GSE270093 -> GSE270nnn
    """
    m = re.match(r"^(GSE)(\d+)$", gse.upper())
    if not m:
        raise ValueError(f"Not a GEO Series accession: {gse}")
    prefix, digits = m.groups()
    return f"{prefix}{digits[:-3]}nnn"


def run_cmd(cmd: list[str], execute: bool, log_file: Path | None = None) -> int:
    """
    Print command. Run it only when execute=True.
    """
    printable = shlex.join(cmd)
    print(f"\n$ {printable}")

    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with log_file.open("a", encoding="utf-8") as fh:
            fh.write(printable + "\n")

    if not execute:
        return 0

    try:
        return subprocess.run(cmd, check=False).returncode
    except FileNotFoundError as exc:
        eprint(f"ERROR: command not found: {cmd[0]}")
        eprint("Install it, for example: conda install -c conda-forge wget curl")
        return 127


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def fetch_json(url: str, params: dict[str, str] | None = None, timeout: int = 120) -> Any:
    if requests is None:
        raise RuntimeError("requests is not installed. Run: pip install requests")
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.json()


def fetch_text(url: str, params: dict[str, str] | None = None, timeout: int = 120) -> str:
    if requests is None:
        raise RuntimeError("requests is not installed. Run: pip install requests")
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.text


def extract_urls(obj: Any) -> list[str]:
    urls: list[str] = []
    if isinstance(obj, dict):
        for v in obj.values():
            urls.extend(extract_urls(v))
    elif isinstance(obj, list):
        for v in obj:
            urls.extend(extract_urls(v))
    elif isinstance(obj, str):
        if obj.startswith(("http://", "https://", "ftp://")):
            urls.append(obj)
        # ENA often returns URLs without scheme, e.g. ftp.sra.ebi.ac.uk/...
        elif obj.startswith("ftp.sra.ebi.ac.uk/") or obj.startswith("ftp.ebi.ac.uk/"):
            urls.append("https://" + obj)
    return sorted(set(urls))


def parse_associated_accessions(value: str) -> dict[str, list[str]]:
    if not value:
        return {}
    try:
        parsed = json.loads(value)
        if isinstance(parsed, dict):
            return {str(k): list(v) if isinstance(v, list) else [str(v)] for k, v in parsed.items()}
    except Exception:
        pass
    return {}


def is_processed_like(url: str) -> bool:
    low = url.lower()
    return low.endswith(PROCESSED_EXTENSIONS)


def is_raw_like(url: str) -> bool:
    low = url.lower()
    return low.endswith(RAW_EXTENSIONS)


def write_url_list(path: Path, urls: Iterable[str]) -> None:
    urls = sorted(set(u for u in urls if u))
    write_text(path, "\n".join(urls) + ("\n" if urls else ""))


def download_url_list(url_file: Path, outdir: Path, execute: bool, log_file: Path) -> None:
    if not url_file.exists() or url_file.stat().st_size == 0:
        print(f"No URLs to download from {url_file}")
        return

    outdir.mkdir(parents=True, exist_ok=True)
    run_cmd(
        [
            "wget",
            "-c",
            "--tries=5",
            "--timeout=60",
            "--waitretry=10",
            "-i",
            str(url_file),
            "-P",
            str(outdir),
        ],
        execute=execute,
        log_file=log_file,
    )


def download_geo(accession: str, row: dict[str, str], outdir: Path, mode: str, execute: bool, log_file: Path) -> None:
    accession = accession.upper()
    dataset_dir = outdir / "geo" / accession
    dataset_dir.mkdir(parents=True, exist_ok=True)

    folder = geo_series_folder(accession)
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{folder}/{accession}"

    # Always save the source metadata page / SOFT text record.
    soft_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&view=full&form=text"
    run_cmd(
        ["curl", "-L", "--retry", "5", "--output", str(dataset_dir / f"{accession}_full_record.soft.txt"), soft_url],
        execute=execute,
        log_file=log_file,
    )

    # Useful GEO directories.
    # matrix = expression matrix when available
    # miniml = XML metadata
    # suppl = processed/supplementary files; can include raw tar archives, so check sizes.
    dirs = ["matrix", "miniml"]
    if mode in {"processed", "all"}:
        dirs.append("suppl")

    for d in dirs:
        url = f"{base}/{d}/"
        target = dataset_dir / d
        target.mkdir(parents=True, exist_ok=True)
        run_cmd(
            [
                "wget",
                "-r",
                "-np",
                "-nH",
                "--cut-dirs=5",
                "-R",
                "index.html*",
                "--tries=5",
                "--timeout=60",
                "-P",
                str(target),
                url,
            ],
            execute=execute,
            log_file=log_file,
        )


def process_ena_project(accession: str, outdir: Path, mode: str, execute: bool, log_file: Path) -> None:
    """
    Create ENA run report and optionally download FASTQ/submitted files.
    For large raw sequence data, use --mode all.
    """
    accession = accession.upper()
    dataset_dir = outdir / "ena" / accession
    dataset_dir.mkdir(parents=True, exist_ok=True)

    fields = ",".join([
        "run_accession",
        "sample_accession",
        "experiment_accession",
        "study_accession",
        "scientific_name",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
        "submitted_ftp",
        "submitted_md5",
        "submitted_bytes",
    ])

    report_url = (
        f"{ENA_FILEREPORT}?accession={quote_plus(accession)}"
        f"&result=read_run&fields={quote_plus(fields)}&format=tsv&download=true"
    )

    report_path = dataset_dir / f"{accession}_ena_run_report.tsv"

    if execute:
        print(f"\nFetching ENA run report: {accession}")
        try:
            text = fetch_text(report_url)
            write_text(report_path, text)
        except Exception as exc:
            eprint(f"WARNING: failed ENA filereport for {accession}: {exc}")
            return
    else:
        print(f"\nWould fetch ENA run report:\n{report_url}")
        print(f"Would save: {report_path}")
        return

    # Extract FASTQ/submitted file URLs.
    fastq_urls: list[str] = []
    submitted_urls: list[str] = []

    with report_path.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            for item in row.get("fastq_ftp", "").split(";"):
                item = item.strip()
                if item:
                    fastq_urls.append("https://" + item if item.startswith("ftp.") else item)
            for item in row.get("submitted_ftp", "").split(";"):
                item = item.strip()
                if item:
                    submitted_urls.append("https://" + item if item.startswith("ftp.") else item)

    write_url_list(dataset_dir / f"{accession}_fastq_urls.txt", fastq_urls)
    write_url_list(dataset_dir / f"{accession}_submitted_urls.txt", submitted_urls)

    if mode == "all":
        download_url_list(dataset_dir / f"{accession}_fastq_urls.txt", dataset_dir / "fastq", execute, log_file)
        download_url_list(dataset_dir / f"{accession}_submitted_urls.txt", dataset_dir / "submitted", execute, log_file)
    else:
        print(f"Saved ENA file reports and URL lists for {accession}.")
        print("FASTQ/submitted raw downloads are skipped unless you run with --mode all.")


def download_metabolights(accession: str, outdir: Path, mode: str, execute: bool, log_file: Path) -> None:
    accession = accession.upper()
    dataset_dir = outdir / "metabolights" / accession
    dataset_dir.mkdir(parents=True, exist_ok=True)

    base_url = f"https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{accession}/"
    page_url = f"https://www.ebi.ac.uk/metabolights/{accession}"

    write_text(dataset_dir / f"{accession}_download_notes.txt", f"Study page: {page_url}\nFTP directory: {base_url}\n")

    cmd = [
        "wget",
        "-r",
        "-np",
        "-nH",
        "--cut-dirs=6",
        "-R",
        "index.html*",
        "--tries=5",
        "--timeout=60",
        "-P",
        str(dataset_dir),
    ]

    if mode in {"metadata", "processed"}:
        cmd.extend([
            "-A",
            "*.txt,*.tsv,*.csv,*.json,*.xml,*.xlsx,*.xls,*.pdf,*.html,*.mztab,*.mzid,*.sdrf",
        ])

    cmd.append(base_url)
    run_cmd(cmd, execute=execute, log_file=log_file)

    if mode != "all":
        print(f"MetaboLights {accession}: metadata/processed mode only.")
        print("Run with --mode all if you also want every raw metabolomics/proteomics file.")


def fetch_pride_file_report(accession: str, outdir: Path, execute: bool) -> tuple[Path, list[str]]:
    accession = accession.upper()
    dataset_dir = outdir / "pride" / accession
    dataset_dir.mkdir(parents=True, exist_ok=True)

    report_json = dataset_dir / f"{accession}_pride_files.json"
    report_tsv = dataset_dir / f"{accession}_pride_file_urls.tsv"

    if not execute:
        print(f"\nWould fetch PRIDE file list for {accession}")
        print(f"Would try endpoint: https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files")
        return report_tsv, []

    all_items: list[Any] = []
    for page in range(0, 200):
        url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files"
        params = {"page": str(page), "size": "100"}
        try:
            data = fetch_json(url, params=params)
        except Exception as exc:
            eprint(f"WARNING: PRIDE API failed for {accession}: {exc}")
            break

        # PRIDE API may return HAL-style _embedded, or a simpler content/list shape.
        page_items: list[Any] = []
        if isinstance(data, dict):
            if "_embedded" in data:
                emb = data.get("_embedded") or {}
                for v in emb.values():
                    if isinstance(v, list):
                        page_items.extend(v)
            if "content" in data and isinstance(data["content"], list):
                page_items.extend(data["content"])
            if "files" in data and isinstance(data["files"], list):
                page_items.extend(data["files"])
            if not page_items and isinstance(data.get("list"), list):
                page_items.extend(data["list"])

        if not page_items:
            break

        all_items.extend(page_items)

        # Stop when API says we reached the last page.
        if isinstance(data, dict):
            page_info = data.get("page") or {}
            if isinstance(page_info, dict):
                total_pages = page_info.get("totalPages")
                if isinstance(total_pages, int) and page + 1 >= total_pages:
                    break

        time.sleep(0.2)

    report_json.write_text(json.dumps(all_items, indent=2), encoding="utf-8")

    rows: list[tuple[str, str, str]] = []
    for item in all_items:
        item_text = json.dumps(item)
        urls = extract_urls(item)
        file_name = ""
        file_type = ""

        if isinstance(item, dict):
            file_name = str(item.get("fileName") or item.get("name") or item.get("fileName", ""))
            file_type = str(item.get("fileType") or item.get("fileCategory") or "")

        # Sometimes the useful URL is in "downloadLink"; extract_urls catches it.
        for url in urls:
            rows.append((file_name, file_type, url))

    with report_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["file_name", "file_type", "url"])
        writer.writerows(rows)

    return report_tsv, [r[2] for r in rows]


def download_pride(accession: str, outdir: Path, mode: str, execute: bool, log_file: Path) -> None:
    accession = accession.upper()
    dataset_dir = outdir / "pride" / accession
    dataset_dir.mkdir(parents=True, exist_ok=True)

    report_tsv, urls = fetch_pride_file_report(accession, outdir, execute)

    if not execute:
        return

    if not urls:
        print(f"No PRIDE URLs found for {accession}; check {report_tsv} and the PRIDE web page.")
        return

    if mode == "all":
        selected = urls
    else:
        selected = [u for u in urls if is_processed_like(u) and not is_raw_like(u)]

    selected_path = dataset_dir / f"{accession}_{mode}_download_urls.txt"
    write_url_list(selected_path, selected)

    print(f"PRIDE {accession}: {len(urls)} total URLs, {len(selected)} selected for mode={mode}")
    download_url_list(selected_path, dataset_dir / mode, execute, log_file)


def download_massive(accession: str, outdir: Path, mode: str, execute: bool, log_file: Path) -> None:
    accession = accession.upper()
    dataset_dir = outdir / "massive" / accession
    dataset_dir.mkdir(parents=True, exist_ok=True)

    page_url = f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession={accession}"
    ftp_url = f"ftp://massive.ucsd.edu/v01/{accession}/"

    write_text(dataset_dir / f"{accession}_download_notes.txt", f"Study page: {page_url}\nLikely FTP directory: {ftp_url}\n")

    print(f"\nMassIVE {accession}")
    print("MassIVE proteomics folders can be very large.")
    print("This script will try the public FTP path. If it fails, open the study page and use MassIVE's file browser.")

    if mode == "metadata":
        run_cmd(
            ["curl", "-L", "--retry", "5", "--output", str(dataset_dir / f"{accession}_study_page.html"), page_url],
            execute=execute,
            log_file=log_file,
        )
        return

    cmd = [
        "wget",
        "-r",
        "-np",
        "-nH",
        "--cut-dirs=2",
        "-R",
        "index.html*",
        "--tries=5",
        "--timeout=60",
        "-P",
        str(dataset_dir),
    ]

    if mode == "processed":
        cmd.extend([
            "-A",
            "*.txt,*.tsv,*.csv,*.json,*.xml,*.xlsx,*.xls,*.pdf,*.html,*.mzid,*.mztab,*.sdrf",
        ])

    cmd.append(ftp_url)
    run_cmd(cmd, execute=execute, log_file=log_file)


def fetch_omicsdi_detail(source: str, accession: str, outdir: Path, execute: bool) -> None:
    """
    Save the OmicsDI normalized metadata JSON as provenance.
    OmicsDI source names sometimes differ from local source labels, so try a few forms.
    """
    candidates = []
    s = source.lower()
    candidates.append(s)

    source_map = {
        "geo": ["geo"],
        "project": ["project", "ena"],
        "metabolights_dataset": ["metabolights_dataset", "metabolights"],
        "pride": ["pride"],
        "massive": ["massive"],
    }
    candidates = source_map.get(s, candidates)

    dataset_dir = outdir / "_omicsdi_json"
    dataset_dir.mkdir(parents=True, exist_ok=True)

    for db in candidates:
        url = f"{OMICS_DI_BASE}/dataset/{db}/{accession}.json"
        path = dataset_dir / f"{safe_name(db + '_' + accession)}.json"

        if not execute:
            print(f"Would fetch OmicsDI JSON: {url}")
            return

        try:
            data = fetch_json(url)
            path.write_text(json.dumps(data, indent=2), encoding="utf-8")
            return
        except Exception:
            continue

    if execute:
        eprint(f"WARNING: could not fetch OmicsDI JSON for {source}:{accession}")


def process_row(row: dict[str, str], outdir: Path, mode: str, execute: bool, log_file: Path) -> None:
    source = row.get("source", "").strip().lower()
    accession = row.get("accession", "").strip()
    title = row.get("title", "").strip()

    if not accession:
        return

    print("\n" + "=" * 100)
    print(f"{source}:{accession}")
    print(title)
    print("=" * 100)

    fetch_omicsdi_detail(source, accession, outdir, execute)

    if source == "geo" or accession.upper().startswith("GSE"):
        download_geo(accession, row, outdir, mode, execute, log_file)

    elif source == "project" or accession.upper().startswith(("PRJ", "SRP", "ERP", "DRP")):
        process_ena_project(accession, outdir, mode, execute, log_file)

    elif source == "metabolights_dataset" or accession.upper().startswith("MTBLS"):
        download_metabolights(accession, outdir, mode, execute, log_file)

    elif source == "pride" or accession.upper().startswith("PXD"):
        download_pride(accession, outdir, mode, execute, log_file)

    elif source == "massive" or accession.upper().startswith("MSV"):
        download_massive(accession, outdir, mode, execute, log_file)

    else:
        print(f"Unknown source type for {source}:{accession}.")
        primary_url = row.get("primary_url", "")
        if primary_url:
            generic_dir = outdir / "other" / safe_name(f"{source}_{accession}")
            generic_dir.mkdir(parents=True, exist_ok=True)
            run_cmd(
                ["curl", "-L", "--retry", "5", "--output", str(generic_dir / "source_page.html"), primary_url],
                execute=execute,
                log_file=log_file,
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Download selected cardio OmicsDI datasets.")
    parser.add_argument("--input", default="cardio_omicsdi_priority.csv", help="Selected TSV/CSV file.")
    parser.add_argument("--outdir", default="downloads_selected", help="Output directory.")
    parser.add_argument(
        "--mode",
        choices=["metadata", "processed", "all"],
        default="processed",
        help=(
            "metadata = source pages, GEO matrix/miniml, API reports; "
            "processed = metadata + GEO supplementary + processed-like repository files; "
            "all = raw files too; can be very large."
        ),
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually download. Without this flag, only prints commands / dry-run.",
    )
    parser.add_argument(
        "--only",
        nargs="*",
        default=None,
        help="Optional accession subset, e.g. --only GSE108157 GSE109096 MTBLS795",
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    log_file = outdir / "_logs" / "download_commands.log"

    if not input_path.exists():
        raise SystemExit(f"Input file not found: {input_path}")

    if shutil.which("wget") is None:
        eprint("WARNING: wget not found. Install with: conda install -c conda-forge wget")
    if shutil.which("curl") is None:
        eprint("WARNING: curl not found. Install with: conda install -c conda-forge curl")

    rows = read_selected_table(input_path)

    if args.only:
        wanted = {x.upper() for x in args.only}
        rows = [r for r in rows if r.get("accession", "").upper() in wanted]

    outdir.mkdir(parents=True, exist_ok=True)

    # Save a copy of the selected manifest.
    manifest_path = outdir / "_manifest_selected_rows.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        if rows:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    print(f"Input rows: {len(rows)}")
    print(f"Output directory: {outdir}")
    print(f"Mode: {args.mode}")
    print(f"Execute: {args.execute}")
    print(f"Manifest copy: {manifest_path}")

    if not args.execute:
        print("\nDRY-RUN ONLY. Add --execute when the printed commands look correct.")

    for row in rows:
        process_row(row, outdir, args.mode, args.execute, log_file)

    print("\nFinished.")
    print(f"Command log: {log_file}")
    if args.mode != "all":
        print("\nRaw FASTQ/proteomics/metabolomics files were not fully downloaded.")
        print("Run again with --mode all --execute only after checking available disk space.")


if __name__ == "__main__":
    main()
