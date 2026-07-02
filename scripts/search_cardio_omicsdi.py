import json
import re
import time
from pathlib import Path
from urllib.parse import unquote

import requests


BASE = "https://www.omicsdi.org/ws"
OUTDIR = Path("metadata")
LINKDIR = Path("links")
OUTDIR.mkdir(exist_ok=True)
LINKDIR.mkdir(exist_ok=True)


CARDIO_QUERIES = [
    'cardiovascular AND "multi-omics"',
    '"coronary artery disease" AND "multi-omics"',
    '"coronary artery disease" AND transcriptomic AND proteomic',
    '"coronary artery disease" AND transcriptomic AND metabolomic',
    'atherosclerosis AND transcriptomic AND proteomic',
    'atherosclerosis AND transcriptomic AND metabolomic',
    '"myocardial infarction" AND transcriptomic AND proteomic',
    '"myocardial infarction" AND transcriptomic AND metabolomic',
    '"heart failure" AND transcriptomic AND proteomic',
    '"heart failure" AND transcriptomic AND metabolomic',
    'cardiomyopathy AND transcriptomic AND proteomic',
    'cardiomyopathy AND transcriptomic AND metabolomic',
    'cardiac AND methylation AND transcriptomic',
    'cardiac AND lipidomic AND transcriptomic',
    'vascular AND proteomic AND transcriptomic',
    'omics_type:"Transcriptomics" AND omics_type:"Proteomics" AND cardiovascular',
    'omics_type:"Transcriptomics" AND omics_type:"Metabolomics" AND cardiovascular',
    'omics_type:"Transcriptomics" AND omics_type:"Proteomics" AND TAXONOMY:"9606"',
    'omics_type:"Transcriptomics" AND omics_type:"Proteomics" AND TAXONOMY:"10090"',
]


ACCESSION_PATTERNS = {
    "GEO_series": r"\bGSE\d+\b",
    "GEO_sample": r"\bGSM\d+\b",
    "PRIDE_ProteomeXchange": r"\bPXD\d+\b",
    "ENA_BioProject": r"\bPRJ[EDN][A-Z]?\d+\b",
    "SRA_study": r"\bSRP\d+\b",
    "SRA_run": r"\bSRR\d+\b",
    "ArrayExpress": r"\bE-[A-Z]+-\d+\b",
    "MetaboLights": r"\bMTBLS\d+\b",
    "MassIVE": r"\bMSV\d+\b",
    "EGA_study": r"\bEGAS\d+\b",
    "EGA_dataset": r"\bEGAD\d+\b",
    "dbGaP": r"\bphs\d+\.v\d+\.p\d+\b",
    "PubMed": r"\bPMID:?\s*\d+\b",
}


def search_omicsdi(query, size=100, max_pages=3):
    """Search OmicsDI and return all dataset summaries."""
    all_datasets = []

    for page in range(max_pages):
        start = page * size
        r = requests.get(
            f"{BASE}/dataset/search",
            params={
                "query": query,
                "start": start,
                "size": size,
                "sort_field": "publication_date",
            },
            timeout=60,
        )
        r.raise_for_status()
        data = r.json()

        datasets = data.get("datasets", [])
        all_datasets.extend(datasets)

        count = int(data.get("count", 0))
        if start + size >= count or not datasets:
            break

        time.sleep(0.5)

    return all_datasets


def get_dataset_detail(source, accession):
    """Fetch normalized OmicsDI JSON for one dataset."""
    source = str(source).lower()
    accession = str(accession)

    url = f"{BASE}/dataset/{source}/{accession}.json"
    r = requests.get(url, timeout=60)

    if r.status_code != 200:
        return None

    return r.json()


def flatten_text(obj):
    """Convert nested JSON into searchable text."""
    if isinstance(obj, dict):
        return " ".join(flatten_text(v) for v in obj.values())
    if isinstance(obj, list):
        return " ".join(flatten_text(v) for v in obj)
    if obj is None:
        return ""
    return str(obj)


def extract_urls(obj):
    urls = []

    if isinstance(obj, dict):
        for v in obj.values():
            urls.extend(extract_urls(v))
    elif isinstance(obj, list):
        for v in obj:
            urls.extend(extract_urls(v))
    elif isinstance(obj, str):
        if obj.startswith(("http://", "https://", "ftp://")):
            urls.append(obj)

    return sorted(set(urls))


def extract_accessions(text):
    found = {}

    for label, pattern in ACCESSION_PATTERNS.items():
        matches = sorted(set(re.findall(pattern, text, flags=re.IGNORECASE)))
        if matches:
            found[label] = matches

    return found


def is_likely_multomics(text):
    text_l = text.lower()

    evidence_terms = [
        "multi-omics",
        "multiomics",
        "multi omics",
        "transcriptomic and proteomic",
        "proteomic and transcriptomic",
        "transcriptomic and metabolomic",
        "metabolomic and transcriptomic",
        "methylation and transcriptomic",
        "lipidomic and transcriptomic",
        "proteogenomic",
        "proteogenomics",
        "integrated transcriptomic",
        "integrated proteomic",
        "integrated metabolomic",
    ]

    omics_hits = 0
    for term in ["transcript", "proteom", "metabolom", "methyl", "genom", "lipidom"]:
        if term in text_l:
            omics_hits += 1

    keyword_hit = any(term in text_l for term in evidence_terms)

    return keyword_hit or omics_hits >= 2


def main():
    all_records = {}
    search_log = {}

    for query in CARDIO_QUERIES:
        print(f"\nSearching: {query}")
        datasets = search_omicsdi(query)
        search_log[query] = datasets
        print(f"  found summaries: {len(datasets)}")

        for d in datasets:
            accession = d.get("id") or d.get("accession")
            source = d.get("source") or d.get("database")

            if not accession or not source:
                continue

            key = f"{source}:{accession}"

            if key not in all_records:
                all_records[key] = {
                    "source": source,
                    "accession": accession,
                    "title": d.get("title") or d.get("name", ""),
                    "description": d.get("description", ""),
                    "publicationDate": d.get("publicationDate", ""),
                    "summary": d,
                    "queries": [],
                }

            all_records[key]["queries"].append(query)

    with open(OUTDIR / "raw_search_results.json", "w", encoding="utf-8") as f:
        json.dump(search_log, f, indent=2)

    print(f"\nUnique candidate datasets: {len(all_records)}")

    enriched = []

    for i, (key, rec) in enumerate(all_records.items(), start=1):
        source = rec["source"]
        accession = rec["accession"]

        print(f"[{i}/{len(all_records)}] Fetching detail: {source} {accession}")

        detail = get_dataset_detail(source, accession)
        text = flatten_text({"summary": rec["summary"], "detail": detail})
        urls = extract_urls(detail) if detail else []
        associated_accessions = extract_accessions(text)
        likely_multi = is_likely_multomics(text)

        rec["detail_available"] = detail is not None
        rec["likely_multiomics"] = likely_multi
        rec["associated_accessions"] = associated_accessions
        rec["urls"] = urls

        if detail:
            safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", key)
            with open(OUTDIR / f"{safe_name}.json", "w", encoding="utf-8") as f:
                json.dump(detail, f, indent=2)

        enriched.append(rec)
        time.sleep(0.3)

    with open(OUTDIR / "cardio_omicsdi_enriched.json", "w", encoding="utf-8") as f:
        json.dump(enriched, f, indent=2)

    tsv_path = OUTDIR / "candidate_datasets.tsv"
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write(
            "source\taccession\tlikely_multiomics\tpublicationDate\ttitle\t"
            "n_urls\tassociated_accessions\tqueries\n"
        )

        for rec in enriched:
            assoc = json.dumps(rec["associated_accessions"], ensure_ascii=False)
            queries = " || ".join(sorted(set(rec["queries"])))
            title = str(rec["title"]).replace("\t", " ").replace("\n", " ")
            f.write(
                f"{rec['source']}\t{rec['accession']}\t{rec['likely_multiomics']}\t"
                f"{rec['publicationDate']}\t{title}\t{len(rec['urls'])}\t"
                f"{assoc}\t{queries}\n"
            )

    links_path = LINKDIR / "all_extracted_urls.txt"
    with open(links_path, "w", encoding="utf-8") as f:
        for rec in enriched:
            for url in rec["urls"]:
                f.write(f"{rec['source']}\t{rec['accession']}\t{url}\n")

    print("\nDone.")
    print(f"Candidate table: {tsv_path}")
    print(f"Extracted URLs:  {links_path}")
    print("\nNext: inspect candidate_datasets.tsv and choose datasets to download.")


if __name__ == "__main__":
    main()
