#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_cardio_downloads.sh dry
#   bash run_cardio_downloads.sh processed
#   bash run_cardio_downloads.sh metadata
#   bash run_cardio_downloads.sh all
#
# Run this from your folder:
#   ~/RSU_work/omicsdi_cardio

INPUT="${INPUT:-cardio_omicsdi_priority.csv}"
OUTDIR="${OUTDIR:-downloads_selected}"
SCRIPT="${SCRIPT:-scripts/download_cardio_omicsdi_selected.py}"

MODE="${1:-dry}"

conda activate omicsdi 2>/dev/null || true

if [[ ! -f "$INPUT" ]]; then
  echo "Input file not found: $INPUT"
  echo "Current folder:"
  pwd
  echo "Files:"
  ls -lah
  exit 1
fi

if [[ ! -f "$SCRIPT" ]]; then
  echo "Script not found: $SCRIPT"
  echo "Copy download_cardio_omicsdi_selected.py to scripts/ first:"
  echo "  mkdir -p scripts"
  echo "  cp download_cardio_omicsdi_selected.py scripts/"
  exit 1
fi

case "$MODE" in
  dry)
    python "$SCRIPT" --input "$INPUT" --outdir "$OUTDIR" --mode processed
    ;;
  metadata)
    python "$SCRIPT" --input "$INPUT" --outdir "$OUTDIR" --mode metadata --execute
    ;;
  processed)
    python "$SCRIPT" --input "$INPUT" --outdir "$OUTDIR" --mode processed --execute
    ;;
  all)
    echo "WARNING: --mode all may download huge FASTQ/proteomics/metabolomics raw files."
    read -r -p "Continue? Type YES: " answer
    if [[ "$answer" != "YES" ]]; then
      echo "Cancelled."
      exit 0
    fi
    python "$SCRIPT" --input "$INPUT" --outdir "$OUTDIR" --mode all --execute
    ;;
  *)
    echo "Unknown mode: $MODE"
    echo "Use one of: dry, metadata, processed, all"
    exit 1
    ;;
esac
