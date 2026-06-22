#!/bin/bash
#SBATCH --job-name=adamix_cv_collect
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --export=ALL
#SBATCH --output=logs/adamix_cv_collect_%j.out
#SBATCH --error=logs/adamix_cv_collect_%j.err

# Title: adamixture-cv-collect.sh
# Description: Aggregation step for the ADAMIXTURE CV array. Scans the per-K CV
#              logs written by adamixture-cv-array.sh, extracts each K's CV error
#              into a single tab-separated table, and reports the K with the
#              lowest CV error (the elbow is best judged by eye from the table).
#              Designed to run as an afterany dependent of the array job, so it
#              still summarizes whatever tasks finished if some K failed.
# Usage:  sbatch --dependency=afterany:<array_job_id> \
#               bash/adamixture-cv-collect.sh <data_path> [name]
#         (normally launched for you by submit-adamixture-cv.sh)
# Developer: Pavel Salazar-Fernandez <pavel.salazar@galatea.bio>
# Dependencies: adamixture-cv-array.sh outputs
# Version: 0.1, 2026-06-22

set -euo pipefail

# ----------------------------------------------------------------------------
# Configuration (must match adamixture-cv-array.sh)
# ----------------------------------------------------------------------------
OUT_SUFFIX=_adamixture_cv   # CV results subfolder
CV_PATTERN='(cv|cross).*error'  # case-insensitive regex for the CV-error line
                                # (matches "CV error" and "Cross-validation error")

# ----------------------------------------------------------------------------
# Arguments — same as the array script, to locate the CV output folder
# ----------------------------------------------------------------------------
if [ $# -lt 1 ]; then
  echo "[ERROR] Missing input. Usage: sbatch $0 <data_path> [name]" >&2
  exit 1
fi

DATA_PATH=$(realpath "$1")
INPUT_DIR=$(dirname "${DATA_PATH}")
INPUT_FILE=$(basename "${DATA_PATH}")
case "${INPUT_FILE}" in
  *.vcf.gz) STEM="${INPUT_FILE%.vcf.gz}" ;;
  *.vcf)    STEM="${INPUT_FILE%.vcf}" ;;
  *.bed)    STEM="${INPUT_FILE%.bed}" ;;
  *.pgen)   STEM="${INPUT_FILE%.pgen}" ;;
  *)        STEM="${INPUT_FILE%.*}" ;;
esac
NAME=${2:-${STEM}}
OUTDIR="${INPUT_DIR}/${NAME}${OUT_SUFFIX}"

if [ ! -d "${OUTDIR}" ]; then
  echo "[ERROR] CV results folder not found: ${OUTDIR}" >&2
  exit 1
fi

# ----------------------------------------------------------------------------
# Collect per-K CV errors
# ----------------------------------------------------------------------------
shopt -s nullglob
LOGS=( "${OUTDIR}/${NAME}.cv."*.log )
if [ ${#LOGS[@]} -eq 0 ]; then
  echo "[ERROR] No per-K CV logs (${NAME}.cv.*.log) in ${OUTDIR}" >&2
  exit 1
fi

TSV="${OUTDIR}/${NAME}.cv_errors.tsv"
printf 'K\tCV_error\n' > "${TSV}"

# Build a K-sorted list from the log filenames, then extract each K's error
while read -r K; do
  LOG="${OUTDIR}/${NAME}.cv.${K}.log"
  # Last number on the last line matching the CV-error pattern (handles sci. notation)
  VAL=$(grep -iE "${CV_PATTERN}" "${LOG}" 2>/dev/null | tail -1 \
        | grep -oE '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?' | tail -1 || true)
  [ -z "${VAL}" ] && VAL=NA
  printf '%s\t%s\n' "${K}" "${VAL}" >> "${TSV}"
done < <(for f in "${LOGS[@]}"; do b=$(basename "${f}"); k=${b#"${NAME}.cv."}; echo "${k%.log}"; done | sort -n)

echo "[!] CV error table -> ${TSV}"
echo "---------------------------------------------"
column -t "${TSV}" 2>/dev/null || cat "${TSV}"
echo "---------------------------------------------"

# ----------------------------------------------------------------------------
# Report the best (lowest CV error) K; flag any K we could not parse
# ----------------------------------------------------------------------------
awk -F'\t' '
  NR>1 && $2=="NA" { missing = missing " " $1 }
  NR>1 && $2!="NA" { if (best=="" || $2+0 < best+0) { best=$2; bestk=$1 } }
  END {
    if (bestk != "") printf "[!] Lowest CV error: K=%s  (CV error = %s)\n", bestk, best
    else             print  "[!] WARNING: no CV errors could be parsed — check CV_PATTERN against the logs."
    if (missing != "") printf "[!] No CV error parsed for K:%s (task failed or different log wording)\n", missing
  }' "${TSV}"

echo "[!] Tip: the optimal K is the minimum OR where the curve flattens (elbow);"
echo "    eyeball the table/plot before committing to a K."
