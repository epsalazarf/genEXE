#!/usr/bin/env bash
# =============================================================================
# check_skipped_variants.sh
# Diagnose the 522/1044 skipped entries warning from --score projection.
#
# HYPOTHESIS: Skipped variants are multiallelic in the ancient dataset —
# the A1 allele in ref_pca.eigenvec.allele does not match any ALT in the
# ancient pvar because the site has additional ALT alleles.
#
# USAGE
#   bash check_skipped_variants.sh [OUTDIR]
#   Default OUTDIR: ./pca_projection_out
#
# OUTPUT (to stdout)
#   - Count of unique IDs in eigenvec.allele
#   - Count of shared variants in intersect_variants.txt
#   - IDs in eigenvec.allele absent from intersect (pre-scoring exclusion)
#   - Shared variants that are multiallelic in the ancient pvar
#   - IDs with >1 row in eigenvec.allele (should be 0 for biallelic-only panels)
# =============================================================================
set -euo pipefail

OUTDIR="${1:-./pca_projection_out}"
EIGENVEC_ALLELE="${OUTDIR}/ref_pca.eigenvec.allele"
INTERSECT="${OUTDIR}/intersect_variants.txt"
ANC_PVAR="${OUTDIR}/ancient_filtered.pvar"

for f in "${EIGENVEC_ALLELE}" "${INTERSECT}" "${ANC_PVAR}"; do
    [[ -f "$f" ]] || { echo "ERROR: required file not found: $f"; exit 1; }
done

echo "=== Skipped variant diagnostic ==="
echo "OUTDIR: ${OUTDIR}"
echo ""

# --- A. Unique IDs in eigenvec.allele ---
awk 'NR>1 {print $2}' "${EIGENVEC_ALLELE}" | sort -u > /tmp/_eigenvec_ids.txt
echo "Unique variant IDs in eigenvec.allele:    $(wc -l < /tmp/_eigenvec_ids.txt)"

# --- B. Intersect list ---
sort "${INTERSECT}" > /tmp/_intersect_sorted.txt
echo "Variants in intersect_variants.txt:       $(wc -l < /tmp/_intersect_sorted.txt)"

# --- C. IDs in eigenvec.allele not in intersect (excluded before --score) ---
comm -23 /tmp/_eigenvec_ids.txt /tmp/_intersect_sorted.txt > /tmp/_missing_from_intersect.txt
echo "eigenvec IDs absent from intersect:       $(wc -l < /tmp/_missing_from_intersect.txt)"
if [[ $(wc -l < /tmp/_missing_from_intersect.txt) -gt 0 ]]; then
    echo "  (first 5):"
    head -5 /tmp/_missing_from_intersect.txt | sed 's/^/    /'
fi
echo ""

# --- D. Multiallelic sites in ancient pvar among shared variants ---
echo "=== Multiallelic check (comma in ALT field of ancient pvar) ==="
# ancient_filtered.pvar col 3=ID, col 5=ALT
awk 'NR>1 {print $3, $5}' "${ANC_PVAR}" | sort -k1,1 > /tmp/_anc_id_alt.txt

join -1 1 -2 1 /tmp/_intersect_sorted.txt /tmp/_anc_id_alt.txt | \
    awk '$2 ~ /,/ {print $1, $2}' > /tmp/_multiallelic_shared.txt

N_MULTI=$(wc -l < /tmp/_multiallelic_shared.txt)
echo "Shared variants multiallelic in ancient pvar: ${N_MULTI}"

if [[ "${N_MULTI}" -gt 0 ]]; then
    echo ""
    echo "First 20 (ID  ALT_field):"
    head -20 /tmp/_multiallelic_shared.txt | column -t | sed 's/^/  /'
fi
echo ""

# --- E. Row-per-ID check in eigenvec.allele ---
# NOTE: eigenvec.allele from --pca allele-wts writes exactly 2 rows per biallelic
# variant (one per allele: REF and ALT). >1 row per ID is therefore expected and
# correct. This check flags anomalies: IDs with exactly 1 row (truncated entry)
# or >2 rows (triallelic or malformed), both of which would indicate a problem.
echo "=== Row-per-ID check in eigenvec.allele ==="
echo "  (Expected: exactly 2 rows per ID for all biallelic SNPs)"
awk 'NR>1 {count[$2]++} END {for (id in count) print id, count[id]}' \
    "${EIGENVEC_ALLELE}" | sort -k2,2n > /tmp/_eigenvec_rowcounts.txt

N_ONE=$(awk '$2 == 1' /tmp/_eigenvec_rowcounts.txt | wc -l)
N_TWO=$(awk '$2 == 2' /tmp/_eigenvec_rowcounts.txt | wc -l)
N_MORE=$(awk '$2 > 2' /tmp/_eigenvec_rowcounts.txt | wc -l)
echo "  IDs with 1 row (anomalous — truncated?):  ${N_ONE}"
echo "  IDs with 2 rows (expected biallelic):      ${N_TWO}"
echo "  IDs with >2 rows (anomalous — triallelic?): ${N_MORE}"

if [[ "${N_ONE}" -gt 0 ]]; then
    echo "  WARNING: IDs with only 1 row (first 5):"
    awk '$2 == 1 {print $1}' /tmp/_eigenvec_rowcounts.txt | head -5 | sed 's/^/    /'
fi
if [[ "${N_MORE}" -gt 0 ]]; then
    echo "  WARNING: IDs with >2 rows (first 5):"
    awk '$2 > 2 {print $1, $2}' /tmp/_eigenvec_rowcounts.txt | head -5 | sed 's/^/    /'
fi
echo ""

# --- F. Accounting: reference-private variants explain the skipped entries ---
echo "=== Accounting ==="
N_ABSENT=$(wc -l < /tmp/_missing_from_intersect.txt)
EXPECTED_SKIPPED=$(( N_ABSENT * 2 ))
echo "Reference-private variants (in eigenvec.allele but not in intersect): ${N_ABSENT}"
echo "  Each has 2 rows in eigenvec.allele → ${EXPECTED_SKIPPED} expected skipped entries."
echo "  PLINK2 reported: 1044 skipped entries (from log warning)."
if [[ "${EXPECTED_SKIPPED}" -eq 1044 ]]; then
    echo "  MATCH: reference-private variants fully account for the warning."
    echo "  These are variants present in the reference panel but absent from the"
    echo "  ancient dataset. They are correctly ignored at scoring — no action needed."
else
    RESIDUAL=$(( 1044 - EXPECTED_SKIPPED ))
    echo "  PARTIAL match: ${RESIDUAL} skipped entries remain unexplained."
    echo "  Possible causes: allele code mismatches, REF/ALT swaps, or variant ID"
    echo "  format inconsistencies between ref and ancient pvar files."
fi

# --- Cleanup ---
rm -f /tmp/_eigenvec_ids.txt /tmp/_intersect_sorted.txt /tmp/_missing_from_intersect.txt \
       /tmp/_anc_id_alt.txt /tmp/_multiallelic_shared.txt /tmp/_eigenvec_rowcounts.txt

echo ""
echo "=== Done ==="
