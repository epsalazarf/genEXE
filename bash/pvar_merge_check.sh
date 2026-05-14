#!/usr/bin/env bash
# pvar_merge_check.sh — Compare two PLINK2 pvar files for merge compatibility.
#
# Usage:
#   pvar_merge_check.sh <ref.pvar> <query.pvar> [output.tsv]
#
# Arguments:
#   ref.pvar    Reference pvar file (REF allele priority for shared sites)
#   query.pvar  Query pvar file to compare against the reference
#   output.tsv  Output TSV path (default: pvar_merge_check.tsv)
#
# Output columns:
#   VARID   CHR-POS-REF-ALT identifier (reference pvar takes priority)
#   CHROM   Chromosome
#   POS     Position
#   REF     Reference allele (from ref.pvar when available)
#   ALT     Alternate allele (from ref.pvar when available)
#   DS1     Present in reference pvar (1/0)
#   DS2     Present in query pvar (1/0)
#   MATCH   Compatibility code:
#             2  Perfect match: same CHROM, POS, REF, ALT
#             1  Semi-match:    same CHROM, POS, REF; different ALT
#            -1  REF mismatch:  same CHROM, POS; different REF allele
#             0  Missing in one of the datasets
#
# Handles both plain pvars and pvars with VCF ##meta-information headers.
# Both files must have CHROM/POS/REF/ALT in columns 1/2/4/5 (standard PLINK2 format).
# Output is sorted by CHROM (version order) then POS (numeric).
#
# Memory: O(1) — uses sort+merge-join; suitable for files with millions of variants.
# Dependencies: awk, sort (GNU coreutils)

set -euo pipefail

REF_PVAR="${1:?Usage: pvar_merge_check.sh <ref.pvar> <query.pvar> [output.tsv]}"
QRY_PVAR="${2:?Usage: pvar_merge_check.sh <ref.pvar> <query.pvar> [output.tsv]}"
OUT_TSV="${3:-pvar_merge_check.tsv}"

[[ -f "$REF_PVAR" ]] || { echo "Error: reference pvar not found: $REF_PVAR" >&2; exit 1; }
[[ -f "$QRY_PVAR" ]] || { echo "Error: query pvar not found: $QRY_PVAR" >&2; exit 1; }

echo "Reference : $REF_PVAR" >&2
echo "Query     : $QRY_PVAR" >&2
echo "Output    : $OUT_TSV"  >&2

# Temp dir cleaned up on exit regardless of success/failure
WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

TMP1="$WORK_DIR/ds1.sorted"
TMP2="$WORK_DIR/ds2.sorted"

# ── Step 1: extract CHROM POS REF ALT from each pvar, sort to disk ───────────
# Sorting to disk keeps peak RAM low regardless of file size.
# A numeric sort key (chrom rank) is prepended so chromosomes sort correctly
# (1-22, X, Y, XY, MT) without relying on locale or version-sort availability.

pvar_extract_sort() {
    local file="$1"
    awk '
    /^##/ { next }
    /^#/  { next }
    {
        c = $1; sub(/^chr/, "", c)
        rank = (c+0 > 0) ? c+0 : (c=="X" ? 23 : c=="Y" ? 24 : c=="XY" ? 25 : \
                                   (c=="MT"||c=="M") ? 26 : 99)
        printf "%02d\t%d\t%s\t%s\t%s\t%s\n", rank, $2+0, $1, $2, $4, $5
    }' "$file" | sort -k1,1n -k2,2n
}

echo "Sorting DS1 (reference)..." >&2
pvar_extract_sort "$REF_PVAR" > "$TMP1"

echo "Sorting DS2 (query)..." >&2
pvar_extract_sort "$QRY_PVAR" > "$TMP2"

# ── Step 2: merge-join the two sorted streams, one line at a time ─────────────
# No arrays — only two lines live in memory at once.
# Output is already sorted; no downstream sort needed.

printf 'VARID\tCHROM\tPOS\tREF\tALT\tDS1\tDS2\tMATCH\n' > "$OUT_TSV"

awk -v f1="$TMP1" -v f2="$TMP2" '
BEGIN {
    OFS = "\t"
    c_perfect = c_semi = c_mismatch = c_ds1only = c_ds2only = 0

    r1 = (getline line1 < f1)
    r2 = (getline line2 < f2)

    while (r1 > 0 || r2 > 0) {

        if (r1 > 0) { split(line1, a, "\t"); rk1=a[1]+0; ps1=a[2]+0; ch1=a[3]; po1=a[4]; re1=a[5]; al1=a[6] }
        if (r2 > 0) { split(line2, a, "\t"); rk2=a[1]+0; ps2=a[2]+0; ch2=a[3]; po2=a[4]; re2=a[5]; al2=a[6] }

        # Determine relative order of the two current records
        if      (r1 <= 0)           cmp =  1   # DS1 exhausted
        else if (r2 <= 0)           cmp = -1   # DS2 exhausted
        else if (rk1 != rk2)        cmp = (rk1 < rk2) ? -1 : 1
        else if (ps1 != ps2)        cmp = (ps1 < ps2) ? -1 : 1
        else                        cmp =  0

        if (cmp < 0) {
            # Present only in DS1
            print ch1 "-" po1 "-" re1 "-" al1, ch1, po1, re1, al1, 1, 0, 0
            c_ds1only++
            r1 = (getline line1 < f1)

        } else if (cmp > 0) {
            # Present only in DS2
            print ch2 "-" po2 "-" re2 "-" al2, ch2, po2, re2, al2, 0, 1, 0
            c_ds2only++
            r2 = (getline line2 < f2)

        } else {
            # Same position — compare alleles; DS1 REF/ALT take priority
            if (re1 == re2) {
                if (al1 == al2) { m = 2; c_perfect++  }
                else             { m = 1; c_semi++     }
            } else               { m = -1; c_mismatch++ }
            print ch1 "-" po1 "-" re1 "-" al1, ch1, po1, re1, al1, 1, 1, m
            r1 = (getline line1 < f1)
            r2 = (getline line2 < f2)
        }
    }

    # Totals
    n_shared = c_perfect + c_semi + c_mismatch
    n_ds1    = n_shared + c_ds1only
    n_ds2    = n_shared + c_ds2only
    total    = n_shared + c_ds1only + c_ds2only

    print ""                                                        > "/dev/stderr"
    print "=== Merge Compatibility Summary ==="                     > "/dev/stderr"
    printf "%-34s %d\n", "Total unique variants:",   total         > "/dev/stderr"
    printf "%-34s %d\n", "DS1 (reference) variants:", n_ds1       > "/dev/stderr"
    printf "%-34s %d\n", "DS2 (query) variants:",     n_ds2       > "/dev/stderr"
    print "---"                                                     > "/dev/stderr"
    printf "%-34s %d\n", "Perfect match     [2]:", c_perfect      > "/dev/stderr"
    printf "%-34s %d\n", "Semi-match        [1]:", c_semi         > "/dev/stderr"
    printf "%-34s %d\n", "REF mismatch     [-1]:", c_mismatch     > "/dev/stderr"
    printf "%-34s %d\n", "DS1 only (missing DS2) [0]:", c_ds1only > "/dev/stderr"
    printf "%-34s %d\n", "DS2 only (missing DS1) [0]:", c_ds2only > "/dev/stderr"
}
' >> "$OUT_TSV"

echo "" >&2
echo "Done. Output written to: $OUT_TSV" >&2
