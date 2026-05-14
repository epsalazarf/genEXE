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
# Dependencies: awk, sort (POSIX standard tools)

set -euo pipefail

REF_PVAR="${1:?Usage: pvar_merge_check.sh <ref.pvar> <query.pvar> [output.tsv]}"
QRY_PVAR="${2:?Usage: pvar_merge_check.sh <ref.pvar> <query.pvar> [output.tsv]}"
OUT_TSV="${3:-pvar_merge_check.tsv}"

[[ -f "$REF_PVAR" ]] || { echo "Error: reference pvar not found: $REF_PVAR" >&2; exit 1; }
[[ -f "$QRY_PVAR" ]] || { echo "Error: query pvar not found: $QRY_PVAR" >&2; exit 1; }

echo "Reference : $REF_PVAR" >&2
echo "Query     : $QRY_PVAR" >&2
echo "Output    : $OUT_TSV"  >&2

printf 'VARID\tCHROM\tPOS\tREF\tALT\tDS1\tDS2\tMATCH\n' > "$OUT_TSV"

awk '
BEGIN { OFS = "\t" }

# Skip VCF meta-information lines (##) and the column header (#CHROM ...)
/^##/ { next }
/^#/  { next }

# ── Pass 1: load reference pvar ──────────────────────────────────────────────
FNR == NR {
    key          = $1 SUBSEP $2
    ref1[key]    = $4
    alt1[key]    = $5
    order1[++n1] = key
    next
}

# ── Pass 2: load query pvar ───────────────────────────────────────────────────
{
    key = $1 SUBSEP $2
    if (!(key in seen2)) {
        seen2[key]    = 1
        ref2[key]     = $4
        alt2[key]     = $5
        order2[++n2]  = key
    }
}

END {
    c_perfect  = 0
    c_semi     = 0
    c_mismatch = 0
    c_ds1only  = 0
    c_ds2only  = 0

    # Emit reference variants (with match codes for shared sites)
    for (i = 1; i <= n1; i++) {
        key = order1[i]
        split(key, a, SUBSEP)
        chrom = a[1]; pos = a[2]
        r1    = ref1[key]
        al1   = alt1[key]
        vid   = chrom "-" pos "-" r1 "-" al1

        if (key in seen2) {
            if (r1 == ref2[key]) {
                if (al1 == alt2[key]) { m = 2;  c_perfect++  }
                else                   { m = 1;  c_semi++     }
            } else                     { m = -1; c_mismatch++ }
            print vid, chrom, pos, r1, al1, 1, 1, m
        } else {
            print vid, chrom, pos, r1, al1, 1, 0, 0
            c_ds1only++
        }
    }

    # Emit query-only variants
    for (j = 1; j <= n2; j++) {
        key = order2[j]
        if (!(key in ref1)) {
            split(key, a, SUBSEP)
            chrom = a[1]; pos = a[2]
            r2    = ref2[key]
            al2   = alt2[key]
            vid   = chrom "-" pos "-" r2 "-" al2
            print vid, chrom, pos, r2, al2, 0, 1, 0
            c_ds2only++
        }
    }

    total = n1 + c_ds2only

    print ""                                                      > "/dev/stderr"
    print "=== Merge Compatibility Summary ==="                   > "/dev/stderr"
    printf "%-34s %d\n", "Total unique variants:",  total        > "/dev/stderr"
    printf "%-34s %d\n", "DS1 (reference) variants:", n1        > "/dev/stderr"
    printf "%-34s %d\n", "DS2 (query) variants:", n2            > "/dev/stderr"
    print "---"                                                   > "/dev/stderr"
    printf "%-34s %d\n", "Perfect match     [2]:", c_perfect    > "/dev/stderr"
    printf "%-34s %d\n", "Semi-match        [1]:", c_semi       > "/dev/stderr"
    printf "%-34s %d\n", "REF mismatch     [-1]:", c_mismatch   > "/dev/stderr"
    printf "%-34s %d\n", "DS1 only (missing DS2) [0]:", c_ds1only > "/dev/stderr"
    printf "%-34s %d\n", "DS2 only (missing DS1) [0]:", c_ds2only > "/dev/stderr"
}
' "$REF_PVAR" "$QRY_PVAR" | sort -t$'\t' -k2,2V -k3,3n >> "$OUT_TSV"

echo "" >&2
echo "Done. Output written to: $OUT_TSV" >&2
