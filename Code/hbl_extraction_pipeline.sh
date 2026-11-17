#!/bin/bash

################################################################################
#                 hbl GENE EXTRACTION PIPELINE - v3.2
#
# Extracts the hblA, hblB, hblC and hblD genes from B. cereus genome assemblies.
#
# Steps:
#   1. Combine isolate assemblies and build a BLAST nucleotide database
#   2. BLAST each hbl reference gene against the database
#   3. Filter hits, collapse overlapping hits within each gene, and resolve
#      hblA vs hblB at overlapping loci
#   4. Extract PLUS-strand hits and reverse-complement MINUS-strand hits
#   5. Write per-gene FASTA files and a strand summary
#
# Requirements: blastn, blastdbcmd, makeblastdb, Rscript (with tidyverse),
#               awk, grep, sed, tr, rev
#
# Inputs (run from the project folder):
#   Refs/       hblA.fasta, hblB.fasta, hblC.fasta, hblD.fasta
#   isolates/   one assembly per isolate, <isolate>_<suffix>.fasta
#               (the isolate name is the file name up to the first underscore)
#
# Usage:
#   bash snp_extraction_pipeline_v3.2.sh
#
# Outputs:
#   <gene>_gene_all.fasta        all extracted hits for the gene
#   <gene>_gene_plus.fasta       PLUS-strand hits
#   <gene>_gene_minus_RC.fasta   MINUS-strand hits, reverse-complemented
#   <gene>_filtered.tsv          filtered and collapsed hits
#   gene_copy_numbers.tsv        hbl gene hits per isolate (0 = not detected)
#   strand_summary.tsv           hit counts per gene and strand
#   contig_name_map.tsv         renamed contigs and their original headers
#   contig_lengths.tsv           contig lengths for the combined database
################################################################################

set -e

################################################################################
# CONFIGURATION
################################################################################

# Inputs are read from the working folder; every output goes into Figure_5_6/,
# with one subfolder per stage:
#   Figure_5_6/1_extraction/            BLAST, filtered hits, hbl gene FASTAs
#   Figure_5_6/2_protein_translation/   translated proteins
#   Figure_5_6/3_alignments/            MAFFT alignments
ROOT_DIR="$(pwd)"
OUT_BASE="$ROOT_DIR/Figure_5_6"
EXTRACT_DIR="$OUT_BASE/1_extraction"
mkdir -p "$EXTRACT_DIR"

REFS_DIR="$ROOT_DIR/Refs"
ISOLATES_DIR="$ROOT_DIR/isolates"

GENES=(hblA hblB hblC hblD)

# BLAST and hit-filtering thresholds
BLAST_PERC_ID=60        # minimum % identity for BLAST (blastn -perc_identity)
BLAST_QCOV=70           # minimum query coverage for BLAST (blastn -qcov_hsp_perc)
PARSE_MIN_PERC_ID=85    # keep hits with % identity greater than this value

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║           hbl GENE EXTRACTION PIPELINE - v3.2                  ║"
echo "║           hblA  hblB  hblC  hblD                               ║"
echo "║ BLAST >= ${BLAST_PERC_ID}% identity, >= ${BLAST_QCOV}% coverag ║"
echo "║           Hits kept at > ${PARSE_MIN_PERC_ID}% identity        ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

################################################################################
# STEP 0: CHECK DEPENDENCIES AND INPUTS
################################################################################

echo "════════════════════════════════════════════════════════════════"
echo "STEP 0: Checking dependencies and input files"
echo "════════════════════════════════════════════════════════════════"
echo ""

for cmd in blastn blastdbcmd makeblastdb Rscript awk grep sed tr rev; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd not found."
        exit 1
    fi
done

echo "✓ All dependencies found"

if [ ! -d "$REFS_DIR" ]; then
    echo "ERROR: $REFS_DIR directory not found"
    exit 1
fi

if [ ! -d "$ISOLATES_DIR" ]; then
    echo "ERROR: $ISOLATES_DIR directory not found"
    exit 1
fi

ref_count=$(find -L "$REFS_DIR" -maxdepth 1 -name "*.fasta" | wc -l)
iso_count=$(find -L "$ISOLATES_DIR" -maxdepth 1 -name "*.fasta" | wc -l)

echo "✓ Found $ref_count reference files in $REFS_DIR"
echo "✓ Found $iso_count isolate files in $ISOLATES_DIR"
echo ""

################################################################################
# STEP 1: COMBINE ISOLATES & BUILD BLAST DATABASE
################################################################################

echo "════════════════════════════════════════════════════════════════"
echo "STEP 1: Building BLAST database"
echo "════════════════════════════════════════════════════════════════"
echo ""

echo "Combining isolate FASTA files..."

cd "$EXTRACT_DIR" || exit 1
> combined_genomes.fasta
printf "contig\tisolate\tfile\toriginal_header\n" > contig_name_map.tsv

# Contigs are renamed <isolate>_<n> (e.g., PS03258_1) so each hit name splits
# cleanly into isolate and contig at the first underscore. Original contig
# headers are kept in contig_name_map.tsv.
seen_isolates=" "

for f in "$ISOLATES_DIR"/*.fasta; do
    [ -e "$f" ] || continue

    iso=$(basename "$f" .fasta)
    iso="${iso%%_*}"

    if [[ "$seen_isolates" == *" $iso "* ]]; then
        echo "ERROR: two isolate files give the same isolate name ($iso): $f"
        exit 1
    fi
    seen_isolates="${seen_isolates}${iso} "

    echo "  Adding $iso..."

    awk -v iso="$iso" -v file="$(basename "$f")" '
        /^>/ {
            n++
            print iso "_" n "\t" iso "\t" file "\t" substr($0, 2) >> "contig_name_map.tsv"
            print ">" iso "_" n
            next
        }
        { print }
    ' "$f" >> combined_genomes.fasta
done

total_seqs=$(grep -c "^>" combined_genomes.fasta || true)

if [ "$total_seqs" -eq 0 ]; then
    echo "ERROR: No sequences found in combined_genomes.fasta"
    exit 1
fi

echo "✓ Combined $total_seqs sequences"
echo ""

echo "Building BLAST database..."

# BLAST splits database paths on spaces, so the database is built in a temporary
# directory without spaces (for example, when the project is in OneDrive).
BLAST_DB_DIR=$(mktemp -d "${TMPDIR:-/tmp}/hbl_blastdb.XXXXXX")
BLAST_DB="$BLAST_DB_DIR/combined_db"
trap 'rm -rf "$BLAST_DB_DIR"' EXIT

if [[ "$BLAST_DB" == *" "* ]]; then
    echo "ERROR: temporary BLAST database path contains spaces: $BLAST_DB"
    echo "       Set TMPDIR to a directory without spaces and rerun."
    exit 1
fi

makeblastdb \
    -in combined_genomes.fasta \
    -dbtype nucl \
    -parse_seqids \
    -out "$BLAST_DB"

echo "✓ BLAST database ready"
echo ""

echo "Recording contig lengths..."

awk '
    /^>/ {
        if (name != "") print name "\t" len
        name = substr($1, 2)
        len = 0
        next
    }
    { len += length($0) }
    END {
        if (name != "") print name "\t" len
    }
' combined_genomes.fasta > contig_lengths.tsv

echo "✓ contig_lengths.tsv created"
echo ""

################################################################################
# STEP 2: RUN BLAST SEARCHES
################################################################################

echo "════════════════════════════════════════════════════════════════"
echo "STEP 2: Running BLAST searches"
echo "════════════════════════════════════════════════════════════════"
echo ""

run_blast() {

    local gene=$1
    local query="$REFS_DIR/${gene}.fasta"

    if [ ! -f "$query" ]; then
        echo "⚠ $query not found, skipping"
        return
    fi

    echo "Searching for $gene (${BLAST_PERC_ID}% identity, ${BLAST_QCOV}% coverage)..."

    blastn \
        -db "$BLAST_DB" \
        -query "$query" \
        -max_target_seqs 10000 \
        -perc_identity "$BLAST_PERC_ID" \
        -qcov_hsp_perc "$BLAST_QCOV" \
        -task megablast \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
        -out "${gene}_blast.tsv"

    hit_count=$(wc -l < "${gene}_blast.tsv" || true)

    echo "  ✓ ${gene}_blast.tsv created: ${hit_count} raw hits"

    plus_count=$(awk '$9 < $10' "${gene}_blast.tsv" | wc -l)
    minus_count=$(awk '$9 > $10' "${gene}_blast.tsv" | wc -l)

    echo "    PLUS-strand raw hits:  $plus_count"
    echo "    MINUS-strand raw hits: $minus_count"
    echo ""
}

for gene in "${GENES[@]}"; do
    run_blast "$gene"
done

################################################################################
# STEP 3: FILTER HITS & RESOLVE OVERLAPPING HITS
################################################################################

echo "════════════════════════════════════════════════════════════════"
echo "STEP 3: Filtering hits and resolving overlaps"
echo "════════════════════════════════════════════════════════════════"
echo ""

# Remove filtered tables from earlier runs so genes with no hits now are not
# extracted from stale files.
rm -f ./*_filtered.tsv

cat > filter_and_deduplicate.R << 'RSCRIPT'

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
MIN_PERC_ID <- as.numeric(args[1])

genes <- c("hblA", "hblB", "hblC", "hblD")

columns <- c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens",
             "q_start", "q_end", "s_start", "s_end", "eval", "score")

# Read BLAST output and keep hits with more than MIN_PERC_ID percent identity.
read_blast <- function(gene) {

  blast_file <- paste0(gene, "_blast.tsv")

  if (!file.exists(blast_file) || file.info(blast_file)$size == 0) {
    return(NULL)
  }

  read_tsv(blast_file, col_names = columns, comment = "#", show_col_types = FALSE) %>%
    separate(acc, c("acc", "contig"), "_") %>%
    drop_na() %>%
    filter(perc_id > MIN_PERC_ID) %>%
    mutate(
      full_name = paste0(acc, "_", contig),
      start     = ifelse(s_start < s_end, s_start, s_end),
      stop      = ifelse(s_start > s_end, s_start, s_end),
      gene      = gene
    )
}

# For one contig with hits sorted by start: when a hit overlaps the hit below it
# (stop >= next start), keep only the hit with the highest % identity from that
# run of overlapping hits. Non-overlapping hits are all kept.
collapse_overlaps <- function(df_test) {

  df_test <- data.frame(df_test)

  if (nrow(df_test) == 1) {
    return(df_test)
  }

  df_top_hits   <- data.frame()
  df_no_overlap <- data.frame()
  scores        <- data.frame()
  n <- 1

  while (n < nrow(df_test)) {

    # First hit of a run, or first hit after a run of overlapping hits
    if (nrow(scores) == 0) {
      scores <- df_test[n, ]
    }

    if (df_test$stop[n] - df_test$start[n + 1] >= 0) {
      # Overlap between hit n and hit n + 1
      scores <- rbind(scores, df_test[n + 1, ])

      if (n + 1 == nrow(df_test)) {
        df_top_hits <- rbind(df_top_hits, scores[which.max(scores$perc_id), ])
      }

    } else {
      # No overlap between hit n and hit n + 1
      if (nrow(scores) > 1) {
        df_top_hits <- rbind(df_top_hits, scores[which.max(scores$perc_id), ])
      } else {
        df_no_overlap <- rbind(df_no_overlap, df_test[n, ])
      }
      scores <- data.frame()

      if (n + 1 == nrow(df_test)) {
        df_no_overlap <- rbind(df_no_overlap, df_test[n + 1, ])
      }
    }

    n <- n + 1
  }

  rbind(df_no_overlap, df_top_hits)
}

# Group hits by isolate and contig, sort by start, and collapse overlaps
collapse_by_contig <- function(df) {

  df %>%
    group_by(full_name) %>%
    arrange(start, .by_group = TRUE) %>%
    group_split() %>%
    map(collapse_overlaps) %>%
    bind_rows()
}

# ---- Collapse overlapping hits within each gene ----

cat(paste0("Keeping hits with > ", MIN_PERC_ID, "% identity\n\n"))

filtered <- list()

for (gene_name in genes) {

  df <- read_blast(gene_name)

  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf("  %-6s no hits\n", gene_name))
    next
  }

  filtered[[gene_name]] <- collapse_by_contig(df)

  cat(sprintf("  %-6s %5d hits -> %4d after collapsing overlapping hits\n",
              gene_name, nrow(df), nrow(filtered[[gene_name]])))
}

# ---- Resolve hblA vs hblB: pool both genes' hits and keep the gene with the
#      higher % identity at each overlapping locus ----

if (!is.null(filtered$hblA) && !is.null(filtered$hblB)) {

  n_hblA <- nrow(filtered$hblA)
  n_hblB <- nrow(filtered$hblB)

  df_AB <- collapse_by_contig(bind_rows(filtered$hblA, filtered$hblB))

  filtered$hblA <- df_AB %>% filter(gene == "hblA")
  filtered$hblB <- df_AB %>% filter(gene == "hblB")

  cat(sprintf("\nhblA/hblB resolution: hblA %d -> %d, hblB %d -> %d\n",
              n_hblA, nrow(filtered$hblA), n_hblB, nrow(filtered$hblB)))
}

# ---- Save filtered results ----

cat("\nFiltered hits by gene:\n")

for (gene_name in names(filtered)) {

  df_out <- filtered[[gene_name]] %>%
    arrange(full_name, start) %>%
    mutate(
      iso     = acc,
      strand  = ifelse(s_start < s_end, "+", "-"),
      acc     = full_name,
      contig  = full_name
    ) %>%
    select(gene, acc, contig, start, stop, strand, perc_id, ali_len,
           mismatches, gap_opens, q_start, q_end, s_start, s_end, eval, score,
           iso, reference_hit = query)

  write_tsv(df_out, paste0(gene_name, "_filtered.tsv"))

  cat(sprintf("  %-8s %d hits in %d isolates (PLUS: %d | MINUS: %d) -> %s\n",
              paste0(gene_name, ":"),
              nrow(df_out),
              n_distinct(df_out$iso),
              sum(df_out$strand == "+"),
              sum(df_out$strand == "-"),
              paste0(gene_name, "_filtered.tsv")))
}

# ---- Gene copy number per isolate (0 = not detected) ----

isolates <- read_tsv("contig_name_map.tsv", show_col_types = FALSE) %>%
  distinct(isolate)

copy_numbers <- bind_rows(filtered) %>%
  count(isolate = acc, gene) %>%
  pivot_wider(names_from = gene, values_from = n)

copy_numbers <- isolates %>%
  left_join(copy_numbers, by = "isolate") %>%
  mutate(across(-isolate, ~ replace_na(.x, 0)))

write_tsv(copy_numbers, "gene_copy_numbers.tsv")

cat("\n✓ gene_copy_numbers.tsv created\n")
cat("✓ Filtering complete\n")

RSCRIPT

Rscript filter_and_deduplicate.R "$PARSE_MIN_PERC_ID"

echo ""

################################################################################
# STEP 4: EXTRACT SEQUENCES
################################################################################

echo "════════════════════════════════════════════════════════════════"
echo "STEP 4: Extracting sequences"
echo "════════════════════════════════════════════════════════════════"
echo ""

# Complement and reverse a sequence from stdin (lowercase is uppercased)
revcomp() {
    tr 'acgtn' 'ACGTN' | tr 'ACGT' 'TGCA' | rev
}

extract_sequences() {

    local gene=$1
    local in_file="${gene}_filtered.tsv"
    local plus_fa="${gene}_gene_plus.fasta"
    local minus_fa="${gene}_gene_minus_RC.fasta"
    local all_fa="${gene}_gene_all.fasta"

    if [ ! -f "$in_file" ]; then
        echo "⚠ $in_file not found, skipping"
        return
    fi

    echo "Extracting $gene..."

    : > "$plus_fa"
    : > "$minus_fa"

    # Columns in the filtered table: gene, acc, contig, start, stop, strand, ...
    # start and stop are already ordered (start < stop).
    tail -n +2 "$in_file" |
    while IFS=$'\t' read -r g acc contig start stop strand _rest
    do
        seq=$(blastdbcmd \
            -db "$BLAST_DB" \
            -entry "$contig" \
            -range "${start}-${stop}" \
            -outfmt "%s" \
            2>/dev/null) || seq=""

        if [ -z "$seq" ]; then
            echo "  ⚠ No sequence returned for ${contig}:${start}-${stop}"
            continue
        fi

        if [[ "$strand" == "+" ]]; then
            printf ">%s:%s-%s|%s|plus\n%s\n" \
                "$contig" "$start" "$stop" "$gene" "$seq" >> "$plus_fa"
        else
            printf ">%s:%s-%s|%s|minus_RC\n%s\n" \
                "$contig" "$start" "$stop" "$gene" "$(printf '%s' "$seq" | revcomp)" \
                >> "$minus_fa"
        fi
    done

    cat "$plus_fa" "$minus_fa" > "$all_fa"

    plus_count=$(grep -c "^>" "$plus_fa" || true)
    minus_count=$(grep -c "^>" "$minus_fa" || true)
    total_count=$(grep -c "^>" "$all_fa" || true)
    expected_count=$(tail -n +2 "$in_file" | wc -l | tr -d ' ')

    echo "  PLUS strand:          $plus_count"
    echo "  MINUS strand:         $minus_count"
    echo "  TOTAL extracted:      $total_count"

    if [ "$total_count" -ne "$expected_count" ] && [ "$expected_count" -gt 0 ]; then
        echo "  ⚠ WARNING: expected $expected_count filtered hits but extracted $total_count sequences"
    else
        echo "  ✓ All filtered hits extracted"
    fi

    echo ""
}

for gene in "${GENES[@]}"; do
    extract_sequences "$gene"
done

################################################################################
# STEP 5: FINAL REPORT
################################################################################

echo "════════════════════════════════════════════════════════════════"
echo "STEP 5: Final results"
echo "════════════════════════════════════════════════════════════════"
echo ""

echo "Summary by gene:"
echo ""

for gene in "${GENES[@]}"; do

    if [ -f "${gene}_gene_all.fasta" ]; then

        total=$(grep -c "^>" "${gene}_gene_all.fasta" || true)
        plus=$(grep -c "^>" "${gene}_gene_plus.fasta" || true)
        minus=$(grep -c "^>" "${gene}_gene_minus_RC.fasta" || true)

        printf "  %-8s Total: %4d | PLUS: %4d | MINUS: %4d\n" \
            "$gene" "$total" "$plus" "$minus"
    else
        printf "  %-8s No output\n" "$gene"
    fi

done

echo ""

echo "Creating strand summary..."

{
    printf "Gene\tTotal\tPlus\tMinus\n"

    for gene in "${GENES[@]}"; do

        if [ -f "${gene}_gene_all.fasta" ]; then

            total=$(grep -c "^>" "${gene}_gene_all.fasta" || true)
            plus=$(grep -c "^>" "${gene}_gene_plus.fasta" || true)
            minus=$(grep -c "^>" "${gene}_gene_minus_RC.fasta" || true)

            printf "%s\t%s\t%s\t%s\n" "$gene" "$total" "$plus" "$minus"
        fi

    done
} > strand_summary.tsv

echo "✓ strand_summary.tsv created"
echo ""

################################################################################
# OUTPUT SUMMARY
################################################################################

echo "════════════════════════════════════════════════════════════════"
echo "✓ PIPELINE COMPLETE (v3.2)"
echo "════════════════════════════════════════════════════════════════"
echo ""

echo "Generated files:"
for gene in "${GENES[@]}"; do
    if [ -f "${gene}_gene_all.fasta" ]; then
        echo "  ${gene}_gene_all.fasta"
    fi
done

echo ""
echo "  Summary file:"
echo "    strand_summary.tsv"
echo ""
echo "  Filtered results by gene (TSV format):"
for gene in "${GENES[@]}"; do
    if [ -f "${gene}_filtered.tsv" ]; then
        echo "    ${gene}_filtered.tsv"
    fi
done
echo ""

################################################################################
# STEPS 6-7: TRANSLATION AND ALIGNMENT (Python, run from this script)
################################################################################

echo "════════════════════════════════════════════════════════════════"
echo "STEP 6-7: Translating and aligning hbl proteins"
echo "════════════════════════════════════════════════════════════════"
cd "$OUT_BASE" || exit 1
python3 - <<'PY_EOF'
# -*- coding: utf-8 -*-
"""
Steps 6-7 of the hbl pipeline: translate the extracted sequences to protein,
then align the hbl protein sequences with MAFFT.

Step 6 (translation)
    Input : <gene>_gene_all.fasta in the working folder (from steps 1-5)
    Output: 2_protein_translation/<gene>_prot.fasta         protein sequences
            2_protein_translation/translation_summary.tsv   per-sequence QC

    The extracted sequences are already strand-corrected (minus-strand hits are
    reverse-complemented), so only the three forward frames are considered. For
    each sequence the frame with the fewest internal stop codons is chosen, using
    the bacterial genetic code (NCBI table 11). A trailing stop is removed;
    internal stops are kept and flagged.

Step 7 (alignment)
    Input : 2_protein_translation/<gene>_prot.fasta   (hblA-hblD only)
    Output: 3_alignments/<gene>_prot_align.fasta
            3_alignments/alignment_summary.tsv

    mafft --amino --op 5 --ep 1 is used. The gap penalties are non-default; the
    MAFFT defaults are --op 1.53 and --ep 0.123, which give more gapped alignments.
"""

import os
import subprocess

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---- Paths and settings (relative to the working folder) ----
# Paths are relative to Figure_5_6/ (the bash script changes into it before this block)
PIPELINE_DIR = "1_extraction"             # where steps 1-5 wrote <gene>_gene_all.fasta
PROT_DIR = "2_protein_translation"
ALIGN_DIR = "3_alignments"

TRANSLATE_GENES = ["hblA", "hblB", "hblC", "hblD", "nheA", "nheB", "nheC", "cytK1", "cytK2"]
ALIGN_GENES = ["hblA", "hblB", "hblC", "hblD"]

TABLE = 11                                # NCBI bacterial genetic code
GAP_ARGS = ["--op", "5", "--ep", "1"]


# ---------------------------------------------------------------------------
# Step 6: translation
# ---------------------------------------------------------------------------

def best_frame(nt):
    """Translation from the forward frame with the fewest internal stops."""
    best = None
    for f in range(3):
        trimmed = nt[f:]
        trimmed = trimmed[: len(trimmed) - len(trimmed) % 3]
        if not trimmed:
            continue
        aa = str(Seq(trimmed).translate(table=TABLE))
        internal = aa.rstrip("*").count("*")
        cand = (internal, f, aa)
        if best is None or cand[:2] < best[:2]:
            best = cand
    return best


def translate_all():
    os.makedirs(PROT_DIR, exist_ok=True)
    rows = []
    for gene in TRANSLATE_GENES:
        src = os.path.join(PIPELINE_DIR, f"{gene}_gene_all.fasta")
        if not os.path.exists(src) or os.path.getsize(src) == 0:
            print(f"{gene}: no input, skipped")
            continue

        out_recs = []
        for rec in SeqIO.parse(src, "fasta"):
            nt = str(rec.seq).upper()
            internal, frame, aa = best_frame(nt)
            aa_trimmed = aa.rstrip("*")
            out_recs.append(SeqRecord(Seq(aa_trimmed), id=rec.id, description=""))
            rows.append({
                "gene": gene, "sequence": rec.id,
                "nt_len": len(nt), "frame": frame, "aa_len": len(aa_trimmed),
                "internal_stops": internal,
                "starts_M": aa_trimmed.startswith("M"),
                "had_terminal_stop": aa.endswith("*"),
            })

        out = os.path.join(PROT_DIR, f"{gene}_prot.fasta")
        SeqIO.write(out_recs, out, "fasta")
        bad = sum(1 for r in rows if r["gene"] == gene and r["internal_stops"] > 0)
        lens = [r["aa_len"] for r in rows if r["gene"] == gene]
        print(f"{gene}: {len(out_recs)} proteins, {min(lens)}-{max(lens)} aa"
              f"{f', {bad} with internal stops' if bad else ''} -> {out}")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(PROT_DIR, "translation_summary.tsv"), sep="\t", index=False)

    print(f"\n{len(df)} sequences translated")
    print(f"frame 0 used for {(df.frame == 0).sum()} of {len(df)}")
    flagged = df[df.internal_stops > 0]
    if len(flagged):
        print(f"\n{len(flagged)} sequences with internal stop codons:")
        print(flagged[["gene", "sequence", "aa_len", "internal_stops"]].to_string(index=False))
    else:
        print("no internal stop codons in any sequence")


# ---------------------------------------------------------------------------
# Step 7: alignment
# ---------------------------------------------------------------------------

def run_mafft(in_file, out_file):
    cmd = ["mafft", "--quiet", "--anysymbol", "--amino"] + GAP_ARGS + [in_file]
    with open(out_file, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


def summarise(region, out_file):
    recs = list(SeqIO.parse(out_file, "fasta"))
    aln_len = len(recs[0].seq)
    ungapped = [len(str(r.seq).replace("-", "")) for r in recs]
    return {
        "region": region,
        "n": len(recs),
        "aln_len": aln_len,
        "ungapped_min": min(ungapped),
        "ungapped_max": max(ungapped),
        "gap_inflation": round(aln_len / max(ungapped), 2),
    }


def align_all():
    os.makedirs(ALIGN_DIR, exist_ok=True)
    rows = []
    for gene in ALIGN_GENES:
        src = os.path.join(PROT_DIR, f"{gene}_prot.fasta")
        if not os.path.exists(src) or os.path.getsize(src) == 0:
            print(f"{gene}: no protein input, skipped")
            continue

        out = os.path.join(ALIGN_DIR, f"{gene}_prot_align.fasta")
        run_mafft(src, out)
        row = summarise(f"{gene}_prot", out)
        rows.append(row)
        print(f"{gene}: {row['n']} seqs -> {out} ({row['aln_len']} columns)")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(ALIGN_DIR, "alignment_summary.tsv"), sep="\t", index=False)

    pd.set_option("display.width", 200)
    print("\n" + df.to_string(index=False))
    print("\ngap_inflation = alignment columns / longest ungapped sequence")
    print("wrote alignment_summary.tsv")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("== Step 6: translation ==")
    translate_all()
    print("\n== Step 7: alignment ==")
    align_all()
PY_EOF
