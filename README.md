# braker_cleaner

Clean BRAKER gene models by selecting one representative transcript per gene, validating sequence quality, integrating start/stop codon annotations into exon/CDS coordinates, and generating cleaned annotation/CDS/AA outputs.

## What This Tool Does

- Parses input gene model (`GFF3` or `GTF`) plus CDS/AA FASTA.
- Integrates codon features into exon/CDS structure when applicable.
- Applies sequence-quality filtering (length, start/stop codon checks, frame/internal stop checks).
- Selects one representative transcript per gene.
- Supports strict non-overlap representative enforcement across genes when `--overlap-threshold 0`.
- Writes cleaned annotation and sequence files.

## Requirements

- Python 3.8+
- Python packages:
  - `pyfaidx`
  - `intervaltree` (recommended)
  - `PyYAML` (optional, only if using YAML config files)

If `intervaltree` is missing, overlap logic falls back to a slower path.

## Installation

```bash
# from repository root
pip install -r requirements.txt  # if available
# or install minimal dependencies directly
pip install pyfaidx intervaltree pyyaml
```

## Input Files

Required:

- `--gene-model`: gene annotation in `GFF3` or `GTF`
- `--cds`: CDS FASTA (IDs must match transcript IDs)
- `--aa`: amino acid FASTA (IDs must match transcript IDs)
- `--genome`: genome FASTA (used for codon-aware sequence reconstruction)
- `--output-dir`: output directory

## Basic Usage

```bash
python pipeline_cli.py \
  --gene-model example_data/braker.gff3 \
  --cds example_data/braker.cds.fa \
  --aa example_data/braker.aa \
  --genome /path/to/genome.fa \
  --output-dir braker_cleaned
```

## Key Parameters

- `--min-length INT` (default: `50`)
  - Minimum amino acid length for valid transcripts.

- `--overlap-threshold FLOAT` (default: `0.5`, allowed range: `0..1`)
  - `0`: strict non-overlap mode. Output representatives are forced to be non-overlapping across genes.
  - `>0`: overlap is assessed/flagged by reciprocal overlap threshold.

- `--config PATH`
  - JSON or YAML config file.

- `--memory-limit INT` (MB)
- `--batch-size INT`
- `--log-level {DEBUG,INFO,WARNING,ERROR}`

## Strict Non-Overlap Mode

Use:

```bash
--overlap-threshold 0
```

Behavior:

- The pipeline tries alternative valid transcripts per gene to avoid genomic overlap with already selected representatives.
- If no valid non-overlapping candidate exists for a gene, that gene is excluded from cleaned output.
- Excluded genes appear in `manual_review_transcripts.txt`.

### Napkin Explanation: Which Gene Is Kept, and Why?

Think of strict mode as a greedy placement process:

1. Build a ranked candidate list for each gene (only quality-passing transcripts).
2. Rank genes globally by their best candidate quality.
3. Place genes one by one:
   - If a candidate does not overlap already-kept representatives, keep it.
   - If it overlaps, try the next candidate from that same gene.
   - If all candidates overlap, mark that gene as `unresolved_spatial_overlap`.

Ranking priority used by the code:

1. Larger AA-hash group size (more redundant support within the gene)
2. Longer AA length
3. Longer CDS length
4. Longer transcript genomic span
5. Lexicographic ID tie-break (`gene_id`, then `transcript_id`)

Practical consequence:

- If two genes overlap and both have only one valid transcript, one is kept and one is dropped.
- The kept one is simply the higher-ranked one by the rules above.
- The dropped one appears in `manual_review_transcripts.txt` with `unresolved_spatial_overlap`.

## Outputs

In `--output-dir`:

- `cleaned.gff3` (or `cleaned.gtf` if input was GTF)
  - Gene + one representative transcript + exon + CDS features.
- `cleaned.cds.fa`
- `cleaned.aa`
- `manual_review_transcripts.txt`
  - Genes with no selected representative.
- `processing_report.txt`

Notes:

- Output annotation does not emit explicit `start_codon` / `stop_codon` features.
- Sequence files may include quality flags in FASTA headers.

## Quality Flags

`manual_review_transcripts.txt` is driven by transcript `quality_flags`.
For each gene without a representative, the report contains the union of flags across all transcripts in that gene.

Common flags:

- `short_transcript`: protein length below `--min-length`.
- `low_quality_no_codons`: both start/stop codon annotations missing.
- `low_quality_no_start_codon`: AA sequence does not start with `M`.
- `low_quality_no_stop_codon`: CDS does not end with valid stop codon (`TAA/TAG/TGA`).
- `internal_stop_codons`: internal `*` detected in AA sequence.
- `cds_frame_error`: CDS length not divisible by 3.
- `empty_aa_sequence`: AA sequence is empty.
- `validated_cds`: CDS passed structural checks.
- `passed_aa_validation`: AA sequence passed checks.
- `representative`: transcript selected as representative.
- `post_selection_overlap`: overlap detected during non-strict overlap assessment (`--overlap-threshold > 0`).
- `unresolved_spatial_overlap`: transcript could not be placed in strict non-overlap mode (`--overlap-threshold 0`).
- `created_exon_for_cds`: exon was created to cover CDS.
- `redundant_start_codon_removed`, `redundant_stop_codon_removed`: codon already covered by existing exon/CDS.
- `merged_start_codon_with_adjacent`, `merged_stop_codon_with_adjacent`: codon merged into adjacent feature.
- `created_start_codon_features`, `created_stop_codon_features`: new exon/CDS created for codon.
- `cds_reconstructed_from_merged_coordinates`: CDS rebuilt from genome using merged coordinates.
- `cds_reconstruction_failed`: genome-based CDS reconstruction failed.

## Example (Bundled Test Inputs)

```bash
python pipeline_cli.py \
  --gene-model example_data/braker.gff3 \
  --cds example_data/braker.cds.fa \
  --aa example_data/braker.aa \
  --genome /path/to/genome.fa \
  --min-length 20 \
  --overlap-threshold 0 \
  --output-dir braker_cleaned
```

`example_data/` does not include a genome FASTA, so `--genome` must point to your own reference genome.

## Quick Validation Commands

Count gene features:

```bash
awk 'BEGIN{FS=OFS="\t"} $3=="gene" {n++} END{print n+0}' cleaned.gff3
```

Check strict pairwise overlaps (should be `0` in strict mode):

```bash
awk 'BEGIN{FS=OFS="\t"} $3=="gene" {
  id=""; if(match($9,/ID=[^;]+/)){id=substr($9,RSTART+3,RLENGTH-3)};
  print $1,$4-1,$5,id,".",$7
}' cleaned.gff3 \
| bedtools sort -i - \
| bedtools intersect -a - -b - -wa -wb \
| awk 'BEGIN{FS=OFS="\t"}
  !($1==$7 && $2==$8 && $3==$9 && $4==$10) {
    a=$1":"$2"-"$3":"$4; b=$7":"$8"-"$9":"$10;
    if (a<b) c++
  }
  END{print c+0}'
```

## Common Failure Causes

- Transcript IDs in `--cds` / `--aa` do not match annotation transcript IDs.
- Missing or unreadable genome FASTA for codon-aware processing.
- Invalid parameter values (for example overlap threshold outside `0..1`).

## Development

Run tests:

```bash
python -m unittest discover -s gene_curation_pipeline/tests -v
```
