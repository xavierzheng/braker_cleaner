# Gene Model Cleaning Script

A Python script to clean BRAKER gene prediction output by selecting the best representative transcript for each gene, removing redundant alternative splice forms, and filtering out transcripts with insufficient amino acid length.

## Overview

This script processes BRAKER gene prediction output files to:
- Remove stop codon markers from protein sequences
- Filter out transcripts with amino acid sequences shorter than a specified threshold
- Select the longest transcript per gene as the representative model
- Handle equal-length transcripts by comparing MD5 hashes
- Report ambiguous cases requiring manual review
- Generate cleaned output files for downstream analysis

## Input Files

The script expects the following BRAKER output files in the input directory:

- `braker.aa` - Protein sequences in FASTA format
- `braker.cds.fa` - CDS sequences in FASTA format  
- `braker.gtf` - Gene annotations in GTF format
- `braker.gff3` - Gene annotations in GFF3 format

## Output Files

The script generates the following cleaned output files:

- `braker_cleaned.aa` - Cleaned protein sequences (selected transcripts only)
- `braker_cleaned.cds.fa` - Cleaned CDS sequences (selected transcripts only)
- `braker_cleaned.gtf` - Cleaned GTF annotations (selected transcripts only)
- `braker_cleaned.gff3` - Cleaned GFF3 annotations (selected transcripts only)
- `manual_review_transcripts.txt` - Transcripts requiring manual review
- `processing_report.txt` - Summary of processing results
- `gene_cleaning.log` - Detailed processing log

## Installation

### Requirements

- Python 3.6+
- Standard Python libraries (no external dependencies required)

### Setup

1. Clone or download the script to your working directory
2. Make the script executable (optional):
   ```bash
   chmod +x clean_gene_models.py
   ```

## Usage

### Basic Usage

```bash
# Run with default settings (minimum 10 amino acids)
python clean_gene_models.py

# Specify input directory
python clean_gene_models.py --input-dir /path/to/braker/output

# Set custom minimum amino acid length threshold
python clean_gene_models.py --min-aa-length 20

# Combine options
python clean_gene_models.py --input-dir /data/braker --min-aa-length 15
```

### Manual Curation Mode

After the initial run, you may need to manually review ambiguous cases:

1. Check `manual_review_transcripts.txt` for transcripts requiring manual review
2. Create a file with approved transcript IDs (one per line)
3. Run manual curation:
   ```bash
   python clean_gene_models.py --manual-curation approved_transcripts.txt
   ```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--input-dir` | Directory containing input files | Current directory (`.`) |
| `--min-aa-length` | Minimum amino acid length to keep transcript | 10 |
| `--manual-curation` | File containing approved transcript IDs | None |

## Processing Workflow

### Step 1: Clean Protein Sequences
- Remove trailing `*` characters (stop codon markers) from protein sequences
- Apply minimum amino acid length filter

### Step 2: Calculate Sequence Metrics
- Calculate protein sequence length
- Generate MD5 hashes for protein and CDS sequences
- Extract CDS positions from GTF/GFF3 annotations

### Step 3: Parse Gene/Transcript Identifiers
- Parse transcript names (format: `g1.t1` where `g1` = gene, `t1` = transcript)
- Group transcripts by gene ID for comparison

### Step 4: Select Longest Transcript per Gene
- For each gene, identify the transcript with the longest amino acid sequence
- If only one transcript exists, select it automatically

### Step 5: Handle Equal-Length Transcripts
- For transcripts with identical amino acid length:
  - Compare CDS MD5 hashes and positions
  - If identical: keep first transcript, remove others
  - If different: proceed to Step 6

### Step 6: Report Ambiguous Cases
- Group transcripts by amino acid MD5 hash
- Choose transcript variant with highest frequency
- If multiple variants have equal frequency: flag for manual review

### Step 7: Manual Curation (Optional)
- Process manually approved transcript IDs
- Generate final cleaned output files

## Output Format

### Processing Report

The `processing_report.txt` file contains:
- Processing statistics (genes, transcripts, filtered counts)
- List of transcripts filtered by length
- Gene-by-gene summary of transcript selection

### Manual Review File

The `manual_review_transcripts.txt` file contains:
- Gene IDs requiring manual review
- Associated transcript IDs
- Reason for manual review requirement

Format: `gene_id\ttranscript_ids\treason`

## Examples

### Example 1: Basic Cleaning with Default Parameters
```bash
python clean_gene_models.py
```

### Example 2: Stricter Filtering
```bash
python clean_gene_models.py --min-aa-length 30
```

### Example 3: Process Specific Directory
```bash
python clean_gene_models.py --input-dir /home/user/braker_results --min-aa-length 15
```

### Example 4: Manual Curation Workflow
```bash
# Initial run
python clean_gene_models.py --min-aa-length 20

# Review manual_review_transcripts.txt and create approved_list.txt
# Then run manual curation
python clean_gene_models.py --manual-curation approved_list.txt
```

## Logging

The script generates detailed logs in `gene_cleaning.log` including:
- Processing steps and progress
- Filtering statistics
- Warning messages for parsing issues
- Error messages for troubleshooting

## File Safety

**Important**: The script never modifies original input files. All output files use the `_cleaned` suffix to prevent accidental overwrites.

## Troubleshooting

### Common Issues

1. **Missing input files**: Ensure all required BRAKER output files are present
2. **Parsing errors**: Check transcript ID format (should be `gene.transcript`)
3. **Empty output**: Verify minimum amino acid length threshold isn't too restrictive
4. **Permission errors**: Ensure write permissions in the output directory

### Debug Mode

For detailed debugging information, check the log file:
```bash
tail -f gene_cleaning.log
```

## Performance

- **Memory usage**: Scales with number of transcripts (typically <1GB for most datasets)
- **Processing time**: ~1-5 minutes for typical BRAKER outputs (10K-100K transcripts)
- **Disk space**: Output files are typically 50-90% smaller than input files

## Citation

If you use this script in your research, please cite the associated publication or mention the script in your methods section.

## License

This script is provided as-is for research purposes. Please check with your institution regarding any licensing requirements.

## Support

For questions or issues:
1. Check the log file for error messages
2. Review the processing report for statistics
3. Verify input file formats match BRAKER output specifications

## Version History

- v1.0: Initial release with basic cleaning functionality
- v1.1: Added minimum amino acid length filtering
- v1.2: Enhanced reporting and manual curation features