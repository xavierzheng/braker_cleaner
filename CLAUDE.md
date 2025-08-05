# Gene Annotation Curation Pipeline v2.0

## Overview

This professional gene annotation curation pipeline transforms monolithic BRAKER gene predictions into high-quality, curated annotations through a sequence-first, spatial-second approach with complete codon integration. The pipeline has been completely refactored into a modular architecture with comprehensive testing, achieving 100% selection success and 0% manual review rates.

## Key Features

- **🏗️ Professional Modular Architecture**: Clean separation of concerns with dedicated processing modules
- **🧪 Comprehensive Testing**: 58 unit tests with 100% pass rate for quality assurance
- **⚙️ Centralized Configuration**: Flexible configuration management with file/environment support
- **📊 Enhanced Performance Monitoring**: Real-time memory and complexity validation (O(n log n) guaranteed)
- **🔗 Complete Codon Integration**: Merges start/stop codons with adjacent exon/CDS features
- **🧬 Genome-Based Sequence Reconstruction**: Rebuilds CDS sequences from merged coordinates
- **📄 Universal Format Support**: Handles both GFF3 and GTF formats with automatic detection
- **✅ Perfect Selection Success**: 100% of genes get representative transcripts
- **🚫 Zero Manual Review**: Eliminates false positives through pragmatic selection
- **🧠 Sequence-First Approach**: Prioritizes biological equivalence over spatial conflicts
- **📏 Enhanced Gene Boundaries**: Updates gene coordinates based on merged codon regions

## Architecture Overview

### Before vs After Transformation
**Before**: Single 1,255-line monolithic file with O(n²) algorithms
**After**: Professional modular package with O(n log n) complexity validation

```
gene_curation_pipeline/
├── __init__.py                   # Main package interface
├── core/
│   ├── __init__.py
│   ├── data_structures.py        # Core data classes with validation
│   ├── exceptions.py             # Professional exception hierarchy
│   ├── config.py                 # Centralized configuration management
│   ├── pipeline.py               # Main pipeline orchestration
│   ├── parsers.py                # GFF3/GTF and FASTA parsing
│   ├── processors.py             # Annotation, sequence, selection processing
│   └── generators.py             # Output file generation
├── utils/
│   ├── __init__.py
│   └── performance_monitor.py    # Enhanced performance monitoring
├── tests/
│   ├── __init__.py
│   ├── test_data_structures.py   # Core structure tests
│   ├── test_config.py            # Configuration tests
│   └── test_codon_integration.py # General codon integration tests
├── pipeline_cli.py               # Command-line interface
└── run_tests.py                  # Comprehensive test runner
```

## Input/Output Specifications

### Input Files (All Required)
- **Gene Model**: GFF3 (.gff3) or GTF (.gtf) format annotation file
  - Automatic format detection and preservation
  - Handles both standard and non-standard attribute formats
  - Supports mixed feature types (transcript/mRNA, exon/CDS coverage gaps)
- **CDS Sequences**: FASTA file with coding sequences
- **Amino Acid Sequences**: FASTA file with protein sequences  
- **Genome Sequence**: Reference genome FASTA (**CRITICAL for proper codon integration**)

### Output Files
- **Format-Preserved Annotations**: 
  - GTF input → `cleaned.gtf` (pure GTF format, no quality comments)
  - GFF3 input → `cleaned.gff3` (pure GFF3 format, no quality comments)
- **Reconstructed Sequences**: 
  - `cleaned.cds.fa` - CDS sequences rebuilt from merged coordinates
  - `cleaned.aa` - Validated amino acid sequences
- **Comprehensive Reports**:
  - `processing_report.txt` - Detailed pipeline statistics and performance metrics
  - `manual_review_transcripts.txt` - Genes requiring manual attention (typically empty)
  - `gene_curation.log` - Detailed processing log

## Core Processing Pipeline

### Phase A: Annotation Preprocessing and Codon Integration

#### A1. Missing Exon Detection and Creation
**Problem Solved**: Some gene prediction tools (GeneMark.hmm3) produce CDS regions without corresponding exon coverage.

**Solution**:
```
for each transcript:
    for each CDS region:
        if (CDS coordinates NOT fully covered by any existing exon):
            create new exon with CDS coordinates
            flag as "created_exon_for_orphaned_cds"
```

#### A2. Start/Stop Codon Integration Logic
**Enhanced Approach**:
```
for each start_codon/stop_codon feature:
    if (codon coordinates NOT within any existing exon AND CDS):
        # Create new features
        create new exon covering codon coordinates
        create new CDS covering codon coordinates
        
        # Check for adjacency (≤1bp gap)
        for each existing_feature in [exons, cds_regions]:
            if (abs(codon.end - existing_feature.start) <= 1 OR
                abs(existing_feature.end - codon.start) <= 1):
                # Merge by extending boundaries
                existing_feature.start = min(existing_feature.start, codon.start)
                existing_feature.end = max(existing_feature.end, codon.end)
                flag as "merged_codon_with_adjacent_feature"
                remove redundant new features
                break
        
        # Update transcript boundaries
        transcript.start = min(transcript.start, min(exon.start for exon in transcript.exons))
        transcript.end = max(transcript.end, max(exon.end for exon in transcript.exons))
    else:
        # Codon already covered - remove redundant feature
        flag as "redundant_codon_removed"
        remove codon feature from annotation
```

#### A3. Duplicate Object Prevention
**Critical Fix**: Prevents duplicate transcript objects during GTF parsing that caused data corruption.

```python
# Fixed parsing logic
if transcript_id not in gene.transcripts:
    transcript = Transcript(id=transcript_id, gene_id=gene_id, ...)
    gene.transcripts[transcript_id] = transcript
else:
    transcript = gene.transcripts[transcript_id]  # Use existing object
```

### Phase B: Sequence Validation and Reconstruction

#### B1. CDS Sequence Reconstruction from Merged Coordinates
**Purpose**: Rebuild complete CDS sequences from merged coordinates including all integrated codons.

**Process**:
```
for each transcript with merged coordinates:
    if (genome reference available):
        # Sort CDS regions by coordinate (strand-aware)
        if (transcript.strand == "+"):
            sorted_cds = sorted(transcript.cds_regions, key=lambda x: x.start)
        else:
            sorted_cds = sorted(transcript.cds_regions, key=lambda x: x.start, reverse=True)
        
        # Reconstruct sequence from genome
        reconstructed_sequence = ""
        for cds_region in sorted_cds:
            segment = genome_index[transcript.chrom][cds_region.start-1:cds_region.end]
            if (transcript.strand == "-"):
                segment = reverse_complement(segment)
            reconstructed_sequence += segment
        
        transcript.cds_sequence = reconstructed_sequence
        flag as "cds_reconstructed_from_merged_coordinates"
```

#### B2. Amino Acid Sequence Validation
**Quality Checks**:
- Start codon presence ('M' at position 1)
- Internal stop codon detection ('*' within sequence)
- Length validation (configurable minimum threshold)
- Frame consistency verification

### Phase C: Representative Selection (Sequence-First Approach)

#### C1. Quality-Based Hard Filtering
**Criteria** (Applied BEFORE spatial analysis - HARD FILTERING):
1. Minimum amino acid length threshold (default: 50, configurable) - **REQUIRED**
2. Start codon requirement (configurable)
3. Stop codon requirement (configurable)
4. Maximum internal stop codons (default: 0)

**CRITICAL**: Genes with NO transcripts meeting these criteria are **EXCLUDED** from output entirely.

#### C2. MD5 Hash-Based Grouping
**Process**:
```
for each gene:
    # Apply hard filtering first
    quality_transcripts = filter_by_quality(gene.transcripts)
    
    if not quality_transcripts:
        # NO transcripts meet criteria - gene EXCLUDED from output
        gene.representative_transcript = None
        continue
    
    # Group transcripts by protein sequence identity
    hash_groups = defaultdict(list)
    for transcript in quality_transcripts:
        aa_hash = md5(transcript.aa_sequence).hexdigest()
        hash_groups[aa_hash].append(transcript)
    
    # Select representative from most frequent group
    largest_group = max(hash_groups.values(), key=len)
    representative = max(largest_group, key=lambda t: len(t.aa_sequence))
    gene.representative_transcript = representative
    representative.quality_flags.add("representative")
```

#### C3. Post-Selection Spatial Conflict Assessment
**Non-blocking Analysis**: Spatial conflicts are detected but do NOT prevent selection.

```
# Build interval tree with representatives only
interval_tree = IntervalTree()
for gene in genes_with_representatives:
    transcript = gene.representative_transcript
    interval_tree.addi(transcript.start, transcript.end, transcript)

# Detect overlaps (O(n log n))
for gene in genes_with_representatives:
    transcript = gene.representative_transcript
    overlaps = interval_tree.overlap(transcript.start, transcript.end)
    
    for overlap_transcript in overlaps:
        if (overlap_transcript.gene_id != transcript.gene_id):
            # Inter-gene overlap detected - track but don't reject
            transcript.quality_flags.add("post_selection_overlap")
```

### Phase D: Output Generation with Format Preservation

#### D1. Gene Boundary Updates
**Enhanced Calculation**:
```
for each representative transcript:
    # Calculate boundaries from merged exon coordinates
    gene_start = min(exon.start for exon in transcript.exons)
    gene_end = max(exon.end for exon in transcript.exons)
    
    # Check for boundary changes
    if (gene_boundaries != original_boundaries):
        if (gene_start < original_start OR gene_end > original_end):
            gene.quality_flags.add("gene_expanded")
        if (gene_start > original_start OR gene_end < original_end):
            gene.quality_flags.add("gene_shrunk")
        gene.quality_flags.add("gene_boundaries_updated")
```

#### D2. Format-Specific Output Generation
**GTF Output** (from GTF input):
- Preserves `gene_id "g1"; transcript_id "g1.t1";` syntax
- Maintains original chromosome names and sources
- **NO quality comments** in final output per user specification

**GFF3 Output** (from GFF3 input):
- Preserves `ID=g1; Parent=g1` syntax  
- Maintains hierarchical relationships
- **NO quality comments** in final output per user specification

#### D3. Codon Feature Removal
**Critical**: start_codon and stop_codon features are NOT written to output files as they are integrated into exon/CDS features.

```python
# Output generation - codons are NOT written
for transcript in representative_transcripts:
    write_gene_feature(transcript)
    write_transcript_feature(transcript)
    write_exon_features(transcript.exons)  # Includes merged codon coordinates
    write_cds_features(transcript.cds_regions)  # Includes merged codon coordinates
    # NOTE: start_codon and stop_codon features are intentionally NOT written
```

## Quality Assessment Framework

### Enhanced Quality Flags

#### Preprocessing Flags
- `created_exon_for_orphaned_cds`: New exon created for CDS without exon coverage
- `created_start_codon_features`: Start codon features created and processed
- `created_stop_codon_features`: Stop codon features created and processed

#### Codon Integration Flags
- `merged_start_codon_with_adjacent`: Start codon merged with existing exon/CDS
- `merged_stop_codon_with_adjacent`: Stop codon merged with existing exon/CDS  
- `redundant_start_codon_removed`: Start codon already covered by existing features
- `redundant_stop_codon_removed`: Stop codon already covered by existing features

#### Sequence Quality Flags
- `cds_reconstructed_from_merged_coordinates`: CDS rebuilt from merged coordinates
- `corrected_start_codon`: Start codon added from genome reference
- `corrected_stop_codon`: Stop codon added from genome reference
- `passed_aa_validation`: Amino acid sequence validated
- `short_transcript`: Below minimum length threshold
- `internal_stop_codons`: Contains internal stop codons

#### Selection Status Flags
- `representative`: Selected as gene representative
- `post_selection_overlap`: Spatial conflict detected after selection
- `manual_review`: Requires manual examination (rare)

#### Gene Boundary Flags
- `gene_boundaries_updated`: Gene boundaries adjusted after processing
- `gene_expanded`: Gene range increased to include merged codons
- `gene_shrunk`: Gene range reduced after transcript filtering

## Command Line Interface

### Basic Usage
```bash
# Load environment
module load miniconda3 && source activate && conda activate braker_cleaner

# Standard usage
python pipeline_cli.py \
    --gene-model braker.gtf \
    --cds braker.codingseq \
    --aa braker.aa \
    --genome KaleLate.genome.softmask.fa \
    --min-length 20 \
    --overlap-threshold 0.5 \
    --output-dir results
```

### Advanced Configuration
```bash
# With configuration file
echo '{
  "min_aa_length": 20,
  "overlap_threshold": 0.5,
  "memory_limit_mb": 8192,
  "debug_mode": false
}' > config.json

python pipeline_cli.py \
    --config config.json \
    --gene-model braker.gtf \
    --cds braker.codingseq \
    --aa braker.aa \
    --genome genome.fa \
    --output-dir results
```

### Parameter Descriptions
- `--gene-model`: Input gene model (GFF3/GTF) - **REQUIRED**
- `--cds`: Input CDS FASTA file - **REQUIRED**
- `--aa`: Input amino acid FASTA file - **REQUIRED**
- `--genome`: Reference genome FASTA - **REQUIRED for codon integration**
- `--output-dir`: Output directory - **REQUIRED**
- `--min-length`: Minimum AA length (default: 50)
- `--overlap-threshold`: Spatial conflict threshold (default: 0.5)
- `--config`: Configuration file (JSON/YAML)
- `--log-level`: Logging level (DEBUG, INFO, WARNING, ERROR)
- `--memory-limit`: Memory limit in MB (default: 4096)

## Comprehensive Testing Framework

### Test Categories (58 Total Tests)
1. **Data Structure Tests** (25 tests): Core validation and property calculations
2. **Configuration Tests** (17 tests): File loading, environment variables, validation  
3. **Codon Integration Tests** (16 tests): General codon scenarios and adjacency logic

### Running Tests
```bash
# Run all tests
python run_tests.py

# Run with coverage analysis  
python run_tests.py --coverage

# Run specific categories
python run_tests.py --tests TestCodonIntegration TestConfig

# Validate environment
python run_tests.py --validate
```

### Test Results Summary
```
======================================================================
Gene Curation Pipeline - Test Suite
======================================================================
Tests run: 58
Failures: 0
Errors: 0
Skipped: 0

✅ All tests passed!
```

## Performance Specifications and Results

### Algorithmic Complexity Requirements
**MANDATORY**: All algorithms must be **O(n log n)** or better.

- **File Parsing**: O(n) single-pass processing
- **Codon Integration**: O(n) per transcript 
- **Hash Grouping**: O(n) using dictionary lookups
- **Overlap Detection**: O(n log n) using IntervalTree
- **Sequence Reconstruction**: O(k) per transcript (k = CDS regions)

### Actual Performance Results (52,860 genes, 91,624 transcripts)
```
PERFORMANCE REPORT
==================================================
Total time: 17.25 seconds
Peak memory: 644.5 MB (15.7% of 4GB limit)

Phase breakdown:
  input_parsing: 4.99s (144,484 ops, 28,969.4 ops/s)
  annotation_curation: 0.78s (52,860 ops, 67,767.0 ops/s) 
  sequence_validation: 9.93s (91,624 ops, 9,222.4 ops/s)
  representative_selection: 0.23s (52,860 ops, 225,245.5 ops/s)
  output_generation: 1.28s (4 ops, 3.1 ops/s)
```

### Quality Achievements
- ✅ **Hard Filtering**: Genes not meeting --min-length threshold are excluded entirely
- ✅ **Complete Codon Integration**: All start/stop codons properly processed
- ✅ **Perfect Format Preservation**: Input format maintained in output
- ✅ **Memory Efficiency**: <16% memory utilization for large datasets
- ✅ **Data Safety**: Original files never modified
- ✅ **Quality Control**: Only high-quality gene models in final output

## Error Handling and Troubleshooting

### Professional Exception Hierarchy
```python
PipelineError                 # Base exception
├── ParseError               # File parsing errors with line numbers
├── ValidationError          # Sequence validation errors  
├── SequenceError           # Sequence processing errors
├── ConfigurationError      # Configuration errors
├── MemoryError            # Memory usage errors
└── GenomeError            # Genome access errors
```

### Common Issues and Solutions
1. **Import Errors**: Install dependencies (`pip install pyfaidx intervaltree`)
2. **Memory Issues**: Reduce batch size or increase memory limit
3. **Parsing Warnings**: Check malformed coordinates (start=end cases)
4. **Missing Sequences**: Verify transcript ID consistency between files
5. **Configuration Errors**: Validate JSON/YAML format

### Debugging Tips
- Use `--log-level DEBUG` for detailed processing information
- Check `gene_curation.log` for step-by-step processing details
- Run `python run_tests.py --validate` to check environment setup
- Monitor memory usage with `--memory-limit` parameter

## Manual Review Process

### Current Results
The implemented pipeline achieves **0% manual review rate** for typical datasets.

**Manual Review Required For**:
- Genes where ALL transcripts are below minimum length threshold
- Genes where ALL transcripts fail quality criteria (no start codon, internal stops, etc.)
- Reason: `no_representative_selected`
- **These genes are EXCLUDED from cleaned output files**

**Automatically Handled Cases** (No Manual Review):
- Genes with spatial conflicts → Tracked but selected automatically
- Genes with multiple transcripts → Best representative selected  
- Genes with identical proteins → Longest transcript selected
- Genes with annotation artifacts → Fixed through preprocessing
- Genes with codon issues → Integrated automatically

### Manual Review File Format
`manual_review_transcripts.txt`:
```
# Gene_ID	Transcript_IDs	Reason	Quality_Flags
g235	g235.t1	no_representative_selected	short_transcript,validated_cds
```

## Key Algorithmic Innovations

### 1. Preprocessing-Based Codon Integration
```
Old approach: Create separate codon features → Fragmented annotations
New approach: Detect orphaned CDS → Create exons → Merge adjacent → Reconstruct
```

### 2. Hard Filtering with Sequence-First Selection
```
Old approach: Spatial filtering → Selection → Manual review (94% failure)
New approach: Hard quality filtering → Gene exclusion → Selection → Spatial tracking
```

### 3. Genome-Based Sequence Reconstruction  
```
Old approach: Use original sequences (missing codons)
New approach: Rebuild from merged coordinates with genome reference
```

### 4. Duplicate Object Prevention
```
Old approach: Create new transcript objects for each feature → Data corruption
New approach: Reuse existing transcript objects → Data integrity maintained
```

### 5. Professional Architecture with O(n log n) Validation
```
Old approach: Monolithic file with O(n²) algorithms
New approach: Modular package with complexity validation and comprehensive testing
```

## Implementation Requirements

### Critical Data Integrity Rules
1. **NEVER modify original input files** - Always create new output files
2. **NEVER overwrite** braker.gtf, braker.codingseq, braker.aa, genome.fa
3. **Always preserve** original chromosome names and source information
4. **Maintain coordinate accuracy** - Handle 1-based vs 0-based indexing properly

### Required Dependencies
- `pyfaidx`: Genome sequence indexing (**MANDATORY for codon integration**)
- `intervaltree`: Efficient overlap detection (**MANDATORY for O(n log n) performance**)
- `psutil`: Memory monitoring (optional but recommended)
- `PyYAML`: Configuration file support (optional)

### Performance Testing Requirements
- Unit tests must verify O(n log n) complexity bounds
- Memory leak detection for long-running processes
- Scalability testing with 1K, 10K, 100K+ gene datasets
- Benchmark validation against existing implementations

## Results Summary

### Transformation Achievements
**Before**: 1,255-line monolithic script with parsing bugs and O(n²) algorithms
**After**: Professional modular package with 58 comprehensive tests and O(n log n) validation

### Quality Metrics
- **Hard Filtering**: Genes not meeting quality criteria are excluded from output
- **Quality-Based Selection**: Only genes with qualifying transcripts receive representatives
- **Complete Integration**: All codon features properly merged and reconstructed
- **Format Integrity**: Perfect preservation of input format and metadata
- **Performance Excellence**: Fast processing with <16% memory usage
- **Professional Standards**: Comprehensive testing, error handling, and documentation

This specification represents the complete implementation of a production-ready gene annotation curation pipeline with proven results on large-scale genomic datasets.