# Gene Annotation Curation Pipeline

A professional, modular gene annotation curation pipeline with comprehensive codon integration, O(n log n) complexity, and zero manual review requirements.

## Overview

This pipeline implements a sequence-first, spatial-second approach to gene annotation curation with **complete start/stop codon integration**. The pipeline has been completely refactored into a professional modular architecture with comprehensive testing and maintains 100% selection success rate with 0% manual review rate.

## Key Features

- **🏗️ Professional Modular Architecture**: Clean separation of concerns with dedicated modules
- **🧪 Comprehensive Testing**: 58 unit tests with 100% pass rate
- **⚙️ Centralized Configuration**: Flexible configuration management with file/environment support
- **📊 Enhanced Performance Monitoring**: Real-time memory and complexity validation
- **🔗 Complete Codon Integration**: Merges start/stop codons with adjacent exon/CDS features
- **🧬 Genome-Based Sequence Reconstruction**: Rebuilds CDS sequences from merged coordinates
- **📄 Universal Format Support**: Handles both GFF3 and GTF formats with automatic detection
- **✅ Perfect Selection Success**: 100% of genes get representative transcripts
- **🚫 Zero Manual Review**: Eliminates false positives through pragmatic selection
- **🧠 Sequence-First Approach**: Prioritizes biological equivalence over spatial conflicts
- **📏 Enhanced Gene Boundaries**: Updates gene coordinates based on merged codon regions
- **⚡ O(n log n) Complexity**: Efficient processing for large datasets

## Quick Start

```bash
# Activate conda environment
module load miniconda3 && source activate && conda activate braker_cleaner

# Run with GTF format and genome reference (RECOMMENDED)
python pipeline_cli.py \
    --gene-model braker.gtf \
    --cds braker.codingseq \
    --aa braker.aa \
    --genome KaleLate.genome.softmask.fa \
    --min-length 20 \
    --overlap-threshold 0.5 \
    --output-dir results

# Run with GFF3 format
python pipeline_cli.py \
    --gene-model braker.gff3 \
    --cds braker.cds.fa \
    --aa braker.aa \
    --genome genome.fa \
    --output-dir results
```

## New Modular Architecture

### Before vs After
**Before**: Single 1,255-line `gene_curation_pipeline.py` file
**After**: Professional modular package with comprehensive testing

```
gene_curation_pipeline/
├── __init__.py                   # Main package interface
├── core/
│   ├── __init__.py
│   ├── data_structures.py        # Core data classes with validation
│   ├── exceptions.py             # Specific exception hierarchy
│   ├── config.py                 # Centralized configuration
│   ├── pipeline.py               # Main pipeline orchestration
│   ├── parsers.py                # GFF3/GTF and FASTA parsing
│   ├── processors.py             # Annotation, sequence, and selection processing
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

## Installation

### Prerequisites
- Python 3.11+
- Conda environment with required dependencies

### Required Dependencies
```bash
# Core dependencies (required)
pip install pyfaidx intervaltree

# Optional dependencies (for enhanced features)
pip install psutil PyYAML coverage
```

### Environment Setup
```bash
# Load conda environment
module load miniconda3
source activate
conda activate braker_cleaner
```

## Input Files

- **Gene Model**: GFF3 or GTF format annotation file (`--gene-model`)
  - Supports both `.gff3` and `.gtf` extensions
  - Automatic format detection and preservation
  - Handles standard and non-standard attribute formats
- **CDS Sequences**: FASTA file with coding sequences (`--cds`)
- **AA Sequences**: FASTA file with amino acid sequences (`--aa`)
- **Genome Sequence**: Reference genome FASTA (`--genome`, **REQUIRED for proper codon integration**)
  - **CRITICAL**: Required for merging start/stop codons with adjacent features
  - **CRITICAL**: Required for accurate CDS sequence reconstruction

## Output Files

- **Format-Preserved Annotations**: 
  - GTF input → `cleaned.gtf` (GTF format with merged codon coordinates)
  - GFF3 input → `cleaned.gff3` (GFF3 format with merged codon coordinates)
- **Reconstructed Sequences**: 
  - `cleaned.cds.fa` - CDS sequences rebuilt from merged coordinates
  - `cleaned.aa` - Validated amino acid sequences
- **Comprehensive Reports**:
  - `processing_report.txt` - Detailed pipeline statistics and performance metrics
  - `manual_review_transcripts.txt` - Genes requiring manual attention (typically empty)
  - `gene_curation.log` - Detailed processing log

## Performance Results

### Actual Test Results (52,860 genes, 91,624 transcripts)
- ✅ **100% Selection Success**: 52,860/52,860 genes get representatives
- ✅ **0% Manual Review**: 0 genes requiring manual review
- ✅ **Fast Processing**: 17.25 seconds total runtime
- ✅ **Memory Efficient**: 644.5 MB peak memory (15.7% of 4GB limit)
- ✅ **Complete Codon Integration**: All start/stop codons properly processed
- ✅ **Accurate Sequence Reconstruction**: CDS sequences rebuilt from merged coordinates

### Phase Performance Breakdown
```
Phase                    | Time    | Operations | Ops/sec   | Memory
-------------------------|---------|------------|-----------|--------
Input Parsing           | 4.99s   | 144,484    | 28,969.4  | 437.4MB
Annotation Curation     | 0.78s   | 52,860     | 67,767.0  | 494.7MB
Sequence Validation     | 9.93s   | 91,624     | 9,222.4   | 619.7MB
Representative Selection| 0.23s   | 52,860     | 225,245.5 | 644.2MB
Output Generation       | 1.28s   | 4          | 3.1       | 644.5MB
```

## Configuration Management

### Configuration Options
The pipeline supports flexible configuration through multiple methods:

#### Default Configuration
```python
memory_limit_mb: 4096
batch_size: 1000
min_aa_length: 50
overlap_threshold: 0.5
enable_codon_integration: True
enable_sequence_reconstruction: True
adjacency_gap_threshold: 1
```

#### File-Based Configuration
Create `config.json` or `config.yaml`:
```json
{
  "min_aa_length": 20,
  "overlap_threshold": 0.5,
  "memory_limit_mb": 2048,
  "debug_mode": true
}
```

#### Environment Variables
```bash
export PIPELINE_MIN_AA_LENGTH=20
export PIPELINE_MEMORY_LIMIT_MB=2048
export PIPELINE_DEBUG_MODE=true
```

#### Command Line Usage
```bash
python pipeline_cli.py \
    --config config.json \
    --gene-model braker.gtf \
    --cds braker.codingseq \
    --aa braker.aa \
    --genome genome.fa \
    --output-dir results
```

## Command Line Options

- `--gene-model`: Input gene model file - supports GFF3 and GTF formats (required)
- `--cds`: Input CDS FASTA file (required)
- `--aa`: Input AA FASTA file (required)
- `--genome`: Reference genome FASTA file (**REQUIRED for proper codon integration**)
- `--output-dir`: Output directory (required)
- `--min-length`: Minimum AA length threshold (default: 50)
- `--overlap-threshold`: Overlap threshold for spatial conflicts (default: 0.5)
- `--config`: Configuration file (JSON or YAML format)
- `--log-level`: Logging level (DEBUG, INFO, WARNING, ERROR)
- `--memory-limit`: Memory limit in MB (default: 4096)
- `--batch-size`: Batch size for processing (default: 1000)

## Testing Framework

### Comprehensive Unit Tests (58 tests)
```bash
# Run all tests
python run_tests.py

# Run with coverage analysis
python run_tests.py --coverage

# Run specific test categories
python run_tests.py --tests TestCodonIntegration TestConfig

# Validate test environment
python run_tests.py --validate
```

### Test Categories
1. **Data Structure Tests** (25 tests): Core validation and property calculations
2. **Configuration Tests** (17 tests): File loading, environment variables, validation
3. **Codon Integration Tests** (16 tests): General codon scenarios and adjacency logic

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

## Error Handling

### Professional Exception Hierarchy
```python
PipelineError                 # Base exception
├── ParseError               # File parsing errors
├── ValidationError          # Sequence validation errors
├── SequenceError           # Sequence processing errors
├── ConfigurationError      # Configuration errors
├── MemoryError            # Memory usage errors
└── GenomeError            # Genome access errors
```

### Error Context Information
- Filename and line numbers for parsing errors
- Coordinate information for sequence errors
- Configuration parameter details for validation errors
- Memory usage details for resource errors

## Examples

### Basic Usage
```bash
# Standard processing with codon integration
python pipeline_cli.py \
    --gene-model braker.gtf \
    --cds braker.codingseq \
    --aa braker.aa \
    --genome KaleLate.genome.softmask.fa \
    --output-dir results
```

### Production Environment
```bash
# Full environment setup with custom parameters
module load miniconda3
source activate
conda activate braker_cleaner

python pipeline_cli.py \
    --gene-model braker.gtf \
    --cds braker.codingseq \
    --aa braker.aa \
    --genome KaleLate.genome.softmask.fa \
    --min-length 20 \
    --overlap-threshold 0.5 \
    --memory-limit 8192 \
    --log-level INFO \
    --output-dir production_results
```

### With Configuration File
```bash
# Create config.json
echo '{
  "min_aa_length": 20,
  "overlap_threshold": 0.5,
  "memory_limit_mb": 8192,
  "debug_mode": false
}' > config.json

# Run with configuration
python pipeline_cli.py \
    --config config.json \
    --gene-model braker.gtf \
    --cds braker.codingseq \
    --aa braker.aa \
    --genome KaleLate.genome.softmask.fa \
    --output-dir configured_results
```

## Implementation Highlights

### Enhanced Codon Integration
1. **Adjacency Detection**: Check ≤1bp gap between codon and existing features
2. **Coordinate Merging**: Extend existing features to include codon coordinates
3. **Sequence Reconstruction**: Rebuild CDS from merged coordinates using genome
4. **Boundary Updates**: Update transcript and gene coordinates after merging

### Sequence-First Selection
1. **Quality Filtering**: Remove transcripts based on sequence quality only
2. **MD5 Hashing**: Group transcripts by protein sequence identity
3. **Representative Selection**: Choose longest from most frequent hash group
4. **Spatial Assessment**: Post-selection overlap detection (non-blocking)

### Professional Architecture Benefits
- **Maintainable**: Clear module separation and single responsibilities
- **Testable**: Comprehensive unit test coverage with isolated components
- **Configurable**: Flexible configuration for different environments
- **Monitorable**: Real-time performance tracking and validation
- **Extensible**: Easy to add new features and processing phases

## Performance Monitoring

### Real-Time Monitoring
- Memory usage tracking with configurable limits
- Processing time measurement for each phase
- Operation counting and throughput calculation
- Algorithmic complexity validation

### Performance Reports
```
PERFORMANCE REPORT
==================================================
Total time: 17.25 seconds
Peak memory: 644.5 MB
Memory limit: 4096 MB
Memory utilization: 15.7%

Phase breakdown:
  input_parsing: 4.99s (144484 ops, 28969.4 ops/s, 437.4MB)
  annotation_curation: 0.78s (52860 ops, 67767.0 ops/s, 494.7MB)
  sequence_validation: 9.93s (91624 ops, 9222.4 ops/s, 619.7MB)
  representative_selection: 0.23s (52860 ops, 225245.5 ops/s, 644.2MB)
  output_generation: 1.28s (4 ops, 3.1 ops/s, 644.5MB)
```

## Quality Assurance

### Enhanced Quality Flags
- **Codon Integration**: `created_start_codon_features`, `merged_stop_codon_with_adjacent`
- **Sequence Validation**: `cds_reconstructed_from_merged_coordinates`, `passed_aa_validation`
- **Selection Status**: `representative`, `post_selection_overlap`
- **Boundary Updates**: `gene_boundaries_updated`, `gene_expanded`

### Data Integrity
- **Original Files Protected**: Never modifies input files
- **Format Preservation**: Maintains original chromosome names and sources
- **Coordinate Accuracy**: Proper handling of 1-based vs 0-based indexing
- **Sequence Consistency**: CDS sequences match merged annotation coordinates

## Troubleshooting

### Common Issues
1. **Import Errors**: Ensure all dependencies are installed (`pip install pyfaidx intervaltree`)
2. **Memory Issues**: Reduce batch size or increase memory limit
3. **Parsing Warnings**: Check for malformed coordinate entries (e.g., start=end)
4. **Missing Sequences**: Verify transcript ID consistency between files
5. **Configuration Errors**: Validate configuration file format (JSON/YAML)

### Performance Issues
- **Slow Processing**: Check memory usage and enable batch processing
- **High Memory Usage**: Reduce batch size or use streaming processing
- **Test Failures**: Run `python run_tests.py --validate` to check environment

### File Format Issues
- **Auto-Detection**: Ensure proper file extensions (.gtf, .gff3)
- **Attribute Parsing**: Pipeline handles both standard and non-standard formats
- **Coordinate Systems**: Verify genome and annotation coordinate consistency

## Algorithm Complexity

All algorithms maintain **O(n log n)** or better complexity:

- **File Parsing**: O(n) single-pass processing
- **Codon Integration**: O(n) adjacency checking per transcript
- **Hash Grouping**: O(n) using dictionary lookups
- **Overlap Detection**: O(n log n) using IntervalTree
- **Sequence Reconstruction**: O(k) per transcript (k = CDS regions)

## Files in This Package

### Core Pipeline Files
- `pipeline_cli.py`: Command-line interface (maintains original API)
- `gene_curation_pipeline/`: Modular package directory
- `run_tests.py`: Comprehensive test runner
- `CLAUDE.md`: Complete specification document

### Generated Documentation
- `IMPROVEMENTS_SUMMARY.md`: Detailed summary of all architectural improvements
- Test coverage reports in `htmlcov/` (when running with `--coverage`)

### Legacy Files
- `gene_curation_pipeline.py`: Original monolithic implementation (for reference)

## Results Summary

### Test Dataset Performance (KaleLate Dataset)
- **Total genes**: 52,860
- **Total transcripts**: 91,624
- **Representatives selected**: 52,860 (100.0%)
- **Manual review needed**: 0 (0.0%)
- **Processing time**: 17.25 seconds
- **Memory usage**: 644.5 MB peak
- **Codon integration**: Complete success
- **Format preservation**: Perfect GTF maintenance

### Quality Achievements
- ✅ **100% Selection Success**: Every gene gets a representative
- ✅ **0% Manual Review Rate**: No genes require manual intervention
- ✅ **Complete Codon Integration**: All start/stop codons properly merged
- ✅ **Professional Architecture**: 58 comprehensive unit tests with 100% pass rate
- ✅ **Performance Excellence**: O(n log n) complexity validated
- ✅ **Memory Efficiency**: <16% memory utilization for large datasets
- ✅ **Format Integrity**: Perfect preservation of original file formats
- ✅ **Data Safety**: Original files never modified

## Citation

If you use this pipeline in your research, please cite:

```
Gene Annotation Curation Pipeline v2.0
Professional modular gene annotation curation with complete codon integration,
comprehensive testing framework, and O(n log n) complexity validation.
```

## License

This pipeline is provided under the terms specified in the project documentation. For questions about usage, please refer to the comprehensive test suite and performance validation results.