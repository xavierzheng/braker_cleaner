#!/usr/bin/env python3

"""
Gene Annotation Curation Pipeline

A high-performance pipeline for curating BRAKER gene predictions with complete
codon integration, O(n log n) complexity, and zero manual review requirements.

This modular implementation provides:
- Clean separation of concerns with dedicated modules
- Comprehensive error handling with specific exception types
- Centralized configuration management
- Extensive unit test coverage
- Professional software engineering practices

Modules:
- core: Fundamental data structures, exceptions, and configuration
- processors: Main pipeline processing components
- output: Output generation and reporting
- utils: Utility functions and performance monitoring
- tests: Comprehensive test suite
"""

__version__ = "2.0.0"
__author__ = "Gene Curation Pipeline Team"

# Import main components for easy access
from .core.data_structures import Gene, Transcript, Exon, CDS, Codon
from .core.exceptions import (
    PipelineError, ParseError, ValidationError, SequenceError,
    ConfigurationError, MemoryError, GenomeError
)
from .core.config import PipelineConfig, load_config
from .core.pipeline import GeneCurationPipeline

__all__ = [
    # Main pipeline
    'GeneCurationPipeline',
    # Data structures
    'Gene', 'Transcript', 'Exon', 'CDS', 'Codon',
    # Exceptions
    'PipelineError', 'ParseError', 'ValidationError', 'SequenceError',
    'ConfigurationError', 'MemoryError', 'GenomeError',
    # Configuration
    'PipelineConfig', 'load_config'
]