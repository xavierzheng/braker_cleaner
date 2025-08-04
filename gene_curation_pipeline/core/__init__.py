#!/usr/bin/env python3

"""
Core module for the gene curation pipeline.

Contains fundamental data structures, exception types, and configuration
management components.
"""

from .data_structures import Gene, Transcript, Exon, CDS, Codon
from .exceptions import (
    PipelineError, ParseError, ValidationError, SequenceError,
    ConfigurationError, MemoryError, GenomeError
)
from .config import PipelineConfig, load_config

__all__ = [
    'Gene', 'Transcript', 'Exon', 'CDS', 'Codon',
    'PipelineError', 'ParseError', 'ValidationError', 'SequenceError',
    'ConfigurationError', 'MemoryError', 'GenomeError',
    'PipelineConfig', 'load_config'
]