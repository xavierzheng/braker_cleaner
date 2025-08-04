#!/usr/bin/env python3

"""
Custom exceptions for the gene curation pipeline.

Provides specific exception types for better error handling and debugging.
"""

class PipelineError(Exception):
    """Base exception for all pipeline-related errors."""
    pass


class ParseError(PipelineError):
    """Error occurred during file parsing."""
    
    def __init__(self, message: str, filename: str = "", line_number: int = 0):
        super().__init__(message)
        self.filename = filename
        self.line_number = line_number
    
    def __str__(self):
        if self.filename and self.line_number:
            return f"Parse error in {self.filename} at line {self.line_number}: {super().__str__()}"
        elif self.filename:
            return f"Parse error in {self.filename}: {super().__str__()}"
        return super().__str__()


class ValidationError(PipelineError):
    """Error occurred during sequence or annotation validation."""
    
    def __init__(self, message: str, transcript_id: str = "", gene_id: str = ""):
        super().__init__(message)
        self.transcript_id = transcript_id
        self.gene_id = gene_id
    
    def __str__(self):
        if self.transcript_id:
            return f"Validation error for transcript {self.transcript_id}: {super().__str__()}"
        elif self.gene_id:
            return f"Validation error for gene {self.gene_id}: {super().__str__()}"
        return super().__str__()


class SequenceError(PipelineError):
    """Error occurred during sequence processing."""
    
    def __init__(self, message: str, sequence_id: str = "", sequence_type: str = ""):
        super().__init__(message)
        self.sequence_id = sequence_id
        self.sequence_type = sequence_type
    
    def __str__(self):
        if self.sequence_id and self.sequence_type:
            return f"Sequence error in {self.sequence_type} {self.sequence_id}: {super().__str__()}"
        elif self.sequence_id:
            return f"Sequence error in {self.sequence_id}: {super().__str__()}"
        return super().__str__()


class ConfigurationError(PipelineError):
    """Error in pipeline configuration."""
    pass


class MemoryError(PipelineError):
    """Memory usage exceeded limits."""
    
    def __init__(self, message: str, current_usage: float, limit: float):
        super().__init__(message)
        self.current_usage = current_usage
        self.limit = limit
    
    def __str__(self):
        return f"Memory error: {super().__str__()} (current: {self.current_usage:.1f}MB, limit: {self.limit:.1f}MB)"


class GenomeError(PipelineError):
    """Error accessing genome reference."""
    
    def __init__(self, message: str, chromosome: str = "", coordinates: str = ""):
        super().__init__(message)
        self.chromosome = chromosome
        self.coordinates = coordinates
    
    def __str__(self):
        if self.chromosome and self.coordinates:
            return f"Genome error at {self.chromosome}:{self.coordinates}: {super().__str__()}"
        elif self.chromosome:
            return f"Genome error at {self.chromosome}: {super().__str__()}"
        return super().__str__()