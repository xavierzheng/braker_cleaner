#!/usr/bin/env python3

"""
Core data structures for the gene curation pipeline.

Defines the fundamental data classes for representing genes, transcripts,
exons, CDS regions, and codons with proper relationships and methods.
"""

import hashlib
from dataclasses import dataclass, field
from typing import Dict, List, Set, Optional


@dataclass
class Codon:
    """Represents a start or stop codon."""
    start: int
    end: int
    strand: str
    codon_type: str  # 'start' or 'stop'
    sequence: str = ""
    
    def __post_init__(self):
        """Validate codon data after initialization."""
        if self.start >= self.end:
            raise ValueError(f"Invalid codon coordinates: {self.start}-{self.end}")
        if self.codon_type not in ('start', 'stop'):
            raise ValueError(f"Invalid codon type: {self.codon_type}")
        if self.strand not in ('+', '-'):
            raise ValueError(f"Invalid strand: {self.strand}")


@dataclass
class Exon:
    """Represents an exon."""
    start: int
    end: int
    strand: str
    id: str = ""
    
    def __post_init__(self):
        """Validate exon data after initialization."""
        if self.start >= self.end:
            raise ValueError(f"Invalid exon coordinates: {self.start}-{self.end}")
        if self.strand not in ('+', '-'):
            raise ValueError(f"Invalid strand: {self.strand}")
    
    @property
    def length(self) -> int:
        """Get exon length."""
        return self.end - self.start + 1
    
    def overlaps_with(self, other: 'Exon') -> bool:
        """Check if this exon overlaps with another."""
        return not (self.end < other.start or self.start > other.end)
    
    def contains_position(self, position: int) -> bool:
        """Check if exon contains a specific position."""
        return self.start <= position <= self.end


@dataclass
class CDS:
    """Represents a CDS region."""
    start: int
    end: int
    strand: str
    phase: int = 0
    id: str = ""
    
    def __post_init__(self):
        """Validate CDS data after initialization."""
        if self.start >= self.end:
            raise ValueError(f"Invalid CDS coordinates: {self.start}-{self.end}")
        if self.strand not in ('+', '-'):
            raise ValueError(f"Invalid strand: {self.strand}")
        if self.phase not in (0, 1, 2):
            raise ValueError(f"Invalid phase: {self.phase}")
    
    @property
    def length(self) -> int:
        """Get CDS length."""
        return self.end - self.start + 1
    
    def overlaps_with(self, other: 'CDS') -> bool:
        """Check if this CDS overlaps with another."""
        return not (self.end < other.start or self.start > other.end)
    
    def contains_position(self, position: int) -> bool:
        """Check if CDS contains a specific position."""
        return self.start <= position <= self.end


@dataclass
class Transcript:
    """Represents a transcript with all its features."""
    id: str
    gene_id: str
    start: int
    end: int
    strand: str
    exons: List[Exon] = field(default_factory=list)
    cds_regions: List[CDS] = field(default_factory=list)
    start_codon: Optional[Codon] = None
    stop_codon: Optional[Codon] = None
    quality_flags: Set[str] = field(default_factory=set)
    cds_sequence: str = ""
    aa_sequence: str = ""
    # Original annotation info
    chrom: str = ""
    source: str = ""
    
    def __post_init__(self):
        """Validate transcript data after initialization."""
        if self.start >= self.end:
            raise ValueError(f"Invalid transcript coordinates: {self.start}-{self.end}")
        if self.strand not in ('+', '-'):
            raise ValueError(f"Invalid strand: {self.strand}")
        if not self.id:
            raise ValueError("Transcript ID cannot be empty")
        if not self.gene_id:
            raise ValueError("Gene ID cannot be empty")
    
    def get_cds_hash(self) -> str:
        """Get MD5 hash of CDS sequence."""
        if not self.cds_sequence:
            return ""
        return hashlib.md5(self.cds_sequence.encode()).hexdigest()
    
    def get_aa_hash(self) -> str:
        """Get MD5 hash of amino acid sequence."""
        if not self.aa_sequence:
            return ""
        return hashlib.md5(self.aa_sequence.encode()).hexdigest()
    
    def get_coordinates_hash(self) -> str:
        """Get MD5 hash of coding region coordinates."""
        if not self.cds_regions:
            return ""
        coord_str = ";".join([f"{cds.start}-{cds.end}" for cds in sorted(self.cds_regions, key=lambda x: x.start)])
        return hashlib.md5(coord_str.encode()).hexdigest()
    
    @property
    def length(self) -> int:
        """Get transcript length."""
        return self.end - self.start + 1
    
    @property
    def exon_count(self) -> int:
        """Get number of exons."""
        return len(self.exons)
    
    @property
    def cds_count(self) -> int:
        """Get number of CDS regions."""
        return len(self.cds_regions)
    
    @property
    def total_exon_length(self) -> int:
        """Get total length of all exons."""
        return sum(exon.length for exon in self.exons)
    
    @property
    def total_cds_length(self) -> int:
        """Get total length of all CDS regions."""
        return sum(cds.length for cds in self.cds_regions)
    
    @property
    def aa_length(self) -> int:
        """Get amino acid sequence length."""
        return len(self.aa_sequence.rstrip('*'))  # Exclude stop codon
    
    def add_quality_flag(self, flag: str) -> None:
        """Add a quality flag."""
        self.quality_flags.add(flag)
    
    def has_quality_flag(self, flag: str) -> bool:
        """Check if transcript has a specific quality flag."""
        return flag in self.quality_flags
    
    def get_sorted_exons(self) -> List[Exon]:
        """Get exons sorted by coordinate (strand-aware)."""
        if self.strand == '+':
            return sorted(self.exons, key=lambda x: x.start)
        else:
            return sorted(self.exons, key=lambda x: x.start, reverse=True)
    
    def get_sorted_cds(self) -> List[CDS]:
        """Get CDS regions sorted by coordinate (strand-aware)."""
        if self.strand == '+':
            return sorted(self.cds_regions, key=lambda x: x.start)
        else:
            return sorted(self.cds_regions, key=lambda x: x.start, reverse=True)
    
    def codon_is_adjacent_to_features(self, codon: Codon, gap_threshold: int = 1) -> bool:
        """Check if codon is adjacent to existing exons/CDS within gap threshold."""
        for exon in self.exons:
            # Check if codon is adjacent to exon (within gap threshold)
            if (abs(codon.end - exon.start) <= gap_threshold or 
                abs(exon.end - codon.start) <= gap_threshold):
                return True
        
        for cds in self.cds_regions:
            # Check if codon is adjacent to CDS (within gap threshold)
            if (abs(codon.end - cds.start) <= gap_threshold or 
                abs(cds.end - codon.start) <= gap_threshold):
                return True
        
        return False
    
    def codon_is_covered_by_features(self, codon: Codon) -> bool:
        """Check if codon is already covered by existing exons AND CDS."""
        exon_covered = any(exon.contains_position(codon.start) and exon.contains_position(codon.end) 
                          for exon in self.exons)
        cds_covered = any(cds.contains_position(codon.start) and cds.contains_position(codon.end) 
                         for cds in self.cds_regions)
        return exon_covered and cds_covered


@dataclass
class Gene:
    """Represents a gene with multiple transcripts."""
    id: str
    transcripts: List[Transcript] = field(default_factory=list)
    representative: Optional[Transcript] = None
    
    def __post_init__(self):
        """Validate gene data after initialization."""
        if not self.id:
            raise ValueError("Gene ID cannot be empty")
    
    def get_transcript_by_id(self, transcript_id: str) -> Optional[Transcript]:
        """Get transcript by ID."""
        for transcript in self.transcripts:
            if transcript.id == transcript_id:
                return transcript
        return None
    
    def add_transcript(self, transcript: Transcript) -> None:
        """Add a transcript to the gene."""
        if transcript.gene_id != self.id:
            raise ValueError(f"Transcript gene_id ({transcript.gene_id}) doesn't match gene id ({self.id})")
        
        # Check for duplicate transcript IDs
        if any(t.id == transcript.id for t in self.transcripts):
            raise ValueError(f"Transcript {transcript.id} already exists in gene {self.id}")
        
        self.transcripts.append(transcript)
    
    @property
    def transcript_count(self) -> int:
        """Get number of transcripts."""
        return len(self.transcripts)
    
    @property
    def has_representative(self) -> bool:
        """Check if gene has a representative transcript."""
        return self.representative is not None
    
    def get_gene_boundaries(self) -> tuple[int, int]:
        """Get gene boundaries based on all transcripts."""
        if not self.transcripts:
            return 0, 0
        
        starts = [t.start for t in self.transcripts]
        ends = [t.end for t in self.transcripts]
        return min(starts), max(ends)
    
    def get_representative_boundaries(self) -> tuple[int, int]:
        """Get gene boundaries based on representative transcript's exons."""
        if not self.representative or not self.representative.exons:
            return self.get_gene_boundaries()
        
        exon_starts = [e.start for e in self.representative.exons]
        exon_ends = [e.end for e in self.representative.exons]
        return min(exon_starts), max(exon_ends)
    
    def set_representative(self, transcript: Transcript) -> None:
        """Set the representative transcript."""
        if transcript not in self.transcripts:
            raise ValueError(f"Transcript {transcript.id} is not part of gene {self.id}")
        
        self.representative = transcript
        transcript.add_quality_flag("representative")
    
    def get_transcript_by_quality_flags(self, required_flags: Set[str]) -> List[Transcript]:
        """Get transcripts that have all required quality flags."""
        return [t for t in self.transcripts 
                if required_flags.issubset(t.quality_flags)]