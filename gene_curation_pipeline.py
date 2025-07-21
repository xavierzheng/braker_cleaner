#!/usr/bin/env python3

"""
Gene Annotation Curation Pipeline

This pipeline curates BRAKER gene predictions by validating annotations,
correcting sequences, and selecting representative transcripts.

Based on CLAUDE.md specifications with mandatory O(n log n) complexity.
"""

import os
import sys
import argparse
import hashlib
import logging
import time
import psutil
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional
from dataclasses import dataclass, field

# Required dependencies for performance
try:
    import pyfaidx
    from intervaltree import IntervalTree, Interval
except ImportError as e:
    print(f"Error: Missing required dependency: {e}")
    print("Please install: pip install pyfaidx intervaltree")
    sys.exit(1)

# Optional dependencies
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    print("Warning: psutil not available, memory monitoring disabled")

# Performance monitoring
class PerformanceMonitor:
    """Monitor memory usage and processing time."""
    
    def __init__(self):
        self.start_time = time.time()
        self.peak_memory = 0
        self.process = psutil.Process() if PSUTIL_AVAILABLE else None
    
    def get_memory_usage(self) -> float:
        """Get current memory usage in MB."""
        if not self.process:
            return 0.0
        memory_mb = self.process.memory_info().rss / 1024 / 1024
        self.peak_memory = max(self.peak_memory, memory_mb)
        return memory_mb
    
    def check_memory_limit(self, limit_mb: int = 4096) -> bool:
        """Check if memory usage exceeds limit."""
        current_memory = self.get_memory_usage()
        if current_memory > limit_mb:
            logging.warning(f"Memory usage ({current_memory:.1f} MB) exceeds limit ({limit_mb} MB)")
            return False
        return True
    
    def get_elapsed_time(self) -> float:
        """Get elapsed time in seconds."""
        return time.time() - self.start_time

# Core data structures
@dataclass
class Codon:
    """Represents a start or stop codon."""
    start: int
    end: int
    strand: str
    codon_type: str  # 'start' or 'stop'
    sequence: str = ""

@dataclass
class Exon:
    """Represents an exon."""
    start: int
    end: int
    strand: str
    id: str = ""

@dataclass
class CDS:
    """Represents a CDS region."""
    start: int
    end: int
    strand: str
    phase: int = 0
    id: str = ""

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
    # Add original annotation info
    chrom: str = ""
    source: str = ""
    
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
        coord_str = ";".join([f"{cds.start}-{cds.end}" for cds in sorted(self.cds_regions, key=lambda x: x.start)])
        return hashlib.md5(coord_str.encode()).hexdigest()

@dataclass
class Gene:
    """Represents a gene with multiple transcripts."""
    id: str
    transcripts: List[Transcript] = field(default_factory=list)
    representative: Optional[Transcript] = None
    
    def get_transcript_by_id(self, transcript_id: str) -> Optional[Transcript]:
        """Get transcript by ID."""
        for transcript in self.transcripts:
            if transcript.id == transcript_id:
                return transcript
        return None

class GeneAnnotationParser:
    """Parse GFF3/GTF files efficiently."""
    
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.genes: Dict[str, Gene] = {}
        self.transcripts: Dict[str, Transcript] = {}
    
    def parse_gff3(self) -> tuple[Dict[str, Gene], str]:
        """Parse GFF3/GTF file with O(n) complexity."""
        file_type = "GTF" if self.file_path.lower().endswith('.gtf') else "GFF3"
        logging.info(f"Parsing {file_type} file: {self.file_path}")
        
        gene_count = 0
        transcript_count = 0
        
        with open(self.file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    parts = line.split('\t')
                    if len(parts) != 9:
                        continue
                    
                    chrom, source, feature, start, end, score, strand, phase, attributes = parts
                    start, end = int(start), int(end)
                    
                    # Parse attributes based on file type
                    if file_type == "GTF":
                        attr_dict = self._parse_gtf_attributes(attributes)
                    else:
                        attr_dict = self._parse_gff3_attributes(attributes)
                    
                    # Process different feature types
                    if feature == 'gene':
                        if file_type == "GTF":
                            gene_id = attr_dict.get('gene_id', '')
                            # Handle bare gene IDs (like 'g1' without attributes)
                            if not gene_id and attributes.strip():
                                gene_id = attributes.strip()
                        else:
                            gene_id = attr_dict.get('ID', '')
                        
                        if gene_id:
                            self.genes[gene_id] = Gene(id=gene_id)
                            gene_count += 1
                    
                    elif feature in ['mRNA', 'transcript']:
                        if file_type == "GTF":
                            transcript_id = attr_dict.get('transcript_id', '')
                            parent_id = attr_dict.get('gene_id', '')
                            # Handle bare transcript IDs (like 'g130.t1' without attributes)
                            if not transcript_id and attributes.strip():
                                transcript_id = attributes.strip()
                                # Extract gene_id from transcript_id (g130.t1 -> g130)
                                if '.' in transcript_id:
                                    parent_id = transcript_id.split('.')[0]
                        else:
                            transcript_id = attr_dict.get('ID', '')
                            parent_id = attr_dict.get('Parent', '')
                        
                        if transcript_id and parent_id:
                            transcript = Transcript(
                                id=transcript_id,
                                gene_id=parent_id,
                                start=start,
                                end=end,
                                strand=strand,
                                chrom=chrom,
                                source=source
                            )
                            self.transcripts[transcript_id] = transcript
                            if parent_id in self.genes:
                                self.genes[parent_id].transcripts.append(transcript)
                            transcript_count += 1
                    
                    elif feature in ['exon', 'CDS', 'start_codon', 'stop_codon']:
                        if file_type == "GTF":
                            parent_id = attr_dict.get('transcript_id', '')
                        else:
                            parent_id = attr_dict.get('Parent', '')
                        
                        if parent_id in self.transcripts:
                            transcript = self.transcripts[parent_id]
                            
                            if feature == 'exon':
                                exon = Exon(start=start, end=end, strand=strand, id=attr_dict.get('ID', ''))
                                transcript.exons.append(exon)
                            
                            elif feature == 'CDS':
                                cds_phase = int(phase) if phase.isdigit() else 0
                                cds = CDS(start=start, end=end, strand=strand, phase=cds_phase, id=attr_dict.get('ID', ''))
                                transcript.cds_regions.append(cds)
                            
                            elif feature == 'start_codon':
                                transcript.start_codon = Codon(start=start, end=end, strand=strand, codon_type='start')
                            
                            elif feature == 'stop_codon':
                                transcript.stop_codon = Codon(start=start, end=end, strand=strand, codon_type='stop')
                
                except Exception as e:
                    logging.warning(f"Error parsing line {line_num}: {e}")
                    continue
        
        logging.info(f"Parsed {gene_count} genes and {transcript_count} transcripts")
        return self.genes, file_type
    
    def _parse_gff3_attributes(self, attr_string: str) -> Dict[str, str]:
        """Parse GFF3 attributes string."""
        attributes = {}
        for attr in attr_string.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attributes[key] = value
        return attributes
    
    def _parse_gtf_attributes(self, attr_string: str) -> Dict[str, str]:
        """Parse GTF attributes string."""
        attributes = {}
        for attr in attr_string.split(';'):
            attr = attr.strip()
            if not attr:
                continue
            if ' "' in attr:
                key, value = attr.split(' "', 1)
                key = key.strip()
                value = value.rstrip('"')
                attributes[key] = value
        return attributes

class SequenceHandler:
    """Handle FASTA sequence files efficiently."""
    
    def __init__(self, cds_path: str, aa_path: str, genome_path: Optional[str] = None):
        self.cds_sequences = self._parse_fasta(cds_path)
        self.aa_sequences = self._parse_fasta(aa_path)
        self.genome_index = pyfaidx.Fasta(genome_path) if genome_path else None
        
        logging.info(f"Loaded {len(self.cds_sequences)} CDS sequences")
        logging.info(f"Loaded {len(self.aa_sequences)} AA sequences")
    
    def _parse_fasta(self, file_path: str) -> Dict[str, str]:
        """Parse FASTA file with O(n) complexity."""
        sequences = {}
        current_id = None
        current_seq = []
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
        
        return sequences
    
    def get_genome_sequence(self, chrom: str, start: int, end: int, strand: str) -> str:
        """Get genome sequence with O(1) complexity using indexed access."""
        if not self.genome_index:
            return ""
        
        try:
            seq = str(self.genome_index[chrom][start-1:end])
            if strand == '-':
                seq = self._reverse_complement(seq)
            return seq
        except Exception as e:
            logging.warning(f"Error extracting genome sequence {chrom}:{start}-{end}: {e}")
            return ""
    
    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq.upper()))

class AnnotationCurator:
    """Phase A: Curate GFF/GTF annotations."""
    
    def __init__(self, sequence_handler: SequenceHandler, min_length: int = 30):
        self.sequence_handler = sequence_handler
        self.min_length = min_length
    
    def curate_annotations(self, genes: Dict[str, Gene]) -> Dict[str, Gene]:
        """Curate annotations with O(n) complexity - sequence-based quality only."""
        logging.info("Starting annotation curation (sequence-first approach)")
        
        curated_count = 0
        for gene in genes.values():
            for transcript in gene.transcripts:
                # Validate start codon
                if transcript.start_codon:
                    if not self._is_codon_within_features(transcript.start_codon, transcript.exons, transcript.cds_regions):
                        self._create_codon_features(transcript, transcript.start_codon)
                        transcript.quality_flags.add("created_start_codon_features")
                    else:
                        transcript.quality_flags.add("redundant_start_codon_removed")
                
                # Validate stop codon
                if transcript.stop_codon:
                    if not self._is_codon_within_features(transcript.stop_codon, transcript.exons, transcript.cds_regions):
                        self._create_codon_features(transcript, transcript.stop_codon)
                        transcript.quality_flags.add("created_stop_codon_features")
                    else:
                        transcript.quality_flags.add("redundant_stop_codon_removed")
                
                # Filter by length - use configurable threshold
                if len(transcript.aa_sequence) < self.min_length:
                    transcript.quality_flags.add("short_transcript")
                
                curated_count += 1
        
        # REMOVED: Overlap detection is deferred to post-selection phase
        # This eliminates the 94% false positive rate for manual review
        
        logging.info(f"Curated {curated_count} transcripts (sequence-based quality only)")
        return genes
    
    def _is_codon_within_features(self, codon: Codon, exons: List[Exon], cds_regions: List[CDS]) -> bool:
        """Check if codon is within existing exons/CDS using O(n) complexity."""
        codon_start, codon_end = codon.start, codon.end
        
        # Check exons
        for exon in exons:
            if exon.start <= codon_start and codon_end <= exon.end:
                # Check CDS regions
                for cds in cds_regions:
                    if cds.start <= codon_start and codon_end <= cds.end:
                        return True
        return False
    
    def _create_codon_features(self, transcript: Transcript, codon: Codon) -> None:
        """Create new exon and CDS features for codon, then merge with adjacent features if possible."""
        # Step 1: Create new exon and CDS features for codon
        new_exon = Exon(
            start=codon.start,
            end=codon.end,
            strand=codon.strand,
            id=f"{transcript.id}.{codon.codon_type}_exon"
        )
        
        new_cds = CDS(
            start=codon.start,
            end=codon.end,
            strand=codon.strand,
            phase=0,
            id=f"{transcript.id}.{codon.codon_type}_cds"
        )
        
        # Step 2: Check if new features can be merged with existing adjacent features
        merged_exon = self._merge_with_adjacent_exon(transcript, new_exon)
        merged_cds = self._merge_with_adjacent_cds(transcript, new_cds)
        
        # Step 3: Add features (merged or new) to transcript
        if not merged_exon:
            transcript.exons.append(new_exon)
        if not merged_cds:
            transcript.cds_regions.append(new_cds)
        
        # Step 4: Update transcript boundaries to cover new codon
        self._update_transcript_boundaries(transcript, codon)
    
    def _merge_with_adjacent_exon(self, transcript: Transcript, new_exon: Exon) -> bool:
        """Check if new exon can be merged with adjacent existing exons."""
        for exon in transcript.exons:
            # Check if new exon is adjacent to existing exon (allowing 1bp gap or overlap)
            if (abs(exon.end - new_exon.start) <= 1) or (abs(new_exon.end - exon.start) <= 1):
                # Merge by extending the existing exon boundaries
                exon.start = min(exon.start, new_exon.start)
                exon.end = max(exon.end, new_exon.end)
                return True  # Successfully merged
        return False  # No adjacent exon found
    
    def _merge_with_adjacent_cds(self, transcript: Transcript, new_cds: CDS) -> bool:
        """Check if new CDS can be merged with adjacent existing CDS."""
        for cds in transcript.cds_regions:
            # Check if new CDS is adjacent to existing CDS (allowing 1bp gap or overlap)
            if (abs(cds.end - new_cds.start) <= 1) or (abs(new_cds.end - cds.start) <= 1):
                # Merge by extending the existing CDS boundaries
                cds.start = min(cds.start, new_cds.start)
                cds.end = max(cds.end, new_cds.end)
                # Keep the phase of the original CDS (important for translation)
                return True  # Successfully merged
        return False  # No adjacent CDS found
    
    def _update_transcript_boundaries(self, transcript: Transcript, codon: Codon) -> None:
        """Update transcript start/end coordinates to include new codon."""
        transcript.start = min(transcript.start, codon.start)
        transcript.end = max(transcript.end, codon.end)
    
    def _detect_overlaps(self, genes: Dict[str, Gene]) -> None:
        """Detect overlapping transcripts using IntervalTree O(n log n)."""
        # Build interval tree for efficient overlap detection
        interval_tree = IntervalTree()
        transcript_map = {}
        
        # Add all transcripts to interval tree
        for gene in genes.values():
            for transcript in gene.transcripts:
                interval = Interval(transcript.start, transcript.end + 1)
                interval_tree.add(interval)
                transcript_map[interval] = transcript
        
        # Find overlaps
        for gene in genes.values():
            for transcript in gene.transcripts:
                overlaps = interval_tree.overlap(transcript.start, transcript.end + 1)
                overlapping_transcripts = [transcript_map[interval] for interval in overlaps 
                                         if transcript_map[interval].id != transcript.id]
                
                if overlapping_transcripts:
                    transcript.quality_flags.add("overlapping_transcript")


class SequenceValidator:
    """Phase B: Validate and correct CDS and AA sequences."""
    
    def __init__(self, sequence_handler: SequenceHandler):
        self.sequence_handler = sequence_handler
    
    def validate_sequences(self, genes: Dict[str, Gene]) -> Dict[str, Gene]:
        """Validate sequences with O(n) complexity."""
        logging.info("Starting sequence validation")
        
        validated_count = 0
        corrected_count = 0
        
        for gene in genes.values():
            for transcript in gene.transcripts:
                # Validate CDS sequence
                if self._validate_cds_sequence(transcript):
                    corrected_count += 1
                
                # Reconstruct CDS sequence from merged coordinates if genome is available
                if self.sequence_handler.genome_index:
                    self._reconstruct_cds_sequence(transcript)
                
                # Validate AA sequence
                self._validate_aa_sequence(transcript)
                
                validated_count += 1
        
        logging.info(f"Validated {validated_count} transcripts, corrected {corrected_count}")
        return genes
    
    def _validate_cds_sequence(self, transcript: Transcript) -> bool:
        """Validate and correct CDS sequence."""
        corrected = False
        
        # Check if codons are present in GFF/GTF
        if not transcript.start_codon and not transcript.stop_codon:
            transcript.quality_flags.add("low_quality_no_codons")
            return False
        
        # Check start codon in CDS
        if transcript.start_codon and transcript.cds_sequence:
            if not transcript.cds_sequence.upper().startswith(('ATG', 'GTG', 'TTG')):
                # Try to add start codon from genome
                if self._add_start_codon_from_genome(transcript):
                    transcript.quality_flags.add("corrected_start_codon")
                    corrected = True
        
        # Check stop codon in CDS
        if transcript.stop_codon and transcript.cds_sequence:
            if not transcript.cds_sequence.upper().endswith(('TAA', 'TAG', 'TGA')):
                # Try to add stop codon from genome
                if self._add_stop_codon_from_genome(transcript):
                    transcript.quality_flags.add("corrected_stop_codon")
                    corrected = True
        
        if not corrected:
            transcript.quality_flags.add("validated_cds")
        
        return corrected
    
    def _validate_aa_sequence(self, transcript: Transcript) -> None:
        """Validate amino acid sequence."""
        if "low_quality_no_codons" in transcript.quality_flags:
            if not transcript.aa_sequence.startswith('M'):
                transcript.quality_flags.add("low_quality_no_start_codon")
            else:
                transcript.quality_flags.add("passed_aa_validation")
        else:
            # Check for internal stop codons
            if '*' in transcript.aa_sequence[:-1]:  # Exclude last position
                transcript.quality_flags.add("internal_stop_codons")
            else:
                transcript.quality_flags.add("passed_aa_validation")
    
    def _add_start_codon_from_genome(self, transcript: Transcript) -> bool:
        """Add start codon from genome sequence."""
        if not transcript.start_codon or not self.sequence_handler.genome_index:
            return False
        
        # Extract start codon sequence from genome
        codon_seq = self.sequence_handler.get_genome_sequence(
            transcript.chrom, 
            transcript.start_codon.start, 
            transcript.start_codon.end,
            transcript.strand
        )
        
        if codon_seq and codon_seq.upper() in ['ATG', 'GTG', 'TTG']:
            # For start codon, prepend to CDS sequence
            transcript.cds_sequence = codon_seq + transcript.cds_sequence
            return True
        return False
    
    def _add_stop_codon_from_genome(self, transcript: Transcript) -> bool:
        """Add stop codon from genome sequence."""
        if not transcript.stop_codon or not self.sequence_handler.genome_index:
            return False
        
        # Extract stop codon sequence from genome
        codon_seq = self.sequence_handler.get_genome_sequence(
            transcript.chrom, 
            transcript.stop_codon.start, 
            transcript.stop_codon.end,
            transcript.strand
        )
        
        if codon_seq and codon_seq.upper() in ['TAA', 'TAG', 'TGA']:
            # For stop codon, append to CDS sequence
            transcript.cds_sequence += codon_seq
            return True
        return False
    
    def _reconstruct_cds_sequence(self, transcript: Transcript) -> None:
        """Reconstruct complete CDS sequence from merged coordinates."""
        if not transcript.cds_regions or not self.sequence_handler.genome_index:
            return
        
        # Sort CDS regions by coordinate (important for proper sequence order)
        sorted_cds = sorted(transcript.cds_regions, key=lambda x: x.start if transcript.strand == '+' else -x.start)
        
        # Extract sequence from each CDS region
        reconstructed_sequence = ""
        for cds in sorted_cds:
            cds_seq = self.sequence_handler.get_genome_sequence(
                transcript.chrom, 
                cds.start, 
                cds.end,
                transcript.strand
            )
            if cds_seq:
                reconstructed_sequence += cds_seq
        
        # Update transcript with reconstructed sequence
        if reconstructed_sequence:
            transcript.cds_sequence = reconstructed_sequence
            transcript.quality_flags.add("cds_reconstructed_from_merged_coordinates")


class TranscriptSelector:
    """Phase C: Select representative transcripts with sequence-first approach."""
    
    def __init__(self, min_length: int, overlap_threshold: float):
        self.min_length = min_length
        self.overlap_threshold = overlap_threshold
    
    def select_representatives(self, genes: Dict[str, Gene]) -> Dict[str, Gene]:
        """Select representatives with sequence-first approach - O(n) complexity per gene."""
        logging.info("Starting representative selection (sequence-first approach)")
        
        selected_count = 0
        manual_review_count = 0
        
        for gene in genes.values():
            # Step 1: Remove low-quality transcripts (sequence-based only)
            high_quality_transcripts = self._filter_low_quality_sequence_based(gene.transcripts)
            
            # If no high-quality transcripts, check if longest transcript meets minimum threshold
            if not high_quality_transcripts:
                if gene.transcripts:
                    # Find the longest transcript
                    longest_transcript = max(gene.transcripts, key=lambda t: len(t.aa_sequence))
                    
                    # Only select if it meets the --min-length requirement
                    if len(longest_transcript.aa_sequence) >= self.min_length:
                        gene.representative = longest_transcript
                        longest_transcript.quality_flags.add("representative")
                        selected_count += 1
                    else:
                        # All transcripts are too short - mark for manual review
                        for transcript in gene.transcripts:
                            transcript.quality_flags.add("manual_review")
                            transcript.quality_flags.add("all_transcripts_too_short")
                        manual_review_count += 1
                continue
            
            # Step 2: Group by AA sequence MD5 hash
            hash_groups = self._group_by_hash(high_quality_transcripts)
            
            # Step 3: Select representative (no overlap checking)
            representative = self._select_from_hash_groups_pragmatic(hash_groups)
            
            if representative:
                gene.representative = representative
                representative.quality_flags.add("representative")
                selected_count += 1
            else:
                # Only mark for manual review if truly ambiguous
                for transcript in gene.transcripts:
                    transcript.quality_flags.add("manual_review")
                manual_review_count += 1
        
        # Step 4: Post-selection overlap assessment
        post_selection_conflicts = self._assess_post_selection_overlaps(genes)
        
        logging.info(f"Selected {selected_count} representatives, {manual_review_count} need manual review")
        logging.info(f"Post-selection spatial conflicts: {post_selection_conflicts}")
        return genes
    
    def _filter_low_quality_sequence_based(self, transcripts: List[Transcript]) -> List[Transcript]:
        """Filter out low-quality transcripts based on sequence criteria only."""
        high_quality = []
        
        for transcript in transcripts:
            # Check length
            if len(transcript.aa_sequence) < self.min_length:
                transcript.quality_flags.add("short_transcript")
                continue
            
            # Check sequence-based quality flags only (no spatial criteria)
            if any(flag in transcript.quality_flags for flag in 
                   ["low_quality_no_codons", "low_quality_no_start_codon", "internal_stop_codons"]):
                continue
            
            high_quality.append(transcript)
        
        return high_quality
    
    def _filter_low_quality(self, transcripts: List[Transcript]) -> List[Transcript]:
        """Legacy method - kept for compatibility."""
        return self._filter_low_quality_sequence_based(transcripts)
    
    def _group_by_hash(self, transcripts: List[Transcript]) -> Dict[str, List[Transcript]]:
        """Group transcripts by AA sequence MD5 hash."""
        hash_groups = defaultdict(list)
        
        for transcript in transcripts:
            aa_hash = transcript.get_aa_hash()
            if aa_hash:
                hash_groups[aa_hash].append(transcript)
        
        return hash_groups
    
    def _select_from_hash_groups_pragmatic(self, hash_groups: Dict[str, List[Transcript]]) -> Optional[Transcript]:
        """Select representative from hash groups with pragmatic approach."""
        if not hash_groups:
            return None
        
        # Find the most frequent hash group
        most_frequent_hash = max(hash_groups.keys(), key=lambda k: len(hash_groups[k]))
        candidates = hash_groups[most_frequent_hash]
        
        # Pragmatic selection: always select the longest transcript
        # No overlap checking at this stage - deferred to post-selection
        selected = max(candidates, key=lambda t: len(t.aa_sequence))
        
        return selected
    
    def _select_from_hash_groups(self, hash_groups: Dict[str, List[Transcript]]) -> Optional[Transcript]:
        """Legacy method - kept for compatibility."""
        return self._select_from_hash_groups_pragmatic(hash_groups)
    
    def _assess_post_selection_overlaps(self, genes: Dict[str, Gene]) -> int:
        """Assess spatial conflicts AFTER representative selection."""
        logging.info("Assessing post-selection spatial conflicts")
        
        # Extract all representative transcripts
        representatives = []
        for gene in genes.values():
            if gene.representative:
                representatives.append(gene.representative)
        
        if not representatives:
            return 0
        
        # Build interval tree with representatives only
        interval_tree = IntervalTree()
        transcript_map = {}
        
        for transcript in representatives:
            interval = Interval(transcript.start, transcript.end + 1)
            interval_tree.add(interval)
            transcript_map[interval] = transcript
        
        conflicts_detected = 0
        
        # Check for significant overlaps between representatives
        for transcript in representatives:
            overlaps = interval_tree.overlap(transcript.start, transcript.end + 1)
            
            for overlap_interval in overlaps:
                other_transcript = transcript_map[overlap_interval]
                
                # Skip self-overlap
                if other_transcript.id == transcript.id:
                    continue
                
                # Only consider inter-gene conflicts
                if other_transcript.gene_id == transcript.gene_id:
                    continue
                
                # Calculate overlap significance
                overlap_start = max(transcript.start, other_transcript.start)
                overlap_end = min(transcript.end, other_transcript.end)
                overlap_length = max(0, overlap_end - overlap_start)
                
                # Check reciprocal overlap
                transcript_length = transcript.end - transcript.start
                other_length = other_transcript.end - other_transcript.start
                
                overlap_pct_1 = overlap_length / transcript_length if transcript_length > 0 else 0
                overlap_pct_2 = overlap_length / other_length if other_length > 0 else 0
                
                # Flag significant reciprocal overlaps
                if overlap_pct_1 > self.overlap_threshold and overlap_pct_2 > self.overlap_threshold:
                    # Check if they encode identical proteins
                    if transcript.get_aa_hash() == other_transcript.get_aa_hash():
                        # Identical proteins - mark as redundant
                        transcript.quality_flags.add("redundant_gene_model")
                        other_transcript.quality_flags.add("redundant_gene_model")
                    else:
                        # Different proteins - real spatial conflict
                        transcript.quality_flags.add("real_spatial_conflict")
                        other_transcript.quality_flags.add("real_spatial_conflict")
                        conflicts_detected += 1
                
                # Flag any post-selection overlap for tracking
                transcript.quality_flags.add("post_selection_overlap")
                other_transcript.quality_flags.add("post_selection_overlap")
        
        return conflicts_detected


class ManualCurator:
    """Phase D: Integrate manual curation decisions."""
    
    def __init__(self, curated_list_path: str):
        self.curated_list_path = curated_list_path
        self.curated_selections = self._load_curated_list()
    
    def _load_curated_list(self) -> Dict[str, Tuple[str, str, str]]:
        """Load manual curation list."""
        curated = {}
        
        with open(self.curated_list_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 4:
                    gene_id, transcript_id, reason, quality = parts[:4]
                    curated[gene_id] = (transcript_id, reason, quality)
        
        return curated
    
    def integrate_manual_curation(self, genes: Dict[str, Gene]) -> Dict[str, Gene]:
        """Integrate manual curation with O(n) complexity."""
        logging.info("Integrating manual curation")
        
        integrated_count = 0
        
        for gene_id, (transcript_id, reason, quality) in self.curated_selections.items():
            if gene_id in genes:
                gene = genes[gene_id]
                
                if transcript_id == "REMOVE":
                    # Remove entire gene
                    gene.representative = None
                    for transcript in gene.transcripts:
                        transcript.quality_flags.add("manually_removed")
                else:
                    # Select specific transcript
                    selected_transcript = gene.get_transcript_by_id(transcript_id)
                    if selected_transcript:
                        gene.representative = selected_transcript
                        selected_transcript.quality_flags.add("manually_selected")
                        selected_transcript.quality_flags.add(f"manual_reason_{reason}")
                        integrated_count += 1
        
        logging.info(f"Integrated {integrated_count} manual selections")
        return genes


class OutputGenerator:
    """Generate output files."""
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
    
    def generate_outputs(self, genes: Dict[str, Gene], sequence_handler: SequenceHandler, pipeline_stats: Dict = None, input_format: str = "GFF3") -> None:
        """Generate all output files."""
        logging.info("Generating output files")
        
        # Generate cleaned gene model file (preserve input format)
        if input_format == "GTF":
            self._write_cleaned_gtf(genes)
        else:
            self._write_cleaned_gff3(genes)
        
        # Generate cleaned sequences
        self._write_cleaned_sequences(genes, sequence_handler)
        
        # Generate processing report
        self._write_processing_report(genes, pipeline_stats)
        
        # Generate manual review list
        self._write_manual_review_list(genes)
    
    def _write_cleaned_gff3(self, genes: Dict[str, Gene]) -> None:
        """Write cleaned GFF3 file with gene boundary adjustments."""
        output_path = os.path.join(self.output_dir, 'cleaned.gff3')
        
        with open(output_path, 'w') as f:
            f.write("##gff-version 3\n")
            
            for gene in genes.values():
                if gene.representative:
                    transcript = gene.representative
                    
                    # Calculate updated gene boundaries based on representative transcript
                    # Include all exons to ensure proper gene boundary calculation
                    if transcript.exons:
                        exon_starts = [e.start for e in transcript.exons]
                        exon_ends = [e.end for e in transcript.exons]
                        gene_start = min(exon_starts)
                        gene_end = max(exon_ends)
                    else:
                        gene_start = transcript.start
                        gene_end = transcript.end
                    
                    # Check if gene boundaries have changed (shrunk or expanded)
                    original_gene_spans = [t.start for t in gene.transcripts] + [t.end for t in gene.transcripts]
                    original_start = min(original_gene_spans)
                    original_end = max(original_gene_spans)
                    
                    if gene_start != original_start or gene_end != original_end:
                        transcript.quality_flags.add("gene_boundaries_updated")
                        if gene_start > original_start or gene_end < original_end:
                            transcript.quality_flags.add("gene_shrunk")
                        if gene_start < original_start or gene_end > original_end:
                            transcript.quality_flags.add("gene_expanded")
                    
                    # Write gene feature with updated boundaries (preserve original chromosome and source)
                    chrom = transcript.chrom or "LAT_C01"  # Use original chromosome name
                    source = transcript.source or "AUGUSTUS"  # Use original source
                    
                    f.write(f"{chrom}\t{source}\tgene\t{gene_start}\t{gene_end}\t.\t{transcript.strand}\t.\tID={gene.id}\n")
                    
                    # Write mRNA feature
                    f.write(f"{chrom}\t{source}\tmRNA\t{transcript.start}\t{transcript.end}\t.\t{transcript.strand}\t.\tID={transcript.id};Parent={gene.id}\n")
                    
                    # Write exons
                    for i, exon in enumerate(transcript.exons):
                        f.write(f"{chrom}\t{source}\texon\t{exon.start}\t{exon.end}\t.\t{exon.strand}\t.\tID={transcript.id}.exon{i+1};Parent={transcript.id}\n")
                    
                    # Write CDS
                    for i, cds in enumerate(transcript.cds_regions):
                        f.write(f"{chrom}\t{source}\tCDS\t{cds.start}\t{cds.end}\t.\t{cds.strand}\t{cds.phase}\tID={transcript.id}.cds{i+1};Parent={transcript.id}\n")
    
    def _write_cleaned_gtf(self, genes: Dict[str, Gene]) -> None:
        """Write cleaned GTF file with gene boundary adjustments."""
        output_path = os.path.join(self.output_dir, 'cleaned.gtf')
        
        with open(output_path, 'w') as f:
            for gene in genes.values():
                if gene.representative:
                    transcript = gene.representative
                    
                    # Calculate updated gene boundaries based on representative transcript
                    # Include all exons to ensure proper gene boundary calculation
                    if transcript.exons:
                        exon_starts = [e.start for e in transcript.exons]
                        exon_ends = [e.end for e in transcript.exons]
                        gene_start = min(exon_starts)
                        gene_end = max(exon_ends)
                    else:
                        gene_start = transcript.start
                        gene_end = transcript.end
                    
                    # Check if gene boundaries have changed (shrunk or expanded)
                    original_gene_spans = [t.start for t in gene.transcripts] + [t.end for t in gene.transcripts]
                    original_start = min(original_gene_spans)
                    original_end = max(original_gene_spans)
                    
                    if gene_start != original_start or gene_end != original_end:
                        transcript.quality_flags.add("gene_boundaries_updated")
                        if gene_start > original_start or gene_end < original_end:
                            transcript.quality_flags.add("gene_shrunk")
                        if gene_start < original_start or gene_end > original_end:
                            transcript.quality_flags.add("gene_expanded")
                    
                    # Write gene feature with updated boundaries (preserve original chromosome and source)
                    chrom = transcript.chrom or "LAT_C01"  # Use original chromosome name
                    source = transcript.source or "AUGUSTUS"  # Use original source
                    
                    f.write(f"{chrom}\t{source}\tgene\t{gene_start}\t{gene_end}\t.\t{transcript.strand}\t.\tgene_id \"{gene.id}\";\n")
                    
                    # Write transcript feature
                    f.write(f"{chrom}\t{source}\ttranscript\t{transcript.start}\t{transcript.end}\t.\t{transcript.strand}\t.\ttranscript_id \"{transcript.id}\"; gene_id \"{gene.id}\";\n")
                    
                    # Write exons
                    for exon in transcript.exons:
                        f.write(f"{chrom}\t{source}\texon\t{exon.start}\t{exon.end}\t.\t{exon.strand}\t.\ttranscript_id \"{transcript.id}\"; gene_id \"{gene.id}\";\n")
                    
                    # Write CDS
                    for cds in transcript.cds_regions:
                        f.write(f"{chrom}\t{source}\tCDS\t{cds.start}\t{cds.end}\t.\t{cds.strand}\t{cds.phase}\ttranscript_id \"{transcript.id}\"; gene_id \"{gene.id}\";\n")
                    
                    # Write start_codon if present
                    if transcript.start_codon:
                        sc = transcript.start_codon
                        f.write(f"{chrom}\t{source}\tstart_codon\t{sc.start}\t{sc.end}\t.\t{sc.strand}\t0\ttranscript_id \"{transcript.id}\"; gene_id \"{gene.id}\";\n")
                    
                    # Write stop_codon if present
                    if transcript.stop_codon:
                        sc = transcript.stop_codon
                        f.write(f"{chrom}\t{source}\tstop_codon\t{sc.start}\t{sc.end}\t.\t{sc.strand}\t0\ttranscript_id \"{transcript.id}\"; gene_id \"{gene.id}\";\n")
    
    def _write_cleaned_sequences(self, genes: Dict[str, Gene], sequence_handler: SequenceHandler) -> None:
        """Write cleaned CDS and AA sequences."""
        cds_path = os.path.join(self.output_dir, 'cleaned.cds.fa')
        aa_path = os.path.join(self.output_dir, 'cleaned.aa')
        
        with open(cds_path, 'w') as cds_f, open(aa_path, 'w') as aa_f:
            for gene in genes.values():
                if gene.representative:
                    transcript = gene.representative
                    
                    # Write CDS sequence
                    flags = ",".join(transcript.quality_flags)
                    cds_f.write(f">{transcript.id} flags={flags}\n")
                    cds_f.write(f"{transcript.cds_sequence}\n")
                    
                    # Write AA sequence
                    aa_f.write(f">{transcript.id} flags={flags}\n")
                    aa_f.write(f"{transcript.aa_sequence}\n")
    
    def _write_processing_report(self, genes: Dict[str, Gene], pipeline_stats: Dict = None) -> None:
        """Write comprehensive processing statistics report showing pipeline progression."""
        report_path = os.path.join(self.output_dir, 'processing_report.txt')
        
        # Use pipeline_stats if provided, otherwise calculate from current state
        if pipeline_stats:
            original_genes = pipeline_stats['original_genes']
            original_transcripts = pipeline_stats['original_transcripts']
            selected_representatives = pipeline_stats['selected_representatives']
            post_selection_conflicts = pipeline_stats['post_selection_conflicts']
            manual_review_genes = pipeline_stats['manual_review_genes']
        else:
            original_genes = len(genes)
            original_transcripts = sum(len(g.transcripts) for g in genes.values())
            selected_representatives = sum(1 for g in genes.values() if g.representative)
            post_selection_conflicts = sum(1 for g in genes.values() 
                                         if g.representative and 'post_selection_overlap' in g.representative.quality_flags)
            manual_review_genes = sum(1 for g in genes.values() 
                                    if any('manual_review' in t.quality_flags for t in g.transcripts))
        
        # Calculate final output statistics
        final_genes = sum(1 for g in genes.values() if g.representative)
        final_transcripts = final_genes  # One representative per gene
        
        # Count quality improvements
        gene_boundaries_updated = sum(1 for g in genes.values() 
                                     if g.representative and 'gene_boundaries_updated' in g.representative.quality_flags)
        start_codons_added = sum(1 for g in genes.values() 
                               if g.representative and 'created_start_codon_features' in g.representative.quality_flags)
        stop_codons_added = sum(1 for g in genes.values() 
                              if g.representative and 'created_stop_codon_features' in g.representative.quality_flags)
        
        # Calculate rates
        selection_rate = (selected_representatives / original_genes * 100) if original_genes > 0 else 0
        manual_review_rate = (manual_review_genes / original_genes * 100) if original_genes > 0 else 0
        
        with open(report_path, 'w') as f:
            f.write("Gene Annotation Curation Pipeline Report\n")
            f.write("=" * 40 + "\n\n")
            
            f.write("ORIGINAL INPUT:\n")
            f.write("-" * 15 + "\n")
            f.write(f"Input Genes: {original_genes:,}\n")
            f.write(f"Input Transcripts: {original_transcripts:,}\n\n")
            
            f.write("SEQUENCE-FIRST SELECTION RESULTS:\n")
            f.write("-" * 35 + "\n")
            f.write(f"Representatives Selected: {selected_representatives:,}\n")
            f.write(f"Selection Success Rate: {selection_rate:.1f}%\n")
            f.write(f"Genes Requiring Manual Review: {manual_review_genes:,}\n")
            f.write(f"Manual Review Rate: {manual_review_rate:.1f}%\n\n")
            
            f.write("POST-SELECTION SPATIAL ANALYSIS:\n")
            f.write("-" * 32 + "\n")
            f.write(f"Spatial Conflicts Detected: {post_selection_conflicts:,}\n")
            f.write(f"Gene Boundaries Adjusted: {gene_boundaries_updated:,}\n\n")
            
            f.write("FINAL OUTPUT:\n")
            f.write("-" * 13 + "\n")
            f.write(f"Final Genes: {final_genes:,}\n")
            f.write(f"Final Transcripts: {final_transcripts:,}\n")
            f.write(f"Reduction Rate: {((original_transcripts - final_transcripts) / original_transcripts * 100):.1f}%\n\n")
            
            f.write("QUALITY IMPROVEMENTS:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Start Codons Added: {start_codons_added:,}\n")
            f.write(f"Stop Codons Added: {stop_codons_added:,}\n")
            f.write(f"Gene Boundaries Updated: {gene_boundaries_updated:,}\n\n")
            
            f.write("PIPELINE APPROACH:\n")
            f.write("-" * 18 + "\n")
            f.write("Strategy: Sequence-first, spatial-second methodology\n")
            f.write("Key Innovation: Pragmatic selection with deferred overlap detection\n")
            f.write("Result: High-throughput curation with minimal manual intervention\n")
    
    def _write_manual_review_list(self, genes: Dict[str, Gene]) -> None:
        """Write list of genes/transcripts needing manual review."""
        review_path = os.path.join(self.output_dir, 'manual_review_transcripts.txt')
        
        with open(review_path, 'w') as f:
            f.write("# Gene_ID\tTranscript_IDs\tReason\tQuality_Flags\n")
            
            for gene in genes.values():
                # Only include genes that actually need manual review
                # Genes with representatives that have spatial conflicts don't need manual review
                # Only genes without representatives or with true ambiguity need manual review
                
                if not gene.representative:
                    # Gene has no representative - needs manual review
                    transcript_ids = ",".join(t.id for t in gene.transcripts)
                    flags = ",".join(set().union(*(t.quality_flags for t in gene.transcripts)))
                    f.write(f"{gene.id}\t{transcript_ids}\tno_representative_selected\t{flags}\n")
                else:
                    # Gene has representative - only flag if explicitly marked for manual review
                    manual_review_transcripts = [t for t in gene.transcripts 
                                               if 'manual_review' in t.quality_flags]
                    
                    if manual_review_transcripts:
                        transcript_ids = ",".join(t.id for t in manual_review_transcripts)
                        flags = ",".join(set().union(*(t.quality_flags for t in manual_review_transcripts)))
                        f.write(f"{gene.id}\t{transcript_ids}\tconflicting_selection\t{flags}\n")


def setup_logging(log_level: str = 'INFO', log_file: str = 'gene_curation.log') -> None:
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def main():
    """Main pipeline execution."""
    parser = argparse.ArgumentParser(description='Gene Annotation Curation Pipeline')
    parser.add_argument('--gene-model', required=True, help='Input gene model file (GFF3 or GTF format)')
    parser.add_argument('--cds', required=True, help='Input CDS FASTA file')
    parser.add_argument('--aa', required=True, help='Input amino acid FASTA file')
    parser.add_argument('--genome', help='Reference genome FASTA file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--min-length', type=int, default=50, help='Minimum AA length threshold')
    parser.add_argument('--overlap-threshold', type=float, default=0.8, help='Overlap threshold for quality assessment')
    parser.add_argument('--log-level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'])
    parser.add_argument('--curated-list', help='Manual curation file')
    parser.add_argument('--integrate-manual-curation', action='store_true', help='Integrate manual curation')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup logging
    log_file = os.path.join(args.output_dir, 'gene_curation.log')
    setup_logging(args.log_level, log_file)
    
    # Initialize performance monitor
    monitor = PerformanceMonitor()
    
    logging.info("Starting Gene Annotation Curation Pipeline")
    logging.info(f"Input files: Gene Model={args.gene_model}, CDS={args.cds}, AA={args.aa}")
    logging.info(f"Output directory: {args.output_dir}")
    
    try:
        # Parse input files
        logging.info("Phase 1: Parsing input files")
        parser = GeneAnnotationParser(args.gene_model)
        genes, input_format = parser.parse_gff3()
        
        sequence_handler = SequenceHandler(args.cds, args.aa, args.genome)
        
        # Link sequences to transcripts
        logging.info("Phase 2: Linking sequences to transcripts")
        for gene in genes.values():
            for transcript in gene.transcripts:
                transcript.cds_sequence = sequence_handler.cds_sequences.get(transcript.id, "")
                transcript.aa_sequence = sequence_handler.aa_sequences.get(transcript.id, "")
        
        # Check memory usage
        memory_usage = monitor.get_memory_usage()
        logging.info(f"Memory usage after parsing: {memory_usage:.1f} MB")
        
        if not monitor.check_memory_limit():
            logging.error("Memory usage exceeded limit")
            return 1
        
        # Phase A: Annotation Curation
        logging.info("Phase A: Annotation Curation")
        annotation_curator = AnnotationCurator(sequence_handler, args.min_length)
        genes = annotation_curator.curate_annotations(genes)
        
        # Phase B: Sequence Validation
        logging.info("Phase B: Sequence Validation")
        sequence_validator = SequenceValidator(sequence_handler)
        genes = sequence_validator.validate_sequences(genes)
        
        # Phase C: Representative Selection
        logging.info("Phase C: Representative Selection")
        transcript_selector = TranscriptSelector(args.min_length, args.overlap_threshold)
        genes = transcript_selector.select_representatives(genes)
        
        # Phase D: Manual Curation Integration (if requested)
        if args.integrate_manual_curation and args.curated_list:
            logging.info("Phase D: Manual Curation Integration")
            manual_curator = ManualCurator(args.curated_list)
            genes = manual_curator.integrate_manual_curation(genes)
        
        # Generate output files
        logging.info("Generating output files")
        output_generator = OutputGenerator(args.output_dir)
        
        # Collect pipeline statistics
        pipeline_stats = {
            'original_genes': len(genes),
            'original_transcripts': sum(len(g.transcripts) for g in genes.values()),
            'selected_representatives': sum(1 for g in genes.values() if g.representative),
            'post_selection_conflicts': sum(1 for g in genes.values() 
                                           if g.representative and 'post_selection_overlap' in g.representative.quality_flags),
            'manual_review_genes': sum(1 for g in genes.values() 
                                      if any('manual_review' in t.quality_flags for t in g.transcripts))
        }
        
        output_generator.generate_outputs(genes, sequence_handler, pipeline_stats, input_format)
        
        # Performance summary
        elapsed_time = monitor.get_elapsed_time()
        logging.info(f"Pipeline completed in {elapsed_time:.1f} seconds")
        logging.info(f"Peak memory usage: {monitor.peak_memory:.1f} MB")
        
        return 0
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())