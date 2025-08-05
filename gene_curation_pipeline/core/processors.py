#!/usr/bin/env python3

"""
Processing classes for annotation curation, sequence validation, and transcript selection.
"""

import logging
import hashlib
from collections import defaultdict
from typing import Dict, List, Set, Optional, Tuple
from .data_structures import Gene, Transcript, Exon, CDS, Codon
from .exceptions import PipelineError, ValidationError, SequenceError

try:
    from intervaltree import IntervalTree, Interval
    INTERVALTREE_AVAILABLE = True
except ImportError:
    INTERVALTREE_AVAILABLE = False


class AnnotationCurator:
    """Phase A: Curate GFF/GTF annotations with codon integration."""
    
    def __init__(self, adjacency_gap_threshold: int = 1, genome=None):
        self.adjacency_gap_threshold = adjacency_gap_threshold
        self.genome = genome
    
    def curate_gene(self, gene: Gene) -> None:
        """Curate a single gene's annotations."""
        for transcript in gene.transcripts:
            # Ensure CDS regions have corresponding exons
            self._ensure_cds_have_exons(transcript)
            
            # Process start codon integration
            if transcript.start_codon:
                self._process_codon_integration(transcript, transcript.start_codon, 'start')
            
            # Process stop codon integration
            if transcript.stop_codon:
                self._process_codon_integration(transcript, transcript.stop_codon, 'stop')
    
    def _ensure_cds_have_exons(self, transcript: Transcript) -> None:
        """Ensure all CDS regions have corresponding exons."""
        if not transcript.cds_regions:
            return
        
        exons_created = 0
        for cds in transcript.cds_regions:
            # Check if this CDS region is covered by any existing exon
            is_covered = False
            for exon in transcript.exons:
                if exon.start <= cds.start and cds.end <= exon.end:
                    is_covered = True
                    break
            
            # If CDS is not covered by any exon, create a matching exon
            if not is_covered:
                new_exon = Exon(
                    start=cds.start,
                    end=cds.end,
                    strand=cds.strand,
                    id=f"{transcript.id}_exon_for_cds_{cds.start}_{cds.end}"
                )
                transcript.exons.append(new_exon)
                exons_created += 1
        
        if exons_created > 0:
            transcript.quality_flags.add("created_exon_for_cds")
    
    def _process_codon_integration(self, transcript: Transcript, codon: Codon, codon_type: str) -> None:
        """Process codon integration with adjacent feature merging."""
        # Check if codon is already covered by existing features
        if self._is_codon_within_features(codon, transcript.exons, transcript.cds_regions):
            transcript.quality_flags.add(f"redundant_{codon_type}_codon_removed")
            return
        
        # Check if codon can be merged with adjacent features
        if self._try_merge_codon_with_adjacent_features(transcript, codon):
            transcript.quality_flags.add(f"merged_{codon_type}_codon_with_adjacent")
        else:
            # Create new features for the codon
            self._create_codon_features(transcript, codon)
            transcript.quality_flags.add(f"created_{codon_type}_codon_features")
    
    def _is_codon_within_features(self, codon: Codon, exons: List[Exon], cds_regions: List[CDS]) -> bool:
        """Check if codon is within existing exons/CDS."""
        codon_start, codon_end = codon.start, codon.end
        
        # Check exons and CDS regions
        for exon in exons:
            if exon.start <= codon_start and codon_end <= exon.end:
                for cds in cds_regions:
                    if cds.start <= codon_start and codon_end <= cds.end:
                        return True
        return False
    
    def _try_merge_codon_with_adjacent_features(self, transcript: Transcript, codon: Codon) -> bool:
        """Try to merge codon with adjacent exons/CDS features."""
        merged = False
        
        # Try to merge with exons
        for exon in transcript.exons:
            if self._is_adjacent(codon, exon, self.adjacency_gap_threshold):
                # Extend exon to include codon
                exon.start = min(exon.start, codon.start)
                exon.end = max(exon.end, codon.end)
                merged = True
                break
        
        # Try to merge with CDS regions
        for cds in transcript.cds_regions:
            if self._is_adjacent(codon, cds, self.adjacency_gap_threshold):
                # Extend CDS to include codon
                cds.start = min(cds.start, codon.start)
                cds.end = max(cds.end, codon.end)
                merged = True
                break
        
        return merged
    
    def _is_adjacent(self, codon: Codon, feature, gap_threshold: int) -> bool:
        """Check if codon is adjacent to a feature within gap threshold."""
        return (abs(codon.end - feature.start) <= gap_threshold or 
                abs(feature.end - codon.start) <= gap_threshold)
    
    def _create_codon_features(self, transcript: Transcript, codon: Codon) -> None:
        """Create new exon and CDS features for codon."""
        # Create new exon
        new_exon = Exon(
            start=codon.start,
            end=codon.end,
            strand=codon.strand,
            id=f"{transcript.id}_{codon.codon_type}_exon"
        )
        transcript.exons.append(new_exon)
        
        # Create new CDS
        new_cds = CDS(
            start=codon.start,
            end=codon.end,
            strand=codon.strand,
            phase=0,
            id=f"{transcript.id}_{codon.codon_type}_cds"
        )
        transcript.cds_regions.append(new_cds)


class SequenceValidator:
    """Phase B: Validate and correct sequences."""
    
    def __init__(self, min_aa_length: int = 50, max_internal_stops: int = 0, genome=None):
        self.min_aa_length = min_aa_length
        self.max_internal_stops = max_internal_stops
        self.genome = genome
    
    def validate_transcript(self, transcript: Transcript) -> None:
        """Validate and correct a single transcript's sequences."""
        # Length validation
        if len(transcript.aa_sequence) < self.min_aa_length:
            transcript.quality_flags.add("short_transcript")
        
        # Codon validation
        if not transcript.start_codon and not transcript.stop_codon:
            transcript.quality_flags.add("low_quality_no_codons")
        
        # AA sequence validation
        self._validate_aa_sequence(transcript)
        
        # CDS sequence reconstruction if genome available
        if self.genome:
            self._reconstruct_cds_sequence(transcript)
        else:
            transcript.quality_flags.add("validated_cds")
    
    def _validate_aa_sequence(self, transcript: Transcript) -> None:
        """Validate amino acid sequence."""
        aa_seq = transcript.aa_sequence.rstrip('*')  # Remove terminal stop
        
        if not aa_seq:
            transcript.quality_flags.add("empty_aa_sequence")
            return
        
        # Check for start codon
        if not aa_seq.startswith('M'):
            transcript.quality_flags.add("low_quality_no_start_codon")
        
        # Check for internal stop codons
        internal_stops = aa_seq.count('*')
        if internal_stops > self.max_internal_stops:
            transcript.quality_flags.add("internal_stop_codons")
        else:
            transcript.quality_flags.add("passed_aa_validation")
    
    def _reconstruct_cds_sequence(self, transcript: Transcript) -> None:
        """Reconstruct CDS sequence from genome coordinates."""
        if not transcript.cds_regions:
            return
        
        try:
            # Sort CDS regions by coordinate
            sorted_cds = sorted(transcript.cds_regions, key=lambda x: x.start)
            if transcript.strand == '-':
                sorted_cds = sorted(sorted_cds, key=lambda x: x.start, reverse=True)
            
            # Reconstruct sequence
            reconstructed_seq = ""
            for cds in sorted_cds:
                cds_segment = self._extract_genome_sequence(
                    transcript.chrom, cds.start, cds.end, transcript.strand
                )
                reconstructed_seq += cds_segment
            
            # Update transcript sequence
            if reconstructed_seq:
                transcript.cds_sequence = reconstructed_seq
                transcript.quality_flags.add("cds_reconstructed_from_merged_coordinates")
            
        except Exception as e:
            logging.warning(f"Failed to reconstruct CDS for {transcript.id}: {e}")
            transcript.quality_flags.add("cds_reconstruction_failed")
    
    def _extract_genome_sequence(self, chrom: str, start: int, end: int, strand: str) -> str:
        """Extract sequence from genome."""
        if not self.genome:
            return ""
        
        try:
            # pyfaidx uses 0-based indexing, GFF uses 1-based
            sequence = str(self.genome[chrom][start-1:end])
            
            if strand == '-':
                sequence = self._reverse_complement(sequence)
            
            return sequence.upper()
        
        except Exception as e:
            logging.warning(f"Failed to extract sequence {chrom}:{start}-{end}: {e}")
            return ""
    
    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement_map.get(base, base) for base in reversed(sequence.upper()))


class TranscriptSelector:
    """Phase C: Select representative transcripts."""
    
    def __init__(self, overlap_threshold: float = 0.5, min_aa_length: int = 50):
        self.overlap_threshold = overlap_threshold
        self.min_aa_length = min_aa_length
        self.spatial_conflicts: List[Dict] = []
    
    def select_representative(self, gene: Gene) -> Optional[Transcript]:
        """Select representative transcript for a gene using sequence-first approach with hard filtering."""
        if not gene.transcripts:
            return None
        
        # Step 1: Filter by sequence quality only (hard filtering - no fallback)
        quality_transcripts = self._filter_by_quality(gene.transcripts)
        
        if not quality_transcripts:
            # No transcripts meet quality criteria - gene will be excluded from output
            return None
        
        # Step 2: Group by sequence similarity (MD5 hash)
        hash_groups = self._group_by_sequence_hash(quality_transcripts)
        
        # Step 3: Select from most frequent hash group within this gene
        best_group = max(hash_groups.values(), key=len)
        representative = max(best_group, key=lambda t: len(t.aa_sequence))
        
        # Mark as representative
        representative.quality_flags.add("representative")
        gene.representative = representative
        
        return representative
    
    def _filter_by_quality(self, transcripts: List[Transcript]) -> List[Transcript]:
        """Filter transcripts by sequence quality only."""
        quality_transcripts = []
        
        for transcript in transcripts:
            # Check quality flags
            if "short_transcript" in transcript.quality_flags:
                continue
            if "low_quality_no_codons" in transcript.quality_flags:
                continue
            if "internal_stop_codons" in transcript.quality_flags:
                continue
            
            quality_transcripts.append(transcript)
        
        return quality_transcripts
    
    def _group_by_sequence_hash(self, transcripts: List[Transcript]) -> Dict[str, List[Transcript]]:
        """Group transcripts by sequence MD5 hash."""
        hash_groups = defaultdict(list)
        
        for transcript in transcripts:
            # Calculate hash for AA sequence
            aa_hash = hashlib.md5(transcript.aa_sequence.encode()).hexdigest()
            hash_groups[aa_hash].append(transcript)
        
        return hash_groups
    
    def assess_spatial_conflicts(self, genes: Dict[str, Gene]) -> None:
        """Assess spatial conflicts between selected representatives (post-selection)."""
        if not INTERVALTREE_AVAILABLE:
            logging.warning("IntervalTree not available, skipping spatial conflict detection")
            return
        
        # Build interval tree with representatives only
        representatives = []
        for gene in genes.values():
            if gene.representative:
                representatives.append(gene.representative)
        
        # Group by chromosome for efficiency
        chrom_trees = defaultdict(lambda: IntervalTree())
        
        for transcript in representatives:
            chrom_trees[transcript.chrom].addi(
                transcript.start, transcript.end, transcript
            )
        
        # Detect overlaps
        conflicts_detected = 0
        for chrom, tree in chrom_trees.items():
            for transcript in representatives:
                if transcript.chrom != chrom:
                    continue
                
                overlaps = tree.overlap(transcript.start, transcript.end)
                for interval in overlaps:
                    other_transcript = interval.data
                    
                    # Skip self-overlaps and same-gene overlaps
                    if (other_transcript.id == transcript.id or 
                        other_transcript.gene_id == transcript.gene_id):
                        continue
                    
                    # Check if overlap is significant
                    if self._is_significant_overlap(transcript, other_transcript):
                        self._record_spatial_conflict(transcript, other_transcript)
                        conflicts_detected += 1
        
        logging.info(f"Detected {conflicts_detected} spatial conflicts (post-selection)")
    
    def _is_significant_overlap(self, t1: Transcript, t2: Transcript) -> bool:
        """Check if overlap between transcripts is significant."""
        overlap_start = max(t1.start, t2.start)
        overlap_end = min(t1.end, t2.end)
        
        if overlap_start >= overlap_end:
            return False
        
        overlap_length = overlap_end - overlap_start
        t1_length = t1.end - t1.start
        t2_length = t2.end - t2.start
        
        # Reciprocal overlap check
        overlap_ratio_1 = overlap_length / t1_length
        overlap_ratio_2 = overlap_length / t2_length
        
        return (overlap_ratio_1 >= self.overlap_threshold and 
                overlap_ratio_2 >= self.overlap_threshold)
    
    def _record_spatial_conflict(self, t1: Transcript, t2: Transcript) -> None:
        """Record a spatial conflict for tracking."""
        # Check if transcripts encode identical proteins
        if t1.aa_sequence == t2.aa_sequence:
            conflict_type = "redundant_gene_model"
        else:
            conflict_type = "real_spatial_conflict"
        
        # Mark transcripts
        t1.quality_flags.add("post_selection_overlap")
        t2.quality_flags.add("post_selection_overlap")
        
        # Record conflict
        self.spatial_conflicts.append({
            'transcript_1': t1.id,
            'transcript_2': t2.id,
            'gene_1': t1.gene_id,
            'gene_2': t2.gene_id,
            'type': conflict_type,
            'identical_proteins': t1.aa_sequence == t2.aa_sequence
        })
    
    def get_spatial_conflict_summary(self) -> Dict:
        """Get summary of spatial conflicts detected."""
        total_conflicts = len(self.spatial_conflicts)
        redundant_models = sum(1 for c in self.spatial_conflicts 
                              if c['type'] == 'redundant_gene_model')
        real_conflicts = total_conflicts - redundant_models
        
        return {
            'total_conflicts': total_conflicts,
            'redundant_gene_models': redundant_models,
            'real_spatial_conflicts': real_conflicts,
            'conflicts': self.spatial_conflicts
        }