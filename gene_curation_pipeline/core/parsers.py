#!/usr/bin/env python3

"""
File parsers for gene annotations and sequences.

Handles GFF3/GTF parsing and FASTA sequence processing.
"""

import logging
from typing import Dict, Tuple, Optional
from .data_structures import Gene, Transcript, Exon, CDS, Codon
from .exceptions import ParseError

try:
    import pyfaidx
    PYFAIDX_AVAILABLE = True
except ImportError:
    PYFAIDX_AVAILABLE = False


class GeneAnnotationParser:
    """Parse GFF3/GTF files efficiently with O(n) complexity."""
    
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.genes: Dict[str, Gene] = {}
        self.transcripts: Dict[str, Transcript] = {}
    
    def parse_gff3(self) -> Tuple[Dict[str, Gene], str]:
        """Parse GFF3/GTF file with O(n) complexity."""
        file_type = "GTF" if self.file_path.lower().endswith('.gtf') else "GFF3"
        logging.info(f"Parsing {file_type} file: {self.file_path}")
        
        gene_count = 0
        transcript_count = 0
        
        try:
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
                            gene_id = self._extract_gene_id(attr_dict, attributes, file_type)
                            if gene_id:
                                self.genes[gene_id] = Gene(id=gene_id)
                                gene_count += 1
                        
                        elif feature in ['mRNA', 'transcript']:
                            transcript_id, parent_id = self._extract_transcript_ids(attr_dict, attributes, file_type)
                            if transcript_id and parent_id:
                                self._add_transcript(transcript_id, parent_id, start, end, strand, chrom, source)
                                transcript_count += 1
                        
                        elif feature in ['exon', 'CDS', 'start_codon', 'stop_codon']:
                            parent_id = self._extract_parent_id(attr_dict, file_type)
                            if parent_id in self.transcripts:
                                self._add_feature(feature, start, end, strand, phase, attr_dict, parent_id)
                    
                    except Exception as e:
                        logging.warning(f"Error parsing line {line_num}: {e}")
                        continue
        
        except FileNotFoundError:
            raise ParseError(f"Gene model file not found: {self.file_path}")
        except Exception as e:
            raise ParseError(f"Failed to parse {file_type} file: {e}", self.file_path)
        
        logging.info(f"Parsed {gene_count} genes and {transcript_count} transcripts")
        return self.genes, file_type
    
    def _extract_gene_id(self, attr_dict: Dict[str, str], attributes: str, file_type: str) -> str:
        """Extract gene ID from attributes."""
        if file_type == "GTF":
            gene_id = attr_dict.get('gene_id', '')
            # Handle bare gene IDs (like 'g1' without attributes)
            if not gene_id and attributes.strip():
                gene_id = attributes.strip()
        else:
            gene_id = attr_dict.get('ID', '')
        return gene_id
    
    def _extract_transcript_ids(self, attr_dict: Dict[str, str], attributes: str, file_type: str) -> Tuple[str, str]:
        """Extract transcript and parent IDs."""
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
        
        return transcript_id, parent_id
    
    def _extract_parent_id(self, attr_dict: Dict[str, str], file_type: str) -> str:
        """Extract parent ID for features."""
        if file_type == "GTF":
            return attr_dict.get('transcript_id', '')
        else:
            return attr_dict.get('Parent', '')
    
    def _add_transcript(self, transcript_id: str, parent_id: str, start: int, end: int, 
                       strand: str, chrom: str, source: str) -> None:
        """Add transcript to the collection, avoiding duplicates."""
        # Check if transcript already exists to avoid duplicates
        if transcript_id not in self.transcripts:
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
        else:
            # Update coordinates if this feature has different boundaries
            existing_transcript = self.transcripts[transcript_id]
            existing_transcript.start = min(existing_transcript.start, start)
            existing_transcript.end = max(existing_transcript.end, end)
    
    def _add_feature(self, feature: str, start: int, end: int, strand: str, 
                    phase: str, attr_dict: Dict[str, str], parent_id: str) -> None:
        """Add feature to the appropriate transcript."""
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
    """Handle FASTA sequence files for CDS and amino acid sequences."""
    
    def __init__(self, cds_file: str, aa_file: str, genome_file: Optional[str] = None):
        self.cds_file = cds_file
        self.aa_file = aa_file
        self.genome_file = genome_file
        
        # Sequence storage
        self.cds_sequences: Dict[str, str] = {}
        self.aa_sequences: Dict[str, str] = {}
        self.genome: Optional['pyfaidx.Fasta'] = None
        
        # Load sequences
        self._load_sequences()
    
    def _load_sequences(self) -> None:
        """Load all sequence files."""
        try:
            # Load CDS sequences
            logging.info(f"Loading CDS sequences from {self.cds_file}")
            self.cds_sequences = self._parse_fasta(self.cds_file)
            logging.info(f"Loaded {len(self.cds_sequences)} CDS sequences")
            
            # Load amino acid sequences
            logging.info(f"Loading AA sequences from {self.aa_file}")
            self.aa_sequences = self._parse_fasta(self.aa_file)
            logging.info(f"Loaded {len(self.aa_sequences)} AA sequences")
            
            # Load genome if provided
            if self.genome_file:
                if not PYFAIDX_AVAILABLE:
                    logging.warning("pyfaidx not available, genome-based operations disabled")
                else:
                    logging.info(f"Loading genome from {self.genome_file}")
                    self.genome = pyfaidx.Fasta(self.genome_file)
                    logging.info("Genome loaded successfully")
            
        except Exception as e:
            raise ParseError(f"Failed to load sequence files: {e}")
    
    def _parse_fasta(self, file_path: str) -> Dict[str, str]:
        """Parse FASTA file into sequence dictionary."""
        sequences = {}
        current_id = None
        current_seq = []
        
        try:
            with open(file_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    if line.startswith('>'):
                        # Save previous sequence
                        if current_id and current_seq:
                            sequences[current_id] = ''.join(current_seq)
                        
                        # Start new sequence
                        current_id = line[1:].split()[0]  # Take first part of header
                        current_seq = []
                    
                    elif line and current_id:
                        current_seq.append(line)
                
                # Save last sequence
                if current_id and current_seq:
                    sequences[current_id] = ''.join(current_seq)
        
        except FileNotFoundError:
            raise ParseError(f"Sequence file not found: {file_path}")
        except Exception as e:
            raise ParseError(f"Failed to parse FASTA file {file_path}: {e}", file_path, line_num)
        
        return sequences
    
    def get_genome_sequence(self, chrom: str, start: int, end: int, strand: str = '+') -> str:
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
            logging.warning(f"Failed to extract genome sequence {chrom}:{start}-{end}: {e}")
            return ""
    
    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement_map.get(base, base) for base in reversed(sequence.upper()))