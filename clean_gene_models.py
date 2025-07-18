#!/usr/bin/env python3
"""
Gene Model Cleaning Script

This script processes BRAKER gene prediction output to select the best representative 
transcript for each gene by analyzing protein sequences and removing redundant 
alternative splice forms.

Based on the workflow defined in CLAUDE.md
"""

import os
import re
import hashlib
import logging
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple, Optional
import argparse


class GeneModelCleaner:
    """Main class for cleaning gene models"""
    
    def __init__(self, input_dir: str = ".", min_aa_length: int = 10):
        self.input_dir = input_dir
        self.min_aa_length = min_aa_length
        self.input_files = {
            'aa': os.path.join(input_dir, 'braker.aa'),
            'cds': os.path.join(input_dir, 'braker.cds.fa'),
            'gtf': os.path.join(input_dir, 'braker.gtf'),
            'gff3': os.path.join(input_dir, 'braker.gff3')
        }
        self.output_files = {
            'aa': os.path.join(input_dir, 'braker_cleaned.aa'),
            'cds': os.path.join(input_dir, 'braker_cleaned.cds.fa'),
            'gtf': os.path.join(input_dir, 'braker_cleaned.gtf'),
            'gff3': os.path.join(input_dir, 'braker_cleaned.gff3'),
            'manual_review': os.path.join(input_dir, 'manual_review_transcripts.txt'),
            'report': os.path.join(input_dir, 'processing_report.txt')
        }
        
        # Data structures
        self.protein_sequences = {}
        self.cds_sequences = {}
        self.gtf_records = defaultdict(list)
        self.gff3_records = defaultdict(list)
        self.gff3_gene_records = {}  # Store gene records
        self.transcript_metrics = {}
        self.gene_transcripts = defaultdict(list)
        self.selected_transcripts = set()
        self.manual_review_cases = []
        self.filtered_transcripts = set()  # Track filtered out transcripts
        
        # Setup logging
        self.setup_logging()
    
    def setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('gene_cleaning.log'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def parse_transcript_id(self, transcript_id: str) -> Tuple[str, str]:
        """Parse transcript ID to extract gene and transcript parts"""
        match = re.match(r'([^.]+)\.([^.]+)', transcript_id)
        if match:
            return match.group(1), match.group(2)
        else:
            self.logger.warning(f"Could not parse transcript ID: {transcript_id}")
            return transcript_id, "t1"
    
    def calculate_md5(self, sequence: str) -> str:
        """Calculate MD5 hash of a sequence"""
        return hashlib.md5(sequence.encode()).hexdigest()
    
    def read_fasta(self, filepath: str) -> Dict[str, str]:
        """Read FASTA file and return sequences"""
        sequences = {}
        current_id = None
        current_seq = []
        
        with open(filepath, 'r') as f:
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
    
    def filter_short_sequences(self):
        """Filter out transcripts with amino acid length < min_aa_length"""
        self.logger.info(f"Filtering transcripts with amino acid length < {self.min_aa_length}")
        
        filtered_sequences = {}
        for transcript_id, sequence in self.protein_sequences.items():
            if len(sequence) >= self.min_aa_length:
                filtered_sequences[transcript_id] = sequence
            else:
                self.filtered_transcripts.add(transcript_id)
                self.logger.debug(f"Filtered out {transcript_id}: length {len(sequence)} < {self.min_aa_length}")
        
        self.protein_sequences = filtered_sequences
        self.logger.info(f"Filtered out {len(self.filtered_transcripts)} transcripts with short amino acid sequences")
        self.logger.info(f"Remaining transcripts: {len(self.protein_sequences)}")
    
    def clean_protein_sequences(self):
        """Step 1: Clean protein sequences by removing stop codon markers"""
        self.logger.info("Step 1: Cleaning protein sequences...")
        
        raw_sequences = self.read_fasta(self.input_files['aa'])
        
        for seq_id, sequence in raw_sequences.items():
            # Remove trailing '*' (stop codon marker)
            cleaned_seq = sequence.rstrip('*')
            self.protein_sequences[seq_id] = cleaned_seq
            
        self.logger.info(f"Cleaned {len(self.protein_sequences)} protein sequences")
        
        # Apply length filter
        self.filter_short_sequences()
    
    def calculate_sequence_metrics(self):
        """Step 2: Calculate sequence metrics and MD5 hashes"""
        self.logger.info("Step 2: Calculating sequence metrics...")
        
        # Read CDS sequences
        all_cds_sequences = self.read_fasta(self.input_files['cds'])
        
        # Filter CDS sequences to match remaining protein sequences
        self.cds_sequences = {tid: seq for tid, seq in all_cds_sequences.items() 
                             if tid in self.protein_sequences}
        
        # Calculate metrics for each transcript
        for transcript_id in self.protein_sequences:
            protein_seq = self.protein_sequences[transcript_id]
            cds_seq = self.cds_sequences.get(transcript_id, "")
            
            # Get CDS positions from GTF/GFF3 (will be updated in parse_annotations)
            cds_positions = []
            
            self.transcript_metrics[transcript_id] = {
                'protein_length': len(protein_seq),
                'protein_md5': self.calculate_md5(protein_seq),
                'cds_md5': self.calculate_md5(cds_seq),
                'cds_positions': cds_positions,
                'cds_positions_md5': ""  # Will be calculated after positions are extracted
            }
        
        self.logger.info(f"Calculated metrics for {len(self.transcript_metrics)} transcripts")
    
    def parse_annotations(self):
        """Parse GTF and GFF3 files to extract transcript annotations"""
        self.logger.info("Step 3: Parsing gene/transcript identifiers...")
        
        # Parse GTF file
        self.parse_gtf()
        
        # Parse GFF3 file
        self.parse_gff3()
        
        # Update CDS positions MD5 after parsing
        for transcript_id in self.transcript_metrics:
            cds_positions = self.transcript_metrics[transcript_id]['cds_positions']
            self.transcript_metrics[transcript_id]['cds_positions_md5'] = self.calculate_md5(str(sorted(cds_positions)))
        
        # Group transcripts by gene (only for remaining transcripts)
        for transcript_id in self.protein_sequences:
            gene_id, _ = self.parse_transcript_id(transcript_id)
            self.gene_transcripts[gene_id].append(transcript_id)
        
        self.logger.info(f"Found {len(self.gene_transcripts)} genes with transcripts")
    
    def parse_gtf(self):
        """Parse GTF file"""
        if not os.path.exists(self.input_files['gtf']):
            self.logger.warning("GTF file not found")
            return
            
        # Store gene records temporarily to add them later
        gene_records = {}
        
        with open(self.input_files['gtf'], 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                    
                feature_type = parts[2]
                start, end = int(parts[3]), int(parts[4])
                attributes = parts[8]
                
                # Handle gene records (no transcript_id)
                if feature_type == 'gene':
                    # Extract gene_id directly from attributes
                    gene_match = re.search(r'gene_id "([^"]+)"', attributes)
                    if not gene_match:
                        # Try alternative format without quotes
                        gene_match = re.search(r'(\w+)', attributes)
                    if gene_match:
                        gene_id = gene_match.group(1)
                        # Store gene record for later processing
                        gene_records[gene_id] = {
                            'feature_type': feature_type,
                            'start': start,
                            'end': end,
                            'line': line.strip()
                        }
                else:
                    # Extract transcript_id for other features
                    transcript_match = re.search(r'transcript_id "([^"]+)"', attributes)
                    if transcript_match:
                        transcript_id = transcript_match.group(1)
                        
                        # Only process if transcript wasn't filtered out
                        if transcript_id in self.protein_sequences:
                            self.gtf_records[transcript_id].append({
                                'feature_type': feature_type,
                                'start': start,
                                'end': end,
                                'line': line.strip()
                            })
                            
                            # Update CDS positions in metrics
                            if feature_type == 'CDS' and transcript_id in self.transcript_metrics:
                                self.transcript_metrics[transcript_id]['cds_positions'].append((start, end))
        
        # Now add gene records to transcripts efficiently
        for transcript_id in self.protein_sequences:
            gene_id, _ = self.parse_transcript_id(transcript_id)
            if gene_id in gene_records:
                self.gtf_records[transcript_id].append(gene_records[gene_id])
    
    def parse_gff3(self):
        """Parse GFF3 file"""
        if not os.path.exists(self.input_files['gff3']):
            self.logger.warning("GFF3 file not found")
            return
            
        with open(self.input_files['gff3'], 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                    
                feature_type = parts[2]
                start, end = int(parts[3]), int(parts[4])
                attributes = parts[8]
                
                # Handle different feature types
                if feature_type == 'gene':
                    gene_id = self.extract_id_from_gff3_attributes(attributes)
                    if gene_id:
                        self.gff3_gene_records[gene_id] = {
                            'feature_type': feature_type,
                            'start': start,
                            'end': end,
                            'line': line.strip()
                        }
                        
                elif feature_type == 'mRNA':
                    transcript_id = self.extract_id_from_gff3_attributes(attributes)
                    if transcript_id and transcript_id in self.protein_sequences:
                        self.gff3_records[transcript_id].append({
                            'feature_type': feature_type,
                            'start': start,
                            'end': end,
                            'line': line.strip()
                        })
                        
                else:
                    # For other features (exon, CDS, start_codon, stop_codon, intron)
                    parent_id = self.extract_parent_from_gff3_attributes(attributes)
                    if parent_id and parent_id in self.protein_sequences:
                        self.gff3_records[parent_id].append({
                            'feature_type': feature_type,
                            'start': start,
                            'end': end,
                            'line': line.strip()
                        })
                        
                        # Update CDS positions for metrics
                        if feature_type == 'CDS' and parent_id in self.transcript_metrics:
                            self.transcript_metrics[parent_id]['cds_positions'].append((start, end))
    
    def extract_id_from_gff3_attributes(self, attributes: str) -> Optional[str]:
        """Extract ID from GFF3 attributes"""
        match = re.search(r'ID=([^;]+)', attributes)
        if match:
            return match.group(1)
        return None
    
    def extract_parent_from_gff3_attributes(self, attributes: str) -> Optional[str]:
        """Extract Parent ID from GFF3 attributes"""
        match = re.search(r'Parent=([^;]+)', attributes)
        if match:
            return match.group(1)
        return None
    
    def select_longest_transcripts(self):
        """Step 4: Select longest transcript per gene"""
        self.logger.info("Step 4: Selecting longest transcript per gene...")
        
        for gene_id, transcripts in self.gene_transcripts.items():
            if len(transcripts) == 1:
                self.selected_transcripts.add(transcripts[0])
                continue
            
            # Find transcript with longest protein sequence
            longest_transcript = max(transcripts, 
                                   key=lambda t: self.transcript_metrics[t]['protein_length'])
            
            # Check if there are multiple transcripts with same length
            max_length = self.transcript_metrics[longest_transcript]['protein_length']
            equal_length_transcripts = [t for t in transcripts 
                                      if self.transcript_metrics[t]['protein_length'] == max_length]
            
            if len(equal_length_transcripts) == 1:
                self.selected_transcripts.add(longest_transcript)
            else:
                # Handle equal-length transcripts
                self.handle_equal_length_transcripts(gene_id, equal_length_transcripts)
    
    def handle_equal_length_transcripts(self, gene_id: str, transcripts: List[str]):
        """Step 5: Handle equal-length transcripts"""
        self.logger.info(f"Step 5: Handling equal-length transcripts for gene {gene_id}")
        
        # Group by CDS MD5 and positions
        identical_groups = defaultdict(list)
        for transcript in transcripts:
            metrics = self.transcript_metrics[transcript]
            key = (metrics['cds_md5'], metrics['cds_positions_md5'])
            identical_groups[key].append(transcript)
        
        # If all transcripts are identical, keep the first one
        if len(identical_groups) == 1:
            self.selected_transcripts.add(transcripts[0])
        else:
            # Handle ambiguous cases
            self.handle_ambiguous_cases(gene_id, transcripts)
    
    def handle_ambiguous_cases(self, gene_id: str, transcripts: List[str]):
        """Step 6: Report ambiguous cases"""
        self.logger.info(f"Step 6: Handling ambiguous cases for gene {gene_id}")
        
        # Group transcripts by amino acid MD5 hash
        protein_groups = defaultdict(list)
        for transcript in transcripts:
            protein_md5 = self.transcript_metrics[transcript]['protein_md5']
            protein_groups[protein_md5].append(transcript)
        
        # Find the group with highest frequency
        max_frequency = max(len(group) for group in protein_groups.values())
        max_frequency_groups = [group for group in protein_groups.values() 
                               if len(group) == max_frequency]
        
        if len(max_frequency_groups) == 1:
            # Single group with highest frequency
            self.selected_transcripts.add(max_frequency_groups[0][0])
        else:
            # Multiple groups with equal highest frequency - needs manual review
            self.manual_review_cases.append({
                'gene_id': gene_id,
                'transcripts': transcripts,
                'reason': 'Multiple protein variants with equal highest frequency'
            })
            # For now, select the first transcript
            self.selected_transcripts.add(transcripts[0])
    
    def write_output_files(self):
        """Create output files with selected transcripts"""
        self.logger.info("Writing output files...")
        
        # Write cleaned protein sequences
        self.write_fasta_output(self.protein_sequences, self.output_files['aa'])
        
        # Write cleaned CDS sequences
        self.write_fasta_output(self.cds_sequences, self.output_files['cds'])
        
        # Write cleaned GTF
        self.write_gtf_output()
        
        # Write cleaned GFF3
        self.write_gff3_output()
        
        # Write manual review cases
        self.write_manual_review_file()
        
        # Write processing report
        self.write_processing_report()
    
    def write_fasta_output(self, sequences: Dict[str, str], output_file: str):
        """Write FASTA sequences for selected transcripts"""
        with open(output_file, 'w') as f:
            for transcript_id in self.selected_transcripts:
                if transcript_id in sequences:
                    f.write(f">{transcript_id}\n")
                    f.write(f"{sequences[transcript_id]}\n")
    
    def write_gtf_output(self):
        """Write GTF records for selected transcripts with proper order"""
        with open(self.output_files['gtf'], 'w') as f:
            # Pre-group selected transcripts by gene (efficient O(n) operation)
            selected_transcripts_by_gene = defaultdict(list)
            for transcript_id in self.selected_transcripts:
                gene_id, _ = self.parse_transcript_id(transcript_id)
                selected_transcripts_by_gene[gene_id].append(transcript_id)
            
            # Define feature order once
            feature_order = {
                'gene': 1,
                'transcript': 2,
                'stop_codon': 3,
                'CDS': 4,
                'exon': 5,
                'start_codon': 6,
                'mRNA': 7
            }
            
            # Write records in proper order
            for gene_id in sorted(selected_transcripts_by_gene.keys()):
                gene_transcripts = selected_transcripts_by_gene[gene_id]
                
                for transcript_id in sorted(gene_transcripts):
                    if transcript_id in self.gtf_records:
                        # Sort records by feature type to maintain proper order
                        records = self.gtf_records[transcript_id]
                        sorted_records = sorted(records, 
                                              key=lambda x: (feature_order.get(x['feature_type'], 99), 
                                                           x['start']))
                        
                        for record in sorted_records:
                            f.write(f"{record['line']}\n")
    
    def write_gff3_output(self):
        """Write GFF3 records for selected transcripts with proper hierarchy"""
        with open(self.output_files['gff3'], 'w') as f:
            # Pre-group selected transcripts by gene (efficient O(n) operation)
            selected_transcripts_by_gene = defaultdict(list)
            for transcript_id in self.selected_transcripts:
                gene_id, _ = self.parse_transcript_id(transcript_id)
                selected_transcripts_by_gene[gene_id].append(transcript_id)
            
            # Define feature order once
            feature_order = {
                'mRNA': 1,
                'exon': 2,
                'CDS': 3,
                'start_codon': 4,
                'stop_codon': 5,
                'intron': 6
            }
            
            # Write records in proper hierarchical order
            for gene_id in sorted(selected_transcripts_by_gene.keys()):
                # Write gene record first
                if gene_id in self.gff3_gene_records:
                    f.write(f"{self.gff3_gene_records[gene_id]['line']}\n")
                
                # Write selected transcripts for this gene
                gene_transcripts = selected_transcripts_by_gene[gene_id]
                
                for transcript_id in sorted(gene_transcripts):
                    if transcript_id in self.gff3_records:
                        # Sort records by feature type to maintain hierarchy
                        records = self.gff3_records[transcript_id]
                        sorted_records = sorted(records, 
                                              key=lambda x: (feature_order.get(x['feature_type'], 99), 
                                                           x['start']))
                        
                        for record in sorted_records:
                            f.write(f"{record['line']}\n")
    
    def write_manual_review_file(self):
        """Write manual review cases to file"""
        with open(self.output_files['manual_review'], 'w') as f:
            f.write("# Transcripts requiring manual review\n")
            f.write("# Format: gene_id\ttranscript_ids\treason\n")
            for case in self.manual_review_cases:
                transcript_list = ','.join(case['transcripts'])
                f.write(f"{case['gene_id']}\t{transcript_list}\t{case['reason']}\n")
    
    def write_processing_report(self):
        """Write processing report"""
        with open(self.output_files['report'], 'w') as f:
            f.write("# Gene Model Cleaning Report\n\n")
            f.write(f"Minimum amino acid length filter: {self.min_aa_length}\n")
            f.write(f"Total genes processed: {len(self.gene_transcripts)}\n")
            f.write(f"Total transcripts in input: {len(self.protein_sequences) + len(self.filtered_transcripts)}\n")
            f.write(f"Transcripts filtered by length: {len(self.filtered_transcripts)}\n")
            f.write(f"Transcripts after filtering: {len(self.protein_sequences)}\n")
            f.write(f"Selected transcripts: {len(self.selected_transcripts)}\n")
            f.write(f"Cases requiring manual review: {len(self.manual_review_cases)}\n\n")
            
            if self.filtered_transcripts:
                f.write("## Filtered Transcripts (too short)\n")
                for transcript_id in sorted(self.filtered_transcripts):
                    f.write(f"{transcript_id}\n")
                f.write("\n")
            
            f.write("## Gene Summary\n")
            for gene_id, transcripts in self.gene_transcripts.items():
                selected = [t for t in transcripts if t in self.selected_transcripts]
                f.write(f"{gene_id}: {len(transcripts)} transcripts -> {len(selected)} selected\n")
    
    def manual_curation(self, approved_transcripts_file: str):
        """Step 7: Manual curation function"""
        self.logger.info("Step 7: Performing manual curation...")
        
        # Read approved transcript IDs
        approved_transcripts = set()
        with open(approved_transcripts_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    approved_transcripts.add(line)
        
        # Update selected transcripts
        self.selected_transcripts = approved_transcripts
        
        # Rewrite output files
        self.write_output_files()
        
        self.logger.info(f"Manual curation completed. {len(approved_transcripts)} transcripts approved.")
    
    def run(self):
        """Run the complete cleaning workflow"""
        self.logger.info("Starting gene model cleaning workflow...")
        self.logger.info(f"Minimum amino acid length filter: {self.min_aa_length}")
        
        # Check input files exist
        for file_type, filepath in self.input_files.items():
            if not os.path.exists(filepath):
                self.logger.error(f"Input file not found: {filepath}")
                return False
        
        try:
            # Execute workflow steps
            self.clean_protein_sequences()
            self.calculate_sequence_metrics()
            self.parse_annotations()
            self.select_longest_transcripts()
            self.write_output_files()
            
            self.logger.info("Gene model cleaning completed successfully!")
            return True
            
        except Exception as e:
            self.logger.error(f"Error during processing: {str(e)}")
            return False


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Clean BRAKER gene models')
    parser.add_argument('--input-dir', default='.', 
                       help='Directory containing input files (default: current directory)')
    parser.add_argument('--min-aa-length', type=int, default=10,
                       help='Minimum amino acid length to keep transcript (default: 10)')
    parser.add_argument('--manual-curation', 
                       help='File containing approved transcript IDs for manual curation')
    
    args = parser.parse_args()
    
    # Create cleaner instance
    cleaner = GeneModelCleaner(args.input_dir, args.min_aa_length)
    
    if args.manual_curation:
        # Run manual curation
        cleaner.manual_curation(args.manual_curation)
    else:
        # Run main workflow
        success = cleaner.run()
        if not success:
            exit(1)


if __name__ == "__main__":
    main()