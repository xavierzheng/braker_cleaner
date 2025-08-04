#!/usr/bin/env python3

"""
Main pipeline class for gene annotation curation.

Integrates all processing phases with the modular architecture.
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Import modular components
from .config import PipelineConfig
from .data_structures import Gene, Transcript
from .exceptions import PipelineError, ParseError, ValidationError, SequenceError, GenomeError
from ..utils.performance_monitor import PerformanceMonitor

# Import processing classes
from .parsers import GeneAnnotationParser, SequenceHandler
from .processors import AnnotationCurator, SequenceValidator, TranscriptSelector
from .generators import OutputGenerator


class GeneCurationPipeline:
    """Main pipeline class that coordinates all processing phases."""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.monitor = PerformanceMonitor(memory_limit_mb=config.memory_limit_mb)
        self.genes: Dict[str, Gene] = {}
        self.input_format: str = ""
        
        # Initialize processors
        self.parser: Optional[GeneAnnotationParser] = None
        self.sequence_handler: Optional[SequenceHandler] = None
        self.curator: Optional[AnnotationCurator] = None
        self.validator: Optional[SequenceValidator] = None
        self.selector: Optional[TranscriptSelector] = None
        self.generator: Optional[OutputGenerator] = None
    
    def run(self, gene_model_file: str, cds_file: str, aa_file: str, 
            genome_file: Optional[str], output_dir: str) -> bool:
        """
        Run the complete gene curation pipeline.
        
        Args:
            gene_model_file: Path to gene model file (GTF/GFF3)
            cds_file: Path to CDS sequences file
            aa_file: Path to amino acid sequences file
            genome_file: Path to genome FASTA file (optional but recommended)
            output_dir: Output directory path
            
        Returns:
            True if pipeline completed successfully
        """
        try:
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Set up logging
            self._setup_pipeline_logging(output_dir)
            
            logging.info("Starting Gene Annotation Curation Pipeline")
            logging.info(f"Configuration: {self.config}")
            logging.info(f"Input files: Gene Model={gene_model_file}, CDS={cds_file}, AA={aa_file}")
            if genome_file:
                logging.info(f"Genome file: {genome_file}")
            logging.info(f"Output directory: {output_dir}")
            
            # Phase 1: Parse input files
            self._parse_inputs(gene_model_file, cds_file, aa_file, genome_file)
            
            # Phase 2: Curate annotations (codon integration)
            self._curate_annotations()
            
            # Phase 3: Validate and correct sequences
            self._validate_sequences()
            
            # Phase 4: Select representative transcripts
            self._select_representatives()
            
            # Phase 5: Generate output files
            self._generate_outputs(output_dir)
            
            # Generate final report
            self._generate_final_report(output_dir)
            
            logging.info("Pipeline completed successfully")
            self.monitor.log_performance_report()
            
            return True
            
        except Exception as e:
            logging.error(f"Pipeline failed: {e}")
            logging.debug("Full traceback:", exc_info=True)
            return False
    
    def _setup_pipeline_logging(self, output_dir: str) -> None:
        """Set up pipeline-specific logging."""
        log_file = Path(output_dir) / 'gene_curation.log'
        
        # Add file handler to root logger
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        ))
        
        root_logger = logging.getLogger()
        root_logger.addHandler(file_handler)
        
        if self.config.debug_mode:
            root_logger.setLevel(logging.DEBUG)
    
    def _parse_inputs(self, gene_model_file: str, cds_file: str, 
                     aa_file: str, genome_file: Optional[str]) -> None:
        """Parse all input files."""
        with self.monitor.phase_context("input_parsing") as metrics:
            try:
                # Parse gene model
                logging.info("Parsing gene model file...")
                self.parser = GeneAnnotationParser(gene_model_file)
                self.genes, self.input_format = self.parser.parse_gff3()
                
                # Parse sequences
                logging.info("Parsing sequence files...")
                self.sequence_handler = SequenceHandler(cds_file, aa_file, genome_file)
                
                # Link sequences to transcripts
                logging.info("Linking sequences to transcripts...")
                self._link_sequences()
                
                # Record operations
                total_transcripts = sum(len(gene.transcripts) for gene in self.genes.values())
                metrics.operations_count = len(self.genes) + total_transcripts
                
                logging.info(f"Parsed {len(self.genes)} genes with {total_transcripts} transcripts")
                
            except Exception as e:
                raise ParseError(f"Failed to parse input files: {e}")
    
    def _link_sequences(self) -> None:
        """Link sequence data to transcript objects."""
        for gene in self.genes.values():
            for transcript in gene.transcripts:
                # Link CDS sequence
                if transcript.id in self.sequence_handler.cds_sequences:
                    transcript.cds_sequence = self.sequence_handler.cds_sequences[transcript.id]
                
                # Link AA sequence
                if transcript.id in self.sequence_handler.aa_sequences:
                    transcript.aa_sequence = self.sequence_handler.aa_sequences[transcript.id]
    
    def _curate_annotations(self) -> None:
        """Curate annotations with codon integration."""
        with self.monitor.phase_context("annotation_curation") as metrics:
            try:
                logging.info("Starting annotation curation...")
                
                if not self.config.enable_codon_integration:
                    logging.info("Codon integration disabled, skipping curation")
                    return
                
                # Initialize curator
                self.curator = AnnotationCurator(
                    self.config.adjacency_gap_threshold,
                    self.sequence_handler.genome if self.sequence_handler else None
                )
                
                # Process each gene
                processed_genes = 0
                for gene in self.genes.values():
                    self.curator.curate_gene(gene)
                    processed_genes += 1
                    
                    # Check memory periodically
                    if processed_genes % self.config.batch_size == 0:
                        self.monitor.check_memory_limit()
                
                metrics.operations_count = processed_genes
                logging.info(f"Curated {processed_genes} genes")
                
            except Exception as e:
                raise PipelineError(f"Failed during annotation curation: {e}")
    
    def _validate_sequences(self) -> None:
        """Validate and correct sequences."""
        with self.monitor.phase_context("sequence_validation") as metrics:
            try:
                logging.info("Starting sequence validation...")
                
                # Initialize validator
                self.validator = SequenceValidator(
                    self.config.min_aa_length,
                    self.config.max_internal_stops,
                    self.sequence_handler.genome if self.sequence_handler else None
                )
                
                # Process transcripts
                processed_transcripts = 0
                for gene in self.genes.values():
                    for transcript in gene.transcripts:
                        self.validator.validate_transcript(transcript)
                        processed_transcripts += 1
                
                metrics.operations_count = processed_transcripts
                logging.info(f"Validated {processed_transcripts} transcripts")
                
            except Exception as e:
                raise ValidationError(f"Failed during sequence validation: {e}")
    
    def _select_representatives(self) -> None:
        """Select representative transcripts."""
        with self.monitor.phase_context("representative_selection") as metrics:
            try:
                logging.info("Starting representative selection...")
                
                # Initialize selector
                self.selector = TranscriptSelector(
                    self.config.overlap_threshold,
                    self.config.min_aa_length
                )
                
                # Select representatives for each gene
                selected_count = 0
                for gene in self.genes.values():
                    selected = self.selector.select_representative(gene)
                    if selected:
                        selected_count += 1
                
                metrics.operations_count = len(self.genes)
                logging.info(f"Selected representatives for {selected_count}/{len(self.genes)} genes")
                
            except Exception as e:
                raise PipelineError(f"Failed during representative selection: {e}")
    
    def _generate_outputs(self, output_dir: str) -> None:
        """Generate output files."""
        with self.monitor.phase_context("output_generation") as metrics:
            try:
                logging.info("Generating output files...")
                
                # Initialize generator
                self.generator = OutputGenerator(
                    self.input_format,
                    self.config.preserve_source_names,
                    self.config.include_quality_flags
                )
                
                # Generate outputs
                output_files = self.generator.generate_outputs(
                    self.genes, 
                    output_dir,
                    self.sequence_handler.cds_sequences if self.sequence_handler else {},
                    self.sequence_handler.aa_sequences if self.sequence_handler else {}
                )
                
                metrics.operations_count = len(output_files)
                logging.info(f"Generated {len(output_files)} output files")
                
                for file_path in output_files:
                    logging.info(f"Created: {file_path}")
                
            except Exception as e:
                raise PipelineError(f"Failed during output generation: {e}")
    
    def _generate_final_report(self, output_dir: str) -> None:
        """Generate comprehensive processing report."""
        try:
            report_file = Path(output_dir) / 'processing_report.txt'
            
            # Collect statistics
            total_genes = len(self.genes)
            total_transcripts = sum(len(gene.transcripts) for gene in self.genes.values())
            selected_genes = sum(1 for gene in self.genes.values() if gene.representative_transcript)
            
            # Performance summary
            performance = self.monitor.get_performance_summary()
            
            # Write report
            with open(report_file, 'w') as f:
                f.write("Gene Curation Pipeline - Processing Report\n")
                f.write("=" * 50 + "\n\n")
                
                f.write("INPUT STATISTICS\n")
                f.write("-" * 20 + "\n")
                f.write(f"Total genes: {total_genes:,}\n")
                f.write(f"Total transcripts: {total_transcripts:,}\n")
                f.write(f"Input format: {self.input_format}\n\n")
                
                f.write("PROCESSING RESULTS\n")
                f.write("-" * 20 + "\n")
                f.write(f"Genes with representatives: {selected_genes:,} ({selected_genes/total_genes*100:.1f}%)\n")
                f.write(f"Selection success rate: {selected_genes/total_genes*100:.1f}%\n\n")
                
                f.write("PERFORMANCE METRICS\n")
                f.write("-" * 20 + "\n")
                f.write(f"Total processing time: {performance['total_elapsed_time']:.2f} seconds\n")
                f.write(f"Peak memory usage: {performance['peak_memory_mb']:.1f} MB\n")
                f.write(f"Memory efficiency: {performance['peak_memory_mb']/self.config.memory_limit_mb*100:.1f}%\n\n")
                
                if performance['phases']:
                    f.write("PHASE BREAKDOWN\n")
                    f.write("-" * 20 + "\n")
                    for phase_name, phase_data in performance['phases'].items():
                        f.write(f"{phase_name}: {phase_data['elapsed_time']:.2f}s ")
                        f.write(f"({phase_data['operations_count']} operations)\n")
                
                f.write(f"\nConfiguration used:\n")
                config_dict = self.config.to_dict()
                for key, value in config_dict.items():
                    f.write(f"  {key}: {value}\n")
            
            logging.info(f"Generated processing report: {report_file}")
            
        except Exception as e:
            logging.warning(f"Failed to generate processing report: {e}")