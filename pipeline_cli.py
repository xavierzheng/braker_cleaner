#!/usr/bin/env python3

"""
Command-line interface for the gene curation pipeline.

Provides the same CLI interface as the original gene_curation_pipeline.py
while using the new modular architecture.
"""

import argparse
import sys
import os
import logging
from pathlib import Path

# Import the modular pipeline components
from gene_curation_pipeline.core.config import PipelineConfig, load_config
from gene_curation_pipeline.core.exceptions import PipelineError

def setup_logging(log_level: str = "INFO") -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def create_argument_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Gene annotation curation pipeline for BRAKER predictions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python pipeline_cli.py --gene-model braker.gtf --cds braker.codingseq --aa braker.aa --genome genome.fa --output-dir results

  # With custom parameters
  python pipeline_cli.py --gene-model braker.gtf --cds braker.codingseq --aa braker.aa --genome genome.fa --min-length 20 --overlap-threshold 0.5 --output-dir results
        """
    )
    
    # Required arguments
    parser.add_argument(
        '--gene-model',
        required=True,
        help='Input gene model file (GFF3 or GTF format)'
    )
    parser.add_argument(
        '--cds',
        required=True,
        help='Input CDS sequences FASTA file'
    )
    parser.add_argument(
        '--aa',
        required=True,
        help='Input amino acid sequences FASTA file'
    )
    parser.add_argument(
        '--genome',
        required=True,
        help='Reference genome FASTA file (REQUIRED for codon integration)'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for cleaned files'
    )
    
    # Optional parameters
    parser.add_argument(
        '--min-length',
        type=int,
        default=50,
        help='Minimum amino acid length threshold (default: 50)'
    )
    parser.add_argument(
        '--overlap-threshold',
        type=float,
        default=0.5,
        help='Overlap threshold for spatial conflict detection (default: 0.5). Use 0 for strict non-overlap enforcement.'
    )
    parser.add_argument(
        '--config',
        help='Configuration file (JSON or YAML)'
    )
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging level (default: INFO)'
    )
    
    # Advanced options
    parser.add_argument(
        '--memory-limit',
        type=int,
        help='Memory limit in MB (default: 4096)'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        help='Batch size for processing (default: 1000)'
    )
    
    return parser

def validate_input_files(args) -> None:
    """Validate that input files exist."""
    input_files = {
        'gene-model': args.gene_model,
        'cds': args.cds,
        'aa': args.aa,
        'genome': args.genome
    }
    
    for file_type, file_path in input_files.items():
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{file_type} file not found: {file_path}")

def main():
    """Main CLI entry point."""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log_level)
    logger = logging.getLogger(__name__)
    
    try:
        # Validate input files
        validate_input_files(args)
        
        # Load configuration
        config = load_config(config_path=args.config, use_env=True)
        
        # Override config with command line arguments
        if args.min_length is not None:
            config.min_aa_length = args.min_length
        if args.overlap_threshold is not None:
            config.overlap_threshold = args.overlap_threshold
        if args.memory_limit is not None:
            config.memory_limit_mb = args.memory_limit
        if args.batch_size is not None:
            config.batch_size = args.batch_size
        
        # Re-validate after CLI overrides.
        config.validate()
        
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("Starting gene curation pipeline...")
        logger.info(f"Gene model: {args.gene_model}")
        logger.info(f"CDS file: {args.cds}")
        logger.info(f"AA file: {args.aa}")
        logger.info(f"Genome file: {args.genome}")
        logger.info(f"Output directory: {args.output_dir}")
        logger.info(f"Min AA length: {config.min_aa_length}")
        logger.info(f"Overlap threshold: {config.overlap_threshold}")
        
        # Initialize and run the pipeline
        from gene_curation_pipeline import GeneCurationPipeline
        
        pipeline = GeneCurationPipeline(config)
        success = pipeline.run(
            gene_model_file=args.gene_model,
            cds_file=args.cds,
            aa_file=args.aa,
            genome_file=args.genome,
            output_dir=args.output_dir
        )
        
        if success:
            logger.info("Pipeline completed successfully!")
            return 0
        else:
            logger.error("Pipeline failed!")
            return 1
        
    except FileNotFoundError as e:
        logger.error(f"File error: {e}")
        return 1
    except PipelineError as e:
        logger.error(f"Pipeline error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        logger.debug("Full traceback:", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())
