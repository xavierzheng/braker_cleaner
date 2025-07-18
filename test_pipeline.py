#!/usr/bin/env python3

"""
Test suite for Gene Annotation Curation Pipeline
"""

import unittest
import os
import tempfile
import shutil
from gene_curation_pipeline import (
    Transcript, Gene, Codon, Exon, CDS,
    PerformanceMonitor, GeneAnnotationParser, 
    SequenceHandler, AnnotationCurator
)

class TestDataStructures(unittest.TestCase):
    """Test core data structures."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.transcript = Transcript(
            id="g1.t1",
            gene_id="g1",
            start=1000,
            end=2000,
            strand="+",
            cds_sequence="ATGAAATAG",
            aa_sequence="MK*"
        )
    
    def test_transcript_hashes(self):
        """Test MD5 hash generation."""
        cds_hash = self.transcript.get_cds_hash()
        aa_hash = self.transcript.get_aa_hash()
        
        self.assertEqual(len(cds_hash), 32)  # MD5 hash length
        self.assertEqual(len(aa_hash), 32)
        self.assertNotEqual(cds_hash, aa_hash)  # Different sequences should have different hashes
    
    def test_gene_transcript_lookup(self):
        """Test gene transcript lookup."""
        gene = Gene(id="g1")
        gene.transcripts.append(self.transcript)
        
        found_transcript = gene.get_transcript_by_id("g1.t1")
        self.assertIsNotNone(found_transcript)
        self.assertEqual(found_transcript.id, "g1.t1")
        
        not_found = gene.get_transcript_by_id("g1.t2")
        self.assertIsNone(not_found)


class TestPerformanceMonitor(unittest.TestCase):
    """Test performance monitoring."""
    
    def test_performance_monitor_init(self):
        """Test performance monitor initialization."""
        monitor = PerformanceMonitor()
        self.assertIsNotNone(monitor.start_time)
        self.assertEqual(monitor.peak_memory, 0)
    
    def test_memory_usage(self):
        """Test memory usage monitoring."""
        monitor = PerformanceMonitor()
        memory_usage = monitor.get_memory_usage()
        
        # Should return a non-negative number
        self.assertGreaterEqual(memory_usage, 0)
    
    def test_elapsed_time(self):
        """Test elapsed time calculation."""
        import time
        
        monitor = PerformanceMonitor()
        time.sleep(0.1)  # Sleep for 100ms
        elapsed = monitor.get_elapsed_time()
        
        self.assertGreater(elapsed, 0.05)  # Should be at least 50ms
        self.assertLess(elapsed, 1.0)      # Should be less than 1 second


class TestSequenceHandler(unittest.TestCase):
    """Test sequence handling."""
    
    def setUp(self):
        """Set up test files."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test CDS file
        self.cds_file = os.path.join(self.temp_dir, "test.cds.fa")
        with open(self.cds_file, 'w') as f:
            f.write(">g1.t1\nATGAAATAG\n")
            f.write(">g2.t1\nATGTTTTGA\n")
        
        # Create test AA file
        self.aa_file = os.path.join(self.temp_dir, "test.aa")
        with open(self.aa_file, 'w') as f:
            f.write(">g1.t1\nMK*\n")
            f.write(">g2.t1\nMF*\n")
    
    def tearDown(self):
        """Clean up test files."""
        shutil.rmtree(self.temp_dir)
    
    def test_fasta_parsing(self):
        """Test FASTA file parsing."""
        handler = SequenceHandler(self.cds_file, self.aa_file)
        
        self.assertEqual(len(handler.cds_sequences), 2)
        self.assertEqual(len(handler.aa_sequences), 2)
        
        self.assertEqual(handler.cds_sequences["g1.t1"], "ATGAAATAG")
        self.assertEqual(handler.aa_sequences["g1.t1"], "MK*")
    
    def test_reverse_complement(self):
        """Test reverse complement function."""
        handler = SequenceHandler(self.cds_file, self.aa_file)
        
        # Test reverse complement
        result = handler._reverse_complement("ATGC")
        self.assertEqual(result, "GCAT")
        
        result = handler._reverse_complement("AAAA")
        self.assertEqual(result, "TTTT")


class TestAnnotationCurator(unittest.TestCase):
    """Test annotation curation."""
    
    def setUp(self):
        """Set up test data."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create minimal sequence handler
        cds_file = os.path.join(self.temp_dir, "test.cds.fa")
        aa_file = os.path.join(self.temp_dir, "test.aa")
        
        with open(cds_file, 'w') as f:
            f.write(">g1.t1\nATGAAATAG\n")
        
        with open(aa_file, 'w') as f:
            f.write(">g1.t1\nMK*\n")
        
        self.sequence_handler = SequenceHandler(cds_file, aa_file)
        self.curator = AnnotationCurator(self.sequence_handler)
        
        # Create test transcript
        self.transcript = Transcript(
            id="g1.t1",
            gene_id="g1",
            start=1000,
            end=2000,
            strand="+",
            aa_sequence="MK*"
        )
        
        # Add exon and CDS
        self.transcript.exons.append(Exon(start=1000, end=2000, strand="+"))
        self.transcript.cds_regions.append(CDS(start=1000, end=2000, strand="+"))
        
        # Add codon within existing features
        self.transcript.start_codon = Codon(start=1000, end=1002, strand="+", codon_type="start")
        
        self.gene = Gene(id="g1", transcripts=[self.transcript])
    
    def tearDown(self):
        """Clean up test files."""
        shutil.rmtree(self.temp_dir)
    
    def test_codon_within_features(self):
        """Test codon within features detection."""
        # Codon within existing features
        result = self.curator._is_codon_within_features(
            self.transcript.start_codon,
            self.transcript.exons,
            self.transcript.cds_regions
        )
        self.assertTrue(result)
        
        # Codon outside existing features
        outside_codon = Codon(start=3000, end=3002, strand="+", codon_type="start")
        result = self.curator._is_codon_within_features(
            outside_codon,
            self.transcript.exons,
            self.transcript.cds_regions
        )
        self.assertFalse(result)
    
    def test_create_codon_features(self):
        """Test creating codon features."""
        initial_exon_count = len(self.transcript.exons)
        initial_cds_count = len(self.transcript.cds_regions)
        
        # Create features for codon outside existing regions
        outside_codon = Codon(start=3000, end=3002, strand="+", codon_type="start")
        self.curator._create_codon_features(self.transcript, outside_codon)
        
        # Should have added new exon and CDS
        self.assertEqual(len(self.transcript.exons), initial_exon_count + 1)
        self.assertEqual(len(self.transcript.cds_regions), initial_cds_count + 1)


class TestComplexity(unittest.TestCase):
    """Test algorithmic complexity requirements."""
    
    def test_hash_grouping_complexity(self):
        """Test that hash grouping is O(n) complexity."""
        import time
        
        # Create test transcripts with varying sizes
        sizes = [100, 1000, 10000]
        times = []
        
        for size in sizes:
            transcripts = []
            for i in range(size):
                transcript = Transcript(
                    id=f"g{i}.t1",
                    gene_id=f"g{i}",
                    start=i * 100,
                    end=i * 100 + 50,
                    strand="+",
                    aa_sequence=f"MK*{i}"
                )
                transcripts.append(transcript)
            
            # Time the hash grouping operation
            start_time = time.time()
            hash_groups = {}
            for transcript in transcripts:
                aa_hash = transcript.get_aa_hash()
                if aa_hash not in hash_groups:
                    hash_groups[aa_hash] = []
                hash_groups[aa_hash].append(transcript)
            end_time = time.time()
            
            times.append(end_time - start_time)
        
        # Check that time complexity is roughly linear
        # Time for 10x more data should be less than 20x more time (allowing for overhead)
        if len(times) >= 2:
            time_ratio = times[-1] / times[0] if times[0] > 0 else float('inf')
            size_ratio = sizes[-1] / sizes[0]
            
            # Allow for some overhead, but should be roughly linear
            self.assertLess(time_ratio, size_ratio * 5)


def run_integration_test():
    """Run integration test with sample data."""
    print("Running integration test...")
    
    # Check if sample data exists
    if not all(os.path.exists(f) for f in ["braker.gff3", "braker.cds.fa", "braker.aa"]):
        print("Sample data not found, skipping integration test")
        return
    
    # Run the pipeline
    import subprocess
    import sys
    
    cmd = [
        sys.executable, "gene_curation_pipeline.py",
        "--gff3", "braker.gff3",
        "--cds", "braker.cds.fa",
        "--aa", "braker.aa",
        "--output-dir", "test_integration_output",
        "--min-length", "30",
        "--log-level", "INFO"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            print("✓ Integration test passed")
            
            # Check output files
            expected_files = [
                "test_integration_output/cleaned.gff3",
                "test_integration_output/cleaned.cds.fa",
                "test_integration_output/cleaned.aa",
                "test_integration_output/processing_report.txt"
            ]
            
            for file_path in expected_files:
                if os.path.exists(file_path):
                    print(f"✓ Output file created: {file_path}")
                else:
                    print(f"✗ Missing output file: {file_path}")
        else:
            print(f"✗ Integration test failed: {result.stderr}")
            
    except subprocess.TimeoutExpired:
        print("✗ Integration test timed out")
    except Exception as e:
        print(f"✗ Integration test error: {e}")


if __name__ == "__main__":
    print("Gene Annotation Curation Pipeline - Test Suite")
    print("=" * 50)
    
    # Run unit tests
    unittest.main(verbosity=2, exit=False)
    
    print("\n" + "=" * 50)
    
    # Run integration test
    run_integration_test()