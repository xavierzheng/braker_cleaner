#!/usr/bin/env python3

"""
Unit tests for core data structures.

Tests the fundamental data classes and their methods for correctness
and error handling.
"""

import unittest
import sys
import os

# Add the parent directory to the path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from gene_curation_pipeline.core.data_structures import Gene, Transcript, Exon, CDS, Codon


class TestCodon(unittest.TestCase):
    """Test the Codon data structure."""
    
    def test_valid_codon_creation(self):
        """Test creating valid codons."""
        start_codon = Codon(start=100, end=102, strand="+", codon_type="start", sequence="ATG")
        self.assertEqual(start_codon.start, 100)
        self.assertEqual(start_codon.end, 102)
        self.assertEqual(start_codon.strand, "+")
        self.assertEqual(start_codon.codon_type, "start")
        self.assertEqual(start_codon.sequence, "ATG")
        
        stop_codon = Codon(start=200, end=202, strand="-", codon_type="stop", sequence="TAA")
        self.assertEqual(stop_codon.codon_type, "stop")
        self.assertEqual(stop_codon.strand, "-")
    
    def test_invalid_codon_coordinates(self):
        """Test that invalid coordinates raise ValueError."""
        # Only reversed coordinates should raise error (start > end)
        with self.assertRaises(ValueError):
            Codon(start=100, end=99, strand="+", codon_type="start")
            
        # Single bp coordinates should be valid (start == end)
        # This should NOT raise an error
        codon = Codon(start=100, end=100, strand="+", codon_type="start")
        self.assertEqual(codon.start, 100)
        self.assertEqual(codon.end, 100)
    
    def test_invalid_codon_type(self):
        """Test that invalid codon types raise ValueError."""
        with self.assertRaises(ValueError):
            Codon(start=100, end=102, strand="+", codon_type="invalid")
    
    def test_invalid_strand(self):
        """Test that invalid strands raise ValueError."""
        with self.assertRaises(ValueError):
            Codon(start=100, end=102, strand="x", codon_type="start")


class TestExon(unittest.TestCase):
    """Test the Exon data structure."""
    
    def test_valid_exon_creation(self):
        """Test creating valid exons."""
        exon = Exon(start=1000, end=1200, strand="+", id="exon1")
        self.assertEqual(exon.start, 1000)
        self.assertEqual(exon.end, 1200)
        self.assertEqual(exon.strand, "+")
        self.assertEqual(exon.id, "exon1")
    
    def test_exon_length(self):
        """Test exon length calculation."""
        exon = Exon(start=1000, end=1200, strand="+")
        self.assertEqual(exon.length, 201)  # 1200 - 1000 + 1
    
    def test_exon_overlaps(self):
        """Test exon overlap detection."""
        exon1 = Exon(start=1000, end=1200, strand="+")
        exon2 = Exon(start=1100, end=1300, strand="+")  # Overlaps
        exon3 = Exon(start=1300, end=1400, strand="+")  # No overlap
        
        self.assertTrue(exon1.overlaps_with(exon2))
        self.assertTrue(exon2.overlaps_with(exon1))
        self.assertFalse(exon1.overlaps_with(exon3))
    
    def test_exon_contains_position(self):
        """Test position containment."""
        exon = Exon(start=1000, end=1200, strand="+")
        self.assertTrue(exon.contains_position(1000))   # Start
        self.assertTrue(exon.contains_position(1100))   # Middle
        self.assertTrue(exon.contains_position(1200))   # End
        self.assertFalse(exon.contains_position(999))   # Before
        self.assertFalse(exon.contains_position(1201))  # After


class TestCDS(unittest.TestCase):
    """Test the CDS data structure."""
    
    def test_valid_cds_creation(self):
        """Test creating valid CDS regions."""
        cds = CDS(start=1000, end=1200, strand="+", phase=0, id="cds1")
        self.assertEqual(cds.start, 1000)
        self.assertEqual(cds.end, 1200)
        self.assertEqual(cds.strand, "+")
        self.assertEqual(cds.phase, 0)
        self.assertEqual(cds.id, "cds1")
    
    def test_cds_phase_validation(self):
        """Test that invalid phases raise ValueError."""
        with self.assertRaises(ValueError):
            CDS(start=1000, end=1200, strand="+", phase=3)
        
        with self.assertRaises(ValueError):
            CDS(start=1000, end=1200, strand="+", phase=-1)
    
    def test_cds_length(self):
        """Test CDS length calculation."""
        cds = CDS(start=1000, end=1200, strand="+")
        self.assertEqual(cds.length, 201)


class TestTranscript(unittest.TestCase):
    """Test the Transcript data structure."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.transcript = Transcript(
            id="t1",
            gene_id="g1",
            start=1000,
            end=2000,
            strand="+",
            chrom="chr1",
            source="TEST"
        )
    
    def test_valid_transcript_creation(self):
        """Test creating valid transcripts."""
        self.assertEqual(self.transcript.id, "t1")
        self.assertEqual(self.transcript.gene_id, "g1")
        self.assertEqual(self.transcript.start, 1000)
        self.assertEqual(self.transcript.end, 2000)
        self.assertEqual(self.transcript.strand, "+")
    
    def test_transcript_validation(self):
        """Test transcript validation."""
        with self.assertRaises(ValueError):
            Transcript(id="", gene_id="g1", start=1000, end=2000, strand="+")
        
        with self.assertRaises(ValueError):
            Transcript(id="t1", gene_id="", start=1000, end=2000, strand="+")
        
        with self.assertRaises(ValueError):
            Transcript(id="t1", gene_id="g1", start=2000, end=1000, strand="+")
    
    def test_transcript_properties(self):
        """Test transcript properties."""
        # Add exons
        exon1 = Exon(start=1000, end=1200, strand="+")
        exon2 = Exon(start=1400, end=1600, strand="+")
        self.transcript.exons.extend([exon1, exon2])
        
        # Add CDS regions
        cds1 = CDS(start=1100, end=1200, strand="+")
        cds2 = CDS(start=1400, end=1500, strand="+")
        self.transcript.cds_regions.extend([cds1, cds2])
        
        self.assertEqual(self.transcript.length, 1001)  # 2000 - 1000 + 1
        self.assertEqual(self.transcript.exon_count, 2)
        self.assertEqual(self.transcript.cds_count, 2)
        self.assertEqual(self.transcript.total_exon_length, 402)  # 201 + 201
        self.assertEqual(self.transcript.total_cds_length, 202)   # 101 + 101
    
    def test_hash_methods(self):
        """Test hash generation methods."""
        self.transcript.cds_sequence = "ATGCGATAA"
        self.transcript.aa_sequence = "MR*"
        
        cds_hash = self.transcript.get_cds_hash()
        aa_hash = self.transcript.get_aa_hash()
        
        self.assertIsInstance(cds_hash, str)
        self.assertEqual(len(cds_hash), 32)  # MD5 hash length
        self.assertIsInstance(aa_hash, str)
        self.assertEqual(len(aa_hash), 32)
        
        # Test empty sequences
        empty_transcript = Transcript(id="t2", gene_id="g1", start=1000, end=2000, strand="+")
        self.assertEqual(empty_transcript.get_cds_hash(), "")
        self.assertEqual(empty_transcript.get_aa_hash(), "")
    
    def test_quality_flags(self):
        """Test quality flag management."""
        self.transcript.add_quality_flag("representative")
        self.transcript.add_quality_flag("validated_cds")
        
        self.assertTrue(self.transcript.has_quality_flag("representative"))
        self.assertTrue(self.transcript.has_quality_flag("validated_cds"))
        self.assertFalse(self.transcript.has_quality_flag("nonexistent"))
        
        self.assertEqual(len(self.transcript.quality_flags), 2)
    
    def test_sorted_features(self):
        """Test strand-aware feature sorting."""
        # Test positive strand
        exon1 = Exon(start=1400, end=1500, strand="+")
        exon2 = Exon(start=1000, end=1200, strand="+")
        exon3 = Exon(start=1600, end=1800, strand="+")
        self.transcript.exons.extend([exon1, exon2, exon3])
        
        sorted_exons = self.transcript.get_sorted_exons()
        self.assertEqual([e.start for e in sorted_exons], [1000, 1400, 1600])
        
        # Test negative strand
        neg_transcript = Transcript(id="t2", gene_id="g1", start=1000, end=2000, strand="-")
        neg_transcript.exons.extend([exon1, exon2, exon3])
        
        sorted_neg_exons = neg_transcript.get_sorted_exons()
        self.assertEqual([e.start for e in sorted_neg_exons], [1600, 1400, 1000])
    
    def test_codon_adjacency_detection(self):
        """Test codon adjacency detection."""
        # Add exon
        exon = Exon(start=1200, end=1300, strand="+")
        self.transcript.exons.append(exon)
        
        # Test adjacent codon (1bp gap)
        adjacent_codon = Codon(start=1301, end=1303, strand="+", codon_type="stop")
        self.assertTrue(self.transcript.codon_is_adjacent_to_features(adjacent_codon, gap_threshold=1))
        
        # Test non-adjacent codon (2bp gap)
        non_adjacent_codon = Codon(start=1302, end=1304, strand="+", codon_type="stop")
        self.assertFalse(self.transcript.codon_is_adjacent_to_features(non_adjacent_codon, gap_threshold=1))
        
        # Test with larger threshold
        self.assertTrue(self.transcript.codon_is_adjacent_to_features(non_adjacent_codon, gap_threshold=5))
    
    def test_codon_coverage_detection(self):
        """Test codon coverage detection."""
        # Add exon and CDS that cover a region
        exon = Exon(start=1200, end=1300, strand="+")
        cds = CDS(start=1200, end=1300, strand="+")
        self.transcript.exons.append(exon)
        self.transcript.cds_regions.append(cds)
        
        # Test covered codon
        covered_codon = Codon(start=1250, end=1252, strand="+", codon_type="start")
        self.assertTrue(self.transcript.codon_is_covered_by_features(covered_codon))
        
        # Test partially covered codon (only by exon, not CDS)
        partial_exon = Exon(start=1400, end=1500, strand="+")
        self.transcript.exons.append(partial_exon)
        
        partial_codon = Codon(start=1450, end=1452, strand="+", codon_type="stop")
        self.assertFalse(self.transcript.codon_is_covered_by_features(partial_codon))  # No CDS coverage
        
        # Test uncovered codon
        uncovered_codon = Codon(start=1600, end=1602, strand="+", codon_type="stop")
        self.assertFalse(self.transcript.codon_is_covered_by_features(uncovered_codon))


class TestGene(unittest.TestCase):
    """Test the Gene data structure."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.gene = Gene(id="g1")
        self.transcript1 = Transcript(id="t1", gene_id="g1", start=1000, end=2000, strand="+")
        self.transcript2 = Transcript(id="t2", gene_id="g1", start=1500, end=2500, strand="+")
    
    def test_valid_gene_creation(self):
        """Test creating valid genes."""
        self.assertEqual(self.gene.id, "g1")
        self.assertEqual(self.gene.transcript_count, 0)
        self.assertFalse(self.gene.has_representative)
    
    def test_gene_validation(self):
        """Test gene validation."""
        with self.assertRaises(ValueError):
            Gene(id="")
    
    def test_transcript_management(self):
        """Test transcript addition and retrieval."""
        self.gene.add_transcript(self.transcript1)
        self.assertEqual(self.gene.transcript_count, 1)
        
        retrieved = self.gene.get_transcript_by_id("t1")
        self.assertEqual(retrieved, self.transcript1)
        
        self.gene.add_transcript(self.transcript2)
        self.assertEqual(self.gene.transcript_count, 2)
        
        # Test duplicate transcript
        with self.assertRaises(ValueError):
            self.gene.add_transcript(self.transcript1)
        
        # Test wrong gene ID
        wrong_transcript = Transcript(id="t3", gene_id="g2", start=1000, end=2000, strand="+")
        with self.assertRaises(ValueError):
            self.gene.add_transcript(wrong_transcript)
    
    def test_representative_management(self):
        """Test representative transcript management."""
        self.gene.add_transcript(self.transcript1)
        
        self.gene.set_representative(self.transcript1)
        self.assertTrue(self.gene.has_representative)
        self.assertEqual(self.gene.representative, self.transcript1)
        self.assertTrue(self.transcript1.has_quality_flag("representative"))
        
        # Test invalid representative
        with self.assertRaises(ValueError):
            self.gene.set_representative(self.transcript2)  # Not in gene
    
    def test_gene_boundaries(self):
        """Test gene boundary calculation."""
        self.gene.add_transcript(self.transcript1)
        self.gene.add_transcript(self.transcript2)
        
        start, end = self.gene.get_gene_boundaries()
        self.assertEqual(start, 1000)  # Min transcript start
        self.assertEqual(end, 2500)    # Max transcript end
        
        # Test representative boundaries
        self.gene.set_representative(self.transcript1)
        exon1 = Exon(start=900, end=1200, strand="+")
        exon2 = Exon(start=1800, end=2100, strand="+")
        self.transcript1.exons.extend([exon1, exon2])
        
        rep_start, rep_end = self.gene.get_representative_boundaries()
        self.assertEqual(rep_start, 900)   # Min representative exon start
        self.assertEqual(rep_end, 2100)    # Max representative exon end
    
    def test_transcript_filtering(self):
        """Test transcript filtering by quality flags."""
        self.transcript1.add_quality_flag("validated")
        self.transcript1.add_quality_flag("representative")
        self.transcript2.add_quality_flag("validated")
        
        self.gene.add_transcript(self.transcript1)
        self.gene.add_transcript(self.transcript2)
        
        # Filter by single flag
        validated = self.gene.get_transcript_by_quality_flags({"validated"})
        self.assertEqual(len(validated), 2)
        
        # Filter by multiple flags
        rep_validated = self.gene.get_transcript_by_quality_flags({"validated", "representative"})
        self.assertEqual(len(rep_validated), 1)
        self.assertEqual(rep_validated[0], self.transcript1)
        
        # Filter by non-existent flag
        missing = self.gene.get_transcript_by_quality_flags({"nonexistent"})
        self.assertEqual(len(missing), 0)


if __name__ == '__main__':
    unittest.main()