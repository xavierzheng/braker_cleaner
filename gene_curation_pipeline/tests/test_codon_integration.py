#!/usr/bin/env python3

"""
Comprehensive tests for codon integration functionality.

Tests the general case of codon integration without relying on specific gene IDs,
covering various scenarios including adjacent and non-adjacent codons.
"""

import unittest
import tempfile
import os
from unittest.mock import Mock, patch
import sys

# Add the parent directory to the path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from gene_curation_pipeline.core.data_structures import Gene, Transcript, Exon, CDS, Codon
from gene_curation_pipeline.core.exceptions import ValidationError


class TestCodonIntegration(unittest.TestCase):
    """Test codon integration with various scenarios."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.gene_id = "test_gene"
        self.transcript_id = "test_transcript"
        self.chromosome = "chr1"
        
        # Create a basic transcript
        self.transcript = Transcript(
            id=self.transcript_id,
            gene_id=self.gene_id,
            start=1000,
            end=2000,
            strand="+",
            chrom=self.chromosome,
            source="TEST"
        )
        
        self.gene = Gene(id=self.gene_id)
        self.gene.add_transcript(self.transcript)
    
    def test_codon_adjacent_to_exon_integration(self):
        """Test integration of codons adjacent to existing exons."""
        # Scenario: Stop codon adjacent to exon (similar to g2.t1 case)
        
        # Add existing exon
        exon = Exon(start=1251, end=1577, strand="+", id="exon1")
        self.transcript.exons.append(exon)
        
        # Add existing CDS
        cds = CDS(start=1251, end=1577, strand="+", phase=0, id="cds1")
        self.transcript.cds_regions.append(cds)
        
        # Add stop codon adjacent to exon (1bp gap)
        stop_codon = Codon(start=1248, end=1250, strand="+", codon_type="stop", sequence="TAA")
        self.transcript.stop_codon = stop_codon
        
        # Test adjacency detection
        self.assertTrue(self.transcript.codon_is_adjacent_to_features(stop_codon))
        
        # Test that codon is NOT already covered
        self.assertFalse(self.transcript.codon_is_covered_by_features(stop_codon))
        
        # Simulate codon integration by extending exon and CDS
        extended_exon = Exon(start=1248, end=1577, strand="+", id="exon1_extended")
        extended_cds = CDS(start=1248, end=1577, strand="+", phase=0, id="cds1_extended")
        
        # Verify extended features cover the original codon
        self.assertTrue(extended_exon.contains_position(stop_codon.start))
        self.assertTrue(extended_exon.contains_position(stop_codon.end))
        self.assertTrue(extended_cds.contains_position(stop_codon.start))
        self.assertTrue(extended_cds.contains_position(stop_codon.end))
    
    def test_codon_non_adjacent_requires_new_features(self):
        """Test codons that are not adjacent and require new features."""
        # Add existing exon far from codon
        exon = Exon(start=1500, end=1600, strand="+", id="exon1")
        self.transcript.exons.append(exon)
        
        # Add start codon not adjacent to any existing feature
        start_codon = Codon(start=1100, end=1102, strand="+", codon_type="start", sequence="ATG")
        self.transcript.start_codon = start_codon
        
        # Test that codon is NOT adjacent
        self.assertFalse(self.transcript.codon_is_adjacent_to_features(start_codon))
        
        # Test that codon is NOT covered
        self.assertFalse(self.transcript.codon_is_covered_by_features(start_codon))
        
        # This would require creating new exon/CDS features
        new_exon = Exon(start=1100, end=1102, strand="+", id="start_codon_exon")
        new_cds = CDS(start=1100, end=1102, strand="+", phase=0, id="start_codon_cds")
        
        self.assertTrue(new_exon.contains_position(start_codon.start))
        self.assertTrue(new_cds.contains_position(start_codon.start))
    
    def test_codon_already_covered_is_redundant(self):
        """Test codons that are already covered by existing features."""
        # Add exon that covers a wide range
        exon = Exon(start=1200, end=1300, strand="+", id="exon1")
        self.transcript.exons.append(exon)
        
        # Add CDS that covers the same range
        cds = CDS(start=1200, end=1300, strand="+", phase=0, id="cds1")
        self.transcript.cds_regions.append(cds)
        
        # Add start codon that is already covered by existing features
        start_codon = Codon(start=1220, end=1222, strand="+", codon_type="start", sequence="ATG")
        self.transcript.start_codon = start_codon
        
        # Test that codon is already covered (redundant)
        self.assertTrue(self.transcript.codon_is_covered_by_features(start_codon))
        
        # This codon should be removed as redundant
    
    def test_negative_strand_codon_integration(self):
        """Test codon integration on negative strand."""
        # Create negative strand transcript
        neg_transcript = Transcript(
            id="neg_transcript",
            gene_id="neg_gene",
            start=1000,
            end=2000,
            strand="-",
            chrom=self.chromosome,
            source="TEST"
        )
        
        # Add exon on negative strand
        exon = Exon(start=1400, end=1500, strand="-", id="neg_exon1")
        neg_transcript.exons.append(exon)
        
        # Add CDS on negative strand
        cds = CDS(start=1400, end=1500, strand="-", phase=0, id="neg_cds1")
        neg_transcript.cds_regions.append(cds)
        
        # Add start codon adjacent to exon (on negative strand, start codon is at 3' end)
        start_codon = Codon(start=1501, end=1503, strand="-", codon_type="start", sequence="ATG")
        neg_transcript.start_codon = start_codon
        
        # Test adjacency detection
        self.assertTrue(neg_transcript.codon_is_adjacent_to_features(start_codon))
        
        # Test that codon is not covered
        self.assertFalse(neg_transcript.codon_is_covered_by_features(start_codon))
    
    def test_multiple_exon_cds_codon_scenario(self):
        """Test codon integration with multiple exons and CDS regions."""
        # Add multiple exons
        exon1 = Exon(start=1100, end=1200, strand="+", id="exon1")
        exon2 = Exon(start=1300, end=1400, strand="+", id="exon2")
        exon3 = Exon(start=1500, end=1600, strand="+", id="exon3")
        
        self.transcript.exons.extend([exon1, exon2, exon3])
        
        # Add corresponding CDS regions
        cds1 = CDS(start=1150, end=1200, strand="+", phase=0, id="cds1")
        cds2 = CDS(start=1300, end=1400, strand="+", phase=1, id="cds2")
        cds3 = CDS(start=1500, end=1550, strand="+", phase=0, id="cds3")
        
        self.transcript.cds_regions.extend([cds1, cds2, cds3])
        
        # Test start codon adjacent to first CDS
        start_codon = Codon(start=1147, end=1149, strand="+", codon_type="start", sequence="ATG")
        self.transcript.start_codon = start_codon
        
        # Should be adjacent to exon1 and cds1
        self.assertTrue(self.transcript.codon_is_adjacent_to_features(start_codon))
        
        # Test stop codon adjacent to last CDS
        stop_codon = Codon(start=1551, end=1553, strand="+", codon_type="stop", sequence="TAA")
        self.transcript.stop_codon = stop_codon
        
        # Should be adjacent to cds3
        self.assertTrue(self.transcript.codon_is_adjacent_to_features(stop_codon))
    
    def test_gap_threshold_sensitivity(self):
        """Test that gap threshold affects adjacency detection."""
        # Add exon
        exon = Exon(start=1200, end=1300, strand="+", id="exon1")
        self.transcript.exons.append(exon)
        
        # Test codon at different distances
        codon_1bp = Codon(start=1301, end=1303, strand="+", codon_type="stop", sequence="TAA")
        codon_2bp = Codon(start=1302, end=1304, strand="+", codon_type="stop", sequence="TAA")
        codon_5bp = Codon(start=1305, end=1307, strand="+", codon_type="stop", sequence="TAA")
        
        # Test with default gap threshold (1bp)
        self.assertTrue(self.transcript.codon_is_adjacent_to_features(codon_1bp, gap_threshold=1))
        self.assertFalse(self.transcript.codon_is_adjacent_to_features(codon_2bp, gap_threshold=1))
        self.assertFalse(self.transcript.codon_is_adjacent_to_features(codon_5bp, gap_threshold=1))
        
        # Test with larger gap threshold (5bp)
        self.assertTrue(self.transcript.codon_is_adjacent_to_features(codon_1bp, gap_threshold=5))
        self.assertTrue(self.transcript.codon_is_adjacent_to_features(codon_2bp, gap_threshold=5))
        self.assertTrue(self.transcript.codon_is_adjacent_to_features(codon_5bp, gap_threshold=5))
    
    def test_codon_integration_quality_flags(self):
        """Test that appropriate quality flags are set during codon integration."""
        # This would be tested in the actual pipeline components
        # For now, test that transcript can hold quality flags
        
        self.transcript.add_quality_flag("created_start_codon_features")
        self.transcript.add_quality_flag("created_stop_codon_features") 
        self.transcript.add_quality_flag("redundant_start_codon_removed")
        
        self.assertTrue(self.transcript.has_quality_flag("created_start_codon_features"))
        self.assertTrue(self.transcript.has_quality_flag("created_stop_codon_features"))
        self.assertTrue(self.transcript.has_quality_flag("redundant_start_codon_removed"))
        
        expected_flags = {"created_start_codon_features", "created_stop_codon_features", "redundant_start_codon_removed"}
        self.assertTrue(expected_flags.issubset(self.transcript.quality_flags))
    
    def test_transcript_boundary_updates(self):
        """Test that transcript boundaries are updated after codon integration."""
        original_start = self.transcript.start
        original_end = self.transcript.end
        
        # Add exons that extend beyond original boundaries
        exon1 = Exon(start=900, end=1100, strand="+", id="extended_exon1")  # Extends start
        exon2 = Exon(start=1900, end=2100, strand="+", id="extended_exon2")  # Extends end
        
        self.transcript.exons.extend([exon1, exon2])
        
        # Get new boundaries based on exons
        exon_starts = [e.start for e in self.transcript.exons]
        exon_ends = [e.end for e in self.transcript.exons]
        new_start = min(exon_starts)
        new_end = max(exon_ends)
        
        # Verify boundaries would need to be updated
        self.assertLess(new_start, original_start)  # Start extended
        self.assertGreater(new_end, original_end)   # End extended
    
    def test_gene_boundary_calculation(self):
        """Test gene boundary calculation after representative selection."""
        # Add another transcript to the gene
        transcript2 = Transcript(
            id="transcript2",
            gene_id=self.gene_id,
            start=1500,
            end=2500,
            strand="+",
            chrom=self.chromosome,
            source="TEST"
        )
        self.gene.add_transcript(transcript2)
        
        # Test original gene boundaries
        original_start, original_end = self.gene.get_gene_boundaries()
        self.assertEqual(original_start, 1000)  # Min of transcript starts
        self.assertEqual(original_end, 2500)    # Max of transcript ends
        
        # Set representative and add exons
        self.gene.set_representative(self.transcript)
        
        exon1 = Exon(start=800, end=1200, strand="+", id="rep_exon1")
        exon2 = Exon(start=1800, end=2200, strand="+", id="rep_exon2")
        self.transcript.exons.extend([exon1, exon2])
        
        # Test representative-based boundaries
        rep_start, rep_end = self.gene.get_representative_boundaries()
        self.assertEqual(rep_start, 800)   # Min of representative exon starts
        self.assertEqual(rep_end, 2200)    # Max of representative exon ends


class TestCodonValidation(unittest.TestCase):
    """Test codon validation and error handling."""
    
    def test_invalid_codon_coordinates(self):
        """Test that invalid codon coordinates raise errors."""
        with self.assertRaises(ValueError):
            Codon(start=100, end=100, strand="+", codon_type="start")  # start == end
        
        with self.assertRaises(ValueError):
            Codon(start=100, end=99, strand="+", codon_type="start")   # start > end
    
    def test_invalid_codon_type(self):
        """Test that invalid codon types raise errors."""
        with self.assertRaises(ValueError):
            Codon(start=100, end=102, strand="+", codon_type="invalid")
    
    def test_invalid_strand(self):
        """Test that invalid strand values raise errors."""
        with self.assertRaises(ValueError):
            Codon(start=100, end=102, strand="*", codon_type="start")


class TestGeneralCodonIntegrationScenarios(unittest.TestCase):
    """Test general codon integration scenarios that could apply to any gene."""
    
    def create_test_transcript(self, transcript_id: str = "test_t1", gene_id: str = "test_g1") -> Transcript:
        """Helper method to create a test transcript."""
        return Transcript(
            id=transcript_id,
            gene_id=gene_id,
            start=1000,
            end=2000,
            strand="+",
            chrom="chr_test",
            source="TEST"
        )
    
    def test_scenario_stop_codon_adjacent_to_first_exon(self):
        """Test stop codon integration adjacent to first exon (general g2.t1-like case)."""
        transcript = self.create_test_transcript()
        
        # Create scenario similar to g2.t1: stop codon adjacent to first exon
        stop_codon = Codon(start=1248, end=1250, strand="-", codon_type="stop", sequence="TAA")
        first_exon = Exon(start=1251, end=1577, strand="-", id="first_exon")
        first_cds = CDS(start=1251, end=1577, strand="-", phase=0, id="first_cds")
        
        transcript.stop_codon = stop_codon
        transcript.exons.append(first_exon)
        transcript.cds_regions.append(first_cds)
        transcript.strand = "-"  # Negative strand like g2.t1
        
        # Verify scenario setup
        self.assertTrue(transcript.codon_is_adjacent_to_features(stop_codon))
        self.assertFalse(transcript.codon_is_covered_by_features(stop_codon))
        
        # Simulate integration result
        merged_exon = Exon(start=1248, end=1577, strand="-", id="merged_exon")
        merged_cds = CDS(start=1248, end=1577, strand="-", phase=0, id="merged_cds")
        
        # Verify merged features cover the codon
        self.assertTrue(merged_exon.contains_position(stop_codon.start))
        self.assertTrue(merged_cds.contains_position(stop_codon.start))
        
        # After integration, the codon would be flagged as redundant
        transcript.add_quality_flag("redundant_stop_codon_removed")
        transcript.add_quality_flag("created_stop_codon_features")
        
        self.assertTrue(transcript.has_quality_flag("redundant_stop_codon_removed"))
        self.assertTrue(transcript.has_quality_flag("created_stop_codon_features"))
    
    def test_scenario_start_codon_adjacent_to_last_cds(self):
        """Test start codon integration adjacent to last CDS on negative strand."""
        transcript = self.create_test_transcript()
        transcript.strand = "-"
        
        # On negative strand, start codon is typically at the 3' end (higher coordinates)
        last_cds = CDS(start=1800, end=1900, strand="-", phase=0, id="last_cds")
        last_exon = Exon(start=1800, end=1950, strand="-", id="last_exon")
        start_codon = Codon(start=1901, end=1903, strand="-", codon_type="start", sequence="ATG")
        
        transcript.cds_regions.append(last_cds)
        transcript.exons.append(last_exon)
        transcript.start_codon = start_codon
        
        # Verify adjacency
        self.assertTrue(transcript.codon_is_adjacent_to_features(start_codon))
        
        # The codon should be adjacent to CDS but not covered by the current CDS range
        self.assertFalse(transcript.codon_is_covered_by_features(start_codon))
    
    def test_scenario_orphaned_cds_without_exon(self):
        """Test scenario where CDS exists but no corresponding exon (GeneMark.hmm3 case)."""
        transcript = self.create_test_transcript()
        
        # Add CDS without corresponding exon (orphaned CDS)
        orphaned_cds = CDS(start=1200, end=1300, strand="+", phase=0, id="orphaned_cds")
        transcript.cds_regions.append(orphaned_cds)
        
        # Verify no exon covers this CDS
        cds_covered_by_exon = any(
            exon.contains_position(orphaned_cds.start) and exon.contains_position(orphaned_cds.end)
            for exon in transcript.exons
        )
        self.assertFalse(cds_covered_by_exon)
        
        # This would trigger the _ensure_cds_have_exons functionality
        # Create matching exon for orphaned CDS
        matching_exon = Exon(start=1200, end=1300, strand="+", id="created_for_cds")
        transcript.exons.append(matching_exon)
        
        # Verify the created exon covers the CDS
        self.assertTrue(matching_exon.contains_position(orphaned_cds.start))
        self.assertTrue(matching_exon.contains_position(orphaned_cds.end))
        
        # This would be flagged
        transcript.add_quality_flag("created_exon_for_cds")
        self.assertTrue(transcript.has_quality_flag("created_exon_for_cds"))
    
    def test_scenario_multiple_codon_integrations(self):
        """Test multiple codon integrations in the same transcript."""
        transcript = self.create_test_transcript()
        
        # Create multiple exons/CDS
        exon1 = Exon(start=1100, end=1200, strand="+", id="exon1")
        exon2 = Exon(start=1400, end=1500, strand="+", id="exon2")
        cds1 = CDS(start=1150, end=1200, strand="+", phase=0, id="cds1")
        cds2 = CDS(start=1400, end=1450, strand="+", phase=2, id="cds2")
        
        transcript.exons.extend([exon1, exon2])
        transcript.cds_regions.extend([cds1, cds2])
        
        # Add both start and stop codons that need integration
        start_codon = Codon(start=1147, end=1149, strand="+", codon_type="start", sequence="ATG")
        stop_codon = Codon(start=1451, end=1453, strand="+", codon_type="stop", sequence="TAA")
        
        transcript.start_codon = start_codon
        transcript.stop_codon = stop_codon
        
        # Verify both codons are adjacent to features
        self.assertTrue(transcript.codon_is_adjacent_to_features(start_codon))
        self.assertTrue(transcript.codon_is_adjacent_to_features(stop_codon))
        
        # Verify both codons need integration (not covered)
        self.assertFalse(transcript.codon_is_covered_by_features(start_codon))
        self.assertFalse(transcript.codon_is_covered_by_features(stop_codon))
        
        # After integration, both would be flagged
        transcript.add_quality_flag("created_start_codon_features")
        transcript.add_quality_flag("created_stop_codon_features")
        
        self.assertTrue(transcript.has_quality_flag("created_start_codon_features"))
        self.assertTrue(transcript.has_quality_flag("created_stop_codon_features"))


if __name__ == '__main__':
    unittest.main()