#!/usr/bin/env python3

"""
Unit tests for strict non-overlap representative enforcement.
"""

import os
import sys
import unittest

# Add the parent directory to the path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from gene_curation_pipeline.core.data_structures import Gene, Transcript
from gene_curation_pipeline.core.processors import TranscriptSelector


class TestStrictOverlapEnforcement(unittest.TestCase):
    """Test strict non-overlap behavior when overlap threshold is zero."""

    def _tx(self, tx_id: str, gene_id: str, start: int, end: int, aa_len: int) -> Transcript:
        return Transcript(
            id=tx_id,
            gene_id=gene_id,
            chrom="chr1",
            source="TEST",
            start=start,
            end=end,
            strand="+",
            aa_sequence="M" * aa_len + "*",
            cds_sequence="ATG" * aa_len + "TAA"
        )

    def test_enforcement_uses_non_overlapping_fallback_candidate(self):
        selector = TranscriptSelector(overlap_threshold=0.0, min_aa_length=20)

        g1 = Gene(id="g1")
        g1.add_transcript(self._tx("g1.t1", "g1", 100, 200, 120))

        g2 = Gene(id="g2")
        # Best-by-length candidate overlaps g1.t1
        g2.add_transcript(self._tx("g2.t1", "g2", 150, 250, 110))
        # Alternative candidate does not overlap g1.t1
        g2.add_transcript(self._tx("g2.t2", "g2", 260, 330, 80))

        genes = {"g1": g1, "g2": g2}

        # Initial per-gene picks select the best local representative.
        selector.select_representative(g1)
        selector.select_representative(g2)
        self.assertEqual(g1.representative.id, "g1.t1")
        self.assertEqual(g2.representative.id, "g2.t1")

        selected_count = selector.enforce_non_overlapping_representatives(genes)

        self.assertEqual(selected_count, 2)
        self.assertEqual(g1.representative.id, "g1.t1")
        self.assertEqual(g2.representative.id, "g2.t2")

    def test_enforcement_drops_gene_when_no_non_overlapping_candidate_exists(self):
        selector = TranscriptSelector(overlap_threshold=0.0, min_aa_length=20)

        g1 = Gene(id="g1")
        g1.add_transcript(self._tx("g1.t1", "g1", 100, 220, 120))

        g2 = Gene(id="g2")
        g2.add_transcript(self._tx("g2.t1", "g2", 150, 250, 115))
        g2.add_transcript(self._tx("g2.t2", "g2", 180, 260, 90))

        genes = {"g1": g1, "g2": g2}
        selector.select_representative(g1)
        selector.select_representative(g2)

        selected_count = selector.enforce_non_overlapping_representatives(genes)

        self.assertEqual(selected_count, 1)
        self.assertIsNotNone(g1.representative)
        self.assertIsNone(g2.representative)
        self.assertTrue(
            any("unresolved_spatial_overlap" in t.quality_flags for t in g2.transcripts)
        )


if __name__ == "__main__":
    unittest.main()
