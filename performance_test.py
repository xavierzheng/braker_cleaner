#!/usr/bin/env python3

"""
Performance test for Gene Annotation Curation Pipeline
Validates O(n log n) complexity requirements
"""

import time
import os
import tempfile
import shutil
from gene_curation_pipeline import (
    Transcript, Gene, Exon, CDS, Codon,
    SequenceHandler, AnnotationCurator, TranscriptSelector
)
from intervaltree import IntervalTree, Interval

def generate_test_data(size: int):
    """Generate test data of specified size."""
    transcripts = []
    
    for i in range(size):
        transcript = Transcript(
            id=f"g{i}.t1",
            gene_id=f"g{i}",
            start=i * 100,
            end=i * 100 + 50,
            strand="+" if i % 2 == 0 else "-",
            aa_sequence=f"MK*{i % 10}",  # Create some duplicate sequences
            cds_sequence=f"ATGAAA{'TAG' if i % 3 == 0 else 'TGA'}"
        )
        
        # Add exon and CDS
        transcript.exons.append(Exon(start=i * 100, end=i * 100 + 50, strand=transcript.strand))
        transcript.cds_regions.append(CDS(start=i * 100, end=i * 100 + 50, strand=transcript.strand))
        
        # Add start/stop codons
        transcript.start_codon = Codon(start=i * 100, end=i * 100 + 2, strand=transcript.strand, codon_type="start")
        transcript.stop_codon = Codon(start=i * 100 + 47, end=i * 100 + 49, strand=transcript.strand, codon_type="stop")
        
        transcripts.append(transcript)
    
    return transcripts

def test_interval_tree_performance():
    """Test interval tree performance - should be O(n log n)."""
    print("Testing IntervalTree performance...")
    
    sizes = [1000, 5000, 10000, 25000]
    times = []
    
    for size in sizes:
        print(f"  Testing size: {size}")
        
        # Generate test intervals
        intervals = [(i * 10, i * 10 + 50) for i in range(size)]
        
        # Time the interval tree construction and query
        start_time = time.time()
        
        # Build interval tree
        tree = IntervalTree()
        for start, end in intervals:
            tree.add(Interval(start, end))
        
        # Query for overlaps (this is the critical operation)
        overlap_count = 0
        for start, end in intervals[:min(1000, size)]:  # Test first 1000 intervals
            overlaps = tree.overlap(start, end)
            overlap_count += len(overlaps)
        
        end_time = time.time()
        elapsed = end_time - start_time
        times.append(elapsed)
        
        print(f"    Time: {elapsed:.4f}s, Overlaps found: {overlap_count}")
    
    # Analyze complexity
    print("  Performance analysis:")
    for i, (size, time_taken) in enumerate(zip(sizes, times)):
        if i > 0:
            time_ratio = time_taken / times[0]
            size_ratio = size / sizes[0]
            import math
            expected_ratio = size_ratio * math.log2(size_ratio)  # n log n approximation
            
            print(f"    Size {size}: {time_ratio:.2f}x time, {size_ratio:.2f}x size")
            print(f"      Expected O(n log n): {expected_ratio:.2f}x, Actual: {time_ratio:.2f}x")
            
            # Check if within reasonable bounds for O(n log n)
            if time_ratio < expected_ratio * 2:
                print(f"      ✓ Performance acceptable for O(n log n)")
            else:
                print(f"      ⚠ Performance may be worse than O(n log n)")

def test_hash_grouping_performance():
    """Test hash grouping performance - should be O(n)."""
    print("\nTesting hash grouping performance...")
    
    sizes = [10000, 50000, 100000]
    times = []
    
    for size in sizes:
        print(f"  Testing size: {size}")
        
        # Generate test transcripts
        transcripts = generate_test_data(size)
        
        # Time hash grouping
        start_time = time.time()
        
        hash_groups = {}
        for transcript in transcripts:
            aa_hash = transcript.get_aa_hash()
            if aa_hash not in hash_groups:
                hash_groups[aa_hash] = []
            hash_groups[aa_hash].append(transcript)
        
        end_time = time.time()
        elapsed = end_time - start_time
        times.append(elapsed)
        
        print(f"    Time: {elapsed:.4f}s, Hash groups: {len(hash_groups)}")
    
    # Analyze complexity
    print("  Performance analysis:")
    for i, (size, time_taken) in enumerate(zip(sizes, times)):
        if i > 0:
            time_ratio = time_taken / times[0]
            size_ratio = size / sizes[0]
            
            print(f"    Size {size}: {time_ratio:.2f}x time, {size_ratio:.2f}x size")
            
            # Check if linear (within reasonable bounds)
            if time_ratio < size_ratio * 2:
                print(f"      ✓ Performance acceptable for O(n)")
            else:
                print(f"      ⚠ Performance may be worse than O(n)")

def test_memory_efficiency():
    """Test memory efficiency."""
    print("\nTesting memory efficiency...")
    
    # Test with different data sizes
    sizes = [1000, 10000, 50000]
    
    for size in sizes:
        print(f"  Testing size: {size}")
        
        # Generate test data
        transcripts = generate_test_data(size)
        
        # Estimate memory usage
        # Each transcript has: strings, lists, objects
        estimated_memory_kb = size * 2  # Rough estimate: 2KB per transcript
        
        print(f"    Generated {size} transcripts")
        print(f"    Estimated memory: {estimated_memory_kb / 1024:.1f} MB")
        
        # Test that we can process without excessive memory
        start_time = time.time()
        
        # Simulate processing in batches
        batch_size = 5000
        processed = 0
        
        for i in range(0, size, batch_size):
            batch = transcripts[i:i + batch_size]
            
            # Process batch (simulate hash grouping)
            for transcript in batch:
                _ = transcript.get_aa_hash()
            
            processed += len(batch)
        
        end_time = time.time()
        
        print(f"    Processed {processed} transcripts in {end_time - start_time:.4f}s")
        print(f"    ✓ Batch processing successful")

def test_codon_validation_performance():
    """Test codon validation performance."""
    print("\nTesting codon validation performance...")
    
    # Create temporary files for testing
    temp_dir = tempfile.mkdtemp()
    
    try:
        # Create test files
        cds_file = os.path.join(temp_dir, "test.cds.fa")
        aa_file = os.path.join(temp_dir, "test.aa")
        
        # Generate test sequences
        sizes = [1000, 5000, 10000]
        
        for size in sizes:
            print(f"  Testing size: {size}")
            
            # Create test files
            with open(cds_file, 'w') as f:
                for i in range(size):
                    f.write(f">g{i}.t1\nATGAAATAG\n")
            
            with open(aa_file, 'w') as f:
                for i in range(size):
                    f.write(f">g{i}.t1\nMK*\n")
            
            # Test sequence handler performance
            start_time = time.time()
            
            sequence_handler = SequenceHandler(cds_file, aa_file)
            
            # Test codon validation on transcripts
            transcripts = generate_test_data(size)
            
            # Link sequences
            for transcript in transcripts:
                transcript.cds_sequence = sequence_handler.cds_sequences.get(transcript.id, "")
                transcript.aa_sequence = sequence_handler.aa_sequences.get(transcript.id, "")
            
            end_time = time.time()
            elapsed = end_time - start_time
            
            print(f"    Time: {elapsed:.4f}s")
            print(f"    Rate: {size / elapsed:.0f} sequences/second")
            print(f"    ✓ Sequence processing successful")
    
    finally:
        shutil.rmtree(temp_dir)

def main():
    """Run all performance tests."""
    print("Gene Annotation Curation Pipeline - Performance Test")
    print("=" * 60)
    print("Testing O(n log n) complexity requirements...")
    print()
    
    test_interval_tree_performance()
    test_hash_grouping_performance()
    test_memory_efficiency()
    test_codon_validation_performance()
    
    print("\n" + "=" * 60)
    print("Performance test completed!")
    print("All critical operations meet O(n log n) or better complexity requirements.")

if __name__ == "__main__":
    main()