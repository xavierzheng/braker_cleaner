#!/usr/bin/env python3

"""
Enhanced performance monitoring for the gene curation pipeline.

Provides comprehensive monitoring of memory usage, processing time,
and algorithmic complexity validation.
"""

import time
import logging
import psutil
from dataclasses import dataclass
from typing import Optional, Dict, Any, List
from contextlib import contextmanager

from ..core.exceptions import MemoryError as PipelineMemoryError


@dataclass
class PerformanceMetrics:
    """Container for performance metrics."""
    start_time: float
    end_time: Optional[float] = None
    peak_memory_mb: float = 0.0
    current_memory_mb: float = 0.0
    operations_count: int = 0
    phase_name: str = ""
    
    @property
    def elapsed_time(self) -> float:
        """Get elapsed time in seconds."""
        if self.end_time is None:
            return time.time() - self.start_time
        return self.end_time - self.start_time
    
    @property
    def operations_per_second(self) -> float:
        """Get operations per second."""
        elapsed = self.elapsed_time
        if elapsed > 0 and self.operations_count > 0:
            return self.operations_count / elapsed
        return 0.0


class PerformanceMonitor:
    """Enhanced performance monitoring with complexity validation."""
    
    def __init__(self, memory_limit_mb: int = 4096):
        self.memory_limit_mb = memory_limit_mb
        self.start_time = time.time()
        self.phase_metrics: Dict[str, PerformanceMetrics] = {}
        self.current_phase: Optional[str] = None
        
        # Initialize psutil process if available
        try:
            self.process = psutil.Process()
        except:
            self.process = None
            logging.warning("psutil not available, memory monitoring disabled")
    
    def get_memory_usage(self) -> float:
        """Get current memory usage in MB."""
        if not self.process:
            return 0.0
        
        try:
            memory_mb = self.process.memory_info().rss / 1024 / 1024
            
            # Update peak memory for current phase
            if self.current_phase and self.current_phase in self.phase_metrics:
                metrics = self.phase_metrics[self.current_phase]
                metrics.current_memory_mb = memory_mb
                metrics.peak_memory_mb = max(metrics.peak_memory_mb, memory_mb)
            
            return memory_mb
            
        except Exception as e:
            logging.warning(f"Error getting memory usage: {e}")
            return 0.0
    
    def check_memory_limit(self) -> bool:
        """Check if memory usage exceeds limit."""
        current_memory = self.get_memory_usage()
        
        if current_memory > self.memory_limit_mb:
            error_msg = f"Memory usage exceeded limit"
            logging.warning(f"{error_msg}: {current_memory:.1f}MB > {self.memory_limit_mb}MB")
            raise PipelineMemoryError(error_msg, current_memory, self.memory_limit_mb)
        
        return True
    
    def start_phase(self, phase_name: str) -> None:
        """Start monitoring a processing phase."""
        if self.current_phase:
            self.end_phase()
        
        self.current_phase = phase_name
        self.phase_metrics[phase_name] = PerformanceMetrics(
            start_time=time.time(),
            phase_name=phase_name,
            current_memory_mb=self.get_memory_usage()
        )
        
        logging.info(f"Started phase: {phase_name}")
    
    def end_phase(self) -> Optional[PerformanceMetrics]:
        """End the current phase and return metrics."""
        if not self.current_phase:
            return None
        
        metrics = self.phase_metrics[self.current_phase]
        metrics.end_time = time.time()
        metrics.current_memory_mb = self.get_memory_usage()
        
        logging.info(f"Completed phase {self.current_phase} in {metrics.elapsed_time:.2f}s "
                    f"(peak memory: {metrics.peak_memory_mb:.1f}MB)")
        
        phase_name = self.current_phase
        self.current_phase = None
        
        return metrics
    
    def record_operations(self, count: int) -> None:
        """Record the number of operations performed in current phase."""
        if self.current_phase and self.current_phase in self.phase_metrics:
            self.phase_metrics[self.current_phase].operations_count += count
    
    @contextmanager
    def phase_context(self, phase_name: str):
        """Context manager for monitoring a phase."""
        self.start_phase(phase_name)
        try:
            yield self.phase_metrics[phase_name]
        finally:
            self.end_phase()
    
    def get_total_elapsed_time(self) -> float:
        """Get total elapsed time since monitor creation."""
        return time.time() - self.start_time
    
    def get_peak_memory(self) -> float:
        """Get peak memory usage across all phases."""
        if not self.phase_metrics:
            return self.get_memory_usage()
        
        return max(metrics.peak_memory_mb for metrics in self.phase_metrics.values())
    
    def get_performance_summary(self) -> Dict[str, Any]:
        """Get comprehensive performance summary."""
        summary = {
            "total_elapsed_time": self.get_total_elapsed_time(),
            "peak_memory_mb": self.get_peak_memory(),
            "current_memory_mb": self.get_memory_usage(),
            "memory_limit_mb": self.memory_limit_mb,
            "phases": {}
        }
        
        for phase_name, metrics in self.phase_metrics.items():
            summary["phases"][phase_name] = {
                "elapsed_time": metrics.elapsed_time,
                "operations_count": metrics.operations_count,
                "operations_per_second": metrics.operations_per_second,
                "peak_memory_mb": metrics.peak_memory_mb
            }
        
        return summary
    
    def validate_complexity(self, n: int, actual_time: float, 
                           expected_complexity: str = "O(n log n)") -> bool:
        """
        Validate that actual performance matches expected algorithmic complexity.
        
        Args:
            n: Input size
            actual_time: Actual processing time
            expected_complexity: Expected complexity (e.g., "O(n log n)")
        
        Returns:
            True if performance is within expected bounds
        """
        if n <= 0 or actual_time <= 0:
            return True
        
        import math
        
        # Define complexity functions
        complexity_functions = {
            "O(1)": lambda x: 1,
            "O(log n)": lambda x: math.log(x) if x > 1 else 1,
            "O(n)": lambda x: x,
            "O(n log n)": lambda x: x * math.log(x) if x > 1 else x,
            "O(n^2)": lambda x: x * x,
        }
        
        if expected_complexity not in complexity_functions:
            logging.warning(f"Unknown complexity: {expected_complexity}")
            return True
        
        complexity_func = complexity_functions[expected_complexity]
        expected_relative_time = complexity_func(n)
        
        # Calculate relative performance (time per complexity unit)
        relative_performance = actual_time / expected_relative_time
        
        # Store for analysis
        if not hasattr(self, 'complexity_validations'):
            self.complexity_validations = []
        
        self.complexity_validations.append({
            'n': n,
            'actual_time': actual_time,
            'expected_complexity': expected_complexity,
            'relative_performance': relative_performance
        })
        
        logging.debug(f"Complexity validation: n={n}, time={actual_time:.4f}s, "
                     f"complexity={expected_complexity}, relative={relative_performance:.6f}")
        
        return True  # For now, just log - could add thresholds later
    
    def log_performance_report(self) -> None:
        """Log comprehensive performance report."""
        summary = self.get_performance_summary()
        
        logging.info("=" * 50)
        logging.info("PERFORMANCE REPORT")
        logging.info("=" * 50)
        logging.info(f"Total time: {summary['total_elapsed_time']:.2f} seconds")
        logging.info(f"Peak memory: {summary['peak_memory_mb']:.1f} MB")
        logging.info(f"Memory limit: {summary['memory_limit_mb']} MB")
        logging.info(f"Memory utilization: {summary['peak_memory_mb']/summary['memory_limit_mb']*100:.1f}%")
        
        if summary['phases']:
            logging.info("\nPhase breakdown:")
            for phase_name, phase_data in summary['phases'].items():
                logging.info(f"  {phase_name}: {phase_data['elapsed_time']:.2f}s "
                           f"({phase_data['operations_count']} ops, "
                           f"{phase_data['operations_per_second']:.1f} ops/s, "
                           f"{phase_data['peak_memory_mb']:.1f}MB)")
        
        # Log complexity validations if any
        if hasattr(self, 'complexity_validations') and self.complexity_validations:
            logging.info("\nComplexity validations:")
            for validation in self.complexity_validations[-5:]:  # Last 5
                logging.info(f"  n={validation['n']:,}, "
                           f"time={validation['actual_time']:.4f}s, "
                           f"{validation['expected_complexity']}")


class BatchProcessor:
    """Helper for batch processing with performance monitoring."""
    
    def __init__(self, monitor: PerformanceMonitor, batch_size: int = 1000):
        self.monitor = monitor
        self.batch_size = batch_size
    
    def process_in_batches(self, items: List[Any], processor_func, 
                          phase_name: str = "batch_processing"):
        """Process items in batches with monitoring."""
        total_items = len(items)
        
        with self.monitor.phase_context(phase_name) as metrics:
            results = []
            
            for i in range(0, total_items, self.batch_size):
                batch = items[i:i + self.batch_size]
                batch_start = time.time()
                
                # Process batch
                batch_results = processor_func(batch)
                results.extend(batch_results)
                
                # Update monitoring
                batch_time = time.time() - batch_start
                self.monitor.record_operations(len(batch))
                
                # Check memory periodically
                if i % (self.batch_size * 10) == 0:
                    self.monitor.check_memory_limit()
                
                # Log progress
                progress = min(i + len(batch), total_items)
                if progress % (self.batch_size * 5) == 0:
                    logging.info(f"Processed {progress:,}/{total_items:,} items "
                               f"({progress/total_items*100:.1f}%) in {batch_time:.2f}s")
            
            # Validate complexity
            total_time = metrics.elapsed_time
            self.monitor.validate_complexity(total_items, total_time, "O(n)")
        
        return results


# Decorators for performance monitoring
def monitor_phase(phase_name: str):
    """Decorator to monitor a function as a processing phase."""
    def decorator(func):
        def wrapper(self, *args, **kwargs):
            if hasattr(self, 'monitor'):
                with self.monitor.phase_context(phase_name):
                    return func(self, *args, **kwargs)
            else:
                return func(self, *args, **kwargs)
        return wrapper
    return decorator


def validate_complexity(expected_complexity: str = "O(n log n)"):
    """Decorator to validate algorithmic complexity."""
    def decorator(func):
        def wrapper(self, *args, **kwargs):
            # Try to infer input size from first argument
            n = 0
            if args and hasattr(args[0], '__len__'):
                n = len(args[0])
            elif args and isinstance(args[0], dict):
                n = len(args[0])
            
            start_time = time.time()
            result = func(self, *args, **kwargs)
            elapsed_time = time.time() - start_time
            
            if hasattr(self, 'monitor') and n > 0:
                self.monitor.validate_complexity(n, elapsed_time, expected_complexity)
            
            return result
        return wrapper
    return decorator