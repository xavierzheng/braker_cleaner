#!/usr/bin/env python3

"""
Unit tests for configuration management.

Tests the configuration loading, validation, and environment
variable handling functionality.
"""

import unittest
import tempfile
import os
import json
import sys

# Add the parent directory to the path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from gene_curation_pipeline.core.config import PipelineConfig, load_config
from gene_curation_pipeline.core.exceptions import ConfigurationError


class TestPipelineConfig(unittest.TestCase):
    """Test PipelineConfig class."""
    
    def test_default_config(self):
        """Test default configuration values."""
        config = PipelineConfig()
        
        # Test default values
        self.assertEqual(config.memory_limit_mb, 4096)
        self.assertEqual(config.batch_size, 1000)
        self.assertEqual(config.min_aa_length, 50)
        self.assertEqual(config.overlap_threshold, 0.5)
        self.assertTrue(config.enable_memory_monitoring)
        self.assertTrue(config.enable_codon_integration)
        self.assertEqual(config.parallel_workers, 1)
        self.assertFalse(config.debug_mode)
    
    def test_config_validation(self):
        """Test configuration validation."""
        # Test valid config
        config = PipelineConfig()
        config.validate()  # Should not raise
        
        # Test invalid values
        with self.assertRaises(ConfigurationError):
            PipelineConfig(min_aa_length=0)
        
        with self.assertRaises(ConfigurationError):
            PipelineConfig(overlap_threshold=-0.01)
        
        with self.assertRaises(ConfigurationError):
            PipelineConfig(overlap_threshold=1.5)
        
        with self.assertRaises(ConfigurationError):
            PipelineConfig(memory_limit_mb=50)
        
        with self.assertRaises(ConfigurationError):
            PipelineConfig(batch_size=0)
        
        with self.assertRaises(ConfigurationError):
            PipelineConfig(parallel_workers=0)
        
        with self.assertRaises(ConfigurationError):
            PipelineConfig(adjacency_gap_threshold=-1)
    
    def test_config_from_dict(self):
        """Test creating config from dictionary."""
        config_dict = {
            "memory_limit_mb": 8192,
            "min_aa_length": 20,
            "overlap_threshold": 0.7,
            "debug_mode": True,
            "unknown_key": "ignored"  # Should be filtered out
        }
        
        config = PipelineConfig.from_dict(config_dict)
        
        self.assertEqual(config.memory_limit_mb, 8192)
        self.assertEqual(config.min_aa_length, 20)
        self.assertEqual(config.overlap_threshold, 0.7)
        self.assertTrue(config.debug_mode)
        # Default values for unspecified parameters
        self.assertEqual(config.batch_size, 1000)
    
    def test_config_to_dict(self):
        """Test converting config to dictionary."""
        config = PipelineConfig(memory_limit_mb=8192, debug_mode=True)
        config_dict = config.to_dict()
        
        self.assertIsInstance(config_dict, dict)
        self.assertEqual(config_dict["memory_limit_mb"], 8192)
        self.assertTrue(config_dict["debug_mode"])
        self.assertIn("batch_size", config_dict)
    
    def test_config_from_json_file(self):
        """Test loading config from JSON file."""
        config_data = {
            "memory_limit_mb": 2048,
            "min_aa_length": 30,
            "parallel_workers": 4,
            "debug_mode": True
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(config_data, f)
            config_path = f.name
        
        try:
            config = PipelineConfig.from_file(config_path)
            
            self.assertEqual(config.memory_limit_mb, 2048)
            self.assertEqual(config.min_aa_length, 30)
            self.assertEqual(config.parallel_workers, 4)
            self.assertTrue(config.debug_mode)
            # Default for unspecified
            self.assertEqual(config.batch_size, 1000)
        finally:
            os.unlink(config_path)
    
    def test_config_from_nonexistent_file(self):
        """Test error handling for nonexistent config file."""
        with self.assertRaises(ConfigurationError):
            PipelineConfig.from_file("/nonexistent/config.json")
    
    def test_config_from_invalid_json(self):
        """Test error handling for invalid JSON."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            f.write("{ invalid json }")
            config_path = f.name
        
        try:
            with self.assertRaises(ConfigurationError):
                PipelineConfig.from_file(config_path)
        finally:
            os.unlink(config_path)
    
    def test_config_save_to_file(self):
        """Test saving config to file."""
        config = PipelineConfig(memory_limit_mb=2048, debug_mode=True)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            config_path = f.name
        
        try:
            config.save_to_file(config_path)
            
            # Load and verify
            loaded_config = PipelineConfig.from_file(config_path)
            self.assertEqual(loaded_config.memory_limit_mb, 2048)
            self.assertTrue(loaded_config.debug_mode)
        finally:
            if os.path.exists(config_path):
                os.unlink(config_path)
    
    def test_config_from_env(self):
        """Test loading config from environment variables."""
        # Set environment variables
        env_vars = {
            'PIPELINE_MEMORY_LIMIT_MB': '2048',
            'PIPELINE_MIN_AA_LENGTH': '25',
            'PIPELINE_OVERLAP_THRESHOLD': '0.6',
            'PIPELINE_DEBUG_MODE': 'true',
            'PIPELINE_PARALLEL_WORKERS': '2'
        }
        
        # Save original env
        original_env = {}
        for key in env_vars:
            original_env[key] = os.environ.get(key)
            os.environ[key] = env_vars[key]
        
        try:
            config = PipelineConfig.from_env()
            
            self.assertEqual(config.memory_limit_mb, 2048)
            self.assertEqual(config.min_aa_length, 25)
            self.assertEqual(config.overlap_threshold, 0.6)
            self.assertTrue(config.debug_mode)
            self.assertEqual(config.parallel_workers, 2)
            # Default for unspecified
            self.assertEqual(config.batch_size, 1000)
            
        finally:
            # Restore original env
            for key, value in original_env.items():
                if value is None:
                    if key in os.environ:
                        del os.environ[key]
                else:
                    os.environ[key] = value
    
    def test_config_from_env_invalid_values(self):
        """Test error handling for invalid environment values."""
        os.environ['PIPELINE_MEMORY_LIMIT_MB'] = 'invalid'
        
        try:
            with self.assertRaises(ConfigurationError):
                PipelineConfig.from_env()
        finally:
            if 'PIPELINE_MEMORY_LIMIT_MB' in os.environ:
                del os.environ['PIPELINE_MEMORY_LIMIT_MB']


class TestLoadConfig(unittest.TestCase):
    """Test the load_config function."""
    
    def test_load_default_config(self):
        """Test loading default configuration."""
        config = load_config()
        
        # Should have default values
        self.assertEqual(config.memory_limit_mb, 4096)
        self.assertEqual(config.min_aa_length, 50)
    
    def test_load_config_with_file(self):
        """Test loading config with file override."""
        config_data = {
            "memory_limit_mb": 1024,
            "min_aa_length": 15
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(config_data, f)
            config_path = f.name
        
        try:
            config = load_config(config_path=config_path)
            
            self.assertEqual(config.memory_limit_mb, 1024)
            self.assertEqual(config.min_aa_length, 15)
            # Default for unspecified
            self.assertEqual(config.batch_size, 1000)
        finally:
            os.unlink(config_path)
    
    def test_load_config_priority(self):
        """Test configuration loading priority: file > env > defaults."""
        # Set environment variable
        os.environ['PIPELINE_MEMORY_LIMIT_MB'] = '2048'
        
        # Create config file with different value
        config_data = {"memory_limit_mb": 1024}
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(config_data, f)
            config_path = f.name
        
        try:
            config = load_config(config_path=config_path, use_env=True)
            
            # File should override environment
            self.assertEqual(config.memory_limit_mb, 1024)
            
        finally:
            os.unlink(config_path)
            if 'PIPELINE_MEMORY_LIMIT_MB' in os.environ:
                del os.environ['PIPELINE_MEMORY_LIMIT_MB']
    
    def test_load_config_no_env(self):
        """Test loading config without environment variables."""
        os.environ['PIPELINE_MEMORY_LIMIT_MB'] = '2048'
        
        try:
            config = load_config(use_env=False)
            
            # Should use default, not environment
            self.assertEqual(config.memory_limit_mb, 4096)
            
        finally:
            if 'PIPELINE_MEMORY_LIMIT_MB' in os.environ:
                del os.environ['PIPELINE_MEMORY_LIMIT_MB']


class TestConfigurationIntegration(unittest.TestCase):
    """Test configuration integration scenarios."""
    
    def test_production_config_scenario(self):
        """Test a realistic production configuration."""
        config_data = {
            "memory_limit_mb": 8192,
            "batch_size": 2000,
            "min_aa_length": 20,
            "overlap_threshold": 0.6,
            "parallel_workers": 8,
            "enable_memory_monitoring": True,
            "enable_codon_integration": True,
            "debug_mode": False,
            "checkpoint_enabled": True
        }
        
        config = PipelineConfig.from_dict(config_data)
        
        # Should validate successfully
        config.validate()
        
        # Verify all values
        self.assertEqual(config.memory_limit_mb, 8192)
        self.assertEqual(config.batch_size, 2000)
        self.assertEqual(config.min_aa_length, 20)
        self.assertEqual(config.overlap_threshold, 0.6)
        self.assertEqual(config.parallel_workers, 8)
        self.assertTrue(config.enable_memory_monitoring)
        self.assertTrue(config.enable_codon_integration)
        self.assertFalse(config.debug_mode)
        self.assertTrue(config.checkpoint_enabled)
    
    def test_minimal_config_scenario(self):
        """Test minimal configuration with only essential parameters."""
        config_data = {
            "min_aa_length": 10,
            "memory_limit_mb": 1024
        }
        
        config = PipelineConfig.from_dict(config_data)
        config.validate()
        
        # Should have specified values and defaults for rest
        self.assertEqual(config.min_aa_length, 10)
        self.assertEqual(config.memory_limit_mb, 1024)
        self.assertEqual(config.batch_size, 1000)  # Default
        self.assertEqual(config.overlap_threshold, 0.5)  # Default
    
    def test_edge_case_values(self):
        """Test configuration with edge case values."""
        # Test minimum valid values
        config = PipelineConfig(
            memory_limit_mb=100,
            batch_size=1,
            min_aa_length=1,
            overlap_threshold=0.1,
            parallel_workers=1,
            adjacency_gap_threshold=0
        )
        
        config.validate()  # Should not raise
        
        # Test maximum reasonable values
        config = PipelineConfig(
            memory_limit_mb=32768,
            batch_size=10000,
            min_aa_length=1000,
            overlap_threshold=1.0,
            parallel_workers=32,
            adjacency_gap_threshold=10
        )
        
        config.validate()  # Should not raise


if __name__ == '__main__':
    unittest.main()
