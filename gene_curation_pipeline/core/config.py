#!/usr/bin/env python3

"""
Configuration management for the gene curation pipeline.

Centralized configuration with support for file-based configuration
and environment variable overrides.
"""

import os
import json
from dataclasses import dataclass, asdict
from typing import Optional, Dict, Any
from .exceptions import ConfigurationError

# Optional YAML support
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


@dataclass
class PipelineConfig:
    """Centralized configuration for the gene curation pipeline."""
    
    # Performance settings
    memory_limit_mb: int = 4096
    batch_size: int = 1000
    enable_memory_monitoring: bool = True
    
    # Biological parameters
    min_aa_length: int = 50
    overlap_threshold: float = 0.5
    
    # Quality control
    max_internal_stops: int = 0
    require_start_codon: bool = True
    require_stop_codon: bool = True
    
    # Processing options
    enable_codon_integration: bool = True
    enable_sequence_reconstruction: bool = True
    adjacency_gap_threshold: int = 1  # bp
    
    # Output settings
    preserve_source_names: bool = True
    include_quality_flags: bool = True
    generate_reports: bool = True
    
    # Advanced settings
    parallel_workers: int = 1
    checkpoint_enabled: bool = False
    debug_mode: bool = False
    
    @classmethod
    def from_file(cls, config_path: str) -> 'PipelineConfig':
        """Load configuration from file (JSON or YAML)."""
        if not os.path.exists(config_path):
            raise ConfigurationError(f"Configuration file not found: {config_path}")
        
        try:
            with open(config_path, 'r') as f:
                if config_path.lower().endswith('.yaml') or config_path.lower().endswith('.yml'):
                    if not YAML_AVAILABLE:
                        raise ConfigurationError("PyYAML not installed but YAML config provided")
                    config_data = yaml.safe_load(f)
                else:
                    config_data = json.load(f)
            
            return cls.from_dict(config_data)
            
        except (json.JSONDecodeError, Exception) as e:
            if YAML_AVAILABLE:
                from yaml import YAMLError
                if isinstance(e, YAMLError):
                    raise ConfigurationError(f"Invalid YAML configuration file format: {e}")
            raise ConfigurationError(f"Invalid configuration file format: {e}")
        except Exception as e:
            raise ConfigurationError(f"Error loading configuration: {e}")
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'PipelineConfig':
        """Create configuration from dictionary."""
        # Filter out unknown keys
        known_keys = set(cls.__dataclass_fields__.keys())
        filtered_dict = {k: v for k, v in config_dict.items() if k in known_keys}
        
        try:
            return cls(**filtered_dict)
        except TypeError as e:
            raise ConfigurationError(f"Invalid configuration parameters: {e}")
    
    @classmethod
    def from_env(cls) -> 'PipelineConfig':
        """Load configuration from environment variables."""
        config = cls()
        
        # Map environment variables to config fields
        env_mappings = {
            'PIPELINE_MEMORY_LIMIT_MB': ('memory_limit_mb', int),
            'PIPELINE_BATCH_SIZE': ('batch_size', int),
            'PIPELINE_MIN_AA_LENGTH': ('min_aa_length', int),
            'PIPELINE_OVERLAP_THRESHOLD': ('overlap_threshold', float),
            'PIPELINE_PARALLEL_WORKERS': ('parallel_workers', int),
            'PIPELINE_DEBUG_MODE': ('debug_mode', lambda x: x.lower() in ('true', '1', 'yes')),
        }
        
        for env_var, (field_name, converter) in env_mappings.items():
            env_value = os.getenv(env_var)
            if env_value:
                try:
                    setattr(config, field_name, converter(env_value))
                except (ValueError, TypeError) as e:
                    raise ConfigurationError(f"Invalid environment variable {env_var}: {e}")
        
        return config
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return asdict(self)
    
    def save_to_file(self, config_path: str) -> None:
        """Save configuration to file."""
        config_dict = self.to_dict()
        
        try:
            with open(config_path, 'w') as f:
                if config_path.lower().endswith('.yaml') or config_path.lower().endswith('.yml'):
                    if not YAML_AVAILABLE:
                        raise ConfigurationError("PyYAML not installed but YAML output requested")
                    yaml.safe_dump(config_dict, f, default_flow_style=False)
                else:
                    json.dump(config_dict, f, indent=2)
        except Exception as e:
            raise ConfigurationError(f"Error saving configuration: {e}")
    
    def validate(self) -> None:
        """Validate configuration parameters."""
        if self.min_aa_length < 1:
            raise ConfigurationError("min_aa_length must be >= 1")
        
        if not 0 <= self.overlap_threshold <= 1:
            raise ConfigurationError("overlap_threshold must be between 0 and 1 (inclusive)")
        
        if self.memory_limit_mb < 100:
            raise ConfigurationError("memory_limit_mb must be >= 100")
        
        if self.batch_size < 1:
            raise ConfigurationError("batch_size must be >= 1")
        
        if self.parallel_workers < 1:
            raise ConfigurationError("parallel_workers must be >= 1")
        
        if self.adjacency_gap_threshold < 0:
            raise ConfigurationError("adjacency_gap_threshold must be >= 0")
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        self.validate()


def load_config(config_path: Optional[str] = None, 
                use_env: bool = True) -> PipelineConfig:
    """
    Load configuration with priority: file > environment > defaults.
    
    Args:
        config_path: Path to configuration file (optional)
        use_env: Whether to load environment variables
    
    Returns:
        PipelineConfig: Loaded configuration
    """
    # Start with defaults
    config = PipelineConfig()
    
    # Override with environment variables if requested
    if use_env:
        try:
            env_config = PipelineConfig.from_env()
            # Merge non-default values from environment
            for field_name, field_def in PipelineConfig.__dataclass_fields__.items():
                env_value = getattr(env_config, field_name)
                default_value = getattr(config, field_name)
                if env_value != default_value:
                    setattr(config, field_name, env_value)
        except ConfigurationError:
            # Environment config is optional
            pass
    
    # Override with file configuration if provided
    if config_path:
        file_config = PipelineConfig.from_file(config_path)
        # Merge all values from file
        for field_name in PipelineConfig.__dataclass_fields__.keys():
            setattr(config, field_name, getattr(file_config, field_name))
    
    return config
