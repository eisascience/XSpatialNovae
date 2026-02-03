#!/usr/bin/env python3
"""
Download and cache Novae models and assets for offline/cluster use.

This script prefetches Novae model weights from Hugging Face and caches them locally.
The cache location respects standard environment variables (HF_HOME, TRANSFORMERS_CACHE,
HUGGINGFACE_HUB_CACHE).

Usage:
    python scripts/download_models.py
    python scripts/download_models.py --novae-model-id MICS-Lab/novae-human-1
    python scripts/download_models.py --cache-dir /path/to/cache
    python scripts/download_models.py --optional-medsam

Environment Variables:
    HF_HOME: Base directory for Hugging Face cache
    TRANSFORMERS_CACHE: Cache directory for transformers models
    HUGGINGFACE_HUB_CACHE: Cache directory for Hugging Face Hub models
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Optional

try:
    from huggingface_hub import hf_hub_download, snapshot_download
    from huggingface_hub.utils import HfHubHTTPError, RepositoryNotFoundError
except ImportError:
    print("Error: huggingface_hub is not installed.")
    print("Install it with: pip install huggingface-hub")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def get_cache_dir(custom_cache_dir: Optional[str] = None) -> Path:
    """
    Get the cache directory for Hugging Face models.
    
    Priority order:
    1. Custom cache_dir argument
    2. HF_HOME environment variable
    3. TRANSFORMERS_CACHE environment variable
    4. HUGGINGFACE_HUB_CACHE environment variable
    5. Default: ~/.cache/huggingface/hub
    
    Parameters
    ----------
    custom_cache_dir : str, optional
        Custom cache directory path.
        
    Returns
    -------
    Path
        Cache directory path.
    """
    if custom_cache_dir:
        return Path(custom_cache_dir)
    
    # Check environment variables
    if "HF_HOME" in os.environ:
        return Path(os.environ["HF_HOME"]) / "hub"
    elif "TRANSFORMERS_CACHE" in os.environ:
        return Path(os.environ["TRANSFORMERS_CACHE"])
    elif "HUGGINGFACE_HUB_CACHE" in os.environ:
        return Path(os.environ["HUGGINGFACE_HUB_CACHE"])
    
    # Default location
    return Path.home() / ".cache" / "huggingface" / "hub"


def download_novae_model(
    model_id: str = "MICS-Lab/novae-human-0",
    cache_dir: Optional[str] = None,
) -> Path:
    """
    Download Novae model from Hugging Face Hub.
    
    Parameters
    ----------
    model_id : str
        Hugging Face model identifier (e.g., "MICS-Lab/novae-human-0").
    cache_dir : str, optional
        Custom cache directory. If None, uses default HF cache.
        
    Returns
    -------
    Path
        Path where the model was cached.
    """
    logger.info(f"Downloading Novae model: {model_id}")
    
    try:
        # Download entire model repository
        model_path = snapshot_download(
            repo_id=model_id,
            cache_dir=cache_dir,
            resume_download=True,
        )
        
        logger.info(f"✓ Model downloaded successfully to: {model_path}")
        return Path(model_path)
        
    except RepositoryNotFoundError:
        logger.error(f"✗ Model repository not found: {model_id}")
        logger.error(f"  Please verify the model ID is correct.")
        logger.error(f"  Available models may include:")
        logger.error(f"    - MICS-Lab/novae-human-0")
        logger.error(f"    - MICS-Lab/novae-human-1")
        raise
        
    except HfHubHTTPError as e:
        logger.error(f"✗ HTTP error downloading model: {e}")
        logger.error(f"  Check your internet connection and Hugging Face access.")
        raise
        
    except Exception as e:
        logger.error(f"✗ Unexpected error downloading model: {e}")
        raise


def download_optional_medsam(cache_dir: Optional[str] = None) -> Optional[Path]:
    """
    Download optional MedSAM model for segmentation workflows.
    
    NOTE: This is a placeholder. Update with actual MedSAM model ID if needed.
    
    Parameters
    ----------
    cache_dir : str, optional
        Custom cache directory.
        
    Returns
    -------
    Path or None
        Path where the model was cached, or None if not available.
    """
    logger.info("MedSAM download is not yet implemented.")
    logger.info("This is reserved for future segmentation features.")
    return None


def print_configuration_info(cache_dir: Path):
    """
    Print information about how to configure the application to use the cache.
    
    Parameters
    ----------
    cache_dir : Path
        Cache directory path.
    """
    logger.info("\n" + "="*70)
    logger.info("Configuration Information")
    logger.info("="*70)
    logger.info(f"\nModels cached at: {cache_dir}")
    logger.info("\nTo use this cache in your application:")
    logger.info("\n1. Set environment variable before running:")
    logger.info("   Unix/Linux/macOS:")
    logger.info(f"     export HF_HOME={cache_dir.parent.parent}")
    logger.info(f"     # OR")
    logger.info(f"     export HUGGINGFACE_HUB_CACHE={cache_dir}")
    logger.info("   Windows (Command Prompt):")
    logger.info(f"     set HF_HOME={cache_dir.parent.parent}")
    logger.info(f"     # OR")
    logger.info(f"     set HUGGINGFACE_HUB_CACHE={cache_dir}")
    logger.info("   Windows (PowerShell):")
    logger.info(f"     $env:HF_HOME='{cache_dir.parent.parent}'")
    logger.info(f"     # OR")
    logger.info(f"     $env:HUGGINGFACE_HUB_CACHE='{cache_dir}'")
    logger.info("\n2. In your code:")
    logger.info("   import os")
    logger.info(f"   os.environ['HF_HOME'] = '{cache_dir.parent.parent}'")
    logger.info("\n3. For Streamlit app (app.py):")
    logger.info("   Add at the top of the file:")
    logger.info(f"   os.environ['HF_HOME'] = '{cache_dir.parent.parent}'")
    logger.info("\n4. For offline/cluster use:")
    logger.info("   The models are now cached and will not be re-downloaded")
    logger.info("   unless the cache is cleared.")
    logger.info("="*70 + "\n")


def main():
    """Main entry point for model download script."""
    parser = argparse.ArgumentParser(
        description="Download and cache Novae models for offline use",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download default Novae model
  python scripts/download_models.py
  
  # Download specific model version
  python scripts/download_models.py --novae-model-id MICS-Lab/novae-human-1
  
  # Use custom cache directory
  python scripts/download_models.py --cache-dir /scratch/models
  
  # Download with optional MedSAM (future feature)
  python scripts/download_models.py --optional-medsam

Environment:
  HF_HOME: Base directory for Hugging Face cache
  TRANSFORMERS_CACHE: Cache directory for transformers
  HUGGINGFACE_HUB_CACHE: Cache directory for HF Hub
        """
    )
    
    parser.add_argument(
        "--novae-model-id",
        type=str,
        default="MICS-Lab/novae-human-0",
        help="Novae model identifier on Hugging Face (default: MICS-Lab/novae-human-0)"
    )
    
    parser.add_argument(
        "--cache-dir",
        type=str,
        default=None,
        help="Custom cache directory (default: uses HF_HOME or ~/.cache/huggingface/hub)"
    )
    
    parser.add_argument(
        "--optional-medsam",
        action="store_true",
        help="Also download optional MedSAM segmentation model (future feature)"
    )
    
    parser.add_argument(
        "--skip-novae",
        action="store_true",
        help="Skip downloading Novae model (useful for testing)"
    )
    
    args = parser.parse_args()
    
    # Get cache directory
    cache_dir = get_cache_dir(args.cache_dir)
    logger.info(f"Using cache directory: {cache_dir}")
    
    # Track success/failure
    success = True
    
    # Download Novae model
    if not args.skip_novae:
        try:
            model_path = download_novae_model(
                model_id=args.novae_model_id,
                cache_dir=str(cache_dir) if args.cache_dir else None
            )
        except Exception as e:
            logger.error(f"Failed to download Novae model: {e}")
            success = False
    else:
        logger.info("Skipping Novae model download (--skip-novae)")
    
    # Download optional MedSAM
    if args.optional_medsam:
        try:
            medsam_path = download_optional_medsam(
                cache_dir=str(cache_dir) if args.cache_dir else None
            )
            if medsam_path:
                logger.info(f"MedSAM model cached at: {medsam_path}")
        except Exception as e:
            logger.warning(f"Failed to download MedSAM model: {e}")
            # Don't fail if optional model fails
    
    # Print configuration info
    if success:
        print_configuration_info(cache_dir)
        logger.info("✓ All downloads completed successfully!")
        return 0
    else:
        logger.error("✗ Some downloads failed. See errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
