"""
Novae-Seurat-GUI: Python-first workflow for Novae spatial foundation model on Seurat objects.

This package provides tools to:
- Load and validate Seurat objects exported as H5AD
- Apply quality control filters
- Run Novae spatial foundation model (zero-shot or train-from-scratch)
- Export results back to R/Seurat format
"""

__version__ = "0.1.0"
__author__ = "EISA Science"

from . import io, qc, modeling, viz, export, cluster_interpretation

__all__ = ["io", "qc", "modeling", "viz", "export", "cluster_interpretation", "__version__"]
