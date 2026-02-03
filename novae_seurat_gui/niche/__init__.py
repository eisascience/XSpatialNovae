"""Niche analysis module (optional secondary runs)."""

from .composition import compute_neighborhood_composition
from .clustering import cluster_niches, assign_niche_labels

__all__ = [
    "compute_neighborhood_composition",
    "cluster_niches",
    "assign_niche_labels",
]
