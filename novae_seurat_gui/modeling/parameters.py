"""Parameter management for Novae runs."""

from dataclasses import dataclass, field
from typing import Literal, Optional


@dataclass
class NovaeParameters:
    """Parameters for Novae model runs."""

    # Data mode
    mode: Literal["cosmx", "phenocycler"] = "cosmx"

    # Preprocessing
    target_sum: Optional[float] = 1e4
    n_top_genes: int = 2000
    n_comps: int = 50
    skip_normalize: bool = False
    skip_log: bool = False

    # Spatial neighbors
    neighbors_method: Literal["radius", "knn"] = "radius"
    neighbors_radius: Optional[float] = 150.0  # in coordinate units
    neighbors_k: Optional[int] = None

    # Model
    model_name: str = "MICS-Lab/novae-human-0"
    use_pretrained: bool = True

    # Domain assignment
    n_domains: int = 10
    hierarchical: bool = False
    n_levels: int = 1

    # Training (if from scratch)
    n_epochs: int = 100
    batch_size: int = 256
    learning_rate: float = 1e-3

    # Misc
    random_state: int = 42
    batch_key: Optional[str] = None  # For batch correction

    def to_dict(self):
        """Convert to dictionary."""
        return {
            "mode": self.mode,
            "target_sum": self.target_sum,
            "n_top_genes": self.n_top_genes,
            "n_comps": self.n_comps,
            "skip_normalize": self.skip_normalize,
            "skip_log": self.skip_log,
            "neighbors_method": self.neighbors_method,
            "neighbors_radius": self.neighbors_radius,
            "neighbors_k": self.neighbors_k,
            "model_name": self.model_name,
            "use_pretrained": self.use_pretrained,
            "n_domains": self.n_domains,
            "hierarchical": self.hierarchical,
            "n_levels": self.n_levels,
            "n_epochs": self.n_epochs,
            "batch_size": self.batch_size,
            "learning_rate": self.learning_rate,
            "random_state": self.random_state,
            "batch_key": self.batch_key,
        }

    @classmethod
    def from_dict(cls, d: dict):
        """Create from dictionary."""
        return cls(**{k: v for k, v in d.items() if k in cls.__annotations__})


def get_default_parameters(mode: Literal["cosmx", "phenocycler"] = "cosmx") -> NovaeParameters:
    """
    Get default parameters for a specific mode.

    Parameters
    ----------
    mode : {'cosmx', 'phenocycler'}
        Data modality.

    Returns
    -------
    NovaeParameters
        Default parameters for the specified mode.
    """
    if mode == "cosmx":
        return NovaeParameters(
            mode="cosmx",
            target_sum=1e4,
            n_top_genes=2000,
            skip_normalize=False,
            skip_log=False,
            neighbors_radius=150.0,
            model_name="MICS-Lab/novae-human-0",
            use_pretrained=True,
        )
    elif mode == "phenocycler":
        return NovaeParameters(
            mode="phenocycler",
            target_sum=None,
            n_top_genes=None,  # Use all proteins
            skip_normalize=True,  # Use quantile scaling instead
            skip_log=True,
            neighbors_radius=50.0,  # Smaller for proteomics
            model_name="phenocycler_custom",
            use_pretrained=False,  # Train from scratch
        )
    else:
        raise ValueError(f"Unknown mode: {mode}")


def validate_parameters(params: NovaeParameters) -> tuple[bool, list[str]]:
    """
    Validate parameters.

    Parameters
    ----------
    params : NovaeParameters
        Parameters to validate.

    Returns
    -------
    tuple of (bool, list)
        (is_valid, list of error messages)
    """
    errors = []

    # Check neighbors
    if params.neighbors_method == "radius" and params.neighbors_radius is None:
        errors.append("neighbors_radius must be specified when neighbors_method='radius'")
    if params.neighbors_method == "knn" and params.neighbors_k is None:
        errors.append("neighbors_k must be specified when neighbors_method='knn'")

    # Check domains
    if params.n_domains < 2:
        errors.append("n_domains must be >= 2")

    # Check hierarchical
    if params.hierarchical and params.n_levels < 2:
        errors.append("n_levels must be >= 2 when hierarchical=True")

    # Check training parameters
    if not params.use_pretrained:
        if params.n_epochs < 1:
            errors.append("n_epochs must be >= 1")
        if params.batch_size < 1:
            errors.append("batch_size must be >= 1")
        if params.learning_rate <= 0:
            errors.append("learning_rate must be > 0")

    is_valid = len(errors) == 0

    return is_valid, errors
