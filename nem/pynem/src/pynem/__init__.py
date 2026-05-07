"""pynem — Neighborhood EM for spatial clustering on graphs."""

from .core import NEM
from . import io, viz, models, spatial, metrics

__version__ = "0.1.0"
__all__ = ["NEM", "io", "viz", "models", "spatial", "metrics"]
