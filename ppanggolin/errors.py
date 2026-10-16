#!/usr/bin/env python3


class PangenomeError(Exception):
    """Base class for PPanGGOLiN-specific errors."""


class AlreadyComputedError(PangenomeError):
    """A step was asked for on a pangenome that already has its results, without --force."""


class PartitioningError(PangenomeError, RuntimeError):
    """NEM could not partition the pangenome."""
