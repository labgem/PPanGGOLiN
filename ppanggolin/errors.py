#!/usr/bin/env python3


class PangenomeError(Exception):
    """Base class of the errors PPanGGOLiN raises on purpose."""


class AlreadyComputedError(PangenomeError):
    """A step was asked for on a pangenome that already has its results, without --force."""


class PartitioningError(PangenomeError, RuntimeError):
    """NEM could not partition the pangenome."""


class WorkerDiedError(PangenomeError, RuntimeError):
    """A worker process died, so its task will never finish."""
