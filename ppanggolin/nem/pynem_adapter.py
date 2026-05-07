"""Adapter that wraps pynem.NEM with a ppanggolin-specific extension.

Extension
---------
``param_file`` init mode — replicates ppanggolin's ``.m`` file initialization:
pass ``init_centers``, ``init_dispersions``, ``init_proportions`` to the
constructor alongside ``init="param_file"``.

Notes
-----
* The ``nansum`` fix for ``0 * -inf = NaN`` in the D criterion is applied
  directly in ``pynem.core.NEM._compute_criteria`` (upstream fix).
* All other algorithm improvements (sparse spatial matrix, vectorised
  log-density) are in pynem core and need no adapter-level overrides.
"""

import numpy as np

from pynem import NEM
from pynem.models import compute_log_density, EPSILON


class PPanGGOLiNEM(NEM):
    """``pynem.NEM`` subclass for use inside PPanGGOLiN.

    Additional constructor parameters
    ----------------------------------
    init_centers : np.ndarray, shape (K, D), optional
        Initial Bernoulli centers (mu).  Required when ``init="param_file"``.
    init_dispersions : np.ndarray, shape (K, D), optional
        Initial dispersions (epsilon).  Required when ``init="param_file"``.
    init_proportions : np.ndarray, shape (K,), optional
        Initial class proportions.  Defaults to uniform when omitted.

    All other parameters are forwarded unchanged to ``pynem.NEM``.
    """

    def __init__(self, *args,
                 init_centers=None,
                 init_dispersions=None,
                 init_proportions=None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self._ppang_centers = (None if init_centers is None
                               else np.asarray(init_centers, dtype=float))
        self._ppang_dispersions = (None if init_dispersions is None
                                   else np.asarray(init_dispersions, dtype=float))
        self._ppang_proportions = (None if init_proportions is None
                                   else np.asarray(init_proportions, dtype=float))

    # ------------------------------------------------------------------
    # param_file initialisation (ppanggolin-specific init mode)
    # ------------------------------------------------------------------

    def _initialize(self, X, ns, K, rng):
        """Override: add ``init="param_file"`` mode on top of upstream modes."""
        if self.init != "param_file":
            return super()._initialize(X, ns, K, rng)

        if self._ppang_centers is None or self._ppang_dispersions is None:
            raise ValueError(
                "PPanGGOLiNEM: init='param_file' requires "
                "init_centers and init_dispersions to be provided."
            )
        N = X.shape[0]
        proportions = (self._ppang_proportions
                       if self._ppang_proportions is not None
                       else np.full(K, 1.0 / K))
        log_pkfki = compute_log_density(
            X,
            self._ppang_centers,
            self._ppang_dispersions,
            proportions,
            self.family,
        )
        return self._normalize_membership(log_pkfki, np.zeros((N, K)))
