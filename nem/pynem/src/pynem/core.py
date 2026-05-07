"""Core NEM algorithm (E-step, M-step, convergence)."""

import numpy as np
import networkx as nx

from .models import (
    Family, Dispersion, Proportion, EPSILON,
    compute_log_density, estimate_parameters,
)
from .spatial import NeighborhoodSystem


class NEM:
    """Neighborhood EM for spatial clustering on graphs.

    Parameters
    ----------
    n_clusters : int
        Number of clusters K.
    beta : float
        Spatial regularization strength. 0 = standard EM.
    algorithm : str
        'nem' (mean field), 'ncem' (ICM hard), 'gem' (Gibbs sampling).
    family : str
        'normal', 'laplace', or 'bernoulli'.
    dispersion : str
        's__', 'sk_', 's_d', or 'skd'.
    proportion : str
        'p_' (equal) or 'pk' (free).
    beta_mode : str
        'fix', 'psgrad', 'heu_d', 'heu_l'.
    init : str
        'sort' or 'random'.
    n_init : int
        Number of random initializations (only for init='random').
    max_iter : int
        Maximum number of EM iterations.
    tol : float
        Convergence threshold.
    convergence : str
        'classification' or 'criterion'.
    missing : str
        'replace' or 'ignore'.
    random_state : int or None
        Random seed.
    verbose : int
        Verbosity level.
    """

    def __init__(self, n_clusters=2, beta=1.0, algorithm="nem",
                 family="normal", dispersion="s__", proportion="pk",
                 beta_mode="fix", init="sort", n_init=1,
                 max_iter=100, tol=1e-3, convergence="classification",
                 missing="replace", random_state=None, verbose=0):
        self.n_clusters = n_clusters
        self.beta = beta
        self.algorithm = algorithm
        self.family = Family(family)
        self.dispersion = Dispersion(dispersion)
        self.proportion = Proportion(proportion)
        self.beta_mode = beta_mode
        self.init = init
        self.n_init = n_init
        self.max_iter = max_iter
        self.tol = tol
        self.convergence = convergence
        self.missing = missing
        self.random_state = random_state
        self.verbose = verbose

    def fit(self, G_or_X, graph=None):
        """Fit the NEM model.

        Parameters
        ----------
        G_or_X : nx.Graph or np.ndarray
            If nx.Graph: nodes must have 'features' attribute.
            If np.ndarray of shape (N, D): feature matrix.
        graph : nx.Graph or None
            Required if G_or_X is an array — provides the graph structure.

        Returns
        -------
        self
        """
        rng = np.random.default_rng(self.random_state)

        if isinstance(G_or_X, nx.Graph):
            G = G_or_X
            n = G.number_of_nodes()
            d = len(G.nodes[0]["features"])
            X = np.array([G.nodes[i]["features"] for i in range(n)],
                         dtype=float)
        else:
            X = np.asarray(G_or_X, dtype=float)
            G = graph
            n, d = X.shape

        K = self.n_clusters
        ns = NeighborhoodSystem(G)

        if self.beta_mode in ("heu_d", "heu_l"):
            self._fit_heuristic(X, G, ns, K, rng)
            return self

        n_runs = self.n_init if self.init == "random" else 1
        best_crit = -np.inf
        best_result = None

        for run in range(n_runs):
            result = self._run_once(X, ns, K, rng)
            crit_val = result["criteria"]["U"]
            if crit_val > best_crit:
                best_crit = crit_val
                best_result = result

        self._store_result(best_result)
        return self

    def predict(self, G_or_X=None):
        """Return hard labels."""
        return self.labels_

    def _run_once(self, X, ns, K, rng):
        """Single NEM run from one initialization."""
        N, D = X.shape
        beta = self.beta

        # Initialize
        C = self._initialize(X, ns, K, rng)

        # Compute sample statistics for random init
        history = []
        old_C = None
        old_crit = None

        for iteration in range(self.max_iter):
            # M-step
            params = estimate_parameters(
                X, C, self.family, self.dispersion, self.proportion,
                miss_mode=self.missing,
                old_centers=params["centers"] if iteration > 0 else None,
                old_dispersions=params["dispersions"] if iteration > 0 else None,
            ) if iteration > 0 else self._first_m_step(X, C)

            # Beta estimation (pseudo-gradient)
            if self.beta_mode == "psgrad":
                beta = self._estimate_beta(C, ns, beta)

            # E-step — also returns cached log_pkfki to avoid recomputation
            old_C = C.copy()
            C, log_pkfki_cached = self._e_step(X, C, params, ns, beta, K, rng)

            # Compute criteria, reusing the log_pkfki already computed in _e_step
            criteria = self._compute_criteria(
                X, C, params, ns, beta, log_pkfki=log_pkfki_cached
            )
            history.append(criteria.copy())

            if self.verbose:
                print(f"  iter {iteration}: U={criteria['U']:.4f} "
                      f"D={criteria['D']:.4f} G={criteria['G']:.4f} "
                      f"beta={beta:.4f}")

            # Convergence test
            if old_crit is not None and self._has_converged(old_C, C,
                                                            old_crit, criteria):
                break
            old_crit = criteria

        return {
            "membership": C,
            "labels": np.argmax(C, axis=1) + 1,
            "centers": params["centers"],
            "dispersions": params["dispersions"],
            "proportions": params["proportions"],
            "beta": beta,
            "criteria": criteria,
            "n_iter": iteration + 1,
            "history": history,
        }

    def _first_m_step(self, X, C):
        """First M-step (no old parameters)."""
        return estimate_parameters(
            X, C, self.family, self.dispersion, self.proportion,
            miss_mode=self.missing,
        )

    def _initialize(self, X, ns, K, rng):
        """Initialize the classification matrix C (N, K)."""
        N, D = X.shape
        C = np.zeros((N, K))

        if self.init == "sort":
            # Sort by first variable, partition into K equal groups
            sorted_var = 0
            vals = X[:, sorted_var].copy()
            # Handle NaN by setting to max for sorting
            nan_mask = np.isnan(vals)
            vals[nan_mask] = np.nanmax(vals) + 1 if not nan_mask.all() else 0
            order = np.argsort(vals)
            for rank, idx in enumerate(order):
                ki = (rank * K) // N
                ki = min(ki, K - 1)
                C[idx, ki] = 1.0

        elif self.init == "random":
            # Pick K distinct data points as centers, then do one E-step
            observed = ~np.isnan(X)
            sample_disp = np.nanvar(X, axis=0)
            sample_disp = np.maximum(sample_disp, EPSILON)

            # Pick K distinct centers
            centers = np.zeros((K, D))
            chosen = []
            for k in range(K):
                for _ in range(100):
                    idx = rng.integers(0, N)
                    if idx not in chosen:
                        chosen.append(idx)
                        break
                for d in range(D):
                    if observed[idx, d]:
                        centers[k, d] = X[idx, d]
                    else:
                        lo = np.nanmin(X[:, d]) if observed[:, d].any() else 0
                        hi = np.nanmax(X[:, d]) if observed[:, d].any() else 1
                        centers[k, d] = rng.uniform(lo, hi)

            dispersions = np.tile(sample_disp / K, (K, 1))
            proportions = np.full(K, 1.0 / K)

            # One E-step to get initial classification
            log_pkfki = compute_log_density(
                X, centers, dispersions, proportions, self.family
            )
            # Without spatial term (beta=0 for init)
            C = self._normalize_membership(log_pkfki, np.zeros((N, K)))

        return C

    def _e_step(self, X, C, params, ns, beta, K, rng):
        """E-step: update classification.

        Returns
        -------
        new_C : (N, K) updated membership matrix
        log_pkfki : (N, K) log joint densities (cached for reuse in criteria)
        """
        N = X.shape[0]

        log_pkfki = compute_log_density(
            X, params["centers"], params["dispersions"],
            params["proportions"], self.family,
        )

        # Sparse matmul: W @ C — no Python loop (NeighborhoodSystem builds CSR matrix)
        contexts = ns.compute_all_contexts(C)  # (N, K)
        new_C = self._normalize_membership(log_pkfki, beta * contexts)

        # Algorithm variant
        if self.algorithm == "ncem":
            hard = np.zeros_like(new_C)
            hard[np.arange(N), np.argmax(new_C, axis=1)] = 1.0
            new_C = hard
        elif self.algorithm == "gem":
            hard = np.zeros_like(new_C)
            for i in range(N):
                k = rng.choice(K, p=new_C[i])
                hard[i, k] = 1.0
            new_C = hard

        return new_C, log_pkfki

    def _normalize_membership(self, log_pkfki, spatial_term):
        """Normalize log-probabilities to membership probabilities."""
        N, K = log_pkfki.shape
        log_num = log_pkfki + spatial_term
        # Log-sum-exp normalization
        max_log = np.max(log_num, axis=1, keepdims=True)
        finite = np.isfinite(max_log.ravel())
        C = np.full((N, K), 1.0 / K)
        if finite.any():
            num = np.exp(log_num[finite] - max_log[finite])
            total = num.sum(axis=1, keepdims=True)
            total = np.maximum(total, EPSILON)
            C[finite] = num / total
        return C

    def _estimate_beta(self, C, ns, beta, max_beta_iter=50, step=0.0):
        """Pseudo-gradient ascent on pseudo-likelihood for beta."""
        N, K = C.shape
        contexts = ns.compute_all_contexts(C)

        for _ in range(max_beta_iter):
            # Gradient and second derivative of pseudo-likelihood
            # Q(beta) = sum_i [beta * sum_k c_ik * con_ik - log(sum_k exp(beta*con_ik))]
            beta_con = beta * contexts  # (N, K)
            max_bc = np.max(beta_con, axis=1, keepdims=True)
            exp_bc = np.exp(beta_con - max_bc)
            Z = exp_bc.sum(axis=1, keepdims=True)  # (N, 1)

            # sum_k c_ik * con_ik
            c_dot_con = (C * contexts).sum(axis=1)  # (N,)

            # sum_k con_ik * exp(beta*con_ik) / Z_i
            mean_con = (contexts * exp_bc / Z).sum(axis=1)  # (N,)

            grad = (c_dot_con - mean_con).sum()

            # Second derivative
            mean_con2 = (contexts ** 2 * exp_bc / Z).sum(axis=1)
            d2 = (mean_con2 - mean_con ** 2).sum()

            if abs(grad) < self.tol * N:
                break

            if step <= 0 and abs(d2) > EPSILON:
                # Newton step with damping factor 4
                beta += grad / (4 * abs(d2))
            else:
                beta += grad * (step / N) if step > 0 else 0.01

            beta = np.clip(beta, -5.0, 5.0)

        return beta

    def _compute_criteria(self, X, C, params, ns, beta, log_pkfki=None):
        """Compute all criteria: U, D, G, L, M.

        Parameters
        ----------
        log_pkfki : (N, K) array or None
            Pre-computed log joint densities.  When provided (e.g. cached from
            ``_e_step``), the expensive ``compute_log_density`` call is skipped.
        """
        N, K = C.shape
        if log_pkfki is None:
            log_pkfki = compute_log_density(
                X,
                params["centers"],
                params["dispersions"],
                params["proportions"],
                self.family,
            )

        # D (Hathaway) — use nansum so that 0 * -inf terms are treated as 0
        # (standard entropy convention: lim_{p->0} p*log(p) = 0)
        log_C = np.log(np.maximum(C, EPSILON))
        D = np.nansum(C * (log_pkfki - log_C))

        # G (geographic cohesion)
        contexts = ns.compute_all_contexts(C)
        G = (C * contexts).sum()

        # U (NEM criterion)
        U = D + 0.5 * beta * G

        # L (mixture log-likelihood)
        max_log = np.max(log_pkfki, axis=1, keepdims=True)
        finite = np.isfinite(max_log.ravel())
        L = 0.0
        if finite.any():
            sum_exp = np.sum(np.exp(log_pkfki[finite] - max_log[finite]),
                             axis=1)
            L = (max_log[finite].ravel() + np.log(np.maximum(sum_exp, EPSILON))).sum()

        # Z and M (Markovian pseudo-likelihood)
        beta_con = beta * contexts
        max_bc = np.max(beta_con, axis=1)
        Z_vals = max_bc + np.log(np.sum(np.exp(beta_con - max_bc[:, None]),
                                        axis=1))
        Z = -Z_vals.sum()
        M = D + beta * G + Z

        return {"U": U, "D": D, "G": G, "L": L, "M": M}

    def _has_converged(self, C_old, C_new, crit_old, crit_new):
        """Check convergence."""
        if self.convergence == "classification":
            return np.max(np.abs(C_new - C_old)) < self.tol
        elif self.convergence == "criterion":
            if abs(crit_new["U"]) < EPSILON:
                return True
            return abs(crit_new["U"] - crit_old["U"]) / abs(crit_new["U"]) < self.tol
        return False

    def _fit_heuristic(self, X, G, ns, K, rng):
        """Heuristic beta estimation: increase beta until criterion drops."""
        # First run with beta=0 (pure EM)
        self.beta = 0.0
        result0 = self._run_once(X, ns, K, rng)
        D0 = result0["criteria"]["D"]
        L0 = result0["criteria"]["L"]

        best_result = result0
        best_U = result0["criteria"]["U"]

        beta_step = 0.1
        max_beta = 2.0

        for trial_beta in np.arange(beta_step, max_beta + beta_step, beta_step):
            self.beta = trial_beta
            result = self._run_once(X, ns, K, rng)

            # Check stopping criterion
            if self.beta_mode == "heu_d":
                if result["criteria"]["D"] < 0.8 * D0:
                    break
            elif self.beta_mode == "heu_l":
                if abs(result["criteria"]["L"]) < 0.02 * abs(L0):
                    break

            if result["criteria"]["U"] > best_U:
                best_U = result["criteria"]["U"]
                best_result = result

        self._store_result(best_result)

    def _store_result(self, result):
        """Store algorithm results as attributes."""
        self.labels_ = result["labels"]
        self.membership_ = result["membership"]
        self.centers_ = result["centers"]
        self.dispersions_ = result["dispersions"]
        self.proportions_ = result["proportions"]
        self.beta_ = result["beta"]
        self.criteria_ = result["criteria"]
        self.n_iter_ = result["n_iter"]
        self.history_ = result["history"]
