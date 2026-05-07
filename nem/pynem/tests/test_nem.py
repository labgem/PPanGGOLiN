"""Tests for pynem, using the essai example data."""

import sys
from pathlib import Path
import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for tests
import matplotlib.pyplot as plt

# Add src to path for editable install
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pynem
from pynem import NEM
from pynem.models import Family, Dispersion, Proportion
from pynem.metrics import adjusted_rand_index

EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples"
BASENAME = str(EXAMPLES_DIR / "essai")


class TestIO:
    def test_read_str(self):
        t, n, d = pynem.io.read_str(BASENAME + ".str")
        assert t == "S"
        assert n == 100
        assert d == 3

    def test_read_str_image(self):
        basename = str(EXAMPLES_DIR / "potts_30x30_3")
        t, n, d = pynem.io.read_str(basename + ".str")
        assert t == "I"
        assert n == 900  # 30 * 30
        assert d == 2

    def test_read_dat(self):
        _, n, d = pynem.io.read_str(BASENAME + ".str")
        X = pynem.io.read_dat(BASENAME + ".dat", n, d)
        assert X.shape == (100, 3)
        assert not np.isnan(X).any()

    def test_read_nei(self):
        _, n, _ = pynem.io.read_str(BASENAME + ".str")
        G = pynem.io.read_nei(BASENAME + ".nei", n)
        assert G.number_of_nodes() == 100
        assert G.number_of_edges() > 0

    def test_read_graph(self):
        G = pynem.io.read_graph(BASENAME)
        assert G.number_of_nodes() == 100
        assert G.graph["d"] == 3
        # Check features are set
        feat = G.nodes[0]["features"]
        assert len(feat) == 3
        assert isinstance(feat, np.ndarray)
        # Check connectivity
        assert G.number_of_edges() > 50


class TestPureEM:
    """NEM with beta=0 reduces to standard EM."""

    def test_pure_em_converges(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.0, family="normal",
                    init="sort", max_iter=100, verbose=0)
        model.fit(G)
        assert model.labels_ is not None
        assert len(model.labels_) == 100
        assert model.n_iter_ > 1
        assert model.n_iter_ <= 100

    def test_pure_em_criteria(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.0, family="normal", init="sort")
        model.fit(G)
        # With beta=0, U = D and G should still be computed
        assert model.criteria_["U"] == pytest.approx(model.criteria_["D"])
        assert np.isfinite(model.criteria_["L"])


class TestNEMSpatial:
    def test_nem_with_beta(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.1, family="normal", init="sort")
        model.fit(G)
        assert len(model.labels_) == 100
        assert set(model.labels_).issubset({1, 2, 3})
        # U should be different from D when beta > 0
        assert model.criteria_["U"] != pytest.approx(model.criteria_["D"])
        assert model.criteria_["G"] > 0


class TestNCEM:
    def test_ncem(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.1, algorithm="ncem",
                    family="normal", init="sort")
        model.fit(G)
        # NCEM produces hard labels; membership should be 0/1
        assert np.all((model.membership_ == 0) | (model.membership_ == 1))
        assert np.all(model.membership_.sum(axis=1) == 1)


class TestFuzzyOutput:
    def test_membership_sums_to_one(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.1, family="normal", init="sort")
        model.fit(G)
        row_sums = model.membership_.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_membership_shape(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.0, family="normal", init="sort")
        model.fit(G)
        assert model.membership_.shape == (100, 3)


class TestSortInit:
    def test_sort_init_balanced(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.0, family="normal", init="sort",
                    max_iter=1)
        model.fit(G)
        # After sort init + 1 iteration, should have all 3 classes
        unique_labels = set(model.labels_)
        assert len(unique_labels) == 3


class TestModelVariants:
    @pytest.mark.parametrize("disp", ["s__", "sk_", "s_d", "skd"])
    def test_dispersion_models(self, disp):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.0, family="normal",
                    dispersion=disp, init="sort", max_iter=20)
        model.fit(G)
        assert model.dispersions_.shape == (3, 3)
        assert np.all(model.dispersions_ > 0)

    @pytest.mark.parametrize("prop", ["p_", "pk"])
    def test_proportion_models(self, prop):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.0, family="normal",
                    proportion=prop, init="sort", max_iter=20)
        model.fit(G)
        np.testing.assert_allclose(model.proportions_.sum(), 1.0, atol=1e-10)
        if prop == "p_":
            np.testing.assert_allclose(model.proportions_,
                                       np.full(3, 1.0 / 3), atol=1e-10)

    def test_laplace_family(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.0, family="laplace",
                    init="sort", max_iter=20)
        model.fit(G)
        assert len(model.labels_) == 100


class TestHeuristicBeta:
    def test_heuristic_d(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta_mode="heu_d", family="normal",
                    init="sort", max_iter=50)
        model.fit(G)
        assert model.beta_ >= 0
        assert len(model.labels_) == 100


class TestRandomInit:
    def test_random_init(self):
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=0.0, family="normal",
                    init="random", n_init=3, random_state=42, max_iter=30)
        model.fit(G)
        assert len(model.labels_) == 100
        assert set(model.labels_).issubset({1, 2, 3})


# ── Metrics tests ────────────────────────────────────────────────────────

class TestAdjustedRandIndex:
    def test_perfect_agreement(self):
        labels = [1, 1, 2, 2, 3, 3]
        assert adjusted_rand_index(labels, labels) == pytest.approx(1.0)

    def test_permuted_labels(self):
        """ARI should be 1.0 even when labels are permuted."""
        true = [1, 1, 2, 2, 3, 3]
        pred = [3, 3, 1, 1, 2, 2]
        assert adjusted_rand_index(true, pred) == pytest.approx(1.0)

    def test_random_partition(self):
        """ARI should be near 0 for unrelated partitions."""
        true = [1, 1, 1, 1, 1, 1]
        pred = [1, 1, 2, 2, 3, 3]
        assert adjusted_rand_index(true, pred) == pytest.approx(0.0)

    def test_two_clusters(self):
        true = [0, 0, 0, 1, 1, 1]
        pred = [0, 0, 1, 1, 1, 1]  # one mismatch
        ari = adjusted_rand_index(true, pred)
        assert -1.0 <= ari < 1.0
        assert ari > 0  # still partially correct

    def test_numpy_arrays(self):
        true = np.array([1, 1, 2, 2])
        pred = np.array([2, 2, 1, 1])
        assert adjusted_rand_index(true, pred) == pytest.approx(1.0)

    def test_with_nem_output(self):
        """ARI works with fitted NEM model labels."""
        G = pynem.io.read_graph(BASENAME)
        model = NEM(n_clusters=3, beta=1.0, family="normal", init="sort")
        model.fit(G)
        # Compare model to itself (should be 1.0)
        ari = adjusted_rand_index(model.labels_, model.labels_)
        assert ari == pytest.approx(1.0)


# ── Visualization tests ─────────────────────────────────────────────────

class TestPlotLabels:
    """Tests for plot_labels on graph data."""

    def setup_method(self):
        self.G = pynem.io.read_graph(BASENAME)
        self.model = NEM(n_clusters=3, beta=1.0, family="normal", init="sort")
        self.model.fit(self.G)

    def test_plot_labels_returns_axes(self):
        ax = pynem.viz.plot_labels(self.G, self.model.labels_)
        assert ax is not None
        plt.close("all")

    def test_plot_labels_with_title(self):
        ax = pynem.viz.plot_labels(self.G, self.model.labels_,
                                   title="Test title")
        assert ax.get_title() == "Test title"
        plt.close("all")

    def test_plot_labels_on_provided_axes(self):
        fig, ax = plt.subplots()
        returned_ax = pynem.viz.plot_labels(self.G, self.model.labels_, ax=ax)
        assert returned_ax is ax
        plt.close("all")


class TestPlotFeature:
    """Tests for plot_feature on graph data."""

    def setup_method(self):
        self.G = pynem.io.read_graph(BASENAME)

    def test_plot_feature_returns_axes(self):
        values = np.array([self.G.nodes[i]["features"][0]
                           for i in range(100)])
        ax = pynem.viz.plot_feature(self.G, values)
        assert ax is not None
        plt.close("all")

    def test_plot_feature_with_custom_cmap(self):
        values = np.array([self.G.nodes[i]["features"][0]
                           for i in range(100)])
        ax = pynem.viz.plot_feature(self.G, values, cmap="coolwarm")
        assert ax is not None
        plt.close("all")


class TestPlotResults:
    """Tests for plot_results composite figure."""

    def setup_method(self):
        self.G = pynem.io.read_graph(BASENAME)
        self.model = NEM(n_clusters=3, beta=1.0, family="normal", init="sort")
        self.model.fit(self.G)

    def test_plot_results_without_true_labels(self):
        fig = pynem.viz.plot_results(self.G, self.model)
        assert fig is not None
        axes = fig.get_axes()
        # Should have: 1 label panel + 3 feature panels (D=3)
        # Plus possibly hidden axes
        assert len(axes) >= 4
        plt.close("all")

    def test_plot_results_with_true_labels(self):
        true_labels = self.model.labels_  # use same as "true" for testing
        fig = pynem.viz.plot_results(self.G, self.model,
                                     true_labels=true_labels)
        assert fig is not None
        plt.close("all")

    def test_plot_results_with_var_names(self):
        fig = pynem.viz.plot_results(self.G, self.model,
                                     var_names=["A", "B", "C"])
        assert fig is not None
        plt.close("all")

    def test_plot_results_2d_features(self):
        """Test with 2D feature data (SBM 100 2D if available)."""
        sbm_basename = str(EXAMPLES_DIR / "sbm_100_2")
        try:
            G = pynem.io.read_graph(sbm_basename)
        except FileNotFoundError:
            pytest.skip("sbm_100_2 example data not available")
        model = NEM(n_clusters=3, beta=1.0, family="normal", init="sort")
        model.fit(G)
        fig = pynem.viz.plot_results(G, model)
        assert fig is not None
        plt.close("all")


class TestPlotLabelsImage:
    """Tests for plot_labels and plot_results on image data."""

    def setup_method(self):
        basename = str(EXAMPLES_DIR / "potts_30x30_3")
        try:
            self.G = pynem.io.read_graph(basename)
        except FileNotFoundError:
            pytest.skip("potts_30x30_3 example data not available")
        self.model = NEM(n_clusters=3, beta=1.0, family="normal", init="sort")
        self.model.fit(self.G)

    def test_image_graph_has_dims(self):
        assert self.G.graph["type"] == "I"
        assert "nl" in self.G.graph
        assert "nc" in self.G.graph
        assert self.G.graph["nl"] == 30
        assert self.G.graph["nc"] == 30

    def test_plot_labels_image(self):
        ax = pynem.viz.plot_labels(self.G, self.model.labels_)
        assert ax is not None
        plt.close("all")

    def test_plot_feature_image(self):
        values = np.array([self.G.nodes[i]["features"][0]
                           for i in range(900)])
        ax = pynem.viz.plot_feature(self.G, values)
        assert ax is not None
        plt.close("all")

    def test_plot_results_image(self):
        fig = pynem.viz.plot_results(self.G, self.model)
        assert fig is not None
        plt.close("all")


class TestExistingViz:
    """Ensure existing viz functions still work."""

    def setup_method(self):
        self.G = pynem.io.read_graph(BASENAME)
        self.model = NEM(n_clusters=3, beta=1.0, family="normal", init="sort")
        self.model.fit(self.G)

    def test_plot_graph_clusters(self):
        ax = pynem.viz.plot_graph_clusters(self.G, self.model.labels_)
        assert ax is not None
        plt.close("all")

    def test_plot_membership(self):
        ax = pynem.viz.plot_membership(self.G, self.model.membership_)
        assert ax is not None
        plt.close("all")

    def test_plot_convergence(self):
        ax = pynem.viz.plot_convergence(self.model.history_)
        assert ax is not None
        plt.close("all")

    def test_plot_cluster_centers(self):
        ax = pynem.viz.plot_cluster_centers(self.model.centers_)
        assert ax is not None
        plt.close("all")
