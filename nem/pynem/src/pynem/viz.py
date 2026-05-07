"""Visualization utilities for NEM results."""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx


def _is_image(G):
    """Check if graph represents an image grid."""
    return G.graph.get("type", "").upper() == "I"


def _get_image_dims(G):
    """Get (nl, nc) for image graphs."""
    return G.graph["nl"], G.graph["nc"]


def _get_features_matrix(G):
    """Extract (N, D) feature matrix from graph."""
    N = G.number_of_nodes()
    D = G.graph.get("d", None)
    if D is None:
        D = len(G.nodes[0]["features"])
    X = np.array([G.nodes[i]["features"] for i in range(N)])
    return X


def plot_graph_clusters(G, labels, pos=None, ax=None, **kwargs):
    """Draw graph with nodes colored by cluster label.

    Parameters
    ----------
    G : nx.Graph
    labels : array-like of shape (N,)
        Cluster labels (1-based).
    pos : dict or None
        Node positions. If None, uses spring layout.
    ax : matplotlib Axes or None
    **kwargs : passed to nx.draw_networkx_nodes.

    Returns
    -------
    ax : matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(8, 6))
    if pos is None:
        pos = nx.spring_layout(G, seed=42)

    labels = np.asarray(labels)
    K = int(labels.max())
    cmap = plt.colormaps.get_cmap("tab10").resampled(K)

    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.2)
    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_color=labels - 1,
        cmap=cmap,
        vmin=0, vmax=K - 1,
        node_size=kwargs.get("node_size", 60),
        alpha=kwargs.get("alpha", 0.9),
    )
    ax.set_title("Graph clusters")
    ax.axis("off")
    return ax


def plot_membership(G, membership, class_idx=0, pos=None, ax=None):
    """Heatmap of soft membership for one class.

    Parameters
    ----------
    G : nx.Graph
    membership : (N, K) array
    class_idx : int
        Which class to display.
    pos : dict or None
    ax : matplotlib Axes or None

    Returns
    -------
    ax : matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(8, 6))
    if pos is None:
        pos = nx.spring_layout(G, seed=42)

    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.1)
    nodes = nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_color=membership[:, class_idx],
        cmap=plt.cm.viridis,
        vmin=0, vmax=1,
        node_size=60,
    )
    plt.colorbar(nodes, ax=ax, label=f"P(class {class_idx + 1})")
    ax.set_title(f"Membership — class {class_idx + 1}")
    ax.axis("off")
    return ax


def plot_convergence(history, criterion="U", ax=None):
    """Plot criterion value vs iteration.

    Parameters
    ----------
    history : list of dict
        Each dict has keys 'U', 'D', 'G', 'L', 'M'.
    criterion : str
        Which criterion to plot.
    ax : matplotlib Axes or None

    Returns
    -------
    ax : matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(6, 4))

    values = [h[criterion] for h in history]
    ax.plot(range(1, len(values) + 1), values, "o-", markersize=3)
    ax.set_xlabel("Iteration")
    ax.set_ylabel(criterion)
    ax.set_title(f"Convergence — {criterion}")
    ax.grid(True, alpha=0.3)
    return ax


def plot_cluster_centers(centers, labels=None, ax=None):
    """Parallel coordinates plot of cluster centers.

    Parameters
    ----------
    centers : (K, D) array
    labels : list of str or None
        Variable names.
    ax : matplotlib Axes or None

    Returns
    -------
    ax : matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(8, 4))

    K, D = centers.shape
    x = np.arange(D)
    cmap = plt.colormaps.get_cmap("tab10").resampled(K)

    for k in range(K):
        ax.plot(x, centers[k], "o-", color=cmap(k), label=f"Class {k + 1}",
                markersize=6)

    if labels is not None:
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
    else:
        ax.set_xticks(x)
        ax.set_xticklabels([f"Var {d + 1}" for d in range(D)])

    ax.set_ylabel("Center value")
    ax.set_title("Cluster centers")
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax


def plot_labels(G, labels, pos=None, ax=None, title="Labels", **kwargs):
    """Display labels on graph or image.

    For image graphs (type='I'), shows a 2D image.
    For other graphs, draws nodes colored by label.

    Parameters
    ----------
    G : nx.Graph
    labels : array-like of shape (N,)
        Cluster labels (1-based).
    pos : dict or None
        Node positions (graph mode only).
    ax : matplotlib Axes or None
    title : str
    **kwargs : passed to nx.draw_networkx_nodes (graph mode).

    Returns
    -------
    ax : matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(6, 5))

    labels = np.asarray(labels)
    K = int(labels.max())
    cmap = plt.colormaps.get_cmap("tab10").resampled(K)

    if _is_image(G):
        nl, nc = _get_image_dims(G)
        img = (labels - 1).reshape(nl, nc)
        ax.imshow(img, cmap=cmap, vmin=0, vmax=K - 1, interpolation="nearest")
    else:
        if pos is None:
            pos = nx.spring_layout(G, seed=42)
        nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.2)
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            node_color=labels - 1,
            cmap=cmap,
            vmin=0, vmax=K - 1,
            node_size=kwargs.get("node_size", 60),
            alpha=kwargs.get("alpha", 0.9),
        )

    ax.set_title(title)
    ax.axis("off")
    return ax


def plot_feature(G, feature_values, pos=None, ax=None, title="Feature",
                 cmap="viridis", **kwargs):
    """Display a single feature (continuous) on graph or image.

    Parameters
    ----------
    G : nx.Graph
    feature_values : array-like of shape (N,)
        Continuous values for each node.
    pos : dict or None
        Node positions (graph mode only).
    ax : matplotlib Axes or None
    title : str
    cmap : str
        Colormap name.
    **kwargs : passed to nx.draw_networkx_nodes (graph mode).

    Returns
    -------
    ax : matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(6, 5))

    values = np.asarray(feature_values)

    if _is_image(G):
        nl, nc = _get_image_dims(G)
        img = values.reshape(nl, nc)
        im = ax.imshow(img, cmap=cmap, interpolation="nearest")
        plt.colorbar(im, ax=ax, shrink=0.8)
    else:
        if pos is None:
            pos = nx.spring_layout(G, seed=42)
        nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.1)
        nodes = nx.draw_networkx_nodes(
            G, pos, ax=ax,
            node_color=values,
            cmap=plt.colormaps.get_cmap(cmap),
            node_size=kwargs.get("node_size", 60),
            alpha=kwargs.get("alpha", 0.9),
        )
        plt.colorbar(nodes, ax=ax, shrink=0.8)

    ax.set_title(title)
    ax.axis("off")
    return ax


def plot_results(G, model, true_labels=None, pos=None, var_names=None):
    """Plot estimated labels and per-dimension feature views.

    Creates a figure with:
    - Row 1: estimated labels (+ true labels if provided)
    - Row 2: one panel per emission dimension

    For image graphs, uses imshow. For general graphs, uses node coloring.

    Parameters
    ----------
    G : nx.Graph
        Graph with 'features' attribute on nodes.
    model : NEM
        Fitted NEM model (must have labels_ attribute).
    true_labels : array-like of shape (N,) or None
        True labels (1-based) for comparison.
    pos : dict or None
        Node positions for graph layout.
    var_names : list of str or None
        Names for each feature dimension.

    Returns
    -------
    fig : matplotlib Figure
    """
    X = _get_features_matrix(G)
    N, D = X.shape

    if var_names is None:
        var_names = [f"Variable {d + 1}" for d in range(D)]

    # Compute layout once for graph mode
    if not _is_image(G) and pos is None:
        pos = nx.spring_layout(G, seed=42)

    # Determine grid layout
    has_true = true_labels is not None
    n_label_cols = 2 if has_true else 1
    n_cols = max(n_label_cols, D)
    n_rows = 2

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4.5 * n_cols, 4 * n_rows))
    if n_cols == 1:
        axes = axes.reshape(n_rows, 1)

    # Row 1: labels
    plot_labels(G, model.labels_, pos=pos, ax=axes[0, 0],
                title="Estimated labels")
    if has_true:
        plot_labels(G, true_labels, pos=pos, ax=axes[0, 1],
                    title="True labels")
    # Hide unused label panels
    for j in range(n_label_cols, n_cols):
        axes[0, j].axis("off")

    # Row 2: feature dimensions
    for d in range(D):
        plot_feature(G, X[:, d], pos=pos, ax=axes[1, d],
                     title=var_names[d])
    # Hide unused feature panels
    for j in range(D, n_cols):
        axes[1, j].axis("off")

    fig.tight_layout()
    return fig
