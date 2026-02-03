"""Embedding plots (PCA, UMAP, etc.)."""

import logging
from typing import Optional

import anndata
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

logger = logging.getLogger(__name__)


def plot_embedding(
    adata: anndata.AnnData,
    basis: str = "X_umap",
    color_by: Optional[str] = None,
    size: float = 3,
    opacity: float = 0.7,
    title: Optional[str] = None,
    color_map: Optional[str] = None,
    width: int = 700,
    height: int = 600,
) -> go.Figure:
    """
    Plot embedding (PCA, UMAP, Novae, etc.) colored by metadata.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    basis : str
        Key in adata.obsm for embedding coordinates.
    color_by : str, optional
        Column in adata.obs to color by.
    size : float
        Marker size.
    opacity : float
        Marker opacity.
    title : str, optional
        Plot title.
    color_map : str, optional
        Colormap name.
    width : int
        Figure width.
    height : int
        Figure height.

    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object.
    """
    if basis not in adata.obsm:
        raise ValueError(f"Basis '{basis}' not found in adata.obsm")

    embedding = adata.obsm[basis]

    # Use first 2 dimensions
    if embedding.shape[1] < 2:
        raise ValueError(f"Embedding '{basis}' has fewer than 2 dimensions")

    x = embedding[:, 0]
    y = embedding[:, 1]

    # Prepare data
    plot_data = pd.DataFrame({"x": x, "y": y})

    if color_by and color_by in adata.obs.columns:
        plot_data[color_by] = adata.obs[color_by].values
        color_col = color_by
    else:
        plot_data["cell"] = "Cell"
        color_col = "cell"

    # Create figure
    if pd.api.types.is_numeric_dtype(plot_data[color_col]):
        fig = px.scatter(
            plot_data,
            x="x",
            y="y",
            color=color_col,
            color_continuous_scale=color_map or "viridis",
            opacity=opacity,
            title=title or f"{basis}: {color_by}" if color_by else basis,
        )
    else:
        fig = px.scatter(
            plot_data,
            x="x",
            y="y",
            color=color_col,
            opacity=opacity,
            title=title or f"{basis}: {color_by}" if color_by else basis,
        )

    fig.update_traces(marker=dict(size=size))
    fig.update_layout(
        width=width,
        height=height,
        xaxis_title=f"{basis}_1",
        yaxis_title=f"{basis}_2",
        plot_bgcolor="white",
        xaxis=dict(showgrid=True, gridcolor="lightgray"),
        yaxis=dict(showgrid=True, gridcolor="lightgray"),
    )

    return fig


def plot_pca(
    adata: anndata.AnnData,
    color_by: Optional[str] = None,
    components: tuple = (0, 1),
    **kwargs,
) -> go.Figure:
    """
    Plot PCA embedding.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with PCA in obsm['X_pca'].
    color_by : str, optional
        Column to color by.
    components : tuple
        Which PC components to plot (0-indexed).
    **kwargs
        Additional arguments passed to plot_embedding.

    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object.
    """
    if "X_pca" not in adata.obsm:
        raise ValueError("PCA not found. Run preprocessing first.")

    # Extract specific components
    pca_subset = adata.obsm["X_pca"][:, list(components)]

    # Temporarily store in obsm
    adata.obsm["X_pca_plot"] = pca_subset

    fig = plot_embedding(
        adata,
        basis="X_pca_plot",
        color_by=color_by,
        title=f"PCA (PC{components[0]+1} vs PC{components[1]+1})",
        **kwargs,
    )

    # Clean up
    del adata.obsm["X_pca_plot"]

    return fig


def plot_umap(
    adata: anndata.AnnData, color_by: Optional[str] = None, **kwargs
) -> go.Figure:
    """
    Plot UMAP embedding.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with UMAP in obsm['X_umap'].
    color_by : str, optional
        Column to color by.
    **kwargs
        Additional arguments passed to plot_embedding.

    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object.
    """
    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP not found. Run UMAP computation first.")

    return plot_embedding(adata, basis="X_umap", color_by=color_by, title="UMAP", **kwargs)


def plot_variance_explained(
    adata: anndata.AnnData, n_comps: int = 50, width: int = 800, height: int = 400
) -> go.Figure:
    """
    Plot PCA variance explained.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with PCA variance in uns['pca']['variance_ratio'].
    n_comps : int
        Number of components to plot.
    width : int
        Figure width.
    height : int
        Figure height.

    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object.
    """
    if "pca" not in adata.uns or "variance_ratio" not in adata.uns["pca"]:
        raise ValueError("PCA variance not found. Run PCA first.")

    variance_ratio = adata.uns["pca"]["variance_ratio"][:n_comps]
    cumsum_variance = np.cumsum(variance_ratio)

    fig = go.Figure()

    # Bar plot of variance explained
    fig.add_trace(
        go.Bar(
            x=np.arange(1, len(variance_ratio) + 1),
            y=variance_ratio * 100,
            name="Variance Explained",
            marker_color="steelblue",
        )
    )

    # Line plot of cumulative variance
    fig.add_trace(
        go.Scatter(
            x=np.arange(1, len(cumsum_variance) + 1),
            y=cumsum_variance * 100,
            name="Cumulative Variance",
            mode="lines+markers",
            line=dict(color="red", width=2),
            marker=dict(size=5),
            yaxis="y2",
        )
    )

    fig.update_layout(
        title="PCA Variance Explained",
        xaxis_title="Principal Component",
        yaxis_title="Variance Explained (%)",
        yaxis2=dict(
            title="Cumulative Variance (%)",
            overlaying="y",
            side="right",
            range=[0, 100],
        ),
        width=width,
        height=height,
        plot_bgcolor="white",
        xaxis=dict(showgrid=True, gridcolor="lightgray"),
        yaxis=dict(showgrid=True, gridcolor="lightgray"),
    )

    return fig
