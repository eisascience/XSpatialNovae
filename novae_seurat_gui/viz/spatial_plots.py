"""Spatial scatter plots."""

import logging
from typing import Optional

import anndata
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

logger = logging.getLogger(__name__)

# Marker size clamps for pixel-based sizing
_MIN_MARKER_SIZE = 1
_MAX_MARKER_SIZE = 20


def _clamp_size(size: float) -> float:
    """Return size clamped to [_MIN_MARKER_SIZE, _MAX_MARKER_SIZE]."""
    return float(max(_MIN_MARKER_SIZE, min(_MAX_MARKER_SIZE, size)))


def plot_spatial_scatter(
    adata: anndata.AnnData,
    color_by: Optional[str] = None,
    spatial_key: str = "spatial",
    size: float = 4,
    opacity: float = 0.7,
    title: Optional[str] = None,
    color_map: Optional[str] = None,
    width: int = 800,
    height: int = 600,
) -> go.Figure:
    """
    Create a spatial scatter plot coloured by metadata or expression.

    Marker size is fixed in pixels (Plotly default) so points remain
    consistently readable regardless of zoom level.  Use the ``size``
    parameter to adjust the apparent size; it is clamped to
    [1, 20] pixels.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    color_by : str, optional
        Column in adata.obs to colour by.  If None, all points are the same
        colour.
    spatial_key : str
        Key in adata.obsm for spatial coordinates.
    size : float
        Marker diameter in pixels (clamped to [1, 20]).
    opacity : float
        Marker opacity (0-1).
    title : str, optional
        Plot title.
    color_map : str, optional
        Colormap name (e.g. ``"viridis"``, ``"Set1"``).
    width : int
        Figure width in pixels.
    height : int
        Figure height in pixels.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")

    coords = adata.obsm[spatial_key]
    x = coords[:, 0]
    y = coords[:, 1]

    plot_data = pd.DataFrame({"x": x, "y": y})

    if color_by and color_by in adata.obs.columns:
        plot_data[color_by] = adata.obs[color_by].values
        color_col = color_by
    else:
        plot_data["cell"] = "Cell"
        color_col = "cell"

    pixel_size = _clamp_size(size)

    if pd.api.types.is_numeric_dtype(plot_data[color_col]):
        fig = px.scatter(
            plot_data,
            x="x",
            y="y",
            color=color_col,
            color_continuous_scale=color_map or "viridis",
            opacity=opacity,
            title=title or (f"Spatial: {color_by}" if color_by else "Spatial"),
        )
    else:
        fig = px.scatter(
            plot_data,
            x="x",
            y="y",
            color=color_col,
            color_discrete_sequence=px.colors.qualitative.Set1 if not color_map else None,
            opacity=opacity,
            title=title or (f"Spatial: {color_by}" if color_by else "Spatial"),
        )

    fig.update_traces(marker=dict(size=pixel_size))
    fig.update_layout(
        width=width,
        height=height,
        xaxis_title="X",
        yaxis_title="Y",
        plot_bgcolor="white",
        xaxis=dict(showgrid=True, gridcolor="lightgray"),
        yaxis=dict(
            showgrid=True,
            gridcolor="lightgray",
            scaleanchor="x",
            scaleratio=1,
        ),
        dragmode="zoom",
    )

    return fig


def plot_qc_spatial(
    adata: anndata.AnnData,
    qc_mask: np.ndarray,
    spatial_key: str = "spatial",
    size: float = 4,
    opacity: float = 0.7,
    title: str = "QC Filtering: Kept vs. Filtered",
    width: int = 800,
    height: int = 600,
) -> go.Figure:
    """
    Plot spatial coordinates coloured by QC pass/fail status.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    qc_mask : np.ndarray
        Boolean mask: True = kept, False = filtered.
    spatial_key : str
        Key in adata.obsm for spatial coordinates.
    size : float
        Marker diameter in pixels (clamped to [1, 20]).
    opacity : float
        Marker opacity (0-1).
    title : str
        Plot title.
    width : int
        Figure width in pixels.
    height : int
        Figure height in pixels.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")

    coords = adata.obsm[spatial_key]
    x = coords[:, 0]
    y = coords[:, 1]

    qc_status = np.where(qc_mask, "Kept", "Filtered")
    plot_data = pd.DataFrame({"x": x, "y": y, "QC Status": qc_status})

    pixel_size = _clamp_size(size)

    fig = px.scatter(
        plot_data,
        x="x",
        y="y",
        color="QC Status",
        color_discrete_map={"Kept": "blue", "Filtered": "red"},
        opacity=opacity,
        title=title,
        category_orders={"QC Status": ["Kept", "Filtered"]},
    )

    fig.update_traces(marker=dict(size=pixel_size))
    fig.update_layout(
        width=width,
        height=height,
        xaxis_title="X",
        yaxis_title="Y",
        plot_bgcolor="white",
        xaxis=dict(showgrid=True, gridcolor="lightgray"),
        yaxis=dict(
            showgrid=True,
            gridcolor="lightgray",
            scaleanchor="x",
            scaleratio=1,
        ),
        dragmode="zoom",
    )

    return fig
