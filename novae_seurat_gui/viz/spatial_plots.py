"""Spatial scatter plots."""

import logging
from typing import Optional, Union

import anndata
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

logger = logging.getLogger(__name__)


def plot_spatial_scatter(
    adata: anndata.AnnData,
    color_by: Optional[str] = None,
    spatial_key: str = "spatial",
    size: float = 3,
    opacity: float = 0.7,
    title: Optional[str] = None,
    color_map: Optional[str] = None,
    width: int = 800,
    height: int = 600,
    scale_with_zoom: bool = True,
) -> go.Figure:
    """
    Create spatial scatter plot colored by metadata or expression.
    
    When scale_with_zoom=True, markers are sized in data coordinates so they
    scale proportionally when zooming in/out.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    color_by : str, optional
        Column in adata.obs to color by. If None, uses default coloring.
    spatial_key : str
        Key in adata.obsm for spatial coordinates.
    size : float
        Marker size. When scale_with_zoom=True, this is interpreted as a
        fraction of the data range (default 3). When False, it's pixels.
    opacity : float
        Marker opacity.
    title : str, optional
        Plot title.
    color_map : str, optional
        Colormap name (e.g., 'viridis', 'Set1').
    width : int
        Figure width in pixels.
    height : int
        Figure height in pixels.
    scale_with_zoom : bool
        If True, markers scale in data coordinates (grow/shrink with zoom).
        If False, markers stay constant pixel size. Default is True.

    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object with zoom-responsive marker sizing.
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")

    coords = adata.obsm[spatial_key]
    x = coords[:, 0]
    y = coords[:, 1]

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
        # Continuous color scale
        fig = px.scatter(
            plot_data,
            x="x",
            y="y",
            color=color_col,
            color_continuous_scale=color_map or "viridis",
            opacity=opacity,
            title=title or f"Spatial: {color_by}" if color_by else "Spatial",
        )
    else:
        # Categorical color scale
        fig = px.scatter(
            plot_data,
            x="x",
            y="y",
            color=color_col,
            color_discrete_sequence=px.colors.qualitative.Set1 if not color_map else None,
            opacity=opacity,
            title=title or f"Spatial: {color_by}" if color_by else "Spatial",
        )

    # Configure marker sizing
    if scale_with_zoom:
        # Scale markers in data coordinates so they grow/shrink with zoom
        # Calculate appropriate size based on data range
        x_range = float(x.max() - x.min())
        y_range = float(y.max() - y.min())
        avg_range = (x_range + y_range) / 2
        
        # Marker size as fraction of data range
        # size parameter acts as a scaling factor
        marker_size_data = avg_range * (size / 300.0)  # 300 is empirical for good default
        
        # Use sizemode='area' with constant size value
        # This makes all markers the same size in DATA coordinates
        fig.update_traces(
            marker=dict(
                size=marker_size_data,
                sizemode='diameter',
                sizeref=1,  # Size is in data coordinates
            )
        )
    else:
        # Use pixel-based sizing (Plotly default)
        fig.update_traces(marker=dict(size=size))
    
    fig.update_layout(
        width=width,
        height=height,
        xaxis_title="X",
        yaxis_title="Y",
        plot_bgcolor="white",
        xaxis=dict(
            showgrid=True, 
            gridcolor="lightgray",
        ),
        yaxis=dict(
            showgrid=True, 
            gridcolor="lightgray", 
            scaleanchor="x", 
            scaleratio=1,
        ),
        # Enable drag mode for zoom/pan
        dragmode='zoom',
    )

    return fig


def plot_qc_spatial(
    adata: anndata.AnnData,
    qc_mask: np.ndarray,
    spatial_key: str = "spatial",
    size: float = 3,
    opacity: float = 0.7,
    title: str = "QC Filtering: Kept vs. Filtered",
    width: int = 800,
    height: int = 600,
    scale_with_zoom: bool = True,
) -> go.Figure:
    """
    Plot spatial coordinates colored by QC pass/fail.
    
    When scale_with_zoom=True, markers scale proportionally with zoom.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    qc_mask : np.ndarray
        Boolean mask indicating cells that pass QC.
    spatial_key : str
        Key in adata.obsm for spatial coordinates.
    size : float
        Marker size.
    opacity : float
        Marker opacity.
    title : str
        Plot title.
    width : int
        Figure width in pixels.
    height : int
        Figure height in pixels.
    scale_with_zoom : bool
        If True, markers scale in data coordinates (grow/shrink with zoom).
        Default is True.

    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object with zoom-responsive marker sizing.
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")

    coords = adata.obsm[spatial_key]
    x = coords[:, 0]
    y = coords[:, 1]

    qc_status = np.where(qc_mask, "Kept", "Filtered")

    plot_data = pd.DataFrame({"x": x, "y": y, "QC Status": qc_status})

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

    # Configure marker sizing
    if scale_with_zoom:
        # Scale markers in data coordinates
        x_range = float(x.max() - x.min())
        y_range = float(y.max() - y.min())
        avg_range = (x_range + y_range) / 2
        marker_size_data = avg_range * (size / 300.0)
        
        fig.update_traces(
            marker=dict(
                size=marker_size_data,
                sizemode='diameter',
                sizeref=1,
            )
        )
    else:
        fig.update_traces(marker=dict(size=size))
    
    fig.update_layout(
        width=width,
        height=height,
        xaxis_title="X",
        yaxis_title="Y",
        plot_bgcolor="white",
        xaxis=dict(showgrid=True, gridcolor="lightgray"),
        yaxis=dict(showgrid=True, gridcolor="lightgray", scaleanchor="x", scaleratio=1),
        dragmode='zoom',
    )

    return fig


def plot_spatial_with_edges(
    adata: anndata.AnnData,
    spatial_key: str = "spatial",
    connectivities_key: str = "spatial_connectivities",
    color_by: Optional[str] = None,
    max_edges: int = 5000,
    size: float = 3,
    opacity: float = 0.7,
    edge_width: float = 0.5,
    edge_color: str = "gray",
    width: int = 800,
    height: int = 600,
) -> go.Figure:
    """
    Plot spatial coordinates with neighbor edges overlaid.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    spatial_key : str
        Key in adata.obsm for spatial coordinates.
    connectivities_key : str
        Key in adata.obsp for connectivity matrix.
    color_by : str, optional
        Column in adata.obs to color points by.
    max_edges : int
        Maximum number of edges to plot (downsampled).
    size : float
        Marker size.
    opacity : float
        Marker opacity.
    edge_width : float
        Edge line width.
    edge_color : str
        Edge color.
    width : int
        Figure width.
    height : int
        Figure height.

    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object.
    """
    from ..spatial.diagnostics import get_neighbor_edges_for_visualization

    # Get edges
    edge_x, edge_y = get_neighbor_edges_for_visualization(
        adata, spatial_key, connectivities_key, max_edges
    )

    # Get coordinates
    coords = adata.obsm[spatial_key]
    x = coords[:, 0]
    y = coords[:, 1]

    # Create figure
    fig = go.Figure()

    # Add edges
    fig.add_trace(
        go.Scatter(
            x=edge_x,
            y=edge_y,
            mode="lines",
            line=dict(color=edge_color, width=edge_width),
            hoverinfo="skip",
            showlegend=False,
            name="Edges",
        )
    )

    # Add points
    if color_by and color_by in adata.obs.columns:
        colors = adata.obs[color_by].values
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="markers",
                marker=dict(size=size, color=colors, opacity=opacity),
                text=adata.obs[color_by].astype(str),
                hovertemplate="<b>%{text}</b><br>X: %{x}<br>Y: %{y}<extra></extra>",
                name="Cells",
            )
        )
    else:
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="markers",
                marker=dict(size=size, opacity=opacity),
                hovertemplate="X: %{x}<br>Y: %{y}<extra></extra>",
                name="Cells",
            )
        )

    fig.update_layout(
        title="Spatial Neighbor Graph",
        width=width,
        height=height,
        xaxis_title="X",
        yaxis_title="Y",
        plot_bgcolor="white",
        xaxis=dict(showgrid=True, gridcolor="lightgray"),
        yaxis=dict(showgrid=True, gridcolor="lightgray", scaleanchor="x", scaleratio=1),
    )

    return fig
