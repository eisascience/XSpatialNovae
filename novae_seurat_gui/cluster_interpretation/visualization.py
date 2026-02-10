"""Visualization functions for cluster interpretation."""

import logging
from typing import Optional

import anndata
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

logger = logging.getLogger(__name__)

# Constants for numerical stability
MIN_PVALUE_FOR_LOG = 1e-300  # Minimum p-value for log transformation to avoid log(0)


def plot_spatial_highlight(
    adata: anndata.AnnData,
    label_col: str,
    group_id: str,
    spatial_key: str = "spatial",
    size: float = 3,
    highlight_color: str = "#e74c3c",
    background_color: str = "#d3d3d3",
    title: Optional[str] = None,
    width: int = 800,
    height: int = 600,
) -> go.Figure:
    """
    Create spatial scatter plot with selected group highlighted.
    
    Background cells are shown in light grey, selected group in a strong color.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    label_col : str
        Column in adata.obs containing group labels.
    group_id : str
        Group ID to highlight.
    spatial_key : str
        Key in adata.obsm for spatial coordinates (default: "spatial").
    size : float
        Marker size.
    highlight_color : str
        Color for highlighted group (default: red).
    background_color : str
        Color for background cells (default: light grey).
    title : str, optional
        Plot title.
    width : int
        Figure width in pixels.
    height : int
        Figure height in pixels.
    
    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object.
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")
    
    if label_col not in adata.obs.columns:
        raise ValueError(f"Label column '{label_col}' not found in adata.obs")
    
    coords = adata.obsm[spatial_key]
    x = coords[:, 0]
    y = coords[:, 1]
    
    # Create highlight mask
    is_selected = (adata.obs[label_col] == group_id).values
    
    # Prepare plot data
    plot_data = pd.DataFrame({
        "x": x,
        "y": y,
        "group": ["Selected" if sel else "Other" for sel in is_selected],
    })
    
    # Create figure with two traces (background and highlighted)
    fig = go.Figure()
    
    # Add background cells first (so they're behind)
    background = plot_data[plot_data["group"] == "Other"]
    fig.add_trace(go.Scatter(
        x=background["x"],
        y=background["y"],
        mode="markers",
        marker=dict(
            size=size,
            color=background_color,
            opacity=0.3,
        ),
        name="Other",
        hovertemplate="<b>Other cells</b><br>x: %{x:.2f}<br>y: %{y:.2f}<extra></extra>",
    ))
    
    # Add highlighted cells
    highlighted = plot_data[plot_data["group"] == "Selected"]
    fig.add_trace(go.Scatter(
        x=highlighted["x"],
        y=highlighted["y"],
        mode="markers",
        marker=dict(
            size=size,
            color=highlight_color,
            opacity=0.8,
        ),
        name=f"{group_id}",
        hovertemplate=f"<b>{group_id}</b><br>x: %{{x:.2f}}<br>y: %{{y:.2f}}<extra></extra>",
    ))
    
    # Update layout
    fig.update_layout(
        title=title or f"Spatial: {group_id} vs Other",
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
        showlegend=True,
        legend=dict(
            x=1.02,
            y=1,
            xanchor="left",
            yanchor="top",
        ),
    )
    
    return fig


def plot_celltype_composition(
    composition_df: pd.DataFrame,
    title: Optional[str] = None,
    width: int = 600,
    height: int = 400,
) -> go.Figure:
    """
    Create stacked bar chart of cell type composition.
    
    Parameters
    ----------
    composition_df : pd.DataFrame
        DataFrame with columns: cell_type, n_cells, percent.
    title : str, optional
        Plot title.
    width : int
        Figure width in pixels.
    height : int
        Figure height in pixels.
    
    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object.
    """
    if len(composition_df) == 0:
        # Return empty figure
        fig = go.Figure()
        fig.update_layout(
            title=title or "Cell Type Composition",
            xaxis_title="Percentage",
            width=width,
            height=height,
            annotations=[
                dict(
                    text="No data available",
                    xref="paper",
                    yref="paper",
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                )
            ],
        )
        return fig
    
    # Create horizontal stacked bar
    fig = go.Figure()
    
    colors = px.colors.qualitative.Set2
    
    for i, row in composition_df.iterrows():
        fig.add_trace(go.Bar(
            x=[row["percent"]],
            y=["Composition"],
            orientation="h",
            name=row["cell_type"],
            marker=dict(color=colors[i % len(colors)]),
            hovertemplate=f"<b>{row['cell_type']}</b><br>"
                         f"{row['percent']:.1f}% ({row['n_cells']} cells)<extra></extra>",
        ))
    
    fig.update_layout(
        title=title or "Cell Type Composition",
        xaxis_title="Percentage (%)",
        barmode="stack",
        width=width,
        height=height,
        showlegend=True,
        legend=dict(
            orientation="v",
            x=1.02,
            y=1,
            xanchor="left",
            yanchor="top",
        ),
        yaxis=dict(showticklabels=False),
    )
    
    return fig


def create_qc_comparison_table(
    comparison_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Format QC comparison table for display.
    
    Parameters
    ----------
    comparison_df : pd.DataFrame
        DataFrame with QC metric comparisons.
    
    Returns
    -------
    pd.DataFrame
        Formatted DataFrame for display.
    """
    if len(comparison_df) == 0:
        return pd.DataFrame(columns=["Metric", "Mean (Group)", "Median (Group)", "Mean (Other)", "Median (Other)"])
    
    # Round numeric columns
    display_df = comparison_df.copy()
    
    numeric_cols = ["mean_in_group", "median_in_group", "mean_other", "median_other"]
    for col in numeric_cols:
        if col in display_df.columns:
            display_df[col] = display_df[col].round(2)
    
    # Rename columns for display
    display_df = display_df.rename(columns={
        "metric": "Metric",
        "mean_in_group": "Mean (Group)",
        "median_in_group": "Median (Group)",
        "mean_other": "Mean (Other)",
        "median_other": "Median (Other)",
    })
    
    return display_df


def plot_marker_genes_volcano(
    markers_df: pd.DataFrame,
    title: Optional[str] = None,
    width: int = 800,
    height: int = 600,
) -> go.Figure:
    """
    Create volcano plot for marker genes.
    
    Parameters
    ----------
    markers_df : pd.DataFrame
        DataFrame with marker gene statistics.
    title : str, optional
        Plot title.
    width : int
        Figure width in pixels.
    height : int
        Figure height in pixels.
    
    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure object.
    """
    if len(markers_df) == 0:
        # Return empty figure
        fig = go.Figure()
        fig.update_layout(
            title=title or "Volcano Plot",
            width=width,
            height=height,
            annotations=[
                dict(
                    text="No marker genes found",
                    xref="paper",
                    yref="paper",
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                )
            ],
        )
        return fig
    
    # Prepare data
    plot_data = markers_df.copy()
    plot_data["-log10(p)"] = -np.log10(plot_data["p_value"] + MIN_PVALUE_FOR_LOG)
    
    # Create scatter plot
    fig = px.scatter(
        plot_data,
        x="logFC",
        y="-log10(p)",
        hover_data=["gene", "adj_p_value"],
        title=title or "Volcano Plot: Marker Genes",
        labels={"logFC": "Log2 Fold Change", "-log10(p)": "-log10(p-value)"},
    )
    
    fig.update_traces(
        marker=dict(size=6, opacity=0.7, color="red"),
    )
    
    fig.update_layout(
        width=width,
        height=height,
        plot_bgcolor="white",
        xaxis=dict(showgrid=True, gridcolor="lightgray"),
        yaxis=dict(showgrid=True, gridcolor="lightgray"),
    )
    
    return fig
