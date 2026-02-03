"""Streamlit GUI for Novae-Seurat workflow."""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

# Add package to path
sys.path.insert(0, str(Path(__file__).parent))

from novae_seurat_gui import io, qc, spatial, modeling, viz, export

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Page config
st.set_page_config(
    page_title="Novae-Seurat GUI",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Session state initialization
if "adata" not in st.session_state:
    st.session_state.adata = None
if "adata_filtered" not in st.session_state:
    st.session_state.adata_filtered = None
if "mappings" not in st.session_state:
    st.session_state.mappings = {}
if "qc_mask" not in st.session_state:
    st.session_state.qc_mask = None
if "preprocessing_done" not in st.session_state:
    st.session_state.preprocessing_done = False
if "novae_done" not in st.session_state:
    st.session_state.novae_done = False


# Title
st.title("🧬 Novae-Seurat GUI")
st.markdown(
    "Python-first workflow for Novae spatial foundation model on Seurat objects"
)

# Sidebar
with st.sidebar:
    st.header("Navigation")
    page = st.radio(
        "Select Page",
        [
            "📁 Load Data",
            "🔍 QC Filtering",
            "🕸️ Spatial Neighbors",
            "🧠 Run Novae",
            "📊 Results",
            "💾 Export",
        ],
    )

    st.divider()
    st.header("Dataset Info")
    if st.session_state.adata is not None:
        st.metric("Total Cells", st.session_state.adata.n_obs)
        st.metric("Features", st.session_state.adata.n_vars)
        if st.session_state.adata_filtered is not None:
            st.metric("Cells (after QC)", st.session_state.adata_filtered.n_obs)


# ==================== Load Data Page ====================
if page == "📁 Load Data":
    st.header("📁 Load Dataset")

    uploaded_file = st.file_uploader("Upload H5AD file", type=["h5ad"])

    if uploaded_file is not None:
        # Save temporarily
        temp_path = Path("/tmp") / uploaded_file.name
        with open(temp_path, "wb") as f:
            f.write(uploaded_file.getbuffer())

        if st.button("Load Dataset"):
            with st.spinner("Loading dataset..."):
                try:
                    adata = io.load_h5ad(str(temp_path))
                    st.session_state.adata = adata

                    # Detect mappings
                    mappings = io.detect_mappings(adata)
                    st.session_state.mappings = mappings

                    st.success(f"Loaded {adata.n_obs} cells × {adata.n_vars} features")

                except Exception as e:
                    st.error(f"Error loading file: {e}")

    if st.session_state.adata is not None:
        st.divider()
        st.subheader("Dataset Summary")

        adata = st.session_state.adata
        col1, col2, col3 = st.columns(3)
        col1.metric("Cells", adata.n_obs)
        col2.metric("Features", adata.n_vars)
        col3.metric("Obs Columns", len(adata.obs.columns))

        st.divider()
        st.subheader("Detected Mappings")

        mappings = st.session_state.mappings

        # Allow user to override
        col1, col2 = st.columns(2)

        with col1:
            st.write("**Spatial Coordinates**")
            x_col = st.selectbox(
                "X column",
                options=adata.obs.columns,
                index=(
                    list(adata.obs.columns).index(mappings["x_col"])
                    if mappings["x_col"]
                    else 0
                ),
            )
            y_col = st.selectbox(
                "Y column",
                options=adata.obs.columns,
                index=(
                    list(adata.obs.columns).index(mappings["y_col"])
                    if mappings["y_col"]
                    else 0
                ),
            )
            units = st.selectbox("Units", options=["mm", "px", "unknown"], index=0)

        with col2:
            st.write("**Metadata**")
            sample_id_col = st.selectbox(
                "Sample ID column",
                options=adata.obs.columns,
                index=(
                    list(adata.obs.columns).index(mappings["sample_id_col"])
                    if mappings["sample_id_col"]
                    else 0
                ),
            )
            cell_type_col = st.selectbox(
                "Cell type column (optional)",
                options=["None"] + list(adata.obs.columns),
                index=(
                    list(adata.obs.columns).index(mappings["cell_type_col"]) + 1
                    if mappings["cell_type_col"]
                    else 0
                ),
            )

        if st.button("Apply Mappings"):
            # Update mappings
            st.session_state.mappings.update(
                {
                    "x_col": x_col,
                    "y_col": y_col,
                    "units": units,
                    "sample_id_col": sample_id_col,
                    "cell_type_col": cell_type_col if cell_type_col != "None" else None,
                }
            )

            # Apply conversions
            adata = io.ensure_spatial_coords(adata, x_col=x_col, y_col=y_col)
            adata = io.normalize_metadata(adata, sample_id_col=sample_id_col)
            st.session_state.adata = adata

            st.success("Mappings applied!")


# ==================== QC Filtering Page ====================
elif page == "🔍 QC Filtering":
    st.header("🔍 Quality Control Filtering")

    if st.session_state.adata is None:
        st.warning("Please load a dataset first")
    else:
        adata = st.session_state.adata

        st.subheader("Filter Criteria")

        # Get numeric columns
        numeric_cols = adata.obs.select_dtypes(include=[np.number]).columns.tolist()

        # Filter selection
        filter_criteria = {}

        for col in numeric_cols[:6]:  # Show top 6 numeric columns
            col1, col2, col3 = st.columns([2, 1, 1])

            with col1:
                use_filter = st.checkbox(f"Filter by {col}", value=False)

            if use_filter:
                with col2:
                    min_val = st.number_input(
                        f"Min {col}",
                        value=float(adata.obs[col].min()),
                        key=f"min_{col}",
                    )
                with col3:
                    max_val = st.number_input(
                        f"Max {col}",
                        value=float(adata.obs[col].max()),
                        key=f"max_{col}",
                    )

                filter_criteria[col] = (min_val, max_val)

        if st.button("Apply Filters"):
            if filter_criteria:
                with st.spinner("Applying filters..."):
                    mask = qc.create_filter_mask(adata, filter_criteria)
                    st.session_state.qc_mask = mask

                    # Compute stats
                    stats = qc.compute_filter_stats(adata, mask)

                    st.success(
                        f"Filters applied: {stats['n_kept']} cells kept "
                        f"({stats['percent_kept']:.1f}%)"
                    )

                    # Create filtered copy
                    adata_filtered = adata[mask, :].copy()
                    st.session_state.adata_filtered = adata_filtered
            else:
                st.warning("No filters selected")

        # Visualize filtering
        if st.session_state.qc_mask is not None:
            st.divider()
            st.subheader("Filtering Results")

            mask = st.session_state.qc_mask

            # Stats
            col1, col2, col3 = st.columns(3)
            n_kept = np.sum(mask)
            n_filtered = np.sum(~mask)
            total = len(mask)

            col1.metric("Kept", n_kept)
            col2.metric("Filtered", n_filtered)
            col3.metric("% Kept", f"{100 * n_kept / total:.1f}%")

            # Spatial plot
            st.subheader("Spatial Distribution")
            if "spatial" in adata.obsm or (
                st.session_state.mappings.get("x_col")
                and st.session_state.mappings.get("y_col")
            ):
                fig = viz.plot_qc_spatial(adata, mask)
                st.plotly_chart(fig, use_container_width=True)


# ==================== Spatial Neighbors Page ====================
elif page == "🕸️ Spatial Neighbors":
    st.header("🕸️ Spatial Neighbor Graph")

    if st.session_state.adata_filtered is None:
        st.warning("Please apply QC filters first")
    else:
        adata = st.session_state.adata_filtered

        st.subheader("Neighbor Parameters")

        col1, col2 = st.columns(2)

        with col1:
            method = st.radio("Method", ["Radius", "K-Nearest Neighbors"])

        with col2:
            if method == "Radius":
                radius = st.number_input("Radius", value=150.0, min_value=1.0)
                n_neighbors = None
            else:
                n_neighbors = st.number_input("K", value=15, min_value=1)
                radius = None

        if st.button("Build Neighbor Graph"):
            with st.spinner("Building spatial neighbor graph..."):
                try:
                    if radius:
                        adata = spatial.compute_neighbors(
                            adata, method="radius", radius=radius
                        )
                    else:
                        adata = spatial.compute_neighbors(
                            adata, method="knn", n_neighbors=n_neighbors
                        )

                    st.session_state.adata_filtered = adata
                    st.success("Neighbor graph built!")

                except Exception as e:
                    st.error(f"Error: {e}")

        # Show diagnostics
        if "spatial_connectivities" in adata.obsp:
            st.divider()
            st.subheader("Graph Diagnostics")

            diag = spatial.graph_diagnostics(adata)

            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Edges", diag["n_edges"])
            col2.metric("Avg Degree", f"{diag['degree']['mean']:.1f}")
            col3.metric("Components", diag["connected_components"]["n_components"])
            col4.metric("Isolated Cells", diag["connected_components"]["isolated_cells"])

            # Degree distribution
            st.subheader("Degree Distribution")
            fig = spatial.plot_degree_distribution(adata)
            st.pyplot(fig)


# ==================== Run Novae Page ====================
elif page == "🧠 Run Novae":
    st.header("🧠 Run Novae Model")

    if st.session_state.adata_filtered is None:
        st.warning("Please complete QC filtering and spatial neighbors first")
    else:
        adata = st.session_state.adata_filtered

        st.subheader("Preprocessing")

        col1, col2 = st.columns(2)
        with col1:
            n_top_genes = st.number_input("Number of HVGs", value=2000, min_value=100)
            n_comps = st.number_input("PCA components", value=50, min_value=10)

        with col2:
            target_sum = st.number_input("Target sum (normalize)", value=10000.0)
            random_state = st.number_input("Random seed", value=42)

        if st.button("Run Preprocessing") and not st.session_state.preprocessing_done:
            with st.spinner("Preprocessing..."):
                try:
                    adata = modeling.preprocess_for_novae(
                        adata,
                        target_sum=target_sum,
                        n_top_genes=n_top_genes,
                        n_comps=n_comps,
                        random_state=random_state,
                    )
                    st.session_state.adata_filtered = adata
                    st.session_state.preprocessing_done = True
                    st.success("Preprocessing complete!")
                except Exception as e:
                    st.error(f"Error: {e}")

        if st.session_state.preprocessing_done:
            st.divider()
            st.subheader("Run Novae")

            col1, col2 = st.columns(2)
            with col1:
                model_name = st.text_input("Model", value="MICS-Lab/novae-human-0")
                n_domains = st.number_input("Number of domains", value=10, min_value=2)

            if st.button("Run Novae Model"):
                with st.spinner("Running Novae..."):
                    try:
                        adata = modeling.run_novae_zeroshot(
                            adata,
                            model_name=model_name,
                            n_domains=n_domains,
                            random_state=random_state,
                        )
                        st.session_state.adata_filtered = adata
                        st.session_state.novae_done = True
                        st.success(
                            f"Novae complete! Assigned {adata.obs['domain'].nunique()} domains"
                        )
                    except Exception as e:
                        st.error(f"Error: {e}")


# ==================== Results Page ====================
elif page == "📊 Results":
    st.header("📊 Results Visualization")

    if not st.session_state.novae_done:
        st.warning("Please run Novae first")
    else:
        adata = st.session_state.adata_filtered

        # Visualization options
        st.subheader("Spatial Plots")

        color_by = st.selectbox(
            "Color by",
            options=["domain"] + list(adata.obs.columns),
            index=0,
        )

        fig = viz.plot_spatial_scatter(adata, color_by=color_by, size=2)
        st.plotly_chart(fig, use_container_width=True)

        st.divider()
        st.subheader("Embeddings")

        # Check available embeddings
        if "X_novae" in adata.obsm:
            fig = viz.plot_embedding(adata, basis="X_novae", color_by=color_by)
            st.plotly_chart(fig, use_container_width=True)

        if "X_umap" in adata.obsm:
            fig = viz.plot_umap(adata, color_by=color_by)
            st.plotly_chart(fig, use_container_width=True)

        st.divider()
        st.subheader("Domain Summary")

        if "domain" in adata.obs:
            domain_counts = adata.obs["domain"].value_counts().sort_index()
            st.bar_chart(domain_counts)

            st.dataframe(
                domain_counts.reset_index().rename(
                    columns={"index": "Domain", "domain": "Cell Count"}
                )
            )


# ==================== Export Page ====================
elif page == "💾 Export":
    st.header("💾 Export Results")

    if not st.session_state.novae_done:
        st.warning("Please run Novae first")
    else:
        adata = st.session_state.adata_filtered

        st.subheader("Export Options")

        output_dir = st.text_input("Output Directory", value="./results")
        embedding_format = st.selectbox("Embedding Format", ["parquet", "csv"])

        if st.button("Export All"):
            with st.spinner("Exporting..."):
                try:
                    # Export
                    exported = export.export_all(
                        adata,
                        output_dir=output_dir,
                        embedding_format=embedding_format,
                    )

                    # Create manifest
                    manifest = export.create_manifest(
                        adata,
                        input_files=[],
                        parameters=adata.uns.get("novae_params", {}),
                    )

                    # Add spatial diagnostics
                    if "spatial_connectivities" in adata.obsp:
                        diag = spatial.graph_diagnostics(adata)
                        manifest = export.add_spatial_summary_to_manifest(manifest, diag)

                    # Save manifest
                    manifest_file = Path(output_dir) / "run_manifest.json"
                    export.save_manifest(manifest, str(manifest_file))

                    st.success("Export complete!")

                    st.subheader("Exported Files")
                    for key, path in exported.items():
                        st.text(f"{key}: {path}")
                    st.text(f"manifest: {manifest_file}")

                except Exception as e:
                    st.error(f"Error: {e}")

        st.divider()
        st.subheader("R Import Instructions")

        st.code(
            """
# In R:
library(Seurat)
library(arrow)

# Load domains
domains <- read.csv("results/domains.csv", row.names = "cell_id")
seurat_obj <- AddMetaData(seurat_obj, domains)

# Load embeddings
embeddings_df <- arrow::read_parquet("results/embeddings.parquet")
embeddings_mat <- as.matrix(embeddings_df[, -1])
rownames(embeddings_mat) <- embeddings_df$cell_id
embeddings_mat <- embeddings_mat[colnames(seurat_obj), ]

novae_dr <- CreateDimReducObject(
  embeddings = embeddings_mat,
  key = "NOVAe_",
  assay = DefaultAssay(seurat_obj)
)
seurat_obj[["novae"]] <- novae_dr

# Visualize
DimPlot(seurat_obj, reduction = "novae", group.by = "domain_level_0")
        """,
            language="r",
        )


# Footer
st.divider()
st.markdown(
    """
    <div style='text-align: center; color: gray;'>
    Novae-Seurat-GUI v0.1.0 | Built with Streamlit
    </div>
    """,
    unsafe_allow_html=True,
)
