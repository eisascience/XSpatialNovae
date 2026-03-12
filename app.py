"""Streamlit GUI for Novae-Seurat workflow."""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

# Add package to path
sys.path.insert(0, str(Path(__file__).parent))

from novae_seurat_gui import io, qc, modeling, viz, export, cluster_interpretation

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Page config
st.set_page_config(
    page_title="Novae-Seurat GUI",
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
st.title("Novae-Seurat GUI")
st.markdown(
    "Python-first workflow for Novae spatial foundation model on Seurat objects"
)

# Sidebar
with st.sidebar:
    st.header("Navigation")
    page = st.radio(
        "Select Page",
        [
            "Load Data",
            "QC Filtering",
            "Run Novae",
            "Results",
            "Domains & Markers",
            "Export",
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
if page == "Load Data":
    st.header("Load Dataset")

    st.markdown("""
    Upload your spatial transcriptomics data in one of the following formats:
    - **H5AD** (AnnData): Direct upload, ready for analysis
    - **RDS** (Seurat): Will be converted to H5AD automatically
    """)

    uploaded_file = st.file_uploader("Upload file", type=["h5ad", "rds"])

    if uploaded_file is not None:
        file_type = uploaded_file.name.split(".")[-1].lower()
        
        # Save temporarily
        temp_path = Path("/tmp") / uploaded_file.name
        with open(temp_path, "wb") as f:
            f.write(uploaded_file.getbuffer())

        st.info(f"File uploaded: {uploaded_file.name} ({file_type.upper()} format)")

        # Handle RDS conversion
        if file_type == "rds":
            st.divider()
            st.subheader("Seurat Conversion Settings")
            
            # Check R availability first
            r_available, r_msg = io.check_r_available()
            
            if not r_available:
                st.error(r_msg)
                st.markdown("""
                **To use .rds files, you need to install:**
                1. R (version >= 4.2)
                2. Required R packages:
                   ```r
                   install.packages(c("Seurat", "hdf5r", "optparse"))
                   remotes::install_github("mojaveazure/seurat-disk")
                   ```
                """)
                uploaded_file = None  # Reset to prevent further processing
            else:
                st.success(r_msg)
                
                # Check R packages
                packages_ok, packages_msg, missing = io.check_r_packages()
                if not packages_ok:
                    st.error("Missing R packages")
                    st.code(packages_msg)
                    uploaded_file = None
                else:
                    st.success("All required R packages installed")
                    
                    # Show conversion options
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        assay = st.text_input("Assay name", value="RNA", 
                                             help="Seurat assay to extract")
                        counts_slot = st.selectbox("Counts slot", 
                                                  options=["counts", "data"],
                                                  help="Which slot to use for count data")
                    
                    with col2:
                        x_col = st.text_input("X coordinate column", value="auto",
                                             help="Column name for x coordinates (auto to detect)")
                        y_col = st.text_input("Y coordinate column", value="auto",
                                             help="Column name for y coordinates (auto to detect)")
                    
                    sample_id_col = st.text_input("Sample ID column", value="auto",
                                                 help="Column name for sample ID (auto to detect)")
                    
                    with st.expander("Advanced Options"):
                        cell_id_col = st.text_input("Cell ID column (optional)", value="",
                                                   help="Leave empty to use rownames")
                        keep_meta_regex = st.text_input("Metadata filter regex (optional)", value="",
                                                       help="Regex to filter metadata columns")
                    
                    # Store conversion params in session state
                    if "rds_conversion_params" not in st.session_state:
                        st.session_state.rds_conversion_params = {}
                    
                    st.session_state.rds_conversion_params = {
                        "assay": assay,
                        "counts_slot": counts_slot,
                        "x_col": x_col,
                        "y_col": y_col,
                        "sample_id_col": sample_id_col,
                        "cell_id_col": cell_id_col if cell_id_col else None,
                        "keep_meta_regex": keep_meta_regex if keep_meta_regex else None,
                    }

        # Load button
        if uploaded_file is not None:
            if st.button("Load Dataset", type="primary"):
                with st.spinner("Loading dataset..."):
                    try:
                        if file_type == "rds":
                            # Convert RDS to H5AD
                            st.info("Converting Seurat .rds to H5AD format...")
                            
                            params = st.session_state.rds_conversion_params
                            output_dir = Path("/tmp/conversions")
                            
                            output_path, adata, info = io.convert_rds_to_h5ad_with_validation(
                                input_rds=str(temp_path),
                                output_dir=str(output_dir),
                                **params,
                                use_cache=True,
                                overwrite=False,
                            )
                            
                            # Show conversion summary
                            st.success("✓ Conversion complete!")
                            
                            with st.expander("Conversion Summary", expanded=True):
                                col1, col2, col3 = st.columns(3)
                                col1.metric("Cells", info["n_cells"])
                                col2.metric("Features", info["n_features"])
                                col3.metric("Assay", info["params"]["assay"])
                                
                                st.write("**Schema validation:**")
                                st.write(f"- Spatial coordinates: {'yes' if info['has_spatial'] else 'no'}")
                                st.write(f"- Cell IDs: {'yes' if info['has_cell_id'] else 'no'}")
                                st.write(f"- Sample IDs: {'yes' if info['has_sample_id'] else 'no'}")
                                
                                if info["cached"]:
                                    st.info("Used cached conversion from previous upload")
                                
                                # Show R script output if available
                                if not info["cached"] and "r_stdout" in info:
                                    with st.expander("Conversion Log", expanded=False):
                                        st.code(info["r_stdout"], language="text")
                            
                            st.session_state.adata = adata
                            st.session_state.conversion_info = info
                            
                        else:
                            # Load H5AD directly
                            adata = io.load_h5ad(str(temp_path))
                            st.session_state.adata = adata
                            st.success(f"Loaded {adata.n_obs} cells × {adata.n_vars} features")
                        
                        # Detect mappings
                        mappings = io.detect_mappings(st.session_state.adata)
                        st.session_state.mappings = mappings

                    except Exception as e:
                        st.error(f"Error loading file: {e}")
                        logger.exception("Error loading file")

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
elif page == "QC Filtering":
    st.header("Quality Control Filtering")

    if st.session_state.adata is None:
        st.warning("Please load a dataset first")
    else:
        adata = st.session_state.adata

        st.subheader("Filter Criteria")

        # Get numeric columns
        numeric_cols = adata.obs.select_dtypes(include=[np.number]).columns.tolist()

        # Filter selection - reactive (no apply button needed)
        filter_criteria = {}

        for col in numeric_cols[:6]:  # Show top 6 numeric columns
            col1, col2, col3 = st.columns([2, 1, 1])

            with col1:
                use_filter = st.checkbox(f"Filter by {col}", value=False, key=f"use_filter_{col}")

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

        # ------ Spatial coordinate filters ------
        st.divider()
        st.subheader("Spatial Coordinate Filters")
        st.caption(
            "Optionally restrict cells to an X/Y bounding box to remove "
            "surrounding whitespace or tissue edges."
        )

        # Detect spatial coordinates
        _has_spatial = "spatial" in adata.obsm
        _has_xy_obs = "x" in adata.obs.columns and "y" in adata.obs.columns

        if _has_spatial or _has_xy_obs:
            if _has_spatial:
                _sx = adata.obsm["spatial"][:, 0]
                _sy = adata.obsm["spatial"][:, 1]
            else:
                _sx = adata.obs["x"].values.astype(float)
                _sy = adata.obs["y"].values.astype(float)

            _x_data_min = float(np.nanmin(_sx))
            _x_data_max = float(np.nanmax(_sx))
            _y_data_min = float(np.nanmin(_sy))
            _y_data_max = float(np.nanmax(_sy))

            use_coord_filter = st.checkbox(
                "Apply spatial bounding-box filter", value=False, key="use_coord_filter"
            )

            _coord_x_min = None
            _coord_x_max = None
            _coord_y_min = None
            _coord_y_max = None

            if use_coord_filter:
                c1, c2, c3, c4 = st.columns(4)
                with c1:
                    _coord_x_min = st.number_input(
                        "X min",
                        value=_x_data_min,
                        key="coord_x_min",
                    )
                with c2:
                    _coord_x_max = st.number_input(
                        "X max",
                        value=_x_data_max,
                        key="coord_x_max",
                    )
                with c3:
                    _coord_y_min = st.number_input(
                        "Y min",
                        value=_y_data_min,
                        key="coord_y_min",
                    )
                with c4:
                    _coord_y_max = st.number_input(
                        "Y max",
                        value=_y_data_max,
                        key="coord_y_max",
                    )
        else:
            st.info(
                "No spatial coordinates found (adata.obsm['spatial'] or obs columns 'x'/'y'). "
                "Spatial bounding-box filter is disabled."
            )
            use_coord_filter = False
            _coord_x_min = _coord_x_max = _coord_y_min = _coord_y_max = None

        # Add reset button
        if st.button("Reset Filters"):
            keys_to_clear = []
            for col in numeric_cols[:6]:
                keys_to_clear += [f"use_filter_{col}", f"min_{col}", f"max_{col}"]
            keys_to_clear += [
                "use_coord_filter", "coord_x_min", "coord_x_max",
                "coord_y_min", "coord_y_max",
            ]
            for key in keys_to_clear:
                if key in st.session_state:
                    del st.session_state[key]
            st.rerun()

        # Build metadata mask reactively
        mask = qc.build_filter_mask(adata, filter_criteria)

        # Apply spatial coordinate mask if requested
        if use_coord_filter and (_has_spatial or _has_xy_obs):
            try:
                coord_mask = qc.create_coordinate_mask(
                    adata,
                    x_min=_coord_x_min,
                    x_max=_coord_x_max,
                    y_min=_coord_y_min,
                    y_max=_coord_y_max,
                )
                mask = mask & coord_mask
            except Exception as _e:
                st.error(f"Error applying coordinate filter: {_e}")

        # Store in session state
        st.session_state.qc_mask = mask

        # Create filtered adata
        if np.all(mask):
            adata_current = adata
        else:
            adata_current = adata[mask, :].copy()

        st.session_state.adata_filtered = adata_current

        # Always show filtering results and spatial plot
        st.divider()
        st.subheader("Filtering Results")

        # Compute stats
        stats = qc.compute_filter_stats(adata, mask)

        # Stats panel - always visible
        col1, col2, col3 = st.columns(3)
        col1.metric("Kept", stats['n_kept'])
        col2.metric("Filtered", stats['n_filtered'])
        col3.metric("% Kept", f"{stats['percent_kept']:.1f}%")

        if stats['n_filtered'] == 0:
            st.info("No filters applied; showing all cells.")
        else:
            st.success(f"Filters applied: {stats['n_kept']} cells kept ({stats['percent_kept']:.1f}%)")

        # Spatial plot - always show
        st.subheader("Spatial Distribution")

        # Point size slider for spatial plots
        _pt_size = st.slider("Point size", min_value=1, max_value=20, value=4, key="qc_point_size")

        if "spatial" in adata.obsm:
            fig = viz.plot_qc_spatial(adata, mask, size=_pt_size)
            st.plotly_chart(fig, use_container_width=True)
        elif (
            st.session_state.mappings.get("x_col")
            and st.session_state.mappings.get("y_col")
        ):
            st.info(
                "Spatial coordinates found in obs columns but not in obsm['spatial']. "
                "Apply mappings on the Load Data page to enable the spatial plot."
            )


# ==================== Run Novae Page ====================
elif page == "Run Novae":
    st.header("Run Novae Model")

    if st.session_state.adata is None:
        st.warning("Please load a dataset first")
    else:
        # Use filtered data if available, otherwise use original
        if st.session_state.adata_filtered is not None:
            adata = st.session_state.adata_filtered
        else:
            adata = st.session_state.adata
            st.info("Using full dataset (no QC filtering applied)")

        # Show mode information (novae availability)
        from novae_seurat_gui.modeling.novae_runner import _NOVAE_AVAILABLE
        if _NOVAE_AVAILABLE:
            st.success("Novae package detected. Running in full Novae mode.")
        else:
            st.info(
                "The optional 'novae' package is not installed. "
                "Running in proxy mode: UMAP + Leiden clustering will be used "
                "as a substitute for the Novae foundation model. "
                "Install 'novae' for the real model."
            )

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
                    # Update the appropriate session state
                    if st.session_state.adata_filtered is not None:
                        st.session_state.adata_filtered = adata
                    else:
                        st.session_state.adata = adata
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
                        # Update the appropriate session state
                        if st.session_state.adata_filtered is not None:
                            st.session_state.adata_filtered = adata
                        else:
                            st.session_state.adata = adata
                        st.session_state.novae_done = True
                        st.success(
                            f"Novae complete! Assigned {adata.obs['domain'].nunique()} domains"
                        )
                    except Exception as e:
                        # Check if it's a dependency error
                        from novae_seurat_gui.utils.deps import MissingDependency
                        if isinstance(e, MissingDependency):
                            st.error("Missing Required Dependencies")
                            st.code(str(e), language="text")
                        else:
                            st.error(f"Error: {e}")
                            logger.exception("Error running Novae")


# ==================== Results Page ====================
elif page == "Results":
    st.header("Results Visualization")

    if not st.session_state.novae_done:
        st.warning("Please run Novae first")
    else:
        # Use whichever dataset has Novae results
        adata = st.session_state.adata_filtered if st.session_state.adata_filtered is not None else st.session_state.adata

        # Visualization options
        st.subheader("Spatial Plots")

        color_by = st.selectbox(
            "Color by",
            options=["domain"] + list(adata.obs.columns),
            index=0,
        )

        _pt_size_results = st.slider("Point size", min_value=1, max_value=20, value=4, key="results_point_size")

        fig = viz.plot_spatial_scatter(adata, color_by=color_by, size=_pt_size_results)
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


# ==================== Domains & Markers Page ====================
elif page == "Domains & Markers":
    st.header("Domains & Markers")
    
    if st.session_state.adata is None:
        st.warning("Please load a dataset first")
    else:
        # Use filtered data if available, otherwise use original
        if st.session_state.adata_filtered is not None:
            adata = st.session_state.adata_filtered
        else:
            adata = st.session_state.adata
        
        st.markdown(
            "Explore your discovered groups: understand what each cluster/domain is, "
            "where it sits in tissue, and what genes/metadata define it."
        )
        
        # ========== 1. Label/Grouping Selector ==========
        st.subheader("1. Select Label Column")
        
        candidate_cols = cluster_interpretation.get_candidate_label_columns(adata)
        
        if len(candidate_cols) == 0:
            st.error("No columns found in adata.obs. Please check your dataset.")
            st.stop()
        
        label_col = st.selectbox(
            "Choose a grouping/label column",
            options=candidate_cols,
            help="Select a column containing cluster, domain, or cell type labels",
        )
        
        # Check for too many unique values
        n_unique = adata.obs[label_col].nunique()
        if n_unique > 200:
            st.warning(f"Column '{label_col}' has {n_unique} unique values. This may be slow.")
        
        # Exclude NA option
        exclude_na = st.checkbox("Exclude NA/missing values", value=True)
        
        # ========== 2. Grouped Summary Table ==========
        st.divider()
        st.subheader("2. Grouped Summary")
        
        # Detect sample column (used as default for cross-tabulation)
        sample_col = cluster_interpretation.utils.get_sample_column(adata)
        
        # Add column selectors in two columns
        col1, col2 = st.columns(2)
        
        with col1:
            st.caption("Primary grouping column selected above")
        
        with col2:
            # Add comparison column selector
            celltype_col_default = cluster_interpretation.utils.get_celltype_column(adata)
            comparison_options = ["None"] + candidate_cols
            
            # Set default index based on celltype column
            default_idx = 0
            if celltype_col_default and celltype_col_default in comparison_options:
                default_idx = comparison_options.index(celltype_col_default)
            
            comparison_col = st.selectbox(
                "Compare against column (optional)",
                options=comparison_options,
                index=default_idx,
                help="Create cross-tabulation with this column"
            )
            
            # Set to None if "None" selected
            if comparison_col == "None":
                comparison_col = None
        
        # Compute summary (with caching)
        @st.cache_data(show_spinner=False)
        def get_cluster_summary(_adata, _label_col, _sample_col, _exclude_na):
            return cluster_interpretation.compute_cluster_summary(
                _adata, _label_col, sample_col=_sample_col, exclude_na=_exclude_na
            )
        
        # Use comparison_col if selected, otherwise use sample_col
        crosstab_col = comparison_col if comparison_col else sample_col
        summary_df = get_cluster_summary(adata, label_col, crosstab_col, exclude_na)
        
        if len(summary_df) == 0:
            st.warning("No groups found in the selected column.")
            st.stop()
        
        st.dataframe(summary_df, use_container_width=True)
        
        # Group selection (using dropdown for simplicity)
        st.subheader("3. Select Group to Analyze")
        
        group_options = summary_df["group_id"].astype(str).tolist()
        selected_group = st.selectbox(
            "Select a group/cluster",
            options=group_options,
            index=0,
            help="Choose a specific group to analyze in detail",
        )
        
        # Convert back to original type
        original_type = adata.obs[label_col].dtype
        if pd.api.types.is_numeric_dtype(original_type):
            try:
                selected_group = type(adata.obs[label_col].iloc[0])(selected_group)
            except (ValueError, TypeError):
                logger.debug(
                    "Could not coerce selected group '%s' to dtype %s; keeping as string.",
                    selected_group,
                    original_type,
                )
        
        st.info(f"Analyzing group: **{selected_group}**")
        
        # ========== 3. Spatial Visualization ==========
        st.divider()
        st.subheader("4. Where in Tissue")
        
        # Check if spatial coordinates exist
        has_spatial, spatial_msg = cluster_interpretation.utils.check_spatial_coords(adata)
        
        if not has_spatial:
            st.warning(spatial_msg)
            st.info("Spatial visualization is not available. Continuing with marker analysis...")
        else:
            st.success(spatial_msg)
            
            # Sample/FOV filter (if applicable)
            if sample_col:
                samples = adata.obs[sample_col].unique().tolist()
                if len(samples) > 1:
                    selected_samples = st.multiselect(
                        "Filter by sample/FOV (optional)",
                        options=samples,
                        default=samples,
                        help="Select which samples to display",
                    )
                    
                    # Filter adata for visualization
                    if len(selected_samples) < len(samples):
                        sample_mask = adata.obs[sample_col].isin(selected_samples)
                        adata_viz = adata[sample_mask, :].copy()
                    else:
                        adata_viz = adata
                else:
                    adata_viz = adata
            else:
                adata_viz = adata
            
            # Create spatial plot
            try:
                fig_spatial = cluster_interpretation.plot_spatial_highlight(
                    adata_viz,
                    label_col=label_col,
                    group_id=selected_group,
                    spatial_key="spatial",
                    size=3,
                )
                st.plotly_chart(fig_spatial, use_container_width=True)
            except Exception as e:
                st.error(f"Error creating spatial plot: {e}")
        
        # ========== 4. Marker Genes ==========
        st.divider()
        st.subheader("5. What Defines It: Marker Genes")
        
        col1, col2 = st.columns(2)
        with col1:
            n_genes = st.number_input("Number of top genes", value=25, min_value=5, max_value=200)
        with col2:
            normalize = st.checkbox("Apply log1p normalization", value=True)
        
        # Compute marker genes (with caching)
        @st.cache_data(show_spinner=False)
        def get_marker_genes(_adata, _label_col, _group_id, _n_genes, _normalize):
            return cluster_interpretation.compute_marker_genes(
                _adata,
                label_col=_label_col,
                group_id=_group_id,
                n_genes=_n_genes,
                normalize=_normalize,
            )
        
        with st.spinner("Computing marker genes..."):
            try:
                markers_df = get_marker_genes(adata, label_col, selected_group, n_genes, normalize)
                
                if len(markers_df) == 0:
                    st.warning("No significant marker genes found for this group.")
                else:
                    st.success(f"Found {len(markers_df)} marker genes")
                    
                    # Display table
                    st.dataframe(markers_df, use_container_width=True)
                    
                    # Download button
                    csv = markers_df.to_csv(index=False)
                    st.download_button(
                        label="Download Marker Genes (CSV)",
                        data=csv,
                        file_name=f"markers_{label_col}_{selected_group}.csv",
                        mime="text/csv",
                    )
                    
                    # Store in adata.uns for export
                    if "xspatialnovae_markers" not in adata.uns:
                        adata.uns["xspatialnovae_markers"] = {}
                    if label_col not in adata.uns["xspatialnovae_markers"]:
                        adata.uns["xspatialnovae_markers"][label_col] = {}
                    
                    adata.uns["xspatialnovae_markers"][label_col][str(selected_group)] = markers_df.to_dict(orient="records")
                    
                    st.info("Marker genes stored in adata.uns['xspatialnovae_markers']")
            
            except Exception as e:
                st.error(f"Error computing marker genes: {e}")
                logger.exception("Error computing marker genes")
        
        # ========== 5. Gene Weights / Domain Drivers ==========
        st.divider()
        st.subheader("6. Gene Weights: Top Drivers for This Domain")
        st.caption(
            "Shows which genes drive this domain/cluster assignment. "
            "Positive logFC = enriched in this group; negative = depleted. "
            "When the 'novae' package is not installed, log fold change vs. "
            "all other cells is used as a proxy weight."
        )

        col1, col2 = st.columns(2)
        with col1:
            n_top_weights = st.number_input(
                "Top N genes (per direction)", value=20, min_value=1, max_value=200,
                key="n_top_weights",
            )
        with col2:
            weight_direction = st.selectbox(
                "Show",
                options=["All (top absolute)", "Positive only", "Negative only"],
                key="weight_direction",
            )

        @st.cache_data(show_spinner=False)
        def get_domain_weights(_adata, _label_col, _group_id, _n_top):
            return cluster_interpretation.compute_domain_weights(
                _adata,
                label_col=_label_col,
                group_id=_group_id,
                n_top=_n_top,
                store_in_uns=True,
            )

        with st.spinner("Computing gene weights..."):
            try:
                weights_df = get_domain_weights(adata, label_col, selected_group, n_top_weights)

                if weight_direction == "Positive only":
                    weights_df = weights_df[weights_df["direction"] == "positive"].reset_index(drop=True)
                elif weight_direction == "Negative only":
                    weights_df = weights_df[weights_df["direction"] == "negative"].reset_index(drop=True)

                if len(weights_df) == 0:
                    st.warning("No gene weights available for this group.")
                else:
                    st.dataframe(weights_df, use_container_width=True)

                    weights_csv = weights_df.to_csv(index=False)
                    st.download_button(
                        label="Export Gene Weights (CSV)",
                        data=weights_csv,
                        file_name=f"weights_{label_col}_{selected_group}.csv",
                        mime="text/csv",
                        key="dl_weights",
                    )

            except Exception as e:
                st.error(f"Error computing gene weights: {e}")
                logger.exception("Error computing domain weights")

        # ========== 6. Metadata Enrichment ==========
        st.divider()
        st.subheader("7. Metadata Enrichment")
        
        # Cell type composition
        celltype_col = cluster_interpretation.utils.get_celltype_column(adata)
        
        if celltype_col:
            st.markdown("**Cell Type Composition**")
            
            # Validate columns exist
            if label_col not in adata.obs.columns:
                st.error(f"Label column '{label_col}' not found in data")
            elif celltype_col not in adata.obs.columns:
                st.error(f"Cell type column '{celltype_col}' not found in data")
            else:
                try:
                    composition_df = cluster_interpretation.summaries.compute_group_composition(
                        adata, label_col, selected_group, celltype_col
                    )
                    
                    if len(composition_df) > 0:
                        col1, col2 = st.columns([1, 1])
                        
                        with col1:
                            st.dataframe(composition_df, use_container_width=True)
                        
                        with col2:
                            fig_comp = cluster_interpretation.plot_celltype_composition(
                                composition_df,
                                title=f"Cell Types in {selected_group}",
                            )
                            st.plotly_chart(fig_comp, use_container_width=True)
                    else:
                        st.info("No cell type data available for this group.")
                
                except Exception as e:
                    st.warning(f"Could not compute cell type composition: {e}")
        else:
            st.info("No cell type column detected. Skipping cell type composition.")
        
        # QC metrics comparison
        st.markdown("**QC Metrics Comparison**")
        
        qc_cols = cluster_interpretation.utils.get_qc_columns(adata)
        
        if len(qc_cols) > 0:
            try:
                qc_comparison_df = cluster_interpretation.summaries.compare_qc_metrics(
                    adata, label_col, selected_group, qc_cols
                )
                
                if len(qc_comparison_df) > 0:
                    display_df = cluster_interpretation.create_qc_comparison_table(qc_comparison_df)
                    st.dataframe(display_df, use_container_width=True)
                else:
                    st.info("No QC metrics available for comparison.")
            
            except Exception as e:
                st.warning(f"Could not compute QC comparison: {e}")
        else:
            st.info("No QC metric columns detected.")
        
        # Store cluster summary in adata.uns
        if "xspatialnovae_cluster_summary" not in adata.uns:
            adata.uns["xspatialnovae_cluster_summary"] = {}
        
        adata.uns["xspatialnovae_cluster_summary"][label_col] = summary_df.to_dict(orient="records")
        
        # Update session state to save changes
        if st.session_state.adata_filtered is not None:
            st.session_state.adata_filtered = adata
        else:
            st.session_state.adata = adata


# ==================== Export Page ====================
elif page == "Export":
    st.header("Export Results")

    if not st.session_state.novae_done:
        st.warning("Please run Novae first")
    else:
        # Use whichever dataset has Novae results
        adata = st.session_state.adata_filtered if st.session_state.adata_filtered is not None else st.session_state.adata

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
