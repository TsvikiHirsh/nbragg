"""
nbragg Interactive Explorer - Streamlit App

Explore neutron Bragg edge transmission spectroscopy analysis interactively!
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import io
import json
from datetime import datetime
import hashlib

# Import nbragg
import nbragg

# Use Agg backend for matplotlib (non-interactive, for conversion to images)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save_user_data(data_content, data_type, user_id):
    """
    Save user-uploaded data to a private directory for analysis improvement.

    Parameters
    ----------
    data_content : str or bytes
        The data content to save
    data_type : str
        Type of data ('signal', 'openbeam', 'transmission')
    user_id : str
        Anonymous user identifier
    """
    try:
        # Create data collection directory if it doesn't exist
        data_dir = Path("user_data_collection")
        data_dir.mkdir(exist_ok=True)

        # Create timestamped filename
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = data_dir / f"{user_id}_{timestamp}_{data_type}.csv"

        # Save data
        if isinstance(data_content, bytes):
            with open(filename, 'wb') as f:
                f.write(data_content)
        else:
            with open(filename, 'w') as f:
                f.write(data_content)

        # Log metadata
        metadata_file = data_dir / "metadata.jsonl"
        metadata = {
            "timestamp": timestamp,
            "user_id": user_id,
            "data_type": data_type,
            "filename": str(filename)
        }
        with open(metadata_file, 'a') as f:
            f.write(json.dumps(metadata) + "\n")

        return True
    except Exception as e:
        st.warning(f"Could not save data for research purposes: {e}")
        return False

def get_anonymous_user_id():
    """Generate anonymous user ID based on session"""
    if 'user_id' not in st.session_state:
        # Use session state to generate a consistent ID per session
        st.session_state.user_id = hashlib.md5(
            str(datetime.now().timestamp()).encode()
        ).hexdigest()[:12]
    return st.session_state.user_id

def mpl_to_plotly(fig):
    """
    Convert matplotlib figure to plotly figure for interactivity.
    """
    axes = fig.get_axes()

    # Create plotly figure
    if len(axes) == 1:
        plotly_fig = go.Figure()
        ax = axes[0]

        # Extract lines and error bars
        for collection in ax.collections:
            if hasattr(collection, 'get_paths') and len(collection.get_paths()) > 0:
                facecolors = collection.get_facecolor()
                if len(facecolors) > 0:
                    for path in collection.get_paths():
                        vertices = path.vertices
                        if len(vertices) > 0:
                            plotly_fig.add_trace(go.Scatter(
                                x=vertices[:, 0],
                                y=vertices[:, 1],
                                mode='lines',
                                fill='toself',
                                fillcolor=matplotlib.colors.to_hex(facecolors[0], keep_alpha=True),
                                line=dict(width=0),
                                opacity=0.2,
                                showlegend=False,
                                hoverinfo='skip'
                            ))

        for line in ax.get_lines():
            xdata = line.get_xdata()
            ydata = line.get_ydata()
            plotly_fig.add_trace(go.Scatter(
                x=xdata, y=ydata,
                mode='lines',
                name=line.get_label(),
                line=dict(color=matplotlib.colors.rgb2hex(line.get_color())),
                showlegend=(line.get_label() and not line.get_label().startswith('_'))
            ))

        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        plotly_fig.update_layout(
            xaxis_title=ax.get_xlabel(),
            yaxis_title=ax.get_ylabel(),
            title=ax.get_title(),
            hovermode='x unified',
            template='plotly_dark',
            xaxis=dict(range=xlim),
            yaxis=dict(range=ylim),
            paper_bgcolor='#0E1117',
            plot_bgcolor='#0E1117'
        )

    elif len(axes) == 2:
        # Two subplots (data + residuals)
        plotly_fig = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.05,
            row_heights=[0.7, 0.3]
        )

        ax_top = axes[0]
        for collection in ax_top.collections:
            if hasattr(collection, 'get_paths') and len(collection.get_paths()) > 0:
                facecolors = collection.get_facecolor()
                if len(facecolors) > 0:
                    for path in collection.get_paths():
                        vertices = path.vertices
                        if len(vertices) > 0:
                            plotly_fig.add_trace(go.Scatter(
                                x=vertices[:, 0],
                                y=vertices[:, 1],
                                mode='lines',
                                fill='toself',
                                fillcolor=matplotlib.colors.to_hex(facecolors[0], keep_alpha=True),
                                line=dict(width=0),
                                opacity=0.2,
                                showlegend=False,
                                hoverinfo='skip'
                            ), row=1, col=1)

        for line in ax_top.get_lines():
            xdata = line.get_xdata()
            ydata = line.get_ydata()
            plotly_fig.add_trace(go.Scatter(
                x=xdata, y=ydata,
                mode='lines',
                name=line.get_label(),
                line=dict(color=matplotlib.colors.rgb2hex(line.get_color())),
                showlegend=(line.get_label() and not line.get_label().startswith('_'))
            ), row=1, col=1)

        ax_bottom = axes[1]
        for line in ax_bottom.get_lines():
            xdata = line.get_xdata()
            ydata = line.get_ydata()
            plotly_fig.add_trace(go.Scatter(
                x=xdata, y=ydata,
                mode='lines',
                name=line.get_label(),
                line=dict(color=matplotlib.colors.rgb2hex(line.get_color())),
                showlegend=False
            ), row=2, col=1)

        top_ylim = ax_top.get_ylim()
        bottom_ylim = ax_bottom.get_ylim()
        xlim = ax_top.get_xlim()

        plotly_fig.update_xaxes(title_text=ax_bottom.get_xlabel(), range=xlim, row=2, col=1)
        plotly_fig.update_xaxes(range=xlim, row=1, col=1)
        plotly_fig.update_yaxes(title_text=ax_top.get_ylabel(), range=top_ylim, row=1, col=1)
        plotly_fig.update_yaxes(title_text=ax_bottom.get_ylabel(), range=bottom_ylim, row=2, col=1)
        plotly_fig.update_layout(
            title=ax_top.get_title(),
            hovermode='x unified',
            template='plotly_dark',
            height=600,
            paper_bgcolor='#0E1117',
            plot_bgcolor='#0E1117'
        )

    else:
        plotly_fig = go.Figure()

    return plotly_fig

# Page config
st.set_page_config(
    page_title="nbragg Explorer",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for dark theme
st.markdown("""
<style>
.main > div {
    padding-top: 2rem;
}
.stPlotly {
    background-color: #0E1117;
}
h1 {
    color: #58D3F7;
}
h2 {
    color: #FFA07A;
}
h3 {
    color: #98D8C8;
}
.stTabs [data-baseweb="tab-list"] {
    gap: 24px;
}
.stTabs [data-baseweb="tab"] {
    height: 50px;
    padding-left: 20px;
    padding-right: 20px;
    background-color: #262730;
}
.stTabs [data-baseweb="tab"][aria-selected="true"] {
    background-color: #3D5A80;
}
.stButton>button {
    width: 100%;
}
</style>
""", unsafe_allow_html=True)

# Title
st.title("⚛️ nbragg Interactive Explorer")
st.markdown("Explore neutron Bragg edge transmission spectroscopy analysis step by step!")

# Initialize session state
if 'fit_result' not in st.session_state:
    st.session_state.fit_result = None
if 'cross_section' not in st.session_state:
    st.session_state.cross_section = None
if 'model' not in st.session_state:
    st.session_state.model = None
if 'data' not in st.session_state:
    st.session_state.data = None
if 'data_consent' not in st.session_state:
    st.session_state.data_consent = False

# Sidebar
st.sidebar.header("⚙️ Analysis Pipeline")

# Run button at the top
run_button = st.sidebar.button("🚀 Run Analysis", type="primary", use_container_width=True)

st.sidebar.markdown("---")
st.sidebar.markdown("Configure each stage of the analysis pipeline:")

# Stage 1: Data Loading
with st.sidebar.expander("📁 1. Data Loading", expanded=True):
    data_source = st.radio(
        "Data Source",
        ["Example Data", "Upload Your Data"],
        help="Use example data or upload your own transmission measurements"
    )

    if data_source == "Example Data":
        example_file = st.selectbox(
            "Select Example",
            ["iron_powder_transmission.csv", "iron_powder.csv + openbeam.csv"],
            help="Choose from available example datasets"
        )

        if example_file == "iron_powder_transmission.csv":
            data_type = "transmission"
            signal_path = "notebooks/iron_powder_transmission.csv"
            openbeam_path = None
        else:
            data_type = "signal_openbeam"
            signal_path = "notebooks/iron_powder.csv"
            openbeam_path = "notebooks/openbeam.csv"

    else:  # Upload Your Data
        st.markdown("**📋 Data Usage Consent**")
        st.info(
            "By uploading data to this application, you consent to its use for "
            "improving the nbragg software and enhancing user experience. "
            "All data is stored securely and used solely for research and development purposes. "
            "Your data will be anonymized and treated confidentially."
        )

        data_consent = st.checkbox(
            "I consent to my uploaded data being used for research purposes",
            value=st.session_state.data_consent
        )
        st.session_state.data_consent = data_consent

        if data_consent:
            upload_type = st.radio(
                "Upload Type",
                ["Transmission", "Signal + Openbeam"],
                help="Upload pre-calculated transmission or separate signal and openbeam files"
            )

            if upload_type == "Transmission":
                uploaded_trans = st.file_uploader(
                    "Upload Transmission CSV",
                    type=['csv'],
                    help="CSV file with columns: wavelength, trans, err"
                )

                if uploaded_trans is not None:
                    # Save uploaded data
                    user_id = get_anonymous_user_id()
                    save_user_data(uploaded_trans.getvalue(), "transmission", user_id)

                    # Save temporarily
                    temp_path = Path(f"/tmp/uploaded_trans_{user_id}.csv")
                    with open(temp_path, 'wb') as f:
                        f.write(uploaded_trans.getvalue())

                    data_type = "transmission"
                    signal_path = str(temp_path)
                    openbeam_path = None
                    st.success("✓ Transmission data uploaded")
                else:
                    signal_path = None
                    openbeam_path = None

            else:  # Signal + Openbeam
                uploaded_signal = st.file_uploader(
                    "Upload Signal CSV",
                    type=['csv'],
                    help="CSV file with sample signal"
                )
                uploaded_openbeam = st.file_uploader(
                    "Upload Openbeam CSV",
                    type=['csv'],
                    help="CSV file with openbeam (no sample)"
                )

                if uploaded_signal is not None and uploaded_openbeam is not None:
                    # Save uploaded data
                    user_id = get_anonymous_user_id()
                    save_user_data(uploaded_signal.getvalue(), "signal", user_id)
                    save_user_data(uploaded_openbeam.getvalue(), "openbeam", user_id)

                    # Save temporarily
                    temp_signal_path = Path(f"/tmp/uploaded_signal_{user_id}.csv")
                    temp_openbeam_path = Path(f"/tmp/uploaded_openbeam_{user_id}.csv")

                    with open(temp_signal_path, 'wb') as f:
                        f.write(uploaded_signal.getvalue())
                    with open(temp_openbeam_path, 'wb') as f:
                        f.write(uploaded_openbeam.getvalue())

                    data_type = "signal_openbeam"
                    signal_path = str(temp_signal_path)
                    openbeam_path = str(temp_openbeam_path)
                    st.success("✓ Signal and openbeam data uploaded")
                else:
                    signal_path = None
                    openbeam_path = None
        else:
            signal_path = None
            openbeam_path = None
            st.warning("Please consent to data usage to upload files")

# Stage 2: Cross-Section Definition
with st.sidebar.expander("⚛️ 2. Cross-Section Definition", expanded=False):
    st.markdown("**Material Selection**")

    # Get list of available materials from the materials dict
    material_list = sorted(list(nbragg.materials.keys()))

    # Number of phases
    n_phases = st.number_input(
        "Number of Phases",
        min_value=1,
        max_value=5,
        value=1,
        help="Number of material phases to combine"
    )

    # Store phase configurations
    phase_configs = {}

    for i in range(n_phases):
        st.markdown(f"**Phase {i+1}**")

        phase_name = st.text_input(
            f"Phase {i+1} Name",
            value=f"phase{i+1}",
            key=f"phase_name_{i}"
        )

        material = st.selectbox(
            f"Material",
            options=material_list,
            index=material_list.index("Fe_sg229_Iron-alpha") if "Fe_sg229_Iron-alpha" in material_list else 0,
            key=f"material_{i}",
            help="Select material from nbragg.materials"
        )

        col1, col2 = st.columns(2)
        with col1:
            weight = st.number_input(
                "Weight",
                min_value=0.0,
                max_value=1.0,
                value=1.0/n_phases,
                step=0.01,
                key=f"weight_{i}"
            )

        with col2:
            temp = st.number_input(
                "Temp (K)",
                min_value=0.0,
                max_value=1000.0,
                value=300.0,
                step=10.0,
                key=f"temp_{i}"
            )

        # Optional parameters
        with st.expander(f"Advanced (Phase {i+1})"):
            enable_sans = st.checkbox(f"Enable SANS", key=f"sans_enable_{i}")
            sans = st.number_input(
                "SANS radius (Å)",
                min_value=0.0,
                max_value=1000.0,
                value=100.0,
                step=10.0,
                key=f"sans_{i}",
                disabled=not enable_sans
            ) if enable_sans else None

            enable_ext = st.checkbox(f"Enable Extinction", key=f"ext_enable_{i}")
            if enable_ext:
                ext_l = st.number_input("ext_l", value=100.0, key=f"ext_l_{i}")
                ext_Gg = st.number_input("ext_Gg", value=100.0, key=f"ext_Gg_{i}")
                ext_L = st.number_input("ext_L", value=100000.0, key=f"ext_L_{i}")
            else:
                ext_l, ext_Gg, ext_L = None, None, None

        # Store configuration
        phase_config = {
            "mat": nbragg.materials[material],
            "weight": weight,
            "temp": temp
        }
        if sans is not None:
            phase_config["sans"] = sans
        if ext_l is not None:
            phase_config["ext_l"] = ext_l
            phase_config["ext_Gg"] = ext_Gg
            phase_config["ext_L"] = ext_L

        phase_configs[phase_name] = phase_config

# Stage 3: Model Configuration
with st.sidebar.expander("🔧 3. Model Configuration", expanded=False):
    st.markdown("**Fitting Options**")

    vary_background = st.checkbox("Vary Background", value=True)
    vary_response = st.checkbox("Vary Response", value=True)
    vary_sans = st.checkbox("Vary SANS", value=False)
    vary_extinction = st.checkbox("Vary Extinction", value=False)
    vary_lattice = st.checkbox("Vary Lattice", value=False)
    vary_orientation = st.checkbox("Vary Orientation", value=False)

    st.markdown("**Wavelength Range**")
    wlmin = st.number_input("Min λ (Å)", value=1.0, step=0.1)
    wlmax = st.number_input("Max λ (Å)", value=6.0, step=0.1)

    st.markdown("**Fitting Method**")
    fit_method = st.selectbox(
        "Method",
        ["rietveld", "staged", "least-squares"],
        help="Rietveld: accumulative stages, Staged: frozen stages, Least-squares: single stage"
    )

    # Optional: Custom stages
    use_custom_stages = st.checkbox("Custom Stages", value=False)
    if use_custom_stages:
        st.text_area(
            "Stages (JSON)",
            value='{"basic": "basic", "background": "background", "response": "response"}',
            help="Define custom fitting stages as JSON"
        )

# Run button at the bottom
st.sidebar.markdown("---")
run_button_bottom = st.sidebar.button("🚀 Run Analysis", type="primary", use_container_width=True, key="run_bottom")

# Main content area
if run_button or run_button_bottom:
    if signal_path is None:
        st.error("Please load data first!")
    else:
        with st.spinner("Running analysis pipeline..."):
            try:
                # Step 1: Load data
                if data_type == "transmission":
                    data = nbragg.Data.from_transmission(signal_path)
                else:
                    data = nbragg.Data(signal_path, openbeam_path)

                st.session_state.data = data
                st.sidebar.success("✓ Data loaded")

                # Step 2: Create cross-section
                xs = nbragg.CrossSection(**phase_configs)
                st.session_state.cross_section = xs
                st.sidebar.success(f"✓ Cross-section created ({n_phases} phase{'s' if n_phases > 1 else ''})")

                # Step 3: Create model
                model = nbragg.TransmissionModel(
                    xs,
                    vary_background=vary_background,
                    vary_response=vary_response,
                    vary_sans=vary_sans,
                    vary_extinction=vary_extinction,
                    vary_lattice=vary_lattice,
                    vary_orientation=vary_orientation
                )
                st.session_state.model = model
                st.sidebar.success("✓ Model created")

                # Step 4: Fit
                result = model.fit(data, wlmin=wlmin, wlmax=wlmax, method=fit_method)
                st.session_state.fit_result = result
                st.sidebar.success(f"✓ Fit completed (χ²/dof: {result.redchi:.2f})")

            except Exception as e:
                st.error(f"Error during analysis: {e}")
                st.exception(e)

# Display results
if st.session_state.fit_result is not None:
    result = st.session_state.fit_result

    # Create tabs
    tab1, tab2, tab3 = st.tabs(["📊 Fit Results", "📈 Cross-Section", "📉 Statistics"])

    with tab1:
        st.header("Fit Results")

        try:
            # Generate plot
            mpl_fig = result.plot()

            # Convert to plotly
            plotly_fig = mpl_to_plotly(mpl_fig)
            st.plotly_chart(plotly_fig, use_container_width=True)
            plt.close(mpl_fig)

        except Exception as e:
            st.error(f"Error generating plot: {e}")

    with tab2:
        st.header("Cross-Section Components")

        try:
            xs = st.session_state.cross_section
            mpl_fig = xs.plot()

            plotly_fig = mpl_to_plotly(mpl_fig)
            st.plotly_chart(plotly_fig, use_container_width=True)
            plt.close(mpl_fig)

        except Exception as e:
            st.error(f"Error generating cross-section plot: {e}")

    with tab3:
        st.header("Fit Statistics")

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Fit Quality")
            st.metric("Reduced χ²", f"{result.redchi:.4f}")
            st.metric("AIC", f"{result.aic:.2f}")
            st.metric("BIC", f"{result.bic:.2f}")

            # Quality indicator
            if result.redchi < 2:
                st.success("✅ Excellent fit")
            elif result.redchi < 5:
                st.info("ℹ️ Good fit")
            else:
                st.warning("⚠️ Poor fit")

        with col2:
            st.subheader("Fit Parameters")

            # Show parameters table
            params_data = []
            for param_name, param in result.params.items():
                stderr_str = f"{param.stderr:.4e}" if param.stderr is not None else "N/A"
                params_data.append({
                    'Parameter': param_name,
                    'Value': f"{param.value:.4e}",
                    'Stderr': stderr_str,
                    'Vary': 'Yes' if param.vary else 'No'
                })

            params_df = pd.DataFrame(params_data)
            st.dataframe(params_df, hide_index=True, use_container_width=True, height=400)

        # Full fit report
        with st.expander("📋 View Full Fit Report"):
            fit_report = result.fit_report()
            st.text(fit_report)

else:
    st.info("👈 Configure the analysis pipeline in the sidebar and click 'Run Analysis' to start!")

    # Show example usage
    with st.expander("📚 How to use this app"):
        st.markdown("""
        ### Quick Start Guide

        1. **Load Data**
           - Choose example data or upload your own
           - For uploads, you must consent to data usage for research purposes

        2. **Define Cross-Section**
           - Select materials from the nbragg database
           - Set phase fractions and temperatures
           - Optionally enable SANS and extinction corrections

        3. **Configure Model**
           - Toggle which parameters to vary during fitting
           - Set wavelength range for analysis
           - Choose fitting method (Rietveld recommended)

        4. **Run Analysis**
           - Click the "Run Analysis" button
           - View interactive plots and fit statistics

        ### Example Workflow

        - **Simple iron sample**: Use "iron_powder_transmission.csv", select "Fe_sg229_Iron-alpha", enable background and response fitting
        - **Complex multi-phase**: Add multiple phases with different weights, enable SANS/extinction for accurate modeling
        """)

# Footer
st.markdown("---")
st.markdown(
    "**nbragg Interactive Explorer** | "
    "[Documentation](https://github.com/yourusername/nbragg) | "
    "[Report Issues](https://github.com/yourusername/nbragg/issues)"
)
