# nbragg Streamlit Web App

Interactive web application for neutron Bragg edge transmission spectroscopy analysis.

![nbragg App](https://img.shields.io/badge/streamlit-app-FF4B4B?style=flat&logo=streamlit)

## Quick Start

### Installation

The streamlit app is **optional** and requires additional dependencies:

```bash
# Install streamlit dependencies
pip install streamlit plotly

# Or install from requirements file
pip install -r requirements-streamlit.txt
```

### Run the App

```bash
# From the nbragg root directory
streamlit run docs/streamlit/streamlit_app.py

# Or navigate to the streamlit directory first
cd docs/streamlit
streamlit run streamlit_app.py
```

The app will open in your browser at `http://localhost:8501`

## Features

### 🎨 Modern Dark Theme Interface
- Professional dark theme optimized for extended use
- Interactive Plotly plots with zoom, pan, and hover
- Collapsible sidebar for workflow management

### 📁 Flexible Data Loading
- **Example datasets**: Pre-loaded iron powder data from notebooks
- **Upload your own**: Support for transmission or signal+openbeam formats
- **Data consent**: Optional anonymous data collection for research (with explicit user consent)

### ⚛️ Advanced Cross-Section Modeling
- Select from all materials in `nbragg.materials` database
- Multi-phase support (up to 5 phases)
- Configure phase fractions, temperatures
- Enable SANS corrections (hard-sphere model)
- Enable extinction corrections (primary/secondary)
- Lattice parameter refinement

### 🔧 Comprehensive Model Configuration
- Toggle fitting parameters:
  - Background correction
  - Instrument response function
  - SANS parameters
  - Extinction parameters
  - Lattice parameters
  - Texture/orientation
- Choose fitting method:
  - **Rietveld**: True Rietveld refinement with accumulative parameters (recommended)
  - **Staged**: Sequential refinement with frozen parameters
  - **Least-squares**: Single-stage fitting
- Custom wavelength range
- Optional custom fitting stages

### 📊 Interactive Results
Three result tabs:
1. **Fit Results**: Data vs fit with residuals plot
2. **Cross-Section**: Individual phase contributions
3. **Statistics**: Fit quality metrics and parameter table with uncertainties

## Screenshots

### Main Interface
```
┌─────────────────────────────────────────────────────────────┐
│  ⚛️ nbragg Interactive Explorer                             │
├──────────────┬──────────────────────────────────────────────┤
│  Sidebar     │                                              │
│              │         Interactive Plotly Plot              │
│ 📁 Data      │                                              │
│ ⚛️ X-Section │                                              │
│ 🔧 Model     │                                              │
│              │                                              │
│ 🚀 Run       ├──────────────────────────────────────────────┤
│              │  📊 Fit | 📈 X-Sect | 📉 Stats              │
└──────────────┴──────────────────────────────────────────────┘
```

## Data Collection (Optional)

The app includes **optional anonymous data collection** to improve nbragg:

### What is Collected
- Uploaded transmission/signal data (CSV files only)
- Anonymous session ID (randomly generated)
- Upload timestamp and data type

### What is NOT Collected
- No IP addresses
- No personal information
- No browsing history
- No system information

### How to Enable
Users must **explicitly consent** by checking the consent box before uploading data.

### Where Data is Stored
By default: `user_data_collection/` directory
- Individual files: `<user_id>_<timestamp>_<type>.csv`
- Metadata log: `metadata.jsonl`

### Alternative Storage Options
See [STREAMLIT_DEPLOYMENT.md](STREAMLIT_DEPLOYMENT.md) for:
- Email notifications
- AWS S3 cloud storage
- PostgreSQL database
- Custom webhooks

### Privacy & GDPR
- Data is anonymized (no personal identifiers)
- Used solely for improving nbragg
- Can be deleted upon request
- See full privacy considerations in deployment guide

## Documentation

- **[STREAMLIT_README.md](STREAMLIT_README.md)** - User guide and features
- **[STREAMLIT_DEPLOYMENT.md](STREAMLIT_DEPLOYMENT.md)** - Deployment options and data collection setup

## Deployment Options

### 1. Streamlit Cloud (Recommended for Public Apps)
- Free hosting
- Automatic deployments from GitHub
- Built-in authentication
- Easy setup at [share.streamlit.io](https://share.streamlit.io)

### 2. Self-Hosted (Recommended for Private/Research Use)
- Full control over data
- Private data collection
- Custom domain
- See deployment guide for setup

### 3. Docker
- Containerized deployment
- Easy scaling
- Cloud-agnostic (AWS, GCP, Azure)

### 4. Heroku
- Free tier available
- Simple git-based deployment
- Automatic HTTPS

Full deployment instructions in [STREAMLIT_DEPLOYMENT.md](STREAMLIT_DEPLOYMENT.md).

## Requirements

### Core Dependencies
```
streamlit>=1.28.0
plotly>=5.17.0
```

These are **not** included in the main nbragg requirements to keep the core package lightweight.

### Installation
```bash
# Create a requirements file for streamlit
cat > requirements-streamlit.txt <<EOF
streamlit>=1.28.0
plotly>=5.17.0
EOF

# Install
pip install -r requirements-streamlit.txt
```

Or add to your existing environment:
```bash
pip install streamlit plotly
```

## Troubleshooting

### Error: "cannot import name streamlit"
**Solution**: Install streamlit: `pip install streamlit plotly`

### Error: "File not found: iron_powder_transmission.csv"
**Solution**: Run from nbragg root directory: `streamlit run docs/streamlit/streamlit_app.py`

### Error: "KeyError: 'clear'"
**Solution**: This has been fixed in the latest version. Update to the latest streamlit_app.py.

### App is slow when fitting
**Solutions**:
- Reduce wavelength range (wlmin/wlmax)
- Use fewer phases
- Disable SANS/extinction if not needed
- Try "least-squares" instead of "rietveld" method

### Cannot upload files
**Solutions**:
- Check you've consented to data usage
- Verify file is CSV format
- Check file size (max 10MB recommended)
- Ensure proper CSV structure (wavelength, trans, err columns)

## Development

### Testing Locally
```bash
# Test imports
python tests/test_streamlit_imports.py

# Run app in development mode
streamlit run docs/streamlit/streamlit_app.py --server.runOnSave true
```

### Customization
The app can be customized by editing `streamlit_app.py`:
- **Theme colors**: Modify CSS in `st.markdown()` section
- **Default parameters**: Change default values in sidebar widgets
- **Available materials**: Controlled by `nbragg.materials`
- **Data collection**: Modify `save_user_data()` function

## Contributing

Contributions welcome! Areas for improvement:
- [ ] Additional plot types (d-spacing, energy)
- [ ] Batch processing multiple files
- [ ] Export results to CSV/Excel
- [ ] Parameter comparison tool
- [ ] Material database browser
- [ ] Advanced fitting diagnostics
- [ ] Real-time fitting progress

Please open an issue or submit a pull request.

## Support

- **GitHub Issues**: [nbragg/issues](https://github.com/your-username/nbragg/issues)
- **Documentation**: [nbragg docs](https://github.com/your-username/nbragg)
- **Tutorial**: See Jupyter notebooks in `notebooks/` directory

## License

Same as nbragg main package (see root LICENSE file).

## Citation

If you use nbragg in your research, please cite:
```bibtex
@software{nbragg,
  title = {nbragg: Neutron Bragg Edge Analysis},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/your-username/nbragg}
}
```
