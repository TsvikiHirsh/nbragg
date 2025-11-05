# Streamlit App Quick Start

Get up and running with the nbragg web app in 5 minutes!

## Prerequisites

- Python 3.9 or higher
- nbragg installed (`pip install nbragg` or `pip install -e .` from repo root)

## Installation

```bash
# Install optional streamlit dependencies
pip install -r requirements-streamlit.txt

# Or install manually
pip install streamlit plotly
```

## Run the App

From the **nbragg repository root** directory:

```bash
streamlit run docs/streamlit/streamlit_app.py
```

The app will open automatically in your browser at `http://localhost:8501`

If it doesn't open automatically, navigate to the URL shown in the terminal.

## First Analysis - Iron Powder Example

1. **Load Example Data** (Sidebar → 📁 1. Data Loading)
   - Keep "Example Data" selected
   - Select "iron_powder_transmission.csv"

2. **Define Cross-Section** (Sidebar → ⚛️ 2. Cross-Section Definition)
   - Keep default: 1 phase
   - Material: Fe_sg229_Iron-alpha (default)
   - Weight: 1.0
   - Temperature: 300 K

3. **Configure Model** (Sidebar → 🔧 3. Model Configuration)
   - ✅ Vary Background
   - ✅ Vary Response
   - Keep defaults: λ from 1.0 to 6.0 Å, Rietveld method

4. **Run Analysis**
   - Click **🚀 Run Analysis** button at top or bottom of sidebar
   - Wait for fitting to complete (~10-30 seconds)

5. **View Results**
   - **📊 Fit Results** tab: See data vs fit with residuals
   - **📈 Cross-Section** tab: View iron cross-section
   - **📉 Statistics** tab: Check χ²/dof and parameters

**Expected result**: χ²/dof ≈ 1-3 (good fit)

## Advanced Example - Multi-Phase with SANS

1. **Load Data**: Same as above

2. **Define Cross-Section**:
   - Number of Phases: **2**

   **Phase 1** (Iron):
   - Name: iron
   - Material: Fe_sg229_Iron-alpha
   - Weight: 0.7
   - Temperature: 300 K
   - Advanced → Enable SANS: ✅
   - SANS radius: 100 Å

   **Phase 2** (Oxide):
   - Name: oxide
   - Material: Al2O3_sg167_Corundum
   - Weight: 0.3
   - Temperature: 300 K

3. **Configure Model**:
   - ✅ Vary Background
   - ✅ Vary Response
   - ✅ Vary SANS

4. **Run Analysis** and view results

## Uploading Your Own Data

1. **Prepare Your Data**:
   - Format: CSV file
   - Required columns:
     - `wavelength` (Å)
     - `trans` (transmission, 0-1)
     - `err` (uncertainty in transmission)

   Example:
   ```csv
   wavelength,trans,err
   1.0,0.95,0.01
   1.1,0.94,0.01
   ...
   ```

2. **Upload**:
   - Sidebar → 📁 1. Data Loading
   - Select "Upload Your Data"
   - ✅ **Check consent box** (required)
   - Upload Transmission: Choose your CSV file

3. **Continue** with cross-section definition and fitting as usual

## Tips

### Faster Fitting
- Reduce wavelength range (e.g., 2.0 to 5.0 Å)
- Use fewer phases
- Try "least-squares" instead of "rietveld"

### Better Fits
- Enable SANS for samples with nanoscale features
- Enable Extinction for large crystallites
- Adjust wavelength range to focus on Bragg edges
- Try different materials from the database

### Exploring Materials
138 materials are available! Search for:
- **Iron**: Fe_sg229_Iron-alpha, Fe_sg225_Iron-gamma
- **Aluminum**: Al_sg225, Al2O3_sg167_Corundum
- **Steel components**: Ni_sg225, Cr_sg229
- **Oxides**: MgO_sg225, ZrO2_sg137
- And many more...

## Troubleshooting

### Port Already in Use
```bash
streamlit run docs/streamlit/streamlit_app.py --server.port 8502
```

### Cannot Find Example Data
Make sure you're running from the nbragg root directory:
```bash
cd /path/to/nbragg
streamlit run docs/streamlit/streamlit_app.py
```

### App Won't Start
Check all dependencies are installed:
```bash
python tests/test_streamlit_imports.py
```

### Fitting Fails
- Check your data format (wavelength, trans, err columns)
- Verify transmission values are between 0 and 1
- Ensure errors are positive
- Try reducing wavelength range

## Next Steps

- **Read full documentation**: [README.md](README.md)
- **Deploy your app**: [STREAMLIT_DEPLOYMENT.md](STREAMLIT_DEPLOYMENT.md)
- **Try Jupyter notebooks**: `notebooks/nbragg_tutorial.ipynb`
- **Explore the API**: [https://nbragg.readthedocs.io](https://nbragg.readthedocs.io)

## Support

- GitHub Issues: [https://github.com/your-username/nbragg/issues](https://github.com/your-username/nbragg/issues)
- Documentation: [https://nbragg.readthedocs.io](https://nbragg.readthedocs.io)

Happy analyzing! 🎉
