# nbragg Streamlit App

Interactive web application for neutron Bragg edge transmission spectroscopy analysis using nbragg.

## Features

- **Dark Theme Interface**: Modern dark theme optimized for long analysis sessions
- **Interactive Data Loading**: Upload your own data or use example datasets
- **Data Collection for Research**: Optional anonymous data collection to improve nbragg (with user consent)
- **Multi-Phase Cross-Section Builder**: Combine multiple materials with custom parameters
- **Advanced Modeling**: Support for SANS, extinction, lattice parameters, and orientation
- **Interactive Plotly Plots**: Zoom, pan, and explore your results interactively
- **Comprehensive Statistics**: View fit quality metrics and parameter uncertainties
- **Rietveld Refinement**: Multi-stage fitting with accumulative parameters

## Installation

```bash
# Install required packages
pip install streamlit plotly nbragg
```

## Running the App

```bash
# From the nbragg directory
streamlit run streamlit_app.py

# Or specify port
streamlit run streamlit_app.py --server.port 8501
```

## Usage

### 1. Data Loading

**Example Data**:
- `iron_powder_transmission.csv`: Pre-calculated transmission data
- `iron_powder.csv + openbeam.csv`: Separate signal and openbeam files

**Upload Your Data**:
- You must consent to data usage for research purposes
- Uploaded data is anonymized and stored securely
- Helps improve nbragg for the community

### 2. Cross-Section Definition

- Select materials from the nbragg database (e.g., Fe_sg229_Iron-alpha)
- Add multiple phases with custom weights
- Set temperature for each phase
- Enable SANS corrections with hard-sphere radius
- Enable extinction corrections with crystallite size parameters

### 3. Model Configuration

**Vary Options**:
- **Background**: Polynomial background correction
- **Response**: Instrument response function
- **SANS**: Small-angle neutron scattering
- **Extinction**: Primary and secondary extinction
- **Lattice**: Lattice parameters (a, b, c)
- **Orientation**: Texture and preferred orientation

**Fitting Method**:
- **Rietveld**: True Rietveld refinement with accumulative parameters (recommended)
- **Staged**: Sequential refinement with frozen parameters
- **Least-squares**: Single-stage fitting

### 4. Results

View three tabs:
- **Fit Results**: Interactive plot showing data, fit, and residuals
- **Cross-Section**: Plot of individual phase contributions
- **Statistics**: Fit quality metrics and parameter table

## Data Collection

The app includes optional anonymous data collection to improve nbragg:

- **What is collected**: Uploaded transmission data (CSV files)
- **How it's used**: To improve fitting algorithms and add new features
- **Privacy**: Data is anonymized with random session IDs
- **Storage**: Saved locally in `user_data_collection/` directory
- **Consent**: Explicit user consent required before any upload

To disable data collection, simply don't check the consent box.

## Customization

### Email Notifications

To receive email notifications when users upload data, you can modify the `save_user_data()` function to send emails:

```python
import smtplib
from email.mime.text import MIMEText

def send_notification_email(metadata):
    msg = MIMEText(f"New data uploaded: {metadata}")
    msg['Subject'] = 'nbragg App: New Data Upload'
    msg['From'] = 'your-email@example.com'
    msg['To'] = 'your-email@example.com'

    with smtplib.SMTP('smtp.gmail.com', 587) as server:
        server.starttls()
        server.login('your-email@example.com', 'your-app-password')
        server.send_message(msg)
```

### Cloud Storage

For production deployment, consider using cloud storage for uploaded data:

```python
# Example using AWS S3
import boto3

def save_to_s3(data_content, filename):
    s3 = boto3.client('s3')
    s3.put_object(
        Bucket='nbragg-user-data',
        Key=filename,
        Body=data_content
    )
```

## Deployment

### Local Development

```bash
streamlit run streamlit_app.py
```

### Streamlit Cloud (Recommended)

1. Push your code to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your GitHub repository
4. Deploy!

### Docker

```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "streamlit_app.py", "--server.port=8501", "--server.address=0.0.0.0"]
```

## Troubleshooting

### Port Already in Use

```bash
streamlit run streamlit_app.py --server.port 8502
```

### Example Data Not Found

Make sure you're running from the nbragg directory where `notebooks/` exists.

### Slow Fitting

- Reduce wavelength range (wlmin, wlmax)
- Use fewer phases
- Disable advanced parameters (SANS, extinction) if not needed

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

Same as nbragg main package.
