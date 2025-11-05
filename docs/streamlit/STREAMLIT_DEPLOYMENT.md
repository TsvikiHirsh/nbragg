# Streamlit App Deployment Guide

## Quick Start

```bash
cd /work/Programs/nbragg
streamlit run streamlit_app.py
```

The app will open in your browser at `http://localhost:8501`

## Data Collection Setup

### Option 1: Local Storage (Default)

Data is saved to `user_data_collection/` directory:
- `user_data_collection/<user_id>_<timestamp>_<type>.csv` - Uploaded data files
- `user_data_collection/metadata.jsonl` - Metadata log

Check collected data:
```bash
ls -lh user_data_collection/
cat user_data_collection/metadata.jsonl
```

### Option 2: Email Notifications

Add this to `streamlit_app.py` after the `save_user_data()` function:

```python
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

def send_email_notification(filepath, metadata):
    """Send email when new data is uploaded"""
    sender_email = "your-email@gmail.com"  # Replace with your email
    receiver_email = "your-email@gmail.com"  # Replace with your email
    password = "your-app-password"  # Use Gmail App Password

    message = MIMEMultipart()
    message["From"] = sender_email
    message["To"] = receiver_email
    message["Subject"] = f"nbragg App: New {metadata['data_type']} Data Uploaded"

    body = f"""
    New data uploaded to nbragg app:

    User ID: {metadata['user_id']}
    Timestamp: {metadata['timestamp']}
    Data Type: {metadata['data_type']}
    Filename: {metadata['filename']}
    """

    message.attach(MIMEText(body, "plain"))

    # Attach the data file
    with open(filepath, "rb") as attachment:
        part = MIMEBase("application", "octet-stream")
        part.set_payload(attachment.read())

    encoders.encode_base64(part)
    part.add_header(
        "Content-Disposition",
        f"attachment; filename= {Path(filepath).name}",
    )
    message.attach(part)

    try:
        with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
            server.login(sender_email, password)
            server.sendmail(sender_email, receiver_email, message.as_string())
        return True
    except Exception as e:
        print(f"Failed to send email: {e}")
        return False

# Then modify save_user_data() to call this:
def save_user_data(data_content, data_type, user_id):
    # ... existing code ...

    if True:  # If save successful
        # Send email notification
        send_email_notification(filename, metadata)

    return True
```

**Gmail App Password Setup:**
1. Go to Google Account settings
2. Security → 2-Step Verification (enable if not enabled)
3. App passwords → Generate new app password
4. Use this password in the code (not your regular Gmail password)

### Option 3: Cloud Storage (AWS S3)

```python
import boto3
from botocore.exceptions import ClientError

def save_to_s3(data_content, data_type, user_id):
    """Save uploaded data to AWS S3"""
    s3_client = boto3.client('s3',
        aws_access_key_id='YOUR_ACCESS_KEY',
        aws_secret_access_key='YOUR_SECRET_KEY',
        region_name='us-east-1'
    )

    bucket_name = 'nbragg-user-uploads'
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    object_name = f"{user_id}_{timestamp}_{data_type}.csv"

    try:
        if isinstance(data_content, bytes):
            s3_client.put_object(
                Bucket=bucket_name,
                Key=object_name,
                Body=data_content,
                Metadata={
                    'user-id': user_id,
                    'data-type': data_type,
                    'timestamp': timestamp
                }
            )
        else:
            s3_client.put_object(
                Bucket=bucket_name,
                Key=object_name,
                Body=data_content.encode('utf-8'),
                Metadata={
                    'user-id': user_id,
                    'data-type': data_type,
                    'timestamp': timestamp
                }
            )
        return True
    except ClientError as e:
        st.warning(f"Could not save to S3: {e}")
        return False
```

**AWS Setup:**
1. Create S3 bucket: `nbragg-user-uploads`
2. Create IAM user with S3 write permissions
3. Get access key and secret key
4. Install: `pip install boto3`

### Option 4: Database Storage (PostgreSQL)

```python
import psycopg2
from psycopg2.extras import Json

def save_to_database(data_content, data_type, user_id):
    """Save uploaded data to PostgreSQL database"""
    try:
        conn = psycopg2.connect(
            host="localhost",
            database="nbragg_data",
            user="your_username",
            password="your_password"
        )
        cur = conn.cursor()

        # Create table if not exists
        cur.execute("""
            CREATE TABLE IF NOT EXISTS uploads (
                id SERIAL PRIMARY KEY,
                user_id VARCHAR(50),
                timestamp TIMESTAMP,
                data_type VARCHAR(50),
                data_content TEXT,
                metadata JSON
            )
        """)

        # Insert data
        cur.execute("""
            INSERT INTO uploads (user_id, timestamp, data_type, data_content, metadata)
            VALUES (%s, %s, %s, %s, %s)
        """, (
            user_id,
            datetime.now(),
            data_type,
            data_content if isinstance(data_content, str) else data_content.decode('utf-8'),
            Json({"source": "streamlit_app"})
        ))

        conn.commit()
        cur.close()
        conn.close()
        return True
    except Exception as e:
        st.warning(f"Could not save to database: {e}")
        return False
```

## Deployment Options

### 1. Streamlit Cloud (Free, Recommended)

**Pros**: Free, easy, managed hosting
**Cons**: Public, limited resources

Steps:
1. Push code to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Sign in with GitHub
4. Click "New app"
5. Select repository: `your-username/nbragg`
6. Main file path: `streamlit_app.py`
7. Click "Deploy"

**Environment Variables** (for secrets):
Add to `.streamlit/secrets.toml`:
```toml
[email]
sender = "your-email@gmail.com"
password = "your-app-password"

[aws]
access_key = "YOUR_ACCESS_KEY"
secret_key = "YOUR_SECRET_KEY"
bucket = "nbragg-user-uploads"
```

Access in code:
```python
sender_email = st.secrets["email"]["sender"]
```

### 2. Heroku

**Pros**: Free tier available, supports background workers
**Cons**: Sleeps after 30 min inactivity

Files needed:

`runtime.txt`:
```
python-3.11.0
```

`Procfile`:
```
web: streamlit run streamlit_app.py --server.port=$PORT --server.address=0.0.0.0
```

`setup.sh`:
```bash
mkdir -p ~/.streamlit/
echo "\
[server]\n\
headless = true\n\
port = $PORT\n\
enableCORS = false\n\
\n\
" > ~/.streamlit/config.toml
```

Deploy:
```bash
heroku create your-app-name
git push heroku main
```

### 3. Docker + Cloud Run

`Dockerfile`:
```dockerfile
FROM python:3.11-slim

WORKDIR /app

# Copy requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy app
COPY . .

# Create data directory
RUN mkdir -p user_data_collection

EXPOSE 8080

CMD streamlit run streamlit_app.py \
    --server.port=8080 \
    --server.address=0.0.0.0 \
    --server.headless=true \
    --server.fileWatcherType=none
```

Build and deploy to Google Cloud Run:
```bash
gcloud builds submit --tag gcr.io/PROJECT_ID/nbragg-app
gcloud run deploy --image gcr.io/PROJECT_ID/nbragg-app --platform managed
```

### 4. Self-Hosted (VPS)

**Best for**: Full control, private data

Setup on Ubuntu server:
```bash
# Install dependencies
sudo apt update
sudo apt install python3.11 python3-pip nginx

# Clone repo
git clone https://github.com/your-username/nbragg.git
cd nbragg

# Install packages
pip install -r requirements.txt

# Run with supervisor (keeps app running)
sudo apt install supervisor

# Create supervisor config
sudo nano /etc/supervisor/conf.d/nbragg.conf
```

`/etc/supervisor/conf.d/nbragg.conf`:
```ini
[program:nbragg]
command=/usr/bin/python3 -m streamlit run streamlit_app.py --server.port=8501
directory=/home/user/nbragg
user=user
autostart=true
autorestart=true
stderr_logfile=/var/log/nbragg.err.log
stdout_logfile=/var/log/nbragg.out.log
```

Start:
```bash
sudo supervisorctl reread
sudo supervisorctl update
sudo supervisorctl start nbragg
```

Setup nginx reverse proxy:
```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://localhost:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
    }
}
```

## Security Considerations

### 1. User Data Privacy

- ✅ Use anonymous user IDs (already implemented)
- ✅ Don't collect IP addresses or personal information
- ✅ Explicit consent required (already implemented)
- ⚠️ Consider GDPR compliance if deploying in EU
- ⚠️ Add privacy policy link

### 2. File Upload Security

Add file validation to `streamlit_app.py`:

```python
def validate_csv_upload(uploaded_file):
    """Validate uploaded CSV file"""
    # Check file size (max 10MB)
    if uploaded_file.size > 10 * 1024 * 1024:
        return False, "File too large (max 10MB)"

    # Check file extension
    if not uploaded_file.name.endswith('.csv'):
        return False, "Only CSV files allowed"

    # Try to parse as CSV
    try:
        df = pd.read_csv(uploaded_file)
        if len(df) > 100000:  # Max 100k rows
            return False, "Too many rows (max 100k)"
        return True, "Valid"
    except Exception as e:
        return False, f"Invalid CSV: {e}"
```

### 3. Rate Limiting

Add basic rate limiting:

```python
from datetime import datetime, timedelta

if 'upload_history' not in st.session_state:
    st.session_state.upload_history = []

def check_rate_limit():
    """Allow max 10 uploads per hour"""
    now = datetime.now()
    # Remove old entries
    st.session_state.upload_history = [
        t for t in st.session_state.upload_history
        if now - t < timedelta(hours=1)
    ]

    if len(st.session_state.upload_history) >= 10:
        return False

    st.session_state.upload_history.append(now)
    return True
```

## Monitoring

### Check Data Collection

```bash
# Count uploaded files
ls user_data_collection/*.csv | wc -l

# View metadata
tail -f user_data_collection/metadata.jsonl | jq '.'

# Disk usage
du -sh user_data_collection/
```

### App Analytics

Add Google Analytics (optional):

```python
# Add to streamlit_app.py
import streamlit.components.v1 as components

def add_google_analytics():
    components.html("""
        <!-- Global site tag (gtag.js) - Google Analytics -->
        <script async src="https://www.googletagmanager.com/gtag/js?id=G-XXXXXXXXXX"></script>
        <script>
          window.dataLayer = window.dataLayer || [];
          function gtag(){dataLayer.push(arguments);}
          gtag('js', new Date());
          gtag('config', 'G-XXXXXXXXXX');
        </script>
    """, height=0)

# Call at app start
add_google_analytics()
```

## Troubleshooting

### App Crashes on Large Files

Increase memory limit:

`.streamlit/config.toml`:
```toml
[server]
maxUploadSize = 200
maxMessageSize = 200
```

### Slow Fitting

Add progress feedback:

```python
with st.spinner("Fitting model (this may take a minute)..."):
    progress_bar = st.progress(0)
    result = model.fit(data, wlmin=wlmin, wlmax=wlmax, progress_bar=False)
    progress_bar.progress(100)
```

### Data Collection Fills Disk

Add auto-cleanup:

```python
def cleanup_old_data(days=30):
    """Remove data older than N days"""
    cutoff = datetime.now() - timedelta(days=days)
    data_dir = Path("user_data_collection")

    for file in data_dir.glob("*.csv"):
        if datetime.fromtimestamp(file.stat().st_mtime) < cutoff:
            file.unlink()
```

## Support

For issues or questions:
- GitHub Issues: https://github.com/your-username/nbragg/issues
- Email: your-email@example.com
