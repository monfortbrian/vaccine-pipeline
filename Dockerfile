# Stable Python version
FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc curl perl make wget tar \
    && rm -rf /var/lib/apt/lists/*


# Install Python deps
COPY requirements_api_lean.txt ./
RUN pip install --no-cache-dir -r requirements_api_lean.txt

# Copy source
COPY src/ src/
COPY data/ data/
COPY api/ api/
COPY run_pipeline.py ./

# Expose port
EXPOSE 8000

# Start FastAPI
CMD ["sh", "-c", "uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8000}"]