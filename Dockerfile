# Stable Python version
FROM python:3.11-slim

WORKDIR /app

# System deps needed for compilation (gcc) or curl if needed
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc curl && rm -rf /var/lib/apt/lists/*

# Only install API dependencies (lean)
COPY requirements_api_lean.txt ./
RUN pip install --no-cache-dir -r requirements_api_lean.txt

# Copy source code
COPY src/ src/
COPY data/ data/
COPY api/ api/
COPY run_pipeline.py ./

# Expose port (Railway uses $PORT)
EXPOSE 8000

# Start FastAPI using dynamic Railway port
CMD ["sh", "-c", "uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8000}"]