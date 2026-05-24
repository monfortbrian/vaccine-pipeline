# TOPE_DEEP API - Dockerfile
# Local safety screening databases downloaded at build time.

FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc curl perl make wget tar libgomp1 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements_api_lean.txt ./
RUN pip install --no-cache-dir -r requirements_api_lean.txt

# MHCflurry models
RUN mhcflurry-downloads fetch \
    && echo "MHCflurry: OK" \
    || echo "WARNING: MHCflurry download failed"

# Copy source
COPY src/ src/
COPY data/ data/
COPY api/ api/
COPY run_pipeline.py ./

# Download safety screening databases
# AllergenOnline (FAO/WHO regulatory allergen database) + Human Swiss-Prot
# Baked into image - zero runtime dependency on external servers
RUN python data/safety_db/download_databases.py \
    && echo "Safety databases: OK" \
    || echo "WARNING: Safety database download failed - N6 will use AllerTOP-only mode"

# Verify IEDB coverage tool
RUN python3 -c "\
import sys; sys.path.insert(0, 'src/tools/population_coverage'); \
from population_calculation import PopulationCoverage; \
print('N7 IEDB tool: OK')" \
    || echo "WARNING: N7 IEDB tool not found - AFND 2020 fallback active"

EXPOSE 8000

CMD ["sh", "-c", "uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8000}"]
