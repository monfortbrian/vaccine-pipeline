# TOPE_DEEP API - Dockerfile
# Python 3.11-slim base
#
# Tools bundled at build time (no runtime downloads needed):
#   MHCflurry 2.0 models       - mhcflurry-downloads fetch (~200MB)
#   IEDB Population Coverage   - bundled in src/tools/population_coverage/
#
# PSORTb note:
#   PSORTb v3.0 requires Docker-in-Docker which Railway does not support.
#   Phobius (phobius.sbc.su.se) is used instead for transmembrane and
#   signal peptide prediction - scientifically equivalent for surface
#   antigen identification.
#   Ref: Kall et al., J Mol Biol 2004; Gardy et al., Bioinformatics 2003.

FROM python:3.11-slim

WORKDIR /app

# System dependencies
# libgomp1: required by MHCflurry (OpenMP threading in numpy/sklearn)
# curl:     used by healthcheck in docker-compose and Railway
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc curl perl make wget tar libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies
COPY requirements_api_lean.txt ./
RUN pip install --no-cache-dir -r requirements_api_lean.txt

# MHCflurry models - downloaded at build time, cached inside the image.
# Avoids runtime dependency on external servers.
# Stored at: /root/.local/share/mhcflurry/
RUN mhcflurry-downloads fetch \
    && echo "MHCflurry models: OK" \
    && ls /root/.local/share/mhcflurry/ \
    || echo "WARNING: MHCflurry model directory not found - N3 MHCflurry fallback will be unavailable"

# Copy application source
COPY src/ src/
COPY data/ data/
COPY api/ api/
COPY run_pipeline.py ./

# Verify IEDB Population Coverage Tool is present.
# Tool lives at src/tools/population_coverage/population_calculation.py.
# If missing, N7 silently uses AFND 2020 static frequency fallback - still valid,
# but IEDB tool v3.0.1 is preferred.
RUN python3 -c "\
import sys; sys.path.insert(0, 'src/tools/population_coverage'); \
from population_calculation import PopulationCoverage; \
print('N7 IEDB coverage tool: OK')" \
    || echo "WARNING: N7 IEDB tool not found - will use AFND 2020 fallback"

# Port 8000 - must match CMD below and docker-compose ports mapping
EXPOSE 8000

# PORT env var is set by Railway at runtime.
# Locally (docker compose) it falls back to 8000.
# --reload is omitted here - hot reload is handled by volume mounts in docker-compose.
CMD ["sh", "-c", "uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8000}"]