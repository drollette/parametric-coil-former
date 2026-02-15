FROM condaforge/mambaforge:latest

WORKDIR /app

# Install conda dependencies (CadQuery with OpenCASCADE)
COPY environment.yml .
RUN mamba env update -n base -f environment.yml && mamba clean -afy

# Install Python dependencies (FastAPI, uvicorn)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY backend/ ./backend/
COPY wasm-coil-former/static/ ./wasm-coil-former/static/

# Create outputs directory
RUN mkdir -p outputs

# Expose port
EXPOSE 8000

# Run the application
CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]
