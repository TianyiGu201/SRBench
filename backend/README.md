# SpatialVis Backend

Python FastAPI backend for the Spatial Transcriptomics Visualization Platform.

## Requirements

- Python 3.9+
- Dependencies listed in `requirements.txt`

## Setup

1. **Create virtual environment:**
```bash
python -m venv venv
```

2. **Activate virtual environment:**

Windows:
```bash
venv\Scripts\activate
```

Linux/Mac:
```bash
source venv/bin/activate
```

3. **Install dependencies:**
```bash
pip install -r requirements.txt
```

4. **Configure environment variables:**

Edit `.env` file to set your data root directory:
```
DATA_ROOT=D:/path/to/your/data_root
PORT=8000
HOST=0.0.0.0
```

## Running the Server

**Development mode:**
```bash
python main.py
```

Or using uvicorn directly:
```bash
uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

The API will be available at `http://localhost:8000`

API documentation: `http://localhost:8000/docs`

## Data Directory Structure

The backend expects data to be organized as follows:

```
data_root/
├── dataset_A/
│   ├── dataset_metaINFO.json
│   ├── metric_vis.json
│   ├── he.tiff
│   ├── he.png
│   ├── reference.h5ad
│   └── prediction/
│       ├── iSTAR/
│       │   ├── out_of_sample_prediction/
│       │   │   ├── 10_20_30_40.csv
│       │   │   └── ...
│       │   └── different_initial_resolution/
│       │       ├── dense/
│       │       │   ├── r_spotSize_55.h5ad
│       │       │   ├── r_spotSize_55.csv
│       │       │   └── ...
│       │       ├── proportional/
│       │       ├── fixed_spotSize/
│       │       └── fixed_spotDistance/
│       ├── XFuse/
│       ├── TESLA/
│       ├── scstGCN/
│       └── SuperST/
└── dataset_B/
    └── ...
```

## API Endpoints

### Dataset Overview
- `GET /api/datasets` - List all datasets
- `GET /api/datasets/{dataset_name}/metadata` - Get dataset metadata
- `GET /api/datasets/{dataset_name}/he-thumbnail` - Get thumbnail image
- `GET /api/datasets/{dataset_name}/he-image` - Download full HE image

### Metric Visualization
- `GET /api/metric-vis/{dataset_name}/options` - Get metric options
- `GET /api/metric-vis/{dataset_name}/{metric}/{value}/genes` - Get similar genes
- `GET /api/spatial-data/{dataset_name}/reference/{gene}` - Get reference spatial data
- `GET /api/spatial-data/{dataset_name}/prediction/{gene}` - Get prediction spatial data

### Out-of-sample Prediction
- `GET /api/out-of-sample/{dataset_name}/{algorithm}/regions` - Get all regions
- `GET /api/out-of-sample/{dataset_name}/{algorithm}/region/{region_id}` - Get region data
- `GET /api/out-of-sample/{dataset_name}/{algorithm}/region/{region_id}/download` - Download region CSV

### Different Initial Resolution
- `GET /api/diff-res/{dataset_name}/{algorithm}/arrangements` - Get spot arrangements
- `GET /api/diff-res/{dataset_name}/{algorithm}/{arrangement}/resolutions` - Get resolutions
- `GET /api/diff-res/{dataset_name}/{algorithm}/{arrangement}/{resolution}/genes` - Get gene list
- `GET /api/diff-res/{dataset_name}/{algorithm}/{arrangement}/{resolution}/gene/{gene}` - Get gene data

### Health Check
- `GET /api/health` - Server health status
- `GET /` - API info

## File Format Notes

**h5ad files:**
- Must contain `obsm['spatial']` - N×2 spatial coordinates
- Must contain `X` - Expression matrix
- Must contain `var_names` - Gene names

**CSV files:**
- Gene metrics must include columns: gene, PCC, SCC, KLD, JSD, SSIM, PSNR
- Out-of-sample prediction CSVs should also include PARnorm column

**JSON files:**
- `dataset_metaINFO.json` - Dataset metadata (name, species, technology, tissue, description, original_data_link)
- `metric_vis.json` - Metric visualization lookup (nested dict: metric -> value -> [genes])
