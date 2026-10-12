# Study Design

A comprehensive overview of the platform's features, evaluated datasets, and visualization tools for benchmarking spatial transcriptomics super-resolution methods.

---

## 🔄 Overview

This platform provides interactive visualization of benchmarking results for computational super-resolution methods applied to spatially resolved transcriptomics (SRT). Rather than presenting a complete study report, this site focuses on making key results explorable and interpretable — offering additional context and visual tools to complement the published findings.

The core benchmarking strategy involves collecting high-resolution SRT datasets, intentionally lowering their resolution to simulate low-resolution input, and using the original high-resolution data as ground truth to evaluate prediction performance. In our study, we evaluated accuracy, scalability and usability of each super-resolution method. Here, in this platform, we focus on the visualization of the prediction accuracy.

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐     ┌───────────────┐
│  High-Res SRT   │     │  LR Simulation  │     │  SR Prediction  │     │   Evaluation  │
├─────────────────┤     ├─────────────────┤     ├─────────────────┤     ├───────────────┤
│   Xenium        │     │ Spot Layout     │     │    XFuse        │     │  Accuracy     │
│   CosMx SMI     │ --> │ Size & Distance │ --> │    TESLA        │ --> │  Scalability  │
│   Visium HD     │     │ Downsampling    │     │    iSTAR        │     │  Usability    │
└─────────────────┘     └─────────────────┘     │    scstGCN      │     └───────────────┘
                                                │    SuperST_t    │    
                                                │    SuperST_t_g  │
                                                └─────────────────┘
```

---

## 🗄️ Evaluated Datasets

Our study utilizes high-resolution spatial transcriptomics datasets from three leading technologies. Detailed metadata for each dataset can be explored in the [Dataset Overview](/dataset-overview).

| Technology | Data Format | Analysis Scope |
|------------|-------------|----------------|
| **Xenium** | Transcript-level coordinates | All genes in panel |
| **CosMx SMI** | Transcript-level coordinates | All genes in panel |
| **Visium HD** | High-resolution bin-like matrix | Spatially variable genes (SVGs) only |

Since different platforms store data in fundamentally different formats (coordinate-based vs. matrix-based), and different super-resolution algorithms require different input formats, we developed a standardized preprocessing pipeline to:
1. Simulate low-resolution SRT data from high-resolution input across various experimental scenarios
2. Generate corresponding ground truth data for consistent evaluation

---

## 🔬 Metric Visualization

To help researchers intuitively understand what different metric values mean in practice, we provide an interactive visualization tool. You can specify a target value for a selected metric, and the platform will display — for a chosen dataset — the H&E image, the ground truth gene expression, and the spatial expression maps of the 3 genes whose prediction metric values are closest to the specified target.

Because gene expression levels vary substantially across genes, two color scale modes are available: an automatic mode that adjusts the scale across all displayed genes simultaneously, and a per-gene mode that allows independent color range specification for each gene.

➡️ [Explore Metric Visualization](/metric-vis)

---

## 🔲 Out-of-Sample Prediction

In this scenario, a super-resolution method is trained using a single selected field of view (FOV) that contains both the H&E image and the corresponding spatial transcriptomics (SRT) gene expression data within that region. After training, the model is applied to predict gene expression across the entire tissue section.

To evaluate how the spatial location of the training region affects prediction performance, each FOV in the dataset is used once as the training region. For each run, predictions across the tissue are compared with the ground-truth SRT data, and the resulting performance is recorded at the center location of the corresponding FOV. This process generates a spatial performance map that reflects how prediction quality varies depending on the location of the training region.

The platform visualizes these results by overlaying a grid of FOVs (each with a side length of 1,500 pixels) on the tissue image. Each square represents one FOV used for training, and its color indicates the prediction performance achieved when that FOV served as the training region.

Hovering over a square highlights the spatial extent of the corresponding FOV on the tissue. Clicking a square reveals detailed per-gene evaluation metrics (e.g., PCC, SCC, JSD, and KLD) for that prediction.

➡️ [View Results](/out-of-sample)

---

## 📚 Initial Resolution Impact

We investigate how the initial resolution of the low-resolution input affects super-resolution performance. Four common spot layout configurations are simulated, each with varying spot diameter or center-to-center distance, to reflect different SRT technology configurations:

| Layout | Description |
|--------|-------------|
| **Proportional** | Spot diameter to center-to-center distance ratio is fixed at 55/100 across all resolutions, consistent with Visium's design. |
| **Dense** | Spots are packed densely such that each spot's diameter equals the center-to-center distance. |
| **Fixed Spot Distance** | Center-to-center distance varies across resolutions, but spot diameter is fixed at 100 µm. |
| **Fixed Spot Size** | Spot diameter varies across resolutions, but center-to-center distance is fixed at 55 µm. |


You can select a super-resolution method and a dataset to visualize how performance changes across different initial spatial resolutions. Results are displayed for a chosen statistical metric, with values representing the mean across 3 replicate predictions. The per-gene metric values for each replicate are also available for detailed inspection.

➡️ [View Results](/diff-resolution)

---

## 📊 Evaluation Metrics

We use a suite of statistical metrics to quantify the similarity between predicted gene expression maps and the ground truth. You can interactively explore these in the [Metric Visualization](/metric-vis) module.

### Metrics Overview

| Metric | Description |
|--------|-------------|
| **PCC** (Pearson Correlation Coefficient) | Measures the linear correlation between predicted and actual expression levels. Range: [−1, 1]. |
| **SCC** (Spearman Correlation Coefficient) | Assesses monotonic relationships, robust to outliers and non-linear associations. Range: [−1, 1]. |
| **JSD** (Jensen–Shannon Divergence) | A symmetric and smoothed version of KLD, providing a stable distance metric between distributions. Range: [0, 1]. |
| **KLD** (Kullback–Leibler Divergence) | Measures information loss when the predicted distribution approximates the true distribution. Range: [0, +∞). |
| **PARnorm** (Normalized Prediction accuracy rank) | Normalized average rank across PCC, SCC, JSD, and KLD; lower values indicate better overall performance. Range: [0, 1]. |
| **SSIM** (Structural Similarity Index) | Evaluates the preservation of spatial patterns and structures in the expression map. Range: [0, 1]. |
| **PSNR** (Peak Signal-to-Noise Ratio) | Quantifies reconstruction quality; higher values indicate lower error relative to signal strength. Range: [0, +∞). |


