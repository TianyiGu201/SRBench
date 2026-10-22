"""Classify genes independently by expression level and spatial variability."""
from __future__ import annotations
import numpy as np
import pandas as pd

def spatial_variability_level(moran_i):
    if pd.isna(moran_i):
        return np.nan
    if moran_i < 0.2:
        return 'Low'
    if moran_i <= 0.5:
        return 'Medium'
    return 'High'

def classify_genes(annotation_csv, prediction_csv, output_csv):
    annotations = pd.read_csv(annotation_csv)
    predictions = pd.read_csv(prediction_csv)
    df = annotations.merge(predictions, on='gene', how='left', validate='one_to_one')
    df['Spatial_Variability_Level'] = df['morans_I'].map(spatial_variability_level)
    q25, q75 = df['gt_counts'].quantile([0.25, 0.75])
    def expression_level(count):
        if pd.isna(count):
            return np.nan
        if count < q25:
            return 'Low'
        if count <= q75:
            return 'Medium'
        return 'High'
    df['Gene_Expression_Level'] = df['gt_counts'].map(expression_level)
    df.to_csv(output_csv, index=False)
    return df
