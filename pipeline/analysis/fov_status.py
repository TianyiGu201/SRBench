"""Keep whole-FOV failures distinct from undefined individual gene metrics."""
import numpy as np
import pandas as pd

ACCURACY_METRICS = ('PCC', 'SCC', 'JSD', 'KLD')
SCORE_COLUMNS = ('PAR', 'PARnorm', 'PAR_all_fov', 'PARnorm_all_fov')

def mask_non_evaluable_fovs(frame, score_columns=SCORE_COLUMNS):
    result = frame.copy()
    missing = set(ACCURACY_METRICS) - set(result.columns)
    if missing:
        raise ValueError(f"FOV table is missing raw accuracy columns: {sorted(missing)}")
    keys = [c for c in ('dataset', 'dataset_id', 'method', 'FOV_size', 'FOV_region') if c in result.columns]
    if 'FOV_region' not in keys or 'FOV_size' not in keys:
        raise ValueError('FOV_size and FOV_region are required to identify prediction attempts')
    values = result[list(ACCURACY_METRICS)].apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan)
    any_finite = values.notna().any(axis=1)
    evaluable = any_finite.groupby([result[k] for k in keys], dropna=False).transform('any')
    for column in score_columns:
        if column in result:
            result.loc[~evaluable, column] = np.nan
    return result
