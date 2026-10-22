"""Plot gene metrics stratified by expression level and spatial variability."""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle


def plot_3x3_metric_heatmap(method: str, dataset: str, metric: str, merged_csv_dir: str='PATH/TO/REQUIRED_INPUT', output_dir: str='PATH/TO/REQUIRED_INPUT', *, cmap_high_better: str='RdBu_r', cmap_low_better: str='RdBu', center: str='median', vmin=None, vmax=None, order=('Low', 'Medium', 'High'), x_tick_rotation: float=90.0, remove_tick_marks: bool=True, grid_linewidth: float=0.45, spine_linewidth: float=0.75, cbar_width_mm: float=3.0, cbar_gap_mm: float=1.5, annot: bool=True, annot_decimals: int=3, auto_text_contrast: bool=True, contrast_threshold: float=0.45, missing_color: str='#F7F7F7', missing_label: str='n.d.', missing_hatch: str='///', missing_hatch_color: str='#C7C7C7', show_colorbar: bool=True, show_tick_labels: bool=True, show_axis_titles: bool=True, title_top: bool=True, hvg_status: str=None, save_pdf: bool=True, show: bool=False, dpi: int=300):
    """Plot expression-by-spatial-variability metric heatmaps with explicit margins."""
    os.makedirs(output_dir, exist_ok=True)
    MM_PER_INCH = 25.4
    FINAL_WIDTH_MM = 58.0
    ORIGINAL_WIDTH_MM = 58.0
    SCALE = FINAL_WIDTH_MM / ORIGINAL_WIDTH_MM
    TARGET_MAIN_PT = 7.0
    TARGET_CBAR_PT = 5.0
    MAIN_FONTSIZE = TARGET_MAIN_PT / SCALE
    CBAR_FONTSIZE = TARGET_CBAR_PT / SCALE
    ORIGINAL_HEIGHT_MM = ORIGINAL_WIDTH_MM
    FIG_WIDTH_IN = ORIGINAL_WIDTH_MM / MM_PER_INCH
    FIG_HEIGHT_IN = ORIGINAL_HEIGHT_MM / MM_PER_INCH
    figsize = (FIG_WIDTH_IN, FIG_HEIGHT_IN)
    left_margin_mm = 14.5
    right_margin_mm = cbar_width_mm + cbar_gap_mm + 2.5
    bottom_margin_mm = 6
    top_margin_mm = 0
    left_frac = left_margin_mm / ORIGINAL_WIDTH_MM
    right_frac = 1.0 - right_margin_mm / ORIGINAL_WIDTH_MM
    bottom_frac = bottom_margin_mm / ORIGINAL_HEIGHT_MM
    top_frac = 1.0 - top_margin_mm / ORIGINAL_HEIGHT_MM
    cbar_width_frac = cbar_width_mm / ORIGINAL_WIDTH_MM
    cbar_gap_frac = cbar_gap_mm / ORIGINAL_WIDTH_MM
    csv_path = os.path.join(merged_csv_dir, f'{dataset}_{method}_classified.csv')
    df = pd.read_csv(csv_path)
    required_cols = {'Spatial_Variability_Level', 'Gene_Expression_Level', metric}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f'CSV is missing columns: {missing_cols}')
    if hvg_status is not None:
        if 'hvg_status' not in df.columns:
            raise ValueError('An hvg_status filter was supplied, but the CSV has no hvg_status column.')
        df = df[df['hvg_status'] == hvg_status].copy()
    order = list(order)
    df = df[df['Spatial_Variability_Level'].isin(order) & df['Gene_Expression_Level'].isin(order)].copy()
    mat = df.pivot_table(index='Spatial_Variability_Level', columns='Gene_Expression_Level', values=metric, aggfunc='mean').reindex(index=order, columns=order)
    data = mat.values.astype(float)
    lower_is_better = metric.upper() in {'JSD', 'KLD'}
    cmap_name = cmap_low_better if lower_is_better else cmap_high_better
    finite_vals = data[np.isfinite(data)]
    if finite_vals.size == 0:
        raise ValueError('No finite heatmap values; check the data and filtering criteria.')
    _vmin = float(np.min(finite_vals)) if vmin is None else vmin
    _vmax = float(np.max(finite_vals)) if vmax is None else vmax
    norm = Normalize(vmin=_vmin, vmax=_vmax)
    data_ma = np.ma.masked_invalid(data)
    cmap_obj = plt.get_cmap(cmap_name)
    try:
        cmap_obj = cmap_obj.copy()
    except Exception:
        pass
    cmap_obj.set_bad(color=missing_color)

    def _luminance(rgb):
        r, g, b = rgb
        return 0.2126 * r + 0.7152 * g + 0.0722 * b

    def _format_no_negzero(x, decimals, eps=None):
        """
        Format number and avoid displaying '-0' caused by rounding.
        eps: values with abs(x) < eps will be treated as 0.
        """
        if eps is None:
            eps = 0.5 * 10 ** (-decimals)
        if np.isfinite(x) and abs(float(x)) < eps:
            x = 0.0
        s = f'{float(x):.{decimals}f}'.rstrip('0').rstrip('.')
        return '0' if s == '-0' else s
    fig = plt.figure(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor('white')
    heatmap_width = right_frac - left_frac - cbar_width_frac - cbar_gap_frac
    heatmap_height = top_frac - bottom_frac
    if heatmap_width < heatmap_height:
        ax_width = heatmap_width
        ax_height = heatmap_width
        ax_bottom = bottom_frac + (heatmap_height - ax_height) / 2
        ax_left = left_frac
    else:
        ax_height = heatmap_height
        ax_width = heatmap_height
        ax_left = left_frac + (heatmap_width - ax_width) / 2
        ax_bottom = bottom_frac
    ax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])
    x = np.arange(4)
    y = np.arange(4)
    quad = ax.pcolormesh(x, y, data_ma, cmap=cmap_obj, norm=norm, shading='flat', antialiased=False, linewidth=0.0, edgecolors='none')
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 3)
    ax.set_aspect('equal', adjustable='box')
    for k in (1, 2):
        ax.axvline(k, color='black', linewidth=grid_linewidth)
        ax.axhline(k, color='black', linewidth=grid_linewidth)
    for s in ax.spines.values():
        s.set_linewidth(spine_linewidth)
        s.set_edgecolor('black')
    ax.set_xticks([0.5, 1.5, 2.5])
    ax.set_yticks([0.5, 1.5, 2.5])
    if show_tick_labels:
        ax.set_xticklabels(order, fontsize=MAIN_FONTSIZE)
        ax.set_yticklabels(order, fontsize=MAIN_FONTSIZE)
        plt.setp(ax.get_xticklabels(), rotation=x_tick_rotation, ha='right', va='center', rotation_mode='anchor')
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    if remove_tick_marks:
        ax.tick_params(axis='both', which='both', length=0, width=0, pad=1.0)
    else:
        ax.tick_params(axis='both', which='major', length=2.0, width=0.6, pad=1.0)
    if show_axis_titles:
        ax.set_xlabel('Gene expression level', fontsize=MAIN_FONTSIZE)
        ax.set_ylabel('Spatial variability level', fontsize=MAIN_FONTSIZE)
    if annot:
        for r in range(3):
            for c in range(3):
                v = data[r, c]
                cx, cy = (c + 0.5, r + 0.5)
                if np.isfinite(v):
                    txt = _format_no_negzero(v, annot_decimals)
                    if auto_text_contrast:
                        rgba = cmap_obj(norm(v))
                        fg = 'white' if _luminance(rgba[:3]) < contrast_threshold else 'black'
                    else:
                        fg = 'black'
                    ax.text(cx, cy, txt, ha='center', va='center', fontsize=MAIN_FONTSIZE, color=fg)
                else:
                    ax.text(cx, cy, missing_label, ha='center', va='center', fontsize=MAIN_FONTSIZE, color='#5A5A5A', fontweight='bold')
                    if missing_hatch:
                        ax.add_patch(Rectangle((c, r), 1.0, 1.0, facecolor='none', edgecolor=missing_hatch_color, hatch=missing_hatch, linewidth=0.0))
    if show_colorbar:
        from matplotlib.ticker import MaxNLocator
        cax_left = ax_left + ax_width + cbar_gap_frac
        cax_bottom = ax_bottom
        cax_width = cbar_width_frac
        cax_height = ax_height
        cax = fig.add_axes([cax_left, cax_bottom, cax_width, cax_height])
        cbar = fig.colorbar(quad, cax=cax)
        direction = 'lower better' if lower_is_better else 'higher better'
        cbar.set_label(f'{metric} ({direction})', fontsize=CBAR_FONTSIZE, labelpad=2)
        vmin_actual = cbar.vmin
        vmax_actual = cbar.vmax
        data_range = vmax_actual - vmin_actual
        if data_range < 0.1:
            n_ticks = 3
        elif data_range < 1.0:
            n_ticks = 4
        elif data_range < 5.0:
            n_ticks = 5
        else:
            n_ticks = 6
        tick_positions = np.linspace(vmin_actual, vmax_actual, n_ticks)
        tick_labels = [_format_no_negzero(val, 2) for val in tick_positions]
        if len(set(tick_labels)) < len(tick_labels):
            n_ticks = max(3, n_ticks - 1)
            tick_positions = np.linspace(vmin_actual, vmax_actual, n_ticks)
            tick_labels = [_format_no_negzero(val, 2) for val in tick_positions]
        cbar.set_ticks(tick_positions)
        cbar.set_ticklabels(tick_labels)
        cbar.ax.tick_params(labelsize=CBAR_FONTSIZE, length=3.0, width=0.6, direction='out', pad=1.0)
        for s in cax.spines.values():
            s.set_linewidth(0.6)
            s.set_edgecolor('black')
    if title_top:
        title_x = ax_left + ax_width / 2
        title_y = ax_bottom + ax_height
        fig.text(title_x, title_y, f'Method:{method}   Metric:{metric}', ha='center', va='bottom', fontsize=MAIN_FONTSIZE)
    out_name = f'{dataset}_{method}_{metric}_heatmap'
    if hvg_status is not None:
        out_name += f'_{hvg_status}'
    out_path = os.path.join(output_dir, out_name + '.pdf')
    if save_pdf:
        fig.savefig(out_path, dpi=dpi, bbox_inches=None, pad_inches=0)
    if show:
        plt.show()
    else:
        plt.close(fig)
    print(f'Heatmap saved to: {out_path}')
    return mat

