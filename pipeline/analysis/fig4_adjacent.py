"""Compare within-slide and adjacent-slide gene predictions."""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pickle
import gzip
from scipy.stats import gaussian_kde, wilcoxon
from PIL import Image

pkl_rngs = {'within_slice': {'pd_row_rng': (104, 532), 'pd_col_rng': (104, 987)}, 'cross_slice': {'pd_row_rng': (532, 960), 'pd_col_rng': (104, 987)}}
need_pad = ['Xe_hsa_heart_h', 'Xe_hsa_liver_h', 'Xe_hsa_skin_h', 'Xe_hsa_kidney_h']

class NumpyCompatUnpickler(pickle.Unpickler):

    def find_class(self, module, name):
        if module.startswith('numpy._core'):
            module = module.replace('numpy._core', 'numpy.core')
        return super().find_class(module, name)

def load_exp_matrix(file_path: str):
    if file_path.endswith('.pkl.gz'):
        with gzip.open(file_path, 'rb') as f:
            return NumpyCompatUnpickler(f).load()
    else:
        with open(file_path, 'rb') as f:
            return NumpyCompatUnpickler(f).load()

def load_he_hw(he_path: str):
    with Image.open(he_path) as im:
        w, h = im.size
    return (h, w)

def build_gt_hist_by_he_grid(coords_yx: np.ndarray, pd_h: int, pd_w: int, he_h: int, he_w: int):
    """
    coords_yx: (N,2) float, columns [y, x] in HE pixel coordinate system
    returns hist: (pd_h, pd_w) counts (zeros where uncovered)
    """
    hist = np.zeros((pd_h, pd_w), dtype=np.float32)
    if coords_yx is None or coords_yx.size == 0:
        return hist
    patch_h = he_h / float(pd_h)
    patch_w = he_w / float(pd_w)
    y = coords_yx[:, 0]
    x = coords_yx[:, 1]
    yi = np.floor(y / patch_h).astype(np.int64)
    xi = np.floor(x / patch_w).astype(np.int64)
    yi = np.clip(yi, 0, pd_h - 1)
    xi = np.clip(xi, 0, pd_w - 1)
    np.add.at(hist, (yi, xi), 1.0)
    return hist

def plot_spatial_gene_prediction(data_path_root: str, method: str, simulation: str, gene_name: str, vmin=0, vmax=20, use_max=False, cmap='RdBu_r', cbar_position='right', fig_name=None):
    """
    Plot spatial gene expression prediction with fixed physical size.

    Illustrator target (after scaling to 28 mm width):
    - Title font ≈ 6 pt
    - Colorbar tick font ≈ 4 pt
    """
    MM_PER_INCH = 25.4
    FIG_WIDTH_MM = 140
    FINAL_WIDTH_MM = 28
    SCALE = FINAL_WIDTH_MM / FIG_WIDTH_MM
    LEFT_MARGIN_MM = 6
    RIGHT_MARGIN_MM = 6
    TOP_MARGIN_FRAC = 0.1
    CBAR_WIDTH_MM_VERTICAL = 6
    CBAR_HEIGHT_MM_HORIZONTAL = 6
    CBAR_GAP_MM = 2
    TARGET_TITLE_PT = 7
    TARGET_CBAR_TICK_PT = 5
    TITLE_FONTSIZE = TARGET_TITLE_PT / SCALE
    CBAR_TICK_FONTSIZE = TARGET_CBAR_TICK_PT / SCALE
    gene_pred = load_exp_matrix(f'{data_path_root}/{method}/{simulation}/{gene_name}.pkl.gz')
    if use_max:
        vmin = float(np.nanmin(gene_pred))
        vmax = float(np.nanmax(gene_pred))
    norm = Normalize(vmin=vmin, vmax=vmax)
    img_h, img_w = gene_pred.shape
    img_aspect = img_w / img_h
    if cbar_position == 'right':
        available_width_mm = FIG_WIDTH_MM - LEFT_MARGIN_MM - RIGHT_MARGIN_MM - CBAR_WIDTH_MM_VERTICAL - CBAR_GAP_MM
    else:
        available_width_mm = FIG_WIDTH_MM - LEFT_MARGIN_MM - RIGHT_MARGIN_MM
    panel_width_mm = available_width_mm
    panel_width_in = panel_width_mm / MM_PER_INCH
    panel_height_in = panel_width_in / img_aspect
    panel_height_mm = panel_height_in * MM_PER_INCH
    if cbar_position == 'bottom':
        FIG_HEIGHT_MM = panel_height_mm * (1 + TOP_MARGIN_FRAC) + CBAR_HEIGHT_MM_HORIZONTAL + CBAR_GAP_MM
    else:
        FIG_HEIGHT_MM = panel_height_mm * (1 + TOP_MARGIN_FRAC)
    FIG_WIDTH_IN = FIG_WIDTH_MM / MM_PER_INCH
    FIG_HEIGHT_IN = FIG_HEIGHT_MM / MM_PER_INCH
    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
    left_frac = LEFT_MARGIN_MM / FIG_WIDTH_MM
    width_frac = panel_width_mm / FIG_WIDTH_MM
    height_frac = panel_height_mm / FIG_HEIGHT_MM
    bottom_frac = 1 - height_frac - TOP_MARGIN_FRAC
    ax = fig.add_axes([left_frac, bottom_frac, width_frac, height_frac])
    ax.imshow(gene_pred, cmap=cmap, norm=norm, origin='upper')
    ax.set_axis_off()
    ax.set_title(f'{method} prediction', fontsize=TITLE_FONTSIZE, pad=10)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    if cbar_position == 'right':
        cbar_left = (LEFT_MARGIN_MM + panel_width_mm + CBAR_GAP_MM) / FIG_WIDTH_MM
        cbar_width = CBAR_WIDTH_MM_VERTICAL / FIG_WIDTH_MM
        cax = fig.add_axes([cbar_left, bottom_frac, cbar_width, height_frac])
        cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
        ticks = np.linspace(vmin, vmax, 6)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f'{t:.2f}' for t in ticks])
        cbar.ax.tick_params(labelsize=CBAR_TICK_FONTSIZE, rotation=0)
        for label in cbar.ax.get_yticklabels():
            label.set_ha('left')
            label.set_va('center')
    else:
        cbar_height = CBAR_HEIGHT_MM_HORIZONTAL / FIG_HEIGHT_MM
        cbar_bottom = bottom_frac - CBAR_GAP_MM / FIG_HEIGHT_MM - cbar_height
        cax = fig.add_axes([left_frac, cbar_bottom, width_frac, cbar_height])
        cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
        ticks = np.linspace(vmin, vmax, 6)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f'{t:.2f}' for t in ticks])
        cbar.ax.tick_params(labelsize=CBAR_TICK_FONTSIZE, rotation=90)
        for label in cbar.ax.get_xticklabels():
            label.set_ha('center')
            label.set_va('top')
    if fig_name:
        plt.savefig(fig_name, format='pdf', dpi=300)
    plt.show()
    return gene_pred.shape

def dump_pkl_gz(obj, path, compresslevel=6):
    with gzip.open(path, 'wb', compresslevel=compresslevel) as f:
        pickle.dump(obj, f)

def pad_pd_by_ratio(pd, ratio, pad_value=0):
    a, b = pd.shape
    top = int(round(ratio['top'] * a))
    bottom = int(round(ratio['bottom'] * a))
    left = int(round(ratio['left'] * b))
    right = int(round(ratio['right'] * b))
    return np.pad(pd, pad_width=((top, bottom), (left, right)), mode='constant', constant_values=pad_value)

def detect_black_border_ratio(he, eps=0):
    h, w, _ = he.shape
    is_black = np.all(he <= eps, axis=2)
    top = 0
    while top < h and np.all(is_black[top, :]):
        top += 1
    bottom = 0
    while bottom < h and np.all(is_black[h - 1 - bottom, :]):
        bottom += 1
    left = 0
    while left < w and np.all(is_black[:, left]):
        left += 1
    right = 0
    while right < w and np.all(is_black[:, w - 1 - right]):
        right += 1
    return {'top': top / h, 'bottom': bottom / h, 'left': left / w, 'right': right / w}

def plot_slice_comparison(method, data_dir, output_dir, box_width=0.5, violin_alpha=0.15, box_fill_alpha=0.4, outlier=False, outlier_size=10, outlier_alpha=0.6, rasterize_outliers=False, dpi=300, show=False, subplot_spacing=0.012, left_margin=0.05, right_margin=0.01, spine_linewidth=2.5, box_linewidth=2.0):
    """Plot six paired within-slide and adjacent-slide metric comparisons."""
    os.makedirs(output_dir, exist_ok=True)
    MM_PER_INCH = 25.4
    FINAL_WIDTH_MM = 85.0
    ORIGINAL_WIDTH_MM = 85.0
    SCALE = FINAL_WIDTH_MM / ORIGINAL_WIDTH_MM
    TARGET_TITLE_PT = 7.0
    TARGET_YLABEL_PT = 5.0
    TARGET_XLABEL_PT = 5.0
    TARGET_YTICK_PT = 5.0
    TITLE_FONTSIZE = TARGET_TITLE_PT / SCALE
    YLABEL_FONTSIZE = TARGET_YLABEL_PT / SCALE
    XLABEL_FONTSIZE = TARGET_XLABEL_PT / SCALE
    YTICK_FONTSIZE = TARGET_YTICK_PT / SCALE
    ORIGINAL_HEIGHT_MM = 30
    FIG_WIDTH_IN = ORIGINAL_WIDTH_MM / MM_PER_INCH
    FIG_HEIGHT_IN = ORIGINAL_HEIGHT_MM / MM_PER_INCH
    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), dpi=dpi)
    within_slice_path = os.path.join(data_dir, method, 'within_slice.csv')
    cross_slice_path = os.path.join(data_dir, method, 'cross_slice.csv')
    df_within = pd.read_csv(within_slice_path)
    df_cross = pd.read_csv(cross_slice_path)
    metrics = ['PCC', 'SCC', 'JSD', 'KLD', 'SSIM', 'PSNR']
    color_within = '#3b85bd'
    color_cross = '#9db648'
    n_subplots = 6
    available_width = 1.0 - left_margin - right_margin
    subplot_width = (available_width - (n_subplots - 1) * subplot_spacing) / n_subplots
    axes = []
    for i in range(n_subplots):
        left = left_margin + i * (subplot_width + subplot_spacing)
        bottom = 0.02
        width = subplot_width
        height = 0.84
        ax = fig.add_axes([left, bottom, width, height])
        axes.append(ax)
    subplot_area_left = left_margin
    subplot_area_right = left_margin + n_subplots * subplot_width + (n_subplots - 1) * subplot_spacing
    subplot_area_center = (subplot_area_left + subplot_area_right) / 2
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        paired = df_within[['gene', metric]].merge(
            df_cross[['gene', metric]], on='gene', suffixes=('_within', '_cross'), validate='one_to_one'
        )
        paired = paired.replace([np.inf, -np.inf], np.nan).dropna(subset=[metric + '_within', metric + '_cross'])
        within_vals = paired[metric + '_within']
        cross_vals = paired[metric + '_cross']
        if len(within_vals) > 0 and len(cross_vals) > 0:
            statistic, p_value = wilcoxon(within_vals, cross_vals, alternative='two-sided')
        else:
            p_value = 1.0
        sig_label_fontsize = XLABEL_FONTSIZE * 1.0
        if p_value < 0.001:
            sig_label = '***'
        elif p_value < 0.01:
            sig_label = '**'
        elif p_value < 0.05:
            sig_label = '*'
        else:
            sig_label = 'ns'
            sig_label_fontsize = sig_label_fontsize * 0.83
        x_positions = [0, 1]
        groups = [within_vals, cross_vals]
        colors = [color_within, color_cross]
        labels = ['Within Slide', 'Across Slide']
        for i, (vals, color, label) in enumerate(zip(groups, colors, labels)):
            if vals.empty:
                continue
            if outlier and len(vals) > 0:
                q1 = np.percentile(vals, 25)
                q3 = np.percentile(vals, 75)
                iqr = q3 - q1
                lower_bound = q1 - 1.5 * iqr
                upper_bound = q3 + 1.5 * iqr
                outliers = vals[(vals < lower_bound) | (vals > upper_bound)]
            if len(vals) > 1 and np.isfinite(vals).all() and np.ptp(vals.to_numpy()) > 0:
                kde = gaussian_kde(vals)
                y = np.linspace(vals.min(), vals.max(), 400)
                d = kde(y)
                d = d / d.max() * (box_width / 2)
                ax.fill_betweenx(y, x_positions[i] - d, x_positions[i] + d, facecolor=color, alpha=violin_alpha, linewidth=0, zorder=1)
            bp = ax.boxplot(vals, positions=[x_positions[i]], widths=box_width, patch_artist=True, showfliers=False, showcaps=False, zorder=2)
            for box in bp['boxes']:
                box.set_facecolor(color)
                box.set_alpha(box_fill_alpha)
                box.set_edgecolor('none')
                box.set_linewidth(0)
                vertices = box.get_path().vertices
                ax.plot(vertices[:, 0], vertices[:, 1], color=color, linewidth=box_linewidth, alpha=1.0, zorder=2.5)
            for whisker in bp['whiskers']:
                whisker.set(color=color, linewidth=box_linewidth, alpha=1.0)
            for median in bp['medians']:
                median.set(color=color, linewidth=box_linewidth, alpha=1.0)
            if outlier and len(outliers) > 0:
                scatter = ax.scatter(np.full(len(outliers), x_positions[i]), outliers, color=color, s=outlier_size, alpha=outlier_alpha, edgecolors='black', linewidths=0.5, zorder=3)
                if rasterize_outliers:
                    scatter.set_rasterized(True)
        y_max = max(within_vals.max() if len(within_vals) > 0 else 0, cross_vals.max() if len(cross_vals) > 0 else 0)
        y_min = min(within_vals.min() if len(within_vals) > 0 else 0, cross_vals.min() if len(cross_vals) > 0 else 0)
        y_range = y_max - y_min
        line_y = y_max + y_range * 0.05
        sig_y = line_y + y_range * 0.01
        ax.plot([0, 1], [line_y, line_y], 'k-', linewidth=0.4, zorder=4)
        ax.plot([0, 0], [line_y, line_y - y_range * 0.02], 'k-', linewidth=0.4, zorder=4)
        ax.plot([1, 1], [line_y, line_y - y_range * 0.02], 'k-', linewidth=0.4, zorder=4)
        ax.text(0.5, sig_y, sig_label, ha='center', va='bottom', fontsize=sig_label_fontsize, zorder=4)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(labels, fontsize=XLABEL_FONTSIZE, rotation=45, ha='right', rotation_mode='anchor')
        ax.set_ylabel(metric, fontsize=YLABEL_FONTSIZE, labelpad=1)
        ax.tick_params(axis='y', labelsize=YTICK_FONTSIZE, length=0.8, width=0.4, direction='out', pad=0.2)
        ax.tick_params(axis='x', length=1, width=0.4, pad=0.5)
        y_upper_limit = sig_y + y_range * 0.1
        ax.set_ylim(y_min - y_range * 0.05, y_upper_limit)
        y_lim_min, y_lim_max = ax.get_ylim()
        y_data_range = y_lim_max - y_lim_min
        if y_data_range <= 0.1:
            tick_interval = 0.02
        elif y_data_range <= 0.5:
            tick_interval = 0.05
        elif y_data_range <= 1.0:
            tick_interval = 0.1
        elif y_data_range <= 5.0:
            tick_interval = 0.5
        else:
            tick_interval = 1.0
        MAX_TICKS = 8
        tick_positions_y = []
        while True:
            first_tick = np.ceil(y_lim_min / tick_interval) * tick_interval
            tick_positions_y = []
            current_tick = first_tick
            while current_tick <= y_lim_max:
                tick_positions_y.append(current_tick)
                current_tick += tick_interval
            if len(tick_positions_y) <= MAX_TICKS:
                break
            if tick_interval < 0.05:
                tick_interval = 0.05
            elif tick_interval < 0.1:
                tick_interval = 0.1
            elif tick_interval < 0.2:
                tick_interval = 0.2
            elif tick_interval < 0.5:
                tick_interval = 0.5
            elif tick_interval < 1.0:
                tick_interval = 1.0
            elif tick_interval < 2.0:
                tick_interval = 2.0
            elif tick_interval < 5.0:
                tick_interval = 5.0
            else:
                tick_interval = tick_interval * 2
        if len(tick_positions_y) < 3:
            if tick_interval >= 1.0:
                tick_interval = 0.5
            elif tick_interval >= 0.5:
                tick_interval = 0.2
            elif tick_interval >= 0.2:
                tick_interval = 0.1
            elif tick_interval >= 0.1:
                tick_interval = 0.05
            else:
                tick_interval = tick_interval / 2
            first_tick = np.ceil(y_lim_min / tick_interval) * tick_interval
            tick_positions_y = []
            current_tick = first_tick
            while current_tick <= y_lim_max:
                tick_positions_y.append(current_tick)
                current_tick += tick_interval
        ax.set_yticks(tick_positions_y)
        tick_labels_y = []
        for tick in tick_positions_y:
            if abs(tick) < 1e-10:
                tick = 0.0
            label = f'{tick:.2f}'.rstrip('0').rstrip('.')
            tick_labels_y.append(label)
        ax.set_yticklabels(tick_labels_y)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_linewidth(spine_linewidth)
            spine.set_edgecolor('black')
        ax.set_xlim(-0.7, 1.7)
    fig.text(left_margin, 0.95, method, fontsize=TITLE_FONTSIZE, ha='left', va='top')
    out_path = os.path.join(output_dir, f'{method}_slice_comparison.pdf')
    plt.savefig(out_path, dpi=dpi, bbox_inches=None, pad_inches=0)
    if show:
        plt.show()
    else:
        plt.close()
    print(f'Figure saved to: {out_path}')

