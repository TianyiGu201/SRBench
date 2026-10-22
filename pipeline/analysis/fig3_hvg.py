"""Compare HVG and non-HVG prediction and structural metrics."""
from __future__ import annotations
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, mannwhitneyu


def plot_hvg_comparison(method, dataset, merged_csv_dir, output_dir, box_width=0.5, violin_alpha=0.15, box_fill_alpha=0.4, outlier=False, outlier_size=10, outlier_alpha=0.6, rasterize_outliers=False, dpi=300, show=False, subplot_spacing=0.012, left_margin=0.05, right_margin=0.01, spine_linewidth=2.5, box_linewidth=2.0):
    """Plot six HVG/non-HVG metric comparisons in a horizontal panel."""
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
    csv_path = os.path.join(merged_csv_dir, f'{dataset}_{method}_classified.csv')
    df = pd.read_csv(csv_path)
    metrics = ['PCC', 'SCC', 'JSD', 'KLD', 'SSIM', 'PSNR']
    color_hvg = '#F28500'
    color_non_hvg = '#2878B5'
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
        hvg_vals = df.loc[df['hvg_status'] == 'HVG', metric].dropna()
        non_hvg_vals = df.loc[df['hvg_status'] == 'non-HVG', metric].dropna()
        if len(hvg_vals) > 0 and len(non_hvg_vals) > 0:
            statistic, p_value = mannwhitneyu(hvg_vals, non_hvg_vals, alternative='two-sided')
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
        groups = [hvg_vals, non_hvg_vals]
        colors = [color_hvg, color_non_hvg]
        labels = ['HVG', 'non-HVG']
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
        y_max = max(hvg_vals.max() if len(hvg_vals) > 0 else 0, non_hvg_vals.max() if len(non_hvg_vals) > 0 else 0)
        y_min = min(hvg_vals.min() if len(hvg_vals) > 0 else 0, non_hvg_vals.min() if len(non_hvg_vals) > 0 else 0)
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
    out_path = os.path.join(output_dir, f'{dataset}_{method}_HVG_comparison.pdf')
    plt.savefig(out_path, dpi=dpi, bbox_inches=None, pad_inches=0)
    if show:
        plt.show()
    else:
        plt.close()
    print(f'Figure saved to: {out_path}')

