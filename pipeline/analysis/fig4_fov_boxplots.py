"""Compare prediction metrics across FOV sizes."""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from sklearn.preprocessing import MinMaxScaler
from .fov_status import mask_non_evaluable_fovs


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'ARIAL'

def preprocess_PARnorm(csv_path, write_back=False):
    df = pd.read_csv(csv_path)
    if 'PARnorm_all_fov' in df.columns:
        return (mask_non_evaluable_fovs(df), 'PARnorm_all_fov')
    rep = df.shape[0]
    df['SCC_rank'] = df['SCC'].rank(ascending=False)
    df['PCC_rank'] = df['PCC'].rank(ascending=False)
    df['KLD_rank'] = df['KLD'].rank()
    df['JSD_rank'] = df['JSD'].rank()
    ranks_df = df[['SCC_rank', 'PCC_rank', 'KLD_rank', 'JSD_rank']].fillna(rep)
    df['PAR_all_fov'] = np.average(ranks_df.values, axis=1)
    scaler = MinMaxScaler()
    df['PARnorm_all_fov'] = scaler.fit_transform(df[['PAR_all_fov']])
    df = mask_non_evaluable_fovs(df)
    if write_back:
        df.to_csv(csv_path, index=False)
    return (df, 'PARnorm_all_fov')

def plot_box_violin_fov(method, csv_path, metrics, output_dir, fov_size_list=(500, 1000, 1500, 2000, 2500), box_width=0.5, violin_alpha=0.15, box_fill_alpha=0.4, outlier=False, outlier_size=10, outlier_alpha=0.6, rasterize_outliers=False, dpi=300, show=False, left_margin=0.3, right_margin=0.02, bottom_margin=0.2, top_margin=0.1, spine_linewidth=0.4, box_linewidth=0.4):
    """Plot metric distributions and median trends across FOV sizes."""
    os.makedirs(output_dir, exist_ok=True)
    MM_PER_INCH = 25.4
    FINAL_WIDTH_MM = 28.0
    ORIGINAL_WIDTH_MM = 28.0
    SCALE = FINAL_WIDTH_MM / ORIGINAL_WIDTH_MM
    TARGET_TITLE_PT = 7.0
    TARGET_LABEL_PT = 7.0
    TARGET_TICK_PT = 5.0
    TITLE_FONTSIZE = TARGET_TITLE_PT / SCALE
    LABEL_FONTSIZE = TARGET_LABEL_PT / SCALE
    TICK_FONTSIZE = TARGET_TICK_PT / SCALE
    ORIGINAL_HEIGHT_MM = 33.6
    FIG_WIDTH_IN = ORIGINAL_WIDTH_MM / MM_PER_INCH
    FIG_HEIGHT_IN = ORIGINAL_HEIGHT_MM / MM_PER_INCH
    color_map = {500: '#EE6677', 1000: '#4477AA', 1500: '#CCBB44', 2000: '#228833', 2500: '#66CCEE'}
    for metric in metrics:
        if metric == 'PARnorm':
            df, metric_col = preprocess_PARnorm(csv_path)
            metric_display = 'PARnorm'
        else:
            df = mask_non_evaluable_fovs(pd.read_csv(csv_path))
            metric_col = metric
            metric_display = metric
        df['FOV_numeric'] = (df['FOV_size'].astype(str).str.lower()
                             .str.replace('fov', '', regex=False)
                             .str.replace('_', '', regex=False).astype(int))
        fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), dpi=dpi)
        fig.patch.set_facecolor('white')
        plot_left = left_margin
        plot_bottom = bottom_margin
        plot_width = 1.0 - left_margin - right_margin
        plot_height = 1.0 - bottom_margin - top_margin
        ax = fig.add_axes([plot_left, plot_bottom, plot_width, plot_height])
        x_positions = np.arange(len(fov_size_list))
        med_x, med_y = ([], [])
        for i, fov in enumerate(fov_size_list):
            vals = df.loc[df['FOV_numeric'] == fov, metric_col].dropna()
            if vals.empty:
                continue
            color = color_map[fov]
            if outlier:
                q1 = np.percentile(vals, 25)
                q3 = np.percentile(vals, 75)
                iqr = q3 - q1
                lower_bound = q1 - 1.5 * iqr
                upper_bound = q3 + 1.5 * iqr
                outliers = vals[(vals < lower_bound) | (vals > upper_bound)]
            if len(vals) > 1 and np.ptp(vals.to_numpy()) > 0:
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
            med = np.median(vals)
            ax.scatter(x_positions[i], med, color='#C82423', s=1, zorder=3)
            med_x.append(x_positions[i])
            med_y.append(med)
        if len(med_x) > 1:
            ax.plot(med_x, med_y, linestyle='--', color='#C82423', linewidth=0.5, zorder=2.5)
        ax.set_xticks(x_positions)
        ax.set_xticklabels([str(f) for f in fov_size_list], fontsize=TICK_FONTSIZE, rotation=45, ha='right', rotation_mode='anchor')
        ax.set_xlabel('FOV size (pixel)', fontsize=LABEL_FONTSIZE)
        ax.set_ylabel(metric_display, fontsize=LABEL_FONTSIZE)
        ax.tick_params(axis='y', labelsize=TICK_FONTSIZE, length=0.8, width=0.4, direction='out', pad=0.2)
        ax.tick_params(axis='x', labelsize=TICK_FONTSIZE, length=1, width=0.4, pad=0.5)
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
        title_x = plot_left + plot_width / 2
        title_y = plot_bottom + plot_height + (1.0 - top_margin - (plot_bottom + plot_height)) / 2
        fig.text(title_x, title_y, method, ha='center', va='bottom', fontsize=TITLE_FONTSIZE)
        out_path = os.path.join(output_dir, f'{method}_{metric_display}_FOV_comparison.pdf')
        fig.savefig(out_path, dpi=dpi, bbox_inches=None, pad_inches=0)
        if show:
            plt.show()
        else:
            plt.close()
        print(f'Figure saved to: {out_path}')
