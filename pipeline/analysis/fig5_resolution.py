"""Aggregate resolution experiments and plot mean performance curves."""
from __future__ import annotations
from typing import Dict, List, Tuple, Iterable
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from matplotlib.pyplot import MultipleLocator
from matplotlib.ticker import FormatStrFormatter


scaler = MinMaxScaler()
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'ARIAL'
DEFAULT_PLOT_CFG = {'spotsize': {'figsize': (10, 5), 'simulations': ['proportional', 'dense', 'fixed_spotDistance'], 'xlabel': 'Spot Size', 'outfile_suffix': 'spot_size', 'x_transform': lambda sim, x: x}, 'spotdistance': {'figsize': (10, 5), 'simulations': ['proportional', 'dense', 'fixed_spotSize'], 'xlabel': 'Distance between Centers of Spots', 'outfile_suffix': 'spot_distance', 'x_transform': lambda sim, x: x * 100 / 55 if sim == 'proportional' else x}}

def process_data(fs, plt_path, method, dataset):
    dfs = []
    for myfile in fs:
        p = Path(myfile)
        try:
            repeat, simulation, res_part = p.parts[-3:]
        except ValueError:
            raise ValueError(f'path too short: {myfile}')
        resolution = int(res_part.split('.csv')[0].split('_')[1])
        df = pd.read_csv(myfile)
        df = df.assign(method=method, dataset=dataset, repeat=repeat, simulation=simulation, resolution=resolution)
        dfs.append(df)
    res = pd.concat(dfs, ignore_index=True)
    res.to_csv(f'{plt_path}/res_all.csv', index=False)
    metric_cols = ['PCC', 'SCC', 'JSD', 'KLD', 'SSIM', 'PSNR']
    res_mean = res.groupby(['method', 'dataset', 'gene', 'simulation', 'resolution'], as_index=False).agg({**{col: 'mean' for col in metric_cols}, 'repeat': lambda x: ','.join(map(str, x)), 'gt_counts': 'first'})
    res_mean.to_csv(f'{plt_path}/res_mean.csv', index=False)
    if method == 'TESLA':
        simulations = ['fixed_spotSize']
    else:
        simulations = ['proportional', 'dense', 'fixed_spotSize', 'fixed_spotDistance']
    res_rank_list = []
    for sim in simulations:
        df = res_mean[res_mean['simulation'] == sim].copy()
        if df.empty:
            raise ValueError(f'Missing required simulation scheme: {sim}; supply the complete method input set')
        df['SCC_rank'] = df['SCC'].rank(ascending=False, method='average', na_option='bottom')
        df['PCC_rank'] = df['PCC'].rank(ascending=False, method='average', na_option='bottom')
        df['KLD_rank'] = df['KLD'].rank(method='average', na_option='bottom')
        df['JSD_rank'] = df['JSD'].rank(method='average', na_option='bottom')
        df['SSIM_rank'] = df['SSIM'].rank(ascending=False, method='average', na_option='bottom')
        df['PSNR_rank'] = df['PSNR'].rank(ascending=False, method='average', na_option='bottom')
        rank_cols_acc = ['SCC_rank', 'PCC_rank', 'KLD_rank', 'JSD_rank']
        rank_cols_img = ['SSIM_rank', 'PSNR_rank']
        df['PAR'] = df[rank_cols_acc].mean(axis=1)
        df['PARnorm'] = scaler.fit_transform(df[['PAR']])
        df['PAR_img'] = df[rank_cols_img].mean(axis=1)
        df['PARnorm_img'] = scaler.fit_transform(df[['PAR_img']])
        res_rank_list.append(df)
    res_rank = pd.concat(res_rank_list, ignore_index=True)
    res_rank.to_csv(f'{plt_path}/res_rank.csv', index=False)
    return res_rank

def set_plot_style(style: str='white') -> None:
    """Apply the shared seaborn plotting style."""
    sns.set_style(style)

def style_axes(ax, x_major_step: int=5, spine_color: str='#CCCCCC', spine_lw: float=1.5) -> None:
    ax.xaxis.set_major_locator(MultipleLocator(x_major_step))
    ax.tick_params(axis='x', labelsize=9, rotation=45)
    ax.tick_params(axis='y', labelsize=10)
    for spine in ax.spines.values():
        spine.set_edgecolor(spine_color)
        spine.set_linewidth(spine_lw)

def summarize_metric_means(res_rank, simulation, metric):
    """Retain the source mean reduction; no uncertainty interval is calculated.

    The input contains gene-level values averaged across runs. With complete
    triplicates, the common repeat label is rep1,rep2,rep3; points summarize genes.
    """
    df = res_rank.loc[res_rank['simulation'] == simulation, ['resolution', 'repeat', metric]].copy()
    per_group = df.groupby(['resolution', 'repeat'], as_index=False)[metric].mean()
    rows = []
    for resolution, group in per_group.groupby('resolution'):
        if len(group) > 1:
            import warnings
            warnings.warn('Different replicate-availability labels occur at one resolution; inspect the retained input groups before interpreting the summary.', RuntimeWarning)
        rows.append({'resolution': resolution, metric: float(np.mean(group[metric].to_numpy(dtype=float)))})
    return pd.DataFrame(rows, columns=['resolution', metric]).sort_values('resolution').reset_index(drop=True)


def plot_curve(*, res_rank: pd.DataFrame, metric_name: str, resolution_type: str, method: str, dataset: str, plt_path: str | Path, color_map: Dict[str, str], simulations: Iterable[str], cfg: Dict=DEFAULT_PLOT_CFG, keep_fig: bool=False, show_legend: bool=True, legend_loc: str='upper right', final_width_mm: float=40.0, left_margin: float=0.18, right_margin: float=0.04, bottom_margin: float=0.3, top_margin: float=0.05, spine_linewidth: float=0.8) -> Tuple[plt.Figure, Path]:
    import warnings
    MM_PER_INCH = 25.4
    FINAL_WIDTH_MM = final_width_mm
    FINAL_HEIGHT_MM = 22.0
    ORIGINAL_WIDTH_MM = FINAL_WIDTH_MM
    ORIGINAL_HEIGHT_MM = FINAL_HEIGHT_MM
    SCALE = FINAL_WIDTH_MM / ORIGINAL_WIDTH_MM
    TARGET_LABEL_PT = 6.0
    TARGET_TICK_PT = 5.0
    LABEL_FONTSIZE = TARGET_LABEL_PT / SCALE
    TICK_FONTSIZE = TARGET_TICK_PT / SCALE
    FIG_WIDTH_IN = ORIGINAL_WIDTH_MM / MM_PER_INCH
    FIG_HEIGHT_IN = ORIGINAL_HEIGHT_MM / MM_PER_INCH
    all_y_values: list[float] = []
    if resolution_type not in cfg:
        raise ValueError(f'resolution_type must be one of {list(cfg.keys())}')
    c = cfg[resolution_type]
    out_dir = Path(plt_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
    plot_left = left_margin
    plot_bottom = bottom_margin
    plot_width = 1.0 - left_margin - right_margin
    plot_height = 1.0 - bottom_margin - top_margin
    ax = fig.add_axes([plot_left, plot_bottom, plot_width, plot_height])
    marker_size = 1 if resolution_type == 'spotsize' else 0.5
    plotted_any = False
    for sim in simulations:
        values = summarize_metric_means(res_rank, sim, metric_name)
        x_raw = c['x_transform'](sim, values['resolution'])
        y = pd.to_numeric(pd.Series(values[metric_name]), errors='coerce').to_numpy(dtype=float)
        x_arr = np.asarray(x_raw)
        if x_arr.ndim == 0:
            try:
                x = np.full(len(y), float(x_arr), dtype=float)
            except Exception:
                x = np.full(len(y), np.nan, dtype=float)
        else:
            x = pd.to_numeric(pd.Series(x_arr.ravel()), errors='coerce').to_numpy(dtype=float)
        if len(x) != len(y):
            raise ValueError(f'x and y length mismatch for sim={sim}: len(x)={len(x)}, len(y)={len(y)}')
        valid = np.isfinite(x) & np.isfinite(y)
        x_plot = x[valid]
        y_plot = y[valid]
        if x_plot.size == 0:
            warnings.warn(f'[plot_curve] sim={sim} has no finite points after filtering; skipped.', RuntimeWarning)
            continue
        plotted_any = True
        all_y_values.extend(y_plot.tolist())
        ax.plot(x_plot, y_plot, marker='o', linestyle='-', ms=marker_size, label=sim, linewidth=0.5, color=color_map.get(sim, None))
    if resolution_type.lower() == 'spotsize':
        xlabel = 'Spot Size (μm)'
    elif resolution_type.lower() == 'spotdistance':
        xlabel = 'Distance between Centers of Spots (μm)'
    else:
        xlabel = c.get('xlabel', resolution_type)
    is_xe = dataset.startswith('Xe')
    is_vihd = dataset.startswith('ViHD')
    is_parnorm = metric_name == 'PARnorm'
    sim_list = list(simulations)
    sim_name = sim_list[0] if len(sim_list) == 1 else ''
    if is_xe:
        if is_parnorm and sim_name in ['proportional', 'dense']:
            if resolution_type.lower() == 'spotsize':
                x_ticks = np.arange(20, 151, 10)
            else:
                x_ticks = np.arange(20, 281, 20)
        elif is_parnorm and sim_name == 'fixed_spotDistance':
            x_ticks = np.arange(20, 101, 10)
        elif is_parnorm and sim_name == 'fixed_spotSize':
            x_ticks = np.arange(55, 210, 15)
        elif resolution_type.lower() == 'spotsize':
            x_ticks = np.arange(20, 151, 10)
        elif method == 'TESLA':
            x_ticks = np.arange(55, 210, 15)
        else:
            x_ticks = np.arange(20, 281, 20)
    elif is_vihd:
        if is_parnorm and sim_name in ['proportional', 'dense']:
            if resolution_type.lower() == 'spotsize':
                x_ticks = np.arange(16, 145, 16)
            else:
                x_ticks = np.arange(16, 167, 32)
        elif is_parnorm and sim_name == 'fixed_spotDistance':
            x_ticks = np.arange(16, 97, 8)
        elif is_parnorm and sim_name == 'fixed_spotSize':
            x_ticks = np.arange(56, 217, 16)
        elif resolution_type.lower() == 'spotsize':
            x_ticks = np.arange(16, 145, 16)
        elif method == 'TESLA':
            x_ticks = np.arange(56, 217, 16)
        else:
            x_ticks = np.arange(16, 273, 32)
    else:
        x_ticks = None
    ax.set_xlabel(xlabel, fontsize=LABEL_FONTSIZE, labelpad=2)
    if x_ticks is not None:
        ax.set_xticks(x_ticks)
        ax.set_xlim(x_ticks[0] - (x_ticks[1] - x_ticks[0]) * 0.3, x_ticks[-1] + (x_ticks[1] - x_ticks[0]) * 0.3)
    ax.set_ylabel(metric_name, fontsize=LABEL_FONTSIZE, labelpad=2)
    y_all = np.asarray(all_y_values, dtype=float)
    y_all = y_all[np.isfinite(y_all)]
    if y_all.size == 0:
        warnings.warn('[plot_curve] No finite y values found. Using default y-limits.', RuntimeWarning)
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks([0.0, 0.5, 1.0])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    else:
        y_min = float(np.min(y_all))
        y_max = float(np.max(y_all))
        y_range = y_max - y_min
        candidate_steps = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
        step = candidate_steps[-1]
        ticks = None
        for s in candidate_steps:
            if not np.isfinite(s) or s <= 0:
                continue
            tick_min = np.floor(y_min / s) * s
            tick_max = np.ceil(y_max / s) * s
            if not (np.isfinite(tick_min) and np.isfinite(tick_max)):
                continue
            if tick_min >= tick_max:
                continue
            candidate_ticks = np.arange(tick_min, tick_max + s * 0.5, s)
            candidate_ticks = candidate_ticks[np.isfinite(candidate_ticks)]
            if 3 <= len(candidate_ticks) <= 8:
                step = float(s)
                ticks = candidate_ticks
                break
        if ticks is None:
            step = float(candidate_steps[0])
            center = (y_min + y_max) / 2.0
            ticks = np.array([center - step, center, center + step], dtype=float)
            ticks = ticks[np.isfinite(ticks)]
            ticks = np.unique(ticks)
            if ticks.size < 3:
                step = float(candidate_steps[1]) if len(candidate_steps) > 1 else step * 2.0
                ticks = np.array([center - step, center, center + step], dtype=float)
                ticks = ticks[np.isfinite(ticks)]
                ticks = np.unique(ticks)
        ylim_lower = float(ticks[0] - 0.3 * step)
        ylim_upper = float(ticks[-1] + 0.3 * step)
        if np.isfinite(ylim_lower) and np.isfinite(ylim_upper) and (ylim_lower < ylim_upper):
            ax.set_ylim(ylim_lower, ylim_upper)
        ax.set_yticks(ticks)
        tick_labels_y = []
        for tick in ticks:
            if abs(tick) < 1e-10:
                tick = 0.0
            label = f'{tick:.2f}'.rstrip('0').rstrip('.')
            tick_labels_y.append(label)
        ax.set_yticklabels(tick_labels_y)
    for spine in ax.spines.values():
        spine.set_linewidth(spine_linewidth)
        spine.set_edgecolor('black')
        spine.set_zorder(1)
    ax.tick_params(axis='y', which='both', labelsize=TICK_FONTSIZE, length=0.8, width=0.4, direction='out', pad=0.2, left=True, right=False, labelleft=True, labelright=False, colors='black', zorder=10)
    ax.tick_params(axis='x', which='both', labelsize=TICK_FONTSIZE, length=1.0, width=0.4, direction='out', pad=0.5, rotation=45, bottom=True, top=False, labelbottom=True, labeltop=False, colors='black', zorder=10)
    for label in ax.get_xticklabels():
        label.set_ha('right')
        label.set_rotation_mode('anchor')
    ax.grid(False)
    if show_legend and plotted_any:
        legend = ax.legend(loc=legend_loc, fontsize=TICK_FONTSIZE, frameon=True, edgecolor='black', fancybox=False)
        legend.get_frame().set_linewidth(0.8)
    if metric_name in PAR_METRICS and len(simulations) == 1:
        sim_suffix = list(simulations)[0]
        out_file = out_dir / f"{metric_name}_{sim_suffix}_{c['outfile_suffix']}.pdf"
    else:
        out_file = out_dir / f"{metric_name}_{c['outfile_suffix']}.pdf"
    fig.savefig(out_file, dpi=300, bbox_inches=None, pad_inches=0)
    if not keep_fig:
        plt.close(fig)
    return (fig, out_file)

def plot_all_metrics(dataset: str, *, res_rank: pd.DataFrame, metrics: Iterable[str], method: str, plt_path: str | Path, color_map: Dict[str, str], cfg: Dict=DEFAULT_PLOT_CFG, resolution_types: Tuple[str, ...]=('spotsize', 'spotdistance'), show_plots: bool=False, **plot_kwargs) -> List[Path]:
    figs: List[plt.Figure] = []
    outs: List[Path] = []
    for m in metrics:
        for rt in resolution_types:
            if method == 'TESLA':
                simulations = ['fixed_spotSize']
            else:
                simulations = cfg[rt]['simulations']
            if m in PAR_METRICS:
                for sim in simulations:
                    fig, out = plot_curve(res_rank=res_rank, metric_name=m, resolution_type=rt, method=method, dataset=dataset, plt_path=plt_path, color_map=color_map, simulations=[sim], cfg=cfg, keep_fig=show_plots, **plot_kwargs)
                    outs.append(out)
                    if show_plots:
                        figs.append(fig)
            else:
                fig, out = plot_curve(res_rank=res_rank, metric_name=m, resolution_type=rt, method=method, dataset=dataset, plt_path=plt_path, color_map=color_map, simulations=simulations, cfg=cfg, keep_fig=show_plots, **plot_kwargs)
                outs.append(out)
                if show_plots:
                    figs.append(fig)
    if show_plots:
        for fig in figs:
            fig.show()
        plt.show()
        for fig in figs:
            plt.close(fig)
    return outs

PAR_METRICS = {'PARnorm', 'PARnorm_img'}
