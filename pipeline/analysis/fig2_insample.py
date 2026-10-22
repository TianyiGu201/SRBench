"""Aggregate repeated in-sample metrics, calculate robustness and rank methods."""
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, Polygon
from matplotlib.colors import TwoSlopeNorm
from glob import glob

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
DEFAULT_METHODS = ['iSTAR', 'scstGCN', 'TESLA', 'XFuse', 'SuperST_t', 'SuperST_t_g']
METRICS = ['PCC', 'SCC', 'KLD', 'JSD', 'SSIM', 'PSNR']
LOWER_BETTER = {'JSD', 'KLD'}
HIGHER_BETTER = set(METRICS) - LOWER_BETTER
ASCENDING = {'PCC': False, 'SCC': False, 'KLD': True, 'JSD': True, 'SSIM': False, 'PSNR': False}

def minmax_norm_series(x: pd.Series) -> pd.Series:
    x = x.astype(float)
    mn, mx = (x.min(), x.max())
    if pd.isna(mn) or pd.isna(mx) or mx == mn:
        return pd.Series(np.zeros(len(x), dtype=float), index=x.index)
    return (x - mn) / (mx - mn)

def get_rep_file_path(base_dir: Path, method: str, dataset: str, rep: str) -> Path:
    base_dir = Path(base_dir)
    if dataset.startswith('Xe'):
        return base_dir / method / dataset / 'raw_resolution' / rep / 'fixed_spotSize' / 'spotDistance_100.csv'
    elif dataset.startswith('ViHD'):
        return base_dir / method / dataset / 'raw_resolution' / rep / 'fixed_spotSize' / 'spotDistance_96.csv'
    elif dataset.startswith('cos'):
        return base_dir / method / dataset / f'{rep}.csv'
    else:
        raise ValueError(f"Unknown dataset prefix for dataset='{dataset}'. Expected it to start with Xe / ViHD / Cos.")

def compute_rep_robustness_overall_sd(long_df: pd.DataFrame, metrics, pa_metrics=('PCC', 'SCC', 'KLD', 'JSD'), sa_metrics=('SSIM', 'PSNR')):
    """Calculate robustness from normalized gene-by-run metrics.

    For each metric, min-max normalize across all gene/run observations, average
    over genes within each run, then calculate the population SD of the run means.
    PArobustness and SArobustness are one minus the mean SD for their metric sets.
    """
    df = long_df.copy()
    metric_sd_over_reps = {}
    for m in metrics:
        norm_col = f'{m}__norm'
        df[norm_col] = minmax_norm_series(df[m])
        rep_means = df.groupby('rep')[norm_col].mean()
        sd_over_reps = float(rep_means.std(ddof=0))
        if np.isnan(sd_over_reps):
            sd_over_reps = 0.0
        metric_sd_over_reps[m] = sd_over_reps

    def safe_mean(metric_list):
        vals = [metric_sd_over_reps[m] for m in metric_list if m in metric_sd_over_reps]
        if len(vals) == 0:
            return np.nan
        return float(np.mean(vals))
    pa_sd_mean = safe_mean(pa_metrics)
    sa_sd_mean = safe_mean(sa_metrics)
    PArobustness = 1 - pa_sd_mean if not np.isnan(pa_sd_mean) else np.nan
    SArobustness = 1 - sa_sd_mean if not np.isnan(sa_sd_mean) else np.nan
    return (PArobustness, SArobustness, metric_sd_over_reps)

def merge_by_method_datasets_and_compute_robustness(base_dir, out_dir, method_datasets: dict, reps=('rep1', 'rep2', 'rep3'), metrics=('PCC', 'SCC', 'KLD', 'JSD', 'SSIM', 'PSNR'), pa_metrics=('PCC', 'SCC', 'KLD', 'JSD'), sa_metrics=('SSIM', 'PSNR'), save_metric_sd_over_reps=True):
    base_dir = Path(base_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    required_cols = ['gene'] + list(metrics)
    all_saved = []
    for method, datasets in method_datasets.items():
        for dataset in datasets:
            rep_dfs = []
            for rep in reps:
                fpath = get_rep_file_path(base_dir, method, dataset, rep)
                if not fpath.exists():
                    raise FileNotFoundError(f'Missing file for method={method}, dataset={dataset}, rep={rep}\nExpected: {fpath}')
                df = pd.read_csv(fpath)
                missing = [c for c in required_cols if c not in df.columns]
                if missing:
                    raise ValueError(f'File {fpath} missing columns: {missing}\nColumns found: {list(df.columns)}')
                use_cols = ['gene']
                if 'gt_counts' in df.columns:
                    use_cols.append('gt_counts')
                use_cols += list(metrics)
                df_use = df[use_cols].copy()
                df_use['rep'] = rep
                rep_dfs.append(df_use)
            long_df = pd.concat(rep_dfs, ignore_index=True)
            PArobustness, SArobustness, metric_sd_over_reps = compute_rep_robustness_overall_sd(long_df=long_df, metrics=metrics, pa_metrics=pa_metrics, sa_metrics=sa_metrics)
            agg_dict = {m: 'mean' for m in metrics}
            if 'gt_counts' in long_df.columns:
                agg_dict['gt_counts'] = 'mean'
            gene_mean = long_df.groupby('gene', as_index=False).agg(agg_dict)
            gene_mean['dataset'] = dataset
            gene_mean['method'] = method
            gene_mean['PArobustness'] = PArobustness
            gene_mean['SArobustness'] = SArobustness
            if save_metric_sd_over_reps:
                for m, sd_val in metric_sd_over_reps.items():
                    gene_mean[f'{m}_repSD'] = sd_val
            out_path = out_dir / f'{dataset}_{method}.csv'
            gene_mean.to_csv(out_path, index=False)
            all_saved.append(out_path)
            print(f'[Saved] {out_path}')
    return all_saved

def find_csvs_for_dataset(input_dir: str, dataset_name: str, methods=None):
    """Find files named <dataset>_<method>.csv in the input directory."""
    if methods is None:
        methods = DEFAULT_METHODS
    found = {}
    for m in methods:
        pattern = os.path.join(input_dir, f'{dataset_name}_{m}.csv')
        hits = glob(pattern)
        if hits:
            found[m] = sorted(hits)[0]
    missing = [m for m in methods if m not in found]
    return (found, missing)

def load_and_merge_dataset(input_dir: str, dataset_name: str, methods=None, gene_col='gene'):
    """Merge method tables and infer dataset/method identifiers from filenames."""
    if methods is None:
        methods = DEFAULT_METHODS
    found, missing = find_csvs_for_dataset(input_dir, dataset_name, methods=methods)
    if missing:
        raise FileNotFoundError(f'[{dataset_name}] Missing method files: {missing}\nExpected filename format: {dataset_name}_<method>.csv\nCheck directory: {input_dir}')
    dfs = []
    for m in methods:
        f = found[m]
        df = pd.read_csv(f)
        if gene_col not in df.columns:
            raise ValueError(f"File {f} does not contain gene column '{gene_col}'")
        df['dataset'] = dataset_name
        df['method'] = m
        for c in METRICS + ['PArobustness', 'SArobustness']:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors='coerce')
        dfs.append(df)
    out = pd.concat(dfs, ignore_index=True)
    return out

def _unique_rank_by_sort(values: pd.Series, methods: pd.Series, method_order):
    """Assign ranks 1..N, breaking ties by the supplied method order."""
    order_map = {m: i for i, m in enumerate(method_order)}
    v = values.astype(float).copy()
    v2 = v.values.copy()
    v2[~np.isfinite(v2)] = np.inf
    mkey = np.array([order_map.get(m, 10 ** 9) for m in methods.values], dtype=np.int64)
    idx = np.lexsort((mkey, v2))
    ranks = np.empty_like(idx, dtype=np.int64)
    ranks[idx] = np.arange(1, len(idx) + 1)
    return ranks

def merge_and_rank_strict(df: pd.DataFrame, gene_col='gene', method_order=None, require_n_methods=6):
    """Add PArank and SArank within each dataset, with lower JSD/KLD ranked better."""
    if method_order is None:
        method_order = DEFAULT_METHODS
    need = set(['dataset', 'method', gene_col] + METRICS)
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f'Input dataframe is missing required columns: {miss}')
    gene_level = df.groupby(['dataset', 'method', gene_col], as_index=False)[METRICS].mean()
    method_rows = []
    for ds, sub in gene_level.groupby('dataset', sort=False):
        methods_in = sub['method'].dropna().unique().tolist()
        if require_n_methods is not None and len(methods_in) != require_n_methods:
            raise ValueError(f'[{ds}] number of methods={len(methods_in)}, required={require_n_methods}.\nobserved methods: {methods_in}')
        tmp = sub.copy()
        for m in METRICS:
            tmp[f'global_rank_{m}'] = tmp[m].rank(ascending=ASCENDING[m], method='average')
        method_mean = tmp.groupby(['dataset', 'method'], as_index=False)[[f'global_rank_{m}' for m in METRICS]].mean().rename(columns={f'global_rank_{m}': f'mean_gene_rank_{m}' for m in METRICS})
        for m in METRICS:
            method_mean[f'rank_{m}'] = _unique_rank_by_sort(values=method_mean[f'mean_gene_rank_{m}'], methods=method_mean['method'], method_order=method_order).astype(float)
        method_mean['PA_rank_mean'] = method_mean[[f'rank_{m}' for m in ['PCC', 'SCC', 'JSD', 'KLD']]].mean(axis=1)
        method_mean['SA_rank_mean'] = method_mean[[f'rank_{m}' for m in ['SSIM', 'PSNR']]].mean(axis=1)
        method_mean['PArank'] = _unique_rank_by_sort(values=method_mean['PA_rank_mean'], methods=method_mean['method'], method_order=method_order).astype(float)
        method_mean['SArank'] = _unique_rank_by_sort(values=method_mean['SA_rank_mean'], methods=method_mean['method'], method_order=method_order).astype(float)
        method_rows.append(method_mean[['dataset', 'method', 'PArank', 'SArank']])
    ranks = pd.concat(method_rows, ignore_index=True)
    out = df.merge(ranks, on=['dataset', 'method'], how='left')
    return out

def _reverse_cmap_name(name: str) -> str:
    """Return the reversed colormap name."""
    if name.endswith('_r'):
        return name[:-2]
    return name + '_r'

def build_df_for_plot(input_dir: str, dataset_names, methods=None, gene_col='gene', require_n_methods=6):
    if methods is None:
        methods = DEFAULT_METHODS
    if isinstance(dataset_names, str):
        dataset_names = [dataset_names]
    dfs = []
    for ds in dataset_names:
        d = load_and_merge_dataset(input_dir, ds, methods=methods, gene_col=gene_col)
        dfs.append(d)
    df_all = pd.concat(dfs, ignore_index=True)
    df_all = merge_and_rank_strict(df_all, gene_col=gene_col, method_order=methods, require_n_methods=require_n_methods)
    return df_all

def _format_tick(val):
    """Format ticks with at most two decimal places and no trailing zeros."""
    if abs(val) < 1e-10:
        return '0'
    s = f'{val:.2f}'
    s = s.rstrip('0').rstrip('.')
    return s

def plot_panel_from_df(df: pd.DataFrame, fig_width_in=16.0, col_widths=None, color_specs=None, xlim_dict=None, xticks_dict=None, row_gap=0.9, dataset_gap=0.35, dataset_heights=None, height_unit_in=0.12, extra_top_in=1.0, extra_bottom_in=0.6, method_order=None, metric_gap=0.24, title_gap=0.012, title_offset=0.02, seed=7, plot_name=None, show=True, font_sizes=None, robust_pa_xlim=(0.8, 1.0), robust_st_xlim=(0.8, 1.0), base_cmap='RdBu_r', autoscale_kde=True, autoscale_pad_ratio=0.0, colorbar_row_ratio=0.9, colorbar_height=0.015, colorbar_y_offset=0.0, rank_col_width_ratio=0.5, per_dataset_xlim=True, colorbar_extend_ratio=0.06):
    """Plot metric distributions, ranks and robustness across datasets.

    xlim_dict accepts metric names or (dataset, metric) keys. per_dataset_xlim
    controls automatic dataset-specific limits. color_specs defines colormaps
    and normalization bounds; robust_pa_xlim and robust_st_xlim set robustness
    limits. colorbar_extend_ratio sets the fraction used for extension arrows.
    """
    if font_sizes is None:
        font_sizes = {}
    fs_method = font_sizes.get('method_labels', 10)
    fs_col_title = font_sizes.get('col_titles', 11)
    fs_group_title = font_sizes.get('group_titles', 13)
    fs_ticks = font_sizes.get('ticks', 9)
    fs_dataset = font_sizes.get('dataset_band', 10)
    fs_cbar_title = font_sizes.get('colorbar_title', 9)
    fs_cbar_ticks = font_sizes.get('colorbar_ticks', 8)
    np.random.seed(seed)
    kde_metrics_pa = ['PCC', 'SCC', 'JSD', 'KLD']
    kde_metrics_st = ['SSIM', 'PSNR']
    kde_metrics_all = kde_metrics_pa + kde_metrics_st
    all_cols = ['PCC', 'SCC', 'JSD', 'KLD', 'ROBUST_PA', 'RANK_PA', 'SSIM', 'PSNR', 'ROBUST_ST', 'RANK_ST']

    def truncate_cmap(cmap, minval=0.55, maxval=0.9, n=256):
        return LinearSegmentedColormap.from_list('trunc_red', cmap(np.linspace(minval, maxval, n)))

    def _get_cmap(cmap_like):
        if isinstance(cmap_like, mpl.colors.Colormap):
            return cmap_like
        if hasattr(mpl, 'colormaps'):
            return mpl.colormaps.get_cmap(cmap_like)
        return plt.get_cmap(cmap_like)

    def _gaussian_kde_1d(samples, grid, bandwidth=None):
        x = np.asarray(samples, dtype=float)
        x = x[np.isfinite(x)]
        n = x.size
        if n == 0:
            return np.zeros_like(grid)
        std = np.std(x, ddof=1) if n > 1 else 0.0
        if bandwidth is None:
            bandwidth = max(1e-09, 1.06 * (std if std > 0 else 1e-09) * n ** (-1 / 5))
        u = (grid[:, None] - x[None, :]) / bandwidth
        phi = np.exp(-0.5 * u * u) / np.sqrt(2 * np.pi)
        dens = phi.mean(axis=1) / bandwidth
        return dens

    def _plot_upper_half_violin(ax, data, y_pos, *, xlim, max_width=0.6, cmap, norm, z=1):
        """Draw a median-colored upper violin with a white median line inside its outline."""
        data = np.asarray(data, dtype=float)
        data = data[np.isfinite(data)]
        lo, hi = xlim
        grid = np.linspace(lo, hi, 240)
        dens = _gaussian_kde_1d(data, grid)
        dens = dens / dens.max() * max_width if dens.max() > 0 else np.zeros_like(dens)
        med = np.median(data) if data.size else (lo + hi) / 2
        color = _get_cmap(cmap)(norm(med))
        y_edge = y_pos - dens
        poly = ax.fill(np.r_[grid, grid[::-1]], np.r_[y_edge, np.full_like(grid, y_pos)[::-1]], facecolor=color, edgecolor='none', alpha=1, zorder=z)[0]
        dens_med = float(np.interp(med, grid, dens)) if dens.size else 0.0
        y_end = y_pos - dens_med
        lw = 1.6
        half_lw_px = lw / 2.0 * (ax.figure.dpi / 72.0)
        x_disp, y_disp = ax.transData.transform((med, y_end))
        _, y_end_inner = ax.transData.inverted().transform((x_disp, y_disp + half_lw_px))
        y_end_inner = min(y_pos, max(y_end_inner, y_end))
        line, = ax.plot([med, med], [y_pos, y_end_inner], color='white', lw=lw, zorder=z + 0.2, solid_capstyle='butt', clip_on=True)
        line.set_clip_path(poly)
        ax.plot(grid, y_edge, color='black', lw=1.2, zorder=z + 0.1)

    def _ds_kind(db):
        if isinstance(db, str):
            if db.startswith('Xe'):
                return 'Xenium'
            if db.startswith('ViHD'):
                return 'VisiumHD'
            if db.startswith('Co'):
                return 'CosMx'
        return 'Other'
    kind_color = {'Xenium': '#E6F2FF', 'VisiumHD': '#E8F7E1', 'CosMx': '#FDECEF', 'Other': '#EEEEEE'}
    if col_widths is None:
        col_widths = [1.2, 0.95, 1.2]
    assert len(col_widths) == 3
    robust_bar_height_ratio = 0.6 * row_gap
    rank_col_width = col_widths[2] * rank_col_width_ratio
    width_ratios = [col_widths[0]] * 4 + [col_widths[2]] + [rank_col_width] + [col_widths[0]] * 2 + [col_widths[2]] + [rank_col_width]
    if method_order is None:
        method_order = list(df['method'].dropna().unique())
    if color_specs is None:
        color_specs = {}
    if xlim_dict is None:
        xlim_dict = {}
    if xticks_dict is None:
        xticks_dict = {}
    db_groups = []
    seen = set()
    for dbn in df['dataset']:
        if pd.isna(dbn) or dbn in seen:
            continue
        seen.add(dbn)
        methods_in = [m for m in method_order if m in set(df.loc[df['dataset'] == dbn, 'method'])]
        if methods_in:
            db_groups.append((dbn, methods_in))
    n_ds = len(db_groups)
    if n_ds == 0:
        raise ValueError('No dataset rows to plot.')

    def _metric_vals(metric: str, dataset: str=None):
        """Return metric values, optionally restricted to one dataset."""
        if dataset is None:
            vals = pd.to_numeric(df[metric], errors='coerce').dropna().values
        else:
            vals = pd.to_numeric(df.loc[df['dataset'] == dataset, metric], errors='coerce').dropna().values
        return vals

    def _auto_range_for_metric(metric: str, dataset: str=None):
        """Calculate the data range for a metric."""
        vals = _metric_vals(metric, dataset)
        if vals.size == 0:
            return (0.0, 1.0)
        lo = float(np.min(vals))
        hi = float(np.max(vals))
        if not np.isfinite(lo) or not np.isfinite(hi):
            return (0.0, 1.0)
        if abs(hi - lo) < 1e-12:
            lo -= 1e-06
            hi += 1e-06
        pad = (hi - lo) * float(autoscale_pad_ratio)
        return (lo - pad, hi + pad)

    def _get_xlim(metric: str, dataset: str):
        """Return the x-axis limits for a dataset and metric."""
        if (dataset, metric) in xlim_dict:
            return xlim_dict[dataset, metric]
        if metric in xlim_dict:
            return xlim_dict[metric]
        if per_dataset_xlim:
            return _auto_range_for_metric(metric, dataset)
        else:
            return _auto_range_for_metric(metric, None)

    def _get_xticks(metric: str, dataset: str, xlim):
        """Return the x-axis ticks for a dataset and metric."""
        if (dataset, metric) in xticks_dict:
            return xticks_dict[dataset, metric]
        if metric in xticks_dict:
            return xticks_dict[metric]
        lo, hi = xlim
        return list(np.linspace(lo, hi, 5))
    data_range_dict = {}
    for m in kde_metrics_all:
        data_range_dict[m] = _auto_range_for_metric(m, None)
        if m in LOWER_BETTER:
            cmap_m = _reverse_cmap_name(base_cmap)
            extend_m = 'min'
        else:
            cmap_m = base_cmap
            extend_m = 'max'
        auto_lo, auto_hi = data_range_dict[m]
        color_specs.setdefault(m, {'cmap': cmap_m, 'vmin': auto_lo, 'vmax': auto_hi, 'extend': extend_m})
    xlim_dict.setdefault('ROBUST_PA', robust_pa_xlim)
    xlim_dict.setdefault('ROBUST_ST', robust_st_xlim)
    xlim_dict.setdefault('RANK_PA', (0.0, 1.0))
    xlim_dict.setdefault('RANK_ST', (0.0, 1.0))
    color_specs.setdefault('PA_TITLE_BAND', '#EAF1FF')
    color_specs.setdefault('ST_TITLE_BAND', '#EAFEEA')
    color_specs.setdefault('RANK_PA', {'cmap': 'RdBu', 'vmin': 1, 'vmax': 6})
    color_specs.setdefault('RANK_ST', {'cmap': 'RdBu', 'vmin': 1, 'vmax': 6})
    if dataset_heights is None:
        heights = [20.0] * n_ds
    elif isinstance(dataset_heights, (int, float)):
        heights = [float(dataset_heights)] * n_ds
    elif isinstance(dataset_heights, (list, tuple)):
        if len(dataset_heights) != n_ds:
            raise ValueError(f'dataset_heights length={len(dataset_heights)},number of datasets={n_ds}')
        heights = [float(x) for x in dataset_heights]
    elif isinstance(dataset_heights, dict):
        heights = [float(dataset_heights.get(dbn, 20.0)) for dbn, _ in db_groups]
    else:
        raise TypeError('dataset_heights must be None, a scalar, a list or a dictionary')
    total_ratio = float(sum(heights) + colorbar_row_ratio)
    fig_height_in = height_unit_in * total_ratio + extra_top_in + extra_bottom_in
    fig_size = (float(fig_width_in), float(fig_height_in))
    tmp = df.copy()
    for c in ['PArank', 'SArank', 'PArobustness', 'SArobustness']:
        if c in tmp.columns:
            tmp[c] = pd.to_numeric(tmp[c], errors='coerce')
    gp = tmp.groupby(['dataset', 'method'], observed=True)
    mean_parank = gp['PArank'].mean().unstack('method') if 'PArank' in tmp.columns else pd.DataFrame()
    mean_sarank = gp['SArank'].mean().unstack('method') if 'SArank' in tmp.columns else pd.DataFrame()
    mean_parob = gp['PArobustness'].mean().unstack('method') if 'PArobustness' in tmp.columns else pd.DataFrame()
    mean_sarob = gp['SArobustness'].mean().unstack('method') if 'SArobustness' in tmp.columns else pd.DataFrame()
    fig = plt.figure(figsize=fig_size)
    gs = GridSpec(nrows=n_ds + 1, ncols=10, width_ratios=width_ratios, height_ratios=heights + [colorbar_row_ratio], wspace=metric_gap, hspace=dataset_gap, figure=fig)
    axes_map = {}

    def add_dataset_band(fig, ax_row0, dbn, facecolor):
        bb = ax_row0.get_position(fig)
        band_w = 0.028
        left = max(0.005, bb.x0 - band_w - 0.06)
        band_ax = fig.add_axes([left, bb.y0, band_w, bb.y1 - bb.y0])
        band_ax.set_facecolor(facecolor)
        band_ax.set_xticks([])
        band_ax.set_yticks([])
        for sp in band_ax.spines.values():
            sp.set_visible(False)
        if '_' in dbn:
            first, rest = dbn.split('_', 1)
            label = f'{first}\n{rest}'
        else:
            label = dbn
        band_ax.text(0.5, 0.5, label, rotation=90, ha='center', va='center', multialignment='left', fontsize=fs_dataset, color='black')
    grid_color = '#BFBFBF'
    grid_lw = 0.6
    for i, (dbn, methods_in) in enumerate(db_groups):
        n_m = len(methods_in)
        y_pos = np.array([1.0 + k * row_gap for k in range(n_m)])
        ymin_local, ymax_local = (y_pos.min() - 0.6 * row_gap, y_pos.max() + 0.6 * row_gap)
        for j, key in enumerate(all_cols):
            ax = fig.add_subplot(gs[i, j])
            axes_map[i, j] = ax
            for sp in ax.spines.values():
                sp.set_visible(True)
            ax.tick_params(axis='both', length=0, labelsize=fs_ticks)
            ax.set_axisbelow(True)
            ax.set_yticks(y_pos)
            if j == 0:
                ax.set_yticklabels(methods_in, fontsize=fs_method)
            else:
                ax.set_yticklabels([])
            ax.set_ylim(ymax_local, ymin_local)
            if i == 0:
                if key in ['RANK_PA', 'RANK_ST']:
                    ax.set_title('Rank', fontsize=fs_col_title)
                elif key in ['ROBUST_PA', 'ROBUST_ST']:
                    ax.set_title('Robustness', fontsize=fs_col_title)
                else:
                    ax.set_title(key, fontsize=fs_col_title)
            if key in kde_metrics_all:
                spec = color_specs[key]
                xlim = _get_xlim(key, dbn)
                ax.set_xlim(*xlim)
                xtks = _get_xticks(key, dbn, xlim)
                ax.set_xticks(xtks)
                ax.set_xticklabels([_format_tick(x) for x in xtks], rotation=90, fontsize=fs_ticks, ha='center', va='top')
                ax.grid(True, axis='x', color=grid_color, linewidth=grid_lw)
                ax.grid(True, axis='y', color=grid_color, linewidth=grid_lw)
                if 'vcenter' in spec:
                    norm = TwoSlopeNorm(vmin=spec['vmin'], vcenter=spec['vcenter'], vmax=spec['vmax'])
                else:
                    norm = mpl.colors.Normalize(vmin=spec['vmin'], vmax=spec['vmax'])
                for yi, mth in zip(y_pos, methods_in):
                    arr = pd.to_numeric(df.loc[(df['dataset'] == dbn) & (df['method'] == mth), key], errors='coerce').dropna().values
                    if arr.size == 0:
                        continue
                    _plot_upper_half_violin(ax, arr, yi, xlim=xlim, cmap=spec['cmap'], norm=norm, z=yi)
            elif key == 'ROBUST_PA':
                xlim = xlim_dict.get('ROBUST_PA', robust_pa_xlim)
                xpad = 0.005
                xlim_view = (xlim[0] - xpad, xlim[1] + xpad)
                ax.set_xlim(*xlim_view)
                xtks = [xlim[0], (xlim[0] + xlim[1]) / 2, xlim[1]]
                ax.set_xticks(xtks)
                ax.set_xticklabels([_format_tick(x) for x in xtks], rotation=90, fontsize=fs_ticks, ha='center', va='top')
                ax.grid(True, axis='x', color=grid_color, linewidth=grid_lw)
                for sp in ax.spines.values():
                    sp.set_zorder(10)
                vals = []
                for mth in methods_in:
                    v = np.nan
                    if not mean_parob.empty and dbn in mean_parob.index and (mth in mean_parob.columns):
                        v = mean_parob.loc[dbn, mth]
                    vals.append(v)
                vals = np.array(vals, dtype=float)
                vals[np.isnan(vals)] = xlim[0]
                bar_lengths = vals - xlim[0]
                bar_lengths[bar_lengths < 0] = 0
                from matplotlib.colors import LinearSegmentedColormap
                cmap = LinearSegmentedColormap.from_list('warm_brick', ['#e7e0db', '#e7d2c8', '#d7a59a', '#c47b6e', '#a9554c'])
                norm = mpl.colors.PowerNorm(gamma=1.8, vmin=xlim[0], vmax=xlim[1])
                colors = cmap(norm(vals))
                bars = ax.barh(y_pos, bar_lengths, left=xlim[0], height=robust_bar_height_ratio, color=colors, edgecolor='none', zorder=3, clip_on=True)
                for b in bars:
                    b.set_clip_path(ax.patch)
            elif key == 'ROBUST_ST':
                xlim = xlim_dict.get('ROBUST_ST', robust_st_xlim)
                xpad = 0.01
                xlim_view = (xlim[0] - xpad, xlim[1] + xpad)
                ax.set_xlim(*xlim_view)
                xtks = [xlim[0], (xlim[0] + xlim[1]) / 2, xlim[1]]
                ax.set_xticks(xtks)
                ax.set_xticklabels([_format_tick(x) for x in xtks], rotation=90, fontsize=fs_ticks, ha='center', va='top')
                ax.grid(True, axis='x', color=grid_color, linewidth=grid_lw)
                for sp in ax.spines.values():
                    sp.set_zorder(10)
                vals = []
                for mth in methods_in:
                    v = np.nan
                    if not mean_sarob.empty and dbn in mean_sarob.index and (mth in mean_sarob.columns):
                        v = mean_sarob.loc[dbn, mth]
                    vals.append(v)
                vals = np.array(vals, dtype=float)
                vals[np.isnan(vals)] = xlim[0]
                bar_lengths = vals - xlim[0]
                bar_lengths[bar_lengths < 0] = 0
                from matplotlib.colors import LinearSegmentedColormap
                cmap = LinearSegmentedColormap.from_list('warm_brick', ['#e7e0db', '#e7d2c8', '#d7a59a', '#c47b6e', '#a9554c'])
                norm = mpl.colors.PowerNorm(gamma=1.8, vmin=xlim[0], vmax=xlim[1])
                colors = cmap(norm(vals))
                bars = ax.barh(y_pos, bar_lengths, left=xlim[0], height=robust_bar_height_ratio, color=colors, edgecolor='none', zorder=3, clip_on=True)
                for b in bars:
                    b.set_clip_path(ax.patch)
            elif key == 'RANK_PA':
                spec = color_specs['RANK_PA']
                rnorm = mpl.colors.Normalize(vmin=spec['vmin'], vmax=spec['vmax'])
                cmap = _get_cmap(spec['cmap'])
                ax.set_xlim(*xlim_dict['RANK_PA'])
                ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                ax._rank_fixed_info = ('PA', dbn, methods_in, rnorm, cmap)
            elif key == 'RANK_ST':
                spec = color_specs['RANK_ST']
                rnorm = mpl.colors.Normalize(vmin=spec['vmin'], vmax=spec['vmax'])
                cmap = _get_cmap(spec['cmap'])
                ax.set_xlim(*xlim_dict['RANK_ST'])
                ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                ax._rank_fixed_info = ('ST', dbn, methods_in, rnorm, cmap)
        add_dataset_band(fig, axes_map[i, 0], dbn, kind_color.get(_ds_kind(dbn), '#EEEEEE'))
    fig.canvas.draw()
    dpi = fig.dpi
    for (i, j), ax in list(axes_map.items()):
        if hasattr(ax, '_rank_fixed_info'):
            kind, dbn, methods_in, rnorm, cmap = ax._rank_fixed_info
            y_ticks = axes_map[i, 0].get_yticks()
            robust_col = 4 if kind == 'PA' else 8
            robust_ax = axes_map[i, robust_col]
            bbox = robust_ax.get_window_extent(renderer=fig.canvas.get_renderer())
            ylim = robust_ax.get_ylim()
            ylim_range = abs(ylim[0] - ylim[1])
            bar_height_px = robust_bar_height_ratio / ylim_range * bbox.height
            bar_height_pts = bar_height_px * 72.0 / dpi
            s_fixed = bar_height_pts ** 2
            rank_values = []
            valid_yticks = []
            for yi, mth in zip(y_ticks, methods_in):
                if kind == 'PA':
                    if not mean_parank.empty and dbn in mean_parank.index and (mth in mean_parank.columns):
                        r = mean_parank.loc[dbn, mth]
                    else:
                        continue
                elif not mean_sarank.empty and dbn in mean_sarank.index and (mth in mean_sarank.columns):
                    r = mean_sarank.loc[dbn, mth]
                else:
                    continue
                if not np.isfinite(r):
                    continue
                ax.scatter([0.5], [yi], s=s_fixed, c=[cmap(rnorm(r))], marker='s', edgecolors='k', linewidths=0.6, zorder=6)
                valid_yticks.append(yi)
                rank_values.append(r)
            if valid_yticks:
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position('right')
                ax.set_yticks(valid_yticks)
                ax.set_yticklabels([f'{int(round(r))}' for r in rank_values])
                ax.tick_params(axis='y', labelsize=fs_cbar_ticks, right=True, left=False, length=2, width=0.8, pad=1)
                for label in ax.get_yticklabels():
                    label.set_horizontalalignment('left')
            else:
                ax.set_yticks([])
    last_row = n_ds

    def _add_slim_colorbar(gs_cell, metric):
        """Draw a colorbar over the observed range using the supplied normalization.

        Mark data endpoints and normalization landmarks. Use extension arrows
        when observations extend beyond the normalization bounds.
        """
        spec = color_specs[metric]
        norm_vmin = spec['vmin']
        norm_vmax = spec['vmax']
        norm_vcenter = spec.get('vcenter', None)
        extend = spec.get('extend', 'neither')
        data_min, data_max = data_range_dict.get(metric, (norm_vmin, norm_vmax))
        if norm_vcenter is not None:
            norm = TwoSlopeNorm(vmin=norm_vmin, vcenter=norm_vcenter, vmax=norm_vmax)
        else:
            norm = mpl.colors.Normalize(vmin=norm_vmin, vmax=norm_vmax)
        cmap = _get_cmap(spec['cmap'])
        host_ax = fig.add_subplot(gs_cell)
        host_ax.axis('off')
        bbox = host_ax.get_position()
        cbar_ax = fig.add_axes([bbox.x0, bbox.y0 + colorbar_y_offset, bbox.width, colorbar_height])
        need_extend_min = extend in ['min', 'both'] and data_min < norm_vmin
        need_extend_max = extend in ['max', 'both'] and data_max > norm_vmax
        arrow_width = colorbar_extend_ratio if need_extend_min or need_extend_max else 0
        display_min = data_min
        display_max = data_max
        n_segments = 256
        x_vals = np.linspace(display_min, display_max, n_segments + 1)
        cbar_ax.set_xlim(display_min, display_max)
        cbar_ax.set_ylim(0, 1)
        for k in range(n_segments):
            x0, x1 = (x_vals[k], x_vals[k + 1])
            mid_val = (x0 + x1) / 2
            color = cmap(norm(mid_val))
            rect = Rectangle((x0, 0), x1 - x0, 1, facecolor=color, edgecolor='none')
            cbar_ax.add_patch(rect)
        arrow_height = 1.0
        arrow_ratio = 0.15
        if need_extend_min:
            arrow_base = display_min + (display_max - display_min) * arrow_width
            arrow_tip = display_min
            arrow_mid = (arrow_base + arrow_tip) / 2
            arrow_color = cmap(norm(norm_vmin))
            arrow_pts = np.array([[arrow_tip, 0.5], [arrow_base, 0], [arrow_base, 1]])
            arrow_patch = Polygon(arrow_pts, closed=True, facecolor=arrow_color, edgecolor='none')
            cbar_ax.add_patch(arrow_patch)
        if need_extend_max:
            arrow_base = display_max - (display_max - display_min) * arrow_width
            arrow_tip = display_max
            arrow_color = cmap(norm(norm_vmax))
            arrow_pts = np.array([[arrow_tip, 0.5], [arrow_base, 0], [arrow_base, 1]])
            arrow_patch = Polygon(arrow_pts, closed=True, facecolor=arrow_color, edgecolor='none')
            cbar_ax.add_patch(arrow_patch)
        ticks = []
        tick_labels = []
        ticks.append(data_min)
        ticks.append(data_max)
        range_span = data_max - data_min
        if norm_vcenter is not None:
            if data_min <= norm_vcenter <= data_max:
                not_too_close = all((abs(norm_vcenter - t) > 0.1 * range_span for t in ticks))
                if not_too_close:
                    ticks.append(norm_vcenter)
        if data_min < norm_vmin < data_max:
            not_too_close = all((abs(norm_vmin - t) > 0.1 * range_span for t in ticks))
            if not_too_close:
                ticks.append(norm_vmin)
        if data_min < norm_vmax < data_max:
            not_too_close = all((abs(norm_vmax - t) > 0.1 * range_span for t in ticks))
            if not_too_close:
                ticks.append(norm_vmax)
        ticks = sorted(set(ticks))
        cbar_ax.set_xticks(ticks)
        cbar_ax.set_xticklabels([_format_tick(t) for t in ticks])
        cbar_ax.set_yticks([])
        cbar_ax.tick_params(axis='x', labelsize=fs_cbar_ticks, pad=2, rotation=90, labeltop=False, labelbottom=True, top=False, bottom=True, direction='out', length=3)
        for spine in ['top', 'left', 'right']:
            cbar_ax.spines[spine].set_visible(False)
        cbar_ax.spines['bottom'].set_visible(True)
        cbar_ax.spines['bottom'].set_linewidth(0.5)
        for label in cbar_ax.get_xticklabels():
            label.set_verticalalignment('top')
            label.set_horizontalalignment('center')
    for j, metric in enumerate(kde_metrics_pa):
        _add_slim_colorbar(gs[last_row, j], metric)
    for j in [4, 5]:
        ax_off = fig.add_subplot(gs[last_row, j])
        ax_off.axis('off')
    for j, metric in zip([6, 7], kde_metrics_st):
        _add_slim_colorbar(gs[last_row, j], metric)
    for j in [8, 9]:
        ax_off = fig.add_subplot(gs[last_row, j])
        ax_off.axis('off')
    fig.canvas.draw()
    ax0_pcc = axes_map[0, 0]
    ax0_rpa = axes_map[0, 5]
    ax0_ssim = axes_map[0, 6]
    ax0_rst = axes_map[0, 9]
    bb_pcc = ax0_pcc.get_position(fig)
    bb_rpa = ax0_rpa.get_position(fig)
    bb_ssim = ax0_ssim.get_position(fig)
    bb_rst = ax0_rst.get_position(fig)
    y_top = max(bb_pcc.y1, bb_rpa.y1, bb_ssim.y1, bb_rst.y1)
    band_h = 0.024
    pa_x0, pa_x1 = (bb_pcc.x0, bb_rpa.x1)
    fig.patches.append(Rectangle((pa_x0, y_top + title_gap + title_offset), pa_x1 - pa_x0, band_h, facecolor=color_specs.get('PA_TITLE_BAND', '#EAF1FF'), edgecolor='black', linewidth=0.8, zorder=50))
    fig.text((pa_x0 + pa_x1) / 2, y_top + title_gap + title_offset + band_h / 2, 'Prediction Accuracy', ha='center', va='center', fontsize=fs_group_title, weight='bold', zorder=51)
    st_x0, st_x1 = (bb_ssim.x0, bb_rst.x1)
    fig.patches.append(Rectangle((st_x0, y_top + title_gap + title_offset), st_x1 - st_x0, band_h, facecolor=color_specs.get('ST_TITLE_BAND', '#EAFEEA'), edgecolor='black', linewidth=0.8, zorder=50))
    fig.text((st_x0 + st_x1) / 2, y_top + title_gap + title_offset + band_h / 2, 'Structural Accuracy', ha='center', va='center', fontsize=fs_group_title, weight='bold', zorder=51)
    if plot_name:
        out_dir = os.path.dirname(plot_name)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        fig.savefig(plot_name, dpi=300, bbox_inches='tight', pad_inches=0.1)
        print(f'Figure saved to: {plot_name}')
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig

