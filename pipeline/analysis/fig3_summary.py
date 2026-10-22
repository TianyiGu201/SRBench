"""Summarize stratified gene-property effects across datasets."""
from __future__ import annotations
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D


METHOD_ORDER = ['XFuse', 'TESLA', 'iSTAR', 'scstGCN', 'SuperST_t', 'SuperST_t_g']
LEVEL_ORDER = ['Low', 'Medium', 'High']

def load_all_csvs(input_dir: str) -> pd.DataFrame:
    paths = sorted(glob.glob(os.path.join(input_dir, '*.csv')))
    if not paths:
        raise ValueError(f'No CSV files found in {input_dir}')
    dfs = [pd.read_csv(p) for p in paths]
    return pd.concat(dfs, ignore_index=True)

def zscore(x: pd.Series) -> np.ndarray:
    x = x.astype(float)
    std = x.std(ddof=0)
    if std == 0 or np.isnan(std):
        return np.zeros(len(x), dtype=float)
    return ((x - x.mean()) / std).to_numpy(dtype=float)

def build_task_scores(df: pd.DataFrame) -> pd.DataFrame:
    """
    Prediction_Score = mean(z(PCC), z(SCC), z(-KLD), z(-JSD))
    Structure_Score  = mean(z(SSIM), z(PSNR))
    Z-scoring is within each dataset–method.
    """
    required_cols = {'dataset', 'method', 'Gene_Expression_Level', 'Spatial_Variability_Level'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f'Missing required columns: {missing}')
    out = df.copy()
    out['Gene_expression_level'] = pd.Categorical(out['Gene_Expression_Level'], categories=LEVEL_ORDER, ordered=True)
    out['Spatial_variability_level'] = pd.Categorical(out['Spatial_Variability_Level'], categories=LEVEL_ORDER, ordered=True)

    def per_group(g: pd.DataFrame) -> pd.DataFrame:
        g = g.copy()
        pred_parts = []
        for c in ['PCC', 'SCC']:
            if c in g.columns:
                pred_parts.append(zscore(g[c]))
        for c in ['KLD', 'JSD']:
            if c in g.columns:
                pred_parts.append(zscore(-g[c]))
        g['Prediction_Score'] = np.mean(np.vstack(pred_parts), axis=0) if pred_parts else np.nan
        struct_parts = []
        for c in ['SSIM', 'PSNR']:
            if c in g.columns:
                struct_parts.append(zscore(g[c]))
        g['Structure_Score'] = np.mean(np.vstack(struct_parts), axis=0) if struct_parts else np.nan
        return g
    return pd.concat([per_group(g) for _, g in out.groupby(['dataset', 'method'], sort=False)], ignore_index=True)

def compute_group_means(df_scores: pd.DataFrame, score_col: str) -> pd.DataFrame:
    g = df_scores.dropna(subset=[score_col]).groupby(['dataset', 'method', 'Gene_expression_level', 'Spatial_variability_level'], observed=True).agg(mean_score=(score_col, 'mean'), n_genes=(score_col, 'size')).reset_index()
    return g

def marginal_contrasts(group_means: pd.DataFrame) -> pd.DataFrame:
    """
    ΔExpression abundance = mean(score | High) - mean(score | Low), averaged across spatial variability levels
    ΔSpatial variability  = mean(score | High) - mean(score | Low), averaged across expression levels
    Missing cells ignored.
    """
    rows = []
    for (dataset, method), g in group_means.groupby(['dataset', 'method'], observed=True):
        pivot = g.pivot(index='Spatial_variability_level', columns='Gene_expression_level', values='mean_score')
        pivot = pivot.reindex(index=LEVEL_ORDER, columns=LEVEL_ORDER)
        gel_marg = pivot.mean(axis=0, skipna=True)
        svl_marg = pivot.mean(axis=1, skipna=True)
        d_gel = np.nan
        d_svl = np.nan
        if pd.notna(gel_marg.get('High')) and pd.notna(gel_marg.get('Low')):
            d_gel = float(gel_marg['High'] - gel_marg['Low'])
        if pd.notna(svl_marg.get('High')) and pd.notna(svl_marg.get('Low')):
            d_svl = float(svl_marg['High'] - svl_marg['Low'])
        rows.append({'dataset': dataset, 'method': method, 'Delta_gene_expression_level': d_gel, 'Delta_spatial_variability': d_svl, 'n_cells_available': int(pivot.notna().sum().sum())})
    return pd.DataFrame(rows)

def _method_colors(methods):
    cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', [f'C{i}' for i in range(10)])
    return {m: cycle[i % len(cycle)] for i, m in enumerate(methods)}

def _set_limits_focus_main_quadrant(ax, x, y, q_low=0, q_high=100, pad_frac=0.08):
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if x.size == 0 or y.size == 0:
        ax.set_xlim((-1, 1))
        ax.set_ylim((-1, 1))
        return
    x_lo, x_hi = np.percentile(x, [q_low, q_high])
    y_lo, y_hi = np.percentile(y, [q_low, q_high])
    x_lo = min(x_lo, 0.0)
    x_hi = max(x_hi, 0.0)
    y_lo = min(y_lo, 0.0)
    y_hi = max(y_hi, 0.0)
    x_pad = (x_hi - x_lo) * pad_frac if x_hi > x_lo else 1.0 * pad_frac
    y_pad = (y_hi - y_lo) * pad_frac if y_hi > y_lo else 1.0 * pad_frac
    ax.set_xlim((x_lo - x_pad, x_hi + x_pad))
    ax.set_ylim((y_lo - y_pad, y_hi + y_pad))

def scatter_with_method_mean_and_iqr(ax, df_task, method_order, focus_quantiles=(0, 100)):
    """
    Small points: dataset–method
    Large point + x/y error bars: across-dataset mean ± IQR (25–75%) for each method
    (Legend handled externally.)
    """
    df_task = df_task[df_task['method'].isin(method_order)].copy()
    df_task['method'] = pd.Categorical(df_task['method'], categories=method_order, ordered=True)
    df_task = df_task.sort_values('method')
    colors = _method_colors(method_order)
    for m in method_order:
        sub = df_task[df_task['method'] == m]
        if sub.empty:
            continue
        x = sub['Delta_gene_expression_level'].to_numpy(float)
        y = sub['Delta_spatial_variability'].to_numpy(float)
        ax.scatter(x, y, s=10, alpha=0.55, linewidths=0, color=colors[m])
        x_fin = x[np.isfinite(x)]
        y_fin = y[np.isfinite(y)]
        if x_fin.size == 0 or y_fin.size == 0:
            continue
        mean_x = float(np.mean(x_fin))
        mean_y = float(np.mean(y_fin))
        q1_x, q3_x = (float(np.percentile(x_fin, 25)), float(np.percentile(x_fin, 75)))
        q1_y, q3_y = (float(np.percentile(y_fin, 25)), float(np.percentile(y_fin, 75)))
        ax.errorbar(mean_x, mean_y, xerr=[[mean_x - q1_x], [q3_x - mean_x]], yerr=[[mean_y - q1_y], [q3_y - mean_y]], fmt='o', markersize=6, markerfacecolor=colors[m], markeredgecolor='black', markeredgewidth=0.5, ecolor=colors[m], elinewidth=0.8, capsize=2.5, capthick=0.8, zorder=6)
    x_all = df_task['Delta_gene_expression_level'].to_numpy(float)
    y_all = df_task['Delta_spatial_variability'].to_numpy(float)
    ql, qh = focus_quantiles
    _set_limits_focus_main_quadrant(ax, x_all, y_all, q_low=ql, q_high=qh, pad_frac=0.08)
    ax.set_aspect('equal', adjustable='box')
    ax.set_box_aspect(1)
    ax.grid(True, linewidth=0.3)
    ax.axhline(0, linewidth=0.5)
    ax.axvline(0, linewidth=0.5)

def build_right_legend_handles(method_order):
    colors = _method_colors(method_order)
    method_handles = [Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=colors[m], markeredgecolor='none', markersize=4, label=m) for m in method_order]
    summary_handle = Line2D([0], [0], marker='o', linestyle='None', markerfacecolor='white', markeredgecolor='black', markersize=5, label='Across-dataset mean ± IQR')
    return method_handles + [summary_handle]

def make_one_page_pdf(df_summary, out_pdf):
    """
    Single panel (a) with two sub-axes:
      left: Prediction Accuracy
      right: Structural Accuracy
    One shared legend on the right.
    """
    width_mm = 174.0
    width_in = width_mm / 25.4
    height_in = 4.6
    pred = df_summary[df_summary['task'] == 'Prediction'].copy()
    st = df_summary[df_summary['task'] == 'Structure'].copy()
    colors = _method_colors(METHOD_ORDER)
    method_handles = [Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=colors[m], markeredgecolor='none', markersize=4, label=m) for m in METHOD_ORDER]
    summary_handle = Line2D([0], [0], marker='o', linestyle='None', markerfacecolor='white', markeredgecolor='black', markersize=5, label='Across-Dataset Mean and IQR')
    handles = method_handles + [summary_handle]
    with PdfPages(out_pdf) as pdf:
        fig = plt.figure(figsize=(width_in, height_in), constrained_layout=True)
        gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.55], wspace=0.15)
        ax_left = fig.add_subplot(gs[0, 0])
        ax_right = fig.add_subplot(gs[0, 1])
        ax_leg = fig.add_subplot(gs[0, 2])
        ax_leg.axis('off')
        fig.suptitle('a  Sensitivity To Gene Characteristics Across Datasets', y=1.02)
        ax_left.set_title('Prediction Accuracy', pad=2)
        ax_right.set_title('Structural Accuracy', pad=2)
        scatter_with_method_mean_and_iqr(ax_left, pred, METHOD_ORDER, focus_quantiles=(0, 100))
        scatter_with_method_mean_and_iqr(ax_right, st, METHOD_ORDER, focus_quantiles=(0, 100))
        for ax in (ax_left, ax_right):
            ax.set_xlabel('Marginal Effect Of Expression Level (High − Low)')
            ax.set_ylabel('Marginal Effect Of Spatial Variability (High − Low)')
            ax.tick_params(axis='both', which='major', pad=1)
        leg = ax_leg.legend(handles=handles, title='Methods', loc='center', frameon=True, facecolor='white', edgecolor='black', framealpha=1.0, ncol=1, borderpad=0.35, labelspacing=0.3, handletextpad=0.5)
        leg.get_frame().set_linewidth(0.6)
        if leg.get_title() is not None:
            leg.get_title().set_fontsize(5)
        pdf.savefig(fig, dpi=300)
        plt.close(fig)

def run_overall_summary(input_dir: str, out_summary_csv: str, out_pdf: str):
    df = load_all_csvs(input_dir)
    df = build_task_scores(df)
    gm_pred = compute_group_means(df, 'Prediction_Score')
    contr_pred = marginal_contrasts(gm_pred)
    contr_pred['task'] = 'Prediction'
    gm_struct = compute_group_means(df, 'Structure_Score')
    contr_struct = marginal_contrasts(gm_struct)
    contr_struct['task'] = 'Structure'
    summary = pd.concat([contr_pred, contr_struct], ignore_index=True)
    summary['method'] = pd.Categorical(summary['method'], categories=METHOD_ORDER, ordered=True)
    summary = summary.sort_values(['task', 'method', 'dataset']).reset_index(drop=True)
    summary.to_csv(out_summary_csv, index=False)
    make_one_page_pdf(summary[['dataset', 'method', 'task', 'Delta_gene_expression_level', 'Delta_spatial_variability']], out_pdf)
    print('Saved summary CSV:', out_summary_csv)
    print('Saved PDF figure:', out_pdf)

